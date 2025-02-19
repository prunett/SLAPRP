module cuts

using JuMP

using ..structures_decomp
using ..cost_functions
using ..BB

abstract type AbstractCut end
abstract type LocGroupCut <:AbstractCut end



include("cuts_og.jl")
include("cuts_LB.jl")

############### Cut data structure ##################
CutCollection = Dict{String,Union{Nothing, Array{LinkCut,3}, Dict{order,Array{OG_cut,1}}}} 


function deactivate_cut!(master_problem::Model,cut::AbstractCut)
# This function deactivate the cut, but does not affect the structures it is stocked in
    # delete cuts
    delete(master_problem, cut.constraint)
    # change its properties
    cut.active = false
    cut.constraint = nothing
    cut.dual_value = nothing
    cut.iter_without_use = 0
end


function update_active_cut_pool!(ins::instance,master_problem::Model,max_iter::Int,xi_vars::Array{VariableRef,2},xi_val::Array{Float64,2},rho_vars::Dict{order,Array{VariableRef,1}},rho_val::Dict{order,Array{Float64,1}},a_coef,cut_pool::CutCollection,EPS::Float64,EPS_cut_violated::Vector{Float64},EPS_vio_SL1::Float64,max_cut_per_order::Int)
# This function performs at each iteration of the CG
# It checks all cuts, for unactive ones it checks if they should be added, for active ones
# it increase counter and check if they should be deactivated
    # keep in memory if we added a cut
    cut_added = false
    # first LinkCut
    
    for l in ins.locList, o in ins.orderList, s in o.sku
        cut = cut_pool["LinkCut"][l.id,o.id,s.id]
        rho_val_o = rho_val[o]
        a_coef_o = a_coef[o]
        if cut.active
            # Update the iter without improvement
            if is_cut_equal(cut,rho_val_o,xi_val,a_coef_o,EPS)
                cut.iter_without_use = 0
            else
                cut.iter_without_use += 1
            end
        else
            # check if we should add it
            if is_cut_violated(cut,rho_val_o,xi_val,a_coef_o,EPS_vio_SL1)
                # add cut
                activate_cut!(master_problem,cut,rho_vars,a_coef,xi_vars)
                cut_added = true
            end
        end
    end
    
    # OG cuts
    for o in ins.orderList
        rho_val_o = rho_val[o]
        a_coef_o = a_coef[o]

        # new parameters for limited cut activation
        nb_active_cuts = 0
        violated_cuts = OG_cut[]
        violated_val = Float64[]

        for cut in cut_pool["OG_cut"][o]
            if cut.active
                nb_active_cuts += 1
                # update iter without imp
                if is_cut_equal(ins,cut,rho_val_o,xi_val,a_coef_o,EPS)
                    cut.iter_without_use = 0
                else
                    cut.iter_without_use += 1
                end
            else
                # check if it should be added
                violated, value = is_cut_violated(ins,cut,rho_val_o,xi_val,a_coef_o,EPS_cut_violated[o.id])
                if violated
                    # register the cut as violated
                    push!(violated_cuts, cut)
                    push!(violated_val, value)
                    # add cut
                    #activate_cut!(ins,master_problem,cut,rho_vars,a_coef,xi_vars)
                    #cut_added = true
                end
            end
        end

        # now we sort the violated cuts (if needed) and activate them
        if length(violated_cuts) + nb_active_cuts <= max_cut_per_order
            for cut in violated_cuts
                activate_cut!(ins,master_problem,cut,rho_vars,a_coef,xi_vars)
                cut_added = true
            end
        else
            # in this case we sort and only activate the best cuts
            perm = partialsortperm(violated_val, 1:(max_cut_per_order - nb_active_cuts))
            for cut in violated_cuts[perm]
                activate_cut!(ins,master_problem,cut,rho_vars,a_coef,xi_vars)
                cut_added = true
            end
        end
    end

    # return false if we added a cut
    return cut_added
end

function deactivate_cut_pool!(ins::instance,master_problem::Model,max_iter::Int,max_iter_SL1::Int,xi_vars::Array{VariableRef,2},xi_val::Array{Float64,2},rho_vars::Dict{order,Array{VariableRef,1}},rho_val::Dict{order,Array{Float64,1}},a_coef,cut_pool::CutCollection,EPS::Float64)
# This function performs at each iteration of the CG
# It deactivates cuts that should be removed
    # first LinkCut
    
    for l in ins.locList, o in ins.orderList, s in o.sku
        cut = cut_pool["LinkCut"][l.id,o.id,s.id]
        if cut.active
            # Deactivate the cut if needed
            if cut.iter_without_use > max_iter_SL1
                deactivate_cut!(master_problem,cut)
            end
        end
    end
    # OG cuts
    for o in ins.orderList, cut in cut_pool["OG_cut"][o]
        if cut.active
            if cut.iter_without_use > max_iter
                deactivate_cut!(master_problem,cut)
            end
        end
    end
end


function initialization_cuts(ins::instance,config)
# initialize cuts
    cut_pool = CutCollection()
    # LinkCut
    initialization_cuts_Link!(ins,cut_pool)
    # OG cuts
    initialization_cuts_OG!(ins,cut_pool,config)

    return cut_pool
end


function collect_dual_variable_cut!(ins::instance,master_problem::Model,cut_pool::CutCollection)
# This function collect the dual variables for the cuts
    #collect activated link cuts dual variables
    for l in ins.locList, o in ins.orderList, s in o.sku
        cut = cut_pool["LinkCut"][l.id,o.id,s.id]
        if cut.active
            cut.dual_value = - dual(cut.constraint)
        end
    end
    # collect OG cuts dual variables
    for o in ins.orderList, cut in cut_pool["OG_cut"][o]
        if cut.active
            cut.dual_value = - dual(cut.constraint)
        end
    end

end

function separate_cut_all!(ins::instance,cut_pool::CutCollection,xi_val::Array{Float64,2},rho_val::Dict{order,Array{Float64,1}},a_coef::Dict{order,Array{route,1}},EPS_cut_separation::Vector{Float64})
# This function separates all cuts
    # separate OG cuts
    violated_cuts, cut_value, added_OG = separate_OG_all!(ins, cut_pool, xi_val, rho_val, a_coef, EPS_cut_separation)

    return violated_cuts, cut_value, added_OG
end

function add_route_to_cut!(r::route,cut_pool::CutCollection,master_problem::Model,touched_constraints::Array{ConstraintRef,1},constraint_coef::Array{Float64,1},route_index::Int)
# This is the main function to add new routes with the active cuts
    add_route_to_cut!(r,cut_pool["LinkCut"],master_problem,touched_constraints,constraint_coef)
    add_route_to_cut!(r,cut_pool["OG_cut"],master_problem,touched_constraints,constraint_coef,route_index)
end

end # end of module
