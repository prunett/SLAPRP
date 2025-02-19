# Link Bis cuts

mutable struct LinkCut <:AbstractCut
    loc::loc
    order::order
    sku::sku
    constraint::Union{Nothing,ConstraintRef}
    active::Bool
    iter_without_use::Int
    dual_value::Union{Nothing,Float64}

    function LinkCut(
        loc::loc,
        order::order,
        sku::sku,
        constraint = nothing,
        active = false,
        iter_without_use = 0,
        dual_value = nothing)
        return new(loc,order,sku,constraint,active,iter_without_use,dual_value)
    end
end

function is_cut_equal(cut::LinkCut,rho_val_o::Array{Float64,1},xi_val::Array{Float64,2},a_coef_o::Array{route,1},EPS::Float64)
# This function checks if an active cut is at equality
    # check if active
    @assert(cut.active)
    #compute the lhs
    lhs = 0.
    for r in 1:length(a_coef_o)
        if cut.loc in a_coef_o[r].loc
            lhs += rho_val_o[r]
        end
    end
    # return the equality value
    return abs(lhs - xi_val[cut.loc.id,cut.sku.id]) < EPS
end

function is_cut_violated(cut::LinkCut,rho_val_o::Array{Float64,1},xi_val::Array{Float64,2},a_coef_o::Array{route,1},EPS::Float64)
# This function checks if an unactive cut is at violated
    # check if active
    @assert(!cut.active)
    #compute the lhs
    lhs = 0.
    for r in 1:length(a_coef_o)
        if cut.loc in a_coef_o[r].loc
            lhs += rho_val_o[r]
        end
    end
    # return the equality value
    return lhs - xi_val[cut.loc.id,cut.sku.id] < - EPS
end

function add_route_to_cut!(cut::LinkCut,master_problem::Model,touched_constraints::Array{ConstraintRef,1},constraints_coef::Array{Float64,1})
# This function add the route to the cut
    # assert the cut is active
    @assert(cut.active)
    # add the good coefficient
    push!(touched_constraints,cut.constraint)
    push!(constraints_coef, -1)
end

function add_route_to_cut!(r::route,cut_array::Array{LinkCut,3},master_problem::Model,touched_constraints::Array{ConstraintRef,1},constraint_coef::Array{Float64,1})
# add the route to all active cut
    for l in r.loc
        for s in r.order.sku
            # only if the cut is active
            if cut_array[l.id,r.order.id,s.id].active
                add_route_to_cut!(cut_array[l.id,r.order.id,s.id],master_problem,touched_constraints,constraint_coef)
            end
        end
    end
end

function activate_cut!(master_problem::Model,cut::LinkCut,rho_vars::Dict{order,Array{VariableRef,1}},a_coef,xi_vars::Array{VariableRef,2})
# activate the cut
    cut.active = true
    cut.iter_without_use = 0
    cut.constraint = @constraint(master_problem, xi_vars[cut.loc.id,cut.sku.id] <= sum((cut.loc in a_coef[cut.order][r].loc ? rho_vars[cut.order][r] : 0) for r=1:length(rho_vars[cut.order])))
    #name
    set_name(cut.constraint,"LinkCut[$(cut.loc.id), $(cut.order.id), $(cut.sku.id)]")
end

function activate_all_link_cut!(ins::instance,master_problem::Model,cut_pool,rho_vars::Dict{order,Array{VariableRef,1}},a_coef,xi_vars::Array{VariableRef,2})
# This function activate all link cut
    for l in ins.locList
        for o in ins.orderList
            for s in o.sku
                if !(s in keys(ins.fixed_positions))
                    activate_cut!(master_problem,cut_pool["LinkCut"][l.id,o.id,s.id],rho_vars,a_coef,xi_vars)
                end
            end
        end
    end
end


function separate_cut_Link()
end

function initialization_cuts_Link!(ins::instance,cut_pool)
# Initialize the cut pool for Link cuts
    cut_pool["LinkCut"] = Array{LinkCut,3}(undef,length(ins.locList),length(ins.orderList),ins.n)
    for l in ins.locList
        for o in ins.orderList
            for s in o.sku
                cut_pool["LinkCut"][l.id,o.id,s.id] = LinkCut(l,o,s)
            end
        end
    end
end
