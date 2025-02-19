
#export OG_cut, separate_OG_all

mutable struct OG_cut <: LocGroupCut
    loc_array::Array{loc,1}
    order::order
    sku::sku
    constraint::Union{Nothing,ConstraintRef}
    active::Bool
    iter_without_use::Int
    dual_value::Union{Nothing,Float64}
    touched_routes_index::Array{Int,1}

    function OG_cut(
        loc_array::Array{loc,1},
        order::order,
        sku::sku,
        constraint = nothing,
        active = false,
        iter_without_use = 0,
        dual_value = nothing;
        touched_routes_index = Int[]) # the index of the touched routes in rho_val or a_coef
        return new(loc_array,order,sku,constraint,active,iter_without_use,dual_value,touched_routes_index)
    end
end

function is_cut_equal(ins::instance,cut::OG_cut,rho_val_o::Array{Float64,1},xi_val::Array{Float64,2},a_coef_o::Array{route,1},EPS::Float64)
# This function checks if the cut is equal
    # check if active
    @assert(cut.active)
    # compute lhs
    lhs = 0.
    for r in cut.touched_routes_index
        lhs += rho_val_o[r]
    end

    return abs(lhs - sum(xi_val[l.id,cut.sku.id] for l in cut.loc_array)) < EPS
end

function is_cut_violated(ins::instance,cut::OG_cut,rho_val_o::Array{Float64,1},xi_val::Array{Float64,2},a_coef_o::Array{route,1},EPS::Float64)
# This function checks if the cut is violated
    # check if unactive
    @assert(!cut.active)
    # compute lhs
    lhs = 0.
    for r in cut.touched_routes_index
        lhs += rho_val_o[r]
    end

    val = lhs - sum(xi_val[l.id,cut.sku.id] for l in cut.loc_array)

    return (val < - EPS), val
end

function add_route_to_cut!(cut::OG_cut,master_problem::Model,touched_constraints::Array{ConstraintRef,1},constraints_coef::Array{Float64,1})
# This function apply the good coefficient
    @assert(cut.active)
    #add the good coefficient
    push!(touched_constraints, cut.constraint)
    push!(constraints_coef, -1)
end

function add_route_to_cut!(r::route,cut_dict::Dict{order,Array{OG_cut,1}},master_problem::Model,touched_constraints::Array{ConstraintRef,1},constraint_coef::Array{Float64,1},route_index::Int)
# add the route to all active cut
    for cut in cut_dict[r.order]
        if intersects_id(r.loc,cut.loc_array)
            # first add the route to the touched routes of the cut 
            push!(cut.touched_routes_index, route_index)
            # then add the route to the constraint of the cut if it is active
            if cut.active
                add_route_to_cut!(cut,master_problem,touched_constraints,constraint_coef)
            end
        end
    end
end



function activate_cut!(ins::instance,master_problem::Model,cut::OG_cut,rho_vars::Dict{order,Array{VariableRef,1}},a_coef,xi_vars::Array{VariableRef,2})
# activate the cut
    cut.active = true
    cut.iter_without_use = 0
    cut.constraint = @constraint(master_problem, sum(xi_vars[l.id, cut.sku.id] for l in cut.loc_array) <= sum(rho_vars[cut.order][r] * intersects_id(a_coef[cut.order][r].loc, cut.loc_array) for r=1:length(rho_vars[cut.order])))
    # name
    set_name(cut.constraint, "OG_cut[$(cut.order.id) ; $(cut.sku.id) ; $([l.id for l in cut.loc_array])]")
end


function initialization_cuts_OG!(ins::instance,cut_pool,config)
# Initialize the empty cut pool for OG
    cut_pool["OG_cut"] = Dict{order,Array{OG_cut,1}}(o => OG_cut[] for o in ins.orderList)
end

function initiate_aisle_cuts!(ins::instance,cut_pool,a_coef)
    # we generate in the pool all the cuts that on an entire aisle
    for a = 1:ins.amax
        # get all the locs of the aisle
        aisle = filter(l -> l.aisle == a, ins.locList)
        for o in ins.orderList
            for s in o.sku
                touched_routes_index = compute_touched_routes_index(ins, o, aisle, a_coef)
                push!(cut_pool["OG_cut"][o], OG_cut(aisle,o,s,touched_routes_index = touched_routes_index))
            end
        end
    end
end


function separate_OG_s(ins::instance,o::order,s::sku,xi_val::Array{Float64,2},rho_val::Dict{order,Array{Float64,1}},a_coef::Dict{order,Array{route,1}},EPS::Float64)
# This function separates the OG cuts for a given s. It is an intermediary function that prepares input to call
# the recursive one that do the actual job.

    # First we start by ordering locList by decreasing order of xi
    locList = copy(ins.locList)
    sort!(locList, by = l -> xi_val[l.id,s.id], rev=true)
    # Create the list of routes available for the cut
    route_av_list = Array{route,1}()
    route_av_value = Float64[]
    for r = 1:length(a_coef[o])
        push!(route_av_list, a_coef[o][r])
        push!(route_av_value, rho_val[o][r])
    end
    # Sort route_av_list in a useful order:
    # not implemented for now

    # Create a dict that indicates, from the id of a loc, its position in locList
    loc_position = Dict{Int,Int}(locList[i].id => i for i=1:length(locList))

    # Create delta array
    delta_xi = Float64[xi_val[l.id,s.id] for l in ins.locList]

    # Call the recursive procedure
    return separate_OG_rec(locList, route_av_list, route_av_value, 0., delta_xi, loc_position, EPS)

end


function separate_OG_all!(ins::instance, cut_pool, xi_val::Array{Float64,2}, rho_val::Dict{order,Array{Float64,1}},a_coef::Dict{order, Array{route,1}},EPS_cut_sep::Vector{Float64})
# Main OG cuts separation procedure
# Maybe need to have the output in a better organization/structure

    violated_OG_cut = OG_cut[]
    violated_cuts_array = Array{Vector{Vector{loc}},2}(undef, length(ins.orderList), length(ins.skuList))
    fill!(violated_cuts_array, [])
    cut_value = Float64[]
    cut_value_array = Array{Vector{Float64},2}(undef, length(ins.orderList), length(ins.skuList))
    fill!(cut_value_array, [])
    added_cuts = false

    order_sku_list = []
    for o in ins.orderList, s in o.sku
        push!(order_sku_list,(o,s))
    end
    for (o,s) in order_sku_list
        cut_list_s, cut_value_s = separate_OG_s(ins,o,s,xi_val,rho_val,a_coef,EPS_cut_sep[o.id])
        violated_cuts_array[o.id,s.id] = cut_list_s
        cut_value_array[o.id,s.id] = cut_value_s
    end
    # post processing
    for o in ins.orderList, s in o.sku, lgind in eachindex(violated_cuts_array[o.id,s.id])
        # compute the touched routes
        touched_routes_index = compute_touched_routes_index(ins, o, violated_cuts_array[o.id,s.id][lgind], a_coef)


        push!(cut_pool["OG_cut"][o], OG_cut(violated_cuts_array[o.id,s.id][lgind],o,s,touched_routes_index = touched_routes_index))
        added_cuts = true
        push!(violated_OG_cut, cut_pool["OG_cut"][o][end])
        push!(cut_value, cut_value_array[o.id,s.id][lgind])
    end

    return violated_OG_cut, cut_value, added_cuts
end

function compute_touched_routes_index(ins::instance, o::order, cut_loc_array::Array{loc,1}, a_coef::Dict{order,Array{route,1}})
# This function computes the list of routes touched by the new generated cut
    touched_routes_index = Int[]
    # iterate over all routes for order o, and register the index of the touched routes
    for r=1:length(a_coef[o])
        if intersects_id(cut_loc_array,a_coef[o][r].loc)
            push!(touched_routes_index, r)
        end
    end

    return touched_routes_index
end


function separate_OG_rec(locList::Array{loc,1},route_av_list::Array{route,1},route_av_value::Array{Float64,1},value_current::Float64,delta_xi::Array{Float64,1},loc_position::Dict{Int,Int},EPS::Float64)
# This function is the recurrence for the separation fo OG cuts. The value_current corresponds to the value
# of the lhs of the inequality: sum (rho) - sum(xi). A cut is violated when this value is < 0

    # Termination criteria
    #if loc list empty
    if isempty(locList)
        return Array{loc,1}[[]], Float64[0.]
    end
    # if xi potential too low to possibly get a violated cut
    if (sum(delta_xi[l.id] for l in locList) <= value_current) && (value_current >= - EPS)
        return Array{loc,1}[[]], Float64[0.]
    end
    # if all remaining xi are equal to 0
    if sum(delta_xi[l.id] for l in locList) <= EPS
        return Array{loc,1}[[]], Float64[0.]
    end
    # if route_av_list is empty or all rho = 0 -> add all nonnegative xi to the cut
    if isempty(route_av_list) || (sum(route_av_value) <= EPS)
        a = Array{loc,1}[[]]
        val = 0.
        for l in locList
            if delta_xi[l.id] >= EPS
                val += delta_xi[l.id]
                push!(a[1], l)
            end
        end
        return a, Float64[- val]
    end


    #=
    # if route_av_list is empty -> add all the remaining locs to the cut
    if isempty(route_av_list)
        return Array{loc,1}[[]], Float64[0.]
        #return Array{loc,1}[locList], Float64[(isempty(locList) ? 0. : - sum(delta_xi[l.id] for l in locList))]
    end
    # this condition not sure about it
    if sum(route_av_value) < EPS
        return Array{loc,1}[[]], Float64[0.]
    end=#


    # Pop the first element of locList
    l = locList[1]
    locList_new = copy(locList[2:end])

    # First branch: without l
    # pop elements of route_av_list and route_av_value that do not contain anymore elemnts of locList
    route_av_list_1 = copy(route_av_list)
    route_av_value_1 = copy(route_av_value)
    current_position = loc_position[l.id]
    ind = 1
    while ind <= length(route_av_list_1)
        if maximum(loc_position[i.id] for i in route_av_list_1[ind].loc) <= current_position
            deleteat!(route_av_list_1, ind)
            deleteat!(route_av_value_1, ind)
            ind -= 1
        end
        ind += 1
    end
    # Call the recursion with the new parameters
    cut_list_1, cut_value_1 = separate_OG_rec(locList_new, route_av_list_1, route_av_value_1, value_current, delta_xi, loc_position, EPS)
    # Filter cuts where the value is >= 0
    ind = 1
    while ind <= length(cut_list_1)
        if cut_value_1[ind] + value_current >= - EPS
            deleteat!(cut_list_1,ind)
            deleteat!(cut_value_1, ind)
            ind -= 1
        end
        ind += 1
    end

    # Second branch: with l
    delta_value = - delta_xi[l.id]
    # pop elements of route_av_list that stop in l + compute delta_value at the same time
    route_av_list_2 = copy(route_av_list)
    route_av_value_2 = copy(route_av_value)
    ind = 1
    while ind <= length(route_av_list_2)
        if (l in route_av_list_2[ind].loc)
            delta_value += route_av_value_2[ind]
            deleteat!(route_av_list_2,ind)
            deleteat!(route_av_value_2, ind)
            ind -= 1
        end
        ind += 1
    end
    # Call the recursion with the new parameters
    cut_list_2, cut_value_2 = separate_OG_rec(locList_new,route_av_list_2,route_av_value_2,value_current + delta_value,delta_xi,loc_position,EPS)
    # Add the value of the addition of l to the cuts
    cut_value_2 .+= delta_value
    # Add l to the cut_list_2
    for c in cut_list_2
        push!(c,l)
    end
    # Filter cuts where the value is >= 0
    ind = 1
    while ind <= length(cut_list_2)
        if cut_value_2[ind] + value_current >= - EPS
            deleteat!(cut_list_2, ind)
            deleteat!(cut_value_2, ind)
            ind -= 1
        end
        ind += 1
    end

    # merge results
    #append!(cut_list_1, cut_list_2)
    #append!(cut_value_1, cut_value_2)

    # return result
    return vcat(cut_list_1,cut_list_2), vcat(cut_value_1,cut_value_2)
end
