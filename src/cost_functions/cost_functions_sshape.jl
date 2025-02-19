module cost_functions_sshape

#using Concorde
using Combinatorics

using ..structures_decomp

#export dist, route_cost, sol_cost, create_solution

################################################################
# IMPORTANT INFO

# If we want to improve this part of the algo, read the literature review
# of the paper on inventory routing in warehouses
# They provide ref for a polynomial algo (via DP) and a MIP formulation
################################################################


#=
isequal(a::loc,b::loc) =
# defines a method for isequal with locations
    if (a.aisle == b.aisle) && (a.col == b.col)
        return true
    else
        return false
    end
end=#

function print_policy()
    println("S-SHAPE routing policy")
end

function dist_locs_by_top(a::loc,b::loc,ins::instance)
    if a.aisle == b.aisle
        return ins.wc*abs(a.col - b.col)
    else
        return ins.wa*abs(a.aisle - b.aisle) + ins.wc*(2*(ins.cmax + 1) - a.col - b.col)
    end
end

function dist_locs_by_bottom(a::loc,b::loc,ins::instance)
    if a.aisle == b.aisle
        return ins.wc*abs(a.col - b.col)
    else
        return ins.wa*(a.aisle - b.aisle) + ins.wc*(a.col + b.col)
    end
end



# This function returns the loc group in the walked order for the midpoint policy
function route_cost(ins::instance,loc_group::Vector{loc},o::order,partial::Bool = false)
    # test if the length matches the order
    if !partial && (length(loc_group) != length(o.sku))
        throw(DimensionMismatch("the length of the arguments are not matching!"))
    end
    if isempty(loc_group)
        return route(o,loc_group,0)
    end

    loc_per_aisle = [[] for a=1:ins.amax]
    for l in loc_group
        push!(loc_per_aisle[l.aisle],l)
    end
    # we remove the empty aisles
    filter!(!isempty,loc_per_aisle)

    # first we compute the cost
    cost = 2*ins.wa*(loc_per_aisle[end][1].aisle - 1)
    if iseven(length(loc_per_aisle))
        # in this case all aisles are completely crossed
        cost += ins.wc*(ins.cmax + 1)*length(loc_per_aisle)
    else
        cost += ins.wc*(ins.cmax + 1)*(length(loc_per_aisle) - 1)
        # return for the last aisle
        cost += 2*ins.wc*maximum([l.col for l in loc_per_aisle[end]])
    end

    # now we do the routing
    route_loc = loc[]
    for i in eachindex(loc_per_aisle)
        if isodd(i)
            sort!(loc_per_aisle[i], by= l-> l.id)
        else
            sort!(loc_per_aisle[i], by= l->l.id, rev=true)
        end
        append!(route_loc,loc_per_aisle[i])
    end

    return route(o,route_loc,cost)
end

    
end # end of module