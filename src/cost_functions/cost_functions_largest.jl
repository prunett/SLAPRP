module cost_functions_largest

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
    println("MIDPOINT routing policy")
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
    first_aisle = loc_per_aisle[1]
    sort!(first_aisle, by=l->l.id)
    # we look at the case where there is only one aisle
    if length(loc_per_aisle) == 1
        cost = 2*ins.wa*(first_aisle[1].aisle-1)
        cost += 2*ins.wc*first_aisle[end].col
        return route(o,first_aisle,cost)
    else
        # now we sort the last aisle
        last_aisle = loc_per_aisle[end]
        sort!(last_aisle, by=l->l.id, rev = true)
        # we initialize the cost with the first and last aisles crossed entirely
        cost = 2*ins.wa*(last_aisle[1].aisle-1)
        cost += 2*ins.wc*(ins.cmax + 1)
        # we iterate on the middle aisles
        locs_from_top = loc[]
        locs_from_bottom = loc[]
        for i in 2:length(loc_per_aisle)-1
            # we sort the loc of the aisle we are looking at
            sort!(loc_per_aisle[i], by=l->l.id)
            # we iterate to find the largest gap
            last_ind_from_bottom = 0
            best_gap = loc_per_aisle[i][1].col
            for last_ind = 1:length(loc_per_aisle[i])
                # we compute the gap
                if last_ind == length(loc_per_aisle[i])
                    gap = ins.cmax + 1 - loc_per_aisle[i][last_ind].col
                else
                    gap = loc_per_aisle[i][last_ind+1].col - loc_per_aisle[i][last_ind].col
                end
                @assert gap >= 0
                # we register if it is the largest gap
                if gap > best_gap
                    best_gap = gap
                    last_ind_from_bottom = last_ind
                end
            end
            # we update the cost
            cost += 2*ins.wc*(ins.cmax+1 - best_gap)
            # we add the locs in the corresponding arrays
            locs_from_bottom = append!(loc_per_aisle[i][1:last_ind_from_bottom],locs_from_bottom)
            append!(locs_from_top, reverse(loc_per_aisle[i][last_ind_from_bottom+1:end]))
        end
        # now we assemble the route
        route_loc = first_aisle
        append!(route_loc,locs_from_top)
        append!(route_loc,last_aisle)
        append!(route_loc,locs_from_bottom)
        
        return route(o, route_loc, cost)
    end
end

    






#=
function sol_cost(ins::instance,stol::Dict{sku,loc}, partial::Bool = false)
#Computes the cost of a solution, using the return policy
    cost = 0
    for o in ins.orderList
        cost += route_cost(ins,stol,o,partial).cost
    end
    return cost
end=#

# New methods
#sol_cost(ins::instance,assignment::assignment) = sol_cost(ins,assignment.stol)

function convert_x_to_dict(ins::instance,x::Array{Int,2})
# This function convert a matrix like the one used for the x variable in pricing plan, to a dist
    stol = Dict{sku,loc}()
    for s in ins.skuList
        for l in ins.locList
            if x[l.id,s.id] == 1
                stol[s] = l
                break
            end
        end
    end
    return stol
end
#=
function create_solution(ins::instance,stol::Dict{sku,loc}, partial::Bool = false)
# Creates a solution type

    routes = Array{route,1}()
    for o in ins.orderList
        push!(routes, route_cost(ins,stol,o,partial))
    end

    return solution(sol_cost(ins,stol,partial),routes,stol)
end=#







end # end of module