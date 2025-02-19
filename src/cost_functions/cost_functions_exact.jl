module cost_functions_exact

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
    println("EXACT routing policy")
end

# The distance between 2 locs
function dist(a::loc,b::loc,ins::instance)
    if a.aisle == b.aisle
        return abs(a.col-b.col)*ins.wc
    else
        return ins.wc * min((a.col + b.col) , (2*(ins.cmax + 1) - a.col - b.col)) + abs(a.aisle - b.aisle)*ins.wa
    end
end

dist(ins::instance,a::loc,b::loc) = dist(a,b,ins)


######################
# Version with concorde
#######################
#=
function route_cost(ins::instance,loc_group::Array{loc,1},o::order, partial::Bool = false)
    # test if the length matches the order
    if !partial && (length(loc_group) != length(o.sku))
        throw(DimensionMismatch("the length of the arguments are not matching!"))
    end
    # Test if the length of the loc_group is 1, in that case concorde will not work and we need to do it by hand
    @assert !isempty(loc_group)
    if length(loc_group) == 1
        cost = 2*dist(ins,ioloc,loc_group[1])
        return route(o,loc_group,cost)
    end
    # build the distance matrix
    loc_with_io = copy(loc_group)
    pushfirst!(loc_with_io, ioloc)
    dist_mat = Int[dist(ins,loc_with_io[i],loc_with_io[j]) for i=1:length(loc_with_io), j=1:length(loc_with_io)]
    # call concorde to solve the tsp
    opt_tour, opt_len = solve_tsp(dist_mat)

    # sort the locs to be in the good order
    loc_with_io = loc_with_io[opt_tour]
    # assert that we indeed pop the ioloc
    @assert loc_with_io[1] == ioloc
    popfirst!(loc_with_io)

    return route(o,loc_with_io,opt_len)
end
=#

######################################## 
# Version with Combinatorics
########################################
function route_cost(ins::instance,loc_group::Array{loc,1},o::order, partial::Bool = false)
    # test if the length matches the order
    if !partial && (length(loc_group) != length(o.sku))
        throw(DimensionMismatch("the length of the arguments are not matching!"))
    end

    # iterate on the permutations of loc_group
    best_cost = typemax(Int)
    best_loc_group = copy(loc_group)
    for loc_group_test in permutations(loc_group)
        cost = 0
        cost += dist(ioloc,loc_group_test[1],ins)
        for i=1:length(loc_group_test)-1
            cost += dist(loc_group_test[i], loc_group_test[i+1],ins)
        end
        cost += dist(loc_group_test[end],ioloc,ins)
        if cost < best_cost
            best_cost = cost
            best_loc_group = copy(loc_group_test)
        end
    end

    return route(o,best_loc_group,best_cost)
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




## Some insertion and delta evaluation function
function delta_insertion_route(ins::instance,route::route,loc::loc)
# Compute the delta insertion cost to insert a stop in loc into the route
    if isempty(route.loc)
        return 2*dist(ins,ioloc,loc)
    end
    # Find the insertion position in the route
    #=i = 1
    while i <= length(route.loc)
        if route.loc[i] > loc
            break
        end
        i += 1
    end=#
    # Now we compute the insertion cost
    delta_best = Inf

    for i=1:length(route.loc)+1
        delta = 0
        if i == 1
            delta += dist(ins,ioloc,loc) + dist(ins,loc,route.loc[1]) - dist(ins,ioloc,route.loc[1])
        elseif i == length(route.loc) + 1
            delta += dist(ins,route.loc[end],loc) + dist(ins,loc,ioloc) - dist(ins,route.loc[end],ioloc)
        else
            delta += dist(ins,route.loc[i-1],loc) + dist(ins,loc,route.loc[i]) - dist(ins,route.loc[i-1],route.loc[i])
        end

        if delta < delta_best
            delta_best = delta
        end
    end

    return convert(Int,delta)
end

function delta_insertion_main(ins::instance, route_array::Array{route,1}, occupation::Dict{loc,Int}, s::sku, l::loc) ::Int
# compute the delta cost for insertion of sku s in loc l
    # is the loc is full then it is 0
    if occupation[l] == ins.n_sym
        return typemax(Int)
    end
    delta = 0
    for o in ins.sku_orderList[s]
        delta += delta_insertion_route(ins, route_array[o.id], l)
    end

    return delta
end





end # end of module
