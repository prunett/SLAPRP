module cost_functions_return

using JuMP, CPLEX

using ..structures_decomp

#export dist, route_cost, sol_cost, create_solution

# Defines a new isless method, to use the sorting functions on locations
Base.isless(a::loc,b::loc) = (a.aisle < b.aisle) || ((a.aisle == b.aisle) && (a.col < b.col))

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
    println("RETURN routing policy")
end

function dist(a::loc,b::loc,ins::instance)
# Computes the walking distance between loc a and b, with return policy
    if a.aisle == b.aisle
        return abs(a.col-b.col)*ins.wc
    else
        return (a.col + b.col)*ins.wc + abs(a.aisle - b.aisle)*ins.wa
    end
end

dist(ins::instance,a::loc,b::loc) = dist(a,b,ins)

#=
function route_cost(ins::instance, stol::Dict{sku,loc}, o::order, partial::Bool = false)
# This function computes the cost of a route with the given assignment, and return policy
# also works with partial routes and partial assignment


    skuList = intersect(collect(keys(stol)), o.sku)
    if !partial && (length(skuList) != length(o.sku))
        throw(DomainError("Length mismatch for non partial routes !"))
    end
    if isempty(skuList)
        return route(o,loc[],0)
    end
    #println("skuList = $skuList")
    #println("stol = $stol")
    routeLoc = [stol[s] for s in skuList]
    routeLoc_ordered = sort(routeLoc)

    cost = 0
    cost += dist(ioloc,routeLoc_ordered[1],ins)
    for i=1:length(routeLoc_ordered)-1
        cost += dist(routeLoc_ordered[i], routeLoc_ordered[i+1],ins)
    end
    cost += dist(routeLoc_ordered[end],ioloc,ins)

    return route(o,routeLoc_ordered,cost)
end=#

# New methods
#route_cost(ins::instance, assignment::assignment,o::order) = route_cost(ins,assignment.stol,o)
#route_cost(ins::instance, assignment::assignment, i::Int) = route_cost(ins,assignment,ins.orderDict[i])



#New method without assignment given
function route_cost(ins::instance,loc_group::Array{loc,1},o::order, partial::Bool = false)
    # test if the length matches the order
    if !partial && (length(loc_group) != length(o.sku))
        throw(DimensionMismatch("the length of the arguments are not matching!"))
    end
    # get cost
    loc_group_ordered = sort(loc_group)
    cost = 0
    cost += dist(ioloc,loc_group_ordered[1],ins)
    for i=1:length(loc_group_ordered)-1
        cost += dist(loc_group_ordered[i], loc_group_ordered[i+1],ins)
    end
    cost += dist(loc_group_ordered[end],ioloc,ins)

    return route(o,loc_group_ordered,cost)
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

# New methods
#create_solution(ins::instance,x::Array{Int,2}) = create_solution(ins, convert_x_to_dict(ins,x))

function create_pricing_problem_pure(ins::instance,o::order)
# This function creates the pricing problem for the return routing policy
    amax = ins.amax
    n = ins.n
    n_loc = length(ins.locList)
    locList = ins.locList
    wa = ins.wa
    wc = ins.wc
    loc_aisle = Array{Array{loc,1},1}()
    for a=1:amax
        push!(loc_aisle, filter(l -> l.aisle == a, locList))
    end

    pricing_route_o = Model(CPLEX.Optimizer)
    # Removing output from CPLEX
    set_optimizer_attribute(pricing_route_o, "CPX_PARAM_SCRIND", 0)
    # No objective for now
    @objective(pricing_route_o, Min, 0)
    # x = 1 if the route passess by l, 0 otherwise
    @variable(pricing_route_o, x[1:n_loc], Bin)
    # y variable
    @variable(pricing_route_o, y[1:amax] >= 0)
    # z variable
    @variable(pricing_route_o, z >= 0)
    # constraint on x
    @constraint(pricing_route_o, sum(x[l] for l=1:n_loc) == length(o.sku))
    # constraint on y
    @constraint(pricing_route_o, [a=1:amax,l=1:length(loc_aisle[a])], y[a] >= loc_aisle[a][l].col*wc* x[loc_aisle[a][l].id])
    # constraint on z
    @constraint(pricing_route_o, [a=1:amax,l=1:length(loc_aisle[a])], z >= (a-1)*wa* x[loc_aisle[a][l].id])

    return pricing_route_o
end

function update_pricing_obj_pure!(ins::instance,pricing_route_o::Model,mu_o::Float64,pi_o::Array{Float64,1})
# This function upates the pricing problem objective with the rteurn routing
    x = pricing_route_o[:x]
    y = pricing_route_o[:y]
    z = pricing_route_o[:z]
    amax = ins.amax
    n_loc = length(ins.locList)
    @objective(pricing_route_o, Min, 2*z + 2*sum(y[a] for a=1:amax) - mu_o - sum(pi_o[l]*x[l] for l=1:n_loc))
end

## Some insertion and delta evaluation function
function delta_insertion_route(ins::instance,route::route,loc::loc)
# Compute the delta insertion cost to insert a stop in loc into the route
    if isempty(route.loc)
        return 2*dist(ins,ioloc,loc)
    end
    # Find the insertion position in the route
    i = 1
    while i <= length(route.loc)
        if route.loc[i] > loc
            break
        end
        i += 1
    end
    # Now we compute the insertion cost
    delta = 0
    if i == 1
        delta += dist(ins,ioloc,loc) + dist(ins,loc,route.loc[1]) - dist(ins,ioloc,route.loc[1])
    elseif i == length(route.loc) + 1
        delta += dist(ins,route.loc[end],loc) + dist(ins,loc,ioloc) - dist(ins,route.loc[end],ioloc)
    else
        delta += dist(ins,route.loc[i-1],loc) + dist(ins,loc,route.loc[i]) - dist(ins,route.loc[i-1],route.loc[i])
    end

    return delta
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
