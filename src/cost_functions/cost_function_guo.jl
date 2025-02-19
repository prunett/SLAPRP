module cost_functions_guo

using JuMP, CPLEX

using ..structures_decomp

#export dist, route_cost, sol_cost, create_solution

# Defines a new isless method, to use the sorting functions on locations
#=function Base.isless(a::loc,b::loc)
    if a.id == 0 || b.id == 0
        return true
    elseif isodd(a.aisle) && isodd(b.aisle)
        return a.id < b.id
    elseif iseven(a.aisle) && iseven(b.aisle)
        return a.aisle > b.aisle || (a.aisle == b.aisle && a.col < b.col)
    elseif isodd(a.aisle) && iseven(b.aisle)
        return true
    else
        return false
    end
end=#
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

dist_matrix = 
[0	4.5	5.5	6.5	7.5	8.5	4.5	5.5	6.5	7.5	8.5	6.5	7.5	8.5	9.5	10.5	6.5	7.5	8.5	9.5	10.5;
4.5	0	1	2	3	4	7	8	9	10	11	9	10	11	12	13	9	10	11	12	13;
5.5	1	0	1	2	3	8	9	10	11	12	10	11	12	13	14	10	11	12	13	14;
6.5	2	1	0	1	2	9	10	11	12	13	11	12	13	14	15	11	12	13	14	15;
7.5	3	2	1	0	1	10	11	12	13	14	12	13	14	15	16	12	13	14	15	16;
8.5	4	3	2	1	0	11	12	13	14	15	13	14	15	16	17	13	14	15	16	17;
4.5	7	8	9	10	11	0	1	2	3	4	9	10	11	12	13	9	10	11	12	13;
5.5	8	9	10	11	12	1	0	1	2	3	10	11	12	13	14	10	11	12	13	14;
6.5	9	10	11	12	13	2	1	0	1	2	11	12	13	14	15	11	12	13	14	15;
7.5	10	11	12	13	14	3	2	1	0	1	12	13	14	15	16	12	13	14	15	16;
8.5	11	12	13	14	15	4	3	2	1	0	13	14	15	16	17	13	14	15	16	17;
6.5	9	10	11	12	13	9	10	11	12	13	0	1	2	3	4	7	8	9	10	11;
7.5	10	11	12	13	14	10	11	12	13	14	1	0	1	2	3	8	9	10	11	12;
8.5	11	12	13	14	15	11	12	13	14	15	2	1	0	1	2	9	10	11	12	13;
9.5	12	13	14	15	16	12	13	14	15	16	3	2	1	0	1	10	11	12	13	14;
10.5	13	14	15	16	17	13	14	15	16	17	4	3	2	1	0	11	12	13	14	15;
6.5	9	10	11	12	13	9	10	11	12	13	7	8	9	10	11	0	1	2	3	4;
7.5	10	11	12	13	14	10	11	12	13	14	8	9	10	11	12	1	0	1	2	3;
8.5	11	12	13	14	15	11	12	13	14	15	9	10	11	12	13	2	1	0	1	2;
9.5	12	13	14	15	16	12	13	14	15	16	10	11	12	13	14	3	2	1	0	1;
10.5	13	14	15	16	17	13	14	15	16	17	11	12	13	14	15	4	3	2	1	0
]








function dist(a::loc,b::loc,ins::instance)
# Computes the walking distance between loc a and b, with return policy

    return dist_matrix[a.id+1,b.id+1]
    
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
