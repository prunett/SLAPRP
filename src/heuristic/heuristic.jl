
module heuristic
#This module provides heuristic tools to solve the problem, mainly to find initial solutions
# Might be interesting to code some kind of regret heuristic

using ..structures_decomp, ..cost_functions, ..BB
using StatsBase, Random

Random.seed!(1234)

include("heuristic_HC.jl")


mutable struct partial_sol
    instance::instance
    fixed_positions::Dict{sku,loc}
    forbidden_positions::Dict{sku,Array{loc,1}}
    cost::Int
    remaining_sku::Array{sku,1}
    occupation::Dict{loc,Int}
    routes::Array{route,1}

    function partial_sol(
        instance,
        fixed_positions = Dict{sku,loc}(),
        forbidden_positions = Dict{sku,Array{loc,1}}(),
        cost = sol_cost(instance,fixed_positions,true),
        remaining_sku = setdiff(instance.skuList, keys(fixed_positions)),
        occupation = count_occupation(instance,fixed_positions),
        routes = cost_functions.create_solution(instance,fixed_positions,true).routes)
        return new(instance, fixed_positions, forbidden_positions, cost, remaining_sku, occupation, routes)
    end
end

function count_occupation(ins::instance,fixed_positions::Dict{sku,loc})
# counts the number of occupation in each loc
    occupation = Dict{loc,Int}(l => 0 for l in ins.locList)
    for s in keys(fixed_positions)
        occupation[fixed_positions[s]] += 1
    end
    return occupation
end

function greedy_heuristic(ins::instance)
# This is a basic greedy heuristic that insert succesively the sku at the best place in a loc plan,
# it considers first the sku with highest demand
    #sort sku list
    skuList = deepcopy(ins.skuList)
    sort!(skuList, by= s -> s.demand, rev=true)

    #copy locList & create empty plan
    locList = copy(ins.locList)
    plan = Dict{sku,loc}()
    loc_number = Dict{loc,Int}(l => 0 for l in locList)

    for s in skuList
        min_cost = Inf
        ind_loc = locList[1]
        for l in locList
            # test if a space is available at loc l
            if loc_number[l] == ins.n_sym
                continue
            end
            # test the insertion
            new_plan = copy(plan)
            new_plan[s] = l
            if sol_cost(ins, new_plan) < min_cost
                min_cost = sol_cost(ins,new_plan)
                ind_loc = l
            end
        end
        # update the plan with the best insertion
        plan[s] = ind_loc
        loc_number[ind_loc] += 1
    end

    return create_solution(ins,plan)
end




function nothing_heuristic(ins::instance)
# This function assigns sku 1 to loc 1, sku 2 to loc 2, etc...
    stol = Dict{sku,loc}(s => ins.locList[s.id] for s in ins.skuList)
    return create_solution(ins,stol)
end

############## K-regret heuristic #############

function main_regret_heuristic(ins::instance, k::Int, fixed_positions::Dict{sku,loc}, forbidden_positions::Dict{sku,Array{loc,1}})
# This is the main procedure for k-regret heuristic
    # Generate inputs
    sol = partial_sol(ins,fixed_positions, forbidden_positions)
    # Compute delta
    delta = create_delta(sol)
    # regret_inside_loop
    while !isempty(sol.remaining_sku)
        regret_insertion!(sol,delta,k)
    end

    # solution checker
    # all sku placed
    @assert length(ins.skuList) == length(keys(sol.fixed_positions))
    @assert length(intersect(ins.skuList,keys(sol.fixed_positions))) == length(ins.skuList)
    # position filling is ok
    @assert maximum(values(count_occupation(ins,sol.fixed_positions))) <= ins.n_sym
    # cost is ok
    @assert sol.cost == cost_functions.sol_cost(ins,sol.fixed_positions)
    #println("aa")
    #println("solution cost = $(sol.cost)")

    return sol

end

main_regret_heuristic(ins::instance,k::Int) = main_regret_heuristic(ins,k,Dict{sku,loc}(),Dict{sku,Array{loc,1}}())
main_regret_heuristic(path::String,k::Int) = main_regret_heuristic(structures_decomp.read_instance_file(path),k)

function regret_insertion!(sol::partial_sol, delta::Array{Int,2}, k::Int)
# This is the inside loop that perform one insertion
    # Compute xsj -> location index with the j-th lowest insertion cost for sku s
    x = fill(typemax(Int),(sol.instance.n, k))
    for s in sol.remaining_sku
        p = sortperm(delta[s.id,:])
        j = 1 # ind to read delta and p
        ind = 1 # ind to read x
        while j <= length(sol.instance.locList) && ind <= k
            u = sol.occupation[sol.instance.locList[p[j]]]
            while u < sol.instance.n_sym && ind <= k
                x[s.id,ind] = p[j]
                u += 1
                ind += 1
            end
            j += 1
        end
        #@assert issorted(delta[s.id,x[s.id,:]])
    end
    # Compute csk
    c = fill(-1,sol.instance.n)
    for s in sol.remaining_sku
        reduced_k = length(sol.instance.locList)*sol.instance.n_sym - length(keys(sol.fixed_positions))
        c[s.id] = sum(delta[s.id, x[s.id,j]] - delta[s.id, x[s.id,1]] for j=1:min(k, reduced_k))
        @assert c[s.id] >= 0
    end
    # Find the sku to isnert in the solution
    ind = 0
    c_max = -1
    cost = Inf
    for s in sol.remaining_sku
        if c[s.id] > c_max
            ind = s.id
            c_max = c[s.id]
            cost = delta[s.id,x[s.id,1]]
        elseif (c[s.id] == c_max) && delta[s.id, x[s.id,1]] < cost
            ind = s.id
            cost = delta[s.id, x[s.id,1]]
        end
    end
    # insert the sku in the partial solution
    insert_sku!(sol, sol.instance.skuList[ind], sol.instance.locList[x[ind,1]])
    # Update delta
    update_delta!(sol, delta, sol.instance.skuList[ind])
end

function create_delta(sol::partial_sol)
# This is the function that compute delta, another function is used for update. delta[s,l]
    delta = Array{Int,2}(undef,(sol.instance.n,length(sol.instance.locList)))
    for s in sol.remaining_sku
        # skip if the sku is already placed
        for l in sol.instance.locList
            # skip if the sku is already placed
            if (s in keys(sol.forbidden_positions)) && (l in forbidden_positions[s])
                delta[s.id,l.id] = Inf
                continue
            end
            delta[s.id,l.id] = cost_functions.delta_insertion_main(sol,s,l)
        end
    end

    return delta
end

function update_delta!(sol::partial_sol,delta::Array{Int,2},sku_added::sku)
# This function updates delta for the relevant sku
    # Construct the list of sku affected
    affected_sku = sku[]
    for s in sol.remaining_sku
        for o in sol.instance.sku_orderList[sku_added]
            if s in o.sku
                push!(affected_sku,s)
                break
            end
        end
    end
    # Modify the delta values accordingly
    for s in affected_sku
        for l in sol.instance.locList
            # if the value is already infinite, skip it
            if delta[s.id,l.id] == Inf
                continue
            end
            delta[s.id,l.id] = cost_functions.delta_insertion_main(sol,s,l)
        end
    end
end



function insert_sku!(sol::partial_sol,sku::sku,loc::loc)
# insert the sku in the partial solution at location loc, and update it
    #Update fixed position
    sol.fixed_positions[sku] = loc
    # update occupation
    sol.occupation[loc] += 1
    # Update remaining sku
    filter!(s -> s.id != sku.id, sol.remaining_sku)
    # Create new solution
    for o in sol.instance.sku_orderList[sku]
        route = cost_functions.route_cost(sol.instance, sol.fixed_positions, sol.instance.orderList[o.id], true)
        sol.cost += route.cost - sol.routes[o.id].cost
        sol.routes[o.id] = route
    end
end

# some new methods
cost_functions.delta_insertion_main(sol::partial_sol,s::sku,l::loc) = cost_functions.delta_insertion_main(sol.instance, sol.routes, sol.occupation, s, l)


####################################################
# Greedy related heuristic
###################################################

function greedy_related_heuristic(ins::instance, weights::Dict{String,Float64}, fixed_positions::Dict{sku,loc}, forbidden_positions::Dict{sku,Array{loc,1}})
# This is the main procedure for the greedy related insertion heuristic
    # Generate input
    sol = partial_sol(ins,fixed_positions, forbidden_positions)
    locList = filter(l -> sol.occupation[l] < ins.n_sym, ins.locList)
    # Order the available loc by distance to ioloc
    sort!(locList, by = l -> cost_functions.dist(ins,ioloc,l))
    # Main loop
    for l in locList
        while sol.occupation[l] < ins.n_sym && !isempty(sol.remaining_sku)
            # Find best sku
            s = find_best_sku_related(sol,weights,l)
            # perform the insertion
            insert_sku!(sol,s,l)
        end
    end

    # solution checker
    # all sku placed
    @assert length(ins.skuList) == length(keys(sol.fixed_positions))
    @assert length(intersect(ins.skuList,keys(sol.fixed_positions))) == length(ins.skuList)
    # position filling is ok
    @assert maximum(values(count_occupation(ins,sol.fixed_positions))) <= ins.n_sym
    # cost is ok
    @assert sol.cost == cost_functions.sol_cost(ins,sol.fixed_positions)
    #println("aa")
    #println("solution cost = $(sol.cost)")

    return sol
end

greedy_related_heuristic(ins::instance,weights::Dict{String,Float64}) = greedy_related_heuristic(ins,weights,Dict{sku,loc}(),Dict{sku,Array{loc,1}}())
greedy_related_heuristic(path::String,weights::Dict{String,Float64}) = greedy_related_heuristic(structures_decomp.read_instance_file(path),weights)


function find_best_sku_related(sol::partial_sol,weights::Dict{String,Float64},loc::loc)
# This function find the best sku to insert in position l, depending on several criteria with weighted sum
    # Initialize
    best_sku = sol.remaining_sku[1]
    max_score = 0.
    for s in sol.remaining_sku
        # demand term
        if s.demand == 0
            demand_term = 0.
        else
            demand_term = s.demand / maximum(sku.demand for sku in sol.remaining_sku)
        end
        # mean order length term
        if s.demand == 0
            order_length_term = 0.
        else
            order_length_term = (maximum(length(o.sku) for o in sol.instance.orderList) - mean(length(o.sku) for o in sol.instance.sku_orderList[s])) / maximum(length(o.sku) for o in sol.instance.orderList)
        end
        # distance to other sku term
        # first we compute the occurences of connected sku
        sku_placed_weights = Dict{sku,Int}(sku => 0 for sku in keys(sol.fixed_positions))
        for o in sol.instance.sku_orderList[s]
            for sbis in o.sku
                if sbis in keys(sku_placed_weights)
                    sku_placed_weights[sbis] += 1
                end
            end
        end
        # mean distance to other sku placed together already placed, divided by largest distance. 0 if first sku to be placed
        if isempty(sku_placed_weights) || maximum(values(sku_placed_weights)) == 0
            distance_term = 0
        else
            distance_term = 0
            for s in keys(sol.fixed_positions)
                distance_term += cost_functions.dist(sol.instance,loc,sol.fixed_positions[s]) * sku_placed_weights[s]
            end
            # make the mean
            distance_term = distance_term / sum(values(sku_placed_weights))
            # divide by the max distance
            distance_term = distance_term / (sol.instance.wa*sol.instance.amax + sol.instance.wc*sol.instance.cmax)
            # transform into the opposite so that the objective is to max it
            distance_term = 1 - distance_term
        end

        # test if it is better than previous best
        score = demand_term*weights["demand"] + order_length_term*weights["order_length"] + distance_term*weights["distance"]
        #println("score = $score, demand = $demand_term, length = $order_length_term, distance = $distance_term")
        if score > max_score
            max_score = score
            best_sku = s
        end
    end

    return best_sku
end



end #end of module
