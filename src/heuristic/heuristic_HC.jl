#using Random

mutable struct HC_sol <:BB.AbstractSolution
    instance::instance
    plan::Dict{sku,loc}
    cost::Int
    routes::Array{route,1}

    function HC_sol(
        instance::instance,
        plan::Dict{sku,loc},
        cost = cost_functions.sol_cost(instance,plan),
        routes = cost_functions.create_solution(instance,plan).routes)
        return new(instance,plan,cost,routes)
    end
end


struct HC_config
    iter_without_imp::Int
    shuffling_prob::Float64 # proba to apply shuffling on an iter
    max_shuffle::Float64 # number max of sku shuffled
    noise_coef::Float64 # coefficient for the noise term
end


function compute_exchange_cost(sol::HC_sol,a::sku,b::sku)
    # Compute the cost of exchanging 2 sku, with delta evaluation
    # Maybe for improvement -> better delta evaluation, with accounting for the routing policy
    # if it is the same loc, then return 0
    if sol.plan[a].id == sol.plan[b].id
        return 0
    end
    # delta
    new_cost = 0
    #=sol_test = deepcopy(sol)
    for o in sol.instance.orderList
        @assert(length(sol.routes[o.id].loc) == length(sol_test.routes[o.id].loc))
        for l=1:length(sol.routes[o.id].loc)
            @assert(sol.routes[o.id].loc[l] == sol_test.routes[o.id].loc[l])
        end
    end
    for o=1:length(sol.instance.orderList)
        @assert(sol.routes[o].order.id == o)
    end=#
    for o in sol.instance.orderList
        # if both sku are in the route, or neither of them -> continue
        a_in_o = a in o.sku
        b_in_o = b in o.sku

        if a_in_o == b_in_o
            new_cost += sol.routes[o.id].cost
            continue
        end
        # get the new loc_group, with a shallow copy
        loc_group = copy(sol.routes[o.id].loc)
        # replace the loc
        to_replace_loc_id = (a_in_o ? sol.plan[a].id : sol.plan[b].id)
        for i=1:length(loc_group)
            # continue until we find one to replace
            if (loc_group[i].id != to_replace_loc_id)
                continue
            end

            # replace the loc
            if a_in_o
                loc_group[i] = sol.plan[b]
            else
                loc_group[i] = sol.plan[a]
            end
            # break the loop not to replace more than one loc
            break
        end

        # compute the delta cost
        new_cost += cost_functions.route_cost(sol.instance,loc_group,o).cost
    end
    return new_cost - sol.cost
end

function compute_exchange_three_cost(sol::HC_sol,a::sku,b::sku,c::sku)
    # Compute the cost of exchanging 3 sku with the permutation
    # a -> loc of b
    # b -> loc of c
    # c -> loc of a
    
    # if it is the same loc, then return 0
    if sol.plan[a].id == sol.plan[b].id
        return 0
    end
    # delta
    new_cost = 0
    
    for o in sol.instance.orderList
        # if both sku are in the route, or neither of them -> continue
        a_in_o = a in o.sku
        b_in_o = b in o.sku
        c_in_o = c in o.sku

        if a_in_o == b_in_o && a_in_o == c_in_o
            new_cost += sol.routes[o.id].cost
            continue
        end
        # get the new loc_group, with a shallow copy
        loc_group = copy(sol.routes[o.id].loc)
        # replace the loc
        loc_to_replace = []
        if a_in_o
            for i in eachindex(loc_group)
                if loc_group[i] == sol.plan[a]
                    loc_group[i] = sol.plan[b]
                    break
                end
            end
        end
        if b_in_o
            for i in eachindex(loc_group)
                if loc_group[i] == sol.plan[b]
                    loc_group[i] = sol.plan[c]
                    break
                end
            end
        end
        if c_in_o
            for i in eachindex(loc_group)
                if loc_group[i] == sol.plan[c]
                    loc_group[i] = sol.plan[a]
                    break
                end
            end
        end
        

        # compute the delta cost
        new_cost += cost_functions.route_cost(sol.instance,loc_group,o).cost
    end
    #prinln("newcost = $new_cost, old cost = $(sol.cost)")
    return new_cost - sol.cost
end




function make_exchange!(sol::HC_sol,a::sku,b::sku)
# This function proceed at the exchange between the two sku
    # compute cost of the exchange
    delta = compute_exchange_cost(sol,a,b)
    # change the routes
    for o in sol.instance.orderList
        # if both sku are in the route, or neither of them -> continue
        a_in_o = a in o.sku
        b_in_o = b in o.sku
        if a_in_o == b_in_o
            continue
        end
        # get the new loc_group, with a shallow copy
        loc_group = copy(sol.routes[o.id].loc)
        # replace the loc
        to_replace_loc_id = (a_in_o ? sol.plan[a].id : sol.plan[b].id)
        for i=1:length(loc_group)
            # continue until we find one to replace
            if (loc_group[i].id != to_replace_loc_id)
                continue
            end
            # replace the loc
            if a_in_o
                loc_group[i] = sol.plan[b]
            else
                loc_group[i] = sol.plan[a]
            end
            # break the loop not to replace more than one loc
            break
        end
        # replace the route
        sol.routes[o.id] = cost_functions.route_cost(sol.instance,loc_group,o)
    end
    # change plan
    old_a_loc = sol.plan[a]
    sol.plan[a] = sol.plan[b]
    sol.plan[b] = old_a_loc
    # change cost
    sol.cost += delta
    #@assert(sol.cost == cost_functions.sol_cost(sol.instance,sol.plan))
end



function make_three_exchange!(sol::HC_sol,a::sku,b::sku,c::sku)
    # This function proceed at the exchange between the three skus
    # compute cost of the exchange
    delta = compute_exchange_three_cost(sol,a,b,c)
    # change the routes
    for o in sol.instance.orderList
        # if both sku are in the route, or neither of them -> continue
        a_in_o = a in o.sku
        b_in_o = b in o.sku
        c_in_o = c in o.sku
        if a_in_o == b_in_o && a_in_o == c_in_o
            continue
        end
    
        # get the new loc_group, with a shallow copy
        loc_group = copy(sol.routes[o.id].loc)
        # replace the loc
        if a_in_o
            for i in eachindex(loc_group)
                if loc_group[i] == sol.plan[a]
                    loc_group[i] = sol.plan[b]
                    break
                end
            end
        end
        if b_in_o
            for i in eachindex(loc_group)
                if loc_group[i] == sol.plan[b]
                    loc_group[i] = sol.plan[c]
                    break
                end
            end
        end
        if c_in_o
            for i in eachindex(loc_group)
                if loc_group[i] == sol.plan[c]
                    loc_group[i] = sol.plan[a]
                    break
                end
            end
        end
        # replace the route
        sol.routes[o.id] = cost_functions.route_cost(sol.instance,loc_group,o)
    end
    # change plan
    old_a_loc = sol.plan[a]
    sol.plan[a] = sol.plan[b]
    sol.plan[b] = sol.plan[c]
    sol.plan[c] = old_a_loc
    # change cost
    sol.cost += delta
    @show sol.cost
    @show delta
    @show cost_functions.sol_cost(sol.instance,sol.plan)
    @assert(sol.cost == cost_functions.sol_cost(sol.instance,sol.plan))
end


function shuffle_solution!(sol::HC_sol,max_shuffle::Float64)
# This function shuffles the solution by doing a permutation of at most max_shuffle % of the solution
# The number of sku permuted is at least 3, at most according to the percentage max
    ins = sol.instance
    # only in the available ones
    available_sku = setdiff(ins.skuList, keys(ins.fixed_positions))
    # number of sku shuffled
    number_shuffled = rand(3:max(3,round(Int,max_shuffle*length(available_sku))))
    # select these sku
    sku_shuffled = sample(available_sku, number_shuffled, replace=false)
    # make a copy of the plan
    new_plan = copy(sol.plan)
    # get a permutation
    perm = randperm(number_shuffled)
    # perform the permutation
    for i=1:number_shuffled
        new_plan[sku_shuffled[i]] = sol.plan[sku_shuffled[perm[i]]]
    end
    # replace the solution
    sol = HC_sol(sol.instance,new_plan)
end


function main_HC_algo(instance::structures_decomp.instance,config::HC_config)
# main HC procedure
    # register fixed loc
    plan = Dict{sku,loc}()
    loc_available = Dict{loc,Int}(l => instance.n_sym for l in instance.locList)
    for s in keys(instance.fixed_positions)
        plan[s] = instance.fixed_positions[s]
        loc_available[instance.fixed_positions[s]] -= 1
    end
    # Create random solution
    loc_group = loc[]
    for l in instance.locList
        for i=1:loc_available[l]
            push!(loc_group,l)
        end
    end
    # shuffle the loc and place the remaining sku
    # CHANGE RANDOM
    #sort!(loc_group, by = l -> dist(instance,structures_decomp.ioloc,l))
    available_sku = setdiff(instance.skuList,keys(instance.fixed_positions))
    for i in eachindex(available_sku)
        plan[available_sku[i]] = loc_group[i]
    end

    # Praparation for the iterations
    best_solution = HC_sol(instance,plan)
    current_solution = deepcopy(best_solution)
    iter_without_imp = 0
    # max time of the heuristic
    max_time = time() + 60
    # Iterations
    # CHANGE RANDOM
    while (iter_without_imp < config.iter_without_imp && time() < max_time)
        break
        # Perform the iteration
        HC_iter!(current_solution,config)
        iter_without_imp += 1
        # test if new best solution
        if current_solution.cost < best_solution.cost
            best_solution = deepcopy(current_solution)
            iter_without_imp = 0
        end
    end

    #checkeur
    ins = instance
    # all sku placed
    @assert length(ins.skuList) == length(keys(best_solution.plan))
    @assert length(intersect(ins.skuList,keys(best_solution.plan))) == length(ins.skuList)
    # position filling is ok
    @assert maximum(values(count_occupation(ins,best_solution.plan))) <= ins.n_sym
    # the fixed positions are respected
    for s in keys(ins.fixed_positions)
        @assert best_solution.plan[s].id == ins.fixed_positions[s].id
    end
    println("best solution cost = $(best_solution.cost), sol cost = $(cost_functions.sol_cost(ins,best_solution.plan))")
    # check cost
    @assert best_solution.cost == cost_functions.sol_cost(ins,best_solution.plan)

    return best_solution
end

main_HC_algo(path::String,config::HC_config) = main_HC_algo(structures_decomp.read_instance_file(path),config)
main_HC_algo(ins::instance,conf::Dict{String,Any}) = main_HC_algo(ins,HC_config(conf["iter"],conf["shuffle_prob"],conf["shuffle_max"],conf["noise"]))

function HC_iter!(current_solution::HC_sol,config::HC_config)
# This function performs an iteration of the Hill climing algorithm
    # Check if we shuffle
    #println("iter !")
    if rand() < config.shuffling_prob
        shuffle_solution!(current_solution,config.max_shuffle)
        return true
    end
    # Pick randomly 2 sku to perform the move
    ins = current_solution.instance
    # only in the available ones
    available_sku = setdiff(ins.skuList, keys(ins.fixed_positions))
    best_delta = 0
    best_a = nothing
    best_b = nothing
    for a_ind = 1:length(available_sku)-1
        for b_ind = a_ind+1:length(available_sku)
            # compute the delta of the move
            delta = compute_exchange_cost(current_solution, available_sku[a_ind], available_sku[b_ind])
            # add noise term
            delta += config.noise_coef*(2*rand() - 1)*(ins.wa*ins.amax + ins.wc*ins.cmax)
            if delta < best_delta
                best_delta = delta
                best_a = available_sku[a_ind]
                best_b = available_sku[b_ind]
            end
        end
    end
    # perform the move if delta < 0
    if best_delta < 0
        make_exchange!(current_solution,best_a,best_b)
        return true
    else
        # in this case we try the exchange of 3 skus
        #=best_delta = 0
        best_a = nothing
        best_b = nothing
        best_c = nothing
        for a in available_sku, b in available_sku, c in available_sku
            if a == b || a == c || b == c
                continue
            end
            # compute the delta of the move
            delta = compute_exchange_three_cost(current_solution, a, b, c)
            # add noise term
            delta += config.noise_coef*(2*rand() - 1)*(ins.wa*ins.amax + ins.wc*ins.cmax)
            if delta < best_delta
                best_delta = delta
                best_a = a
                best_b = b
                best_c = c
                break
            end
        end
        if best_delta < 0
            make_three_exchange!(current_solution,best_a,best_b,best_c)
            return true
        end=#
    end

    return false
end
