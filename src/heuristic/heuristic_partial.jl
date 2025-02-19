module heuristic_partial

using ..structures_decomp, ..cost_functions, ..BB
using Random, StatsBase, JuMP

##############################################
# This module contains the code for a greedy heuristic that produce 
# primal solutions from a partial assignment (at a node from the BandB tree)
# it is launched after generating a dual bound, and does not handle infeasible MP
# This procedure does not account for forbidden position, since its aim is to find
# improved primal bound: even if the solution is not valid for the current node
# this is not a problem since it is a valid solution for the problem
# description of the weights for the heuristic:
# 1 = weight for the best xi
# noise_term = weight for the noise term
####################################################

mutable struct part_sol
    ins::instance
    xi_val::Array{Float64,2}
    remaining_sku::Vector{sku}
    occupation::Dict{loc,Int}
    scores::Vector{Float64}
    current_score_table::Array{Float64,2}
    best_l::Vector{loc}
    assignment::Dict{sku,loc}
    noise_term::Float64

    function part_sol(node::BB.JumpNode, tree::BB.SearchTree)
        sol = new()
        sol.ins = tree.auxiliary_data["instance"]
        sol.xi_val = JuMP.value.(tree.auxiliary_data["xi_vars"])
        # compute the fixed positions from branches
        fixed_branches = collect_fixed_to_one_branches(node,sol.ins)
        # remaining skus
        sol.remaining_sku = filter(s -> s.demand > 0, sol.ins.skuList)
        filter!(s -> !(s in keys(fixed_branches)) , sol.remaining_sku)
        filter!(s -> !(s in keys(sol.ins.fixed_positions)) , sol.remaining_sku)
        # occupation
        sol.occupation = Dict{loc,Int}(l => 0 for l in sol.ins.locList)
        for l in values(fixed_branches)
            sol.occupation[l] += 1
        end
        for l in values(sol.ins.fixed_positions)
            sol.occupation[l] += 1
        end
        @assert maximum(values(sol.occupation)) <= sol.ins.n_sym 
        # the rest
        sol.scores = [0. for s in sol.ins.skuList]
        sol.current_score_table = [0. for l in sol.ins.locList, s in sol.ins.skuList]
        sol.best_l = [sol.ins.locList[1] for s in sol.ins.skuList] 
        # sol assignment
        sol.assignment = Dict{sku,loc}()
        for s in keys(fixed_branches)
            sol.assignment[s] = fixed_branches[s]
        end
        for s in keys(sol.ins.fixed_positions)
            sol.assignment[s] = sol.ins.fixed_positions[s]
        end
        sol.noise_term = tree.auxiliary_data["config"].partial_heuristic_noise
        return sol
    end
end

function main_partial_heuristic(node::BB.JumpNode, tree::BB.SearchTree)
    # This is the main function
    #  initialization
    sol = part_sol(node,tree)
    compute_scores!(sol)
    while !isempty(sol.remaining_sku)
        insert_best_sku_and_update!(sol)
    end
    # post processing of the solution
    # UNIT testint
    # first test the assignment
    key_list = [s for s in keys(sol.assignment)]
    for s in sol.ins.skuList
        @assert s in key_list || s.demand == 0
    end
    loc_count = [0 for l in sol.ins.locList]
    for s in key_list
        loc_count[sol.assignment[s].id] += 1
    end
    @assert maximum(loc_count) <= sol.ins.n_sym
    # end UNIT test
    # we add the skus with a demand of 0 to the solution, so that all sku are placed
    for s in sol.ins.skuList
        if !(s in key_list)
            for l in sol.ins.locList 
                if loc_count[l.id] < sol.ins.n_sym
                    sol.assignment[s] = l
                    loc_count[l.id] += 1
                    break
                end
            end
        end
    end

    # now we create a list of routes and compute the cost of the solution
    route_list = route[]
    cost = 0
    for o in sol.ins.orderList
        push!(route_list,cost_functions.route_cost(sol.ins,[sol.assignment[s] for s in o.sku],o))
        cost += route_list[end].cost
    end
    return cost, route_list, sol.assignment
end

function compute_scores!(sol::part_sol)
    # This function computes the initial score table
    for s in sol.remaining_sku
        for l in sol.ins.locList
            sol.current_score_table[l.id,s.id] = sol.occupation[l] < sol.ins.n_sym ? sol.xi_val[l.id,s.id] + sol.noise_term*rand() : -Inf
        end
    end
    # now we compute the best l and the scores
    for s in sol.remaining_sku
        sol.scores[s.id] = maximum(sol.current_score_table[:,s.id])
        sol.best_l[s.id] = sol.ins.locList[findmax(sol.current_score_table[:,s.id])[2]]
    end
end

function insert_best_sku_and_update!(sol::part_sol)
    # This function insert the best sku, and update the scores and occupation
    best_sku_ind = findmax(sol.scores[[s.id for s in sol.remaining_sku]])[2]
    best_sku_id = sol.remaining_sku[best_sku_ind].id
    insertion_loc = sol.best_l[best_sku_id]
    # insert it in the solution
    sol.assignment[sol.ins.skuList[best_sku_id]] = insertion_loc
    # update occupation
    sol.occupation[insertion_loc] += 1
    #remove the sku from the remaining ones
    deleteat!(sol.remaining_sku, best_sku_ind)
    #filter!(s -> s.id != best_sku_id, sol.remaining_sku)
    # update the scores
    if sol.occupation[insertion_loc] == sol.ins.n_sym
        for s in sol.remaining_sku
            sol.current_score_table[insertion_loc.id,s.id] = - Inf
        end
    end
    for s in sol.remaining_sku
        sol.scores[s.id] = maximum(sol.current_score_table[:,s.id])
        sol.best_l[s.id] = sol.ins.locList[findmax(sol.current_score_table[:,s.id])[2]]
    end
end

function collect_fixed_to_one_branches(node::BB.JumpNode,ins::instance)
    # This function collects the sku fixed to some locs with branching constraints
    # in order to convexify these constrainst and pass them to the subproblems
    fixed_branches = Dict{sku,loc}()
    pnode = node
    if !isnothing(node.branch) && !isnothing(node.branch.fixed_sku_id)
        fixed_branches[ins.skuList[node.branch.fixed_sku_id]] = ins.locList[node.branch.fixed_loc_id]
    end
    while !isnothing(pnode.parent)
        pnode = pnode.parent
        if !isnothing(pnode.branch) && !isnothing(pnode.branch.fixed_sku_id)
            fixed_branches[ins.skuList[pnode.branch.fixed_sku_id]] = ins.locList[pnode.branch.fixed_loc_id]
        end
    end
    return fixed_branches
end
end