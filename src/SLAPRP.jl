module SLAPRP

#adding BB related modules
include("BandB/mainBB.jl")

# adding CG related modules
include("structures/structures_decomp.jl")

include("cost_functions/cost_functions_exact.jl")
include("cost_functions/cost_functions_return.jl")
include("cost_functions/cost_functions_midpoint.jl")
include("cost_functions/cost_functions_sshape.jl")
include("cost_functions/cost_functions_largest.jl")
include("cost_functions/cost_function_guo.jl")
include("cost_functions/cost_functions.jl")

include("cuts/cuts_main.jl")
include("heuristic/heuristic.jl")
include("heuristic/heuristic_partial.jl")
include("Column_generation/CG_utils.jl")
include("pricing/pricing_exact.jl")
include("pricing/pricing_return.jl")
include("pricing/pricing_general.jl")
include("pricing/pricing.jl")
include("Column_generation/CG_main.jl")




using JuMP, MathOptInterface

const MOI = MathOptInterface
const EPS = 0.001

struct bcp_config
    time_limit::Int
    cpu_limit::Int
    EPS::Float64
    instance_type::String
    max_iter_cut::Int
    EPS_cut_separation_initial::Float64
    EPS_cut_violated_initial::Float64
    EPS_delta_separation::Float64
    EPS_delta_violated::Float64
    EPS_cut_obj_improvement::Float64
    EPS_cut_vio_SL1::Float64
    max_iter_SL1::Int
    max_round_OG::Int
    min_nodes_with_cuts::Int
    max_depth_cut_generation::Int
    search_strategy::String
    branching_strategy::String
    branching_weight_gub_demand::Float64
    branching_weight_gub_xinumber::Float64
    branching_weight_gub_balance::Float64
    branching_weight_gub_closeness::Float64
    branching_weight_gl_distance::Float64
    branching_weight_gl_demand::Float64
    branching_weight_gl_xinumber::Float64
    branching_weight_gl_closeness::Float64
    nb_candidates::Int
    bound_method::String
    partial_heuristic_noise::Float64
    fake_route::Bool
    routing_policy::String
    bi_directional::Bool
    activate_all_link_cuts::Bool # if they are activated from the start
    initiate_aisle_cuts_in_pool::Bool # if they are generated and placed in the cut pool
    simplex_algo::Int
    max_cut_per_order::Int
end

##############################
## Adapt B and B functions to the SLAP
##################################

# Branching
include("BandB/branching.jl")


    

function BB.dual_bound_node!(node::BB.JumpNode,tree::BB.SearchTree,max_time)
# This function computes the dual bound using column generation
    # solve the CG
    valid_dual, obj_master, print_message = SLAP_CG_main.CG_iteration!(tree,node,EPS,max_time)
    # Query solution status
    node.solution_status = JuMP.termination_status(node.model)
    # If it is infeasible then launch the primal bound procedure and restart the CG
    if node.solution_status == MOI.INFEASIBLE
        #JuMP.write_to_file(node.model, "master_problem.lp")
        try
            BB.print_solution(tree,node)
        catch
            println("INFEASIBLE MP with no solution in the model")
        end
        n = node

        # print the branching decisions
        branches = BB.collect_all_branches(node)
        for b in branches 
            println(b.ub)
            println(b.lb)
        end
        xi_vars = node.model[:xi]
        #println("lower bound for xi[11,5] = $(JuMP.lower_bound(xi_vars[11,5]))")
        
        #throw(MethodError("Infeasible MP !!"))
        return "Pruned by infeasability"
        #BB.primal_bound_node!(node,tree)
        #valid_dual, obj_master, print_message = SLAP_CG_main.CG_iteration!(tree,node,EPS,max_time)
    end
    @assert node.solution_status == MOI.OPTIMAL
    # start comment
    #println("Variables fixed to 0 => $fixed_0, fixed to 1 => $fixed_1")
    # update the node infos
    if valid_dual
        node.dual_bound = obj_master
    end

    # test if the value for the fake routes is < EPS
    if tree.auxiliary_data["config"].fake_route
        for o in tree.auxiliary_data["instance"].orderList
            if !(JuMP.value(tree.auxiliary_data["rho_vars"][o][1]) < tree.auxiliary_data["config"].EPS)
                println("o = $(o.id) value = $(JuMP.value(tree.auxiliary_data["rho_vars"][o][1]))")
                BB.print_solution(tree,node)
                #JuMP.print(node.model)
                nodet = node
                while !isnothing(nodet.parent)
                    @show nodet.branch
                    nodet = nodet.parent
                end
                ins = tree.auxiliary_data["instance"]
                fixed_branches = SLAP_CG_main.collect_fixed_to_one_branches(node,ins)
                #println("fixed_branches function = $fixed_branches")
                for r in tree.auxiliary_data["a_coef"][o]
                    println("route r locs = $([l.id for l in r.loc])")
                end
                # call for the pricing with return
                ins = tree.auxiliary_data["instance"]
                master_problem = node.model
                mu = dual.(master_problem[:rho_sum]) + dual.(master_problem[:rho_sum2])
                pi = - dual.(master_problem[:link_constraint])
                pi_o = pi[:,o.id]
                # Compute the sigma dual costs
                sigma_o = [0. for l in ins.locList]
                # first with the Linkcuts
                cut_pool = tree.auxiliary_data["cut_pool"]
                for l in ins.locList, s in o.sku
                    cut = cut_pool["LinkCut"][l.id,o.id,s.id]
                    if cut.active
                        sigma_o[cut.loc.id] += cut.dual_value
                    end
                end
                
                println("mu_o = $(mu[o.id]), 
                pi_o = $(pi_o)
                sigma_o = $(sigma_o)")
                println("node id  = $(node.id)")
            end
            
            @assert JuMP.value(tree.auxiliary_data["rho_vars"][o][1]) < tree.auxiliary_data["config"].EPS
        end
    end


    return print_message
end

function transform_branches_into_dict(tree::BB.SearchTree,branch_objects::Array{BB.VariableBranch,1})
    # transform the input for the primal heuristic
    fixed_positions = Dict{structures_decomp.sku,structures_decomp.loc}()
    forbidden_positions = Dict{structures_decomp.sku,Array{structures_decomp.loc,1}}()
    for branch in branch_objects
        # begin by lower bounds => equal to 1 => fixed_positions
        for var_lb in keys(branch.lb)
            # test the value
            @assert(branch.lb[var_lb] == 1)
            # retrieve the corresponding sku and loc
            lid , sid = CG_utils.read_xi_name(var_lb)
            # test if already existing
            if tree.auxiliary_data["instance"].skuList[sid] in keys(fixed_positions)
                println(branch_objects)
            end
            @assert(!(tree.auxiliary_data["instance"].skuList[sid] in keys(fixed_positions)))
            # register the info
            fixed_positions[tree.auxiliary_data["instance"].skuList[sid]] = tree.auxiliary_data["instance"].locList[lid]
        end
        # upper bound => equal to 0 => forbidden positions
        for var_ub in keys(branch.ub)
            # test the value
            @assert(branch.ub[var_ub] == 0)
            # retrieve the corresponding sku and loc
            lid,sid = CG_utils.read_xi_name(var_ub)
            # register the info
            if tree.auxiliary_data["instance"].skuList[sid] in keys(forbidden_positions)
                push!(forbidden_positions[tree.auxiliary_data["instance"].skuList[sid]], tree.auxiliary_data["instance"].locList[lid])
            else
                forbidden_positions[tree.auxiliary_data["instance"].skuList[sid]] = structures_decomp.loc[tree.auxiliary_data["instance"].locList[lid]]
            end
        end
    end

    return fixed_positions, forbidden_positions
end

# new methods
function transform_branches_into_dict(tree::BB.SearchTree,branch_objects::Array{BB.SumBranch,1})
    # transform the input for the primal heuristic
    fixed_positions = Dict{structures_decomp.sku,structures_decomp.loc}()
    forbidden_positions = Dict{structures_decomp.sku,Array{structures_decomp.loc,1}}()
    for branch in branch_objects
        # Begin with fixed positions
        if branch.ub == 0
            for var_ub in branch.sum_variables
                # retrieve the corresponding sku and loc
                lid , sid = CG_utils.read_xi_name(var_ub)
                # register the info
                if !(tree.auxiliary_data["instance"].skuList[sid] in keys(forbidden_positions))
                    forbidden_positions[tree.auxiliary_data["instance"].skuList[sid]] = tree.auxiliary_data["instance"].locList[lid]
                end
            end
        else
            #throw(DomainError("Think again if we cannot add some informations"))
            # Maybe add another attribute about its class to the SumBranch structure
        end
    end

    return fixed_positions, forbidden_positions
end

function BB.primal_bound_node!(node::BB.JumpNode,tree::BB.SearchTree)
# This function computes the primal bound using the partial greedy heuristic
    if node.solution_status == MOI.OPTIMAL
        cost, route_list, sol_assignment = heuristic_partial.main_partial_heuristic(node,tree)
        # Update primal bound 
        node.primal_bound = cost
        if cost < tree.current_primal_bound
            tree.best_primal_solution = sol_assignment
        end
    end
end

function BB.prune_bound(tree::BB.SearchTree,node::BB.JumpNode)
# This function checks if we should prune by bound
    # prune by bounds -> accounting for the integer part of the obj
    ins = tree.auxiliary_data["instance"]
    if node.dual_bound >= tree.current_primal_bound - 2*gcd(ins.wa,ins.wc) + EPS
        return true
    end

    return false
end

function BB.prune_node(tree::BB.SearchTree,node::BB.JumpNode)
# This function check if the node should be pruned, or branched on
    # prune by infeasability
    if node.solution_status in [MOI.INFEASIBLE]
        return true
    end

    @assert node.solution_status == MOI.OPTIMAL

    # prune by bounds -> accounting for the integer part of the obj
    if BB.prune_bound(tree,node)
        return true
    end

    # test to see if the solution is integer
    ins = tree.auxiliary_data["instance"]
    xi_val = value.(node.model[:xi])
    integer = true
    for l in ins.locList, s in ins.skuList
        # we don't care about integrality of sku that have no demand
        if s.demand == 0
            continue
        end
        # This part was for the symmetry
        #=
        sum_value = sum(xi_val[l.id,s.id] for s in sym.skus)
        sum_value -= floor(sum_value)
        if (sum_value > EPS) && (1 - sum_value > EPS)
            integer = false
            break
        end=#
        if (xi_val[l.id,s.id] > EPS) && (1 - xi_val[l.id,s.id] > EPS)
            integer = false
            break
        end
    end



    ########## If the solution is integer, prune by integrity
    if integer
        # check if the rho variables are integer
        for o in keys(tree.auxiliary_data["rho_vars"])
            for val in value.(tree.auxiliary_data["rho_vars"][o])
                if !(val < EPS || val > 1 - EPS)
                    println("val = $val")
                    BB.print_solution(tree,node)
                end
                @assert(val < EPS || val > 1 - EPS)
            end
        end
        # check if it improves the best primal bound
        if tree.current_primal_bound > node.dual_bound + EPS
            # update primal bound node
            node.primal_bound = round(node.dual_bound)
            # update best primal bound
            tree.current_primal_bound = node.primal_bound
            # update best primal solution
            tree.best_primal_solution = collect_sol_from_xi_val(ins,xi_val,EPS)
        end
        # tests for dual bound
        if node.dual_bound < tree.auxiliary_data["min_bound_close_nodes"]
            tree.auxiliary_data["min_bound_close_nodes"] = node.dual_bound
        end

        return true
    end

    return false
end

function collect_sol_from_xi_val(ins::structures_decomp.instance,xi_val::Array{Float64,2},EPS::Float64)
    # collect stol and remaining l
    stol = Dict{structures_decomp.sku,structures_decomp.loc}()
    remaining_loc = structures_decomp.loc[]
    for l in ins.locList
        i = 0
        while i < ins.n_sym
            push!(remaining_loc,l)
            i += 1
        end
    end
    remaining_s = structures_decomp.sku[]
    for s in ins.skuList
        # we don't care about integrality of sku that have no demand
        if s.demand == 0
            push!(remaining_s,s)
            continue
        end
        for l in ins.locList
            if xi_val[l.id,s.id] > 1 - EPS
                stol[s] = l
                deleteat!(remaining_loc, findfirst(l2 -> l2.id == l.id, remaining_loc))
                break
            end
        end
    end

    # add the remaining s to the remaining l
    for ind in eachindex(remaining_s)
        stol[remaining_s[ind]] = remaining_loc[ind]
    end

    # return
    return stol
end

function insert_node_queue!(tree::BB.SearchTree,node::BB.JumpNode)
    if tree.auxiliary_data["config"].search_strategy == "priority_queue"
        insert_node_priority_queue!(tree,node)
    elseif tree.auxiliary_data["config"].search_strategy == "fifo"
        insert_node_fifo!(tree,node)
    elseif tree.auxiliary_data["config"].search_strategy == "lifo"
        insert_node_lifo!(tree,node)
    else
        throw(MethodError("Unknown search strategy"))
    end
end
    

function insert_node_priority_queue!(tree::BB.SearchTree,node::BB.JumpNode)
# this function add the new node at the good position in the priority queue
    # first find the position to insert the node
    ind = 1
    while ind <= length(tree.nodes) && tree.nodes[ind].dual_bound < node.dual_bound
        ind += 1
    end
    # update infos
    tree.node_counter += 1
    node.id = tree.node_counter
    # insert the node in the good position
    insert!(tree.nodes, ind, node)
end

function insert_node_fifo!(tree::BB.SearchTree,node::BB.JumpNode)
    tree.node_counter += 1
    node.id = tree.node_counter
    push!(tree.nodes, node)
end

function insert_node_lifo!(tree::BB.SearchTree,node::BB.JumpNode)
    tree.node_counter += 1
    node.id = tree.node_counter
    pushfirst!(tree.nodes, node)
end



function BB.initialize_tree(ins::structures_decomp.instance,main_config::bcp_config) ::BB.SearchTree
# This function is the set up of the BB tree for the SLAP

    # test the instance
    #structures_decomp.save_instance(ins, "test_silva_instance.txt")

    # Find an incumbent
    incumbent = heuristic.main_HC_algo(ins,Dict{String,Any}("iter" => 5000, "shuffle_prob" => .2, "shuffle_max" => .4, "noise" => 0.005))
    # old values 5000, 0.05, 0.2, 5
    # Set up the CG model
    model, rho_vars, a_coef, cut_pool = SLAP_CG_main.set_up_CG(incumbent,main_config.fake_route,main_config.cpu_limit,main_config)
    # Create the search tree with branching strategy
    if main_config.branching_strategy in ["variable", "symmetry", "strong", "aisle", "aisle_without"] 
        tree = BB.initialize_tree(model,BB.VariableBranch)
    elseif main_config.branching_strategy in ["sum_gub", "sum_gl"]
        tree = BB.initialize_tree(model,BB.SumBranch)
    else
        throw(DomainError("Unknown branching strategy"))
    end
    # Update some infos
    tree.current_primal_bound = incumbent.cost
    # Add elements to it
    tree.auxiliary_data["instance"] = incumbent.instance
    tree.auxiliary_data["rho_vars"] = rho_vars
    tree.auxiliary_data["xi_vars"] = model[:xi]
    tree.auxiliary_data["a_coef"] = a_coef
    tree.auxiliary_data["cut_pool"] = cut_pool
    tree.auxiliary_data["branching_strategy"] = main_config.branching_strategy
    tree.auxiliary_data["config"] = main_config

    ############################################################
    ## Parameters that influence the CG
    ############################################################
    tree.auxiliary_data["max_iter_cut"] =  main_config.max_iter_cut
    #tree.auxiliary_data["EPS_cut_separation"] = main_config.EPS_cut_separation_initial
    tree.auxiliary_data["EPS_cut_sep"] = [main_config.EPS_cut_separation_initial for o in incumbent.instance.orderList]
    tree.auxiliary_data["EPS_cut_violated"] = [main_config.EPS_cut_violated_initial for o in incumbent.instance.orderList]
    #tree.auxiliary_data["EPS_cut_violated"] = main_config.EPS_cut_violated_initial
    tree.auxiliary_data["routing_policy"] = main_config.routing_policy

    # Dual bound parameter
    tree.auxiliary_data["min_bound_close_nodes"] = Inf
    tree.auxiliary_data["node_explored"] = 0

    # Best primal solution
    tree.best_primal_solution = incumbent.plan

    return tree
end

function main_BP_function(ins::structures_decomp.instance,main_config::bcp_config)
    # This is the main branch and price function
    # initialize B&B tree and timer
    timer = time()
    tree = BB.initialize_tree(ins,main_config)
    # run the tree
    BB.run(tree, timer + main_config.time_limit)
    # print some stuff
    println("BP finished, completed: $(BB.isempty(tree) && !isempty(tree.processed)), time = $(time() - timer), primal = $(tree.current_primal_bound), dual = $(tree.current_dual_bound), dual root = $(isempty(tree.processed) ? 0 : tree.processed[1].dual_bound)")

    return tree
end
#new method
function main_BP_function(path::String,main_config::bcp_config) 
    if main_config.instance_type == "classic"
        return main_BP_function(structures_decomp.read_instance_file(path,main_config.routing_policy),main_config)
    elseif main_config.instance_type == "silva"
        return main_BP_function(structures_decomp.read_silva_instance_file(path,main_config.routing_policy),main_config)
    else
        throw(DomainError("Unknown instance type"))
    end
end

function BB.print_iteration_summary(tree::BB.SearchTree,node::BB.AbstractNode,message::Dict,start_time)
# This function prints the summary of an iteration
    #println("##### Iteration node $(node.id) #####")
    #println("Summary CG")
    #println(message["dual"])
    #println("## Summary of the node ##")
    #println("Primal : $(node.primal_bound), Dual = $(node.dual_bound)")
    #println(message["branch"])
    #println("## Summary tree ##")
    gap = isinf(tree.current_primal_bound) ? Inf : (tree.current_primal_bound - tree.current_dual_bound) / tree.current_primal_bound
    #println("Node processed = $(length(tree.processed)) counter = $(tree.node_counter), current primal = $(tree.current_primal_bound), current dual = $(tree.current_dual_bound), gap = $gap, time = $(time() - start_time)")
end

end
