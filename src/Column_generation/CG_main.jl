
module SLAP_CG_main

using JuMP, CPLEX, MathOptInterface

using ..structures_decomp
using ..cost_functions
using ..BB
using ..cuts
using ..heuristic
using ..CG_utils
using ..pricing

MOI = MathOptInterface


function CG_iteration!(tree::BB.SearchTree,node::BB.JumpNode,EPS::Float64,max_time)
# This is the main function called by the branch and bound to get a dual bound
# there is a max time
# the function returns true if the bound is usable (ie a valid dual bound), false otherwise
    # Set up some parameters
    ins = tree.auxiliary_data["instance"]
    master_problem = node.model
    rho_vars = tree.auxiliary_data["rho_vars"]
    xi_vars = tree.auxiliary_data["xi_vars"]
    a_coef = tree.auxiliary_data["a_coef"]
    cut_pool = tree.auxiliary_data["cut_pool"]
    config = tree.auxiliary_data["config"]

    xi_val = Array{Float64,2}(undef,length(ins.locList),ins.n)
    rho_val = Dict{order,Array{Float64,1}}()

    iteration = 1
    n_generated = 0
    time_begin = time()
    done = false
    best_dual = 0
    round_OG = 0
    tree.auxiliary_data["node_explored"] += 1

    last_obj_cut_OG = 0

    time_cg = time()

    # we collect the sku fixed to locs via branching
    fixed_branches = collect_fixed_to_one_branches(node,ins)

    # UNIT
    for o in ins.orderList
        @assert length(a_coef[o]) == length(rho_vars[o])
    end

    while !done #&& iteration < 100 
        # Solve the master problem
        time_before_mp = time()
        #CPLEX.CPXpresolve(unsafe_backend(master_problem).env, unsafe_backend(master_problem).lp, CPX_ALG_DUAL)
        optimize!(master_problem)
        # UNIT counting the number of cuts
        ncut = 0
        for o in ins.orderList
            for c in cut_pool["OG_cut"][o]
                if c.active
                    ncut += 1
                end
            end
        end

        
        # Check time
        if time() > max_time
            return !(best_dual == 0), best_dual, ""
        end
        # check if it is unfeasible
        if JuMP.termination_status(master_problem) == MOI.INFEASIBLE
            return false, best_dual, "INFEASIBLE"
        end
        # Get the value of the different variables
        xi_val = value.(xi_vars)
        rho_val = Dict{order,Array{Float64,1}}()
        for o in ins.orderList
            rho_val[o] = value.(rho_vars[o])
        end
        
        # Try to add cuts from the cut pool
        cut_activated = cuts.update_active_cut_pool!(ins,master_problem,tree.auxiliary_data["max_iter_cut"],xi_vars,xi_val,rho_vars,rho_val,a_coef,cut_pool,EPS,tree.auxiliary_data["EPS_cut_violated"],tree.auxiliary_data["config"].EPS_cut_vio_SL1,tree.auxiliary_data["config"].max_cut_per_order)
        # If successful then start again the iteration
        if cut_activated
            iteration += 1
            #println("cut activated")
            continue
        end

        # Collect the dual variables
        mu = dual.(master_problem[:rho_sum]) + dual.(master_problem[:rho_sum2])
        pi = - dual.(master_problem[:link_constraint])
        # Collect dual variables from cuts
        cuts.collect_dual_variable_cut!(ins,master_problem,cut_pool)
        obj_master = objective_value(master_problem)

        # Now computing the pricing problem for routes because the pricing for plans failed
        reduced_cost_route = Dict{order,Float64}()
        done_route = Bool[]

        # UNIT
        number_generated_routes = Dict{order,Int}()
        routes_generated_all = Dict{order,Vector{route}}(o => [] for o in ins.orderList)

        time_before_pricing = time()
        for o in ins.orderList
            if time() > max_time
                push!(done_route,false)
                break
            end
            # Aggregate pi dual variables
            pi_o = pi[:,o.id]
            # Compute the sigma dual costs
            sigma_o = [0. for l in ins.locList]
            # first with the Linkcuts
            for l in ins.locList, s in o.sku
                cut = cut_pool["LinkCut"][l.id,o.id,s.id]
                if cut.active
                    sigma_o[cut.loc.id] += cut.dual_value
                end
            end
            # Then with the landp cuts
            #cuts.collect_dual_variables_landp!(o,sigma_o,tree)
            # Solve the pricing
            routes_generated, rcost_generated = solve_pricing(ins,o,mu[o.id],pi_o,sigma_o,cut_pool,tree.auxiliary_data["config"].bi_directional,fixed_branches,EPS)

            reduced_cost_route[o] = (length(rcost_generated) > 0 ? minimum(rcost_generated) : 0)
            push!(done_route, (reduced_cost_route[o] >= - EPS))
            routes_generated_all[o] = routes_generated

            # Verify if all reduced cost are <= 0
            if (length(rcost_generated) > 0) &&  (maximum(rcost_generated) >= 0) #&& false
                println("Error in the reduced cost computation")
            end
            # Count the number of routes
            n_generated += length(routes_generated)
            number_generated_routes[o] = length(routes_generated)
        end
        for o in ins.orderList
            # Add the routes to the master problem and register them
            for r in routes_generated_all[o]
                add_route_to_master!(r,master_problem,rho_vars[o],cut_pool)
                push!(a_coef[o],r)
            end
        end

        iteration += 1
        done = prod(done_route)
        


        # Generate cuts
        # UNIT for only looking for cuts if depth <= 4
        depth = 1
        n = node
        while !isnothing(n.parent)
            depth += 1
            n = n.parent
        end
        if done 
            # First we have a valid dual bound
            if obj_master > best_dual
                best_dual = obj_master
            end
            #println("before cut separation")
            # UNIT cut only with depth <= 4
            if (depth <= config.max_depth_cut_generation || tree.auxiliary_data["node_explored"] <= config.min_nodes_with_cuts) && (round_OG < config.max_round_OG) && (obj_master - last_obj_cut_OG > config.EPS_cut_obj_improvement) && ((node.id == 1 && time() - time_cg < 400) || time() - time_cg < 120)
                test, test_val, cuts_found = cuts.separate_cut_all!(ins,cut_pool,xi_val,rho_val,a_coef,tree.auxiliary_data["EPS_cut_sep"])
                round_OG += 1
                #println("after sep OG, found = $cuts_found, number = $(length(test)), cut_pool = $([length(cut_pool["OG_cut"][o]) for o in ins.orderList])")
                # register the value of the objective function
                last_obj_cut_OG = obj_master
                if cuts_found
                    #activate the cuts
                    activated = cuts.update_active_cut_pool!(ins,master_problem,tree.auxiliary_data["max_iter_cut"],xi_vars,xi_val,rho_vars,rho_val,a_coef,cut_pool,EPS,tree.auxiliary_data["EPS_cut_violated"],tree.auxiliary_data["config"].EPS_cut_vio_SL1,tree.auxiliary_data["config"].max_cut_per_order)
                    done = !activated
                end
            end
        end
            

        # Deactivate cuts
        if !done
            cuts.deactivate_cut_pool!(ins,master_problem,tree.auxiliary_data["max_iter_cut"],tree.auxiliary_data["config"].max_iter_SL1,xi_vars,xi_val,rho_vars,rho_val,a_coef,cut_pool,EPS)
        end
    end

    @assert CG_utils.checker(ins,master_problem,xi_val,rho_val,a_coef,EPS)
    return true, best_dual, "The CG performed $(iteration) iterations, in $(time() - time_begin) seconds, and generated $(n_generated) routes."

end # end of function



function set_up_CG(incumbent,add_fake_route,cpu_limit,main_config)
# This function set up the CG
    ins = incumbent.instance
    # Defining empty arrays for rho variables, one array for each order.
    rho_vars = Dict{order, Array{VariableRef,1}}()
    for o in ins.orderList
        rho_vars[o] = VariableRef[]
    end
    # To keep track of columns
    a_coef = Dict{order,Array{route,1}}(o => route[] for o in ins.orderList)

    # Set up cuts
    cut_pool = cuts.initialization_cuts(ins,main_config)

    # master
    master_problem = CG_utils.create_master_problem(ins,cpu_limit,main_config.simplex_algo)

    # Create fake routes
    if add_fake_route
        fake_route = CG_utils.build_fake_routes(ins)
        # Keep in memory
        add_route_to_master!(fake_route, master_problem, rho_vars, cut_pool)
        for r in fake_route
            push!(a_coef[r.order], r)
        end
    end
    
    # routing variables
    add_route_to_master!(incumbent.routes, master_problem, rho_vars,cut_pool)
    
    #keep in memory
    for r in incumbent.routes
        push!(a_coef[r.order], r)
    end

    # UNIT
    # activate all the link cuts
    if main_config.activate_all_link_cuts
        cuts.activate_all_link_cut!(ins,master_problem,cut_pool,rho_vars,a_coef,master_problem[:xi])
    end
    # initiate SL cut on aisle
    if main_config.initiate_aisle_cuts_in_pool
        cuts.initiate_aisle_cuts!(ins,cut_pool,a_coef)
    end


    return master_problem, rho_vars, a_coef, cut_pool
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

end#end module

#SLAP_ColGen_compact.launch()
