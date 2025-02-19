
# This is the branching functions

# classical branching
function branching_variable!(tree::BB.SearchTree,node::BB.JumpNode{BB.VariableBranch})
    # First branching function, branch on the most fractionate variable with a coefficient for the demand
    ins = tree.auxiliary_data["instance"]
    xi = node.model[:xi]
    # Find the best variable to branch on
    l_best = 0
    s_best = 0
    best_val = 0
    for l=1:size(xi,1)
        for s=1:size(xi,2)
            if min(value(xi[l,s]), 1 - value(xi[l,s])) * ins.skuList[s].demand > best_val
                l_best = l
                s_best = s
                best_val = min(value(xi[l,s]), 1 - value(xi[l,s])) * ins.skuList[s].demand
            end
        end
    end
    
    ### Tests
    branch_objects = BB.collect_all_branches(node)
    @assert(best_val > EPS)
    
    # Now create the two branches
    # Set the variable to zero
    child1 = BB.create_child_node_with_ub(node, xi[l_best,s_best], 0)
    insert_node_queue!(tree,child1)
    # Set the variable to 1
    child2 = BB.create_child_node_with_lb(node, xi[l_best,s_best], 1)
    insert_node_queue!(tree,child2)

    return "Branching on variable xi[$l_best,$s_best], value = $(value(xi[l_best,s_best]))"
end


function BB.branch!(tree::BB.SearchTree,node::BB.JumpNode{T}) where { T <: BB.AbstractBranch}
# This function branch with the chosen strategy
    # First check the strategy
    if tree.auxiliary_data["branching_strategy"] == "symmetry"
        return branching_symmetry!(tree,node)
    elseif tree.auxiliary_data["branching_strategy"] == "aisle"
        return branching_aisle!(tree,node)
    elseif tree.auxiliary_data["branching_strategy"] == "variable"
        return branching_variable!(tree,node)
    elseif tree.auxiliary_data["branching_strategy"] == "sum_gub"
        return branching_sum_gub!(tree,node)
    elseif tree.auxiliary_data["branching_strategy"] == "sum_gl"
        return branching_sum_gl!(tree,node)
    elseif tree.auxiliary_data["branching_strategy"] == "strong"
        return strong_branching!(tree,node)
    elseif tree.auxiliary_data["branching_strategy"] == "aisle_without"
        return branching_aisle_without!(tree,node)
    else
        throw(DomainError("Unknown branching strategy"))
    end
end

function heuristic_branching_partition(set::Array{T,1}, values::Array{Float64,1}, offset::Float64) where T <: Union{structures_decomp.loc,structures_decomp.sku}
# This function creates a partition of a given set, trying to create a partition of the given set while
# trying to balance the values equally between the 2 sets, with an offset for the first one
    @assert length(set) == length(values)
    # first order the set
    perm = sortperm(values,rev=true)
    # preparation
    seta = T[]
    valuea = 0.
    valueb = offset
    for i=1:length(set)
        if valuea < valueb
            push!(seta, set[perm[i]])
            valuea += values[perm[i]]
        else
            valueb += values[perm[i]]
        end
    end

    return seta, valuea
end


#######################################################
# Sum loc branching - GUB branching
#######################################################

function branching_sum_gub!(tree,node)
# This is the main function for branching when the strategy is GUB branching
    # setup
    ins = tree.auxiliary_data["instance"]
    xi = node.model[:xi]
    xi_val = value.(xi)
    best_s = nothing::Union{Nothing,structures_decomp.sku}
    best_score = 0.
    set_value_registered = 0. # only for debugging purpose
    best_set = structures_decomp.loc[]
    # main loop, for each s compute the partition and check if it's better than the current one
    for s in ins.skuList
        # only if the demand is greater than 0
        if s.demand == 0
            continue
        end
        # compute the heuristic partition for this s
        current_set, set_value = heuristic_branching_partition(ins.locList,xi_val[:,s.id],0.)
        if set_value < EPS || set_value > 1 - EPS
            continue
        end
        # compute the score of the partition
        demand_term = s.demand / maximum(t.demand for t in ins.skuList)
        xi_number_term = min(count(l -> xi_val[l,s.id] > EPS, 1:size(xi_val,1)) / 10 , 1)
        balance_term = min(min(count(l -> xi_val[l.id,s.id] > EPS, current_set) , xi_number_term - count(l -> xi_val[l.id,s.id] > EPS, current_set)) / 5 , 1)
        closeness_term = min(set_value, 1 - set_value)
        current_score = tree.auxiliary_data["config"].branching_weight_gub_demand*demand_term + tree.auxiliary_data["config"].branching_weight_gub_xinumber*xi_number_term + tree.auxiliary_data["config"].branching_weight_gub_balance*balance_term + tree.auxiliary_data["config"].branching_weight_gub_closeness*closeness_term
        # register the score if best
        if current_score > best_score
            best_score = current_score
            best_set = current_set
            best_s = s
            set_value_registered = set_value
        end
    end
    # do the branhing 
    branch_objects = BB.collect_all_branches(node)
    
    # Now create the two branches
    @assert !isnothing(best_s)
    # Set the variable to zero
    variable_sum1 = JuMP.VariableRef[xi[l.id,best_s.id] for l in best_set]
    child1 = BB.create_child_node(node, variable_sum1, nothing, 0)
    insert_node_queue!(tree,child1)
    # Set the variable to 1
    variable_sum2 = JuMP.VariableRef[xi[l.id,best_s.id] for l in ins.locList if !(structures_decomp.intersects_id([l],best_set))]
    child2 = BB.create_child_node(node, variable_sum2, nothing, 0)
    insert_node_queue!(tree,child2)

    return "Branching on GUB for sku s = $(best_s.id), score = $best_score, sum value $(set_value_registered), length best set = $(length(best_set))"
end

#############################################################################
# Generalized loc branching GL
#############################################################################

function branching_sum_gl!(tree,node)
    # This is the main function for branching when the strategy is GUB branching
        # setup
        ins = tree.auxiliary_data["instance"]
        xi = node.model[:xi]
        xi_val = value.(xi)
        best_l = nothing::Union{Nothing,structures_decomp.sku}
        best_score = 0.
        set_value_registered = 0. # only for debugging purpose
        best_set = structures_decomp.sku[]
        # main loop, for each s compute the partition and check if it's better than the current one
        for l in ins.locList
            # compute the heuristic partition for this s
            current_set, set_value = heuristic_branching_partition(ins.skuList,xi_val[l.id,:],-1.)
            if set_value > 1 + EPS
                println("set value = $set_value")
            end
            @assert set_value <= 1 + EPS
            if set_value < EPS || set_value > 1 - EPS
                continue
            end
            # compute the score of the partition
            demand_term = sum(s.demand*xi_val[l.id,s.id] for s in ins.skuList) / ( 2*maximum(s.demand for s in ins.skuList))
            xi_number_term = min(count(s -> xi_val[l.id,s] > EPS, 1:size(xi_val,2)) / 10 , 1)
            distance_term = cost_functions.dist(ins,structures_decomp.ioloc,l) / (ins.wa*ins.amax + ins.wc*ins.cmax)
            closeness_term = min(set_value, 1 - set_value)
            current_score = tree.auxiliary_data["config"].branching_weight_gl_demand*demand_term + tree.auxiliary_data["config"].branching_weight_gl_xinumber*xi_number_term + tree.auxiliary_data["config"].branching_weight_gl_distance*distance_term + tree.auxiliary_data["config"].branching_weight_gl_closeness*closeness_term
            # register the score if best
            if current_score > best_score
                best_score = current_score
                best_set = current_set
                best_l = l
                set_value_registered = set_value
            end
        end
    # do the branhing 
    branch_objects = BB.collect_all_branches(node)
        
        # Now create the two branches
    @assert !isnothing(best_l)
    # Set the variable to zero
    variable_sum1 = JuMP.VariableRef[xi[best_l.id,s.id] for s in best_set]
    child1 = BB.create_child_node(node, variable_sum1, nothing, 0)
    insert_node_queue!(tree,child1)
    # Set the variable to 1
    variable_sum2 = copy(variable_sum1)
    child2 = BB.create_child_node(node, variable_sum2, 1, nothing)
    insert_node_queue!(tree,child2)
    
    return "Branching on GL for loc l = $(best_l.id), score = $best_score, sum value $(set_value_registered), length best set = $(length(best_set))"
end


function branching_symmetry!(tree,node)
    # This is the branching function symmetry branching
    ins = tree.auxiliary_data["instance"]
    xi = node.model[:xi]

    l_best = 0
    s_best = 0
    best_val = 0
    for l=1:size(xi,1)
        for s=1:size(xi,2)
            if min(value(xi[l,s]), 1 - value(xi[l,s])) * ins.skuList[s].demand > best_val
                l_best = l
                s_best = s
                best_val = min(value(xi[l,s]), 1 - value(xi[l,s])) * ins.skuList[s].demand
            end
        end
    end

    # we find the corresponding symmetry group
    best_sku = ins.skuList[s_best]
    best_sym = nothing
    for sym in ins.sku_symmetry_list
        if best_sku in sym.skus
            best_sym = sym
            break
        end
    end
    @assert !isnothing(best_sym)

    # now we create the group of variables to set to 0
    if best_sym.size == 1
        zero_var_list = JuMP.VariableRef[xi[l_best,s_best]]
    else
        #println("node id == $(node.id)")
        zero_var_list = create_zero_var_list(node,l_best, best_sym)
    end
        


    # UNIT 
    @assert best_val > EPS
    # now we check that the value is good
    #@assert value(xi[l_best,s_best]) > 0.05
    #@assert value(xi[l_best,s_best]) < 0.95
    # we assert that the s sku has not been fixed in the previous branches
    for branch in BB.collect_all_branches(node)
        #@assert !(branch.fixed_sku_id in [s.id for s in best_sym.skus])
        @assert branch.fixed_sku_id != s_best
    end

    # first the one with the group of variables set to 0
    child1 = BB.create_child_node_multiple0(node, zero_var_list)
    insert_node_queue!(tree,child1)
    # then the branch with the variable set to 1
    child2 = BB.create_child_node_fixed_one(node, xi[l_best,s_best], s_best, l_best)
    insert_node_queue!(tree,child2)

    return "node id = $(node.id) Branching on symmetry id = $(best_sku.id), loc = $(l_best) size = $(best_sym.size)"
end


#####################################################################
# Aisle branching
#####################################################################

function branching_aisle!(tree,node)
    # This is the branching function for branching on aisle
    # first we test if we should branch on aisle, then if the score is too low we branch on a single variable

    ins = tree.auxiliary_data["instance"]
    xi = node.model[:xi]
    xi_val = JuMP.value.(xi)
    aisle_list = [structures_decomp.loc[] for a=1:ins.amax]
    for l in ins.locList
        push!(aisle_list[l.aisle], l)
    end

    # first we look at aisle branching
    a_best = 0
    s_best = 0
    best_val = 0
    for a = 1:ins.amax
        for s=1:size(xi,2)
            val = sum(xi_val[l.id,s] for l in aisle_list[a]) 
            if min(val, 1 - val) * ins.skuList[s].demand > best_val
                a_best, s_best, best_val = a, s, min(val, 1 - val) * ins.skuList[s].demand
            end
        end
    end

    # if the score is too low, we branch on a single variable instead
    if best_val < 0.25
        return branching_symmetry!(tree,node)
    end

    # we find the corresponding symmetry group
    best_sku = ins.skuList[s_best]
    best_sym = nothing
    for sym in ins.sku_symmetry_list
        if best_sku in sym.skus
            best_sym = sym
            break
        end
    end
    @assert !isnothing(best_sym)

    # first we create the branch where the sku is NOT in the chosen aisle
    zero_var_list1 = create_zero_var_list(node, aisle_list[a_best], best_sym)
    child1 = BB.create_child_node_multiple0(node, zero_var_list1)
    insert_node_queue!(tree,child1)

    # then we create the branch where the sku IS FIXED in the chosen aisle
    zero_var_list2 = JuMP.VariableRef[]
    for l in ins.locList
        if l.aisle != a_best
            push!(zero_var_list2, xi[l.id,s_best])
        end
    end
    child2 = BB.create_child_node_multiple0(node, zero_var_list2, s_best)
    insert_node_queue!(tree,child2)

    return "node id = $(node.id) Branching on sku id = $(best_sku.id), aisle = $(a_best) size = $(best_sym.size)"
end


#####################################################################
# aisle branching without sym
#####################################################################

function branching_aisle_without!(tree,node)
    # This is the branching function for branching on aisle
    # first we test if we should branch on aisle, then if the score is too low we branch on a single variable

    ins = tree.auxiliary_data["instance"]
    xi = node.model[:xi]
    xi_val = JuMP.value.(xi)
    aisle_list = [structures_decomp.loc[] for a=1:ins.amax]
    for l in ins.locList
        push!(aisle_list[l.aisle], l)
    end

    # first we look at aisle branching
    a_best = 0
    s_best = 0
    best_val = 0
    for a = 1:ins.amax
        for s=1:size(xi,2)
            val = sum(xi_val[l.id,s] for l in aisle_list[a]) 
            if min(val, 1 - val) * ins.skuList[s].demand > best_val
                a_best, s_best, best_val = a, s, min(val, 1 - val) * ins.skuList[s].demand
            end
        end
    end

    # if the score is too low, we branch on a single variable instead
    if best_val < 0.25
        return branching_variable!(tree,node)
    end

    # we find the corresponding symmetry group
    best_sku = ins.skuList[s_best]

    # first we create the branch where the sku is NOT in the chosen aisle
    zero_var_list1 = JuMP.VariableRef[]
    for l in ins.locList
        if l.aisle == a_best
            push!(zero_var_list1, xi[l.id,s_best])
        end
    end
    child1 = BB.create_child_node_multiple0(node, zero_var_list1)
    insert_node_queue!(tree,child1)

    # then we create the branch where the sku IS FIXED in the chosen aisle
    zero_var_list2 = JuMP.VariableRef[]
    for l in ins.locList
        if l.aisle != a_best
            push!(zero_var_list2, xi[l.id,s_best])
        end
    end
    child2 = BB.create_child_node_multiple0(node, zero_var_list2)
    insert_node_queue!(tree,child2)

    return "node id = $(node.id) Branching on sku id = $(best_sku.id), aisle = $(a_best)"
end

####################################################################
# Strong branching
####################################################################


function strong_branching!(tree,node)
    # This is the branching function symmetry branching
    ins = tree.auxiliary_data["instance"]
    config = tree.auxiliary_data["config"]
    xi = node.model[:xi]
    xi_val = value.(xi)

    # first we order the potential choices by order of the value score, in decreasing order
    candidate_list = [] 
    for l=1:size(xi,1), s=1:size(xi,2)
        push!(candidate_list,(l,s))
    end
    sort!(candidate_list, by= c -> min(xi_val[c[1],c[2]], 1 - xi_val[c[1],c[2]]) * ins.skuList[c[2]].demand, rev=true)

    # then we look for the xxx first ones and compute an approximate dual bound for each
    candidate_bound = []
    candidate_bound_zero = []
    candidate_bound_one = []
    for (l,s) in candidate_list[1:config.nb_candidates]
        if xi_val[l,s] < EPS || xi_val[l,s] > 1 - EPS
            push!(candidate_bound,0)
            continue
        end
        # first we add the constraint equals to 0
        # we find the corresponding symmetry group
        best_sym = nothing
        for sym in ins.sku_symmetry_list
            if ins.skuList[s] in sym.skus
                best_sym = sym
                break
            end
        end
        zero_var_list = create_zero_var_list(node,l,best_sym)
        for v in zero_var_list
            JuMP.set_upper_bound(v, 0.)
        end
        dual_bound_zero = approximate_dual_bound(tree,node,nothing,nothing)
        for v in zero_var_list
            JuMP.set_upper_bound(v, 1.)
        end

        # then we add the constraint equals to 1
        JuMP.set_lower_bound(xi[l,s], 1.)
        dual_bound_one = approximate_dual_bound(tree,node,ins.skuList[s],ins.locList[l])
        JuMP.set_lower_bound(xi[l,s], 0.)
        
        # register the bound
        push!(candidate_bound, 0.99*min(dual_bound_zero, dual_bound_one) + 0.01* max(dual_bound_zero, dual_bound_one))
        push!(candidate_bound_zero, dual_bound_zero)
        push!(candidate_bound_one, dual_bound_one)
    end

    # now we take the best candidate
    @show candidate_list[1:10]
    #@show candidate_bound_zero
    #@show candidate_bound_one
    @show findmax(candidate_bound)
    best_ind = findmax(candidate_bound)[2]
    best_loc, best_sku = candidate_list[best_ind]

    # we find the corresponding symmetry group
    best_sym = nothing
    for sym in ins.sku_symmetry_list
        if ins.skuList[best_sku] in sym.skus
            best_sym = sym
            break
        end
    end
    @assert !isnothing(best_sym)

    # now we create the group of variables to set to 0
    if best_sym.size == 1
        zero_var_list = JuMP.VariableRef[xi[best_loc,best_sku]]
    else
        zero_var_list = create_zero_var_list(node,best_loc, best_sym)
    end
        


    # UNIT 
    # now we check that the value is good
    #@assert value(xi[l_best,s_best]) > 0.05
    #@assert value(xi[l_best,s_best]) < 0.95
    # we assert that the s sku has not been fixed in the previous branches
    for branch in BB.collect_all_branches(node)
        #@assert !(branch.fixed_sku_id in [s.id for s in best_sym.skus])
        @assert branch.fixed_sku_id != best_sku
    end

    # first the one with the group of variables set to 0
    child1 = BB.create_child_node_multiple0(node, zero_var_list)
    insert_node_queue!(tree,child1)
    # then the branch with the variable set to 1
    child2 = BB.create_child_node_fixed_one(node, xi[best_loc,best_sku], best_sku, best_loc)
    insert_node_queue!(tree,child2)

    return "node id = $(node.id) Branching on symmetry id = $(best_sku), loc = $(best_loc) size = $(best_sym.size)"
end


function approximate_dual_bound(tree,node,fixed_sku::Union{Nothing,structures_decomp.sku},fixed_loc::Union{Nothing,structures_decomp.loc})
    # This function returns an estimate of the dual bound for strong branching
    # first there is only one possible option -> basic where we just check the bound without the addition of columns or constraints
    if tree.auxiliary_data["config"].bound_method == "basic"
        return approximate_dual_bound_basic(tree,node)
    elseif tree.auxiliary_data["config"].bound_method == "cg"
        return approximate_dual_bound_cg(tree,node,fixed_sku,fixed_loc)
    else
        error("Unknown bound method")
    end
end

function approximate_dual_bound_basic(tree,node)
    # This function is the basic evaluation, it just solves the master problem and return the result
    mp = node.model
    optimize!(mp)
    return objective_value(mp)
end


function approximate_dual_bound_cg(tree,node,fixed_sku::Union{Nothing,structures_decomp.sku},fixed_loc::Union{Nothing,structures_decomp.loc})
    # This function is an advanced function to evaluate the dual bound (it actually provide a valid bound)
    # we generate new columns and SL1 cuts until convergence

    # Set up some parameters
    ins = tree.auxiliary_data["instance"]
    master_problem = node.model
    rho_vars = tree.auxiliary_data["rho_vars"]
    xi_vars = tree.auxiliary_data["xi_vars"]
    a_coef = tree.auxiliary_data["a_coef"]
    cut_pool = tree.auxiliary_data["cut_pool"]
    config = tree.auxiliary_data["config"]

    # we collect the sku fixed to locs via branching
    fixed_branches = SLAP_CG_main.collect_fixed_to_one_branches(node,ins)
    if !isnothing(fixed_sku)
        fixed_branches[fixed_sku] = fixed_loc
    end

    done = false
    while !done 
        optimize!(master_problem)
        # Get the value of the different variables
        xi_val = value.(xi_vars)
        rho_val = Dict{structures_decomp.order,Array{Float64,1}}()
        for o in ins.orderList
            rho_val[o] = value.(rho_vars[o])
        end
        
        # Try to add cuts from the cut pool
        cut_activated = cuts.update_active_cut_pool!(ins,master_problem,tree.auxiliary_data["max_iter_cut"],xi_vars,xi_val,rho_vars,rho_val,a_coef,cut_pool,EPS,tree.auxiliary_data["EPS_cut_violated"])
        # If successful then start again the iteration
        if cut_activated
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
        reduced_cost_route = Dict{structures_decomp.order,Float64}()
        done_route = Bool[]

        # UNIT
        number_generated_routes = Dict{structures_decomp.order,Int}()
        routes_generated_all = Dict{structures_decomp.order,Vector{structures_decomp.route}}(o => [] for o in ins.orderList)

        for o in ins.orderList
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
            # Solve the pricing
            routes_generated, rcost_generated = pricing.solve_pricing(ins,o,mu[o.id],pi_o,sigma_o,cut_pool,tree.auxiliary_data["config"].bi_directional,fixed_branches,EPS)

            reduced_cost_route[o] = (length(rcost_generated) > 0 ? minimum(rcost_generated) : 0)
            push!(done_route, (reduced_cost_route[o] >= - EPS))
            routes_generated_all[o] = routes_generated

            # Verify if all reduced cost are <= 0
            if (length(rcost_generated) > 0) &&  (maximum(rcost_generated) >= 0) #&& false
                println("Error in the reduced cost computation")
            end
        end
        for o in ins.orderList
            # Add the routes to the master problem and register them
            for r in routes_generated_all[o]
                CG_utils.add_route_to_master!(r,master_problem,rho_vars[o],cut_pool)
                push!(a_coef[o],r)
            end
        end

        done = prod(done_route)
    end

    return objective_value(master_problem)

end



function create_zero_var_list(node,l_best::Int,best_sym)
    # This function creates the list of variables to be set to zero in the case of symmetry branching
    xi = node.model[:xi]
    zero_var_list = JuMP.VariableRef[xi[l_best,s.id] for s in best_sym.skus]
    # now explore the other branches of to see if one sku of the sym group has already been set to 1
    branches = BB.collect_all_branches(node)
    for b in branches
        # we check that there is a lb in the branch (i.e. a variable is put at one)
        # and if this variable (because it is a single one) is in the list of variables to put to zero
        if !isempty(b.lb) 
            @assert length(b.lb) == 1
            for branch_variable in keys(b.lb)
                if branch_variable in zero_var_list
                    filter!(var -> var != branch_variable, zero_var_list)
                end
            end
        end
    end

    return zero_var_list

end

function create_zero_var_list(node,l_best::Vector{structures_decomp.loc},best_sym)
    # This function creates the list of variables to be set to zero in the case of aisle branching, new method
    # with multiple locs
    xi = node.model[:xi]
    sku_list_sym = copy(best_sym.skus)
    # first we check if there is a sku that has already been fixed to an aisle in the symmetry group
    branches = BB.collect_all_branches(node)
    for b in branches
        if !isnothing(b.fixed_sku_aisle)
            filter!(s -> s.id != b.fixed_sku_aisle, sku_list_sym)
        end
    end

    # we create the list of variables, with the updated sku_list_sym that accounts for the previously fixed to aisle skus
    zero_var_list = JuMP.VariableRef[]
    for l in l_best, s in sku_list_sym
        push!(zero_var_list, xi[l.id,s.id])
    end

    # now we can remove the variables that have been fixed to one as usual
    for b in branches
        # first we check is there is a lb in the branch (i.e. the variable has been fixed to 1 by a variable branching)
        if !isempty(b.lb) 
            @assert length(b.lb) == 1
            for branch_variable in keys(b.lb)
                if branch_variable in zero_var_list
                    filter!(var -> var != branch_variable, zero_var_list)
                end
            end
        end

    end

    return zero_var_list

end


            