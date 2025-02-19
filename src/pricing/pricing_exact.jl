module pricing_exact_module

using ..structures_decomp
using ..cost_functions
using ..cuts

###############################################################
# Structures definition
###############################################################

abstract type AbstractNode end
abstract type AbstractPath end


mutable struct Pricing_problem{T<:AbstractNode, P<:AbstractPath}
    ins::instance
    order::order
    node_list::Vector{T}
    ionode::T
    available_cuts::Vector{cuts.OG_cut}
    distance::Array{Float64,2}
    labels::Vector{Vector{P}} #Vector of labels, classified by node they are attached to
    closed_labels::Vector{Vector{P}} # Vector of labels with the good length(s)
    current_k::Int # we keep these info in memory to know we we are in the search
    # in order to know which label should be treated next, and when to stop
    current_node::T # same, which node are we exploring
    mu_o::Float64
    pi_o::Vector{Float64}
    sigma_o::Vector{Float64}
    EPS::Float64
    bi_directional::Bool

    function Pricing_problem{T,P}(ins::instance,o::structures_decomp.order,mu_o::Float64,pi_o::Vector{Float64},sigma_o::Vector{Float64},cut_pool::cuts.CutCollection, bi_directional::Bool, EPS::Float64) where {T <: AbstractNode, P <: AbstractPath}
        pricing_problem = new()
        pricing_problem.ins = ins
        pricing_problem.order = o
        pricing_problem.node_list = create_nodelist(ins,pricing_problem)
        pricing_problem.ionode = T(0,ioloc,0,pricing_problem)
        pricing_problem.available_cuts = cuts.OG_cut[]
        for cut in cut_pool["OG_cut"][o]
            if cut.active && abs(cut.dual_value) > EPS
                push!(pricing_problem.available_cuts, cut)
            end
        end
        pricing_problem.distance = create_distance_matrix(ins,pricing_problem.node_list,pi_o,sigma_o)
        #pricing_problem.open_labels = P[Partial_path{T}(ins,pricing_problem,mu_o,pricing_problem.ionode)]
        pricing_problem.labels = Vector{P}[[] for n1 in pricing_problem.node_list]
        pricing_problem.closed_labels = Vector{P}[[] for n1 in pricing_problem.node_list]
        pricing_problem.current_k = 0
        pricing_problem.current_node = T(0,ioloc,0,pricing_problem)
        pricing_problem.mu_o = mu_o
        pricing_problem.pi_o = pi_o
        pricing_problem.sigma_o = sigma_o
        pricing_problem.EPS = EPS
        if length(o.sku) <= 2
            pricing_problem.bi_directional = false
        else
            pricing_problem.bi_directional = bi_directional
        end

        return pricing_problem
    end
end

mutable struct Partial_path{T<:AbstractNode} <:AbstractPath
    node::T
    rcost::Float64
    k::UInt8
    visited_nodes::Vector{T}
    collected_cuts::Vector{Bool}
    loc_reachability::Vector{Bool}
    cut_reachability::Vector{Bool}

end




# Second constructor that build the first (empty) label associated with the iopoint
function Partial_path{T}(ins::instance,pricing::Pricing_problem, mu_o::Float64,ionode::T) where {T<:AbstractNode}
    node = ionode
    rcost = - mu_o
    k = 0
    visited_nodes = T[node]
    collected_cuts = Bool[false for c in pricing.available_cuts]
    loc_reachability = Bool[true for l in ins.locList]
    cut_reachability = Bool[true for c in pricing.available_cuts]

    return Partial_path{T}(node,rcost,k,visited_nodes,collected_cuts,loc_reachability,cut_reachability)
end

struct Node <:AbstractNode
    id::Int
    loc::loc
    i::Int
    problem::Pricing_problem

    function Node(id::Int,loc::loc,i::Int,pricing_problem::Pricing_problem)
        return new(id,loc,i,pricing_problem)
    end
end

Base.:(==)(a::Node,b::Node) = (a.id == b.id) 
Base.isequal(a::Node,b::Node) = isequal(a.id, b.id)
Base.hash(a::Node, h::UInt) = hash(a.id, hash(:Node,h))

# function that creates the node list
function create_nodelist(ins::instance,pricing_problem::Pricing_problem{T,P}) where {T<:AbstractNode, P <:AbstractPath}
    node_list = Vector{T}()
    id = 1
    for l in ins.locList
        for i = 1:ins.n_sym
            push!(node_list, Node(id,l,i,pricing_problem))
            id += 1
        end
    end

    # UNIT Tests
    @assert length(node_list) == ins.n_sym*length(ins.locList)
    for i in eachindex(node_list)
        n = node_list[i]
        @assert i == n.id
        @assert n.id == (n.loc.id - 1)*ins.n_sym + n.i 
    end

    return node_list
end

# function that compute the distance matrix
function create_distance_matrix(ins::instance, nodelist::Vector{T},pi_o::Vector{Float64},sigma_o::Vector{Float64}) where {T<:AbstractNode}
    distance_matrix = fill(Inf,(length(nodelist)+1,length(nodelist)+1))
    # normal case
    for n1 in nodelist, n2 in nodelist
        if n1 == n2
            continue
        end
        if (n2.i == 1) && (n1.loc != n2.loc)
            distance_matrix[n1.id + 1, n2.id + 1] = distance_locs(n1.loc,n2.loc,ins) - pi_o[n2.loc.id] - sigma_o[n2.loc.id]
        elseif (n1.loc == n2.loc) && (n2.i == n1.i + 1)
            distance_matrix[n1.id + 1, n2.id + 1] = distance_locs(n1.loc,n2.loc,ins) - pi_o[n2.loc.id]
        end
    end
    # with iopoint
    for n in nodelist
        if n.i == 1
            distance_matrix[1,n.id+1] = distance_locs(ioloc,n.loc,ins) - pi_o[n.loc.id] - sigma_o[n.loc.id]
        end
        distance_matrix[n.id+1,1] = distance_locs(n.loc,ioloc,ins)
    end

    # UNIT Tests
    @assert size(distance_matrix,1) == size(distance_matrix,2) == length(nodelist) + 1
    for n1 in nodelist, n2 in nodelist
        # we test if the dist is Inf (not reachable), or n2 is reachable from other nodes (i = 1, and not the same loc)
        # or if they are the same loc with different i
        @assert isinf(distance_matrix[n1.id+1,n2.id+1]) || ((n2.i == 1) && (n1.loc != n2.loc)) || ((n1.loc == n2.loc) && (n2.i == n1.i+1)) 
    end

    return distance_matrix
end


# The distance between 2 locs
function distance_locs(a::loc,b::loc,ins::instance)
    if a.aisle == b.aisle
        return abs(a.col-b.col)*ins.wc
    else
        return ins.wc * min((a.col + b.col) , (2*(ins.cmax + 1) - a.col - b.col)) + abs(a.aisle - b.aisle)*ins.wa
    end
end

# the distance between 2 nodes
function dist(a::Node,b::Node)
    return a.problem.distance[a.id + 1,b.id + 1]
end

#################################################################################
# New methods on base functions with new structs
#################################################################################

# Some base functions on pricing_problem
#Base.isempty(pricing::Pricing_problem{T,P}) where {T<:AbstractNode, P<:AbstractPath} = Base.isempty(pricing.open_labels)
# Define a copy method for partial paths
Base.copy(path::Partial_path{Node}) = Partial_path{Node}(path.node,path.rcost,path.k,copy(path.visited_nodes),copy(path.collected_cuts),copy(path.loc_reachability),copy(path.cut_reachability))


###########################################################################################################
# Main function
###########################################################################################################

function pricing_exact(ins::instance,o::order,mu_o::Float64,pi_o::Array{Float64,1},sigma_o::Array{Float64,1},cut_pool::cuts.CutCollection,bi_directional::Bool,EPS::Float64)
# This is the main function for pricing with exact routing
    # Initialization
    pricing = Pricing_problem{Node,Partial_path{Node}}(ins,o,mu_o,pi_o,sigma_o,cut_pool,bi_directional,EPS)
    # Main loop for label extension
    while !isfinished(pricing)
        # pop the first label of the list
        path = pop_next_label!(pricing)
        # Check if the label is dominated
        if !isdominated!(path,pricing) 
            # Check the length of the label, to put it in closed_labels
            store_label = (!pricing.bi_directional && path.k >= length(o.sku)) || (pricing.bi_directional && path.k >= floor((length(o.sku)) / 2))
            if store_label
                push!(pricing.closed_labels[path.node.id], path)
            end
            # Now extend the label if it is still extandable
            if !done_extending(path)
                for new_node in pricing.node_list
                    # Only if it is a valid successor
                    if feasible_extension(path,new_node)
                        new_path = extend(path,new_node)
                        # insert the label
                        push!(pricing.labels[new_path.node.id], new_path)
                    end
                end
            end
        end
    end
    # UNIT testing -> print the resulting labels
    #=
    for n in pricing.node_list
        k_max1 = isempty(pricing.labels[n.id]) ? 0 : maximum([path.k for path in pricing.labels[n.id]])
        number_labels_1 = k_max1 == 0 ? Int[] : [count(path -> path.k == i, pricing.labels[n.id]) for i = 1:k_max1]
        k_max2 = isempty(pricing.closed_labels[n.id]) ? 0 : maximum([path.k for path in pricing.closed_labels[n.id]])
        number_labels_2 = k_max2 == 0 ? Int[] : [count(path -> path.k == i, pricing.closed_labels[n.id]) for i = 1:k_max2]
        println("node n= $(n.id), number of open labels depending on k = $(number_labels_1), closed labels = $(number_labels_2)")
    end=#
    # Merge the labels of return directly the routes
    if pricing.bi_directional
        routes, rcosts = merge_labels(pricing)
    else
        routes, rcosts = create_routes(pricing)
    end
    # UNIT test if 2 routes are idnetical
    #=
    for rida = 1:(length(routes)-1)
        for ridb = rida+1:length(routes)
            if same_route_test(routes[rida], routes[ridb], pricing.ins)
                println("routea = $([l.id for l in routes[rida].loc])")
                println("routeb = $([l.id for l in routes[ridb].loc])")
            end
            # UNIT we removed this part
            @assert !same_route_test(routes[rida], routes[ridb], pricing.ins)
        end
    end=#
    return routes, rcosts
end

########################################################################
# Basic utility functions
########################################################################

function isfinished(pricing::Pricing_problem{Node,Partial_path{Node}}) 
# This function return true if we are done extending labels (no more open label to extend)
# UNIT extension until the end
    if pricing.bi_directional
        return pricing.current_k > ceil((length(pricing.order.sku)) / 2)
    else
        return pricing.current_k > length(pricing.order.sku) 
    end
end

function pop_next_label!(pricing::Pricing_problem{Node,Partial_path{Node}}) 
    # This function returns the next label to explore, and remove it from the list
    # It also update (potentially) current label list to explore
    # first we code the case of the current node is the ionode
    if pricing.current_node == pricing.ionode
        # pass to the first k at the first node
        pricing.current_k = 1
        pricing.current_node = pricing.node_list[1]
        # return an empty initiated label
        return Partial_path{Node}(pricing.ins,pricing,pricing.mu_o,pricing.ionode)
    end
    # termination condition on k
    if pricing.current_k > (pricing.bi_directional ? ceil((length(pricing.order.sku)) / 2) : length(pricing.order.sku))
        return missing
    end
    # then we check if the current list of label is empty, to avoid bugs
    if isempty(pricing.labels[pricing.current_node.id])
        if pricing.current_node == pricing.node_list[end]
            pricing.current_k += 1
            pricing.current_node = pricing.node_list[1]
            return pop_next_label!(pricing)
        else
            pricing.current_node = pricing.node_list[pricing.current_node.id + 1]
            return pop_next_label!(pricing)
        end
    end 
    # then we check if there is still a label to explore in the current node of length k
    if pricing.labels[pricing.current_node.id][1].k == pricing.current_k
        return popfirst!(pricing.labels[pricing.current_node.id])
    end
    # If that's not the case we look for the next node, or the next k if we are the end of nodelist
    if pricing.current_node == pricing.node_list[end]
        pricing.current_k += 1
        pricing.current_node = pricing.node_list[1]
        # we call the procedure recursively
        return pop_next_label!(pricing)
    else
        pricing.current_node = pricing.node_list[pricing.current_node.id + 1]
        # we call the procedure recursively
        return pop_next_label!(pricing)
    end
end

#=function insert_label(path::Partial_path{T}) where T<:AbstractNode
# This function insert the label in the corresponding array in the pricing problem
    node = path.node
    if node.id != 0
        push!(node.problem.labels[node.id], path)
    end
end=#

function done_extending(path::Partial_path{T}) where T<:AbstractNode
# This function returns true if the label should not be extended further (max length)
    if path.node.problem.bi_directional
        return path.k >= ceil((length(path.node.problem.order.sku)) / 2)
    else
        return path.k >= length(path.node.problem.order.sku) 
    end
end


#######################################################################
# Label extension
#######################################################################

function feasible_extension(path::Partial_path{Node}, new_node::Node)
    # returns true if the extension is feasible
    if new_node.id == 0
        return false
    elseif isinf(dist(path.node, new_node))
        return false
    elseif !path.loc_reachability[new_node.loc.id] && (path.node.loc != new_node.loc)
        return false
    else
        return true
    end
end

function update_loc_reachability(path::Partial_path{Node},old_node::Node)
    # This function updates the loc reachability vector
    # !!!!! IMPORTANT, the node of path should be the new node
    node = path.node
    ins = node.problem.ins
    # first we manage the case of the first extension (from the io node), so that the function doesn't bug
    if old_node.id == 0
        for id = (node.loc.id - node.loc.col + 1) : node.loc.id
            # update it
            path.loc_reachability[id] = false
        end
    # Now we go for the general case
    elseif node.loc.aisle == old_node.loc.aisle
        # put on unreachable the nodes between the two
        for id= min(node.loc.id, old_node.loc.id) : max(node.loc.id, old_node.loc.id)
            path.loc_reachability[id] = false
        end
    elseif (old_node.loc.col + node.loc.col) <= (2*ins.cmax - node.loc.col - old_node.loc.col)
        # In this case we pass by the front cross aisle
        # first the locs in the aisle of old_node
        for id = (old_node.loc.id - old_node.loc.col + 1) : old_node.loc.id
            # update it
            path.loc_reachability[id] = false
        end
        # Then locs in the aisle of new_node
        for id = (node.loc.id - node.loc.col + 1) : node.loc.id
            # update it
            path.loc_reachability[id] = false
        end
    else
        # In this case we pass by the back cross aisle
        # first the locs in the aisle of old_node
        for id = (old_node.loc.id) : (old_node.loc.id + ins.cmax - old_node.loc.col)
            # update it
            path.loc_reachability[id] = false
        end
        # Then locs in the aisle of new_node
        for id = (node.loc.id) : (node.loc.id + ins.cmax - node.loc.col)
            # update it
            path.loc_reachability[id] = false
        end
    end
end

function update_cut_path(path::Partial_path{Node})
    # This function update the fields related to cuts
    available_cuts = path.node.problem.available_cuts
    for i in eachindex(available_cuts)
        # if the cut is already not reachable, nothing to do
        if !path.cut_reachability[i]
            continue
        end
        # test if we add the cut
        if path.node.loc in available_cuts[i].loc_array
            # add the dual cost
            path.rcost -= available_cuts[i].dual_value
            # mark the cut as collected
            path.collected_cuts[i] = true
            # consume the resource
            path.cut_reachability[i] = false
            continue
        end
        # if the cut was not added, test if it is still reachable
        reachable = false
        for loc in available_cuts[i].loc_array
            # for each loc defining the cut, we test if it is still reachable
            if path.loc_reachability[loc.id]
                reachable = true
                break
            end
        end
        if !reachable
            path.cut_reachability[i] = false
        end
    end
end


function extend(path::Partial_path{Node},new_node::Node)
    # This function extend a label to the new_node
    # shallow copy of the new path
    new_path = copy(path)
    # update node and k of new_path
    new_path.node = new_node
    new_path.k += 1
    # update visited locs and loc reachability
    push!(new_path.visited_nodes,new_node)
    update_loc_reachability(new_path,path.node)
    # update cost without cuts
    new_path.rcost += dist(path.node,new_node)
    # update all related to cuts
    update_cut_path(new_path)

    # UNIT testing
    #=
    println("deepcopy: node = $(old_path_deepcopy.node.id), rcost = $(old_path_deepcopy.rcost),
    k = $(old_path_deepcopy.k), visited_nodes = $([n.id for n in old_path_deepcopy.visited_nodes])
    collected_cuts = $(old_path_deepcopy.collected_cuts) 
    loc_reachability = $(old_path_deepcopy.loc_reachability)
    cut_reachability = $(old_path_deepcopy.cut_reachability)")
    println("normal: node = $(path.node.id), rcost = $(path.rcost),
    k = $(path.k), visited_nodes = $([n.id for n in path.visited_nodes])
    collected_cuts = $(path.collected_cuts) 
    loc_reachability = $(path.loc_reachability)
    cut_reachability = $(path.cut_reachability)")=#
    #@assert old_path_deepcopy == path

    return new_path
end

#################################################################
## Label dominance
#################################################################

function isdominated!(path::Partial_path{Node}, pricing::Pricing_problem{Node,Partial_path{Node}})
    # main dominance function
    # first check that the node is not the io point
    if path.node.id == 0
        return false
    end
    # then we look at the dominance in the open labels
    i = 1
    # we removed the second part of the test for UNIT testing
    while i <= length(pricing.labels[path.node.id]) #&& (pricing.labels[1].k <= path.k)
        # we check if the label dominates path, and the opposite
        if dominance_test(path, pricing.labels[path.node.id][i])
            # in this case path dominates the other label, we need to remove it
            deleteat!(pricing.labels[path.node.id], i)
            i -= 1
        elseif dominance_test(pricing.labels[path.node.id][i], path)
            # in this case path is dominated
            return true
        end
        i += 1
    end

    return false
end

# New method for missing argument
isdominated!(path::Missing, pricing::Pricing_problem{Node,Partial_path{Node}}) = true

function dominance_test(patha::Partial_path{Node}, pathb::Partial_path{Node})
    # returns true if path a dominates path b
    @assert patha.node == pathb.node
    # first test if we have the same k
    if patha.k != pathb.k
        return false
    end
    # Then we test the loc reachability
    for i in eachindex(patha.loc_reachability)
        if patha.loc_reachability[i] < pathb.loc_reachability[i]
            return false
        end
    end
    # We test the cost with the improved criteria from Diego for the cuts
    costb = pathb.rcost
    for i in eachindex(patha.cut_reachability)
        if patha.cut_reachability[i] < pathb.cut_reachability[i]
            costb -= patha.node.problem.available_cuts[i].dual_value
        end
    end
    # UNIT testing for cut dominance management
    #=if !prod(patha.cut_reachability .>= pathb.cut_reachability)
        return false
    end=#
    # so we test the dominance on the cost
    if patha.rcost <= costb
        return true
    else
        return false
    end
end

##########################################################################
# Label concatenation
##########################################################################

function merge_labels(pricing::Pricing_problem{Node,Partial_path{Node}})
# This function merges the different half labels in complete paths
    resulting_routes = route[]
    resulting_rcosts = Float64[]

    # we regroup the list of candidate labels for concatenation
    first_candidates, second_candidates = candidate_list(pricing)
            
    # we iterate on the first list of labels
    for first_label_ind in eachindex(first_candidates)
        first_label = first_candidates[first_label_ind]
        for second_label_ind in (iseven(length(pricing.order.sku)) ? (first_label_ind:length(second_candidates)) : eachindex(second_candidates))
            second_label = second_candidates[second_label_ind]
            # we check if the concatenation is 1. feasible and 2. has a rcost < 0
            ok, rcost = try_concatenation(first_label,second_label)
            if ok 
                # we concatenate them
                new_route = concatenation(first_label,second_label)
                # if the rcost if the good one, and the route was not already found we insert the route
                test_and_insert_route!(new_route,rcost,resulting_routes,resulting_rcosts,pricing)
                

                # UNIT test
                # check if the rcost is computed correctly
                #=
                if abs(rcost - rcost_computation_test(new_route,pricing)) > pricing.EPS
                    println("rcost = $(rcost), computed = $(rcost_computation_test(new_route,pricing))")
                    println("locs = $([l.id for l in new_route.loc]) , cost = $(new_route.cost), mu_o = $(pricing.mu_o) 
                    pi_o = $([pricing.pi_o[l] for l in eachindex(pricing.pi_o)]), sigma_o = $([pricing.sigma_o[lid] for lid in eachindex(pricing.sigma_o)]) ")
                end
                @assert abs(rcost - rcost_computation_test(new_route,pricing)) < pricing.EPS=#

            end
        end
    end

    # UNIT test
    # we test that there is no double route
    #=
    for rida = 1:(length(resulting_routes)-1)
        for ridb = rida+1:length(resulting_routes)
            # UNIT we removed this part
            @assert !same_route_test(resulting_routes[rida], resulting_routes[ridb], pricing.ins)
        end
    end=#

    return resulting_routes, resulting_rcosts
end

function candidate_list(pricing::Pricing_problem{Node,Partial_path{Node}})
# Creates the list of candidate labels for concatenation
# The two lists correspond to the two potential lengths of labels
    first_candidates = Partial_path{Node}[]
    second_candidates = Partial_path{Node}[]
    # compute the lengths corresponding to the first list, and the second one
    first_length = ceil((length(pricing.order.sku)) / 2)
    second_length = floor((length(pricing.order.sku)) / 2)
    for node in pricing.node_list
        for label in pricing.closed_labels[node.id]
            # then we check if the label has a good length to be a candidate
            if label.k == first_length
                push!(first_candidates,label)
            end
            if label.k == second_length
                push!(second_candidates, label)
            end
        end 
    end

    # UNIT
    #println("length first = $(length(first_candidates)), second = $(length(second_candidates))")
    #@assert !isempty(first_candidates)
    #@assert !isempty(second_candidates)

    return first_candidates, second_candidates
end

function try_concatenation(first_label::Partial_path{Node},second_label::Partial_path{Node})
# This function provides a quick procedure to see if we try the concatenation of two labels
    pricing = first_label.node.problem
    ins = pricing.ins
    first_locs = [n.loc for n in first_label.visited_nodes[2:end]]
    second_locs = [n.loc for n in second_label.visited_nodes[2:end]]
    # we keep in memory the sigma cost counted twice
    sigma_twice = 0
    for l in intersect(first_locs,second_locs)
        sigma_twice += pricing.sigma_o[l.id]
        count = 0
        for node in first_label.visited_nodes[2:end]
            if node.loc == l
                count +=1
            end
        end
        for node in second_label.visited_nodes[2:end]
            if node.loc == l
                count +=1
            end
        end
        if count > ins.n_sym
            return false, missing
        end
    end
    # Finally, we check if the rcost < 0
    rcost = first_label.rcost + second_label.rcost
    # we add the distance between the two last locations
    rcost += distance_locs(first_label.visited_nodes[end].loc, second_label.visited_nodes[end].loc, ins)
    # we remove the dual costs counted twice
    # first we add the cost associated with mu_o
    rcost += pricing.mu_o
    # then the cost associated with the unique visit of the loc
    rcost += sigma_twice
    # finally we look if we did not count a cut twice 
    for i in eachindex(pricing.available_cuts)
        if first_label.collected_cuts[i] && second_label.collected_cuts[i]
            rcost += pricing.available_cuts[i].dual_value
        end
    end
    
    # we return true if the rcost is < 0
    if rcost < - first_label.node.problem.EPS
        return true, rcost
    else
        return false, missing
    end
end

function concatenation(first_path::Partial_path{Node},second_path::Partial_path{Node})
# This function performs the concatenation
    # we construct the list of locs
    loc_group = loc[]
    for node in first_path.visited_nodes[2:end]
        push!(loc_group, node.loc)
    end
    # then we add the locs of second_path, we do not count the last node, already in first_path
    for node in second_path.visited_nodes[end : -1 : 2]
        push!(loc_group, node.loc)
    end
    # generate the route
    new_route = cost_functions.route_cost(first_path.node.problem.ins, loc_group, first_path.node.problem.order)

    # UNIT tests


    return new_route
end

function test_and_insert_route!(new_route::route,rcost_found::Float64,resulting_routes::Vector{route},resulting_rcosts::Vector{Float64},pricing::Pricing_problem{Node,Partial_path{Node}})
# This function is called when we just generated a new route, to test if we add it in the list
# of resulting routes, and in that case insert it in the good position
# We presuppose the list resulting_routes is already ordered by increasing value of rcost
    # first we compute the rcost
    rcost = rcost_computation_test(new_route,pricing)
    EPS = pricing.EPS
    #UNIT
    @assert rcost < - EPS
    # First case, the rcost was not the computed one
    if abs(rcost_found - rcost) > EPS
        return
    end
    # Second case if the vector of resulting routes is empty
    if isempty(resulting_routes)
        push!(resulting_routes, new_route)
        push!(resulting_rcosts, rcost)
        return
    end
    
    # Normal case we go through the vector of resulting routes until we have a comparable cost
    i = 1
    insertion = true
    insertion_index = 0
    for i = 1:length(resulting_routes)
        if resulting_rcosts[i] < rcost - EPS
            continue
        end
        # we check if it is a good place for insertion
        if (resulting_rcosts[i] >= rcost) && ((i == 1 ? -Inf : resulting_rcosts[i-1]) < rcost)
            insertion_index = i
        end
        # then we check if the route are the same
        if same_route_test(new_route,resulting_routes[i],pricing.ins)
            insertion = false
            break
        end
        # then we have a condition to break the loop at the end
        if resulting_rcosts[i] > rcost + EPS
            break
        end
    end

    # we perform the insertion if need be
    if insertion
        # we check if the insertion index == 0
        if insertion_index == 0
            # UNIT
            @assert resulting_rcosts[end] <= rcost 
            insertion_index = length(resulting_routes) + 1
        end
        insert!(resulting_routes, insertion_index, new_route)
        insert!(resulting_rcosts, insertion_index, rcost)
    end
end


##################################################################
# Test functions
##################################################################

function rcost_computation_test(r::route, pricing_problem::Pricing_problem{Node,Partial_path{Node}})
    # This function computes the rcost of a given route
    rcost = r.cost - pricing_problem.mu_o
    counted_sigma = fill(false,length(pricing_problem.sigma_o))
    # first pi_o
    for id in eachindex(r.loc)
        lid = r.loc[id].id
        # pi_o
        rcost -= pricing_problem.pi_o[lid]
        # sigma_o
        if !counted_sigma[lid]
            rcost -= pricing_problem.sigma_o[lid]
            counted_sigma[lid] = true
        end
    end
    # cuts
    for c in pricing_problem.available_cuts
        if structures_decomp.intersects(c.loc_array,r.loc)
            rcost -= c.dual_value
        end
    end

    return rcost
end

function same_route_test(ra::route,rb::route, ins::instance)
# This function returns true if both routes are identical
    # test order
    if ra.order != rb.order 
        return false
    end
    # test locs
    for l in ins.locList
        # we count if they have the same number of occurences in both routes
        if count(h -> h == l, ra.loc) != count(h -> h == l, rb.loc)
            return false
        end
    end

    return true
end


#################################################################################
# Utility function for UNIT testing
#################################################################################
function Base.:(==)(a::Partial_path{Node},b::Partial_path{Node}) 
    if a.node != b.node 
        return false
    elseif a.rcost != b.rcost
        return false
    elseif a.k != b.k
        return false
    elseif a.visited_nodes != b.visited_nodes
        return false
    elseif a.loc_reachability != b.loc_reachability
        return false
    elseif a.cut_reachability != b.cut_reachability
        return false
    end

    return true
end

function transform_label_into_route(path::Partial_path{Node})
# This function transform the label in a route
    ins = path.node.problem.ins
    return route_cost(ins,[n.loc for n in path.visited_nodes[2:end]], path.node.problem.order)
end

function transform_label_into_route(path_vector::Vector{Partial_path{Node}})
    resulting_routes = structures_decomp.route[]
    for p in path_vector
        route = transform_label_into_route(p)
        push!(resulting_routes,route)
    end
    return resulting_routes
end

function create_routes(pricing::Pricing_problem{Node,Partial_path{Node}})
    # This function takes a pricing problem that is done with the label extension, and return
    # the list of routes, WITHOUT concatenation implemented (i.e. the labels should be extended up to the end)

    route_list = []
    rcost_list = []
    for n in pricing.node_list
        for label in pricing.closed_labels[n.id]
            # we chekc that the label is extended up to the end, and the rcost is < 0
            if done_extending(label) && label.rcost + dist(label.visited_nodes[end],pricing.ionode) < - pricing.EPS
                # first we check that the route is not already present in the list
                new_route = transform_label_into_route(label)
                insert_route = true
                for r in route_list
                    if same_route_test(r, new_route, pricing.ins)
                        insert_route = false
                        break
                    end
                end
                if insert_route
                    push!(route_list, transform_label_into_route(label))
                    push!(rcost_list, rcost_computation_test(route_list[end],pricing))
                end
            end
        end
    end

    return route_list, rcost_list
end



end # end of module