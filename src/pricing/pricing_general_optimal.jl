##############################################
# Distance matrix
###############################################
function create_distance_matrix(ins::instance, nodelist::Vector{T},pi_o::Vector{Float64},sigma_o::Vector{Float64},number_blocks_layout::Int) where {T<:AbstractNode}
    if number_blocks_layout == 1
        return create_distance_matrix_single(ins,nodelist,pi_o,sigma_o)
    else
        throw(MethodError)
    end
end

# function that compute the distance matrix
function create_distance_matrix_single(ins::instance, nodelist::Vector{T},pi_o::Vector{Float64},sigma_o::Vector{Float64}) where {T<:AbstractNode}
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

#########################################################################################
# Extension
##########################################################################################

function feasible_extension(path::Partial_path_optimal, new_node::Node)
    # returns true if the extension is feasible
    if new_node.id == 0
        return false
    elseif isinf(dist(path.node, new_node))
        return false
    elseif !path.loc_reachability[new_node.loc.id] 
        # in this case the node is unreachable, but it can still be allowed for a second stop
        # but in this case we must ensure there is still free stops, or it is a remaining
        # mandatory stop
        if ((path.node.loc == new_node.loc) 
            && (new_node.i <= new_node.problem.initial_loc_availability[new_node.loc.id]) 
            && (path.remaining_free_stops > 0 || path.remaining_stops[new_node.loc.id] > 0))
            return true
        else
            return false
        end
    else
        return true
    end
end

function update_loc_reachability(path::Partial_path_optimal,old_node::Node)
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
    # Update the loc reachability in case of no more free stops
    if path.remaining_free_stops == 0
        # then we remove all the extensions apart from the remaining mandatory stops
        for id in eachindex(path.loc_reachability)
            if path.remaining_stops[id] == 0
                path.loc_reachability[id] = false
            end
        end
    end


end

function extend(path::Partial_path_optimal,new_node::Node)
    # This function extend a label to the new_node
    # shallow copy of the new path
    new_path = copy(path)
    # update node and k of new_path
    new_path.node = new_node
    new_path.k += 1
    # update the parameters related to partial
    if new_path.remaining_stops[new_node.loc.id] == 0
        new_path.remaining_free_stops -= 1
    else
        new_path.remaining_stops[new_node.loc.id] -= 1
    end
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

    # now we test if a mandatory stop has become unreachable
    for id in eachindex(new_path.remaining_stops)
        if new_path.remaining_stops[id] != 0 && !new_path.loc_reachability[id] && id != new_path.node.loc.id
            return nothing
        end
    end

    return new_path
end