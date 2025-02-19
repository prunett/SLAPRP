##############################################
# Distance matrix
###############################################

# function that compute the distance matrix
function create_distance_matrix_sshape(ins::instance, nodelist::Vector{T},pi_o::Vector{Float64},sigma_o::Vector{Float64}) where {T<:AbstractNode}
    distance_matrix_from_bottom = fill(Inf,(length(nodelist)+1,length(nodelist)+1))
    distance_matrix_from_top = fill(Inf,(length(nodelist)+1,length(nodelist)+1))
    # normal case
    for n1 in nodelist, n2 in nodelist
        if n1 == n2
            continue
        end
        if (n2.i == 1) && (n1.loc != n2.loc)
            # case same aisle
            if n1.loc.aisle == n2.loc.aisle
                if n1.loc.col < n2.loc.col
                    distance_matrix_from_bottom[n1.id + 1,n2.id + 1] = distance_locs(n1.loc,n2.loc,ins) - pi_o[n2.loc.id] - sigma_o[n2.loc.id]
                else
                    distance_matrix_from_top[n1.id + 1,n2.id + 1] = distance_locs(n1.loc,n2.loc,ins) - pi_o[n2.loc.id] - sigma_o[n2.loc.id]
                end
            # case where we go forward in the aisles (the other is Inf)
            elseif n1.loc.aisle < n2.loc.aisle
                distance_matrix_from_bottom[n1.id + 1,n2.id + 1] = dist_locs_by_top(n1.loc,n2.loc,ins) - pi_o[n2.loc.id] - sigma_o[n2.loc.id]
                distance_matrix_from_top[n1.id + 1,n2.id + 1] = dist_locs_by_bottom(n1.loc,n2.loc,ins) - pi_o[n2.loc.id] - sigma_o[n2.loc.id]
            end
        elseif (n1.loc == n2.loc) && (n2.i == n1.i + 1)
            distance_matrix_from_bottom[n1.id + 1, n2.id + 1] = - pi_o[n2.loc.id]
            distance_matrix_from_top[n1.id + 1, n2.id + 1] = - pi_o[n2.loc.id]
        end
    end
    # with iopoint
    for n in nodelist
        if n.i == 1
            distance_matrix_from_bottom[1,n.id+1] = distance_locs(ioloc,n.loc,ins) - pi_o[n.loc.id] - sigma_o[n.loc.id]
        end
        distance_matrix_from_bottom[n.id+1,1] = distance_locs(n.loc,ioloc,ins)
        distance_matrix_from_top[n.id+1,1] = distance_locs(n.loc,ioloc,ins)
    end

    # UNIT Tests
    #@assert size(distance_matrix,1) == size(distance_matrix,2) == length(nodelist) + 1
    #=for n1 in nodelist, n2 in nodelist
        # we test if the dist is Inf (not reachable), or n2 is reachable from other nodes (i = 1, and not the same loc)
        # or if they are the same loc with different i
        @assert isinf(distance_matrix[n1.id+1,n2.id+1]) || ((n2.i == 1) && (n1.loc != n2.loc)) || ((n1.loc == n2.loc) && (n2.i == n1.i+1)) 
    end=#

    return distance_matrix_from_bottom, distance_matrix_from_top
end

# the distance between 2 nodes
function dist(p::Partial_path_sshape,b::Node)
    if p.entered_by_bottom
        return b.problem.distance_from_bottom[p.node.id + 1,b.id + 1]
    else
        return b.problem.distance_from_top[p.node.id + 1,b.id + 1]
    end
end

#########################################################################################
# Extension
##########################################################################################

function feasible_extension(path::Partial_path_sshape, new_node::Node)
    # returns true if the extension is feasible
    if new_node.id == 0
        return false
    elseif isinf(dist(path, new_node))
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

function update_loc_reachability(path::Partial_path_sshape,old_node::Node)
    # This function updates the loc reachability vector
    # !!!!! IMPORTANT, the node of path should be the new node
    node = path.node
    ins = node.problem.ins
    # Like return, no need for reachability
    if old_node.id == 0
        for id = 1 : node.loc.id
            # update it
            path.loc_reachability[id] = false
        end
    elseif node.loc.aisle == old_node.loc.aisle
        # put on unreachable the nodes between the two
        for id= min(node.loc.id, old_node.loc.id) : max(node.loc.id, old_node.loc.id)
            path.loc_reachability[id] = false
        end
    elseif path.entered_by_bottom
        # we put unreachable the locs from bottom
        for id = (old_node.loc.id - old_node.loc.col + 1) : node.loc.id
            path.loc_reachability[id] = false
        end
    else
        # we put unreachable the locs from the top
        for id = old_node.loc.id : (node.loc.id - node.loc.col)
            path.loc_reachability[id] = false
        end
        for id = (node.loc.id) : (node.loc.id + ins.cmax - node.loc.col)
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

function extend(path::Partial_path_sshape,new_node::Node)
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
    # update the entered by bottom
    if (path.node.id != 0) && (new_node.loc.aisle != path.node.loc.aisle)
        new_path.entered_by_bottom = !path.entered_by_bottom
    end
    # update visited locs and loc reachability
    push!(new_path.visited_nodes,new_node)
    update_loc_reachability(new_path,path.node)
    # update cost without cuts
    new_path.rcost += dist(path,new_node)
    # update all related to cuts
    update_cut_path(new_path)

    # now we test if a mandatory stop has become unreachable
    for id in eachindex(new_path.remaining_stops)
        if new_path.remaining_stops[id] != 0 && !new_path.loc_reachability[id] && id != new_path.node.loc.id
            return nothing
        end
    end

    return new_path
end