##########################################
# Notes
##########################################
# if the number of rows is odd, the part below the midpoint is larger than the part
# above the midpoint



function dist(path::Partial_path_largest,new_path::Partial_path_largest)
    ins = path.node.problem.ins
    pi_o = path.node.problem.pi_o
    sigma_o = path.node.problem.sigma_o
    # first the case of the ionode
    if path.node.id == 0
        return dist(path.node,new_path.node)
    elseif path.node.loc.aisle == new_path.node.loc.aisle
        return dist(path.node,new_path.node)
    elseif path.node.loc.aisle < new_path.node.loc.aisle
        # in this case we pass by the top
        di =  dist_locs_by_top(path.node.loc,new_path.node.loc,ins)
        di -= pi_o[new_path.node.loc.id]
        di -= sigma_o[new_path.node.loc.id]
        return di
    else
        di =  dist_locs_by_bottom(path.node.loc,new_path.node.loc,ins)
        di -= pi_o[new_path.node.loc.id]
        di -= sigma_o[new_path.node.loc.id]
        return di
    end
end

#########################################################################################
# Extension
##########################################################################################

function feasible_extension(path::Partial_path_largest, new_node::Node)
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
    elseif path.node.loc.id == 0
        return true
    elseif (new_node.loc.aisle == path.first_aisle) || (new_node.loc.aisle == path.last_aisle)
        # in this case return true
        # if the loc was unreachable the case is already managed
        return true
    elseif (path.node.loc.aisle == new_node.loc.aisle) && (path.node.loc.col < new_node.loc.col)
        # in this case we must not reduce the middle gap too much
        # we compute the new middle gap
        a = path.node.loc.aisle
        new_middle_gap = path.last_col_from_top[a] - new_node.loc.col
        # we compute the new minimum
        new_min = max(path.min_middle_gap[a] , new_node.loc.col - path.node.loc.col)
        if new_middle_gap < new_min
            return false
        else
            return true
        end
    elseif (path.node.loc.aisle > new_node.loc.aisle)
        # in this case we must not reduce the middle gap too much
        # we compute the new middle gap
        a = new_node.loc.aisle
        new_middle_gap = path.last_col_from_top[a] - new_node.loc.col
        # we compute the new minimum
        new_min = max(path.min_middle_gap[a] , new_node.loc.col)
        if new_middle_gap < new_min
            return false
        else
            return true
        end
    elseif (path.node.loc.aisle < new_node.loc.aisle) && path.last_aisle != 0
        return false
    else
        return true
    end
end

function update_loc_reachability(path::Partial_path_largest,old_node::Node)
    # This function updates the loc reachability vector
    # !!!!! IMPORTANT, the node of path should be the new node
    node = path.node
    ins = node.problem.ins
    
    # first we manage the case of the first extension (from the io node), so that the function doesn't bug
    if old_node.id == 0
        for id = 1 : node.loc.id
            path.loc_reachability[id] = false
        end
        # we also consume the first_aisle resource
        path.first_aisle = node.loc.aisle
    # now the case where the new loc is in the first aisle
    elseif node.loc.aisle == path.first_aisle
        for id = old_node.loc.id : node.loc.id
            path.loc_reachability[id] = false
        end
    # now the case where it is the last aisle
    elseif node.loc.aisle == path.last_aisle
        #we put unreachable everything above
        for id in (node.loc.id) : (node.loc.id - node.loc.col + ins.cmax)
            path.loc_reachability[id] = false
        end
    # now the case where they are in the same aisle
    elseif old_node.loc.aisle == node.loc.aisle
        for id in min((node.loc.id),(old_node.loc.id)) : max((node.loc.id),(old_node.loc.id))
            path.loc_reachability[id] = false
        end
        # we have to check if the last aisle resource is consumed
        if path.last_aisle == 0 
            # first we update the values
            a = node.loc.aisle
            path.last_col_from_top[a] = node.loc.col
            path.middle_gap[a] = node.loc.col
            if old_node.loc.col - node.loc.col > path.min_middle_gap[a]
                path.min_middle_gap[a] = old_node.loc.col - node.loc.col
            end
            # now we make the check
            if path.middle_gap[a] < path.min_middle_gap[a]
                path.last_aisle = node.loc.aisle
                for id = node.loc.id:length(ins.locList)
                    path.loc_reachability[id] = false
                end
            end
        else
            # in this case we are already below, the only thing to do is to update the gap
            # the check is done in the feasible extension function
            a = node.loc.aisle
            path.middle_gap[a] = path.last_col_from_top[a] - node.loc.col
            if node.loc.col - old_node.loc.col > path.min_middle_gap[a]
                path.min_middle_gap[a] = node.loc.col - old_node.loc.col
            end
        end
    # now the case where they are not in the same aisle
    # first the case coming from top, before the consumption of last aisle
    elseif node.loc.aisle > old_node.loc.aisle
        @assert path.last_aisle == 0
        # we manage the reachability in the new aisle
        for id in node.loc.id : (node.loc.id - node.loc.col + ins.cmax)
            path.loc_reachability[id] = false
        end
        # we have to manage the reachability in the old aisle if it's the first one
        if old_node.loc.aisle == path.first_aisle
            for id = old_node.loc.id : (old_node.loc.id - old_node.loc.col + ins.cmax)
                path.loc_reachability[id] = false
            end
        end
        # first we update the values
        a = node.loc.aisle
        path.last_col_from_top[a] = node.loc.col
        path.middle_gap[a] = node.loc.col
        if ins.cmax + 1 - node.loc.col > path.min_middle_gap[a]
            path.min_middle_gap[a] = ins.cmax + 1 - node.loc.col
        end
        # now we make the check
        if path.middle_gap[a] < path.min_middle_gap[a]
            path.last_aisle = node.loc.aisle
            for id = node.loc.id:length(ins.locList)
                path.loc_reachability[id] = false
            end
        end
    # the case of different aisles coming from bottom, i.e. after last aisle consumption
    else
        for id in (node.loc.id - node.loc.col + 1) : node.loc.id
            path.loc_reachability[id] = false
        end
        # we check if we have to consume the last aisle
        if path.last_aisle == 0
            path.last_aisle == old_node.loc.aisle
            # all the aisle after are unreachable
            for id = (node.loc.id - node.loc.col + ins.cmax) + 1 : length(ins.locList)
                path.loc_reachability[id] = false
            end
        else
            # in this case we must also put some loc unreachable
            for id = (node.loc.id - node.loc.col + ins.cmax) + 1 : (path.last_aisle*ins.cmax)
                path.loc_reachability[id] = false
            end
        end
        # we also update the gap
        a = node.loc.aisle
        path.middle_gap[a] = path.last_col_from_top[a] - node.loc.col
        if node.loc.col > path.min_middle_gap[a]
            path.min_middle_gap[a] = node.loc.col 
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

function extend(path::Partial_path_largest,new_node::Node)
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
    new_path.rcost += dist(path,new_path)
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