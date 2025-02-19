##########################################
# Notes
##########################################
# if the number of rows is odd, the part below the midpoint is larger than the part
# above the midpoint



function dist(path::Partial_path_midpoint,new_path::Partial_path_midpoint)
    ins = path.node.problem.ins
    mp = ceil(path.node.problem.ins.cmax / 2)
    pi_o = path.node.problem.pi_o
    sigma_o = path.node.problem.sigma_o
    # first the case of the ionode
    if path.node.id == 0
        return dist(path.node,new_path.node)
    elseif path.node.loc.aisle == new_path.node.loc.aisle
        return dist(path.node,new_path.node)
    elseif path.node.loc.aisle == new_path.first_aisle
        # in this case we pass by the top
        di =  dist_locs_by_top(path.node.loc,new_path.node.loc,ins)
        di -= pi_o[new_path.node.loc.id]
        di -= sigma_o[new_path.node.loc.id]
        return di
    elseif (path.node.loc.col <= mp && new_path.node.loc.col <= mp) || (path.node.loc.col >= mp + 1 && new_path.node.loc.col >= mp + 1)
        # this is the case where they are both above mp or both below
        return dist(path.node,new_path.node)
    elseif new_path.node.loc.aisle == new_path.last_aisle
        # in this case the last aisle is the one from new_node
        di =  dist_locs_by_top(path.node.loc,new_path.node.loc,ins)
        di -= pi_o[new_path.node.loc.id]
        di -= sigma_o[new_path.node.loc.id]
        return di
    elseif path.node.loc.aisle == new_path.last_aisle
        di =  dist_locs_by_bottom(path.node.loc,new_path.node.loc,ins)
        di -= pi_o[new_path.node.loc.id]
        di -= sigma_o[new_path.node.loc.id]
        return di
    else
        println("last aisle = $(new_path.last_aisle)")
        println("new_path = $([l.loc.id for l in new_path.visited_nodes])")
        error("This case should not appear")
    end
end

#########################################################################################
# Extension
##########################################################################################

function feasible_extension(path::Partial_path_midpoint, new_node::Node)
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

function update_loc_reachability(path::Partial_path_midpoint,old_node::Node)
    # This function updates the loc reachability vector
    # !!!!! IMPORTANT, the node of path should be the new node
    node = path.node
    ins = node.problem.ins
    mp = Int(ceil(ins.cmax / 2)) # midpoint
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
    # now the case where the loc is in the upper part but another aisle than the first one
    elseif node.loc.col >= mp + 1
        # first the case where they are in the same aisle
        if node.loc.aisle == old_node.loc.aisle
            for id in (node.loc.id) : (old_node.loc.id)
                path.loc_reachability[id] = false
            end
        # now the case where it is not the same aisle
        else
            # we put unreachable the locs above mp from the aisle of old node to the aisle just before 
            # the new node
            for a in old_node.loc.aisle : node.loc.aisle - 1 
                for id in ((a-1)*ins.cmax + mp + 1) : ((a-1)*ins.cmax + ins.cmax)
                    path.loc_reachability[id] = false
                end
            end
            # then above the old node (in case of it is in the first aisle)
            for id in (old_node.loc.id) : (old_node.loc.id - old_node.loc.col + ins.cmax)
                path.loc_reachability[id] = false
            end
            # then the locs of the new aisle
            for id in (node.loc.id) : (node.loc.id - node.loc.col + ins.cmax)
                path.loc_reachability[id] = false
            end
        end
    # now the case where it is below the mp and the resource has not be consumed yet
    elseif path.last_aisle == 0
        # we consume the resource 
        path.last_aisle = max(old_node.loc.aisle,node.loc.aisle)
        # we put unreachable the first aisle
        for id in ((path.first_aisle - 1)*ins.cmax + 1) : (path.first_aisle * ins.cmax)
            path.loc_reachability[id] = false
        end
        # then we put unreachable the locs after the aisle of new node (may not be the last)
        for id in (node.loc.aisle*ins.cmax + 1) : length(ins.locList)
            path.loc_reachability[id] = false
        end
        # then we put all the locs above mp
        for id in (path.first_aisle*ins.cmax + 1) : ((path.last_aisle - 1)*ins.cmax)
            if ins.locList[id].col >= mp + 1
                path.loc_reachability[id] = false
            end
        end
        # now there are two case depending of if we visit the last aisle from above of bottom
        # first the case from above
        if node.loc.aisle == path.last_aisle
            # then we put the locs from the top of the last aisle to the new node
            for id in (node.loc.id) : (node.loc.id - node.loc.col + ins.cmax)
                path.loc_reachability[id] = false
            end
        else
            # in this case it is visited from below
            for id in (node.loc.id - node.loc.col + 1) : node.loc.id
                path.loc_reachability[id] = false
            end
        end
    # now the case where it is the last aisle
    elseif node.loc.aisle == path.last_aisle
        #we put unreachable everything above
        for id in (node.loc.id) : (node.loc.id - node.loc.col + ins.cmax)
            path.loc_reachability[id] = false
        end
    # now the last case, a loc below mp, not the first, not the last
    else
        #first the case where they are in the same aisle
        if node.loc.aisle == old_node.loc.aisle
            for id in old_node.loc.id : node.loc.id
                path.loc_reachability[id] = false
            end
        else
            # first the locs below mp in the aisles strictly after the one of new node
            # are put to unreachable
            for a in (node.loc.aisle + 1) : old_node.loc.aisle
                for id in ((a-1)*ins.cmax + 1) : ((a-1)*ins.cmax + mp)
                    path.loc_reachability[id] = false
                end
            end
            # then the locs below the new node in this aisle are unreachable
            for id in (node.loc.id - node.loc.col + 1) : node.loc.id
                path.loc_reachability[id] = false
            end
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

function extend(path::Partial_path_midpoint,new_node::Node)
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