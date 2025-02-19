##############################################
# Distance matrix
###############################################

# function that compute the distance matrix
function create_distance_matrix_guo(ins::instance, nodelist::Vector{T},pi_o::Vector{Float64},sigma_o::Vector{Float64}) where {T<:AbstractNode}
    distance_matrix = fill(Inf,(length(nodelist)+1,length(nodelist)+1))
    # normal case
    for n1 in nodelist, n2 in nodelist
        if n1 == n2
            continue
        end
        # for return
        if n1.loc.id > n2.loc.id
            continue
        end
        if (n2.i == 1) && (n1.loc != n2.loc)
            distance_matrix[n1.id + 1, n2.id + 1] = distance_locs_guo(n1.loc,n2.loc,ins) - pi_o[n2.loc.id] - sigma_o[n2.loc.id]
        elseif (n1.loc == n2.loc) && (n2.i == n1.i + 1)
            distance_matrix[n1.id + 1, n2.id + 1] = distance_locs_guo(n1.loc,n2.loc,ins) - pi_o[n2.loc.id]
        end
    end
    # with iopoint
    for n in nodelist
        if n.i == 1
            distance_matrix[1,n.id+1] = distance_locs_guo(ioloc,n.loc,ins) - pi_o[n.loc.id] - sigma_o[n.loc.id]
        end
        distance_matrix[n.id+1,1] = distance_locs_guo(n.loc,ioloc,ins)
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

dist_locs_guo = 
[0	4.5	5.5	6.5	7.5	8.5	4.5	5.5	6.5	7.5	8.5	6.5	7.5	8.5	9.5	10.5	6.5	7.5	8.5	9.5	10.5;
4.5	0	1	2	3	4	7	8	9	10	11	9	10	11	12	13	9	10	11	12	13;
5.5	1	0	1	2	3	8	9	10	11	12	10	11	12	13	14	10	11	12	13	14;
6.5	2	1	0	1	2	9	10	11	12	13	11	12	13	14	15	11	12	13	14	15;
7.5	3	2	1	0	1	10	11	12	13	14	12	13	14	15	16	12	13	14	15	16;
8.5	4	3	2	1	0	11	12	13	14	15	13	14	15	16	17	13	14	15	16	17;
4.5	7	8	9	10	11	0	1	2	3	4	9	10	11	12	13	9	10	11	12	13;
5.5	8	9	10	11	12	1	0	1	2	3	10	11	12	13	14	10	11	12	13	14;
6.5	9	10	11	12	13	2	1	0	1	2	11	12	13	14	15	11	12	13	14	15;
7.5	10	11	12	13	14	3	2	1	0	1	12	13	14	15	16	12	13	14	15	16;
8.5	11	12	13	14	15	4	3	2	1	0	13	14	15	16	17	13	14	15	16	17;
6.5	9	10	11	12	13	9	10	11	12	13	0	1	2	3	4	7	8	9	10	11;
7.5	10	11	12	13	14	10	11	12	13	14	1	0	1	2	3	8	9	10	11	12;
8.5	11	12	13	14	15	11	12	13	14	15	2	1	0	1	2	9	10	11	12	13;
9.5	12	13	14	15	16	12	13	14	15	16	3	2	1	0	1	10	11	12	13	14;
10.5	13	14	15	16	17	13	14	15	16	17	4	3	2	1	0	11	12	13	14	15;
6.5	9	10	11	12	13	9	10	11	12	13	7	8	9	10	11	0	1	2	3	4;
7.5	10	11	12	13	14	10	11	12	13	14	8	9	10	11	12	1	0	1	2	3;
8.5	11	12	13	14	15	11	12	13	14	15	9	10	11	12	13	2	1	0	1	2;
9.5	12	13	14	15	16	12	13	14	15	16	10	11	12	13	14	3	2	1	0	1;
10.5	13	14	15	16	17	13	14	15	16	17	11	12	13	14	15	4	3	2	1	0
]


# The distance between 2 locs
function distance_locs_guo(a::loc,b::loc,ins::instance)
    return dist_locs_guo[a.id+1,b.id+1]
end

# the distance between 2 nodes
function dist(a::Node,b::Node)
    return a.problem.distance[a.id + 1,b.id + 1]
end

#########################################################################################
# Extension
##########################################################################################

#=function Base.isless(a::loc,b::loc)
    if a.id == 0 || b.id == 0
        return true
    elseif isodd(a.aisle) && isodd(b.aisle)
        return a.id < b.id
    elseif iseven(a.aisle) && iseven(b.aisle)
        return a.aisle > b.aisle || (a.aisle == b.aisle && a.col < b.col)
    elseif isodd(a.aisle) && iseven(b.aisle)
        return true
    else
        return false
    end
end=#
Base.isless(a::loc,b::loc) = (a.aisle < b.aisle) || ((a.aisle == b.aisle) && (a.col < b.col))

function feasible_extension(path::Partial_path_guo, new_node::Node)
    # returns true if the extension is feasible
    if new_node.id == 0
        return false
    elseif isinf(dist(path.node, new_node))
        return false
    elseif new_node.loc < path.node.loc
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

function update_loc_reachability(path::Partial_path_guo,old_node::Node)
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
    else
        for id = (old_node.loc.id) : node.loc.id
            path.loc_reachability[id] = false
        end
    end
    # RETURN
    # for return no need for loc reachability in general, the distance matrix is enough

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

function extend(path::Partial_path_guo,new_node::Node)
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