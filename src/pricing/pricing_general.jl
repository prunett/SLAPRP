module pricing_general

using ..structures_decomp
using ..cost_functions
using ..cuts

#################################################################
# Structures
#################################################################



include("pricing_utilities.jl")
include("pricing_general_optimal.jl")
include("pricing_general_midpoint.jl")
include("pricing_general_return.jl")
include("pricing_general_sshape.jl")
include("pricing_general_largest.jl")
include("pricing_guo_return.jl")




##############################################################
# Main function
##############################################################

function pricing(ins::instance,o::order,mu_o::Float64,pi_o::Array{Float64,1},sigma_o::Array{Float64,1},cut_pool::cuts.CutCollection,fixed_branches::Dict{sku,loc},EPS::Float64,policy::String)
    # This is the main function for pricing 
    # Initialization
    pricing = create_pricing_problem(ins,o,mu_o,pi_o,sigma_o,cut_pool,EPS,policy,fixed_branches)
    # Main loop for label extension
    while !isfinished(pricing)
        # pop the first label of the list
        path = pop_next_label!(pricing)
        # Check if the label is dominated
        if !isdominated!(path,pricing) 
            # Check the length of the label, to put it in closed_labels
            if path.k >= length(o.sku)
                push!(pricing.closed_labels[path.node.id], path)
            end
            # Now extend the label if it is still extandable
            if !done_extending(path)
                for new_node in pricing.node_list
                    # Only if it is a valid successor
                    # UNIT
                    #=if o.id == 7 && new_node.loc.id == 4 && length(path.visited_nodes) == 1# && path.visited_nodes[2].loc.id == 4
                        println("locs of path = $([v.loc.id for v in path.visited_nodes])")
                        println("loc reachability = $(path.loc_reachability)")
                        println("feasible extension = $(feasible_extension(path,new_node))")
                        println("first aisle = $(path.first_aisle), last = $(path.last_aisle)")
                        @show path.last_col_from_top
                        @show path.middle_gap
                        @show path.min_middle_gap
                    end=#
                    if feasible_extension(path,new_node)
                        new_path = extend(path,new_node)
                        # insert the label
                        if !isnothing(new_path)
                            push!(pricing.labels[new_path.node.id], new_path)
                        end
                    end
                end
            end
        end
    end
       
    
    routes, rcosts = create_routes(pricing)

    return routes, rcosts
end

########################################################################
# Basic utility functions
########################################################################

function isfinished(pricing::AbstractPricing) 
    # This function return true if we are done extending labels (no more open label to extend)
    return pricing.current_k > length(pricing.order.sku) 
end
    
function pop_next_label!(pricing::AbstractPricing) 
    # This function returns the next label to explore, and remove it from the list
    # It also update (potentially) current label list to explore
    # first we code the case of the current node is the ionode
    if pricing.current_node == pricing.ionode
        # pass to the first k at the first node
        pricing.current_k = 1
        pricing.current_node = pricing.node_list[1]
        # return an empty initiated label
        return create_empty_label(pricing.ins,pricing,pricing.mu_o,pricing.ionode)
    end
    # termination condition on k
    if pricing.current_k > length(pricing.order.sku)
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

function done_extending(path::AbstractPath)
    return path.k >= length(path.node.problem.order.sku)
end

function collect_mandatory_stops(o::order,ins::instance,fixed_branches::Dict{sku,loc})
    # This function collects the mandatory stops
    mandatory_stops = zeros(Int,length(ins.locList))
    # we look for each sku of the order
    for s in o.sku
        if s in keys(ins.fixed_positions)
            mandatory_stops[ins.fixed_positions[s].id] += 1
        elseif s in keys(fixed_branches)
            mandatory_stops[fixed_branches[s].id] += 1
        end
    end

    return mandatory_stops
end

function compute_limited_loc(ins::instance,o::order,fixed_branches::Dict{sku,loc})
    # This function computes the limited loc (i.e. the ones with only one stop allowed)
    # it returns an array with boolean for each loc, where true = max one stop
    initial_loc_availability = fill(2,length(ins.locList))
    for s in keys(fixed_branches)
        if !(s in o.sku)
            initial_loc_availability[fixed_branches[s].id] -= 1
        end
    end
    for s in keys(ins.fixed_positions)
        if !(s in o.sku)
            initial_loc_availability[ins.fixed_positions[s].id] -= 1
        end
    end
    
    return initial_loc_availability
end


##############################################################
# Label extension
##############################################################

function feasible_extension(path::AbstractPath,new_node::Node) end
function update_loc_reachability(path::AbstractPath,old_node::Node) end
function extend(path::AbstractPath,new_node::Node) end


function update_cut_path(path::AbstractPath)
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

#################################################################
## Label dominance
#################################################################

function isdominated!(path::AbstractPath, pricing::AbstractPricing)
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
isdominated!(path::Missing, pricing::AbstractPricing) = true

function dominance_test(patha::AbstractPath, pathb::AbstractPath)
    # returns true if path a dominates path b
    @assert patha.node == pathb.node
    # first test if we have the same k
    if patha.k != pathb.k
        return false
    end
    # Then test the consumption of the mandatory stops
    for i in eachindex(patha.remaining_stops)
        if patha.remaining_stops[i] != pathb.remaining_stops[i]
            return false
        end
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

############################################################
# Route creation
############################################################

function create_routes(pricing::AbstractPricing)
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
                # UNIT
                #=if abs(rcost_computation_test(new_route,pricing) - label.rcost - dist(label.visited_nodes[end],pricing.ionode)) > pricing.EPS
                    println("rcost_computation_test = $(rcost_computation_test(new_route,pricing)), label = $(label.rcost + dist(label.visited_nodes[end],pricing.ionode))")
                    println("pi_o = $(pricing.pi_o)")
                    println("sigma = $(pricing.sigma_o)")
                    println("mu = $(pricing.mu_o)")
                    println("label locs = $([n.loc.id for n in label.visited_nodes])")
                    println("locs = $([l.id for l in new_route.loc]) route cost = $(new_route.cost)")
                    #error("crong rcost")
                end=#
                if insert_route
                    push!(route_list, new_route)
                    push!(rcost_list, rcost_computation_test(route_list[end],pricing))
                end
            end
        end
    end

    return route_list, rcost_list
end

function transform_label_into_route(path::AbstractPath)
    # This function transform the label in a route
    ins = path.node.problem.ins
    return route_cost(ins,[n.loc for n in path.visited_nodes[2:end]], path.node.problem.order)
end
    
function transform_label_into_route(path_vector::Vector{T}) where {T<:AbstractPath}
    resulting_routes = structures_decomp.route[]
    for p in path_vector
        route = transform_label_into_route(p)
        push!(resulting_routes,route)
    end
    return resulting_routes
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

function rcost_computation_test(r::route, pricing_problem::AbstractPricing)
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

end # end fo module