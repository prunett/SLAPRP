module pricing_return_module

using ..structures_decomp
using ..cost_functions
using ..cuts

#export pricing_main, pricing_improved



mutable struct partial_path
    rcost::Float64
    k::UInt8
    locs::Array{loc,1}
    og_cut_array::Array{Bool,1}
    n_current_loc::UInt8 # number of stops at the last location up to now (should be > n_sym)
end

# Defines a new isless method, to use the sorting functions on locations
function is_move_feasible(a::loc,b::loc,same_loc_ok::Bool)
# returns true if the move from a to b is possible
    if a.aisle < b.aisle
        return true
    elseif a.aisle == b.aisle
        if a.col < b.col
            return true
        elseif (a.col == b.col) && same_loc_ok
            return true
        end
    end
    return false
end
#=
function pricing_main(ins::instance,o::order,mu_o::Float64,pi_o::Array{Float64,1}, EPS::Float64)
# This is the main function for the label setting algorithm

    # setup
    locList = ins.locList
    #dist_o, successor_list = pricing_setup(ins,pi_o)
    K = length(o.sku)
    U = partial_path[partial_path(- mu_o ,0,loc[ioloc],1)]
    S = partial_path[]

    # Main loop where we extend a partial path at each iteration
    while !isempty(U)
        # Pick a path from U (and remove it)
        p = popfirst!(U)
        #for lid in (length(p.locs) > 1 ? successor_list[p.locs[end].id] : 1:length(locList))
        for lid = 1:length(locList)
            # Check if it is possible to go to lid, OR if it's the same loc and still allowed
            same_loc_ok = (p.n_current_loc < ins.n_sym)
            if !is_move_feasible(p.locs[end],locList[lid],same_loc_ok)
                continue
            end
            # Extend the path p to p' and create label
            #p_new = partial_path(p.rcost + (length(p.locs) > 1 ? dist_o[p.locs[end].id,lid] : dist(ioloc, locList[lid], ins) - pi_o[lid]) , p.k + 1, vcat(p.locs,[locList[lid]]))
            p_new = partial_path(p.rcost + dist(p.locs[end],locList[lid],ins) - pi_o[lid] , p.k + 1, vcat(p.locs,[locList[lid]]), (p.locs[end].id == lid ? p.n_current_loc + 1 : 1))
            # Check cuts

            #test where to add it + check dominance
            if p_new.k < K
                keep = basic_dominance!(U,p_new)
                if keep
                    push!(U,p_new)
                end
            else
                p_new.rcost += dist(p_new.locs[end],ioloc,ins)
                if p_new.rcost < - EPS && basic_dominance!(S,p_new)
                    push!(S,p_new)
                end
            end
            # Apply dominance rules
        end
    end

    # maybe some post treatment on S
    routes_generated = route_transform(ins,o,S)
    rcost_generated = [p.rcost for p in S]
    for i=1:length(S)
        if abs(rcost_generated[i] - routes_generated[i].cost + mu_o + sum(pi_o[l.id] for l in routes_generated[i].loc)) > 0.001
            println("problem !!")
        end
    end
    return routes_generated, rcost_generated
end=#

is_less_loc(a::loc,b::loc) = (a.aisle < b.aisle) || ((a.aisle == b.aisle) && (a.col < b.col))

function is_og_cut_reachable(current_loc::loc, cut::cuts.OG_cut)
# return true if the cut is still reachable starting from this loc
    for l in cut.loc_array
        if is_less_loc(current_loc, l)
            return true
        end
    end
    return false
end



function pricing_improved(ins::instance,o::order,mu_o::Float64,pi_o::Array{Float64,1},sigma_o::Array{Float64,1},cut_pool::cuts.CutCollection,EPS::Float64)
# This function is the pricing problem with an improved data structure.

    # setup
    K = length(o.sku)
    locList = ins.locList
    n_sym = ins.n_sym
    # OG_cuts available
    OG_available = cuts.OG_cut[]
    for cut in cut_pool["OG_cut"][o]
        if cut.active
            push!(OG_available, cut)
        end
    end
    #array with the labels
    #U = fill(partial_path[], (length(locList), n_sym, K))
    #U = fill(partial_path[partial_path(Inf,0,loc[],Bool[],0)], (length(locList), n_sym, K))
    U = Array{Array{partial_path,1},3}(undef, (length(locList),n_sym,K))
    for l=1:length(locList), s=1:n_sym, k=1:K
        U[l,s,k] = partial_path[]
    end

    # Main loop
    done = false
    l = 0
    s = 1
    k = 0
    # the current index of the label in the array in U
    current_label_index = 1
    while !done
        #println("start, l = $l, s = $s, k = $k, current = $current_label_index")
        # Get the current partial path to extend
        if l == 0
            p = partial_path(-mu_o,0,loc[ioloc],[true for c in OG_available],1)
            l = 1
        elseif isempty(U[l,s,k])
            # continue if the label array is empty (mostly due to s and k being too low)
            done,l,s,k, current_label_index = update_current_label!(ins,K,l,s,k,U,current_label_index)
            continue
        else
            p = U[l,s,k][current_label_index]
            #@assert(p.k == k)
            #@assert(p.n_current_loc == s)
            #@assert(p.locs[end].id == l)
        end


        # Extend the label
        for lid = 1:length(locList)
            # Check move feasability
            same_loc_ok = (p.n_current_loc < n_sym)
            if !is_move_feasible(p.locs[end],locList[lid],same_loc_ok)
                continue
            end

            # Extend the path
            new_cost = p.rcost + dist(p.locs[end],locList[lid],ins) - pi_o[lid]
            new_s = (p.locs[end].id == lid ? p.n_current_loc + 1 : 1)
            # add the sigma cost if relevant
            new_cost -= (new_s == 1 ? sigma_o[lid] : 0)

            # add the OG_cuts cost + manage available cuts for extended path
            new_available = copy(p.og_cut_array)
            for i in eachindex(new_available)
                if !new_available[i]
                    # cut non available anymore
                    continue
                end
                # test if we add the cut
                if locList[lid] in OG_available[i].loc_array
                    # add dual cost
                    new_cost -= OG_available[i].dual_value
                    # remove it from available
                    new_available[i] = false
                elseif !is_og_cut_reachable(locList[lid], OG_available[i])
                    # if the cut is not reachable anymore,
                    new_available[i] = false
                end
            end

            # Check if the dominance
            dominated = false
            ind = 1
            while ind <= length(U[lid,new_s,p.k+1])
                # if the new cut is dominated => break
                if dominate_OG(U[lid,new_s,p.k+1][ind], new_cost,new_available,OG_available)
                    #println("dominated !!")
                    dominated = true
                    break
                end
                # if the old cut is dominated => delete it
                if dominate_OG(new_cost,new_available, U[lid,new_s,p.k+1][ind],OG_available)
                    #println("dominate !!")
                    deleteat!(U[lid,new_s,p.k+1], ind)
                    ind -= 1
                end
                ind += 1
            end
            # if the new cut is not dominated => add it to the list
            if !dominated
                new_locs = copy(p.locs)
                push!(new_locs, locList[lid])
                #new_locs = vcat(p.locs, locList[lid])
                #@assert(lid == locList[lid].id)
                push!(U[lid,new_s,p.k+1], partial_path(new_cost, p.k + 1, new_locs, new_available, new_s))
            end

        end

        # get to next label to look at
        done,l,s,k,current_label_index = update_current_label!(ins,K,l,s,k,U,current_label_index)
    end

    # Retrieve the non dominated routes if their reduced cost is < 0
    return retrieve_routes(ins,U,K,o,EPS)

end

function dominate_OG(p::partial_path,new_cost::Float64, new_available::Array{Bool,1},OG_available::Vector{cuts.OG_cut})
# checks if p dominates the new path
    cost_test = new_cost
    for i in eachindex(new_available)
        if p.og_cut_array[i] < new_available[i]
            cost_test -= OG_available[i].dual_value
        end
    end
    return p.rcost <= cost_test
end

function dominate_OG(new_cost::Float64, new_available::Array{Bool,1},p::partial_path,OG_available::Vector{cuts.OG_cut})
    # checks if new_path dominates p
        cost_test = p.rcost
        for i in eachindex(new_available)
            if p.og_cut_array[i] > new_available[i]
                cost_test -= OG_available[i].dual_value
            end
        end
        return new_cost <= cost_test
    end
#=
function dominate_OG(new_cost::Float64,new_available::Array{Bool,1},p::partial_path)
# This is the opposite function, it checks if the new path dominates p
    return new_cost <= p.rcost && prod(p.og_cut_array .<= new_available)
end=#

function retrieve_routes(ins::instance,U,K,o,EPS)
# This function returns the routes with negative reduced cost
    retrieved_routes = route[]
    retrieved_rcost = Float64[]
    locList = ins.locList
    for l=1:length(locList)
        for s=1:ins.n_sym
            for i in eachindex(U[l,s,K])
            # test if reduced cost negative
                if U[l,s,K][i].rcost + dist(locList[l],ioloc,ins) < - EPS #|| true
                    push!(retrieved_routes, route_cost(ins,U[l,s,K][i].locs[2:end],o))
                    push!(retrieved_rcost, U[l,s,K][i].rcost + dist(locList[l],ioloc,ins))
                end
            end
        end
    end
    return retrieved_routes, retrieved_rcost
end

function update_current_label!(ins::instance,K,l,s,k,U,current_label_index)
# This function updates the values of l,s and k
    # special case if there is only one stop possible
    if K == 1
        return true,l,s,k,current_label_index
    end
    # if k = 0 then k+1
    if k == 0
        return false, l, s, k+1, 1
    end
    # first check if we are at the last label of the array in U
    if current_label_index < length(U[l,s,k])
        return false, l, s, k, current_label_index+1
    end
    if (k + 1 < K) #&& (k + 1 <= (l-1)*ins.n_sym + s)
        return false, l, s, k+1, 1
    end
    if s + 1 <= ins.n_sym
        return false, l, s+1, 1, 1
    end
    if l + 1 <= length(ins.locList)
        return false, l+1, 1, 1, 1
    end
    # The current label was the last one to check
    return true,l,s,k,current_label_index
end

function pricing_setup(ins::instance,pi_o::Array{Float64,1})
# This is a setup function that prepares the input for the pricing algorithm
    locList = ins.locList

    # Create a distance matrix
    dist_o = zeros(length(locList),length(locList))
    successor_list = Array{Array{Int,1},1}()
    for l1 = 1:length(locList)
        push!(successor_list, Int[])
        for l2 = 1:length(locList)
            if locList[l1] < locList[l2]
                dist_o[l1,l2] = dist(locList[l1],locList[l2],ins) - pi_o[l2]
                push!(successor_list[l1], l2)
            else
                dist_o[l1,l2] = Inf
            end
        end
    end

    return dist_o, successor_list
end

function route_transform(ins,o,S)
# This function transforms the resulting partial paths into routes
    # Array of routes
    routes_generated = route[]
    # Generation in the same order as in S
    for p in S
        push!(routes_generated, route_cost(ins,p.locs[2:end],o))
    end

    return routes_generated
end


function basic_dominance!(U,p_new)
# This function applies the basic dominance rules
# It can be improved in terms of complexity by a better data structure

    for pind = 1:length(U)
        # continue if the rules do not apply
        if (U[pind].locs[end].id != p_new.locs[end].id) || (U[pind].k != p_new.k) || (U[pind].n_current_loc != p_new.n_current_loc)
            continue
        end
        # check rcost
        if U[pind].rcost <= p_new.rcost
            return false
        else
            deleteat!(U,pind)
            return true
        end
    end

    return true

end



end # end module
