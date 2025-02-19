

abstract type AbstractNode end
abstract type AbstractPath end
abstract type AbstractPricing end

struct Node <:AbstractNode
    id::Int
    loc::loc
    i::Int
    problem::AbstractPricing

    function Node(id::Int,loc::loc,i::Int,pricing_problem::AbstractPricing)
        return new(id,loc,i,pricing_problem)
    end
end

Base.:(==)(a::Node,b::Node) = (a.id == b.id) 
Base.isequal(a::Node,b::Node) = isequal(a.id, b.id)
Base.hash(a::Node, h::UInt) = hash(a.id, hash(:Node,h))


###################################################
# Data structures optimal
###################################################

mutable struct Pricing_problem_optimal{T<:AbstractNode, P<:AbstractPath} <:AbstractPricing
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
    mandatory_stops::Vector{Int}
    fixed_branches::Dict{sku,loc}
    initial_loc_availability::Vector{UInt8}

    function Pricing_problem_optimal{T,P}(ins::instance,o::structures_decomp.order,mu_o::Float64,pi_o::Vector{Float64},sigma_o::Vector{Float64},cut_pool::cuts.CutCollection, number_blocks_layout::Int, EPS::Float64,fixed_branches::Dict{sku,loc}) where {T <: AbstractNode, P <: AbstractPath}
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
        pricing_problem.distance = create_distance_matrix(ins,pricing_problem.node_list,pi_o,sigma_o,ins.number_blocks_layout)
        #pricing_problem.open_labels = P[Partial_path{T}(ins,pricing_problem,mu_o,pricing_problem.ionode)]
        pricing_problem.labels = Vector{P}[[] for n1 in pricing_problem.node_list]
        pricing_problem.closed_labels = Vector{P}[[] for n1 in pricing_problem.node_list]
        pricing_problem.current_k = 0
        pricing_problem.current_node = T(0,ioloc,0,pricing_problem)
        pricing_problem.mu_o = mu_o
        pricing_problem.pi_o = pi_o
        pricing_problem.sigma_o = sigma_o
        pricing_problem.EPS = EPS

        pricing_problem.mandatory_stops = collect_mandatory_stops(o,ins,fixed_branches)
        pricing_problem.fixed_branches = fixed_branches
        pricing_problem.initial_loc_availability = compute_limited_loc(ins,o,fixed_branches)
        
        return pricing_problem
    end
end

mutable struct Partial_path_optimal{T<:AbstractNode} <:AbstractPath
    node::T
    rcost::Float64
    k::UInt8
    remaining_free_stops::UInt8
    visited_nodes::Vector{T}
    collected_cuts::Vector{Bool}
    loc_reachability::Vector{Bool}
    cut_reachability::Vector{Bool}
    remaining_stops::Vector{UInt8}
end



function create_empty_label(ins::instance,pricing::Pricing_problem_optimal, mu_o::Float64,ionode::AbstractNode) 
    node = ionode
    rcost = - mu_o
    k = 0
    remaining_free = length(pricing.order.sku) - sum(pricing.mandatory_stops)
    visited_nodes = [node]
    collected_cuts = Bool[false for c in pricing.available_cuts]
    loc_reachability = Bool[true for l in ins.locList]
    cut_reachability = Bool[true for c in pricing.available_cuts]
    remaining_stops = copy(pricing.mandatory_stops)

    # we check the case where there is no free stop, in this case we update the reachability
    if remaining_free == 0
        for i in eachindex(loc_reachability)
            if remaining_stops[i] == 0
                loc_reachability[i] = false
            end
        end
    end

    # Then we put also to unreachable the locs where there is no initial availability
    for i in eachindex(pricing.initial_loc_availability)
        if pricing.initial_loc_availability[i] == 0
            loc_reachability[i] = false
        end
    end

    return Partial_path_optimal{typeof(ionode)}(node,rcost,k,remaining_free,visited_nodes,collected_cuts,loc_reachability,cut_reachability,remaining_stops)
end

# copy methods
Base.copy(path::Partial_path_optimal{Node}) = Partial_path_optimal{Node}(path.node,path.rcost,path.k,
    path.remaining_free_stops,copy(path.visited_nodes),copy(path.collected_cuts),copy(path.loc_reachability),
    copy(path.cut_reachability),copy(path.remaining_stops))


###################################################
# Data structures return
###################################################

mutable struct Pricing_problem_return{T<:AbstractNode, P<:AbstractPath} <:AbstractPricing
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
    mandatory_stops::Vector{Int}
    fixed_branches::Dict{sku,loc}
    initial_loc_availability::Vector{UInt8}

    function Pricing_problem_return{T,P}(ins::instance,o::structures_decomp.order,mu_o::Float64,pi_o::Vector{Float64},sigma_o::Vector{Float64},cut_pool::cuts.CutCollection, number_blocks_layout::Int, EPS::Float64,fixed_branches::Dict{sku,loc}) where {T <: AbstractNode, P <: AbstractPath}
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
        pricing_problem.distance = create_distance_matrix_return(ins,pricing_problem.node_list,pi_o,sigma_o)
        #pricing_problem.open_labels = P[Partial_path{T}(ins,pricing_problem,mu_o,pricing_problem.ionode)]
        pricing_problem.labels = Vector{P}[[] for n1 in pricing_problem.node_list]
        pricing_problem.closed_labels = Vector{P}[[] for n1 in pricing_problem.node_list]
        pricing_problem.current_k = 0
        pricing_problem.current_node = T(0,ioloc,0,pricing_problem)
        pricing_problem.mu_o = mu_o
        pricing_problem.pi_o = pi_o
        pricing_problem.sigma_o = sigma_o
        pricing_problem.EPS = EPS

        pricing_problem.mandatory_stops = collect_mandatory_stops(o,ins,fixed_branches)
        pricing_problem.fixed_branches = fixed_branches
        pricing_problem.initial_loc_availability = compute_limited_loc(ins,o,fixed_branches)
        
        return pricing_problem
    end
end

mutable struct Partial_path_return{T<:AbstractNode} <:AbstractPath
    node::T
    rcost::Float64
    k::UInt8
    remaining_free_stops::UInt8
    visited_nodes::Vector{T}
    collected_cuts::Vector{Bool}
    loc_reachability::Vector{Bool}
    cut_reachability::Vector{Bool}
    remaining_stops::Vector{UInt8}
end



function create_empty_label(ins::instance,pricing::Pricing_problem_return, mu_o::Float64,ionode::AbstractNode) 
    node = ionode
    rcost = - mu_o
    k = 0
    remaining_free = length(pricing.order.sku) - sum(pricing.mandatory_stops)
    visited_nodes = [node]
    collected_cuts = Bool[false for c in pricing.available_cuts]
    loc_reachability = Bool[true for l in ins.locList]
    cut_reachability = Bool[true for c in pricing.available_cuts]
    remaining_stops = copy(pricing.mandatory_stops)

    # we check the case where there is no free stop, in this case we update the reachability
    if remaining_free == 0
        for i in eachindex(loc_reachability)
            if remaining_stops[i] == 0
                loc_reachability[i] = false
            end
        end
    end

    # Then we put also to unreachable the locs where there is no initial availability
    for i in eachindex(pricing.initial_loc_availability)
        if pricing.initial_loc_availability[i] == 0
            loc_reachability[i] = false
        end
    end

    return Partial_path_return{typeof(ionode)}(node,rcost,k,remaining_free,visited_nodes,collected_cuts,loc_reachability,cut_reachability,remaining_stops)
end

# copy methods
Base.copy(path::Partial_path_return{Node}) = Partial_path_return{Node}(path.node,path.rcost,path.k,
    path.remaining_free_stops,copy(path.visited_nodes),copy(path.collected_cuts),copy(path.loc_reachability),
    copy(path.cut_reachability),copy(path.remaining_stops))


    ###################################################
# Data structures return guo
###################################################

mutable struct Pricing_problem_guo{T<:AbstractNode, P<:AbstractPath} <:AbstractPricing
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
    mandatory_stops::Vector{Int}
    fixed_branches::Dict{sku,loc}
    initial_loc_availability::Vector{UInt8}

    function Pricing_problem_guo{T,P}(ins::instance,o::structures_decomp.order,mu_o::Float64,pi_o::Vector{Float64},sigma_o::Vector{Float64},cut_pool::cuts.CutCollection, number_blocks_layout::Int, EPS::Float64,fixed_branches::Dict{sku,loc}) where {T <: AbstractNode, P <: AbstractPath}
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
        pricing_problem.distance = create_distance_matrix_guo(ins,pricing_problem.node_list,pi_o,sigma_o)
        #pricing_problem.open_labels = P[Partial_path{T}(ins,pricing_problem,mu_o,pricing_problem.ionode)]
        pricing_problem.labels = Vector{P}[[] for n1 in pricing_problem.node_list]
        pricing_problem.closed_labels = Vector{P}[[] for n1 in pricing_problem.node_list]
        pricing_problem.current_k = 0
        pricing_problem.current_node = T(0,ioloc,0,pricing_problem)
        pricing_problem.mu_o = mu_o
        pricing_problem.pi_o = pi_o
        pricing_problem.sigma_o = sigma_o
        pricing_problem.EPS = EPS

        pricing_problem.mandatory_stops = collect_mandatory_stops(o,ins,fixed_branches)
        pricing_problem.fixed_branches = fixed_branches
        pricing_problem.initial_loc_availability = compute_limited_loc(ins,o,fixed_branches)
        
        return pricing_problem
    end
end

mutable struct Partial_path_guo{T<:AbstractNode} <:AbstractPath
    node::T
    rcost::Float64
    k::UInt8
    remaining_free_stops::UInt8
    visited_nodes::Vector{T}
    collected_cuts::Vector{Bool}
    loc_reachability::Vector{Bool}
    cut_reachability::Vector{Bool}
    remaining_stops::Vector{UInt8}
end



function create_empty_label(ins::instance,pricing::Pricing_problem_guo, mu_o::Float64,ionode::AbstractNode) 
    node = ionode
    rcost = - mu_o
    k = 0
    remaining_free = length(pricing.order.sku) - sum(pricing.mandatory_stops)
    visited_nodes = [node]
    collected_cuts = Bool[false for c in pricing.available_cuts]
    loc_reachability = Bool[true for l in ins.locList]
    cut_reachability = Bool[true for c in pricing.available_cuts]
    remaining_stops = copy(pricing.mandatory_stops)

    # we check the case where there is no free stop, in this case we update the reachability
    if remaining_free == 0
        for i in eachindex(loc_reachability)
            if remaining_stops[i] == 0
                loc_reachability[i] = false
            end
        end
    end

    # Then we put also to unreachable the locs where there is no initial availability
    for i in eachindex(pricing.initial_loc_availability)
        if pricing.initial_loc_availability[i] == 0
            loc_reachability[i] = false
        end
    end

    return Partial_path_guo{typeof(ionode)}(node,rcost,k,remaining_free,visited_nodes,collected_cuts,loc_reachability,cut_reachability,remaining_stops)
end

# copy methods
Base.copy(path::Partial_path_guo{Node}) = Partial_path_guo{Node}(path.node,path.rcost,path.k,
    path.remaining_free_stops,copy(path.visited_nodes),copy(path.collected_cuts),copy(path.loc_reachability),
    copy(path.cut_reachability),copy(path.remaining_stops))


###################################################
# Data structures sshape
###################################################

mutable struct Pricing_problem_sshape{T<:AbstractNode, P<:AbstractPath} <:AbstractPricing
    ins::instance
    order::order
    node_list::Vector{T}
    ionode::T
    available_cuts::Vector{cuts.OG_cut}
    distance::Array{Float64,2}
    distance_from_bottom::Array{Float64,2}
    distance_from_top::Array{Float64,2}
    labels::Vector{Vector{P}} #Vector of labels, classified by node they are attached to
    closed_labels::Vector{Vector{P}} # Vector of labels with the good length(s)
    current_k::Int # we keep these info in memory to know we we are in the search
    # in order to know which label should be treated next, and when to stop
    current_node::T # same, which node are we exploring
    mu_o::Float64
    pi_o::Vector{Float64}
    sigma_o::Vector{Float64}
    EPS::Float64
    mandatory_stops::Vector{Int}
    fixed_branches::Dict{sku,loc}
    initial_loc_availability::Vector{UInt8}

    function Pricing_problem_sshape{T,P}(ins::instance,o::structures_decomp.order,mu_o::Float64,pi_o::Vector{Float64},sigma_o::Vector{Float64},cut_pool::cuts.CutCollection, number_blocks_layout::Int, EPS::Float64,fixed_branches::Dict{sku,loc}) where {T <: AbstractNode, P <: AbstractPath}
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
        pricing_problem.distance = create_distance_matrix(ins,pricing_problem.node_list,pi_o,sigma_o,1)
        pricing_problem.distance_from_bottom, pricing_problem.distance_from_top = create_distance_matrix_sshape(ins,pricing_problem.node_list,pi_o,sigma_o)
        #pricing_problem.open_labels = P[Partial_path{T}(ins,pricing_problem,mu_o,pricing_problem.ionode)]
        pricing_problem.labels = Vector{P}[[] for n1 in pricing_problem.node_list]
        pricing_problem.closed_labels = Vector{P}[[] for n1 in pricing_problem.node_list]
        pricing_problem.current_k = 0
        pricing_problem.current_node = T(0,ioloc,0,pricing_problem)
        pricing_problem.mu_o = mu_o
        pricing_problem.pi_o = pi_o
        pricing_problem.sigma_o = sigma_o
        pricing_problem.EPS = EPS

        pricing_problem.mandatory_stops = collect_mandatory_stops(o,ins,fixed_branches)
        pricing_problem.fixed_branches = fixed_branches
        pricing_problem.initial_loc_availability = compute_limited_loc(ins,o,fixed_branches)
        
        return pricing_problem
    end
end

mutable struct Partial_path_sshape{T<:AbstractNode} <:AbstractPath
    node::T
    rcost::Float64
    k::UInt8
    remaining_free_stops::UInt8
    visited_nodes::Vector{T}
    collected_cuts::Vector{Bool}
    loc_reachability::Vector{Bool}
    cut_reachability::Vector{Bool}
    remaining_stops::Vector{UInt8}
    entered_by_bottom::Bool
end



function create_empty_label(ins::instance,pricing::Pricing_problem_sshape, mu_o::Float64,ionode::AbstractNode) 
    node = ionode
    rcost = - mu_o
    k = 0
    remaining_free = length(pricing.order.sku) - sum(pricing.mandatory_stops)
    visited_nodes = [node]
    collected_cuts = Bool[false for c in pricing.available_cuts]
    loc_reachability = Bool[true for l in ins.locList]
    cut_reachability = Bool[true for c in pricing.available_cuts]
    remaining_stops = copy(pricing.mandatory_stops)
    entered_by_bottom = true

    # we check the case where there is no free stop, in this case we update the reachability
    if remaining_free == 0
        for i in eachindex(loc_reachability)
            if remaining_stops[i] == 0
                loc_reachability[i] = false
            end
        end
    end

    # Then we put also to unreachable the locs where there is no initial availability
    for i in eachindex(pricing.initial_loc_availability)
        if pricing.initial_loc_availability[i] == 0
            loc_reachability[i] = false
        end
    end

    return Partial_path_sshape{typeof(ionode)}(node,rcost,k,remaining_free,visited_nodes,collected_cuts,loc_reachability,cut_reachability,remaining_stops,entered_by_bottom)
end

# copy methods
Base.copy(path::Partial_path_sshape{Node}) = Partial_path_sshape{Node}(path.node,path.rcost,path.k,
    path.remaining_free_stops,copy(path.visited_nodes),copy(path.collected_cuts),copy(path.loc_reachability),
    copy(path.cut_reachability),copy(path.remaining_stops),copy(path.entered_by_bottom))


##################################################################
# Data structures midpoint
##################################################################

mutable struct Pricing_problem_midpoint{T<:AbstractNode, P<:AbstractPath} <:AbstractPricing
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
    mandatory_stops::Vector{Int}
    fixed_branches::Dict{sku,loc}
    initial_loc_availability::Vector{UInt8}

    function Pricing_problem_midpoint{T,P}(ins::instance,o::structures_decomp.order,mu_o::Float64,pi_o::Vector{Float64},sigma_o::Vector{Float64},cut_pool::cuts.CutCollection, number_blocks_layout::Int, EPS::Float64,fixed_branches::Dict{sku,loc}) where {T <: AbstractNode, P <: AbstractPath}
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
        pricing_problem.distance = create_distance_matrix(ins,pricing_problem.node_list,pi_o,sigma_o,ins.number_blocks_layout)
        #pricing_problem.open_labels = P[Partial_path{T}(ins,pricing_problem,mu_o,pricing_problem.ionode)]
        pricing_problem.labels = Vector{P}[[] for n1 in pricing_problem.node_list]
        pricing_problem.closed_labels = Vector{P}[[] for n1 in pricing_problem.node_list]
        pricing_problem.current_k = 0
        pricing_problem.current_node = T(0,ioloc,0,pricing_problem)
        pricing_problem.mu_o = mu_o
        pricing_problem.pi_o = pi_o
        pricing_problem.sigma_o = sigma_o
        pricing_problem.EPS = EPS

        pricing_problem.mandatory_stops = collect_mandatory_stops(o,ins,fixed_branches)
        pricing_problem.fixed_branches = fixed_branches
        pricing_problem.initial_loc_availability = compute_limited_loc(ins,o,fixed_branches)
        
        return pricing_problem
    end
end

mutable struct Partial_path_midpoint{T<:AbstractNode} <:AbstractPath
    node::T
    rcost::Float64
    k::UInt8
    remaining_free_stops::UInt8
    visited_nodes::Vector{T}
    collected_cuts::Vector{Bool}
    loc_reachability::Vector{Bool}
    cut_reachability::Vector{Bool}
    remaining_stops::Vector{UInt8}
    first_aisle::UInt8
    last_aisle::UInt8
end



function create_empty_label(ins::instance,pricing::Pricing_problem_midpoint, mu_o::Float64,ionode::AbstractNode) 
    node = ionode
    rcost = - mu_o
    k = 0
    remaining_free = length(pricing.order.sku) - sum(pricing.mandatory_stops)
    visited_nodes = [node]
    collected_cuts = Bool[false for c in pricing.available_cuts]
    loc_reachability = Bool[true for l in ins.locList]
    cut_reachability = Bool[true for c in pricing.available_cuts]
    remaining_stops = copy(pricing.mandatory_stops)

    # we check the case where there is no free stop, in this case we update the reachability
    if remaining_free == 0
        for i in eachindex(loc_reachability)
            if remaining_stops[i] == 0
                loc_reachability[i] = false
            end
        end
    end

    # Then we put also to unreachable the locs where there is no initial availability
    for i in eachindex(pricing.initial_loc_availability)
        if pricing.initial_loc_availability[i] == 0
            loc_reachability[i] = false
        end
    end

    # first and last aisle for midpoint
    first_aisle = 0
    last_aisle = 0

    return Partial_path_midpoint{typeof(ionode)}(node,rcost,k,remaining_free,visited_nodes,collected_cuts,loc_reachability,cut_reachability,remaining_stops,first_aisle,last_aisle)
end

# copy methods
Base.copy(path::Partial_path_midpoint{Node}) = Partial_path_midpoint{Node}(path.node,path.rcost,path.k,
    path.remaining_free_stops,copy(path.visited_nodes),copy(path.collected_cuts),copy(path.loc_reachability),
    copy(path.cut_reachability),copy(path.remaining_stops),path.first_aisle,path.last_aisle)


##################################################################
# Data structures largest
##################################################################

mutable struct Pricing_problem_largest{T<:AbstractNode, P<:AbstractPath} <:AbstractPricing
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
    mandatory_stops::Vector{Int}
    fixed_branches::Dict{sku,loc}
    initial_loc_availability::Vector{UInt8}

    function Pricing_problem_largest{T,P}(ins::instance,o::structures_decomp.order,mu_o::Float64,pi_o::Vector{Float64},sigma_o::Vector{Float64},cut_pool::cuts.CutCollection, number_blocks_layout::Int, EPS::Float64,fixed_branches::Dict{sku,loc}) where {T <: AbstractNode, P <: AbstractPath}
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
        pricing_problem.distance = create_distance_matrix(ins,pricing_problem.node_list,pi_o,sigma_o,ins.number_blocks_layout)
        #pricing_problem.open_labels = P[Partial_path{T}(ins,pricing_problem,mu_o,pricing_problem.ionode)]
        pricing_problem.labels = Vector{P}[[] for n1 in pricing_problem.node_list]
        pricing_problem.closed_labels = Vector{P}[[] for n1 in pricing_problem.node_list]
        pricing_problem.current_k = 0
        pricing_problem.current_node = T(0,ioloc,0,pricing_problem)
        pricing_problem.mu_o = mu_o
        pricing_problem.pi_o = pi_o
        pricing_problem.sigma_o = sigma_o
        pricing_problem.EPS = EPS

        pricing_problem.mandatory_stops = collect_mandatory_stops(o,ins,fixed_branches)
        pricing_problem.fixed_branches = fixed_branches
        pricing_problem.initial_loc_availability = compute_limited_loc(ins,o,fixed_branches)
        
        return pricing_problem
    end
end

mutable struct Partial_path_largest{T<:AbstractNode} <:AbstractPath
    node::T
    rcost::Float64
    k::UInt8
    remaining_free_stops::UInt8
    visited_nodes::Vector{T}
    collected_cuts::Vector{Bool}
    loc_reachability::Vector{Bool}
    cut_reachability::Vector{Bool}
    remaining_stops::Vector{UInt8}
    first_aisle::UInt8
    last_aisle::UInt8
    # for the following, it is one paramtere per aisle
    last_col_from_top::Vector{UInt8} # the col of the last visited loc from the top
    middle_gap::Vector{UInt8} # the current value of the middle gap
    min_middle_gap::Vector{UInt8} # the minimum value for the middle gap, otherwise another is the
    # largest
end



function create_empty_label(ins::instance,pricing::Pricing_problem_largest, mu_o::Float64,ionode::AbstractNode) 
    node = ionode
    rcost = - mu_o
    k = 0
    remaining_free = length(pricing.order.sku) - sum(pricing.mandatory_stops)
    visited_nodes = [node]
    collected_cuts = Bool[false for c in pricing.available_cuts]
    loc_reachability = Bool[true for l in ins.locList]
    cut_reachability = Bool[true for c in pricing.available_cuts]
    remaining_stops = copy(pricing.mandatory_stops)

    # we check the case where there is no free stop, in this case we update the reachability
    if remaining_free == 0
        for i in eachindex(loc_reachability)
            if remaining_stops[i] == 0
                loc_reachability[i] = false
            end
        end
    end

    # Then we put also to unreachable the locs where there is no initial availability
    for i in eachindex(pricing.initial_loc_availability)
        if pricing.initial_loc_availability[i] == 0
            loc_reachability[i] = false
        end
    end

    # first and last aisle for midpoint
    first_aisle = 0
    last_aisle = 0
    last_col_from_top = UInt8[ins.cmax+1 for a=1:ins.amax]
    middle_gap = UInt8[ins.cmax+1 for a=1:ins.amax]
    min_middle_gap = UInt8[0 for a=1:ins.amax]

    return Partial_path_largest{typeof(ionode)}(node,rcost,k,remaining_free,visited_nodes,collected_cuts,loc_reachability,cut_reachability,remaining_stops,first_aisle,last_aisle,last_col_from_top,middle_gap,min_middle_gap)
end

# copy methods
Base.copy(path::Partial_path_largest{Node}) = Partial_path_largest{Node}(path.node,path.rcost,path.k,
    path.remaining_free_stops,copy(path.visited_nodes),copy(path.collected_cuts),copy(path.loc_reachability),
    copy(path.cut_reachability),copy(path.remaining_stops),path.first_aisle,path.last_aisle,
    copy(path.last_col_from_top),copy(path.middle_gap),copy(path.min_middle_gap))


################################################################
# Main interface functions
################################################################




# function that creates the node list
function create_nodelist(ins::instance,pricing_problem::AbstractPricing) 
    node_list = Vector{Node}()
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

function create_pricing_problem(ins::instance,o::order,mu_o::Float64,pi_o::Vector{Float64},sigma_o::Vector{Float64},cut_pool::cuts.CutCollection,EPS::Float64,policy::String,fixed_branches::Dict{sku,loc})
    # Interface function for creating the pricing problem
    if policy == "optimal"
        return Pricing_problem_optimal{Node,Partial_path_optimal}(ins,o,mu_o,pi_o,sigma_o,cut_pool,ins.number_blocks_layout,EPS,fixed_branches)
    elseif policy == "return"
        return Pricing_problem_return{Node,Partial_path_return}(ins,o,mu_o,pi_o,sigma_o,cut_pool,ins.number_blocks_layout,EPS,fixed_branches)
    elseif policy == "midpoint"
        return Pricing_problem_midpoint{Node,Partial_path_midpoint}(ins,o,mu_o,pi_o,sigma_o,cut_pool,ins.number_blocks_layout,EPS,fixed_branches)
    elseif policy == "sshape"
        return Pricing_problem_sshape{Node,Partial_path_sshape}(ins,o,mu_o,pi_o,sigma_o,cut_pool,ins.number_blocks_layout,EPS,fixed_branches)
    elseif policy == "largest"
        return Pricing_problem_largest{Node,Partial_path_largest}(ins,o,mu_o,pi_o,sigma_o,cut_pool,ins.number_blocks_layout,EPS,fixed_branches)
    elseif policy == "guo_return"
        return Pricing_problem_guo{Node,Partial_path_guo}(ins,o,mu_o,pi_o,sigma_o,cut_pool,ins.number_blocks_layout,EPS,fixed_branches)
    else
        throw(DomainError("Unknown routing policy"))
    end
end

function create_empty_label(ins::instance,pricing::AbstractPricing,mu_o::Float64,ionode::AbstractNode) end

#########################################################
# Cost functions
#########################################################

# This function gives the distance by passing by the top cross aisle
function dist_locs_by_top(a::loc,b::loc,ins::instance)
    if a.aisle == b.aisle
        return ins.wc*abs(a.col - b.col)
    else
        return ins.wa*abs(a.aisle - b.aisle) + ins.wc*(2*(ins.cmax + 1) - a.col - b.col)
    end
end

function dist_locs_by_bottom(a::loc,b::loc,ins::instance)
    if a.aisle == b.aisle
        return ins.wc*abs(a.col - b.col)
    else
        return ins.wa*abs(a.aisle - b.aisle) + ins.wc*(a.col + b.col)
    end
end