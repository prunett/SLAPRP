mutable struct SearchTree
    node_counter::Int
    nodes::Vector{AbstractNode}
    processed::Vector{AbstractNode}
    current_dual_bound::Real
    current_primal_bound::Real
    best_primal_solution::Union{Nothing,Dict}
    auxiliary_data::Dict # used to store the instance for example

    # Constructor, in the case of MINIMIZATION
    function SearchTree(node_counter::Int = 0)
        return new(node_counter, [], [], - Inf, Inf, nothing, Dict())
    end
end

# main runing algorithm
function run(tree::SearchTree,max_time)
    # UNIT we stop at the root node
    start_time = time()
    while !termination(tree,max_time) #&& isempty(tree.processed)
        # Array of the messages to print
        message = Dict()
        # Get the node to process according to the search strategy
        node = next_node!(tree)
        # Check if the algorithm should terminate with priority queue
        if prune_bound(tree,node)
            if tree.auxiliary_data["config"].search_strategy == "priority_queue"
                # delete all unexplored node (for checking optimality)
                tree.nodes = []
                # Then break the procedure
                break
            elseif tree.auxiliary_data["config"].search_strategy in ["fifo","lifo"]
                continue
            end
        end
        # Apply branches choices to the model
        apply_changes!(node)
        # Compute dual bound
        message["dual"] = dual_bound_node!(node,tree,max_time)
        # Compute primal bound
        primal_bound_node!(node,tree)
        # Update best primal bound
        update_best_primal_bound!(tree,node)
        # Mark node as processed
        processed!(tree,node)
        # Check if it is the root node. In that case we print the solution & obj
        if node.id == 1
            #println("node 159 !!")
            #break
            #print_solution(tree,node)
            # UNIT
            #JuMP.write_to_file(node.model,"master_problem.lp")
        end
        # Check if the node is pruned, or if we branch on it
        prune = prune_node(tree,node)
        # branch to create child if the node has not been pruned
        if !prune
            message["branch"] = branch!(tree,node)
        else
            message["branch"] = "The node has been pruned !"
        end
        # Update current dual bound
        update_current_dual_bound!(tree)
        # Reverse branching choices
        reverse_changes!(node)
        # Print summary
        print_iteration_summary(tree,node,message,start_time)
        # UNIT that cut after the first iterationn
        #@assert node.id <= 3
        #print_solution(tree,node)
    end
end

# basic termination function
termination(tree::SearchTree,max_time) = isempty(tree) || time() > max_time
Base.isempty(tree::SearchTree) = Base.isempty(tree.nodes)

# Add node to the tree
function push!(tree::SearchTree, node::AbstractNode)
    tree.node_counter += 1
    node.id = tree.node_counter
    Base.push!(tree.nodes,node)
    # sort!(tree)
end

# Add several nodes
function push!(tree::SearchTree, nodes::Vector{T}) where T<:AbstractNode
    for node in nodes
        push!(tree,node)
    end
end

# Mark node as processed
processed!(tree::SearchTree, node::AbstractNode) = Base.push!(tree.processed, node)

# update the current dual bound
function update_current_dual_bound!(tree::SearchTree)
    if !isempty(tree)
        tree.current_dual_bound = min(Base.minimum([node.dual_bound for node in tree.nodes]), tree.auxiliary_data["min_bound_close_nodes"])
    elseif length(tree.processed) == 1
        tree.current_dual_bound = tree.processed[1].dual_bound
    end
end
# update the best primal bound
function update_best_primal_bound!(tree::SearchTree,node::AbstractNode)
    if node.primal_bound < tree.current_primal_bound
        tree.current_primal_bound = node.primal_bound
    end
    # UNIT
    #tree.current_primal_bound = 1000
end

# branch
function branch!(tree::SearchTree, node::AbstractNode)
    children = branch(node)
    push!(tree,children)
end

# empty summary function
function print_iteration_summary(tree::SearchTree,node::AbstractNode,message::Dict)
    println("Iteration summary not defined")
end

#empty branching function
branch(node::AbstractNode)::Vector{AbstractNode} = []

#empty bouding function
function dual_bound_node!(node::AbstractNode,tree::SearchTree) end

#empty primal solution function
function primal_bound_node!(node::AbstractNode,tree::SearchTree) end

# empty pruning function
prune_node(tree,node::AbstractNode) = nothing
function prune_bound(tree::SearchTree,node::AbstractNode) end

#empty printing function
function print_solution(node::AbstractNode,tree::SearchTree) end
