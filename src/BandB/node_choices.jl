using JuMP

abstract type AbstractNode end
abstract type AbstractBranch end
abstract type AbstractSolution end

## General abstract Node definitions ##

macro abstract_node_fields()
    return esc(quote
        id::Int
        parent::Union{Nothing,AbstractNode}
        depth::Int
        branch::Union{Nothing,T}
        solution_status::MOI.TerminationStatusCode
        primal_bound::Real
        dual_bound::Real
        solution::Dict{Any,Real}
    end)
end



# apply choices (branches) from ancestor to child
function apply_changes!(node::AbstractNode, ancestor::AbstractNode) end

#function reverse_changes!(node::AbstractNode) end

apply_changes!(node::AbstractNode) = apply_changes!(node,node)

##############################
# JumpNode definition & functions
##############################

mutable struct JumpNode{T<:AbstractBranch} <:AbstractNode
    @abstract_node_fields

    model::Union{Nothing,JuMP.Model}
    auxiliary_data::Dict

    # Constructor
    function JumpNode{T}(
        model = nothing,
        parent = nothing,
        depth = 0,
        dual_bound = -Inf,
        branch = nothing,
        id = -1,
        solution_status = MOI.OPTIMIZE_NOT_CALLED,
        primal_bound = Inf,
        solution = Dict{Any,Real}(),
        auxiliary_data = Dict()) where T <:AbstractBranch
        return new{T}(id, parent, depth, branch, solution_status, primal_bound, dual_bound, solution,model, auxiliary_data)
    end
end

JumpNode{T}(parent::JumpNode{T}) where T<:AbstractBranch = JumpNode{T}(nothing, parent, parent.depth+1, parent.dual_bound)
JumpNode(parent::JumpNode{T}) where T <:AbstractBranch = JumpNode{T}(parent)


## Branches definitions ##

# Basic branching on one variable structure
struct VariableBranch <: AbstractBranch
    lb::Dict{JuMP.VariableRef,Real}
    ub::Dict{JuMP.VariableRef,Real}
    fixed_sku_id::Union{Nothing,Int}
    fixed_loc_id::Union{Nothing,Int}
    fixed_sku_aisle::Union{Nothing,Int}
end
# Sum of variables branch
mutable struct SumBranch <: AbstractBranch
    sum_variables::Array{JuMP.VariableRef,1}
    lb::Union{Nothing,Real}
    ub::Union{Nothing,Real}
    constraint_ub::Union{JuMP.ConstraintRef,Nothing}
    constraint_lb::Union{JuMP.ConstraintRef,Nothing}
end


empty_variable_bound() = Dict{JuMP.VariableRef,Real}()


# Create child node with branch object
function create_child_node(current_node::JumpNode{T}, branch::AbstractBranch) where T<:AbstractBranch
    node = JumpNode(current_node)
    node.branch = branch
    return node
end

# Create child node with variable bounds
create_child_node(current_node::JumpNode{T}, variable::JuMP.VariableRef, lb::Real, ub::Real) where T<:AbstractBranch = create_child_node(current_node,VariableBranch(Dict(variable=>lb),Dict(variable=>ub),nothing,nothing,nothing))
create_child_node_with_lb(current_node::JumpNode{T}, variable::JuMP.VariableRef, lb::Real) where T<:AbstractBranch = create_child_node(current_node,VariableBranch(Dict(variable=>lb), empty_variable_bound(),nothing,nothing,nothing))
create_child_node_with_ub(current_node::JumpNode{T}, variable::JuMP.VariableRef, ub::Real) where T<:AbstractBranch = create_child_node(current_node,VariableBranch(empty_variable_bound(), Dict(variable=>ub),nothing,nothing,nothing))
# method for symmetry => multiple variables with 
create_child_node_multiple0(current_node::JumpNode{VariableBranch}, variable_zeros::Vector{JuMP.VariableRef},fixed_sku_aisle) = create_child_node(current_node,VariableBranch(empty_variable_bound(), Dict(v => 0 for v in variable_zeros),nothing,nothing,fixed_sku_aisle))
create_child_node_multiple0(current_node::JumpNode{VariableBranch}, variable_zeros::Vector{JuMP.VariableRef}) = create_child_node_multiple0(current_node, variable_zeros, nothing)
# method for partial replenishment
create_child_node_fixed_one(current_node::JumpNode{T}, variable::JuMP.VariableRef, sku_id::Int,loc_id::Int) where T<:AbstractBranch = create_child_node(current_node,VariableBranch(Dict(variable=>1), empty_variable_bound(),sku_id,loc_id,nothing))

# Create child node with sum bound branching constraint
create_child_node(current_node::JumpNode{T}, sum_variables::Array{JuMP.VariableRef,1}, lb::Union{Nothing,Real}, ub::Union{Nothing,Real}) where T<:AbstractBranch = create_child_node(current_node,SumBranch(sum_variables,lb,ub,nothing,nothing))

Base.push!(a::Array{T,1},b::Nothing) where T<:AbstractBranch = a

# Collect all branches objects
function collect_all_branches(node::JumpNode{T}) where T<:AbstractBranch
    branch_objects = Array{T,1}()

    # backtrack the tree to collect all the branches
    pnode = node
    Base.push!(branch_objects, pnode.branch)
    while !isnothing(pnode.parent)
        pnode = pnode.parent
        Base.push!(branch_objects, pnode.branch)
    end
    @assert isnothing(pnode.parent)
    @assert !isnothing(pnode.model)

    return branch_objects
end

#set_normalized_rhs
# apply_changes function that set the constraints into the jump model
function apply_changes!(node::JumpNode{T}) where T<:AbstractBranch

    branch_objects = collect_all_branches(node)

    # Shallow copy of the root model
    pnode = node
    while !isnothing(pnode.parent)
        pnode = pnode.parent
    end
    node.model = pnode.model

    # mark bounds to be changed to reverse to previous state
    mark_bound_changes!(node, branch_objects)

    # Apply branch objects
    apply_branch_objects(branch_objects)
    #println("length of branch objects applied = $(length(branch_objects))")
end

# This function enforce the variable branches into the model
function apply_branch_objects(branch_objects::Array{VariableBranch,1})
    for branch in branch_objects
        for (j,v) in branch.lb
            JuMP.set_lower_bound(j, v)
        end
        for (j,v) in branch.ub
            if !isinf(v)
                JuMP.set_upper_bound(j, v)
            else
                JuMP.delete_upper_bound(j)
            end
        end
    end
end

function apply_branch_objects(branch_objects::Array{SumBranch,1})
    for branch in branch_objects
        # first test if the branch is empty
        if isempty(branch.sum_variables)
            continue
        end
        # create the lb constraint
        if !isnothing(branch.lb)
            branch.constraint_lb = @constraint(JuMP.owner_model(branch.sum_variables[1]), sum(v for v in branch.sum_variables) >= branch.lb)
        end
        # create the ub constraint
        if !isnothing(branch.ub)
            branch.constraint_ub = @constraint(JuMP.owner_model(branch.sum_variables[1]), sum(v for v in branch.sum_variables) <= branch.ub)
        end
    end
end



# This function makes a list of branch to be applied to reverse all the choices of branch_objects
function mark_bound_changes!(node::JumpNode, branch_objects::Array{VariableBranch,1})
    node.auxiliary_data["bounds_changed"] = VariableBranch[]
    for branch in branch_objects
        lb_changed = Dict{JuMP.VariableRef,Real}()
        ub_changed = Dict{JuMP.VariableRef,Real}()
        for (j,v) in branch.lb
            lb_changed[j] = JuMP.has_lower_bound(j) ? JuMP.lower_bound(j) : -Inf
        end
        for (j,v) in branch.ub
            ub_changed[j] = JuMP.has_upper_bound(j) ? JuMP.upper_bound(j) : Inf
        end
        Base.push!(node.auxiliary_data["bounds_changed"], VariableBranch(lb_changed, ub_changed,nothing,nothing,nothing))
    end
end
# in case of sum branch it is useless
function mark_bound_changes!(node::JumpNode, branch_objects::Array{SumBranch,1}) end

# Reverse changes in the model to get back at the root model
function reverse_changes!(node::JumpNode{VariableBranch})
    apply_branch_objects(node.auxiliary_data["bounds_changed"])
end
# In case of sum branch, we need to delete the constraints
function reverse_changes!(node::JumpNode{SumBranch})
    branch_objects = collect_all_branches(node)
    for branch in branch_objects
        if !isnothing(branch.constraint_lb)
            JuMP.delete(JuMP.owner_model(branch.constraint_lb), branch.constraint_lb)
            branch.constraint_lb = nothing
        end
        if !isnothing(branch.constraint_ub)
            JuMP.delete(JuMP.owner_model(branch.constraint_ub), branch.constraint_ub)
            branch.constraint_ub = nothing
        end
    end
end



#########################
# Tree initialization
##########################

function initialize_tree(root::JumpNode{T}) ::SearchTree where T <:AbstractBranch
    tree = SearchTree()
    push!(tree,root)
    return tree
end
initialize_tree(model::JuMP.Model, T = AbstractBranch)::SearchTree = initialize_tree(JumpNode{T}(model))
