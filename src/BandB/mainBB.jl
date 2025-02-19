module BB

using JuMP
using MathOptInterface

const MOI = MathOptInterface

include("node_choices.jl")
include("tree.jl")

function next_node!(tree::SearchTree)
# This function is the search strategy, for now it is LIFO
    node = popfirst!(tree.nodes)
    return node
end

end # end module
