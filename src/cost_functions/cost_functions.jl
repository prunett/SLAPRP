module cost_functions

using ..structures_decomp, ..cost_functions_exact, ..cost_functions_return, ..cost_functions_guo, ..cost_functions_midpoint, ..cost_functions_largest, ..cost_functions_sshape

export dist, route_cost, sol_cost, create_solution

function dist(a::loc,b::loc,ins::instance, policy::String = ins.policy)
    if policy == "optimal" || policy == "exact_old"
        return cost_functions_exact.dist(a,b,ins)
    elseif policy == "return" || policy == "return_old"
        return cost_functions_return.dist(a,b,ins)
    else
        throw(DomainError("Unknown routing policy"))
    end
end

dist(ins::instance,a::loc,b::loc,policy::String = ins.policy) = dist(a,b,ins,policy)

function route_cost(ins::instance,loc_group::Vector{loc},o::order,partial::Bool = false, policy::String = ins.policy)
    if policy == "optimal" || policy == "exact_old"
        return cost_functions_exact.route_cost(ins,loc_group,o,partial)
    elseif policy == "return" || policy == "return_old"
        return cost_functions_return.route_cost(ins,loc_group,o,partial)
    elseif policy == "midpoint"
        return cost_functions_midpoint.route_cost(ins,loc_group,o,partial)
    elseif policy == "sshape"
        return cost_functions_sshape.route_cost(ins,loc_group,o,partial)
    elseif policy == "largest"
        return cost_functions_largest.route_cost(ins,loc_group,o,partial)
    elseif policy == "guo_return"
        return cost_functions_guo.route_cost(ins,loc_group,o,partial)
    else
        throw(DomainError("Unknown routing policy"))
    end
end

# new method
function route_cost(ins::instance, stol::Dict{sku,loc}, o::order, partial::Bool = false, policy::String = ins.policy)
    loc_group = loc[stol[s] for s in o.sku]
    return route_cost(ins,loc_group,o,partial,policy)
end

function sol_cost(ins::instance,stol::Dict{sku,loc},partial::Bool = false, policy::String = ins.policy)
    cost = 0
    for o in ins.orderList
        cost += route_cost(ins,stol,o,partial,policy).cost
    end
    return cost
end


function create_solution(ins::instance,stol::Dict{sku,loc}, partial::Bool = false, policy::String = ins.policy)
# Creates a solution type
    
    routes = Array{route,1}()
    for o in ins.orderList
        push!(routes, route_cost(ins,stol,o,partial,policy))
    end
    
    return solution(sol_cost(ins,stol,partial,policy),routes,stol)
end

function delta_insertion_main(ins::instance, route_array::Array{route,1}, occupation::Dict{loc,Int}, s::sku, l::loc,policy::String = ins.policy)
    if policy == "optimal" || policy == "exact_old"
        return cost_functions_exact.delta_insertion_main(ins,route_array,occupation,s,l)
    elseif policy == "return" || policy == "return_old"
        return cost_functions_return.delta_insertion_main(ins,route_array,occupation,s,l)
    else
        throw(DomainError("Unknown routing policy"))
    end
end


end