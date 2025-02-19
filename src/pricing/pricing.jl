module pricing

using ..structures_decomp, ..cuts, ..pricing_return_module, ..pricing_exact_module, ..pricing_general

export solve_pricing

function solve_pricing(ins::instance,o::order,mu_o::Float64,pi_o::Array{Float64,1},sigma_o::Array{Float64,1},cut_pool::cuts.CutCollection,bi_directional::Bool,fixed_branches::Dict{sku,loc},EPS::Float64, policy::String = ins.policy)
    if policy == "exact_old"
        if length(o.sku) >= 4
            return pricing_exact_module.pricing_exact(ins,o,mu_o,pi_o,sigma_o,cut_pool,bi_directional,EPS)
        else
            return pricing_exact_module.pricing_exact(ins,o,mu_o,pi_o,sigma_o,cut_pool,false,EPS)
        end
    elseif policy == "return_old"
        return pricing_return_module.pricing_improved(ins,o,mu_o,pi_o,sigma_o,cut_pool,EPS)
    else
        return pricing_general.pricing(ins,o,mu_o,pi_o,sigma_o,cut_pool,fixed_branches,EPS,policy)
    end
end


end