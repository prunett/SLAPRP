module CG_utils

using JuMP, CPLEX
using ..structures_decomp
using ..cost_functions
using ..cuts
using ..BB


export visualisation_plan, add_route_to_master!, write_model_to_file

EPS = .001


function read_xi_name(xi::VariableRef)
# This function reads the name of a xi variable
    name = JuMP.name(xi)
    lid = 0
    sid = 0
    ind = 1
    while lid == 0 || sid == 0
        if name[ind] == ','
            ind += 1
            ind2 = ind+1
            while name[ind2] != ']'
                ind2 += 1
            end
            sid = parse(Int,name[ind:ind2-1])
        end
        if name[ind] == '['
            ind += 1
            ind2 = ind+1
            while name[ind2] != ','
                ind2 += 1
            end
            lid = parse(Int,name[ind:ind2-1])
        else
            ind += 1
        end
    end
    return lid, sid
end

function checker(ins::instance,master_problem::JuMP.Model,xi_val,rho_val,a_coef,EPS)
# This function is a checker and ensure that a given solution does not violate any constraint

    check = true

    n = size(xi_val,2)
    n_loc = size(xi_val,1)

    for s = 1:n
        if abs(sum(xi_val[:,s]) - 1) > EPS
            println("xi problem on s")
            check = false
        end
    end

    for l = 1:n_loc
        if sum(xi_val[l,:]) > ins.n_sym + EPS
            println("xi problem on l")
            check = false
        end
    end

    for o in ins.orderList
        if sum(rho_val[o]) < 1 - EPS
            println("problem for rho sum order $(o.id) = $(sum(rho_val[o]))")
            check = false
        end
    end

    for l in ins.locList, o in ins.orderList
        if sum( count(l2 -> l2.id == l.id, a_coef[o][r].loc)*rho_val[o][r] for r=1:length(a_coef[o])) < (isempty(o.sku) ? 0 : sum(xi_val[l.id,s.id] for s in o.sku)) - EPS
            println("problem for link constraint l = $(l.id), o= $(o.id), sum 1 = $(sum( count(l2 -> l2.id == l.id, a_coef[o][r].loc)*rho_val[o][r] for r=1:length(a_coef[o]))), sum 2 = $((isempty(o.sku) ? 0 : sum(xi_val[l.id,s.id] for s in o.sku)))")
            check = false
        end
    end

    # look at fixed positions
    for s in keys(ins.fixed_positions)
        if xi_val[ins.fixed_positions[s].id,s.id] < 1 - EPS
            check = false
        end
    end

    return check
end


function create_master_problem(ins::instance,cpu_limit::Int,simplex_algo::Int)
# This function creates the mastre problem for CG stable
    n = ins.n
    m = ins.nOrd
    n_loc = length(ins.locList)
    orderList = ins.orderList

    # Create model with CPLEX
    
    master_problem = Model(CPLEX.Optimizer)
    # Remove output
    set_optimizer_attribute(master_problem, "CPX_PARAM_SCRIND", 0)
    # Only one thread
    set_optimizer_attribute(master_problem, "CPXPARAM_Threads", cpu_limit)
    # Changing the lp solver of cplex
    # 1 = primal simplex, 2 = dual simplex, 3 = network simplex, 4 = barrier, 5 = sifting
    set_optimizer_attribute(master_problem, "CPXPARAM_LPMethod", simplex_algo)
    # TEST
    #set_optimizer_attribute(master_problem, "CPXPARAM_Preprocessing_Reduce", 3)
    # dependancy checker
    #set_optimizer_attribute(master_problem, "CPXPARAM_Preprocessing_Dependency", 3)
    #master_problem = Model(Gurobi.Optimizer)
    #set_optimizer_attribute(master_problem, "OutputFlag", 0)
    #set_optimizer_attribute(master_problem, "Threads", 1)


    # Xi variables
    @variable(master_problem, xi[1:n_loc,1:n], lower_bound = 0., upper_bound = 1.)#, upper_bound = 1.)
    # Constraints
    @constraint(master_problem, xi_sum_1[l=1:n_loc], sum(xi[l,s] for s=1:n) <= ins.n_sym) # <=
    @constraint(master_problem, xi_sum_2[s=1:n], sum(xi[l,s] for l=1:n_loc) == 1)
    @constraint(master_problem, rho_sum[1:m], 0 >= 1) # >=
    @constraint(master_problem, rho_sum2[1:m], 0 <= 1) # <=
    #@constraint(master_problem, link_constraint[l=1:n_loc,o=1:m,s in [s2.id for s2 in orderList[o].sku]], xi[l,s] <= 0) # <=
    @constraint(master_problem, link_constraint[l=1:n_loc,o=1:m], sum(xi[l,s.id] for s in orderList[o].sku) <= 0) # <=
    #@constraint(master_problem, link_constraint_bis[l=1:n_loc,o=1:m, s in [s.id for s in orderList[o].sku]], xi[l,s] <= 0)
    #Defining Objective
    @objective(master_problem, Min, 0)

    # set up the fixed positions constraints
    for s in keys(ins.fixed_positions)
        @constraint(master_problem, xi[ins.fixed_positions[s].id,s.id] >= 1)
    end

    return master_problem
end

function create_pricing_problem(ins::instance,o::order)
    return cost_functions.create_pricing_problem_pure(ins,o)
end

function update_pricing_obj!(ins::instance,pricing_route_o::Model,mu_o::Float64,pi_o::Array{Float64,1})
# This function updates the pricing problem objective

    cost_functions.update_pricing_obj_pure!(ins,pricing_route_o,mu_o,pi_o)
end


function add_route_to_master!(r::route,master_problem::Model,rho_vars_o::Array{VariableRef,1},cut_pool::cuts.CutCollection)
# This function takes a route into argument and add to the master problem the corresponding variable

    # Retrieveing data from inputs
    rho_sum = master_problem[:rho_sum]
    rho_sum2 = master_problem[:rho_sum2]
    link_constraint = master_problem[:link_constraint]
    #link_constraint_bis = master_problem[:link_constraint_bis]
    o = r.order
    # Creation of arrays to register modified constraints
    touched_constraints = ConstraintRef[]
    constraint_coefs = Float64[]
    #sum constraint
    push!(touched_constraints, rho_sum[o.id])
    push!(constraint_coefs, 1)
    # sum 2 constraint
    push!(touched_constraints, rho_sum2[o.id])
    push!(constraint_coefs, 1)

    for l in r.loc
        push!(touched_constraints, link_constraint[l.id,o.id])
        push!(constraint_coefs, - count(l2 -> (l2.id == l.id), r.loc))
    end
    #link constraints bis
    #=for l in r.loc
        for s in o.sku
            push!(touched_constraints, link_constraint_bis[l.id,o.id,s.id])
            push!(constraint_coefs, -1)
        end
    end=#

    # Add to active cuts
    cuts.add_route_to_cut!(r,cut_pool,master_problem,touched_constraints,constraint_coefs, length(rho_vars_o)+1)


    # Add the variable to model
    name = string("rho[",o.id,",", size(rho_vars_o,1) + 1,"]")
    rho_jump = @variable(master_problem, lower_bound = 0.)
    set_name(rho_jump, name)
    # Add coefficient to objective
    set_objective_coefficient(master_problem, rho_jump, Float64(r.cost))
    # Add constraints coefficients
    for i=1:length(touched_constraints)
        set_normalized_coefficient(touched_constraints[i], rho_jump, constraint_coefs[i])
    end
    # Pushing the new variable to the array
    push!(rho_vars_o, rho_jump)
end

function add_route_to_master!(routes::Array{route,1},master_problem::Model,rho_vars::Dict{order,Array{VariableRef,1}},cut_pool::cuts.CutCollection)
    for r in routes
        add_route_to_master!(r, master_problem, rho_vars[r.order],cut_pool)
    end
end
#=
function add_route_to_master!(sol::solution,master_problem::Model,rho_vars::Dict{order,Array{VariableRef,1}},ins::instance)
# Add all routes
    for r in sol.routes
        add_route_to_master!(r,master_problem,rho_vars[r.order])
        ### add
        #add_route_to_master!(reverse_route(ins,r),master_problem,rho_vars[r.order])
    end
end=#





function visualisation_plan(ins::instance,ltos::Dict{loc,sku})
# Gives a visualisation of the plan
    loc_aisle = Array{Array{loc,1},1}()
    for a=1:ins.amax
        push!(loc_aisle, filter(l -> l.aisle == a, ins.locList))
    end
    map = []
    for a=1:ins.amax
        aisle_g = []
        aisle_d = []
        for l in loc_aisle[a]
            for c in StepRange(ins.cmax,-1,Int8(1))
                if l.col == c
                    if l.side == 1
                        if l in keys(ltos)
                            push!(aisle_g,ltos[l].id)
                        else
                            push!(aisle_g,0)
                        end
                    else
                        if l in keys(ltos)
                            push!(aisle_d,ltos[l].id)
                        else
                            push!(aisle_d,0)
                        end
                    end
                end
            end
        end
        push!(map,aisle_g)
        push!(map,aisle_d)
        push!(map, "///////////")
    end

    for i in map
        println(i)
    end
end

visualisation_plan(ins::instance, sol::solution) = visualisation_plan(ins, sol.assignment.ltos)

function write_model_to_file(model::Model, file_name::String)

    # Write pricing problem to file
    f = open(file_name, "w")
    print(f, model)
    close(f)

end

function build_fake_routes(ins::instance)
# This function builds fake routes to ensure feasability of the master problem at all time
    route_list = route[]
    for o in ins.orderList
        loc_group = loc[]
        for l in ins.locList
            push!(loc_group,l)
            push!(loc_group,l)
        end
        r = route(o,loc_group,1000*(ins.wa*ins.amax + ins.wc*ins.cmax*ins.amax))
        push!(route_list,r)
    end
    return route_list
end




end # end of module
