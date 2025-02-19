include("structures_decomp.jl")
include("main.jl")

module launch

using ..structures_decomp
using ..SLAP_main_module
using DataFrames
using CSV
using StatsBase

function run()

###################################
# Parameters
###################################
time_limit = 2000
cpu_limit = Threads.nthreads()
EPS = 0.001

max_iter_cut =  1 # before deactivation
EPS_cut_separation_initial = 0.1
EPS_cut_violated_initial = 0.1
EPS_cut_vio_SL1 = 0.025
max_iter_SL1 = 10
EPS_delta_separation = 0. # unused
EPS_delta_violated = 0. # unused
max_round_OG = 6
max_cut_per_order = 500 # max number of active cut for an order

max_depth_cut_generation = 3
min_nodes_with_cuts = 5 # the first xxx nodes searched will have cuts generated
EPS_cut_obj_improvement = 0.01

# cut intiation
activate_all_link_cuts = true
initiate_aisle_cuts_in_pool = true

# lift and project
landp_main = "none" # can be none or classic_z
EPS_xi_landp = 0.1
EPS_landp_cut_violated = 0.01
landp_max_depth = 3
landp_max_rounds = 5
landp_strenghtening = true
landp_separate_zol_disjunction = true

# branching
branching_strategy = "aisle_without" # can be variable, sum_gub, sum_gl or symmetry or strong or aisle or aisle_without
# GUB greedy heuristic parameters
branching_weight_gub_demand = 20#parameter_list[parse(Int,ARGS[2])][1]
branching_weight_gub_xinumber = 0#parameter_list[parse(Int,ARGS[2])][2]
branching_weight_gub_balance = 0
branching_weight_gub_closeness = 5
# GL greedy heuristic parameters
branching_weight_gl_distance = 0#parameter_list[parse(Int,ARGS[2])][2]
branching_weight_gl_demand = 20#parameter_list[parse(Int,ARGS[2])][1]
branching_weight_gl_xinumber = 0
branching_weight_gl_closeness = 5

# parameters for strong branching
nb_candidates = 10
bound_method = "cg" # can be basic or cg

# partial heuristic noise
partial_heuristic_noise = 0.1

#search strategy
search_strategy = "priority_queue" # can be priority_queue, lifo or fifo

# fake route in the model
fake_route = true

# routing policy can be exact_old, return_old, optimal, midpoint, return, sshape, largest, guo_return
routing_policy = "optimal"
bi_directional = false

# instance type can be classic or silva
instance_type = "silva"

# simplex algo 2 = dual simplex, 4 = barrier
simplex_algo = 2



main_config = SLAP_main_module.bcp_config(
    time_limit,
    cpu_limit,
    EPS,
    instance_type,
    max_iter_cut,
    EPS_cut_separation_initial,
    EPS_cut_violated_initial,
    EPS_delta_separation,
    EPS_delta_violated,
    EPS_cut_obj_improvement,
    EPS_cut_vio_SL1,
    max_iter_SL1,
    max_round_OG,
    min_nodes_with_cuts,
    max_depth_cut_generation,
    search_strategy,
    branching_strategy,
    branching_weight_gub_demand,
    branching_weight_gub_xinumber,
    branching_weight_gub_balance,
    branching_weight_gub_closeness,
    branching_weight_gl_distance,
    branching_weight_gl_demand,
    branching_weight_gl_xinumber,
    branching_weight_gl_closeness,
    nb_candidates,
    bound_method,
    partial_heuristic_noise,
    fake_route,
    routing_policy,
    bi_directional,
    landp_main,
    EPS_xi_landp,
    EPS_landp_cut_violated,
    landp_max_depth,
    landp_max_rounds,
    landp_strenghtening,
    landp_separate_zol_disjunction,
    activate_all_link_cuts,
    initiate_aisle_cuts_in_pool,
    simplex_algo,
    max_cut_per_order
)


path = joinpath("C:\\Users\\thibault.prunet\\Documents\\Code\\decomposition\\test_instances","test_CG_1.txt")
#path = joinpath("C:\\Users\\thibault.prunet\\Documents\\Code\\decomposition\\instances_second_set","instance_T2_100_50_2.txt")
#path = joinpath("C:\\Users\\thibault.prunet\\Documents\\Code\\decomposition\\instances_LG1","instance_TG1_50_20_0.5_1.txt")
#path = joinpath("C:\\Users\\thibault.prunet\\Documents\\Code\\decomposition\\SLAP_CG","instance_test_complex_pricing.txt")

#silva instances
path = joinpath("C:\\Users\\thibault.prunet\\Documents\\Experiments\\silva_instances\\Instances\\Instances\\Small",
"SLAP-PRP_A5_B10_O10_I5_v2.txt")
#path = joinpath("C:\\Users\\thibault.prunet\\Documents\\Experiments\\silva_instances\\Instances\\Instances\\Regular",
#"SLAP-PRP_A5_B10_O10_I5_D3_v2.txt")
# new odysseus instances
#path = joinpath("C:\\Users\\thibault.prunet\\Documents\\Code\\decomposition\\instances_odysseus_classic",
#"classic_A1_B10_O20_I3_Z0.737_v1.txt")
#path = joinpath("C:\\Users\\thibault.prunet\\Documents\\Experiments\\silva_instances\\Instances\\Instances\\partial_instances",
#"SLAP-PRP_A10_B10_O50_I5_D1_a0.3_v1.txt")
#"SLAP-PRP_A5_B10_O30_I20_D1_a0.3_v1.txt")

# extended instances
#path = joinpath("C:\\Users\\thibault.prunet\\Documents\\Code\\decomposition\\SLAP_CG\\instances_silva_extended",
#"SLAPRP_A5_B10_O20_I10_Z33_v2.txt")

# guo instances
#path = joinpath("C:\\Users\\thibault.prunet\\Documents\\Experiments\\Partial_instances_paper\\small_instances_guo",
#"SLAPRP_Guo_small_O50_alpha0.2_v1.txt")

println("#######")
println("########")
println("thread = $(Threads.nthreads())")
tree = SLAP_main_module.main_BP_function(path,main_config)

#println("processed = $(length(tree.processed))")
println("dual bound root = $((isempty(tree.processed) ? 0 : tree.processed[1].dual_bound))")





#=


instance_list = instance[]
path_dir = joinpath(pwd(),ARGS[1])
file_name_list = readdir(path_dir)
for name in file_name_list
    #path = joinpath("C:\\Users\\thibault.prunet\\Documents\\Code\\decomposition","test_instances\\test_CG_$i.txt")
    path = joinpath(path_dir,name)
    push!(instance_list, structures_decomp.read_instance_file(path))
end


# Compute the avg gap of best known solutions for each config
gap = Dict{instance,Float64}()
optimality = Dict{instance,Bool}()
primal_bound = Dict{instance,Float64}()
dual_bound = Dict{instance,Float64}()
dual_bound_root = Dict{instance,Float64}()
node_processed = Dict{instance,Float64}()
time_run = Dict{instance,Float64}()
cut_generated = Dict{instance,Float64}()
# for each instance perform the tests and register the solution
for i=1:length(instance_list)
    ins = instance_list[i]
    time_start = time()
    tree = SLAP_main_module.main_BP_function(joinpath(path_dir,file_name_list[i]),main_config)
    time_run[ins] = time() - time_start
    optimality[ins] = isempty(tree.nodes)
    primal_bound[ins] = tree.current_primal_bound
    dual_bound[ins] = tree.current_dual_bound
    gap[ins] = (primal_bound[ins] - dual_bound[ins]) / primal_bound[ins]
    node_processed[ins] = length(tree.processed)
    dual_bound_root[ins] = (isempty(tree.processed) ? 0 : tree.processed[1].dual_bound)
    cut_generated[ins] = sum(length(tree.auxiliary_data["cut_pool"]["OG_cut"][s]) for s in keys(tree.auxiliary_data["cut_pool"]["OG_cut"]))
end

# Save results
df = DataFrame(
    name = file_name_list,
    optimality = [optimality[ins] for ins in instance_list],
    time = [time_run[ins] for ins in instance_list],
    gap = [gap[ins] for ins in instance_list],
    tree_size = [node_processed[ins] for ins in instance_list],
    cut_generated = [cut_generated[ins] for ins in instance_list],
    root_dual_bound = [dual_bound_root[ins] for ins in instance_list],
    primal_bound = [primal_bound[ins] for ins in instance_list],
    dual_bound = [dual_bound[ins] for ins in instance_list]
)

CSV.write("test_bcp_baseline$(ARGS[1]).csv", df)

# Here to do/undo comments
=#
end #end function run

end #end module

launch.run() 

#=
structures_decomp.save_instance(structures_decomp.createInstance(n,amax,cmax,hmax,wa,wc,nOrd,nPick,zipfParameter,fixed),"instance_test_complex_pricing.txt")

n = Int16(60)
amax = Int8(5)
cmax = Int8(6)
hmax = Int8(1)
wa = Int16(2)
wc = Int16(1)
nOrd = Int16(50)
nPick = Int16(500)
zipfParameter = 0.737
fixed = 0.

ins = structures_decomp.read_instance_file(joinpath(pwd(),"instance_TG1_50_20_0.5_1.txt"))=#
