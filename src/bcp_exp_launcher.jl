########################################################
# Preambule
########################################################

# The first argument passed to this file is the dir to the instances
# the second arguments is the type of instances
# the third is the routing policy
# the forth is the number of the instance to be read


########################################################
# Main
########################################################

include("structures/structures_decomp.jl")
include("SLAPRP.jl")

module launch

using ..structures_decomp, ..SLAP_main_module

using CSV,DataFrames

function run()

    ###################################
    # Parameters
    ###################################
    # Parameters
    time_limit = 7200
    cpu_limit = Threads.nthreads()
    routing_policy = "return" #ARGS[3]
    instance_type = "silva" #ARGS[2]
    
    EPS = 0.001
    max_iter_cut =  1
    EPS_cut_separation_initial = 0.1
    EPS_cut_violated_initial = 0.1
    EPS_cut_vio_SL1 = 0.025
    max_iter_SL1 = 1
    EPS_delta_separation = 0.1 # unused
    EPS_delta_violated = 0.1 # unused
    max_round_OG = 6
    max_cut_per_order = 500

    EPS_cut_obj_improvement = 0.01
    min_nodes_with_cuts = 5 # the first xxx nodes searched will have cuts generated
    max_depth_cut_generation = 3

    # cut intiation
    activate_all_link_cuts = false
    initiate_aisle_cuts_in_pool = true
    

    # branching
    branching_strategy = "aisle" # can be variable, sum_gub, sum_gl or symmetry
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
    nb_candidates = 20
    bound_method = "basic"

    partial_heuristic_noise = 0.1

    #search strategy
    search_strategy = "priority_queue" # can be priority_queue, lifo or fifo

    # fake route in the model
    fake_route = true

    # routing policy
    bi_directional = false

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
        activate_all_link_cuts,
        initiate_aisle_cuts_in_pool,
        simplex_algo,
        max_cut_per_order
    )

    # read the instance
    #path_dir = joinpath(pwd(), ARGS[1])
    #file_name_list = readdir(path_dir)

    # first launch on the first instance for compilation
    tree_test = SLAP_main_module.main_BP_function(joinpath(pwd(),"instances/silva_small/SLAP-PRP_A1_B5_O1_I3_v1.txt"),main_config)

    for id in [parse(Int,ARGS[4])]

        path = joinpath(path_dir,file_name_list[id])
        if instance_type == "silva"
            ins = structures_decomp.read_silva_instance_file(path,routing_policy)
        elseif instance_type == "classic"
            ins = structures_decomp.read_instance_file(path,routing_policy)
        end

        println(path)

        # run the model
        time_start = time()
        tree = SLAP_main_module.main_BP_function(path,main_config)
        run_time = time() - time_start
        run_time = round(run_time, digits=1)
        optimality = isempty(tree.nodes)
        primal_bound = tree.current_primal_bound
        dual_bound = 2 * ceil((tree.current_dual_bound - 0.0001) / 2)
        node_processed = length(tree.processed)
        dual_bound_root = (isempty(tree.processed) ? 0 : tree.processed[1].dual_bound)
        cut_generated = sum(length(tree.auxiliary_data["cut_pool"]["OG_cut"][s]) for s in keys(tree.auxiliary_data["cut_pool"]["OG_cut"]))

        # register the solution
        ins = tree.auxiliary_data["instance"]
        solution = []
        sol_dict = tree.best_primal_solution
        for l in ins.locList
            for s in ins.skuList
                if sol_dict[s] == l
                    push!(solution,s.id)
                end
            end
        end

        # register the results
        file_path = joinpath(pwd(),joinpath("results","result_$(id).csv"))
        file_name = file_name_list[id]
        open(file_path, "w") do f
            write(f, "Instance,Aisles,Rows,Orders,OrderLines,Slots,Picks,Skewness,Version,UB,LB,Gap,Time,Optimality,Nodes,Root_bound,Cuts,Policy,Solution\n")
            write(f,file_name)
            write(f,",$(ins.amax),$(ins.cmax),$(ins.nOrd),$(Int(ins.nPick/ins.nOrd))")
            write(f, ",$(2*ins.amax*ins.cmax),$(ins.nPick),$(ins.zipfParameter)")
            write(f, ",$(file_name[end-4])") # version
            gap = optimality ? 0 : round((primal_bound - dual_bound)*100/primal_bound,digits=2)
            write(f, ",$(round(primal_bound, digits=2)),$(round(dual_bound, digits=2)),$(gap)")
            write(f, ",$(run_time),$(optimality),$(node_processed),$(round(dual_bound_root, digits=2)),$(cut_generated),$(routing_policy),")
            for s in solution
                write(f, "$s/")
            end
        end
    end
end



end

launch.run()