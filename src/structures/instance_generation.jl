include("structures_decomp.jl")

module launch

using Random, StatsBase

using ..structures_decomp

function run()
weighta = 80
weightb = 15
weightc = 5
for a in [1,3,5]
    for b in [5,10]
        for o in [5,10,20,30,50]
            for qo in [3,5]
                for zipf in [0]
                    for v in [1,2,3]
                        ins = structures_decomp.createInstance(Int16(2*a*b),Int8(a),Int8(b),Int8(1),Int16(1),Int16(1),Int16(o),Int16(qo),0.,"return",weighta,weightb,weightc)
                        structures_decomp.save_instance_silva(ins,"SLAPRP_classic_A$(a)_B$(b)_O$(o)_I$(qo)_Z$(weighta)_v$(v).txt")
                    end
                end
            end
        end
    end
end
#=
path_dir = joinpath(pwd(), "regular_test")
file_name_list = readdir(path_dir)
for id in eachindex(file_name_list)
    path = joinpath(path_dir,file_name_list[id])
    ins = structures_decomp.read_silva_instance_file(path,"return")
    @assert isempty(ins.fixed_positions)

    alpha_list = [0.1,0.2,0.3]
    for alpha in alpha_list
        new_ins = deepcopy(ins)
        fixed_sku_list = sample(ins.skuList, Int(floor((1-alpha)*length(ins.skuList))), replace=false)
        loc_list_sample = [l.id for l in ins.locList]
        append!(loc_list_sample,loc_list_sample)
        fixed_loc_id = sample(loc_list_sample, Int(floor((1-alpha)*length(ins.skuList))), replace=false)
        shuffle!(fixed_loc_id)
        for s in eachindex(fixed_sku_list)
            new_ins.fixed_positions[new_ins.skuList[fixed_sku_list[s].id]] = new_ins.locList[fixed_loc_id[s]]
        end
        # now we save the instance
        name = "$(file_name_list[id][1:end-7])_a$(alpha)$(file_name_list[id][end-6:end])"
        structures_decomp.save_instance_silva(new_ins, name)
    end
end=#


end

end

launch.run()