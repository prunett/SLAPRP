module structures_decomp

using Random, Distributions, StatsBase

export sku, order, loc, ioloc, route, solution, instance, intersects, intersects_id

struct sku 
    id::Int32
    demand::Int16 #demand
end

Base.:(==)(a::sku,b::sku) = (a.id == b.id) 
Base.isequal(a::sku,b::sku) = isequal(a.id, b.id)
Base.hash(a::sku, h::UInt) = hash(a.id, hash(:sku,h))

struct order 
    id::Int32
    sku::Array{sku,1} #list of sku
end

Base.:(==)(a::order,b::order) = (a.id == b.id) 
Base.isequal(a::order,b::order) = isequal(a.id, b.id)
Base.hash(a::order, h::UInt) = hash(a.id, hash(:order,h))

struct loc
    id::Int8
    aisle::Int8
    col::Int16
    height::Int8 #height
    side::Int8
end

Base.:(==)(a::loc,b::loc) = (a.id == b.id) 
Base.isequal(a::loc,b::loc) = isequal(a.id, b.id)
Base.hash(a::loc, h::UInt) = hash(a.id, hash(:loc,h))

ioloc = loc(0,1,0,0,0)

struct route
    order::order # corresponding order
    loc::Array{loc,1}
    cost::Int
end


struct solution
    cost::Int
    routes::Array{route,1}
    stol::Dict{sku,loc}
end

struct sku_symmetry
    id::Int
    skus::Vector{sku}
    size::Int
end

struct instance
    n::Int16 #number of sku
    #m::Int16 # number of storage locations
    amax::Int8 #number of aisles
    cmax::Int8 #number of columns in one aisle
    hmax::Int8
    wa::Int16
    wc::Int16
    nOrd::Int16 #number of orders
    nPick::Int16
    zipfParameter::Union{Float64,Nothing} #parameter for the pareto distribution
    skuList::Array{sku,1}
    orderList::Array{order,1}
    orderDict::Dict{Int32,order}
    locList::Array{loc,1}
    sku_orderList::Dict{sku,Array{order,1}}
    n_sym::Int8 # number of equivalent positions
    fixed_positions::Dict{sku,loc}
    policy::String
    #sku_from_code_to_file::Dict{sku,Int} # this is the correspondance between the sku from the instance and the ones used
    # in the code, the keys are the sku from the code, and the values the sku from file (the id)
    #sku_from_file_to_code::Dict{Int,Union{sku,Nothing}} # this is the opposite
    sku_symmetry_list::Vector{sku_symmetry}
    number_blocks_layout::Int
end



Base.:(==)(a::sku_symmetry,b::sku_symmetry) = (a.id == b.id) 
Base.isequal(a::sku_symmetry,b::sku_symmetry) = isequal(a.id, b.id)
Base.hash(a::sku_symmetry, h::UInt) = hash(a.id, hash(:sku_symmetry,h))





function createInstance(n::Int16,amax::Int8,cmax::Int8,hmax::Int8,wa::Int16,wc::Int16,nOrd::Int16,pick_per_order::Int16,fixed::Float64,policy::String,weighta::Int,weightb::Int,weightc::Int)
    #This function create a random instance according to the specified parameters
    ### !!!!!!!!!!! Problem with Pareto distribution
    # Basic definition for e, simple layout

    skuIdList = Int32[i for i=1:n]
    orderIdList = Int32[i for i=1:nOrd]

    #First step: number of items per order
    orderNumber = fill(1,nOrd)
    orderPick = nOrd
    while orderPick < nOrd*pick_per_order
        new_pick_order_id = rand(1:nOrd)
        if orderNumber[new_pick_order_id] < n
            orderNumber[new_pick_order_id] += 1
            orderPick += 1
        end
    end

    # New method for order number
    #orderNumber = fill(pick_per_order,nOrd)

    #skuWeightDist = Pareto(pareto) ########### Maybe to change, problem with pareto distribution
    #skuWeight = Weights(rand(skuWeightDist, n))

    #create the demand array for each order
    skuDemand = Int16[0 for i=1:n]
    weight_list = [0 for k=1:n]
    ab_limit = floor(Int,n/3)
    bc_limit = floor(Int,2*n/3)
    for k=1:ab_limit
        weight_list[k] = weighta
    end
    for k=ab_limit+1:bc_limit
        weight_list[k] = weightb
    end
    for k=bc_limit+1:n
        weight_list[k] = weightc
    end

    skuWeight = Weights(weight_list)

    #for each order generate the list of SKU ID to include
    orderSkuIdList = []
    for o=1:nOrd
        i = 1
        currentSkuIdList = sample(1:n, skuWeight, orderNumber[o], replace = false)
        for i in currentSkuIdList
            skuDemand[i] += 1
        end

        push!(orderSkuIdList, deepcopy(currentSkuIdList))
    end

    #create the list of sku with demand
    skuList = sku[]
    for i=1:n
        push!(skuList, sku(skuIdList[i], skuDemand[i]))
    end

    #create the list of order
    orderList = order[]
    for o=1:nOrd
        listOfSku = sku[skuList[orderSkuIdList[o][i]] for i=1:length(orderSkuIdList[o])]
        push!(orderList, order(orderIdList[o], deepcopy(listOfSku)))
    end
    orderDict = Dict{Int32,order}(o.id => o for o in orderList)
    #sku order list
    sku_orderList = Dict{sku, Array{order,1}}(s => deepcopy(filter(o -> s in o.sku, orderList)) for s in skuList)

    #locations list
    locList = loc[]
    idloc = 1
    for a=1:amax
        for c=1:cmax
            for h=1:hmax
                push!(locList,loc(idloc,a,c,h,1))
                idloc += 1
            end
        end
    end

    # fixed_positions setup
    fixed_positions = Dict{sku,loc}()
    available_loc = loc[]
    for l in locList
        for i=1:2*hmax
            push!(available_loc,l)
        end
    end
    # fix the positions
    if fixed > 0
        sku_fixed = sample(skuList,floor(Int,fixed*length(skuList)),replace=false)
        loc_fixed = sample(available_loc,floor(Int,fixed*length(skuList)),replace=false)
        shuffle!(loc_fixed)
        for i=1:length(sku_fixed)
            fixed_positions[sku_fixed[i]] = loc_fixed[i]
        end
    end

    return instance(n,amax,cmax,hmax,wa,wc,nOrd,nOrd*pick_per_order,weighta,skuList,orderList,orderDict,locList,sku_orderList,2*hmax,fixed_positions,"exact",Vector{sku_symmetry}(),1)
end

createInstance(n::Int16,amax::Int8,cmax::Int8,hmax::Int8,wa::Int16,wc::Int16,nOrd::Int16,nPick::Int16,zipfParameter::Float64) = createInstance(n,amax,cmax,hmax,wa,wc,nOrd,nPick,zipfParameter,0.,"exact")

function save_instance(ins::instance,name::String)
#This function save a given instance to a text file with the given name
#Need to be checked for ergo load
    #save it to file
    filePath = joinpath(pwd(),name)
    open(filePath,"w") do f
        write(f,"///// Instance $name /////\n\n")
        write(f,"/// Layout a_max   c_max   h_max\n")
        write(f,"$(ins.amax)   $(ins.cmax)   $(ins.hmax)\n")
        write(f,"Distance w_a   w_c\n")
        write(f,"$(ins.wa)   $(ins.wc)\n\n")
        write(f,"/// Number of SKU n = \n$(ins.n)\n\n")
        write(f,"/// Number of orders nOrd = \n$(ins.nOrd)\n\n")
        write(f,"/// Number of picks =\n$(ins.nPick)\n\n")
        write(f,"/// Number of fixed positions = \n$(length(keys(ins.fixed_positions)))\n\n")
        write(f,"/// Parameter of the zipf distribution used for sampling =\n$(ins.zipfParameter)\n\n")
        write(f,"/// Different ergonomic load depending on height\nh = ")
        for i=1:ins.hmax
            write(f,"$i   ")
        end
        write(f,"\n")
        #=for i=1:ins.hmax
            write(f,"$(ins.locList[i].e)   ")
        end=#
        write(f,"1   0")
        write(f,"\n\n/// List of SKU\nID   Demand\n")
        for i=1:ins.n
            write(f,"$(ins.skuList[i].id)   $(ins.skuList[i].demand)\n")
        end
        write(f,"\n\n/// List of Orders\nID   Number of order lines   (List of SKU ID)\n")
        for i=1:ins.nOrd
            write(f,"$(ins.orderList[i].id)   $(length(ins.orderList[i].sku))   ")
            for j=1:length(ins.orderList[i].sku)
                write(f,"$(ins.orderList[i].sku[j].id) ")
            end
            write(f,"\n")
        end
        # fixed positions
        write(f,"\n /// List of fixed positions\n SKU ID => loc ID\n")
        for s in keys(ins.fixed_positions)
            write(f,"$(s.id)   $(ins.fixed_positions[s].id)\n")
        end

    end
end

function save_instance_silva(ins::instance,name::String)
    # this function save the instance in silva format
    filePath = joinpath(pwd(),name)
    open(filePath,"w") do f
        write(f, "$(ins.amax) $(ins.cmax)\n")
        write(f, "$(ins.wa) $(ins.wc)\n")
        write(f, "$(ins.amax*ins.cmax*2)\n")
        write(f, "$(ins.nOrd)\n")
        for o in ins.orderList
            write(f, "$(length(o.sku)) ")
        end
        write(f, "\n")
        for o in ins.orderList
            for s in o.sku
                write(f, "$(s.id) ")
            end
            write(f,"\n")
        end
        for s in keys(ins.fixed_positions)
            write(f, "$(s.id) $(ins.fixed_positions[s].id)\n")
        end
    end
end

function read_instance_file(path::String,policy::String)
    # This function retrieve the data from an instance file to create an instanceFirst type
    
    amax = Int8(0)
    cmax = Int8(0)
    hmax = Int8(0)
    wa = Int16(0)
    wc = Int16(0)
    n = Int16(0)
    nOrd = Int16(0)
    nPick = Int16(0)
    zipf = 0.
    skuList = Array{sku,1}(undef,0)
    orderList = Array{order,1}(undef,0)
    locList = Array{loc,1}(undef,0)
    fixed_positions = Dict{sku,loc}()
    #sku_from_code_to_file = Dict{sku,Int}()
    #sku_from_file_to_code = Dict{Int,Union{sku,Nothing}}()

    open(path) do f
        readline(f)
        readline(f)
        readline(f)
        temp = split(readline(f))
        amax = parse(Int8,temp[1])
        cmax = parse(Int8,temp[2])
        hmax = parse(Int8,temp[3])
        readline(f)
        wtemp = split(readline(f))
        wa = parse(Int16,wtemp[1])
        wc = parse(Int16,wtemp[2])
        readline(f)
        readline(f)
        n_file = parse(Int16,readline(f))
        readline(f)
        readline(f)
        nOrd = parse(Int16,readline(f))
        readline(f)
        readline(f)
        nPick = parse(Int16,readline(f))
        readline(f)
        readline(f)
        nfixed = parse(Int16,readline(f))
        readline(f)
        readline(f)
        zipf = parse(Float64,readline(f))
        readline(f)
        readline(f)
        readline(f)
        ergString = split(readline(f))
        erg = Int16[parse(Int16,ergString[i]) for i=1:length(ergString)]
        readline(f)
        readline(f)
        readline(f)

        #here we begin to read SKU
        skuList = sku[]
        for i=1:n_file
            sktemp = split(readline(f))
            id_temp = parse(Int32,sktemp[1])
            demand = parse(Int16,sktemp[2])
            #=if demand == 0 # case where we do not keep this dummy sku
                sku_from_file_to_code[id_temp] = nothing
            else
                push!(skuList, sku(length(skuList)+1,demand))
                sku_from_code_to_file[skuList[end]] = id_temp
                sku_from_file_to_code[id_temp] = skuList[end]
            end=#
            push!(skuList, sku(length(skuList)+1,demand))
        end
        n = n_file
        readline(f)
        readline(f)
        readline(f)
        readline(f)

        #here we begin to read orders
        skuDic = Dict{Int32,sku}(s.id => s for s in skuList) #Dict for sku
        orderList = Array{order,1}(undef,nOrd)
        for i=1:nOrd
            orderTemp = split(readline(f))
            id = parse(Int32,orderTemp[1])
            curNumb = parse(Int,orderTemp[2])
            curSkuList = Array{sku,1}(undef,curNumb)
            for j=1:curNumb
                #curSkuList[j] = sku_from_file_to_code[parse(Int,orderTemp[j+2])]
                curSkuList[j] = skuList[parse(Int,orderTemp[j+2])]
                @assert !isnothing(curSkuList[j])
            end
            orderList[i] = order(id,copy(curSkuList))
        end

        #creation of locList
        i=1
        for a=1:amax
            for c=1:cmax
                push!(locList, loc(Int8(i),Int8(a),Int16(c),Int8(1),Int8(1)))
                i += 1
            end
        end

        # here we read fixed positions
        readline(f)
        readline(f)
        readline(f)
        for i=1:nfixed
            split_string = split(readline(f))
            sku_id = parse(Int32,split_string[1])
            loc_id = parse(Int32,split_string[2])
            #@assert(sku_id == skuList[sku_id].id)
            @assert(loc_id == locList[loc_id].id)
            fixed_positions[skuList[sku_id]] = locList[loc_id]
        end

    end
    orderDict = Dict{Int32,order}(o.id => o for o in orderList)
    sku_orderList = Dict{sku, Array{order,1}}(s => filter(o -> s in o.sku, orderList) for s in skuList)

    sku_symmetry_list = create_symmetry_list(skuList,orderList)
    
    return instance(n,amax,cmax,hmax,wa,wc,nOrd,nPick,zipf,skuList,orderList,orderDict,
    locList,sku_orderList,2*hmax,fixed_positions,policy,
    sku_symmetry_list,1)
end

#function read_silva_instance_file
function read_silva_instance_file(path::String,policy::String)
    # This function retrieve the data from an instance file to create an instanceFirst type

    amax = Int8(0)
    cmax = Int8(0)
    hmax = Int8(0)
    wa = Int16(0)
    wc = Int16(0)
    n = Int16(0)
    nOrd = Int16(0)
    nPick = Int16(0)
    zipf = 0.
    skuList = Array{sku,1}(undef,0)
    orderList = Array{order,1}(undef,0)
    locList = Array{loc,1}(undef,0)
    fixed_positions = Dict{sku,loc}()
    #sku_from_code_to_file = Dict{sku,Int}()
    #sku_from_file_to_code = Dict{Int,Union{sku,Nothing}}()

    open(path) do f
        temp = split(readline(f))
        amax = parse(Int8,temp[1])
        cmax = parse(Int8,temp[2])
        hmax = 1
        wtemp = split(readline(f))
        wa = parse(Int16,wtemp[1])
        wc = parse(Int16,wtemp[2])
        n_temp = parse(Int16,readline(f))
        nOrd = parse(Int16,readline(f))
        pick_per_order_list = parse.(Int,split(readline(f)))
        #pick_per_order = parse(Int,split(readline(f))[1])
        pick_per_order = pick_per_order_list[1]
        nPick = nOrd * pick_per_order
        nfixed = 0
        zipf = nothing

        order_lines = []
        # read the order list
        for i=1:nOrd
            push!(order_lines,parse.(Int,split(readline(f))))
        end

        #here we begin to read SKU
        skuList = sku[]
        for i=1:n_temp
            sku_count = 0
            for u=1:nOrd
                for v in eachindex(order_lines[u])
                    if i == order_lines[u][v]
                        sku_count += 1
                    end
                end
            end

            #=if sku_count == 0
                sku_from_file_to_code[i] = nothing
            else
                push!(skuList, sku(length(skuList)+1, sku_count))
                sku_from_code_to_file[skuList[end]] = i
                sku_from_file_to_code[i] = skuList[end]
            end=#
            push!(skuList, sku(length(skuList)+1, sku_count))
        end
        n = length(skuList)
        

        #here we begin to read orders
        skuDic = Dict{Int32,sku}(s.id => s for s in skuList) #Dict for sku
        orderList = Array{order,1}(undef,nOrd)
        for i=1:nOrd
            id = i
            curNumb = pick_per_order_list[i]
            curSkuList = Array{sku,1}(undef,curNumb)
            for j=1:curNumb
                curSkuList[j] = skuList[order_lines[i][j]]
            end
            orderList[i] = order(id,curSkuList)
        end

        #creation of locList
        i=1
        for a=1:amax
            for c=1:cmax
                push!(locList, loc(Int8(i),Int8(a),Int16(c),Int8(1),Int8(1)))
                i += 1
            end
        end

        # fixed positions
        line = readline(f)
        while !isempty(line)
            fixed_sku = parse(Int,split(line)[1])
            fixed_loc = parse(Int,split(line)[2])
            @assert skuList[fixed_sku].id == fixed_sku
            @assert locList[fixed_loc].id == fixed_loc
            fixed_positions[skuList[fixed_sku]] = locList[fixed_loc]
            line = readline(f)
        end


    end
    orderDict = Dict{Int32,order}(o.id => o for o in orderList)
    sku_orderList = Dict{sku, Array{order,1}}(s => filter(o -> s in o.sku, orderList) for s in skuList)

    sku_symmetry_list = create_symmetry_list(skuList,orderList)
    
    return instance(n,amax,cmax,hmax,wa,wc,nOrd,nPick,zipf,skuList,orderList,orderDict,
    locList,sku_orderList,2*hmax,fixed_positions,policy,
    sku_symmetry_list,1)
end

function read_guo_instance_file(path::String,policy::String)
    # This function retrieve the data from an instance file to create an instanceFirst type

    amax = Int8(0)
    cmax = Int8(0)
    hmax = Int8(0)
    wa = Int16(0)
    wc = Int16(0)
    n = Int16(0)
    nOrd = Int16(0)
    nPick = Int16(0)
    zipf = 0.
    skuList = Array{sku,1}(undef,0)
    orderList = Array{order,1}(undef,0)
    locList = Array{loc,1}(undef,0)
    fixed_positions = Dict{sku,loc}()
    #sku_from_code_to_file = Dict{sku,Int}()
    #sku_from_file_to_code = Dict{Int,Union{sku,Nothing}}()

    open(path) do f
        temp = split(readline(f))
        amax = parse(Int8,temp[1])
        cmax = parse(Int8,temp[2])
        hmax = 1
        wtemp = split(readline(f))
        
        n_temp = parse(Int16,readline(f))
        nOrd = parse(Int16,readline(f))
        pick_per_order = 0
        nPick = 0
        nfixed = 0
        zipf = nothing

        order_lines = []
        # read the order list
        for i=1:nOrd
            push!(order_lines,parse.(Int,split(readline(f))))
        end

        #here we begin to read SKU
        skuList = sku[]
        for i=1:n_temp
            sku_count = 0
            for u=1:nOrd
                for v = 1:pick_per_order
                    if i == order_lines[u][v]
                        sku_count += 1
                    end
                end
            end

            #=if sku_count == 0
                sku_from_file_to_code[i] = nothing
            else
                push!(skuList, sku(length(skuList)+1, sku_count))
                sku_from_code_to_file[skuList[end]] = i
                sku_from_file_to_code[i] = skuList[end]
            end=#
            push!(skuList, sku(length(skuList)+1, sku_count))
        end
        n = length(skuList)
        

        #here we begin to read orders
        skuDic = Dict{Int32,sku}(s.id => s for s in skuList) #Dict for sku
        orderList = Array{order,1}(undef,nOrd)
        for i=1:nOrd
            id = i
            curNumb = pick_per_order
            curSkuList = Array{sku,1}(undef,curNumb)
            for j=1:curNumb
                curSkuList[j] = skuList[order_lines[i][j]]
            end
            orderList[i] = order(id,curSkuList)
        end

        #creation of locList
        i=1
        for a=1:amax
            for c=1:cmax
                push!(locList, loc(Int8(i),Int8(a),Int16(c),Int8(1),Int8(1)))
                i += 1
            end
        end

        # fixed positions
        line = readline(f)
        while !isempty(line)
            fixed_sku = parse(Int,split(line)[1])
            fixed_loc = parse(Int,split(line)[2])
            @assert skuList[fixed_sku].id == fixed_sku
            @assert locList[fixed_loc].id == fixed_loc
            fixed_positions[skuList[fixed_sku]] = locList[fixed_loc]
            line = readline(f)
        end


    end
    orderDict = Dict{Int32,order}(o.id => o for o in orderList)
    sku_orderList = Dict{sku, Array{order,1}}(s => filter(o -> s in o.sku, orderList) for s in skuList)

    sku_symmetry_list = create_symmetry_list(skuList,orderList)
    
    return instance(n,amax,cmax,hmax,0,0,nOrd,nPick,zipf,skuList,orderList,orderDict,
    locList,sku_orderList,2*hmax,fixed_positions,policy,
    sku_symmetry_list,1)
end

function isequivalent(a::sku,b::sku,orderList::Vector{order})
    # this function returns true if a and b are equivalent
    for o in orderList
        if (a in o.sku) && !(b in o.sku)
            return false
        elseif !(a in o.sku) && (b in o.sku)
            return false
        end
    end
    
    return true
end
        
    
function create_symmetry_list(skuList::Vector{sku},orderList::Vector{order})
# This function creates the list of sku symmetries. a sky symmetry is a list of equivalent skus
    symmetry_list = sku_symmetry[]
    sku_available = copy(skuList)
    # we iterate until they are all placed somewhere
    while !isempty(sku_available)
        first_sku = popfirst!(sku_available)
        sku_array = [first_sku]
        size = 1
        sind = 1
        # we look for skus equivalent to first_sku
        while sind <= length(sku_available)
            s = sku_available[sind]
            if isequivalent(first_sku,s,orderList)
                push!(sku_array,s)
                deleteat!(sku_available,sind)
                size += 1
                sind -= 1
            end
            sind += 1
        end
        # now we have a symmetry
        push!(symmetry_list, sku_symmetry(length(symmetry_list)+1,copy(sku_array),size))
    end

    return symmetry_list
end

function intersects(u, v)
    for x in u
        if x in v
            return true
        end
    end
    false
end

#=
function intersects_id(u, v)
    for x in u
        for y in v
            if x.id == y.id
                return true
            end
        end
    end
    false
end=#
function intersects_id(u,v)
    return !isdisjoint(u,v)
end

function aggregate_instance(ins::instance)
# This function aggregates the equivalent locs of one instance
# It returns the new instance file and the number of equivalent positions
    # Computes number of equivalent position, or throw error
    if ins.n_sym == 1
        throw(DomainError("Problem !! Trying to aggregate an instance already aggregated !"))
    else
        n_sym = 2*ins.hmax
    end

    # Duplicate locList and filter unwanted elements
    locList = filter(l -> (l.side == 1) && (l.height == 1), ins.locList)

    return instance(ins.n,ins.amax,ins.cmax,1,ins.wa,ins.wc,ins.nOrd,ins.nPick,ins.zipf,ins.skuList,ins.orderList,ins.orderDict,locList,ins.sku_orderList,n_sym,ins.policy)
end

end #end of module
