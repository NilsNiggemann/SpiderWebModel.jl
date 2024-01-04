using OrderedCollections
const P1 =  SA[
    -1.0  1.0   1.0;
    -1.0  0.0  -1.0;
    1.0  1.0  -1.0]
const P2 = -P1

function generateAllFluctuations(StartConfig,operator = findOps(StartConfig)[1];maxiter = 500)

    Configs = [[copy(StartConfig)]]
    uniqeConfs = Set([copy(StartConfig)])
    for iter in 1:maxiter-1
        Confs = Configs[iter]
        newconf = empty(Confs)
        for Conf in Confs
            plaqs = getApplicablePlaquettes(Conf,operator)
            # isempty(plaqs) && (@warn "no fluctuations possible"; break)
            
            for p in plaqs
                Conf2 = copy(Conf)
                pij = getPlaquette(Conf2,p...)
                pij .+= operator
                if Conf2 ∉ uniqeConfs
                    push!(newconf,Conf2)
                    push!(uniqeConfs,Conf2)
                end
            end
            if isempty(newconf) 
                @info "" length(uniqeConfs)
                return Configs,uniqeConfs
            end
            # @assert i1 == i2 "i1 = $i1, i2 = $i2"
        end
        # @info "" length(newconf) length(uniqeConfs)
        push!(Configs,newconf)
    end
    # filter!(x -> fulFillsConstraint(x,verbose = false),Configs)
    # @assert all(fulFillsConstraint.(Configs,verbose = true))
    @warn "max iterations reached" length(uniqeConfs)
    return Configs,uniqeConfs
end

function getStateNum!(AllConfigs::Dict{T,Int},Conf::T) where T
    num = get(AllConfigs,Conf,0)
    if num == 0
        push!(AllConfigs,Conf => length(AllConfigs)+1)
        return length(AllConfigs)+1
    end
    return num
end

function getNeighborStates!(AllConfigs,StartConfig,operator)

    NeighborStates = Int[]

    plaqs = getApplicablePlaquettes(StartConfig,operator)
    
    for p in plaqs
        Conf2 = copy(StartConfig)
        pij = getPlaquette(Conf2,p...)
        pij .+= operator
        
        num = getStateNum!(AllConfigs,Conf2)
        push!(NeighborStates,num)
    end
    
    return NeighborStates
end

function getAllNeighborStates(StartConfig)
    AllConfigs = Dict(StartConfig => 1)
    
    Nplus = Dict(1 => Int[])
    Nminus = Dict(1 => Int[])

    for (Conf,num) in AllConfigs
        neighbors = getNeighborStates!(AllConfigs,Conf,P1)
        Nplus[num] = neighbors
        
        neighbors = getNeighborStates!(AllConfigs,Conf,P2)
        Nminus[num] = neighbors
        println(num)
    end
    @info "" length(AllConfigs) length(Nplus) length(Nminus) 
    
    return (;AllConfigs,Nplus,Nminus)

    # return (;
    # AllConfigs = extractVectors(AllConfigs),
    # Nplus = extractVectors(Nplus),
    # Nminus = extractVectors(Nminus),
    # )
end

function extractVectors(ConfDict)
    AllConfigsVec = collect(keys(ConfDict))
    vals = collect(values(ConfDict))
    return AllConfigsVec[sortperm(vals)]
end

function invertDict(D)
    D2 = Dict(v => k for (k,v) in D)
end