
function initialize_GFMC(L,Spin = 1)
    Config = stencilConfig(zeros(Int8, L, L), Int8(Spin)) 
end


getOperatorRep(i,j,opNum) = i,j,opNum
# getOperatorRep(i,j,opNum) = CartesianIndex(i,j,opNum)
# getOperatorNumber(i,j,opNum,L) = LinearIndices((L,L,2))[i,j,opNum]

# function getOperatorFromNumber(Config, opNumber, operators)

# end

function getMoves!(
    moves,
    Conf
)
    empty!(moves)

    for I in plaquetteIterator(Conf)
        applPlus, applMinus = P_applicable(Conf, I)
        i,j = I
        if applMinus
            push!(moves, getOperatorRep(i,j,-1))
        end
        if applPlus
            push!(moves, getOperatorRep(i,j,1))
        end
    end
    return moves
end



import StatsBase
function performMarkovStep!(weights,moves,Config,weightfunc::T,Λ) where T
    empty!(weights)
    # moves = getMoves!(moves,Config)
    weights = getWeightList!(weights,moves,Config,weightfunc,Λ)

    # TotalWeight = sum(weights)
    if isempty(weights)
        @info "No moves available" 
    end
    moveidx = StatsBase.sample(StatsBase.Weights(weights))
    move = moves[moveidx]
    if move != (0,0,0)
        applyPlaquette!(Config, move[1], move[2], move[3]) 
    end
    return move,weights
end

# function findAffectedPlaquettes!(plaqs,Config,i,j)
#     for I in plaquetteIterator(Config)
#         if !plaquettesAreSeparated(I,(i,j))
#             push!(plaqs,CartesianIndex(I))
#         end
#     end
#     return plaqs
# end

function getWeightList!(weights,moves,Config,weightfunc::T,Λ) where T
    empty!(weights)
    moves = getMoves!(moves,Config)
    ψₓ = weightfunc(Config)

    for operator in moves
        i,j,opNum = operator
        applyPlaquette!(Config, i, j, opNum)
        ψₓ´ = weightfunc(Config)
        push!(weights,ψₓ´/ψₓ)
        applyPlaquette!(Config, i, j, -opNum)
    end
    if Λ != 0
        push!(moves, (0,0,0))
        push!(weights,Λ)
    end
    return weights
end

function getLocalEnergy(weights,Λ)
    return -sum(weights) + Λ
end


function NPlaquettes(Conf)
    moves = 0
    for I in plaquetteIterator(Conf)
        @inbounds applPlus, applMinus = P_applicable(Conf, I)
        moves += applPlus + applMinus
    end
    return moves
end

function varitationalFunc(α,NPlaq::Integer,NPlaqEstimate)
    return exp(α*(NPlaq-NPlaqEstimate))
end

function ConstructVaritationalFunc(α,ConfEx=nothing)
    NPlaqEst = 0
    if ConfEx !== nothing
        NPlaqEst = NPlaquettes(ConfEx)
    end
    ψ = let α = α, NPlaqEst = NPlaqEst
        Conf -> varitationalFunc(α,NPlaquettes(Conf),NPlaqEst)
    end
    # return Conf -> varitationalFunc(α,NPlaquettes(Conf),NPlaqEst)
end

function startSingleWalkerGFMC(InitialState,NSteps,weightfunc::T,Λ) where T
    # Config = initialize_GFMC(InitialState)
    # Config = copy(InitialState)
    Config = copy(InitialState)
    weights = Float64[]
    Allmoves = fill((Int8(0),Int8(0),Int8(0)),NSteps)
    LocalMoves = empty(Allmoves)
    energies = zeros(NSteps)
    TotalWeights = zeros(NSteps)
    
    for i in 1:NSteps
        move,weights = performMarkovStep!(weights,LocalMoves,Config,weightfunc,Λ)
        
        bx_TotalWeight = sum(weights)
        
        TotalWeights[i] = bx_TotalWeight
        
        Allmoves[i] = move
        getWeightList!(weights,LocalMoves,Config,weightfunc,Λ)
        energies[i] = getLocalEnergy(weights,Λ)
        
    end
    return (;Allmoves,energies,TotalWeights)
end

function normalizedAccWeight(weights,n,p)
    meanweight = mean(weights)
    # meanweight = 1
    prod(weights[n-j]/meanweight for j in 1:p)
end

function precomputeNormalizedAccWeight(weights,nThermal,PMax)
    bn = @view weights[nThermal:end]
    meanweight = mean(bn)
    # meanweight = mean(weights)
    
    Gnp = zeros(length(bn),PMax)
    for n in axes(Gnp,1)
        Gnp[n,1] = bn[n]/meanweight 
    end
    for p in 2:PMax
        for n in p:length(bn)
            # Gnp[n,p] = Gnp[n,p-1]*bn[n-p]/meanweight
            Gnp[n,p] = Gnp[n-1,p-1]*Gnp[n,1]
        end
    end
    return Gnp
end


function iterateProjector!(Gnp_new,Gnp,Gn1,p)
    
    for n in axes(Gnp,1)[p:end]
        Gnp_new[n] = Gnp[n]*Gn1[n-p+1]
    end
    return Gnp_new
end


function getEnergy(weights,localEnergies,p,nthermalization)
    meanweight = mean(weights)
    num = 0.
    denom = 0.
    
    for n in nthermalization:length(localEnergies)
        Gnp = 1.
        for j in 1:p
            Gnp *= weights[n-j]/meanweight
        end
        num += Gnp*localEnergies[n]
        denom += Gnp
    end
    return num/denom
end

function getEnergies(weights,localEnergies,nthermalization,PMax;
    Gnp = precomputeNormalizedAccWeight(weights,nthermalization,PMax)
    )
    
    EL_thermalized = @view localEnergies[nthermalization:end]
    N = lastindex(EL_thermalized)
    num = zeros(PMax)
    denom = zeros(PMax)
    for p in 1:PMax
        for n in p:N
            num[p] += Gnp[n,p]*EL_thermalized[n]
            denom[p] += Gnp[n,p]
        end
    end
    # num = fetch.([Threads.@spawn sum(Gnp[n+p,p]*EL_thermalized[n] for n in axes(Gnp,1)[begin:end-p]) for p in 1:PMax])
    # denom = fetch.([Threads.@spawn sum(Gnp[n+p,p] for n in axes(Gnp,1)[begin:end-p]) for p in 1:PMax])
    return num ./denom
end


function setupProjector(weights,nThermal)
    bn = @view weights[nThermal:end]
    meanweight = mean(bn)
    Gn1 = bn./meanweight
end


function computeGnpForward(weights,p,N=0)
    meanweight = mean(weights)
    # meanweight = mean(weights)
    
    Gnp = zeros(length(weights),p+N)

    for n in axes(Gnp,1)
        Gnp[n,1] = weights[n]/meanweight 
    end
    for p in 2:PMax
        for n in p:length(weights)
            # Gnp[n,p] = Gnp[n,p-1]*weights[n-p]/meanweight
            Gnp[n,p] = Gnp[n-1,p-1]*Gnp[n+N,1]*Gnp[n-p+1,1]
        end
    end
    return Gnp
end


function getObs(Gnp,AllConfigs,reconfigurationTable,ObsFunc,m=size(Gnp,2)÷2)
    N = lastindex(AllConfigs)
    Obs = ObsFunc(AllConfigs[1][1])
    num = zero(typeof(Obs))
    denom = zero(typeof(Obs))
    # fill!(Obs,zero(eltype(Obs)))

    Nw = size(reconfigurationTable,1)
    p = size(Gnp,2)
    for n in m+1:N
        Gn = Gnp[n,p]
        denom += Gn*Nw
        for α in 1:Nw
            α´ = α
            for i_m in 1:m
                α´ = reconfigurationTable[α´,n-i_m]
            end
            num += Gn*ObsFunc(AllConfigs[n-m][α´])
        end
    end
    return num/denom
end


function getObservablesOld(result,StartConf,ObsFunc,nthermalization,PMax)
    
    weights = result.TotalWeights
    localEnergies = result.energies
    
    Gn1 = setupProjector(weights,nthermalization)
    Gnp = zero(Gn1)
    Gnp_new = zero(Gn1)
    
    EL_thermalized = @view localEnergies[nthermalization:end]
    
    Energy_num = zeros(PMax)
    Denom = zeros(PMax)
    
    Gnp .= Gn1
    
    Energy_num[1] = Gn1'*EL_thermalized
    Denom[1] = sum(Gn1)


    for p in 2:PMax
        iterateProjector!(Gnp_new,Gnp,Gn1,p)
        Gnp .= Gnp_new

        Energy_num[p] = Gnp' * EL_thermalized
        Denom[p] = sum(Gnp)
    end

    Conf = copy(StartConf)
    moves = result.Allmoves
    for i in 1:nthermalization
        i,j,type = moves[i]
        type != 0 && applyPlaquette!(Conf, i,j,type)
    end

    Newmoves = @view moves[nthermalization+1:end]
    
    Obs_num = zero(ObsFunc(Conf))

    for n in eachindex(Gnp)[PMax÷2+1:end-PMax]
        i,j,type = Newmoves[n]
        type != 0 && applyPlaquette!(Conf, i,j,type)

        Obs_num .+= Gnp[n+PMax]*ObsFunc(Conf)
    end
    
    E0 = Energy_num ./Denom
    Obs = Obs_num /Denom[end]
    
    return (;E0,Obs)
end

function splitIntoBins(array,binsize)
    Iterators.partition(array,binsize)
end

function startManyWalkerGFMC(InitialState,Nwalkers,NSteps,nBranch,weightfunc::T,Λ) where T
    Walkers = [copy(InitialState) for _ in 1:Nwalkers]

    bx_List = [Float64[] for _ in 1:Nwalkers] # List of weights for each walker at each step
    
    weights = ones(Nwalkers)
    SaveConfigs = Vector{Vector{typeof(InitialState)}}()

    LocalMoves = [ Vector{Tuple{Int8,Int8,Int8}}() for _ in 1:Nwalkers]
    energies = zeros(NSteps)
    TotalWeights = zeros(NSteps)
    reconfigurationList = zeros(Int,Nwalkers)

    reconfTable = zeros(Int,Nwalkers,NSteps)

    for i in 1:NSteps
        # for (α,Config) in enumerate(Walkers)
         for α in eachindex(Walkers)
            Config = Walkers[α]
            moveList = LocalMoves[α]
            weightList = bx_List[α]
            
            # weights[α] = one(weights[α])
            w = 1
            for step in 1:nBranch
                move,weightList = performMarkovStep!(weightList,moveList,Config,weightfunc,Λ)
                bx = sum(weightList)/size(Config,1)
                w *= bx
            end
            weights[α] = w
            getWeightList!(weightList,moveList,Config,weightfunc,Λ)
        end
        energies[i] = getLocalEnergyWalkers_before(weights,bx_List,Λ)

        TotalWeights[i] = mean(weights)
        reconfiguration!(Walkers,reconfigurationList,weights)
        reconfTable[:,i] .= reconfigurationList

        SaveWalkers = empty(Walkers)
        # for (α,α´) in enumerate(reconfigurationList)
        #     if α == α´ # convention: surviving walkers always remain a copy at their index
        #         push!(SaveWalkers,copy(Walkers[α]))
        #     end
        # end

        for Walker in Walkers
            push!(SaveWalkers,copy(Walker))
        end

        push!(SaveConfigs,SaveWalkers)
    end
    return (;TotalWeights, energies, SaveConfigs,reconfTable)
end

function reconfiguration!(Walkers,reconfigurationList,weights)
    StatsBase.sample!(eachindex(Walkers),StatsBase.Weights(weights),reconfigurationList,ordered = true)
    minimizeReconfiguration!(reconfigurationList)
    for (α,α´) in enumerate(reconfigurationList)
        if α´ != α
            Walkers[α] .= Walkers[α´]
        end
    end
end

function getLocalEnergyWalkers_before(weights,bx_List,Λ)
    Nw = length(weights)
    num = 0.
    denom = 0.

    for α in eachindex(weights,bx_List)
        eloc = Λ - sum(bx_List[α])
        num += weights[α]*eloc
        denom += weights[α]
    end

    return num/denom
end

function getLocalEnergyWalkers_after(weights,bx_List,Λ)
    Nw = length(weights)
    num = 0.

    for α in eachindex(weights)
        num += getLocalEnergy(bx_List[α],Λ)
    end
    return num/Nw
end

"""given a sorted list reconfiguration indices, minimizes the number of reconfigurations by swapping elements in the list. Each walker that survives a reconfiguration step remains unchanged while walkers that are killed get assigned to a new index."""
function minimizeReconfiguration!(list)
    for i in eachindex(list)
        jRange = searchsorted(list,i)
        if length(jRange) != 0
            li = list[i]
            list[i] = i
            list[jRange[1]] = li
        end
    end
    return list
end
