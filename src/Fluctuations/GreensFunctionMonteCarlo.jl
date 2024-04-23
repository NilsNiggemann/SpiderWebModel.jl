
abstract type AbstractWalker end

struct SpiderWebWalker{C} <: AbstractWalker
    Config::C
    moves::Vector{Tuple{Int8,Int8,Int8}}
    weights::Vector{Float64}
    indices::Vector{Int}
    n_x::Vector{Int}
    n_x´::Vector{Int}
end
function spiderWebWalker(S)
    Config = copy(S)
    moves = Vector{Tuple{Int8,Int8,Int8}}()
    weights = Vector{Float64}()
    indices = Vector{Int}()
    n_x = zeros(Int,length(collect(plaquetteIterator(Config))))
    n_x´ = Vector{Int}()
    return SpiderWebWalker(Config,moves,weights,indices,n_x,n_x´)
end

get_config(Walker::SpiderWebWalker) = Walker.Config
get_weights(Walker::SpiderWebWalker) = Walker.weights

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

function findAffectedPlaquettes!(Plaq_indices,Config,i,j)
    empty!(Plaq_indices)
    for (index,I) in enumerate(plaquetteIterator(Config))
        if !plaquettesAreSeparated(I,(i,j))
            push!(Plaq_indices,index)
        end
    end
    return Plaq_indices
end

# Assumes that length(n_plaq) == length(plaquetteIterator(Config))
function getNPlaq!(n_plaq, Config)
    for (i,I) in enumerate(plaquetteIterator(Config))
        @inbounds applPlus, applMinus = P_applicable(Config, I)
        n = applPlus + applMinus
        n_plaq[i] = n
    end
    return n_plaq
end

function getNPlaq!(n_plaq, Config,Plaq_indices)
    ind = 1
    # resize!(n_plaq,length(Plaq_indices))
    # println(Plaq_indices)
    empty!(n_plaq)
    for (i,I) in enumerate(plaquetteIterator(Config))
        if i == Plaq_indices[ind]
            # println("\t",i," ",I, " ",ind)
            # i ∉ Plaq_indices && continue
            # ind += 1
            @inbounds applPlus, applMinus = P_applicable(Config, I)
            n = applPlus + applMinus
            push!(n_plaq,n)
            ind += 1
            ind > length(Plaq_indices) && break
        end
    end
    return n_plaq
end

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


function getNPlaq_difference(nPlaq_x,nPlaq_x´,affectedPlaquettes)
    N□ = 0
    for (i,plaqIndex) in enumerate(affectedPlaquettes)
        N□ += nPlaq_x´[i] - nPlaq_x[plaqIndex]
    end
    return N□
end

function getWeightList!(Walker::SpiderWebWalker,weightfunc::T,Λ) where T
    (;Config,moves,weights,indices,n_x,n_x´) = Walker
    empty!(weights)

    moves = getMoves!(moves,Config)

    getNPlaq!(n_x,Config)
    
    for operator in moves
        i,j,opNum = operator
        findAffectedPlaquettes!(indices,Config,i,j)
        applyPlaquette!(Config, i, j, opNum)
        
        getNPlaq!(n_x´,Config,indices)
        
        N□ = getNPlaq_difference(n_x,n_x´,indices) 
        weight = weightfunc(N□)

        push!(weights,weight)
        applyPlaquette!(Config, i, j, -opNum)
    end
    if Λ != 0
        push!(moves, (0,0,0))
        push!(weights,Λ)
    end
    return weights
end

function performMarkovStep!(Walker::SpiderWebWalker,weightfunc::T,Λ) where T
    # moves = getMoves!(moves,Config)
    weights = getWeightList!(Walker,weightfunc,Λ)

    if isempty(weights)
        @info "No moves available" 
    end
    moveidx = StatsBase.sample(StatsBase.Weights(weights))
    move = Walker.moves[moveidx]
    if move != (0,0,0)
        applyPlaquette!(Walker.Config, move[1], move[2], move[3]) 
    end
    return move,weights
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

add_elementwise!(x::AbstractArray,y) = (x .+= y)
add_elementwise!(x::Number,y::Number) = x + y

mult_elementwise!(x::Number,y::Number) = x * y
mult_elementwise!(x::AbstractArray,y::Number) = (x .*= y)

divide_elementwise!(x::Number,y::Number) = x / y
divide_elementwise!(x::AbstractArray,y::AbstractArray) = (x ./= y)

function getObs(Gnp,AllConfigs,reconfigurationTable,ObsFunc,m=size(Gnp,2)÷2)
    N = lastindex(AllConfigs)
    Obs = ObsFunc(AllConfigs[1][1])
    num = zero(Obs)
    denom = zero(Obs)
    # fill!(Obs,zero(eltype(Obs)))

    Nw = size(reconfigurationTable,1)
    p = size(Gnp,2)
    for n in m+1:N
        Gn = Gnp[n,p]
        # denom += Gn*Nw
        denom = add_elementwise!(denom,Gn*Nw)
        for α in 1:Nw
            α´ = α
            for i_m in 1:m
                α´ = reconfigurationTable[α´,n-i_m]
            end
            O = ObsFunc(AllConfigs[n-m][α´])
            GnO = mult_elementwise!(O,Gn)
            # @. num += Gn*O
            num = add_elementwise!(num,GnO)
        end
    end
    return divide_elementwise!(num,denom)
end

function splitIntoBins(array,binsize)
    Iterators.partition(array,binsize)
end

function startManyWalkerGFMC(InitialState,Nwalkers,NSteps,nBranch,weightfunc::T,Λ) where T
    # Walkers = fetch.([Threads.@spawn spiderWebWalker(InitialState) for _ in 1:Nwalkers])# use threads to initialize walkers on correct NUMA domains (hopefully)
    Walkers = [spiderWebWalker(InitialState) for _ in 1:Nwalkers]

    weights = ones(Nwalkers)
    energies = zeros(NSteps)
    TotalWeights = zeros(NSteps)
    reconfigurationList = zeros(Int,Nwalkers)
    
    SaveConfigs = Vector{Vector{typeof(InitialState)}}()
    reconfTable = zeros(Int,Nwalkers,NSteps)

    for i in 1:NSteps
        # for (α,Config) in enumerate(Walkers)
        # for α in eachindex(Walkers)
        Threads.@threads for α in eachindex(Walkers)
            Walker = Walkers[α]

            w = 1
            for step in 1:nBranch
                move,weightList = performMarkovStep!(Walker,weightfunc,Λ)
                bx = sum(weightList)/size(InitialState,1)
                w *= bx
            end
            weights[α] = w
            getWeightList!(Walker,weightfunc,Λ)
        end
        energies[i] = getLocalEnergyWalkers_before(weights,Walkers,Λ)

        TotalWeights[i] = mean(weights)
        reconfiguration!(Walkers,reconfigurationList,weights)
        reconfTable[:,i] .= reconfigurationList

        # for (α,α´) in enumerate(reconfigurationList)
        #     if α == α´ # convention: surviving walkers always remain a copy at their index
        #         push!(SaveWalkers,copy(Walkers[α]))
        #     end
        # end

        push!(SaveConfigs,saveConfigs(Walkers))
    end
    return (;TotalWeights, energies, SaveConfigs,reconfTable)
end

function reconfiguration!(Walkers::AbstractVector{<:AbstractWalker},reconfigurationList,weights)
    StatsBase.sample!(eachindex(Walkers),StatsBase.Weights(weights),reconfigurationList,ordered = true)
    minimizeReconfiguration!(reconfigurationList)
    for (α,α´) in enumerate(reconfigurationList)
        if α´ != α
            get_config(Walkers[α]) .= get_config(Walkers[α´])
        end
    end
end

function getLocalEnergyWalkers_before(weights,Walkers::AbstractVector{<:AbstractWalker},Λ)
    Nw = length(weights)
    num = 0.
    denom = 0.

    for α in eachindex(weights,Walkers)
        bx = sum(get_weights(Walkers[α]))
        eloc = Λ - bx
        num += weights[α]*eloc
        denom += weights[α]
    end

    return num/denom
end

function saveConfigs(Walkers::AbstractVector{<:AbstractWalker})
    [copy(get_config(Walker)) for Walker in Walkers]
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
