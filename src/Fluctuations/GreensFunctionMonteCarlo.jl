"""Abstract supertype for a walker in a GFMC simulation. A walker needs to have a method `get_config` which returns the configuration (for example the spin configuration), and a method `get_weights` which returns a vector of weights for each possible move."""
abstract type AbstractWalker end

struct SpiderWebWalker{C} <: AbstractWalker
    Config::C
    moves::Vector{Tuple{Int8,Int8,Int8}}
    weights::Vector{Float64}
    Plaquette_positions::Vector{Tuple{Int,Int}}
    n_x::Vector{Int8}
    n_x´::Vector{Int8}
end
function SpiderWebWalker(Config,Plaquette_positions)
    moves = Vector{Tuple{Int8,Int8,Int8}}()
    weights = Vector{Float64}()
    # Plaquette_positions = collect(plaquetteIterator(Config))
    n_x = zeros(Int8,length(Plaquette_positions))
    n_x´ = Vector{Int8}()
    return SpiderWebWalker(copy(Config),moves,weights,Plaquette_positions,n_x,n_x´)
end

get_config(Walker::SpiderWebWalker) = Walker.Config
get_weights(Walker::SpiderWebWalker) = Walker.weights
plaquetteIterator(Walker::SpiderWebWalker) = Walker.Plaquette_positions

getOperatorRep(i,j,opNum) = i,j,opNum
# getOperatorRep(i,j,opNum) = CartesianIndex(i,j,opNum)
# getOperatorNumber(i,j,opNum,L) = LinearIndices((L,L,2))[i,j,opNum]

# function getOperatorFromNumber(Config, opNumber, operators)

# end

function getMoves!(
    Walker::SpiderWebWalker
)

    Conf = get_config(Walker)
    moves = Walker.moves
    empty!(moves)
    for I in plaquetteIterator(Walker)
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
# getMoves!(Walker::SpiderWebWalker) = getMoves!(Walker.moves,Walker)

import StatsBase

function findAffectedPlaquettes!(Plaq_indices,S,i,j) 
    A = parent(S)
    bound = Stencils.boundary(A)
    pad = Stencils.padding(A)
    findAffectedPlaquettes!(Plaq_indices,S,i,j,bound,pad)
end
function findAffectedPlaquettes!(Plaq_indices,Config,i,j,::Stencils.Remove,::Stencils.Conditional)
    empty!(Plaq_indices)
    for (index,I) in enumerate(plaquetteIterator(Config))
        if !plaquettesAreSeparated(I,(i,j))
            push!(Plaq_indices,index)
        end
    end
    return Plaq_indices
end
function findAffectedPlaquettes!(Plaq_indices,Config,i,j,::Stencils.Wrap,::Stencils.Conditional)
    empty!(Plaq_indices)
    AllPlaqs = collect(plaquetteIterator(Config))
    A = parent(Config)
    for i´ in i-2:i+2
        for j´ in j-2:j+2
            isodd(i´+j´) || continue
            (i´,j´) = get_wrappend_inds(A,(i´,j´))
            index = findfirst(isequal((i´,j´)),AllPlaqs)
            push!(Plaq_indices,index)
        end
    end
    return Plaq_indices
end

function precomputeAffectedPlaquettes(Config)
    AffPlaqMatrix = [Vector{Int}() for i in 1:size(Config,1), j in 1:size(Config,2)]

    for (index,I) in enumerate(plaquetteIterator(Config))
        i,j = I
        findAffectedPlaquettes!(AffPlaqMatrix[i,j],Config,i,j)
    end
    return AffPlaqMatrix
end
# Assumes that length(n_plaq) == length(plaquetteIterator(Config))
function getNPlaq!(Walker::SpiderWebWalker)
    (;n_x,Config) = Walker
    for (i,I) in enumerate(plaquetteIterator(Walker))
        @inbounds applPlus, applMinus = P_applicable(Config, I)
        n = applPlus + applMinus
        n_x[i] = n
    end
    return n_x
end

function getNPlaq!(Walker::SpiderWebWalker,affected_indices)
    ind = 1
    (;n_x´,Config,Plaquette_positions) = Walker
    resize!(n_x´,length(affected_indices))
    # println(affected_indices)
    # empty!(n_x´)
    for (i,PlaqIndex) in enumerate(affected_indices)
        I = Plaquette_positions[PlaqIndex]
        @inbounds applPlus, applMinus = P_applicable(Config, I)
        n = applPlus + applMinus
        n_x´[i] = n
    end
    return n_x´
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

function getWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,weightfunc::T,Λ) where T
    (;Config,moves,weights,n_x,n_x´) = Walker
    empty!(weights)

    getNPlaq!(Walker)
    
    for operator in moves
        i,j,opNum = operator
        indices = AffectedPlaquetteList[i,j]
        applyPlaquette!(Config, i, j, opNum)
        
        getNPlaq!(Walker,indices)
        
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

function performMarkovStep!(Walker::SpiderWebWalker,AffectedPlaquetteList,weightfunc::T,Λ) where T
    moves = getMoves!(Walker)
    for (i,j,opNum) in moves
        @assert i != 0 && j != 0 "err"
    end
    weights = getWeightList!(Walker,AffectedPlaquetteList,weightfunc,Λ)

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
        moves = getMoves!(moves,Config)
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
    return num ./denom
end


function setupProjector(weights,nThermal)
    bn = @view weights[nThermal:end]
    meanweight = mean(bn)
    Gn1 = bn./meanweight
end

add_elementwise!(x::AbstractArray,y) = (x .+= y)
add_elementwise!(x::Number,y::Number) = x + y

mult_elementwise!(x::Number,y::Number) = x * y
mult_elementwise!(x::AbstractArray,y::Number) = (x .*= y)

divide_elementwise!(x::Number,y::Number) = x / y
divide_elementwise!(x::AbstractArray,y::AbstractArray) = (x ./= y)

function getObs(Gnp,AllConfigs,reconfigurationTable,ObsFunc,m=size(Gnp,2)÷2)
    N = lastindex(AllConfigs,4)
    exampleConf = @view AllConfigs[:,:,begin,begin]
    Obs = ObsFunc(exampleConf)
    num = zero(Obs)
    denom = zero(Obs)
    # fill!(Obs,zero(eltype(Obs)))

    Nw = size(reconfigurationTable,1)
    p = size(Gnp,2)
    # surviving_walker_mapping_list = zeros(Int,Nw)

    for n in m+1:N
        Gn = Gnp[n,p]
        # denom += Gn*Nw
        denom = add_elementwise!(denom,Gn*Nw)
        # reconfigList = @view reconfigurationTable[:,n-m]
        # surviving_walker_mapping!(surviving_walker_mapping_list,reconfigList)
        for α in 1:Nw
            α´ = α
            for i_m in 1:m
                α´ = reconfigurationTable[α´,n-i_m]
            end
            conf = @view AllConfigs[:,:,α´,n-m]
            O = ObsFunc(conf)
            # surviving_index = surviving_walker_mapping_list[α´]
            # O = ObsFunc(AllConfigs[n-m][surviving_index])
            GnO = mult_elementwise!(O,Gn)
            # @. num += Gn*O
            num = add_elementwise!(num,GnO)
        end
    end
    return divide_elementwise!(num,denom)
end

function surviving_walker_mapping!(mappingarr,reconfigList)
    fill!(mappingarr,0)
    surviving_walkers = 0
    for (α,α´) in enumerate(reconfigList)
        if α == α´
            surviving_walkers += 1
            mappingarr[α´] = surviving_walkers
        end
    end
    # println(reconfigList)
    # println(mappingarr)
    for α in eachindex(reconfigList,mappingarr)
        if mappingarr[α] == 0
            α´ = reconfigList[α]
            # println((α,α´,mappingarr[α´]))
            mappingarr[α] = mappingarr[α´]
        end
    end
    return mappingarr
end

function splitIntoBins(array,binsize)
    Iterators.partition(array,binsize)
end

function setup_many_walker_GFMC(InitialState::ConfType,Nwalkers,NSteps) where {ConfType <: StencilSpinConfig}
    AffectedPlaquetteList = precomputeAffectedPlaquettes(InitialState)
    plaquettePositions = collect(plaquetteIterator(InitialState))
    Walkers = Vector{SpiderWebWalker{ConfType}}(undef,Nwalkers)
    Threads.@threads for α in eachindex(Walkers)
        Walkers[α] = SpiderWebWalker(InitialState,plaquettePositions)
    end
    weights = ones(Nwalkers)
    TotalWeights = zeros(NSteps)
    reconfiguration_buffer = zeros(Nwalkers)
    reconfigurationList = zeros(Int,Nwalkers)
    
    return (;AffectedPlaquetteList,Walkers,weights,TotalWeights,reconfiguration_buffer,reconfigurationList)
end

function setupObservables(InitialState,Nwalkers,NSteps)
    energies = zeros(NSteps)
    ConfigDims = size(InitialState)
    SaveConfigs = zeros(eltype(InitialState),ConfigDims...,Nwalkers,NSteps)
    reconfTable = zeros(Int,Nwalkers,NSteps)
    return (;energies,SaveConfigs,reconfTable)
end

function startManyWalkerGFMC(InitialState::ConfType,Nwalkers,NSteps,nBranch,weightfunc::Fun,Λ,saveObs=true) where {T,ConfType <: StencilSpinConfig{T},Fun}
    
    (;AffectedPlaquetteList,Walkers,weights,TotalWeights,reconfiguration_buffer,reconfigurationList) = setup_many_walker_GFMC(InitialState,Nwalkers,NSteps)
    (;energies,SaveConfigs,reconfTable) = setupObservables(InitialState,Nwalkers,NSteps)

    for i in 1:NSteps
        # for (α,Config) in enumerate(Walkers)
        propagateWalkers!(Walkers,weights,AffectedPlaquetteList,weightfunc,Λ,nBranch)

        energies[i] = getLocalEnergyWalkers_before(weights,Walkers,Λ)

        TotalWeights[i] = mean(weights)
        reconfiguration!(Walkers,reconfigurationList,reconfiguration_buffer,weights)
        reconfTable[:,i] .= reconfigurationList

        CurrentConfs = @view SaveConfigs[:,:,:,i]
        saveConfigs!(CurrentConfs,Walkers)
    end
    return (;TotalWeights, energies, SaveConfigs,reconfTable)
end

function propagateWalkers!(Walkers,weights,AffectedPlaquetteList,weightfunc::Fun,Λ,nBranch) where {Fun}
    L = size(get_config(first(Walkers)),1)
    # for α in eachindex(Walkers)
    Threads.@threads for α in eachindex(Walkers)
        Walker = Walkers[α]
        w = 1
        for step in 1:nBranch
            move,weightList = performMarkovStep!(Walker,AffectedPlaquetteList,weightfunc,Λ)
            bx = sum(weightList)/L
            w *= bx
        end
        weights[α] = w
        getMoves!(Walker)
        getWeightList!(Walker,AffectedPlaquetteList,weightfunc,Λ)
    end
end
# function reconfiguration!(Walkers::AbstractVector{<:AbstractWalker},reconfigurationList,weights)
#     StatsBase.sample!(eachindex(Walkers),StatsBase.Weights(weights),reconfigurationList,ordered = true)
#     minimizeReconfiguration!(reconfigurationList)
#     for (α,α´) in enumerate(reconfigurationList)
#         if α´ != α
#             get_config(Walkers[α]) .= get_config(Walkers[α´])
#         end
#     end
# end
"""Performs an efficient reconfiguration of walkers. This reconfiguration will not remove walkers if they all have the same weight, which increases the efficiency as more walkers can contribute to the average.

Matteo Calandra Buonaura and Sandro Sorella
Phys. Rev. B 57, 11446 (1998)
"""
function reconfiguration!(Walkers::AbstractVector{<:AbstractWalker},reconfigurationList,reconfiguration_buffer,weights)
    Nw = length(Walkers)
    reconfiguration_buffer = cumsum!(reconfiguration_buffer,weights) 
    wTotal = sum(weights)
    reconfiguration_buffer ./= wTotal
    for α in eachindex(Walkers)
        ξα = rand()
        zα = (α + ξα - 1)/Nw
        α´ = searchsortedfirst(reconfiguration_buffer,zα)
        reconfigurationList[α] = α´
    end
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

function saveConfigs!(SaveConfigs,Walkers::AbstractVector{<:AbstractWalker})
    for (α,Config) in enumerate(Walkers)
        SaveConfigs[:,:,α] .= get_config(Config)
    end

end

"""given a list of reconfiguration indices, minimizes the number of reconfigurations by swapping elements in the list. Each walker that survives a reconfiguration step remains unchanged while walkers that are killed get assigned to a new index."""
function minimizeReconfiguration!(list)
    for (α,α´) in enumerate(list)
        if α´ != α
            otherIndex = findfirst(isequal(α),list)
            isnothing(otherIndex) && continue
            swapIndices!(list,α,otherIndex)
        end
    end
    return list
end

function swapIndices!(list,i,j)
    list[i],list[j] = list[j],list[i]
    return list
end