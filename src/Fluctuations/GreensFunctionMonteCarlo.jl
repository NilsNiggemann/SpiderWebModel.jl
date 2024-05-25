"""Abstract supertype for a walker in a GFMC simulation. A walker needs to have a method `get_config` which returns the configuration (for example the spin configuration), and a method `get_weights` which returns a vector of weights for each possible move."""
abstract type AbstractWalker end

abstract type AbstractGuidingFunction end

struct SpiderWebWalker{C} <: AbstractWalker
    Config::C
    moves::Vector{Tuple{Int,Int,Int}}
    weights::Vector{Float64}
    Plaquette_positions::Vector{Tuple{Int,Int}}
    n_x::Vector{Float64}
    n_x´::Vector{Float64}
end
function spiderWebWalker(Config,Plaquette_positions)
    moves = Vector{Tuple{Int,Int,Int}}()
    weights = Vector{Float64}()
    # Plaquette_positions = collect(plaquetteIterator(Config))
    n_x = zeros(Float64,length(Plaquette_positions))
    n_x´ = zeros(Float64,length(Plaquette_positions))
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

struct PlaquetteNumberGuidingFunction <: AbstractGuidingFunction
    α::Float64
end
(ψG::PlaquetteNumberGuidingFunction)(ΔNPlaq::Real) = exp(ψG.α*ΔNPlaq)

guidingfunc_name(F::Function) = string(typeof(F))
guidingfunc_name(F::AbstractGuidingFunction) = string(typeof(F))
guidingfunc_name(F::PlaquetteNumberGuidingFunction) = "PlaquetteNumberGuidingFunction"

variational_parameters(P::PlaquetteNumberGuidingFunction) = Dict([:alpha=>P.α])
variational_parameters(P::Function) = Dict([x => getproperty(P,x) for x in propertynames(P)])

function varitationalFunc(α,NPlaq::Real,NPlaqEstimate)
    return exp(α*(NPlaq-NPlaqEstimate))
end

struct VariationalGuidingFunction{A<:AbstractArray} <: AbstractGuidingFunction
    params::A
end

function constructVariationalFunction(State,α::Real=0.1)
    plaqs = collect(plaquetteIterator(State))
    N = length(plaqs)
    params = zeros(N,N+1)
    ψ = VariationalGuidingFunction(params)
    get_alpha_i(ψ) .= α
    return ψ
end

function get_alpha_i(ψG::VariationalGuidingFunction)
    return @view ψG.params[:,begin]
end

function get_beta_ij(ψG::VariationalGuidingFunction)
    return @view ψG.params[:,begin+1:end]
end

function (ψG::VariationalGuidingFunction)(N□::AbstractArray) 
    return exp(guidingfunc_exponent(ψG,N□))
end

function guidingfunc_exponent(ψG::VariationalGuidingFunction,N□::AbstractArray) 
    α = get_alpha_i(ψG)
    β = get_beta_ij(ψG)
    exponent = α' * N□ + dot(N□,β,N□)
    return exponent
end

# function (ψG::VariationalGuidingFunction)(N□´::AbstractArray,N□::AbstractArray) 
#     α = get_alpha_i(ψG)
#     β = get_beta_ij(ψG)
#     exponent = zero(eltype(α))

#     for i in eachindex(N□´,N□)
#         exponent += α[i]*(N□´[i] - N□[i])
#         for j in eachindex(N□´,N□)  # can be optimized using j < i
#             exponent += β[i,j]*(N□´[i]*N□´[j] - N□[i]*N□[j])
#         end
#     end
#     return exp(exponent)
# end

function guidingfuncRatio(ψG::VariationalGuidingFunction,n::AbstractArray,n´::AbstractArray) 
    return exp(guidingfuncRatio_exponent(ψG,n,n´))
end

function guidingfuncRatio_exponent(ψG::VariationalGuidingFunction,n::AbstractArray,n´::AbstractArray) 
    α = get_alpha_i(ψG)
    β = get_beta_ij(ψG)

    exponent_α = sum(ai*(n´i - ni) for (ai,n´i,ni) in zip(α,n´,n))

    exponent_β = zero(eltype(α))
    exponent_β = dot(n´,β,n´) - dot(n,β,n)

    return (exponent_α + exponent_β)
end

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
    IndexDict = Dict(AllPlaqs[i] => i for i in eachindex(AllPlaqs))
    A = parent(Config)
    for i_vic in i-2:i+2
        for j_vic in j-2:j+2
            (i´,j´) = get_wrappend_inds(A,(i_vic,j_vic))
            isodd(i´+j´) || continue
            index = get(IndexDict,(i´,j´),0)
            index !== 0 && push!(Plaq_indices,index)
        end
    end
    return Plaq_indices
end

function precomputeAffectedPlaquettes(Config)
    AffPlaqMatrix = Matrix{Vector{Int}}(undef,size(Config))

    for I in plaquetteIterator(Config)
        i,j = I
        AffPlaqs = Vector{Int}()
        findAffectedPlaquettes!(AffPlaqs,Config,i,j)
        AffPlaqMatrix[i,j] = AffPlaqs
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
    # println(n_x´,affected_indices)
    # println(Walker.n_x)
    # error("")
    return n_x´
end


function getNPlaqfilled!(Walker::SpiderWebWalker,affected_indices)
    (;n_x,n_x´,Config,Plaquette_positions) = Walker
    # println(affected_indices)
    # empty!(n_x´)
    # resize!(n_x´,length(n_x))

    n_x´ .= n_x

    for PlaqIndex in affected_indices
        I = Plaquette_positions[PlaqIndex]
        @inbounds applPlus, applMinus = P_applicable(Config, I)
        n = applPlus + applMinus
        n_x´[PlaqIndex] = n
    end
    # println(n_x´)
    # error("")
    return n_x´
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

function getWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,weightfunc::VariationalGuidingFunction,Λ)
    (;Config,moves,weights) = Walker
    empty!(weights)

    n_x = getNPlaq!(Walker)
    
    for operator in moves
        i,j,opNum = operator
        indices = AffectedPlaquetteList[i,j]
        applyPlaquette!(Config, i, j, opNum)
         
        n_x´ = getNPlaqfilled!(Walker,indices)
        
        weight = guidingfuncRatio(weightfunc,n_x,n_x´)
        # ndiff = sum(n_x´ .- n_x)
        # weight = weightfunc(ndiff)
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
        Walkers[α] = spiderWebWalker(InitialState,plaquettePositions)
    end
    weights = ones(Nwalkers)
    TotalWeights = zeros(NSteps)
    reconfiguration_buffer = zeros(Nwalkers)
    reconfigurationTable = zeros(Int,Nwalkers,NSteps)
    
    return (;AffectedPlaquetteList,Walkers,weights,TotalWeights,reconfiguration_buffer,reconfigurationTable)
end

function setupObservables(InitConfig,NWalkers,NSteps,outfile::Nothing)
    energies = zeros(NSteps)
    Lx,Ly = size(InitConfig)
    SaveConfigs = zeros(eltype(InitConfig),Lx,Ly,NWalkers,NSteps)
    return (;energies,SaveConfigs)
end

function createMMapArray(file::HDF5.File,datasetname::String,type,dims)
    SaveConfigs_dset = create_dataset(file,datasetname,datatype(type),dataspace(dims);alloc_time = HDF5.H5D_ALLOC_TIME_EARLY)
    @assert HDF5.ismmappable(SaveConfigs_dset) "Dataset is not mappable for given type $(eltype(InitConfig))"
    return HDF5.readmmap(SaveConfigs_dset)
end

function readMMapArray(filename::AbstractString,datasetname::String)
    h5open(filename,"r") do file
        SaveConfigs_dset = file[datasetname]
        @assert HDF5.ismmappable(SaveConfigs_dset) "Dataset is not mappable for given type $(eltype(InitConfig))"
        return HDF5.readmmap(SaveConfigs_dset)
    end
end

function setupObservables(InitConfig,NWalkers,NSteps,filename::String)
    energies = zeros(NSteps)
    Lx,Ly = size(InitConfig)
    SaveConfigs = h5open(filename,"cw") do file
        createMMapArray(file,"SaveConfigs",eltype(InitConfig),(Lx,Ly,NWalkers,NSteps))
    end
    return (;energies,SaveConfigs)
end

function saveParameters(filename::String,Λ,equilibration_steps,nBranch,weightfunc,w_avg_estimate)
    h5open(filename,"cw") do file
        file["Λ"] = Λ
        file["equilibration_steps"] = equilibration_steps
        file["nBranch"] = nBranch
        file["w_avg_estimate"] = w_avg_estimate
        saveVariationalParameter(file,weightfunc)
    end
end
saveParameters(::Nothing,args...) = nothing

function saveVariationalParameter(file::HDF5.File,weightfunc)
    pars = variational_parameters(weightfunc)
    funcName = guidingfunc_name(weightfunc)
    for (key,val) in pars
        file[string(funcName,"/",key)] = val
    end
end

function startManyWalkerGFMC(InitialState::ConfType,outfile,Nwalkers::Int,NSteps::Int,equilibration_steps::Int,nBranch::Int,weightfunc::Fun,Λ::Real,w_avg_estimate) where {T,ConfType <: StencilSpinConfig{T},Fun}
    
    (;AffectedPlaquetteList,Walkers,weights,TotalWeights,reconfiguration_buffer,reconfigurationTable) = setup_many_walker_GFMC(InitialState,Nwalkers,NSteps)
    (;energies,SaveConfigs) = setupObservables(InitialState,Nwalkers,NSteps,outfile)
    saveParameters(outfile,Λ,equilibration_steps,nBranch,weightfunc,w_avg_estimate)

    for _ in 1:equilibration_steps
        propagateWalkers!(Walkers,weights,AffectedPlaquetteList,weightfunc,Λ,nBranch,w_avg_estimate)
        reconfigurationList = @view reconfigurationTable[:,1]
        reconfiguration!(Walkers,reconfigurationList,reconfiguration_buffer,weights)
    end

    for i in 1:NSteps
        # for (α,Config) in enumerate(Walkers)
        propagateWalkers!(Walkers,weights,AffectedPlaquetteList,weightfunc,Λ,nBranch,w_avg_estimate)

        energies[i] = getLocalEnergyWalkers_before(weights,Walkers,Λ)
        TotalWeights[i] = mean(weights)
        
        reconfigurationList = @view reconfigurationTable[:,i]
        reconfiguration!(Walkers,reconfigurationList,reconfiguration_buffer,weights)

        saveConfigs!(SaveConfigs,i,Walkers)            
    end
    saveObservables(outfile,TotalWeights,energies,reconfigurationTable)
    return (;TotalWeights, energies, SaveConfigs, reconfigurationTable)
end

function startManyWalkerGFMC(InitialState,Nwalkers,NSteps,nBranch,weightfunc,Λ;equilibration_steps = 0,outfile=nothing,w_avg_estimate=size(InitialState,1))
    startManyWalkerGFMC(InitialState,outfile,Nwalkers,NSteps,equilibration_steps,nBranch,weightfunc,Λ,w_avg_estimate)
end


function propagateWalkers!(Walkers,weights,AffectedPlaquetteList,weightfunc::Fun,Λ,nBranch,w_avg_estimate) where {Fun}
    # L = size(get_config(first(Walkers)),1)
    # for α in eachindex(Walkers)
    w_avg_estimate⁻¹ = 1. / w_avg_estimate
    Threads.@threads for α in eachindex(Walkers)
        Walker = Walkers[α]
        w = 1.
        for step in 1:nBranch
            move,weightList = performMarkovStep!(Walker,AffectedPlaquetteList,weightfunc,Λ)
            bx = sum(weightList)*w_avg_estimate⁻¹
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

function saveObservables(outfile,TotalWeights,energies,reconfigurationTable)
    h5open(outfile,"cw") do file
        file["energies"] = energies
        file["TotalWeights"] = TotalWeights
        file["reconfigurationTable"] = reconfigurationTable
    end
end

saveObservables(::Nothing,args...;kwargs...) = nothing

function saveConfigs!(SaveConfigs,i,Walkers::AbstractVector{<:AbstractWalker})
    for (α,Config) in enumerate(Walkers)
        SaveConfigs[:,:,α,i] .= get_config(Config)
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

