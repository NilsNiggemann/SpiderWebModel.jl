function _setup_GFMC_problem(InitialState::StencilSpinConfig,method::AbstractGFMCMethod,Nwalkers::Integer,NSteps::Integer,nThreads,ψG,outfile)
    setup = setup_many_walker_GFMC(InitialState,Nwalkers,nThreads)
    Guiding_function_buffer = allocate_GWF_buffers_threads(ψG,InitialState,Nwalkers)
    
    (;Walkers,weights,reconfiguration_buffer) = setup
    ObsSetup = setupObservables(InitialState,Nwalkers,NSteps,outfile)
    (;energies,SaveConfigs,TotalWeights,reconfigurationTable) = ObsSetup

    Observables = GFMCObservables(energies,SaveConfigs,TotalWeights,reconfigurationTable,outfile)

    return SpiderwebGFMCProblem(method,InitialState,ψG,Walkers,weights,Guiding_function_buffer,reconfiguration_buffer,Observables)
end

function setup_GFMC_problem(InitialState::ConfType, method::AbstractGFMCMethod, Nwalkers::Integer, NSteps::Integer,nThreads::Integer, ψG;
    outfile = nothing
    ) where {ConfType <: StencilSpinConfig}
    _setup_GFMC_problem(InitialState,method,Nwalkers,NSteps,nThreads,ψG,outfile)
end

function getMoves!(
    Walker::SpiderWebWalker
)
    Conf = get_config(Walker)
    moves = Walker.moves
    getMoves!(moves,Conf)
end

function getMoves!(
    moves::AbstractVector,
    Conf::StencilSpinConfig
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
# getMoves!(Walker::SpiderWebWalker) = getMoves!(Walker.moves,Walker)

function plaquettesAreSeparated(P1, P2)
    separation = P2 .- P1

    return any(x -> abs(x) > 2, separation)
end

function findAffectedPlaquettes!(Plaq_indices,S,i,j) 
    A = parent(S)
    bound = Stencils.boundary(A)
    pad = Stencils.padding(A)
    findAffectedPlaquettes!(Plaq_indices,S,i,j,bound,pad)
end
function findAffectedPlaquettes!(Plaq_indices,Config,i,j,::Stencils.Remove,::Any)
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
    # AffPlaqMatrix = Matrix{OrderedSet{Int}}(undef,size(Config))
    AffPlaqMatrix = Matrix{Vector{Int}}(undef,size(Config))

    for I in plaquetteIterator(Config)
        i,j = I
        AffPlaqs = OrderedSet{Int}()
        findAffectedPlaquettes!(AffPlaqs,Config,i,j)
        AffPlaqMatrix[i,j] = collect(AffPlaqs)
    end
    return AffPlaqMatrix
end

# Assumes that length(n_plaq) == length(plaquetteIterator(Config))
function getNPlaq!(n_x::AbstractVector,Config::StencilSpinConfig)
    for (i,I) in enumerate(plaquetteIterator(Config))
        @inbounds applPlus, applMinus = P_applicable(Config, I)
        n = applPlus + applMinus
        n_x[i] = n
    end
    return n_x
end
getNPlaq!(Walker::SpiderWebWalker) = getNPlaq!(Walker.n_x,Walker.Config)
getNPlaq(Config::StencilSpinConfig) = getNPlaq!(zeros(length(collect(plaquetteIterator(Config)))),Config)

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
    n_x´ .= n_x

    for PlaqIndex in affected_indices
        I = Plaquette_positions[PlaqIndex]
        @inbounds applPlus, applMinus = P_applicable(Config, I)
        n = applPlus + applMinus
        n_x´[PlaqIndex] = n
    end
    return n_x´
end

function getNPlaq_difference(nPlaq_x,nPlaq_x´,affectedPlaquettes)
    N□ = 0
    for (i,plaqIndex) in enumerate(affectedPlaquettes)
        N□ += nPlaq_x´[i] - nPlaq_x[plaqIndex]
    end
    return N□
end

const DIAGONAL_MOVE_ID = (0,0,0)

function performMarkovStep!(Walker::SpiderWebWalker)
    (;moves,weights) = Walker
    moveidx = StatsBase.sample(StatsBase.Weights(weights))
    move = Walker.moves[moveidx]
    if move != DIAGONAL_MOVE_ID
        applyPlaquette!(Walker.Config, move[1], move[2], move[3]) 
    end
    empty!(weights)
    empty!(moves)
    return move
end

function getLocalEnergy(weights,Λ=0)
    return -sum(weights) + Λ
end

function getLocalEnergy(Walker::SpiderWebWalker,method::DiscreteTimeMethod)
    return -sum(Walker.weights) + method.Λ
end

function getLocalEnergy(Walker::SpiderWebWalker,method::ContinuousTimeMethod)
    return -sum(Walker.weights) + method.Hxx(Walker)
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
        Gnp[n,1] = 1 #zero projection order
        Gnp[n,2] = bn[n]/meanweight # first projection order
    end
    for p in 3:PMax
        for n in p:length(bn)
            # Gnp[n,p] = Gnp[n,p-1]*bn[n-p]/meanweight
            Gnp[n,p] = Gnp[n-1,p-1]*Gnp[n,2]
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

function getEnergies(weights,localEnergies,nthermalization,PMax;
    Gnp = precomputeNormalizedAccWeight(weights,nthermalization,PMax)
    )
    
    EL_thermalized = @view localEnergies[nthermalization:end]
    N = lastindex(EL_thermalized)
    num = zeros(PMax)
    denom = zeros(PMax)
    for p in 1:PMax
        for n in p+1:N
            # G_np = p != 0 ? Gnp[n,p] : 1
            num[p] += Gnp[n,p]*EL_thermalized[n]
            denom[p] += Gnp[n,p]
        end
    end
    return num ./denom
end

function getEnergies(results,nthermalization,PMax);
    return getEnergies(results.TotalWeights,results.energies,nthermalization,PMax)
end


function setupProjector(weights,nThermal)
    bn = @view weights[nThermal:end]
    meanweight = mean(bn)
    Gn1 = bn./meanweight
end

"""
Allocates walkers in a multithreaded way such that the memory of each walker is associated with the correct memory domain
...hopefully...
"""
function _allocateWalkers(InitialState::ConfType,Nwalkers,nThreads,plaquettePositions) where {ConfType <: StencilSpinConfig}
    Walkers = Vector{SpiderWebWalker{ConfType}}(undef,Nwalkers)
    chunks = ChunkSplitters.chunks(eachindex(Walkers), n = nThreads)
    Threads.@threads for (i_chunk,αinds) in enumerate(chunks)
        for α in αinds
            Walkers[α] = spiderWebWalker(InitialState,plaquettePositions)
        end
    end
    return Walkers
end

function setup_many_walker_GFMC(InitialState::ConfType,Nwalkers::Integer,nThreads::Integer) where {ConfType <: StencilSpinConfig}
    plaquettePositions = collect(plaquetteIterator(InitialState))
    Walkers = _allocateWalkers(InitialState,Nwalkers,nThreads,plaquettePositions)
    weights = ones(Nwalkers)
    reconfiguration_buffer = zeros(Nwalkers)
    
    return (;Walkers,weights,reconfiguration_buffer)
end

function setupObservables(InitConfig,NWalkers,NSteps,outfile::Nothing)
    energies = zeros(NSteps)
    Lx,Ly = size(InitConfig)
    SaveConfigs = zeros(eltype(InitConfig),Lx,Ly,NWalkers,NSteps)

    TotalWeights = zeros(NSteps)
    reconfigurationTable = zeros(Int,NWalkers,NSteps)

    return (;energies,SaveConfigs,TotalWeights,reconfigurationTable)
end
function setupObservables(InitConfig,NWalkers,NSteps,filename::String)
    Lx,Ly = size(InitConfig)
    h5open(filename,"cw") do file
        energies = createMMapArray(file,"energies",Float64,(NSteps,))
        SaveConfigs = createMMapArray(file,"SaveConfigs",eltype(InitConfig),(Lx,Ly,NWalkers,NSteps))
        TotalWeights = createMMapArray(file,"TotalWeights",Float64,(NSteps,))
        reconfigurationTable = createMMapArray(file,"reconfigurationTable",Int,(NWalkers,NSteps))
        return (;energies,SaveConfigs,TotalWeights,reconfigurationTable)
    end
end

function createMMapArray(file::HDF5.File,datasetname::String,type,dims)
    SaveConfigs_dset = create_dataset(file,datasetname,datatype(type),dataspace(dims);alloc_time = HDF5.H5D_ALLOC_TIME_EARLY)
    @assert HDF5.ismmappable(SaveConfigs_dset) "Dataset is not mappable for given type $(eltype(InitConfig))"
    return HDF5.readmmap(SaveConfigs_dset)
end

function readMMapArray(filename::AbstractString,datasetname::String)
    h5open(filename,"r") do file
        SaveConfigs_dset = file[datasetname]
        return HDF5.readmmap(SaveConfigs_dset)
    end
end

function saveParameters(filename::String,Λ,equilibration_steps,ψG,w_avg_estimate;kwargs...)
    h5open(filename,"cw") do file
        file["Λ"] = Λ
        file["equilibration_steps"] = equilibration_steps
        file["w_avg_estimate"] = w_avg_estimate
        for (k,v) in kwargs
            file[String(k)] = v
        end
        saveVariationalParameters(file,ψG)
    end
end
function saveParameters(filename::String,equilibration_steps,method::DiscreteTimeMethod,ψG)
    (;Λ,nBranch,w_avg_estimate) = method
    saveParameters(filename,Λ,equilibration_steps,ψG,w_avg_estimate;nBranch=nBranch)
end
function saveParameters(filename::String,equilibration_steps,method::ContinuousTimeMethod,ψG)
    w_avg_estimate = method.w_avg_estimate
    saveParameters(filename,Inf,equilibration_steps,ψG,w_avg_estimate;τ=method.τ)
end

saveParameters(::Nothing,args...) = nothing

function saveVariationalParameters(file::HDF5.File,ψG)
    pars = get_params(ψG)
    funcName = guidingfunc_name(ψG)
    _save_h5_array(file,string(funcName,"/params"),pars)
end
function _save_h5_array(file,datasetname,arr::RecursiveArrayTools.ArrayPartition)
    for i in eachindex(arr.x)
        file[datasetname*"/"*string(i)] = arr.x[i]
    end
end
function _save_h5_array(file,datasetname,arr)
    file[datasetname] = arr
end

abstract type AbstractGFMCInitializer end

struct UnguidedWalkInitializer <: AbstractGFMCInitializer 
    pre_equilibration_steps::Int
    scatter_fraction::Float64
end

function initialize!(Walkers::AbstractVector{<:SpiderWebWalker},I::UnguidedWalkInitializer)
    if I.pre_equilibration_steps > 0 && I.scatter_fraction > 0
        random_init_walkers!(Walkers,I.pre_equilibration_steps,I.scatter_fraction)
    end
    return
end

"""draw random configurations with given pre-given weights to accelerate equilibration times"""
struct WeightedConfigsInitializers{T1<:AbstractVector{<:AbstractMatrix},T2<:AbstractVector} <: AbstractGFMCInitializer
    configs::T1
    weights::T2
end

function WeightedConfigsInitializers(SaveConfigs::AbstractArray{<:Number,4},TotalWeights::AbstractVector)
    configs = collect(eachslice(SaveConfigs,dims=(3,4)))
    configsVec = reshape(configs,length(configs))
    weights = [w for w in TotalWeights for _ in axes(SaveConfigs,3)]
    return WeightedConfigsInitializers(configsVec,weights)
end

function WeightedConfigsInitializers(resultsArr::AbstractVector,weight::Symbol=:TotalWeights)

    w1 = WeightedConfigsInitializers(resultsArr[begin].SaveConfigs,getproperty(resultsArr[begin],weight))
    for i in eachindex(resultsArr)[2:end]
        w2 = WeightedConfigsInitializers(resultsArr[i].SaveConfigs,getproperty(resultsArr[i],weight))
        append!(w1.configs,w2.configs)
        append!(w1.weights,w2.weights)
    end
    return w1
end

function initialize!(Walkers::AbstractVector{<:SpiderWebWalker},I::WeightedConfigsInitializers)
    for Walker in Walkers
        rand_conf = StatsBase.sample(I.configs,StatsBase.Weights(I.weights))
        get_config(Walker) .= rand_conf
    end
end

struct CombinedInitializer{I1 <: AbstractGFMCInitializer,I2 <: AbstractGFMCInitializer} <: AbstractGFMCInitializer
    I1::I1
    I2::I2
end

function initialize!(Walkers::AbstractVector{<:SpiderWebWalker},I::CombinedInitializer)
    initialize!(Walkers,I.I1)
    initialize!(Walkers,I.I2)
end

function initializeGFMC!(prob::AbstractGFMCProblem,equilibration_steps=0,initializer = UnguidedWalkInitializer(equilibration_steps ÷ 5,0.8))
    
    (;Walkers,Observables,method,ψG) = prob
    initialize!(Walkers,initializer)

    return prob,Observables
end

function equilibrate!(prob::AbstractGFMCProblem,equilibration_steps;nThreads=Threads.nthreads(),reconfigure=true)
    runGFMC!(prob,equilibration_steps,nThreads,reconfigure,false,false)
end

function pre_estimate_energies!(Observables,weights,Walkers,method,i,equilibration_steps)
    energies = get_energies(Observables)
    TotalWeights = get_TotalWeights(Observables)
    
    N_steps = length(energies)
    isave = i - equilibration_steps + N_steps
    if isave ∉ eachindex(energies,TotalWeights)
        return
    end
    energies[isave] = getLocalEnergyWalkers_before(weights,Walkers,method)
    TotalWeights[isave] = mean(weights)
    return
end

function startManyWalkerGFMC!(prob::AbstractGFMCProblem,NStepsRange,nThreads::Int,equilibration_steps::Int,initializer = UnguidedWalkInitializer(equilibration_steps ÷ 5,0.8))
    initializeGFMC!(prob,nThreads,initializer)
    equilibrate!(prob,equilibration_steps;nThreads)
    outfile = get_outfile(prob.Observables)
    saveParameters(outfile,equilibration_steps,prob.method,prob.ψG)

    fill_all_Buffers!(prob,nThreads)
    runGFMC!(prob,NStepsRange;nThreads)
end
function startManyWalkerGFMC(InitialState::StencilSpinConfig,method::AbstractGFMCMethod,Nwalkers::Integer,nSteps::Integer,ψG; equilibration_steps = 0, pre_equilibration_steps = equilibration_steps ÷ 5, scatter_fraction = 0.8,initializer = UnguidedWalkInitializer(pre_equilibration_steps,scatter_fraction),nThreads=2*Threads.nthreads(),kwargs...)
    prob = setup_GFMC_problem(InitialState,method,Nwalkers,nSteps,nThreads,ψG;kwargs...)
    startManyWalkerGFMC!(prob,nSteps,nThreads,equilibration_steps,initializer)
end

function fill_all_Buffers!(prob::AbstractGFMCProblem,nThreads)
    (;Walkers,Guiding_function_buffer,ψG) = prob
    batches = ChunkSplitters.chunks(eachindex(Walkers,Guiding_function_buffer), n = nThreads)
    Threads.@threads for (i_chunk,αinds) in enumerate(batches)
        for α in αinds
            GWFBuffer = Guiding_function_buffer[α]
            Walker = Walkers[α]
            compute_GWF_buffer!(GWFBuffer,ψG,Walker)
        end
    end
    return nothing
end
# function startManyWalkerGFMC(InitialState,Nwalkers,NSteps,G::AbstractPropagator;equilibration_steps = 0,outfile=nothing,pre_equilibration_steps=5*equilibration_steps)
#     startManyWalkerGFMC(InitialState,outfile,Nwalkers,NSteps,equilibration_steps,pre_equilibration_steps,G)
# end


function runGFMC!(prob::AbstractGFMCProblem,range::UnitRange,nThreads,reconfigure::Bool = true,save_energies::Bool = true,saveObservables::Bool = true)
    (;Walkers,weights,Guiding_function_buffer,reconfiguration_buffer,Observables,ψG,method) = prob
    reconfigurationTable = get_reconfigurationTable(Observables)
    iter = 0
    for i in range
        iter += 1
        # for (α,Config) in enumerate(Walkers)
        propagateWalkers!(Walkers,weights,Guiding_function_buffer,nThreads,ψG,method)
        save_energies && updateEnergies!(Observables,i,Walkers,weights,method)
        if reconfigure
            reconfigurationList = @view reconfigurationTable[:,i]
            reconfiguration!(Walkers,Guiding_function_buffer,reconfigurationList,reconfiguration_buffer,weights)
            # fill_all_Buffers!(prob,nThreads)

        end
        saveObservables && saveObservables!(Observables,i,Walkers)

        if iter%1000 == 0 # recompute buffers only occasionally to avoid accumulation of floating point errors 
            fill_all_Buffers!(prob,nThreads)
        end
    end
    return Observables
end
runGFMC!(prob::AbstractGFMCProblem,range::Integer,nThreads,reconfigure::Bool = true,save_energies::Bool = true,saveObservables::Bool = true) = runGFMC!(prob,1:range,nThreads,reconfigure,save_energies,saveObservables)
runGFMC!(prob::AbstractGFMCProblem,Nsteps::Int;nThreads = 2*Threads.nthreads(),reconfigure=true,save_energies=true,saveObservables = true) = runGFMC!(prob,1:Nsteps,nThreads,reconfigure,save_energies,saveObservables)

function propagateWalkers!(Walkers,weights,Guiding_function_buffer,nThreads,ψG,method::DiscreteTimeMethod)
    (;Λ,nBranch) = method
    w_avg_estimate = get_w_avg_estimate(method)
    w_avg_estimate⁻¹ = 1. / w_avg_estimate

    batches = ChunkSplitters.chunks(eachindex(Walkers), n = nThreads)

    Threads.@threads for (i_chunk,αinds) in enumerate(batches)
        for α in αinds
            GWFBuffer = Guiding_function_buffer[α]
            Walker = Walkers[α]
            w = 1.
            for step in 1:nBranch
                weightList = updateWeightList!(Walker,GWFBuffer,ψG,Λ)
                bx = sum(weightList)*w_avg_estimate⁻¹
                w *= bx
                performMarkovStep!(Walker)
            end
            weights[α] = w
            updateWeightList!(Walker,GWFBuffer,ψG,Λ)
        end
    end
end

function propagateWalkers!(Walkers,weights,Guiding_function_buffer,nThreads,ψG,method::ContinuousTimeMethod)
    (;Hxx,τ,w_avg_estimate) = method
    batches = ChunkSplitters.chunks(eachindex(Walkers), n = nThreads)

    Threads.@threads for (i_chunk,αinds) in enumerate(batches)
        for α in αinds
            Walker = Walkers[α]
            # compute_GWF_buffer!(GWFBuffer,ψG,Walker)
            GWFBuffer = Guiding_function_buffer[α]
            log_w = 0.
            weightList = updateWeightList!(Walker,GWFBuffer,ψG)
            H_xx = Hxx(Walker)
            el_x = H_xx + getLocalEnergy(weightList)
            βleft = τ
            while βleft > 0
                ξ = rand()
                # dτ = log(1-ξ)/(el_x - H_xx)
                dτ = min(βleft,log(1-ξ)/(el_x - H_xx))
                βleft -= dτ
                if isinf(βleft)
                    @info "" dτ el_x H_xx βleft maximum(weightList)
                    error("Infinite propagation time encountered. Check for too large values in guiding wavefunction or its Buffers!")
                end
                log_w += -dτ*el_x
                if βleft > 0 
                    last_move = performMarkovStep!(Walker)
                    post_move_update_GWF_buffer!(GWFBuffer,ψG,Walker,last_move)
                    updateWeightList!(Walker,GWFBuffer,ψG)

                    H_xx = Hxx(Walker)
                    el_x = H_xx + getLocalEnergy(weightList)
                end
            end
            w = exp(log_w - τ* w_avg_estimate)
            # w = exp(log_w)
            weights[α] = w
        end
    end
end

"""Performs an efficient reconfiguration of walkers. This reconfiguration will not remove walkers if they all have the same weight, which increases the efficiency as more walkers can contribute to the average.

Matteo Calandra Buonaura and Sandro Sorella
Phys. Rev. B 57, 11446 (1998)
"""
function reconfiguration!(Walkers::AbstractVector{<:AbstractWalker},Guiding_function_buffer,reconfigurationList,reconfiguration_buffer,weights)
    Nw = length(Walkers)
    reconfiguration_buffer = cumsum!(reconfiguration_buffer,weights) 
    wTotal = sum(weights)
    reconfiguration_buffer ./= wTotal
    for α in eachindex(Walkers)
        ξα = rand()
        zα = (α + ξα - 1)/Nw
        α´ = searchsortedfirst(reconfiguration_buffer,zα)
        reconfigurationList[α] = α´
        empty!(Walkers[α].weights)
        empty!(Walkers[α].moves)
    end
    minimizeReconfiguration!(reconfigurationList)
    for (α,α´) in enumerate(reconfigurationList)
        if α´ != α
            get_config(Walkers[α]) .= get_config(Walkers[α´])
            Guiding_function_buffer[α] = GWFBuffer_set_to!(Guiding_function_buffer[α],Guiding_function_buffer[α´])
        end
    end
end

function getLocalEnergyWalkers_before(weights,Walkers::AbstractVector{<:AbstractWalker},method::AbstractGFMCMethod)
    Nw = length(weights)
    num = 0.
    denom = 0.
    for α in eachindex(weights,Walkers)
        eloc = getLocalEnergy(Walkers[α],method)
        num += weights[α]*eloc
        denom += weights[α]
    end

    return num/denom
end

function updateEnergies!(Observables::AbstractGFMCObservables,i,Walkers::AbstractVector{<:AbstractWalker},weights,method)
    energies = get_energies(Observables)
    TotalWeights = get_TotalWeights(Observables)
    energies[i] = getLocalEnergyWalkers_before(weights,Walkers,method)
    TotalWeights[i] = mean(weights)
    return nothing
end

function saveObservables!(Observables::GFMCObservables,i,Walkers::AbstractVector{<:AbstractWalker})
    SaveConfigs = Observables.SaveConfigs
    for (α,Config) in enumerate(Walkers)
        SaveConfigs[:,:,α,i] .= get_config(Config)
    end
    return nothing
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

function random_init_walkers!(Walkers::AbstractVector{<:SpiderWebWalker},equilibration_steps,fraction=1.0)
    Nw = length(Walkers)

    Threads.@threads for α in eachindex(Walkers)[1:round(Int,fraction*Nw)]
        Walker = Walkers[α]
        
        for _ in 1:equilibration_steps
            movepos = Tuple(rand(Walker.Plaquette_positions))
            movesgn = rand(1:2)
            P_applicable(Walker.Config, movepos)[movesgn] || continue
            applyPlaquette!(Walker.Config, movepos[1], movepos[2], (1,-1)[movesgn])
        end
    end
end