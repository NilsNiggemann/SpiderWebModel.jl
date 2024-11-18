"""Abstract supertype for a walker in a GFMC simulation. A walker needs to have a method `get_config` which returns the configuration (for example the spin configuration), and a method `get_weights` which returns a vector of weights for each possible move."""
abstract type AbstractWalker end
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

abstract type AbstractGFMCMethod end

struct DiscreteTimeMethod <: AbstractGFMCMethod 
    Λ::Float64
    nBranch::Int
    w_avg_estimate::Float64
end

DiscreteTimeMethod(;Λ=1.,nBranch,w_avg_estimate=1.) = DiscreteTimeMethod(Λ,nBranch,w_avg_estimate)

struct ContinuousTimeMethod{F2} <: AbstractGFMCMethod 
    τ::Float64
    nBranch::Int
    w_avg_estimate::Float64
    Hxx::F2
end
ContinuousTimeMethod(τ,nBranch::Integer,w_avg_estimate=1.,Hxx=Hxx_zero()) = ContinuousTimeMethod(float(τ),nBranch,float(w_avg_estimate),Hxx)
ContinuousTimeMethod(τ;nBranch=1,w_avg_estimate=1.,Hxx=Hxx_zero()) = ContinuousTimeMethod(τ,nBranch,float(w_avg_estimate),Hxx)

abstract type AbstractGFMCProblem end
struct SpiderwebGFMCProblem{MethodType<:AbstractGFMCMethod,T<:AbstractFloat,BufferType,C,F,W,O} <: AbstractGFMCProblem
    method::MethodType
    InitialState::C
    ψG::F
    Walkers::Vector{W}
    weights::Vector{T}
    Guiding_function_buffer::BufferType
    reconfiguration_buffer::Vector{T}
    Observables::O
end

function get_Guiding_function_buffer(problem::SpiderwebGFMCProblem)
    return problem.Guiding_function_buffer
end

abstract type AbstractGFMCObservables end
struct GFMCObservables{DT<:AbstractFloat,T,T2} <: AbstractGFMCObservables
    energies::Vector{DT}
    SaveConfigs::T
    TotalWeights::Vector{DT}
    reconfigurationTable::Matrix{Int}
    outfile::T2
end

function _setup_GFMC_problem(InitialState::StencilSpinConfig,method::AbstractGFMCMethod,Nwalkers::Integer,NSteps::Integer,ψG,outfile)
    setup = setup_many_walker_GFMC(InitialState,Nwalkers)
    Guiding_function_buffer = allocate_GWF_buffers_threads(ψG,InitialState)
    
    (;Walkers,weights,reconfiguration_buffer) = setup
    ObsSetup = setupObservables(InitialState,Nwalkers,NSteps,outfile)
    (;energies,SaveConfigs,TotalWeights,reconfigurationTable) = ObsSetup

    Observables = GFMCObservables(energies,SaveConfigs,TotalWeights,reconfigurationTable,outfile)

    return SpiderwebGFMCProblem(method,InitialState,ψG,Walkers,weights,Guiding_function_buffer,reconfiguration_buffer,Observables)
end

function setup_GFMC_problem(InitialState::ConfType, method::AbstractGFMCMethod, Nwalkers::Integer, NSteps::Integer, ψG;
    outfile = nothing
    ) where {ConfType <: StencilSpinConfig}
    _setup_GFMC_problem(InitialState,method,Nwalkers,NSteps,ψG,outfile)
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

add_elementwise!(x::AbstractArray,y) = (x .+= y)
add_elementwise!(x::Number,y::Number) = x + y

mult_elementwise!(x::Number,y::Number) = x * y
function mult_elementwise!(x::AbstractArray,y::Number)
    
    for i in eachindex(x)
        x[i] *= y
    end
    return x
end

divide_elementwise!(x::Number,y::Number) = x / y
divide_elementwise!(x::AbstractArray,y) = (x ./= y)

function getObs(Gnp,AllConfigs,reconfigurationTable,ObsFunc,m::Integer=size(Gnp,2)÷2)
    N = lastindex(AllConfigs,4)
    exampleConf = @view AllConfigs[:,:,begin,begin]
    Obs = ObsFunc(exampleConf)
    num = float(zero(Obs))
    denom = 0.
    GnO = float(zero(Obs))
    # fill!(Obs,zero(eltype(Obs)))

    Nw = size(reconfigurationTable,1)
    p = size(Gnp,2)
    # surviving_walker_mapping_list = zeros(Int,Nw)
    WalkerMultiplicities = zeros(Int,Nw)

    for n in m+1:N
        Gn = Gnp[n,p]
        # denom += Gn*Nw
        denom += Gn*Nw
        # reconfigList = @view reconfigurationTable[:,n-m]
        # surviving_walker_mapping!(surviving_walker_mapping_list,reconfigList)
        WalkerMultiplicities .= 0
        for α in 1:Nw
            α´ = α
            for i_m in 1:m
                α´ = reconfigurationTable[α´,n-i_m]
            end
            WalkerMultiplicities[α´] += 1
        end

        for α in 1:Nw
            mult = WalkerMultiplicities[α]
            mult == 0 && continue
            conf = @view AllConfigs[:,:,α,n-m]
            O = ObsFunc(conf)
            # surviving_index = surviving_walker_mapping_list[α´]
            # O = ObsFunc(AllConfigs[n-m][surviving_index])
            GnO = _set_to!(GnO,O)
            GnO = mult_elementwise!(GnO,Gn*mult)
            # @. num += Gn*O
            num = add_elementwise!(num,GnO)
        end
    end
    return divide_elementwise!(num,denom)
end
function getObs(result,ObsFunc,p::Integer)
    Gnp = precomputeNormalizedAccWeight(result.TotalWeights,1,p)
    return getObs(Gnp,result.SaveConfigs,result.reconfigurationTable,ObsFunc,p÷2)
end

function getObs(Gnp,AllConfigs,reconfigurationTable,ObsFunc::AbstractObservable,m_values::ABSTRACTCOLLECTION)
    N = lastindex(AllConfigs,4)

    pMax = maximum(m_values)

    Obs = obs(ObsFunc)
    num_m = [zeros(size(Obs)) for _ in m_values]
    denom = 0.

    Nw = size(reconfigurationTable,1)
    p = size(Gnp,2)
    WalkerMultiplicities = zeros(Int,Nw)
    ObsBuffer = [similar(Obs) for α in 1:Nw, n in 1:(pMax)]
    
    wrap_idx(n) = (n-1) % (pMax) + 1
    obsBuffer(α,n) = ObsBuffer[α,wrap_idx(n)]

    _fill_obs_buffer!(ObsBuffer,1:pMax,ObsFunc,AllConfigs,pMax)

    for n in pMax+1:N
        Gn = Gnp[n,p]
        denom += Gn*Nw
        _fill_obs_buffer!(ObsBuffer,n-1,ObsFunc,AllConfigs,pMax)
        for (i_m,m) in enumerate(m_values)
            WalkerMultiplicities .= 0
            for α in 1:Nw
                α´ = α
                for i_m in 1:m
                    α´ = reconfigurationTable[α´,n-i_m]
                end
                WalkerMultiplicities[α´] += 1
            end


            for α in 1:Nw
                mult = WalkerMultiplicities[α]
                mult == 0 && continue

                O = obsBuffer(α,n-m)
                @. num_m[i_m] += O*Gn*mult
            end
        end
    end
    for i in eachindex(num_m)
        divide_elementwise!(num_m[i],denom)
    end
    return num_m
end

function _fill_obs_buffer!(ObsBuffer,nRange,ObsFunc!,AllConfigs,pMax)
    wrap_idx(n) = (n-1) % (pMax) + 1
    obsBuffer(α,n) = ObsBuffer[α,wrap_idx(n)]

    for n in nRange, α in axes(AllConfigs,3)
        conf = @view AllConfigs[:,:,α,n]
        ObsFunc!(obsBuffer(α,n),conf)
    end
    return
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

function setup_many_walker_GFMC(InitialState::ConfType,Nwalkers::Integer) where {ConfType <: StencilSpinConfig}
    plaquettePositions = collect(plaquetteIterator(InitialState))
    Walkers = Vector{SpiderWebWalker{ConfType}}(undef,Nwalkers)
    Threads.@threads for α in eachindex(Walkers)
        Walkers[α] = spiderWebWalker(InitialState,plaquettePositions)
    end
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

function saveParameters(filename::String,Λ,equilibration_steps,nBranch,ψG,w_avg_estimate;kwargs...)
    h5open(filename,"cw") do file
        file["Λ"] = Λ
        file["equilibration_steps"] = equilibration_steps
        file["nBranch"] = nBranch
        file["w_avg_estimate"] = w_avg_estimate
        for (k,v) in kwargs
            file[String(k)] = v
        end
        saveVariationalParameter(file,ψG)
    end
end
function saveParameters(filename::String,equilibration_steps,method::DiscreteTimeMethod,ψG)
    (;Λ,nBranch,w_avg_estimate) = method
    saveParameters(filename,Λ,equilibration_steps,nBranch,ψG,w_avg_estimate)
end
function saveParameters(filename::String,equilibration_steps,method::ContinuousTimeMethod,ψG)
    (;τ,nBranch,w_avg_estimate) = method
    saveParameters(filename,Inf,equilibration_steps,nBranch,ψG,w_avg_estimate;τ=τ)
end

saveParameters(::Nothing,args...) = nothing

function saveVariationalParameter(file::HDF5.File,ψG)
    pars = variational_parameters(ψG)
    funcName = guidingfunc_name(ψG)
    for (key,val) in pars
        file[string(funcName,"/",key)] = val
    end
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
    
    (;Guiding_function_buffer,Walkers,weights,reconfiguration_buffer,Observables,method,ψG) = prob
    (;outfile,reconfigurationTable) = Observables

    saveParameters(outfile,equilibration_steps,method,ψG)

    initialize!(Walkers,initializer)

    #fill buffers for available steps and weights
    reconfigurationList = @view reconfigurationTable[:,1]
    for _ in 1:equilibration_steps
        propagateWalkers!(Walkers,weights,Guiding_function_buffer,ψG,method)
        reconfiguration!(Walkers,reconfigurationList,reconfiguration_buffer,weights)
    end

    return prob,Observables
end

function startManyWalkerGFMC!(prob::AbstractGFMCProblem,equilibration_steps::Int,initializer = UnguidedWalkInitializer(equilibration_steps ÷ 5,0.8))
    initializeGFMC!(prob,equilibration_steps,initializer)
    runGFMC!(prob)
end

# function startManyWalkerGFMC(InitialState,Nwalkers,NSteps,G::AbstractPropagator;equilibration_steps = 0,outfile=nothing,pre_equilibration_steps=5*equilibration_steps)
#     startManyWalkerGFMC(InitialState,outfile,Nwalkers,NSteps,equilibration_steps,pre_equilibration_steps,G)
# end

function startManyWalkerGFMC(InitialState::StencilSpinConfig,method::AbstractGFMCMethod,Nwalkers::Integer,nSteps::Integer,ψG; equilibration_steps = 0, pre_equilibration_steps = equilibration_steps ÷ 5, scatter_fraction = 0.8,initializer = UnguidedWalkInitializer(pre_equilibration_steps,scatter_fraction),kwargs...)
    prob = setup_GFMC_problem(InitialState,method,Nwalkers,nSteps,ψG;kwargs...)
    startManyWalkerGFMC!(prob,equilibration_steps,initializer)
end

function runGFMC!(prob::AbstractGFMCProblem,range,reconfigure::Bool=true)
    (;Walkers,weights,Guiding_function_buffer,reconfiguration_buffer,Observables,ψG,method) = prob
    (;energies,SaveConfigs,outfile,TotalWeights, reconfigurationTable) = Observables
    
    for i in range
        # for (α,Config) in enumerate(Walkers)
        propagateWalkers!(Walkers,weights,Guiding_function_buffer,ψG,method)

        energies[i] = getLocalEnergyWalkers_before(weights,Walkers,method)
        TotalWeights[i] = mean(weights)
        
        if reconfigure
            reconfigurationList = @view reconfigurationTable[:,i]
            reconfiguration!(Walkers,reconfigurationList,reconfiguration_buffer,weights)
        end
        saveConfigs!(SaveConfigs,i,Walkers)
    end
    return Observables
end
runGFMC!(prob::AbstractGFMCProblem;reconfigure=true) = runGFMC!(prob,eachindex(prob.Observables.TotalWeights),reconfigure)
runGFMC!(prob::AbstractGFMCProblem,Nsteps::Int;reconfigure=true) = runGFMC!(prob,eachindex(prob.Observables.TotalWeights)[1:Nsteps],reconfigure)

function propagateWalkers!(Walkers,weights,Guiding_function_buffer,ψG,method::DiscreteTimeMethod)
    (;Λ,nBranch,w_avg_estimate) = method

    w_avg_estimate⁻¹ = 1. / w_avg_estimate

    batches = ChunkSplitters.chunks(eachindex(Walkers), n = Threads.nthreads())

    Threads.@threads for (i_chunk,αinds) in enumerate(batches)
        GWFBuffer = Guiding_function_buffer[i_chunk]
        for α in αinds
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

function propagateWalkers!(Walkers,weights,Guiding_function_buffer,ψG,method::ContinuousTimeMethod)
    (;Hxx,nBranch,τ,w_avg_estimate) = method
    
    batches = ChunkSplitters.chunks(eachindex(Walkers), n = Threads.nthreads())

    Threads.@threads for (i_chunk,αinds) in enumerate(batches)
        GWFBuffer = Guiding_function_buffer[i_chunk]
        for α in αinds
            Walker = Walkers[α]
            compute_GWF_buffer!(GWFBuffer,ψG,Walker)
            log_w = 0.
            weightList = updateWeightList!(Walker,GWFBuffer,ψG)
            H_xx = Hxx(Walker)
            el_x = H_xx + getLocalEnergy(weightList)
            for _ in 1:nBranch
                βleft = τ
                while βleft > 0
                    ξ = rand()
                    # dτ = log(1-ξ)/(el_x - H_xx)
                    dτ = min(βleft,log(1-ξ)/(el_x - H_xx))
                    βleft -= dτ
                    log_w += -dτ*el_x
                    if βleft > 0 
                        last_move = performMarkovStep!(Walker)
                        post_move_update_GWF_buffer!(GWFBuffer,ψG,Walker,last_move)
                        updateWeightList!(Walker,GWFBuffer,ψG)

                        H_xx = Hxx(Walker)
                        el_x = H_xx + getLocalEnergy(weightList)
                    end
                end
            end
            w = exp(log_w - nBranch*τ* w_avg_estimate)
            # w = exp(log_w)
            weights[α] = w
        end
    end
end

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
        empty!(Walkers[α].weights)
        empty!(Walkers[α].moves)
    end
    minimizeReconfiguration!(reconfigurationList)
    for (α,α´) in enumerate(reconfigurationList)
        if α´ != α
            get_config(Walkers[α]) .= get_config(Walkers[α´])
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


function getImagTimeCorr(Gnp,reconfigurationTable,ObsFunc::T,mtau=size(Gnp,2)÷4, m=size(Gnp,2)÷2) where {T}
    N = lastindex(reconfigurationTable,2)
    Obs = [float(zero(ObsFunc(1,2m))) for i in 1:mtau]
    O0 = float(zero(ObsFunc(1,2m)))
    # num = zero(Obs)
    denom = 0.

    Nw = size(reconfigurationTable,1)
    p = size(Gnp,2)

    BranchingMatrix = zeros(Int,Nw,m+1)

    WalkerMultiplicities = zeros(Int,Nw,m+1)
    for n in m+1:N
        Gn = Gnp[n,p]
        denom += Gn*Nw
        
        getBranchingMatrix!(BranchingMatrix,WalkerMultiplicities,reconfigurationTable,n,m)
        for α in 1:Nw
            # O0 = ObsFunc(BranchingMatrix[α,m],n-m)
            O0 = _set_to!(O0,ObsFunc(BranchingMatrix[α,m],n-m))
            for ntau in 0:mtau-1
                mult = WalkerMultiplicities[α,m-ntau]
                mult == 0 && continue
                Otau = ObsFunc(BranchingMatrix[α,m-ntau],n-m+ntau)
                Obs[ntau+1] = _add_numerator!(Obs[ntau+1],Gn*mult,O0,Otau)
                # Obs[ntau+1] += Gn*mult*O0*Otau
            end
        end
    end
    
    for i in eachindex(Obs)
        Obs[i] /= denom
    end
    return (Obs)
end
_add_numerator!(Obsn::AbstractArray,Gnmult,O0::AbstractArray,Otau::AbstractArray) = (@. Obsn += Gnmult*O0*Otau)
_add_numerator!(Obsn::Number,Gnmult,O0::Number,Otau::Number) = (Obsn += Gnmult*O0*Otau)
_set_to!(a,b) = (a=b)
_set_to!(a::AbstractArray,b::AbstractArray) = (a.=b)

function getBranchingMatrix!(BranchingMatrix::AbstractMatrix,PopulationMatrix,reconfigurationTable::AbstractMatrix,n,projectionLength)
    # BranchingMatrix[:,begin] .= @view reconfigurationTable[:,begin]
    PopulationMatrix .= 0 
    for α in axes(reconfigurationTable,1)
        α´ = α
        for i_m in 0:projectionLength
            α´ = reconfigurationTable[α´,n-i_m]
            # println((; α,α´,i_m))
            BranchingMatrix[α,i_m+1] = α´
            PopulationMatrix[α´,i_m+1] += 1
        end
    end
    return (;BranchingMatrix,PopulationMatrix)
end

function getBranchingMatrix(reconfigurationTable,n,projectionLength) 
    BranchingMatrix = zeros(Int,size(reconfigurationTable,1),projectionLength+1)
    PopulationMatrix = zeros(Int,size(reconfigurationTable,1),projectionLength+1)
    getBranchingMatrix!(BranchingMatrix,PopulationMatrix,reconfigurationTable,n,projectionLength)
end

function constructSWF_operator(AllConfigs,OpFunc::T) where T
    function Obsfunc(α,n)
        conf = @view AllConfigs[:,:,α,n]
        return OpFunc(conf)
    end
end