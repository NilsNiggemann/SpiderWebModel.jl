struct GFMCObservables_StructureFac{T_high<:AbstractFloat,T_low<:AbstractFloat,FFTType<:SqFFT,T2} <: AbstractGFMCObservables
    TotalWeights::Vector{T_high}
    energies::Vector{T_high}
    FFTBuffers::Vector{FFTType}
    SqBuffers::Array{T_low,4}
    Sq_numerator::Array{T_high,3}
    Sq_denominator::Vector{T_high}
    Gnps::Matrix{T_high}
    reconfigurationTable::Matrix{Int}
    PopulationMatrix::Matrix{Int}
    outfile::T2
end

get_pMax(O::GFMCObservables_StructureFac) = size(O.SqBuffers,4)

function setup_Sq_Observables(InitConfig,NWalkers,NSteps,m_proj,outfile::Nothing)
    energies = zeros(NSteps)
    Lx,Ly = size(InitConfig)

    FFTBuffers = fetch.([Threads.@spawn SqFFT((Lx,Ly)) for i in 1:NWalkers])
    SqBuffers = zeros(Float32,Lx,Ly,NWalkers,m_proj)
    Sq_numerator = zeros(Lx,Ly,m_proj)
    Sq_denominator = zeros(Float64,m_proj)

    Gnps = zeros(Float64,NSteps,2m_proj)
    
    TotalWeights = zeros(NSteps)
    reconfigurationTable = zeros(Int,NWalkers,NSteps)
    PopulationMatrix = zeros(Int,NWalkers,m_proj)
    return GFMCObservables_StructureFac(TotalWeights,energies,FFTBuffers,SqBuffers,Sq_numerator,Sq_denominator,Gnps,reconfigurationTable,PopulationMatrix,outfile)
end
function setup_Sq_Observables(InitConfig,NWalkers,NSteps,m_proj,filename::String)
    Lx,Ly = size(InitConfig)
    h5open(filename,"cw") do file
        energies = createMMapArray(file,"energies",Float64,(NSteps,))
        TotalWeights = createMMapArray(file,"TotalWeights",Float64,(NSteps,))
        reconfigurationTable = createMMapArray(file,"reconfigurationTable",Int,(NWalkers,NSteps))
        FFTBuffers = fetch.([Threads.@spawn SqFFT((Lx,Ly)) for i in 1:NWalkers])
        SqBuffers = zeros(Float32,Lx,Ly,NWalkers,m_proj)
        Sq_numerator = createMMapArray(file,"Sq_numerator",Float64,(Lx,Ly,m_proj))
        Sq_denominator = createMMapArray(file,"Sq_denominator",Float64,(m_proj,))
        Gnps = createMMapArray(file,"Gnps",Float64,(NSteps,2m_proj))
        PopulationMatrix = zeros(Int,NWalkers,m_proj)
        return GFMCObservables_StructureFac(TotalWeights,energies,FFTBuffers,SqBuffers,Sq_numerator,Sq_denominator,Gnps,reconfigurationTable,PopulationMatrix,filename)
    end
end

function setup_Sq_problem(InitialState::StencilSpinConfig,method::AbstractGFMCMethod,Nwalkers::Integer,NSteps::Integer,m_proj,nThreads,ψG,outfile)
    setup = setup_many_walker_GFMC(InitialState,Nwalkers,nThreads)
    Guiding_function_buffer = allocate_GWF_buffers_threads(ψG,InitialState,Nwalkers)
    
    (;Walkers,weights,reconfiguration_buffer) = setup
    Observables = setup_Sq_Observables(InitialState,Nwalkers,NSteps,m_proj,outfile)

    return SpiderwebGFMCProblem(method,InitialState,ψG,Walkers,weights,Guiding_function_buffer,reconfiguration_buffer,Observables)
end

wrap_idx(n,pMax) = (n-1) % (pMax) + 1

function compute_Sq_Walkers!(SqBuffers::Array{T,4},Walkers,n,FFTBuffers::AbstractVector{<:AbstractObservable}) where T
    @assert axes(SqBuffers,3) == eachindex(FFTBuffers) == eachindex(Walkers)
    pMax = size(SqBuffers,4)

    Threads.@threads for α in eachindex(FFTBuffers)
        conf = get_config(Walkers[α])
        Sq_view = @view SqBuffers[:,:,α,wrap_idx(n,pMax)]

        ObsFunc! = FFTBuffers[α]
        Sq = obs(ObsFunc!)
        ObsFunc!(Sq,conf)
        Sq_view .= Sq
    end
    return
end

function updateGnp!(Gnp,TotalWeights,n)
    pMax = size(Gnp,2)
    for p in 1:pMax
        if n-p < 1
            Gnp[n,p] = 0
            continue
        end
        Gnp[n,p] = prod(@view TotalWeights[n-p:n])
    end
    return
end

function saveObservables!(Observables::GFMCObservables_StructureFac,n,Walkers::AbstractVector{<:SpiderWebWalker})

    (;reconfigurationTable,SqBuffers,Sq_numerator,Sq_denominator,PopulationMatrix) = Observables

    compute_Sq_Walkers!(Observables.SqBuffers,Walkers,n,Observables.FFTBuffers)
    updateGnp!(Observables.Gnps,Observables.TotalWeights,n)
    
    Nw = length(Walkers)
    
    pMax = get_pMax(Observables)
    # n <= pMax && return
    
    getPopulationMatrix!(PopulationMatrix,reconfigurationTable,n,pMax-1)
    
    Nw⁻¹ = 1/Nw

    m_values = 0:pMax-1
    Threads.@threads for m_index in eachindex(m_values)
        m = m_values[m_index]
        Gnp = Observables.Gnps[n,1+2m]
        Sq_denominator[m_index] += Gnp
        # Sq_denominator[m_index] += Gnp*Nw
        @views for α in 1:Nw
            mult = PopulationMatrix[α,m_index]
            mult == 0 && continue
            mult *= Nw⁻¹
            O = SqBuffers[:,:,α,wrap_idx(n-m,pMax)]
            @. Sq_numerator[:,:,m_index] += O*Gnp*mult
        end
    end

end

function get_alpha_prime_table(reconfigurationTable::Matrix{Int},n::Integer,pMax::Integer)
    Nw = size(reconfigurationTable,1)
    alpha_prime_table = zeros(Int,Nw)
    for α in 1:Nw
        α´ = α
        for i in 1:pMax
            α´ = reconfigurationTable[α´,n-i]
        end
        alpha_prime_table[α] = α´
    end
    return alpha_prime_table
end
function measure_Sq_GFMC(InitialState::StencilSpinConfig,method::AbstractGFMCMethod,Nwalkers::Integer,nSteps::Integer,mProj,ψG; equilibration_steps = 0, pre_equilibration_steps = equilibration_steps ÷ 5, scatter_fraction = 0.8,initializer = UnguidedWalkInitializer(pre_equilibration_steps,scatter_fraction),nThreads=2*Threads.nthreads(),outfile = nothing,kwargs...)
    prob = setup_Sq_problem(InitialState,method,Nwalkers,nSteps,mProj,nThreads,ψG,outfile)
    startManyWalkerGFMC!(prob,nThreads,equilibration_steps,initializer)
    normalize_numerator!(prob.Observables)
    return prob.Observables    
end

function normalize_numerator!(Observables::GFMCObservables_StructureFac)
    numerator = Observables.Sq_numerator
    denominator = Observables.Sq_denominator

    for i in eachindex(Observables.Sq_denominator)
        @. numerator[:,:,i] ./= denominator[i]
        # denominator[i] = 1
    end
    return
end