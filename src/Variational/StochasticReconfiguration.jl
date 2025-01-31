import LinearMaps, IterativeSolvers
struct QuantumMetric{T,T2<:AbstractMatrix{T}} <: LinearMaps.LinearMap{T}
    O_km::T2
    size::Dims{2}
end
Base.size(S::QuantumMetric) = S.size
LinearAlgebra.issymmetric(S::QuantumMetric) = true
LinearAlgebra.ishermitian(S::QuantumMetric) = true
LinearAlgebra.isposdef(S::QuantumMetric) = true
Base.:(==)(S1::QuantumMetric, S2::QuantumMetric) = S1.O_km == S2.O_km

function QuantumMetric(O_km)
    O_k_avg = reshape(mean(O_km,dims=1),size(O_km,2))
    @boundscheck checkbounds(O_km,1,lastindex(O_k_avg))

    LoopVectorization.@turbo for i in eachindex(O_k_avg)
        for j in axes(O_km,1)
            O_km[j,i] -= O_k_avg[i]
        end
        # O_km[:,i] .-= O_k_avg[i]
    end
    Nparams = size(O_km,2)
    return QuantumMetric(O_km,Dims((Nparams,Nparams)))
end

function LinearMaps._unsafe_mul!(y, S::QuantumMetric, v::AbstractVector)
    N_MC = size(S.O_km,2)
    O = S.O_km
    
    zk1 = O' * (O * v) / N_MC

    y .= zk1 .+1e-6.*v# .- zk2

end

function reconf_obs(InitialState::StencilSpinConfig,method::AbstractGFMCMethod,configs,ψG,
    Ok_i = zeros(Float32,length(configs),Nparams),
    E_i = zeros(Float32,length(configs)),
    inequivParams=eachindex(get_params(ψG))
    )
    plaqs = collect(plaquetteIterator(InitialState))
    Guiding_function_buffers = allocate_GWF_buffers_threads(ψG,InitialState,Threads.nthreads())

    Nparams = length(inequivParams)


    WorkChunks = ChunkSplitters.chunks(eachindex(IndexLinear(),configs),n=length(Guiding_function_buffers))
    
    Threads.@threads for (ichunk,chunkinds) in enumerate(WorkChunks)
        Guiding_function_buffer = Guiding_function_buffers[ichunk]
        Walker = spiderWebWalker(InitialState,plaqs)

        for (i,iconf) in enumerate(chunkinds)
            get_config(Walker) .= configs[iconf]
            compute_GWF_buffer!(Guiding_function_buffer,ψG,Walker)
            updateWeightList!(Walker,Guiding_function_buffer,ψG)
            elocal = getLocalEnergy(Walker,method)
            @inbounds for (ik,k) in enumerate(inequivParams)
                O_xk = getOx_k(ψG,Walker,k)
                Ok_i[iconf,ik] = O_xk
            end
            E_i[iconf] = elocal
        end

    end

    N = length(configs)

    EL_avg = mean(E_i)
    EL_err = std(E_i)

    return (;EL_avg,EL_err,E_i,Ok_i)
end


getSOperator(O_ki,::Type{QuantumMetric}) = QuantumMetric(O_ki)

abstract type AbstractSRSolver end

struct ExplicitSRSolver <: AbstractSRSolver end
struct IterativeSRSolver{T} <: AbstractSRSolver 
    kwargs::T
end
IterativeSRSolver(;kwargs...) = IterativeSRSolver(kwargs)

function stochastic_reconfiguration_step(E_i::AbstractVector,Ok_i::AbstractMatrix,solver::ExplicitSRSolver = ExplicitSRSolver())
    println("computing cov")

    @time S = cov(Ok_i)
    S += 1e-6I
    GC.gc()

    @time F = StatsBase.cov(Ok_i,E_i)
    
    @time SChol = cholesky!(Symmetric(S),check = false)
    if !issuccess(SChol)
        @warn "Cholesky factorization failed"
        SChol = cholesky!(Symmetric(S+1e-3I),check = false)
        !issuccess(SChol) && return zeros(length(F))
    end
    @time δα = SChol \ F
    for i in axes(δα,1)
        δα[i] = -δα[i]
    end
    return δα
end

function stochastic_reconfiguration_step(E_i::AbstractVector,Ok_i::AbstractMatrix,solver::IterativeSRSolver;kwargs...)
    S = QuantumMetric(Ok_i)
    F = reshape(StatsBase.cov(Ok_i,E_i),size(Ok_i,2))
    res = -IterativeSolvers.cg(S,F;solver.kwargs...,kwargs...)
    return res
end
stochastic_reconfiguration_step(E_i,Ok_i,::AbstractSRSolver) = error("solver not implemented")

function _stochastic_reconfiguration(InitialState,method::AbstractGFMCMethod,solver::AbstractSRSolver,NSteps::AbstractVector,GWF::SymmetryReducedWaveFunction,n,dt::AbstractVector,equilibration_steps=1000,rel_tolerance=1e-2,Nwalkers = Threads.nthreads(),nThreads=n_threads_default(Nwalkers),outfile=nothing,initializer = UnguidedWalkInitializer(equilibration_steps,0.8);verbose=true,report_steps=1,reset = false)
    
    (;psi,indicesMapping,uniqueInds) = GWF
    ψG = deepcopy(psi)
    params = get_params(ψG)

    convergedSteps = 0

    normDelta = Inf
    # _, = getDistReduction(InitialState,ψG)
    maxNSteps = maximum(NSteps)
    prob = setup_GFMC_problem(InitialState,method,Nwalkers,maxNSteps,nThreads,ψG)
    initializeGFMC!(prob,equilibration_steps,initializer)
    equilibrate!(prob,equilibration_steps;nThreads,reconfigure=false)
    results = get_stoch_rec_Observables(n,ψG,outfile)

    ind = 0
    Nparams = length(uniqueInds)

    All_Ok_i = zeros(Float32,maxNSteps*Nwalkers,Nparams)
    All_E_i = zeros(Float32,maxNSteps*Nwalkers)

    Guiding_function_buffers = allocate_GWF_buffers_threads(ψG,InitialState,Nwalkers)
    
    for i in 1:n
        range = eachindex(prob.Observables.TotalWeights)[1:NSteps[i]]
        # for w in prob.Walkers
        #     get_config(w) .= InitialState
        # end
        fill_all_Buffers!(prob,nThreads)
        if reset
            initializeGFMC!(prob,equilibration_steps,initializer)
            equilibrate!(prob,equilibration_steps;nThreads,reconfigure=false)
        end
        # @time res = runGFMC!(prob,range,reconfigure)
        res = runGFMC!(prob,range,nThreads,false)

        resSlice = @view res.SaveConfigs[:,:,:,range]
        confs = eachslice( resSlice,dims=(3,4))
        Ok_i = @view All_Ok_i[1:NSteps[i]*Nwalkers,:]
        E_i = @view All_E_i[1:NSteps[i]*Nwalkers]

        observables = reconf_obs(InitialState,method,confs,ψG,Ok_i,E_i,uniqueInds)
        # (;EL_avg,EL_err ) = observables
        # return observables
        δα = stochastic_reconfiguration_step(observables.E_i,observables.Ok_i,solver)
        normDelta = norm(δα)/norm(params)

        add_reconstructedFullParams!(ψG,indicesMapping,δα .*dt[i])
        
        E0 = mean(E_i)
        results.E0[i] = E0
        ΔE0 = std(E_i)
        results.ΔE[i] = ΔE0
        
        # push!(params_steps, copy(params))
        selectdim(results.params_steps,arraydim(results.params_steps),i) .= params

        # α = get_alpha_i(ψG)
        w_avg = mean(res.TotalWeights[range])
        if verbose && i % report_steps == 0
            @info "optimization step $i" dt[i] numConfigs = NSteps[i]*Nwalkers length(uniqueInds) length(params) "max(|α|)" = maximum(abs,params) "||δα||" = normDelta δα[1] E0 ΔE0 = ΔE0 convergedSteps
        end
        ind = i
        if normDelta < rel_tolerance
            convergedSteps += 1

            convergedSteps > 4 && break
        else
            convergedSteps = 0
        end

        # fill_all_Buffers!(prob,nThreads)
        Accessors.@reset prob.Guiding_function_buffer = allocate_GWF_buffers_threads(ψG,InitialState,Nwalkers) # recompute GWF buffers in case they change
    end

    params = selectdim(results.params_steps,arraydim(results.params_steps),ind)

    return (;results...,params)
end
arraydim(a::AbstractArray{T,N}) where {T,N} = N

function stochastic_reconfiguration(InitialState,method::AbstractGFMCMethod,NSteps,ψG,n,dt,solver::AbstractSRSolver = IterativeSRSolver();
    equilibration_steps=1000,
    rel_tolerance=1e-2,
    Nwalkers = 6,
    outfile=nothing,
    pre_equilibration_steps=5*equilibration_steps,
    scatter_fraction=0.8,
    nThreads=n_threads_default(Nwalkers),
    initializer = UnguidedWalkInitializer(pre_equilibration_steps,scatter_fraction),
    kwargs...)
    
    NStepsVec = makeVec(NSteps,n)
    dtVec = makeVec(dt,n)
    GWF = _default_symmetry(InitialState,ψG)

    return _stochastic_reconfiguration(InitialState,method,solver,NStepsVec,GWF,n,dtVec,equilibration_steps,rel_tolerance,Nwalkers,nThreads,outfile,initializer;kwargs...)
end

_default_symmetry(InitialState,ψG::AbstractGuidingFunction) = getDistReduction(InitialState,ψG)
_default_symmetry(InitialState,ψG::SymmetryReducedWaveFunction) = ψG

makeVec(x::AbstractVector,len) = x
makeVec(x::Number,len) = fill(x,len)
makeVec(f::Function,len) = f.(1:len)

function get_stoch_rec_Observables(Nsteps,ψG,::Nothing)
    E0 = zeros(Nsteps)
    ΔE = zeros(Nsteps)
    params = zeros(size(get_params(ψG))...,Nsteps)
    return (;E0 = E0,ΔE = ΔE,params_steps = params)
end

function get_stoch_rec_Observables(Nsteps,ψG,outfile::AbstractString)
    h5open(outfile,"cw") do file
        E0 = createMMapArray(file,"E0",Float64,(Nsteps,))
        ΔE = createMMapArray(file,"ΔE",Float64,(Nsteps,))
        params = createMMapArray(file,"params_steps",Float32,(size(get_params(ψG))...,Nsteps))
        E0 .= 0
        ΔE .= 0
        params .= 0
        return (;E0 = E0,ΔE = ΔE,params_steps = params)
    end
end