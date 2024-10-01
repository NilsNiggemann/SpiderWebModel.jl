import LinearMaps, IterativeSolvers
struct QuantumMetric{T} <: LinearMaps.LinearMap{T}
    O_km::Matrix{T}
    size::Dims{2}
end
Base.size(S::QuantumMetric) = S.size
LinearAlgebra.issymmetric(S::QuantumMetric) = true
LinearAlgebra.ishermitian(S::QuantumMetric) = true
LinearAlgebra.isposdef(S::QuantumMetric) = true
Base.:(==)(S1::QuantumMetric, S2::QuantumMetric) = S1.O_km == S2.O_km

function QuantumMetric(O_km)
    O_k_avg = reshape(mean(O_km,dims=1),size(O_km,2))
    for i in eachindex(O_k_avg)
        O_km[:,i] .-= O_k_avg[i]
    end
    Nparams = size(O_km,2)
    return QuantumMetric(O_km,Dims((Nparams,Nparams)))
end

function LinearMaps._unsafe_mul!(y, S::QuantumMetric, v::AbstractVector)
    N_MC = size(S.O_km,2)
    O = S.O_km
    # for k in axes(S,1)
    #     z_k1 = zero(eltype(v))
    #     z_k2 = zero(eltype(v))
    #     for m in axes(S,1)
    #         z_k2 +=  S.O_k_avg[m] * v[m]
    #     end
    #     z_k2 *= S.O_k_avg[k]

    #     for μ in axes(S,2)
    #         z_k1_m = zero(eltype(v))
    #         for m in axes(S,1)
    #             z_k1_m += O[m,μ] * v[m]
    #         end
    #         z_k1 += O[k,μ] * z_k1_m / N_MC
    #     end
    # end

    zk1 = O' * (O * v) / N_MC
    
    # zk2 = (Ō' * v) .* Ō
    
    y .= zk1 .+1e-6.*v# .- zk2

end

function getOx_k(ψG,Walker::SpiderWebWalker,k)
    return getOx_k(ψG,get_config(Walker),k)
end

function getOx_k(ψG::Union{FullVariationalGuidingFunction,LocalPlaquetteGuidingFunction},Walker::SpiderWebWalker,k)
    return getOx_k_plaqs(ψG,Walker.n_x,k)
end

function getOx_k(ψG::OrderGuidingFunction,Walker::SpiderWebWalker,k)
    α = get_alpha_i(ψG)

    if k in eachindex(α)
        return Walker.n_x[k]
    end
    # m = get_m_i(ψG)

    Ok = get_config(Walker)[k-lastindex(α)]

    return Ok
end

function getOx_k_plaqs(ψG::FullVariationalGuidingFunction,n::AbstractArray,k)
    par = ψG.params

    α = get_alpha_i(ψG)

    if k in eachindex(α)
        return n[k]
    end
    β = get_beta_ij(ψG)

    (i,j) = Tuple(CartesianIndices(β)[k-lastindex(α)])

    Ok = n[i]*n[j]

    return Ok
end
getOx_k_plaqs(::LocalPlaquetteGuidingFunction,n::AbstractArray,k) = n[k]

function getOx_k(ψG::RBM,x::AbstractMatrix,k)
    par = ψG.params

    paramtype = getParameterType(ψG,k)

    if paramtype == 1 
        return x[k]
    end
    bj = get_b_j(ψG)
    Wij = get_W_ij(ψG)
    
    if paramtype == 2
        j = k - ψG.N
        θj = get_theta_j(x,j,bj,Wij)
        return tanh(θj)
    elseif paramtype == 3
        i,j = Tuple(CartesianIndices(Wij)[k-ψG.N-length(bj)])
        θj = get_theta_j(x,j,bj,Wij)
        return x[i] * tanh(θj)
    end
    return Ok
end

function nearbyInt(x1,x2,x_size)
    x_rsize = 1.0 / x_size

    dx = x1 - x2
    dx -= x_size * round(Int,dx * x_rsize)
end

isperiodic(S::Stencils.StencilArray) = Stencils.boundary(S) == Stencils.Wrap()
isperiodic(S::StencilSpinConfig) = isperiodic(parent(S))

function getDistReduction(S,ψG::AbstractGuidingFunction) 
    AllDists = Dict{SVector{2,Int},Int}()
    indicesMapping = collect(eachindex(ψG.params))
    uniqueInds = collect(indicesMapping)
    return (;AllDists,indicesMapping,uniqueInds)
end

function getDistReduction(S,ψG::FullVariationalGuidingFunction)
    
    AllDists = Dict{SVector{2,Int},Int}()
    if !isperiodic(S)
        indicesMapping = collect(eachindex(ψG.params))
        uniqueInds = collect(indicesMapping)
        return (;AllDists,indicesMapping,uniqueInds)
    end

    α = get_alpha_i(ψG)
    Allplaqs = collect(plaquetteIterator(S))

    # AllDists = Dict{Tuple{Int,Int},Int}()
    
    betaIndex = lastindex(α)
    indicesMapping = ones(Int,betaIndex)
    uniqueInds = [1]
    # indicesMapping = collect(eachindex(α))
    # uniqueInds = collect(eachindex(α))
    Lx,Ly = size(S)

    for (i,ri) in enumerate(Allplaqs)
        for (j,rj) in enumerate(Allplaqs)
            Rij = getReducedDist(ri,rj,Lx,Ly)
            # Rij = SVector(0,0)
            # Rij = (i,j)
            betaIndex += 1
            if Rij ∉ keys(AllDists)
                uniqueInds = push!(uniqueInds,betaIndex)
                AllDists[Rij] = length(uniqueInds)
            end
            push!(indicesMapping,AllDists[Rij])
        end
    end

    return (;AllDists,indicesMapping,uniqueInds)

end
function getReducedDist(ri,rj,Lx,Ly) 
    rpr = abs.(SVector(nearbyInt.(ri, rj,(Lx,Ly))))
    return sort(rpr)
end

centralPos(Lx,Ly) = ((Lx+1)//2,(Ly+1)//2)
centralPos(S::AbstractMatrix) = centralPos(size(S)...)
function getCentralPlaquette(S)
    allplaqs = collect(plaquetteIterator(S))
    central = centralPos(S)
    return allplaqs[findmin([norm(central .- r) for r in allplaqs])[2]]
end
function symmetryReducePlaquettes(S,R_ref)
    
    AllDists = Dict{SVector{2,Int},Int}()

    Allplaqs = collect(plaquetteIterator(S))

    indicesMapping = Int[]
    uniqueInds = Int[]
    Lx,Ly = size(S)
    ri = R_ref
    for (j,rj) in enumerate(Allplaqs)
        Rij = getReducedDist(ri,rj,Lx,Ly)
        if Rij ∉ keys(AllDists)
            uniqueInds = push!(uniqueInds,j)
            AllDists[Rij] = length(uniqueInds)
        end
        push!(indicesMapping,AllDists[Rij])
    end
    return (;AllDists,indicesMapping,uniqueInds)
end

function getDistReduction(S,ψG::LocalPlaquetteGuidingFunction)
    AllDists = Dict{SVector{2,Rational{Int}},Int}()
    if isperiodic(S)
        indicesMapping = ones(Int,length(ψG.params))
        uniqueInds = [1]
        return AllDists,indicesMapping,uniqueInds
    end
    
    α = get_alpha_i(ψG)
    Allplaqs = collect(plaquetteIterator(S))

    indicesMapping = Int[]
    uniqueInds = Int[]
    LxLy = size(S)
    r_Central = centralPos(S)
    for (i,ri) in enumerate(Allplaqs)
        x,y = ri .- r_Central
        if y < -x
            x,y = -y,-x
        end
        if y>x
            x,y = y,x
        end

        symMapped = SVector(x,y)
        if symMapped ∉ keys(AllDists)
            uniqueInds = push!(uniqueInds,i)
            AllDists[symMapped] = length(uniqueInds)
        end
        push!(indicesMapping,AllDists[symMapped])
    end
    
    return (;AllDists,indicesMapping,uniqueInds)

end

function add_reconstructedFullParams!(ψG,indicesMapping,trimmedparams)
    for (i,k) in enumerate(indicesMapping)
        ψG.params[i] += trimmedparams[k]
    end
    return ψG
end

function reconf_obs(InitialState::ConfType,method::AbstractGFMCMethod,configs,ψG,inequivParams=eachindex(ψG.params)) where {ConfType}
    plaqs = collect(plaquetteIterator(InitialState))
    AffectedPlaquetteList = precomputeAffectedPlaquettes(InitialState)
    
    NThreads = Threads.nthreads()

    Nparams = length(inequivParams)
    Ok_i = zeros(Float32,length(configs),Nparams)
    E_i = zeros(Float32,length(configs))

    WorkChunks = ChunkSplitters.chunks(eachindex(configs),n=NThreads)
    
    # println("collecting obs")
    # @time Threads.@threads for (ichunk,chunkinds) in enumerate(WorkChunks)
    Threads.@threads for (ichunk,chunkinds) in enumerate(WorkChunks)

        Walker = spiderWebWalker(InitialState,plaqs)

        for (i,iconf) in enumerate(chunkinds)
            get_config(Walker) .= configs[iconf]
            updateWeightList!(Walker,AffectedPlaquetteList,ψG)
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
struct IterativeSRSolver <: AbstractSRSolver end

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
    res = -IterativeSolvers.cg(S,F;kwargs...)
    return res
end
stochastic_reconfiguration_step(E_i,Ok_i,::AbstractSRSolver) = error("solver not implemented")

function _stochastic_reconfiguration(InitialState,method::AbstractGFMCMethod,solver::AbstractSRSolver,NSteps::AbstractVector,ψG,n,dt::AbstractVector,equilibration_steps=1000,rel_tolerance=1e-2,Nwalkers = 6,outfile=nothing,reconfigure=true,initializer = UnguidedWalkInitializer(equilibration_steps,0.8);verbose=true,report_steps=1,reset = true)
    
    ψG = deepcopy(ψG)
    params = ψG.params

    convergedSteps = 0

    normDelta = Inf

    _,indicesMapping,uniqueInds = getDistReduction(InitialState,ψG)
    maxNSteps = maximum(NSteps)
    prob = setup_GFMC_problem(InitialState,method,Nwalkers,maxNSteps,ψG) 
    initializeGFMC!(prob,equilibration_steps,initializer)

    results = get_stoch_rec_Observables(n,ψG,outfile)

    ind = 0
    for i in 1:n
        
        range = eachindex(prob.Observables.TotalWeights)[1:NSteps[i]]
        # for w in prob.Walkers
        #     get_config(w) .= InitialState
        # end
        reset && initializeGFMC!(prob,equilibration_steps,initializer)

        # @time res = runGFMC!(prob,range,reconfigure)
        res = runGFMC!(prob,range,reconfigure)

        resSlice = @view res.SaveConfigs[:,:,:,range]
        confs = eachslice( resSlice,dims=(3,4))
        
        observables = reconf_obs(InitialState,method,confs,ψG,uniqueInds)
        # (;EL_avg,EL_err ) = observables
        # return observables
        δα = stochastic_reconfiguration_step(observables.E_i,observables.Ok_i,solver)
        normDelta = norm(δα)/norm(params)

        add_reconstructedFullParams!(ψG,indicesMapping,δα .*dt[i])
        
        E0 = mean(res.energies[range])
        results.E0[i] = E0
        ΔE0 = sqrt(var(res.energies[range]))
        results.ΔE[i] = ΔE0
        
        # push!(params_steps, copy(params))
        selectdim(results.params_steps,arraydim(results.params_steps),i) .= params

        α = get_alpha_i(ψG)
        w_avg = mean(res.TotalWeights[range])
        if verbose && i % report_steps == 0
            @info "optimization step $i" dt[i] NSteps[i] "|δα|" = normDelta δα[1] E0 ΔE0 = ΔE0 convergedSteps mean(α) w_avg
        end
        ind = i
        if normDelta < rel_tolerance
            convergedSteps += 1

            convergedSteps > 4 && break
        else
            convergedSteps = 0
        end
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
    reconfigure=true,
    initializer = UnguidedWalkInitializer(pre_equilibration_steps,scatter_fraction),
    kwargs...)
    
    NStepsVec = makeVec(NSteps,n)
    dtVec = makeVec(dt,n)
    return _stochastic_reconfiguration(InitialState,method,solver,NStepsVec,ψG,n,dtVec,equilibration_steps,rel_tolerance,Nwalkers,outfile,reconfigure,initializer;kwargs...)
end

makeVec(x::AbstractVector,len) = x
makeVec(x::Number,len) = fill(x,len)
makeVec(f::Function,len) = f.(1:len)

function get_stoch_rec_Observables(Nsteps,ψG,::Nothing)
    E0 = zeros(Nsteps)
    ΔE = zeros(Nsteps)
    params = zeros(size(ψG.params)...,Nsteps)
    return (;E0 = E0,ΔE = ΔE,params_steps = params)
end

function get_stoch_rec_Observables(Nsteps,ψG,outfile::AbstractString)
    h5open(outfile,"cw") do file
        E0 = createMMapArray(file,"E0",Float64,(Nsteps,))
        ΔE = createMMapArray(file,"ΔE",Float64,(Nsteps,))
        params = createMMapArray(file,"params_steps",Float32,(size(ψG.params)...,Nsteps))
        E0 .= 0
        ΔE .= 0
        params .= 0
        return (;E0 = E0,ΔE = ΔE,params_steps = params)
    end
end