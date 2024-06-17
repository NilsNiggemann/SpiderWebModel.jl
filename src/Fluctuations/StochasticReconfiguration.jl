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
    
    y .= zk1# .- zk2

end

function getOx_k(ψG::FullVariationalGuidingFunction,n::AbstractArray,k)
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

getOx_k(::LocalPlaquetteGuidingFunction,n::AbstractArray,k) = n[k]

function nearbyInt(x1,x2,x_size)
    x_rsize = 1.0 / x_size

    dx = x1 - x2
    dx -= x_size * round(Int,dx * x_rsize)
end

isperiodic(S::Stencils.StencilArray) = Stencils.boundary(S) == Stencils.Wrap()
isperiodic(S::StencilSpinConfig) = isperiodic(parent(S))

function getDistReduction(S,ψG::FullVariationalGuidingFunction)
    
    AllDists = Dict{SVector{2,Int},Int}()
    if !isperiodic(S)
        indicesMapping = collect(eachindex(ψG.params))
        uniqueInds = collect(indicesMapping)
        return AllDists,indicesMapping,uniqueInds
    end

    α = get_alpha_i(ψG)
    Allplaqs = collect(plaquetteIterator(S))

    # AllDists = Dict{Tuple{Int,Int},Int}()
    
    betaIndex = lastindex(α)
    indicesMapping = ones(Int,betaIndex)
    uniqueInds = [1]
    # indicesMapping = collect(eachindex(α))
    # uniqueInds = collect(eachindex(α))
    LxLy = size(S)

    for (i,ri) in enumerate(Allplaqs)
        for (j,rj) in enumerate(Allplaqs)
            Rij = abs.(SVector(nearbyInt.(ri, rj,LxLy)))
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

    return AllDists,indicesMapping,uniqueInds

end

function getDistReduction(S,ψG::LocalPlaquetteGuidingFunction)
    AllDists = Dict{SVector{2,Rational{Int}},Int}()
    # if !isperiodic(S)
    #     indicesMapping = collect(eachindex(ψG.params))
    #     uniqueInds = collect(indicesMapping)
    #     return AllDists,indicesMapping,uniqueInds
    # end
    
    α = get_alpha_i(ψG)
    Allplaqs = collect(plaquetteIterator(S))

    indicesMapping = Int[]
    uniqueInds = Int[]
    LxLy = size(S)
    r_Central = (LxLy .+1) .//2 
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
    
    return AllDists,indicesMapping,uniqueInds

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
    
    println("collecting obs")
    @time Threads.@threads for (ichunk,chunkinds) in enumerate(WorkChunks)

        Walker = spiderWebWalker(InitialState,plaqs)

        for (i,iconf) in enumerate(chunkinds)
            get_config(Walker) .= configs[iconf]
            updateWeightList!(Walker,AffectedPlaquetteList,ψG)
            elocal = getLocalEnergy(Walker,method)
            @inbounds for (ik,k) in enumerate(inequivParams)
                O_xk = getOx_k(ψG,Walker.n_x,k)
                Ok_i[iconf,ik] = O_xk
            end
            E_i[iconf] = elocal
        end

    end
    N = length(configs)

    EL_avg = mean(E_i)
    EL_err = sqrt(var(E_i))

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
        !issuccess(SChol) && return zeros(length(ψG.params))
    end
    @time δα = SChol \ F
    for i in axes(δα,1)
        δα[i] = -δα[i]
    end
    return δα
end
function stochastic_reconfiguration_step(E_i::AbstractVector,Ok_i::AbstractMatrix,solver::IterativeSRSolver;kwargs...)
    S = QuantumMetric(Ok_i)
    F = StatsBase.cov(Ok_i,E_i)[:]
    res = -IterativeSolvers.cg(S,F;kwargs...)
    return res
end
stochastic_reconfiguration_step(E_i,Ok_i,::AbstractSRSolver) = error("solver not implemented")

function _stochastic_reconfiguration(InitialState,method::AbstractGFMCMethod,solver::AbstractSRSolver,NSteps::AbstractVector,ψG,n,dt::AbstractVector,equilibration_steps=1000,rel_tolerance=1e-2,Nwalkers = 6,outfile=nothing,pre_equilibration_steps=5*equilibration_steps)
    
    ψG = deepcopy(ψG)
    params = ψG.params
    
    convergedSteps = 0

    E0s = Vector{Float64}()
    ΔE = Vector{Float64}()
    params_steps = Vector{typeof(params)}()
    normDelta = Inf

    _,indicesMapping,uniqueInds = getDistReduction(InitialState,ψG)
    maxNSteps = maximum(NSteps)
    prob = setup_GFMC_problem(InitialState,method,Nwalkers,maxNSteps,ψG) 
    initializeGFMC!(prob,equilibration_steps,pre_equilibration_steps)
    for i in 1:n
        
        range = eachindex(prob.TotalWeights)[1:NSteps[i]]
        @time res = runGFMC!(prob,range)

        resSlice = @view res.SaveConfigs[:,:,:,range]
        confs = eachslice( resSlice,dims=(3,4))
        
        observables = reconf_obs(InitialState,method,confs,ψG,uniqueInds)
        # (;EL_avg,EL_err ) = observables
        # return observables
        δα = stochastic_reconfiguration_step(observables.E_i,observables.Ok_i,solver)
        normDelta = norm(δα)/norm(params)

        add_reconstructedFullParams!(ψG,indicesMapping,δα .*dt[i])
        
        E0 = mean(res.energies[range])
        push!(E0s, E0)
        ΔE0i = sqrt(var(res.energies[range]))
        push!(ΔE,ΔE0i)

        push!(params_steps, copy(params))

        α = get_alpha_i(ψG)
        @info "optimization step $i" dt[i] NSteps[i] "|δα|" = normDelta δα[1] E0 ΔE0 = ΔE0i convergedSteps mean(α) 

        if normDelta < rel_tolerance
            convergedSteps += 1

            convergedSteps > 4 && break
        else
            convergedSteps = 0
        end
    end
    params_steps_arr = stack(params_steps)
    if !isnothing(outfile)
        h5write(outfile,"params",ψG.params)
        h5write(outfile,"E0s",E0s)
        h5write(outfile,"ΔE",ΔE)
        h5write(outfile,"params_steps",params_steps_arr)
    end

    return (;params = params,E0_i = E0s,ΔE_i = ΔE,params_steps = params_steps_arr)
end

function stochastic_reconfiguration(InitialState,method::AbstractGFMCMethod,NSteps,ψG,n,dt,solver::AbstractSRSolver = IterativeSRSolver();equilibration_steps=1000,rel_tolerance=1e-2,Nwalkers = 6,outfile=nothing,pre_equilibration_steps=5*equilibration_steps)
    
    NStepsVec = makeVec(NSteps,n)
    dtVec = makeVec(dt,n)

    return _stochastic_reconfiguration(InitialState,method,solver,NStepsVec,ψG,n,dtVec,equilibration_steps,rel_tolerance,Nwalkers,outfile,pre_equilibration_steps)
end

makeVec(x::AbstractVector,len) = x
makeVec(x::Number,len) = fill(x,len)
makeVec(f::Function,len) = f.(1:len)