function getOx_k(ψG::VariationalGuidingFunction,n::AbstractArray,k)
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

function getDistReduction(S,ψG)

    α = get_alpha_i(ψG)
    Allplaqs = collect(plaquetteIterator(S))
    # AllDists = Dict{SVector{2,Int},Int}()
    AllDists = Dict{Tuple{Int,Int},Int}()
    
    betaIndex = lastindex(α)
    indicesMapping = ones(Int,betaIndex)
    uniqueInds = [1]
    return AllDists,indicesMapping,uniqueInds
    for (i,ri) in enumerate(Allplaqs)
        for (j,rj) in enumerate(Allplaqs)
            # Rij = SVector(ri .- rj) #.% (size(S) .÷2)
            Rij = (i,j) #.% (size(S) .÷2)
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

function add_reconstructedFullParams!(ψG,indicesMapping,trimmedparams)
    for (i,k) in enumerate(indicesMapping)
        ψG.params[i] += trimmedparams[k]
    end
    return ψG
end
import CovarianceEstimation

function reconf_obs(InitialState::ConfType,configs,ψG,Λ,inequivParams=eachindex(ψG.params)) where {ConfType}
    plaqs = collect(plaquetteIterator(InitialState))
    AffectedPlaquetteList = precomputeAffectedPlaquettes(InitialState)
    
    NThreads = Threads.nthreads()

    Nparams = length(inequivParams)
    Ok_i = zeros(length(configs),Nparams)
    E_i = zeros(length(configs))

    WorkChunks = ChunkSplitters.chunks(eachindex(configs),n=NThreads)
    
    println("collecting obs")
    @time Threads.@threads for (ichunk,chunkinds) in enumerate(WorkChunks)

        Walker = spiderWebWalker(InitialState,plaqs)

        for (i,iconf) in enumerate(chunkinds)
            get_config(Walker) .= configs[iconf]
            moves = getMoves!(Walker)
            weights = getWeightList!(Walker,AffectedPlaquetteList,ψG,Λ)
            elocal = getLocalEnergy(weights,Λ)
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

    println("computing cov")
    # @time S = cov(Ok_i) #+ 1e-12I
    method = CovarianceEstimation.LinearShrinkage(CovarianceEstimation.ConstantCorrelation())
    @time S = cov(method,(Ok_i)) #+ 1e-12I
    for i in axes(S,1)
        S[i,i] += 1e-5
    end
    @time F = StatsBase.cov(Ok_i,E_i)
    
    return (;EL_avg,EL_err,S,F,Ok_i)
end

function stochastic_reconfiguration_step(InitialState,configs,ψG,Λ,inequivalentIndices=eachindex(ψG.params))
    (;EL_avg,EL_err,F,S) = reconf_obs(InitialState,configs,ψG,Λ,inequivalentIndices)
    println("constructing cholesky")
    @time SChol = cholesky!(Symmetric(S),check = false)
    if !issuccess(SChol)
        @warn "Cholesky factorization failed"
        SChol = cholesky!(Symmetric(S+1e-3I),check = false)
        !issuccess(SChol) && return (;δα = zeros(length(ψG.params)),EL_avg,EL_err)
    end
    @time δα = SChol \ F
    for i in axes(δα,1)
        δα[i] = -δα[i]
    end
    return (;δα,EL_avg,EL_err)
end


function stochastic_reconfiguration(InitialState,nSteps,ψG,n,dt=1e-3;equilibration_steps=1000,Λ=1,rel_tolerance=1e-2,Nwalkers = 6,nbra =10)
    
    ψG = deepcopy(ψG)
    params = ψG.params
    
    convergedSteps = 0

    E0s = fill(NaN,n)
    ΔE = fill(NaN,n)
    AllParams = zeros(size(params)...,n)
    # αs = fill(NaN,n)
    normDelta = Inf

    AllDists,indicesMapping,uniqueInds = getDistReduction(InitialState,ψG)
    # inequivalentIndices = unique(indicesMapping)

    for i in 1:n
        
        @time res = startManyWalkerGFMC(InitialState,Nwalkers,nSteps,nbra,ψG,Λ;equilibration_steps,pre_equilibration_steps=5equilibration_steps)

        # αs[i] = params

        # conf_dims = size(res.SaveConfigs)
        # confs = reshape(res.SaveConfigs,conf_dims[1],conf_dims[2],conf_dims[3]*conf_dims[4])
        confs = eachslice(res.SaveConfigs,dims=(3,4))

        @time (;δα,EL_avg,EL_err) = stochastic_reconfiguration_step(InitialState,confs,ψG,Λ,uniqueInds)
        
        normDelta = norm(δα)/norm(params)

        # params[eachindex(δα)] .+= dt*δα
        add_reconstructedFullParams!(ψG,indicesMapping,δα .*dt)
        E0s[i] = EL_avg
        ΔE[i] = EL_err
        AllParams[:,:,i] .= ψG.params
        α = get_alpha_i(ψG) 
        # @views α .+= dt*δα[begin:length(α)]
        # @info "optimization step $i" δα params E0 = mean(res.energies) σ = sqrt(var(res.energies)) convergedSteps
        @info "optimization step $i" "δα" = normDelta E0 = mean(res.energies) ΔE0 = sqrt(var(res.energies)) convergedSteps mean(α)

        # ψG = VariationalGuidingFunction(params) 
        if normDelta < rel_tolerance
            convergedSteps += 1

            convergedSteps > 4 && break
        else
            convergedSteps = 0
        end
    end
    # return (;params = params, E0_i = E0s,ΔE_i = σs,α_i = αs)
    return (;params = params,E0_i = E0s,ΔE_i = ΔE,AllParams = AllParams)
end