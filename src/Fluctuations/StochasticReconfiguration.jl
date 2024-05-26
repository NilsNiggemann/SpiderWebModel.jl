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

function reconf_obs(InitialState::ConfType,configs,ψG,Λ) where {ConfType}
    plaqs = collect(plaquetteIterator(InitialState))
    AffectedPlaquetteList = precomputeAffectedPlaquettes(InitialState)
    
    NThreads = Threads.nthreads()

    Nparams = length(ψG.params)
    Oik = zeros(length(configs),Nparams)
    E_i = zeros(length(configs))

    WorkChunks = ChunkSplitters.chunks(eachindex(configs),n=NThreads)
    
    Threads.@threads for (ichunk,chunkinds) in enumerate(WorkChunks)

        Walker = spiderWebWalker(InitialState,plaqs)

        @inbounds for (i,iconf) in enumerate(chunkinds)
            get_config(Walker) .= configs[iconf]
            moves = getMoves!(Walker)
            weights = getWeightList!(Walker,AffectedPlaquetteList,ψG,Λ)
            elocal = getLocalEnergy(weights,Λ)
            for k in 1:Nparams
                O_xk = getOx_k(ψG,Walker.n_x,k)
                Oik[iconf,k] = O_xk
            end
            E_i[iconf] = elocal
        end

    end
    N = length(configs)

    EL_avg = mean(E_i)
    EL_err = sqrt(var(E_i))

    S = cov(Oik,dims=1) #+ 1e-12I
    for i in axes(S,1)
        S[i,i] += 1e-12
    end
    F = cov(Oik,E_i)
    
    return (;EL_avg,EL_err,S,F)
end

function stochastic_reconfiguration_step(InitialState,configs,ψG,Λ)
    (;EL_avg,EL_err,F,S) = reconf_obs(InitialState,configs,ψG,Λ)
    δα = - S \ F
    return (;δα,EL_avg,EL_err)
end

function stochastic_reconfiguration(InitialState,nSteps,ψG,n,dt=1e-3;equilibration_steps=1000,Λ=1,rel_tolerance=1e-2,Nwalkers = 6,nbra =10)
    
    ψG = deepcopy(ψG)
    params = ψG.params
    
    convergedSteps = 0

    E0s = fill(NaN,n)
    ΔE = fill(NaN,n)
    # αs = fill(NaN,n)
    normDelta = Inf
    for i in 1:n
        
        @time res = startManyWalkerGFMC(InitialState,Nwalkers,nSteps,nbra,ψG,Λ;equilibration_steps)

        # αs[i] = params

        # conf_dims = size(res.SaveConfigs)
        # confs = reshape(res.SaveConfigs,conf_dims[1],conf_dims[2],conf_dims[3]*conf_dims[4])
        confs = eachslice(res.SaveConfigs,dims=(3,4))

        @time (;δα,EL_avg,EL_err) = stochastic_reconfiguration_step(InitialState,confs,ψG,Λ)
        
        normDelta = norm(δα)/norm(params)

        params[:] .+= dt*δα

        E0s[i] = EL_avg
        ΔE[i] = EL_err
        # α = get_alpha_i(ψG) 
        # @views α .+= dt*δα[begin:length(α)]
        # @info "optimization step $i" δα params E0 = mean(res.energies) σ = sqrt(var(res.energies)) convergedSteps
        @info "optimization step $i" "δα" = normDelta E0 = mean(res.energies) ΔE0 = sqrt(var(res.energies)) convergedSteps

        # ψG = VariationalGuidingFunction(params) 
        if normDelta < rel_tolerance
            convergedSteps += 1

            convergedSteps > 4 && break
        else
            convergedSteps = 0
        end
    end
    # return (;params = params, E0_i = E0s,ΔE_i = σs,α_i = αs)
    return (;params = params,E0_i = E0s,ΔE_i = ΔE)
end