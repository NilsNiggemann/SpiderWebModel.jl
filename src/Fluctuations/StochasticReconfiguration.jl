function reconf_obs(InitialState::ConfType,configs,ψG,Λ) where {ConfType}
    plaqs = collect(plaquetteIterator(InitialState))
    AffectedPlaquetteList = precomputeAffectedPlaquettes(InitialState)
    
    NThreads = Threads.nthreads()
    # NThreads = 1
    O_sum_buffer = zeros(NThreads)
    O²_sum_buffer = zeros(NThreads)
    EL_sum_buffer = zeros(NThreads)
    OEL_sum_buffer = zeros(NThreads)

    Walkers = Vector{SpiderWebWalker{ConfType}}(undef,NThreads)
    Threads.@threads for α in eachindex(Walkers)
        Walkers[α] = spiderWebWalker(InitialState,plaqs)
    end

    WorkChunks = ChunkSplitters.chunks(eachindex(configs),n=NThreads)

    Threads.@threads for (ichunk,chunkinds) in enumerate(WorkChunks)
        Walker = Walkers[ichunk]
        O_sum = 0.
        O²_sum = 0.
        EL_sum = 0.
        OEL_sum = 0.
        for i in chunkinds
            get_config(Walker) .= configs[i]
            moves = getMoves!(Walker)
            weights = getWeightList!(Walker,AffectedPlaquetteList,ψG,Λ)
            elocal = getLocalEnergy(weights,Λ)
            O_x = ψG.α * length(moves)
            O_sum += O_x
            O²_sum += O_x^2
            EL_sum += elocal
            OEL_sum += O_x*elocal
        end
        O_sum_buffer[ichunk] = O_sum
        O²_sum_buffer[ichunk] = O²_sum
        EL_sum_buffer[ichunk] = EL_sum
        OEL_sum_buffer[ichunk] = OEL_sum
    end
    N = length(configs)
    O_avg = sum(O_sum_buffer)/N
    O²_avg = sum(O²_sum_buffer)/N
    EL_avg = sum(EL_sum_buffer)/N
    OEL_avg = sum(OEL_sum_buffer)/N
    return (;O_avg,O²_avg,EL_avg,OEL_avg)
end

function stochastic_reconfiguration(InitialState,configs,ψG,Λ)
    (;O_avg,O²_avg,EL_avg,OEL_avg) = reconf_obs(InitialState,configs,ψG,Λ)
    F = OEL_avg - O_avg*EL_avg
    InitialState = O²_avg - O_avg^2
    δα = F/InitialState
end

function repeatStochReconf(InitialState,nSteps,ψG,n,dt=1e-3;equilibration_steps=1000,Λ=1,error_threshold=1e-2,Nwalkers = 6,nbra =10)
    α = ψG.α
    convergedSteps = 0

    E0s = fill(NaN,n)
    σs = fill(NaN,n)
    αs = fill(NaN,n)

    for i in 1:n
        
        res = startManyWalkerGFMC(InitialState,Nwalkers,nSteps,nbra,ψG,Λ;equilibration_steps)
        E0s[i] = mean(res.energies)
        σs[i] = sqrt(var(res.energies))
        αs[i] = α

        # conf_dims = size(res.SaveConfigs)
        # confs = reshape(res.SaveConfigs,conf_dims[1],conf_dims[2],conf_dims[3]*conf_dims[4])
        confs = eachslice(res.SaveConfigs,dims=(3,4))
        δα = -stochastic_reconfiguration(InitialState,confs,ψG,Λ)
        α += dt*δα
        @info "optimization step $i" δα α E0 = mean(res.energies) σ = sqrt(var(res.energies)) convergedSteps
        ψG = PlaquetteNumberGuidingFunction(α) 
        if abs(δα) < error_threshold
            convergedSteps += 1

            convergedSteps > 4 && break
        else
            convergedSteps = 0
        end
    end
    return (;α = α, E0_i = E0s,ΔE_i = σs,α_i = αs)
end