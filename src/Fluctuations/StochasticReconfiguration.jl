function reconf_obs(S,configs,ψG,Λ)
    plaqs = collect(plaquetteIterator(S))
    AffectedPlaquetteList = precomputeAffectedPlaquettes(S)
    
    Walker = SpiderWebWalker(S,plaqs)

    O_avg = 0.
    O²_avg = 0.
    EL_avg = 0.
    OEL_avg = 0.
    
    for conf in configs
        get_config(Walker) .= conf
        moves = getMoves!(Walker)
        weights = getWeightList!(Walker,AffectedPlaquetteList,ψG,Λ)
        elocal = getLocalEnergy(weights,Λ)
        O_x = ψG.α * length(moves)
        O_avg += O_x
        O²_avg += O_x^2
        EL_avg += elocal
        OEL_avg += O_x*elocal
    end
    O_avg /= length(configs)
    O²_avg /= length(configs)
    EL_avg /= length(configs)
    OEL_avg /= length(configs)
    return (;O_avg,O²_avg,EL_avg,OEL_avg)
end

function stochastic_reconfiguration(S,configs,ψG,Λ)
    (;O_avg,O²_avg,EL_avg,OEL_avg) = reconf_obs(S,configs,ψG,Λ)
    F = OEL_avg - O_avg*EL_avg
    S = O²_avg - O_avg^2
    δα = F/S
end

function repeatStochReconf(S,nSteps,ψG,n,dt=1e-3;equilibration_steps=1000,Λ=1,error_threshold=1e-2,Nwalkers = 6,nbra =10)
    α = ψG.α
    convergedSteps = 0

    E0s = fill(NaN,n)
    σs = fill(NaN,n)
    αs = fill(NaN,n)

    for i in 1:n
        
        res = startManyWalkerGFMC(S,Nwalkers,nSteps,nbra,ψG,Λ;equilibration_steps)
        E0s[i] = mean(res.energies)
        σs[i] = sqrt(var(res.energies))
        αs[i] = α

        confs = eachslice(res.SaveConfigs,dims=(3,4))
        δα = -stochastic_reconfiguration(S,confs,ψG,Λ)
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