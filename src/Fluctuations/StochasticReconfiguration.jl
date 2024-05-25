
function singleElementmatrix(size,idx,val)
    arr = zeros(size)
    arr[idx] = val
    return sparse(arr)
end

function reconf_obs(InitialState::ConfType,configs,ψG,Λ) where {ConfType}
    plaqs = collect(plaquetteIterator(InitialState))
    AffectedPlaquetteList = precomputeAffectedPlaquettes(InitialState)
    
    NThreads = Threads.nthreads()
    # NThreads = 1

    Nparams = length(ψG.params)
    
    O_sum_buffer = [zeros(Nparams) for _ in 1:NThreads]
    O²_sum_buffer = [zeros(Nparams,Nparams) for _ in 1:NThreads]
    EL_sum_buffer = zeros(NThreads)
    OEL_sum_buffer = [zeros(Nparams) for _ in 1:NThreads]

    Walkers = Vector{SpiderWebWalker{ConfType}}(undef,NThreads)
    Threads.@threads for α in eachindex(Walkers)
        Walkers[α] = spiderWebWalker(InitialState,plaqs)
    end

    PseudoParametersMatrices = [
        singleElementmatrix(size(ψG.params),i,1.) for i in 1:Nparams
    ]

    WorkChunks = ChunkSplitters.chunks(eachindex(configs),n=NThreads)
    
    Threads.@threads for (ichunk,chunkinds) in enumerate(WorkChunks)

        Walker = Walkers[ichunk]
        O_sum = O_sum_buffer[ichunk]
        O²_sum = O²_sum_buffer[ichunk]
        OEL_sum = OEL_sum_buffer[ichunk]

        O_x = zeros(length(chunkinds),Nparams)

        @inbounds for (i,iconf) in enumerate(chunkinds)
            get_config(Walker) .= configs[iconf]
            moves = getMoves!(Walker)
            weights = getWeightList!(Walker,AffectedPlaquetteList,ψG,Λ)
            elocal = getLocalEnergy(weights,Λ)
            for k in 1:Nparams
                PseudoGuidingWF = VariationalGuidingFunction(PseudoParametersMatrices[k])
                O_xk = guidingfunc_exponent(PseudoGuidingWF,Walker.n_x)
                O_x[i,k] = O_xk
                OEL_sum[k] += O_xk*elocal
                O_sum[k] += O_xk
            end
            EL_sum_buffer[ichunk] += elocal
        end

        @inbounds for k in 1:Nparams
            for j in 1:Nparams
                O²_sumval = 0.
                for (i,iconf) in enumerate(chunkinds)
                    O²_sumval += O_x[i,k] * O_x[i,j]
                end
                O²_sum[j,k] = O²_sumval
            end
        end
        # O²_sum .+= O²_sum'
        # O²cov_sum .= naiveCov(O_x)
        # @inbounds for k in 1:Nparams
        #     for j in k+1:Nparams
        #         O²cov_sum[k,j] = O²cov_sum[j,k]
        #     end
        # end
    end
    N = length(configs)
    O_avg = sum(O_sum_buffer)./N
    O²_avg = sum(O²_sum_buffer)./N
    EL_avg = sum(EL_sum_buffer)./N
    OEL_avg = sum(OEL_sum_buffer)./N
    return (;O_avg,O²_avg,EL_avg,OEL_avg)
end

function stochastic_reconfiguration(InitialState,configs,ψG,Λ)
    (;O_avg,O²_avg,EL_avg,OEL_avg) = reconf_obs(InitialState,configs,ψG,Λ)
    F = OEL_avg .- O_avg .* EL_avg
    S = O²_avg - O_avg * O_avg' + 1e-5I
    # S = O²_avg + 1e-5I
    return δα = S \ F
end

function repeatStochReconf(InitialState,nSteps,ψG,n,dt=1e-3;equilibration_steps=1000,Λ=1,error_threshold=1e-2,Nwalkers = 6,nbra =10)
    
    ψG = deepcopy(ψG)
    params = ψG.params
    
    convergedSteps = 0

    E0s = fill(NaN,n)
    ΔE = fill(NaN,n)
    # αs = fill(NaN,n)
    normDelta = Inf
    for i in 1:n
        
        @time res = startManyWalkerGFMC(InitialState,Nwalkers,nSteps,nbra,ψG,Λ;equilibration_steps)
        E0s[i] = mean(res.energies)
        ΔE[i] = sqrt(var(res.energies))
        # αs[i] = params

        # conf_dims = size(res.SaveConfigs)
        # confs = reshape(res.SaveConfigs,conf_dims[1],conf_dims[2],conf_dims[3]*conf_dims[4])
        confs = eachslice(res.SaveConfigs,dims=(3,4))
        @time for _ in 1:1
            δα = -stochastic_reconfiguration(InitialState,confs,ψG,Λ)
            normDelta = norm(δα)
            params[:] .+= dt*δα
        end
        # α = get_alpha_i(ψG) 
        # @views α .+= dt*δα[begin:length(α)]
        # @info "optimization step $i" δα params E0 = mean(res.energies) σ = sqrt(var(res.energies)) convergedSteps
        @info "optimization step $i" "δα" = normDelta E0 = mean(res.energies) ΔE0 = sqrt(var(res.energies)) convergedSteps

        # ψG = VariationalGuidingFunction(params) 
        if normDelta < error_threshold
            convergedSteps += 1

            convergedSteps > 4 && break
        else
            convergedSteps = 0
        end
    end
    # return (;params = params, E0_i = E0s,ΔE_i = σs,α_i = αs)
    return (;params = params,E0_i = E0s,ΔE_i = ΔE)
end