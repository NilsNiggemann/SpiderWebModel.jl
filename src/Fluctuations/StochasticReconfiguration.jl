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

function getDistReduction(S,ψG::FullVariationalGuidingFunction)

    α = get_alpha_i(ψG)
    Allplaqs = collect(plaquetteIterator(S))
    AllDists = Dict{SVector{2,Int},Int}()
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

    α = get_alpha_i(ψG)
    Allplaqs = collect(plaquetteIterator(S))
    AllDists = Dict{SVector{2,Rational{Int}},Int}()

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

function reconf_obs(InitialState::ConfType,configs,ψG,Λ,inequivParams=eachindex(ψG.params)) where {ConfType}
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
    # target = CovarianceEstimation.PerfectPositiveCorrelation()
    # shrinkage = :lw # Ledoit-Wolf optimal shrinkage
    # method = CovarianceEstimation.LinearShrinkage(target,0.)


    @time S = cov(Ok_i) #+ 1e-12I
    GC.gc()
    for i in axes(S,1)
        S[i,i] += 1e-5
    end
    @time F = StatsBase.cov(Ok_i,E_i)
    
    return (;EL_avg,EL_err,S,F,Ok_i)
end

function stochastic_reconfiguration_step(InitialState,configs,ψG,Λ,inequivalentIndices=eachindex(ψG.params))
    (;EL_avg,EL_err,F,S) = reconf_obs(InitialState,configs,ψG,Λ,inequivalentIndices)
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


function stochastic_reconfiguration(InitialState,nSteps,ψG,n,dt=1e-3;equilibration_steps=1000,Λ=1,rel_tolerance=1e-2,Nwalkers = 6,nbra =10,outfile=nothing)
    
    ψG = deepcopy(ψG)
    params = ψG.params
    
    convergedSteps = 0

    E0s = fill(NaN,n)
    ΔE = fill(NaN,n)
    params_steps = [similar(params) for _ in 1:n]
    normDelta = Inf

    _,indicesMapping,uniqueInds = getDistReduction(InitialState,ψG)

    for (i,ParamsSlice) in enumerate(params_steps)
        
        @time res = startManyWalkerGFMC(InitialState,Nwalkers,nSteps,nbra,ψG,Λ;equilibration_steps,pre_equilibration_steps=5equilibration_steps)

        confs = eachslice(res.SaveConfigs,dims=(3,4))

        (;δα,EL_avg,EL_err) = stochastic_reconfiguration_step(InitialState,confs,ψG,Λ,uniqueInds)
        
        normDelta = norm(δα)/norm(params)

        add_reconstructedFullParams!(ψG,indicesMapping,δα .*dt)
        E0s[i] = mean(res.energies)
        ΔE[i] = sqrt(var(res.energies))
        ParamsSlice .= ψG.params
        α = get_alpha_i(ψG) 
        @info "optimization step $i" "|δα|" = normDelta E0 = mean(res.energies) ΔE0 = sqrt(var(res.energies)) convergedSteps mean(α) δα[1]

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