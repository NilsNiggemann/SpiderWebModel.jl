function getRandConfs(InitialState::ConfType, N, Nwalkers; equilibration_steps = 1000,samplingRate=1e-6) where ConfType
    plaquettePositions = collect(plaquetteIterator(InitialState))

    Walkers = Vector{SpiderWebWalker{ConfType}}(undef,Nwalkers)
    Threads.@threads for α in eachindex(Walkers)
        Walkers[α] = spiderWebWalker(InitialState,plaquettePositions)
    end
    Lx,Ly = size(InitialState)
    SaveConfigs = zeros(eltype(InitialState),Lx,Ly,N,Nwalkers)

    random_init_walkers!(Walkers,equilibration_steps)
    Threads.@threads for α in eachindex(Walkers)
        Walker = Walkers[α]
        
        n = 1
        randPlaqs = rand(plaquettePositions,10N)
        randmoves = rand(Bool,10N)
        iter = 1
        while n<=N
            # movepos = Tuple(rand(plaquettePositions))
            iter +=1
            if iter > lastindex(randPlaqs)
                iter = 1
                rand!(plaquettePositions,randPlaqs)
                rand!(randmoves)
            end
                
            movepos = randPlaqs[iter]
            movesgn = rand(1:2)
            # movesgn = randmoves[iter] + 1
            P_applicable(Walker.Config, movepos)[movesgn] || continue
            applyPlaquette!(Walker.Config, movepos[1], movepos[2], (1,-1)[movesgn])

            if rand() < samplingRate
                SaveConfigs[:,:,n,α] .= get_config(Walkers[α])
                n +=1
            end
        end
    end
    return SaveConfigs
end