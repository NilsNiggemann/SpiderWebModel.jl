function plaquettesAreSeparated(P1, P2)
    separation = P2 .- P1

    return any(x -> abs(x) > 2, separation)
end

function getApplicablePlaquettes!(plaqsPos, Conf::SpinConfig)
    empty!(plaqsPos)

    for i in axes(Conf.Mat, 1)
        for j in axes(Conf.Mat, 2)
            if canFlipPlaquette(Conf, i, j)
                push!(plaqsPos, (i, j))
            end
        end
    end
    return plaqsPos
end

function getRandomSeparatedPlaquettes!(plaquettes, sep_plaquettes, Conf)
    getApplicablePlaquettes!(plaquettes, Conf)
    shuffle!(plaquettes)
    empty!(sep_plaquettes)
    for p1 in plaquettes
        compatible = true
        for p2 in sep_plaquettes
            plaquettesAreSeparated(p1, p2) || (compatible = false; break)
        end
        compatible && push!(sep_plaquettes, p1)
    end
    return plaquettes
end

function generateRandomPaths(InitialState, N, flipDepth; acceptanceRate = 0.5)
    startpath = empty(BitSet(1))

    nThreads = Threads.nthreads()

    # Conf_buffer = [copy(InitialState) for _ in 1:nThreads]
    AppendPaths_buffer = [empty([startpath]) for _ = 1:nThreads]
    # CurrentPaths_buffer = [empty([startpath]) for _ in 1:nThreads]

    # FlipPlaq_buffer = [empty([(0,0)]) for _ in 1:nThreads]

    batches = ChunkSplitters.chunks(1:N, n = nThreads, split = :batch)
    Threads.@threads for (iChunk, inds) in enumerate(batches)
        Conf = copy(InitialState)
        AppendPaths = AppendPaths_buffer[iChunk]

        FlipPlaq = empty([(0, 0)])
        path = copy(startpath)
        sepPlaquettes = Set(FlipPlaq)

        LI = LinearIndices(Conf)

        for i in inds
            for _ = 1:flipDepth
                Conf = spinConfig!(Conf, path, InitialState)
                getRandomSeparatedPlaquettes!(FlipPlaq, sepPlaquettes, Conf)
                for p in sepPlaquettes
                    if rand() < acceptanceRate
                        # path = flipPlaquette!(Conf,p)
                        pInt = LI[CartesianIndex(p)]

                        updatePath!(path, pInt)
                    end
                end
            end
            push!(AppendPaths, copy(path))
        end
    end
    return [
        LazyConfig(InitialState, paths) for appendPaths in AppendPaths_buffer for
        paths in appendPaths
    ]
end
