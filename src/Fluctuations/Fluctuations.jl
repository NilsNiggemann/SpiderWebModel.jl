using LinearAlgebra, SparseArrays, Arpack
import KrylovKit
const P1 = SA[-1.0 1.0 1.0;
              -1.0 0.0 -1.0;
              1.0 1.0 -1.0]
const P2 = -P1

const P1_SITES = -SVector(getSitesFromPlaquette(P1)) / 2
const P2_SITES = -SVector(getSitesFromPlaquette(P2)) / 2

const P1_SITES_BOOL = P1_SITES .== 1 / 2
const P2_SITES_BOOL = P2_SITES .== 1 / 2

function plaquetteFlippable(plaqSites::SVector{8})
    return plaqSites == P1_SITES || plaqSites == P2_SITES
end

function plaquetteFlippable(plaqSites::SVector{8, Union{Integer, Bool}})
    return plaqSites == P1_SITES_BOOL || plaqSites == P2_SITES_BOOL
end

function flipPlaquette!(Conf::AbstractMatrix, i, j)
    P = getPlaquette(Conf, i, j)

    P[1, 1] = -P[1, 1]
    P[1, 2] = -P[1, 2]
    P[1, 3] = -P[1, 3]
    P[2, 1] = -P[2, 1]
    # P[2,2] = -P[2,2]
    P[2, 3] = -P[2, 3]
    P[3, 1] = -P[3, 1]
    P[3, 2] = -P[3, 2]
    P[3, 3] = -P[3, 3]

    return Conf
end

function flipPlaquette!(Conf::AbstractMatrix{Bool}, i, j)
    P = getPlaquette(Conf, i, j)

    P[1, 1] = !P[1, 1]
    P[1, 2] = !P[1, 2]
    P[1, 3] = !P[1, 3]
    P[2, 1] = !P[2, 1]
    # P[2,2] = !P[2,2]
    P[2, 3] = !P[2, 3]
    P[3, 1] = !P[3, 1]
    P[3, 2] = !P[3, 2]
    P[3, 3] = !P[3, 3]

    return Conf
end

function flipPlaquette!(Conf::AbstractMatrix, pos::Tuple)
    i, j = pos
    flipPlaquette!(Conf, i, j)
end

function flipPlaquette!(Conf::AbstractMatrix, pos::CartesianIndex{2})
    i, j = Tuple(pos)
    flipPlaquette!(Conf, i, j)
end

function flipPlaquette!(Conf::AbstractMatrix, pos::Int)
    ij = CartesianIndices(Conf)[pos]
    flipPlaquette!(Conf, ij)
end

struct LazyConfig{T}
    parent::T
    path::BitSet
end

Base.hash(L::LazyConfig) = hash(L.path)
Base.isequal(L1::LazyConfig, L2::LazyConfig) = L1.path == L2.path
Base.:(==)(L1::LazyConfig, L2::LazyConfig) = L1.path == L2.path

Base.show(io::IO, ::MIME"text/plain", L::LazyConfig) = print(io, "LazyConfig: ", L.path)

function updatePath!(path, pos)
    if pos ∈ path
        delete!(path, pos)
    else
        push!(path, pos)
    end
    return path
end

function spinConfig!(Conf::AbstractMatrix, path, InitConf::AbstractMatrix)
    Conf .= InitConf
    for op in path
        flipPlaquette!(Conf, op)
    end
    return Conf
end

function spinConfig!(Conf::AbstractMatrix, L::LazyConfig)
    spinConfig!(Conf, L.path, L.parent)
end

function spinConfig(InitConf::AbstractMatrix, path)
    NewConf = copy(InitConf)
    spinConfig!(NewConf, path, InitConf)
end

function spinConfig(L::LazyConfig)
    spinConfig(L.parent, L.path)
end

function plotSpinConfig!(ax, L::LazyConfig, args...; kwargs...)
    plotSpinConfig(spinConfig(L), args...; kwargs...)
end
plotSpinConfig(L::LazyConfig; kwargs...) = plotSpinConfig(spinConfig(L), kwargs...)
function plotApplPlaquettes!(ax, L::LazyConfig; kwargs...)
    plotApplPlaquettes!(ax, spinConfig(L); kwargs...)
end
plotApplPlaquettes(L::LazyConfig; kwargs...) = plotApplPlaquettes(spinConfig(L); kwargs...)

function generateAllPaths(InitialState)
    startpath = empty(BitSet(1))
    # AllPaths = RobinDict(
    #    startpath => 1
    # )
    AllPaths = DataStructures.SwissDict(startpath => 1)

    # AllPaths = Dictionary(
    #    [startpath] , [1]
    # )
    NewPaths = [startpath]

    nThreads = Threads.nthreads()
    Conf_buffer = [copy(InitialState) for _ in 1:nThreads]
    AppendPaths_buffer = [empty([startpath]) for _ in 1:nThreads]

    while length(NewPaths) > 0
        batches = ChunkSplitters.chunks(NewPaths, n = nThreads, split = :batch)

        Threads.@threads for (iChunk, inds) in enumerate(batches)
            # for (path) in NewPaths
            Conf = Conf_buffer[iChunk]
            AppendPaths = AppendPaths_buffer[iChunk]
            empty!(AppendPaths)

            for i in inds
                path = NewPaths[i]
                Conf = spinConfig!(Conf, path, InitialState)
                getNewPaths!(AppendPaths, AllPaths, Conf, path)
            end
            unique!(AppendPaths)
        end
        empty!(NewPaths)

        for AppendPaths in AppendPaths_buffer
            appendPaths!(AllPaths, NewPaths, AppendPaths)
        end
    end
    # return Dict(spinConfig(InitialState,path) => num for (path,num) in AllPaths)
    return AllPaths
end

function appendPaths!(AllPaths::Dictionaries.Dictionary, NewPaths, AppendPaths)
    for path in AppendPaths
        haskey, token = gettoken!(AllPaths, path)
        if !haskey
            # AllPaths[path] = length(AllPaths)+1
            settokenvalue!(AllPaths, token, length(AllPaths))
            push!(NewPaths, path)
        end
    end
    return AllPaths, NewPaths
end

function appendPaths!(AllPaths, NewPaths, AppendPaths)
    for path in AppendPaths
        ind = get(AllPaths, path, 0)
        if ind == 0
            AllPaths[path] = length(AllPaths) + 1
            push!(NewPaths, path)
        end
    end
    return AllPaths, NewPaths
end

function getNewPaths!(AllNewPaths, AllPaths, Conf, path)
    plaqs = getApplicablePlaquettes(Conf)
    LI = LinearIndices(Conf)
    newpath = copy(path)
    for p in plaqs
        pInt = LI[CartesianIndex(p)]
        updatePath!(newpath, pInt)

        ind = get(AllPaths, newpath, 0)
        if ind == 0
            newpath2 = copy(newpath)
            push!(AllNewPaths, newpath2)
        end
        updatePath!(newpath, pInt) #undo the path change
    end

    return AllNewPaths
end

function getNeighbors(Conf, AllPaths, path)
    Neighbors = Int[]
    plaqs = getApplicablePlaquettes(Conf)
    LI = LinearIndices(Conf)
    newpath = copy(path)
    for p in plaqs
        pInt = LI[CartesianIndex(p)]
        updatePath!(newpath, pInt)

        ind = get(AllPaths, newpath, 0)
        updatePath!(newpath, pInt) # undo the path change
        if ind == 0
            error("new Config found")
        end

        push!(Neighbors, ind)
    end

    return Neighbors
end

function getAllNeighborStates(AllPaths, AllPaths_array, InitialState)
    Neighbors = [Int[] for i in eachindex(AllPaths)]
    nThreads = Threads.nthreads()
    Conf_buffer = [copy(InitialState) for _ in 1:nThreads]

    batches = ChunkSplitters.chunks(AllPaths_array, n = nThreads, split = :batch)
    Threads.@threads for (iChunk, inds) in enumerate(batches)
        # for (path) in NewPaths
        Conf = Conf_buffer[iChunk]
        for i in inds
            path = AllPaths_array[i]
            Conf = spinConfig!(Conf, path, InitialState)
            neighs = getNeighbors(Conf, AllPaths, path)
            Neighbors[i] = neighs
        end
    end

    return Neighbors
end

function getAllNeighborStates(StartConfig)
    AllConfigs = generateAllPaths(StartConfig)
    AllConfigsList = sortByValueOrder(AllConfigs)
    Neighbors = getAllNeighborStates(AllConfigs, AllConfigsList, StartConfig)

    AllStates = LazyConfig.(Ref(StartConfig), AllConfigsList)

    return (; AllStates, Neighbors)
end

function sortByValueOrder(D)
    ks = collect(keys(D))
    vals = collect(values(D))
    return ks[sortperm(vals)]
end

function sortbyKeyOrder(D)
    ks = collect(keys(D))
    vals = collect(values(D))
    return vals[sortperm(ks)]
end

function invertDict(D)
    D2 = Dict(v => k for (k, v) in D)
end

function H(AllStates, neighbors, mu::T = 0.0) where {T <: Number}
    dim = length(AllStates)
    rows = Int[]
    cols = Int[]
    vals = T[]
    function addTerm!(n, m, val)
        push!(rows, n)
        push!(cols, m)
        push!(vals, val)
    end

    for n in 1:dim
        for m in neighbors[n]
            addTerm!(n, m, -one(T))
        end

        val = mu * (length(neighbors[n]))
        addTerm!(n, n, val)
    end

    return Hermitian(sparse(rows, cols, vals))
end

SolveH(H, range = 1:1) = eigen(H, range)
const SparseMat = Union{SparseMatrixCSC, Hermitian{<:Number, <:SparseMatrixCSC}}

function SolveH(H::SparseMat; kwargs...)
    if size(H) == (1, 1)
        return (; values = [float(real(only(H)))], vectors = [1.0])
    end
    values, vectors = eigs(H, nev = 1, which = :SR, explicittransform = :none; kwargs...)
    return (; values, vectors)
end

function SolveHKrylov(H; kwargs...)
    values, vectors, _ = KrylovKit.eigsolve(H, 1, :SR; kwargs...)
    return (; values, vectors)
end

function flipSpinsAlongLine!(Conf, org, slope)
    slope ∈ (-Inf, Inf) && return flipSpinsAlongRow!(Conf, org[2])
    for i in axes(Conf.Mat, 1), j in axes(Conf.Mat, 2)
        if slope * (i - org[1]) == j - org[2]
            Conf[i, j] *= -1
        end
    end
    return Conf
end

function flipSpinsAlongDiagonal!(Conf, org, slope)
    j = org
    for i in axes(Conf.Mat, 1)
        j += slope
        if checkbounds(Bool, Conf, i, j)
            Conf[i, j] *= -1
        end
    end
    return Conf
end

function flipSpinsAlongRow!(Conf, i, skip = 2)
    Conf[1:2:end, i] .*= -1
    return Conf
end

function flipSpinsAlongCol!(Conf, i, skip = 2)
    Conf[i, 1:2:end] .*= -1
    return Conf
end

function CanApply(Conf::SpinConfig, Op::AbstractMatrix, i, j)
    plaquetteIsInBounds(Conf, i, j) || return false
    P = getPlaquette(Conf, i, j)

    P .+= Op

    applicable = fulFillsConstraint(Conf)

    P .-= Op

    return applicable
end

function canFlipPlaquette(Conf::SpinConfig, i, j)
    isodd(i + j) || return false
    plaquetteIsInBounds(Conf, i, j) || return false

    # P = getPlaquette(Conf,i,j)
    # sites = SVector(Mat[2,3],Mat[1,3],Mat[1,2],Mat[1,1],Mat[2,1],Mat[3,1],Mat[3,2],Mat[3,3])
    sites = SVector(Conf[i, j + 1],
        Conf[i - 1, j + 1],
        Conf[i - 1, j],
        Conf[i - 1, j - 1],
        Conf[i, j - 1],
        Conf[i + 1, j - 1],
        Conf[i + 1, j],
        Conf[i + 1, j + 1])
    return plaquetteFlippable(sites)
end

function CanApplyNonStrict(Conf::SpinConfig, Op, i, j)
    plaquetteIsInBounds(Conf, i, j) || return false
    isodd(i + j) || return false
    P = getPlaquette(Conf, i, j)

    P .+= Op

    applicable = all(x -> abs(x) <= Conf.S, P)

    P .-= Op

    return applicable
end

function CanApplyAnywhere(Conf::SpinConfig, Op::AbstractMatrix)
    a1 = axes(Conf.Mat, 1)
    a2 = axes(Conf.Mat, 2)

    Opx, Opy = size(Op)
    for i in a1, j in a2
        firstindex(a1) + Opx <= i <= lastindex(a1) - Opx || continue
        firstindex(a2) + Opy <= j <= lastindex(a2) - Opy || continue

        if CanApply(Conf, Op, i, j)
            return true
        end
    end
    return false
end

"""assumes that Op is already an allowed operator"""
function getApplicablePlaquettes(Conf::SpinConfig, Op)
    plaqPos = [(i, j) for i in axes(Conf.Mat, 1)
               for j in axes(Conf.Mat, 2) if CanApplyNonStrict(Conf, Op, i, j)]
    return plaqPos
end

function getApplicablePlaquettes(Conf::SpinConfig)
    plaqPos = [(i, j) for i in axes(Conf.Mat, 1)
               for j in axes(Conf.Mat, 2) if canFlipPlaquette(Conf, i, j)]
    return plaqPos
end
