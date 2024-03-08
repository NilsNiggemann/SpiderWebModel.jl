struct TilingMatrix{T,MatType<:AbstractMatrix{T}} <: AbstractMatrix{T}
    BasisStates::Vector{MatType}
    Tiles::Matrix{Int}
    TileSize::Tuple{Int,Int}
    size::Tuple{Int,Int}
end

Base.size(T::TilingMatrix) = T.size
Base.copy(T::TilingMatrix) = TilingMatrix(T.BasisStates, copy(T.Tiles), T.TileSize, T.size)

function TilingMatrix(BasisStates, Tiles::AbstractMatrix, shape::Tuple{Int,Int})
    TilingMatrix(BasisStates, Tiles, size(BasisStates[begin]), shape)
end

@inline function Base.getindex(T::TilingMatrix, i, j)
    Lx, Ly = T.TileSize

    iChunk = (i - 1) ÷ Lx + 1
    jChunk = (j - 1) ÷ Ly + 1
    i = (i - 1) % Lx + 1
    j = (j - 1) % Ly + 1

    Tile = T.BasisStates[T.Tiles[iChunk, jChunk]]

    return Tile[i, j]
end

function constructRandomTiling!(Mat, Tiles, Lx, Ly)
    rand!(Mat, eachindex(Tiles))
    return TilingMatrix(Tiles, Mat, size(Tiles[begin]), (Lx, Ly))
end

function constructGSFromTiles(Tiles, Lx, Ly; numTries = 100000)
    Tile_Lx, Tile_Ly = size(Tiles[begin])

    TileMat = zeros(Int, cld(Lx, Tile_Lx), cld(Ly, Tile_Ly))

    Tiling = constructRandomTiling!(TileMat, Tiles, Lx, Ly)

    Conf = SpinConfig(Tiling, 1 / 2)
    GS = typeof(Conf)[]

    for i = 1:numTries
        Tiling = constructRandomTiling!(TileMat, Tiles, Lx, Ly)

        if fulFillsConstraint_nonStrict(Conf)
            push!(GS, copy(Conf))
        end
    end
    return GS
end

function constructGSFromTiles_Threads(Tiles, Lx, Ly; numTries = 100000)
    StatesExample = constructGSFromTiles(Tiles, Lx, Ly, numTries = 1)

    nThreads = Threads.nthreads()
    StatesCollection = Vector{typeof(StatesExample)}(undef, nThreads)

    Threads.@threads for i = 1:nThreads
        StatesCollection[i] =
            constructGSFromTiles(Tiles, Lx, Ly, numTries = numTries ÷ nThreads)
    end
    return append!(StatesCollection...)
end
