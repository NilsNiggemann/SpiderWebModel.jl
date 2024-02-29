struct PeriodicMatrix{T, Mat <: AbstractMatrix{T}} <: AbstractMatrix{T}
    UC::Mat
    Lx::Int
    Ly::Int
    offset::Int
end

PeriodicMatrix(Mat, Lx, Ly) = PeriodicMatrix(Mat, Lx, Ly, 0)
PeriodicMatrix(Mat) = PeriodicMatrix(Mat, size(Mat)..., 0)

function remapIndex(i, L)
    i = (i - 1) % L + 1
    if i <= 0
        i += L
    end
    return i
end

function remapIndices(x, y, Lx, Ly, offset)
    xRegion = (x - 1) ÷ Lx

    y = y - offset * xRegion

    x = remapIndex(x, Lx)
    y = remapIndex(y, Ly)
    return x, y
end

function remapIndices(x::UnitRange, y::UnitRange, Lx, Ly, offset)
    Inds = (remapIndices(i, j, Lx, Ly, offset) for i in x for j in y)
end

function Base.getindex(P::PeriodicMatrix, i::Integer, j::Integer)
    Lx, Ly = size(P.UC)
    i, j = remapIndices(i, j, Lx, Ly, P.offset)

    getindex(P.UC, i, j)
end

function Base.getindex(P::PeriodicMatrix, i, j)
    Lx, Ly = size(P.UC)
    checkbounds(Bool, P.UC, i, j) && return getindex(P.UC, i, j)
    inds = [CartesianIndex(remapIndices(ii, jj, Lx, Ly, P.offset)) for ii in i, jj in j]

    getindex(P.UC, inds)
end

function Base.setindex!(P::PeriodicMatrix, x, i::Integer)
    setindex!(P.UC, x, CartesianIndices(P)[i])
end
function Base.setindex!(P::PeriodicMatrix, x, i, j)
    checkbounds(Bool, P.UC, i, j) && return setindex!(P.UC,x, i, j)
    setindex!(P.UC, x, remapIndices(i, j, size(P.UC)..., P.offset)...)
end

Base.size(P::PeriodicMatrix) = (P.Lx, P.Ly)
Base.copy(P::PeriodicMatrix) = PeriodicMatrix(copy(P.UC), P.Lx, P.Ly, P.offset)

@inline function Base.checkbounds(::Type{Bool}, arr::PeriodicMatrix, I...)
    # i, j = remapIndices(I..., size(arr)..., arr.offset)
    # checkbounds(Bool, arr, i,j) || throw_boundserror(arr, ij)
    
    true
end

plaquetteIsInBounds(Conf::PeriodicMatrix, iCenter::Integer, jCenter::Integer) = true

function getPeriodicState(UC, Lx, Ly, offset = 0)
    Mat = PeriodicMatrix(UC, Lx, Ly, offset) |> Array
    S = SpinConfig(Mat, 1 / 2)
end
"""converts a state to a state with periodic boundary conditions"""
function getPeriodic(parent::SpinConfig)
    state = parent |> Array
    SpinConfig(PeriodicMatrix(state),parent.S)
end

function addFullTile!(Pij::SpinConfig, Tile::SpinConfig)
    addTile!(Pij, Tile)
    t22 = Tile[2, 2]
    isnan(t22) && return Pij
    Pij[2, 2] = Tile[2, 2]
    return Pij
end

function getXPath(Lx, Ly)
    correctPath!(xdirecPath(cld(Lx, 2), cld(Ly, 2)), SpinConfig(fill(NaN, Lx, Ly), 1 / 2))
end

function constructAllConfigs(Lx, Ly, PlaquetteList)
    LPx = cld(Lx, 2)
    LPy = cld(Ly, 2)
    path = xdirecPath(LPx, LPy)

    Mat = fill(NaN, Lx, Ly)

    El = PlaquetteList[begin]
    P = SpinConfig(Mat, El.S)

    correctPath!(path, P)

    AllConfigs_current = [zeros(Int, length(path))]

    AllConfigs_next = empty(AllConfigs_current)
    iter = 0
    TileListBuffer = collect(eachindex(PlaquetteList))

    while iter < lastindex(path)
        iter += 1
        i, j = path[iter]
        for history_buff in AllConfigs_current
            history = @view history_buff[1:(iter - 1)]
            P = reconstructTiling!(P, history, PlaquetteList, path)
            Pij = getPlaquette(P, i, j)
            Tiles = getFittingTiles!(TileListBuffer, Pij, PlaquetteList)
            for Tile in Tiles
                newhistory = copy(history_buff)
                newhistory[iter] = Tile
                push!(AllConfigs_next, newhistory)
            end
        end
        display(stack(AllConfigs_next))
        AllConfigs_current = AllConfigs_next
        AllConfigs_next = empty(AllConfigs_current)
        println("progress: ",
            iter,
            "/",
            length(path),
            " = ",
            round(iter * 100 / length(path), digits = 1),
            "% \tnumber of Configs: ",
            length(AllConfigs_current))
    end
    return AllConfigs_current
end

function applyStep!(P, tilingHistory, PlaquetteList, path, n)
    i, j = path[n]
    T = PlaquetteList[tilingHistory[n]]
    Pij = getPlaquette(P, i, j)
    addFullTile!(Pij, T)
    return P
end

function reconstructTiling!(P, tilingHistory, PlaquetteList, path)
    fill!(P, NaN)
    for i in eachindex(tilingHistory)
        P = applyStep!(P, tilingHistory, PlaquetteList, path, i)
    end
    return P
end

function reconstructTiling(tilingHistory, PlaquetteList, path)
    Lx = maximum(p[1] for p in path) + 1
    Ly = maximum(p[2] for p in path) + 1

    P = SpinConfig(fill(NaN, Lx, Ly), 1 / 2)
    path = correctPath(path, P)
    reconstructTiling!(P, tilingHistory, PlaquetteList, path)
end

function reconstructTiling_xDirec(Lx::Int, Ly::Int, tilingHistory, PlaquetteList)
    P = SpinConfig(fill(NaN, Lx, Ly), 1 / 2)
    LPx = cld(Lx, 2)
    LPy = cld(Ly, 2)
    path = correctPath!(xdirecPath(LPx, LPy), P)
    reconstructTiling!(P, tilingHistory, PlaquetteList, path)
end

function fillEmptyState(State::SpinConfig)
    NaNPos = findall(isnan, State)
    AllPossibilities = Iterators.product(((-State.S):(State.S) for i in NaNPos)...)
    newStates = [copy(State) for _ in 1:length(AllPossibilities)]
    for (i, conf) in enumerate(AllPossibilities)
        State = newStates[i]
        for (j, pos) in enumerate(NaNPos)
            State[pos] = conf[j]
        end
    end
    return newStates
end

function fillEmptyStates(States, Lx, Ly, PlaquetteList)
    newStates = empty([reconstructTiling_xDirec(Lx, Ly, first(States), PlaquetteList)])
    for State in States
        spinState = reconstructTiling_xDirec(Lx, Ly, State, PlaquetteList)
        FilledStates = fillEmptyState(spinState)
        append!(newStates, FilledStates)
    end
    return newStates
end

function getPeriodicState_history(Lx,
        Ly,
        ULx,
        ULy,
        tilingHistory,
        PlaquetteList,
        offset = 0)
    UC = reconstructTiling_xDirec(ULx, ULy, tilingHistory, PlaquetteList)
    State = getPeriodicState(UC, Lx, Ly, offset)
end

function isPeriodicTiling(Lx, Ly, tilingHistory, PlaquetteList, offset = 0)
    State = getPeriodicState_history((2 + 1 * offset) * Lx,
        (2 + 1 * offset) * Ly,
        Lx,
        Ly,
        tilingHistory,
        PlaquetteList,
        offset)
    return fulFillsConstraint_nonStrict(State)
end

function isPeriodicTiling(UC::AbstractMatrix, Lx, Ly, offset = 0)
    State = getPeriodicState(UC, Lx, Ly, offset)
    return fulFillsConstraint_nonStrict(State)
end

function getAllPeriodicStates(Lx, Ly)
    AllStates = constructAllConfigs(Lx, Ly, PlaquetteList)
end
