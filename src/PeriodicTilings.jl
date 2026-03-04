function getPeriodicState(UC, Lx, Ly, offset = 0)
    Mat = PeriodicMatrix(UC, Lx, Ly, offset) #|> Array
    S = SpinConfig(Mat, 1 / 2)
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
            history = @view history_buff[1:(iter-1)]
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
        println(
            "progress: ",
            iter,
            "/",
            length(path),
            " = ",
            round(iter * 100 / length(path), digits = 1),
            "% \tnumber of Configs: ",
            length(AllConfigs_current),
        )
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

function _getSpinFromPlaquetteList(PlaquetteList)
    return maximum(y -> maximum(x->isnan(x) ? 0 : abs(x), y), PlaquetteList)
end
function reconstructTiling(tilingHistory, PlaquetteList, path)
    Lx = maximum(p[1] for p in path) + 1
    Ly = maximum(p[2] for p in path) + 1
    Spin = _getSpinFromPlaquetteList(PlaquetteList)
    P = SpinConfig(fill(NaN, Lx, Ly), Spin)
    path = correctPath(path, P)
    reconstructTiling!(P, tilingHistory, PlaquetteList, path)
end

function reconstructTiling_xDirec(Lx::Int, Ly::Int, tilingHistory, PlaquetteList)
    Spin = _getSpinFromPlaquetteList(PlaquetteList)
    P = SpinConfig(fill(NaN, Lx, Ly), Spin)
    LPx = cld(Lx, 2)
    LPy = cld(Ly, 2)
    path = correctPath!(xdirecPath(LPx, LPy), P)
    reconstructTiling!(P, tilingHistory, PlaquetteList, path)
end

function fillEmptyState(State::SpinConfig)
    NaNPos = findall(isnan, State)
    AllPossibilities = Iterators.product(((-State.S):(State.S) for i in NaNPos)...)
    newStates = [copy(State) for _ = 1:length(AllPossibilities)]
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

function getPeriodicState_history(
    Lx,
    Ly,
    ULx,
    ULy,
    tilingHistory,
    PlaquetteList,
    offset = 0,
)
    UC = reconstructTiling_xDirec(ULx, ULy, tilingHistory, PlaquetteList)
    State = getPeriodicState(UC, Lx, Ly, offset)
end

function isPeriodicTiling(Lx, Ly, tilingHistory, PlaquetteList, offset = 0)
    State = getPeriodicState_history(
        (2 + 1 * offset) * Lx,
        (2 + 1 * offset) * Ly,
        Lx,
        Ly,
        tilingHistory,
        PlaquetteList,
        offset,
    )
    return fulFillsConstraint_nonStrict(State)
end

function isPeriodicTiling(UC::AbstractMatrix, Lx, Ly, offset = 0)
    State = getPeriodicState(UC, Lx, Ly, offset)
    return fulFillsConstraint_nonStrict(State)
end