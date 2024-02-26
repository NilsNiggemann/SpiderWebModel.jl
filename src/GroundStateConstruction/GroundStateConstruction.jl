_isZeroOrNaN(x) = x == 0 || isnan(x)

function addTile!(Pij::SpinConfig, Tile::SpinConfig)
    Pij[1, 1] = Tile[1, 1]
    Pij[1, 2] = Tile[1, 2]
    Pij[1, 3] = Tile[1, 3]

    Pij[2, 1] = Tile[2, 1]
    #middle of tile is empty
    Pij[2, 3] = Tile[2, 3]

    Pij[3, 1] = Tile[3, 1]
    Pij[3, 2] = Tile[3, 2]
    Pij[3, 3] = Tile[3, 3]
    return Pij
end

function canPlaceTile(P::SpinConfig, T1)
    for (p, t) in zip(P, T1)
        x = p - t
        if !(_isZeroOrNaN(x))
            return false
        end
    end
    return true
end

function canPlaceTile(P::Tuple, T1)
    Tsites = getSitesFromPlaquette(T1)
    return all(_isZeroOrNaN, p - t for (p, t) in zip(P, Tsites))
end

function getFittingTiles(P, PlaquetteList)
    [i for (i, t) in enumerate(PlaquetteList) if canPlaceTile(P, t)]
end

function getFittingTiles!(TileList, P, PlaquetteList)
    len = 0
    for (i, t) in enumerate(PlaquetteList)
        if canPlaceTile(P, t)
            len += 1
            TileList[len] = i
        end
    end
    return @view TileList[1:len]
end

function getFittingTilesDict!(fittingTilesDict, P, PlaquetteList)::Vector{Int}
    tilenums = get(fittingTilesDict, P, nothing)
    if tilenums !== nothing
        return tilenums
    else
        newtiles = getFittingTiles(P, PlaquetteList)
        fittingTilesDict[P] = newtiles
        return newtiles
    end
end

function constructConfigPath(P::SpinConfig, PlaquetteList, setup;
        maxiter = 10000,
        deleteSteps = getStepDeleter(size(P, 2) + 2, 1, 15),
        verbose = true,
        plotSteps = false)
    path = setup.path
    emptyTilesList = setup.emptyTilesList

    # path = spiralPath(LPx)
    tilingHistory = Int[]
    sizehint!(tilingHistory, length(path))
    iter = 1
    TotIter = 1

    function applyStep!(P, n)
        i, j = path[n]
        T = PlaquetteList[tilingHistory[n]]
        Pij = getPlaquette(P, i, j)
        addTile!(Pij, T)
        return P
    end

    function resetToStep!(P, iterNum)
        emptyTiles = emptyTilesList[iterNum + 1]
        for i in eachindex(emptyTiles)
            I = emptyTiles[i]
            P[I] = NaN
        end
    end

    P_init = getSitesFromPlaquette(getPlaquette(P, path[lastindex(tilingHistory) + 1]...))
    fittingTilesDict = Dict(P_init => getFittingTiles(P_init, PlaquetteList))

    counter = Ref(0)
    successfulplacements = 0

    while iter < lastindex(path) - 1
        iter = lastindex(tilingHistory)

        i, j = path[iter + 1]
        TotIter += 1
        if TotIter > maxiter
            if verbose
                @warn "maxiter reached"
            end
            break
        end
        Pij = getPlaquette(P, i, j)
        PijSites = getSitesFromPlaquette(Pij)

        tileList = getFittingTilesDict!(fittingTilesDict, PijSites, PlaquetteList)
        if isempty(tileList)
            deleteat!(tilingHistory, max(2, iter - deleteSteps(iter, counter)):iter)
            iter = lastindex(tilingHistory)
            resetToStep!(P, iter)
            successfulplacements = 0
            continue
        end
        iT = rand(tileList)
        T = PlaquetteList[iT]
        addTile!(Pij, T)
        push!(tilingHistory, iT)
        if successfulplacements > 1
            counter[] = 0
        end
        successfulplacements += 1
        plotSteps && display(plotSpinConfig(P))
    end
    return P
end

function constructConfigPath(LPx, LPy, PlaquetteList,
        setup = setupCalc!(xdirecPath(LPx, LPy), LPx, LPy, PlaquetteList);
        kwargs...)
    Lx = 2 * LPx + 1
    Ly = 2 * LPy + 1
    Mat = fill(NaN, Lx, Ly)

    El = PlaquetteList[begin]
    P = SpinConfig(Mat, El.S)
    return constructConfigPath(P, PlaquetteList, setup; kwargs...)
end

function constructConfigPath(LPx, LPy, PlaquetteList, path::AbstractVector; kwargs...)
    setup = setupCalc!(copy(path), LPx, LPy, PlaquetteList)
    constructConfigPath(LPx, LPy, PlaquetteList, setup; kwargs...)
end

function constructConfigPath(LPx, LPy, PlaquetteList, pathFunc::Function; kwargs...)
    path = pathFunc(LPx, LPy)
    setup = setupCalc!(path, LPx, LPy, PlaquetteList)
    constructConfigPath(LPx, LPy, PlaquetteList, setup; kwargs...)
end
