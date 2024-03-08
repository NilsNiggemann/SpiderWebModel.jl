function xdirecPath(LPx, LPy = LPx)
    [(i, j) for j = 1:(2LPy+1) for i = 1:(2LPx+1) if iseven(i + j)]
end
function xdirecPathReverse(LPx, LPy = LPx)
    [(i, j) for j = 1:(2LPy+1) for i = (2LPx+1):-1:1 if iseven(i + j)]
end
function ydirecPath(LPx, LPy = LPx)
    [(i, j) for i = 1:(2LPx+1) for j = 1:(2LPy+1) if iseven(i + j)]
end
function ydirecPathReverse(LPx, LPy = LPx)
    [(i, j) for i = 1:(2LPx+1) for j = (2LPy+1):-1:1 if iseven(i + j)]
end

function correctPath!(path, Config)
    filter!(x -> plaquetteIsInBounds(Config, x...), path)
end

function correctPath(path, Config)
    filter(x -> plaquetteIsInBounds(Config, x...), path)
end

function spiralPath(L)
    num_points = 4L^2
    coords = [(0, 0)]
    x, y = 0, 0
    dx, dy = 1, 0
    side_length = 1
    steps_in_side = 0

    for i = 2:num_points
        x += dx
        y += dy
        push!(coords, (x, y))
        steps_in_side += 1

        if steps_in_side == side_length
            steps_in_side = 0
            if dx == 1 && dy == 0
                dx, dy = 0, 1
            elseif dx == 0 && dy == 1
                dx, dy = -1, 0
                side_length += 1
            elseif dx == -1 && dy == 0
                dx, dy = 0, -1
            elseif dx == 0 && dy == -1
                dx, dy = 1, 0
                side_length += 1
            end
        end
    end
    filter!(x -> iseven(x[1] + x[2]), coords)
    for i in eachindex(coords)
        coords[i] = coords[i] .+ (L, L)
    end
    return coords
end

function spiralPath(Lx, Ly)
    @assert Lx == Ly "Lx must be equal to Ly"
    spiralPath(Lx)
end

function setupCalc!(path, LPx, LPy, PlaquetteList)
    Lx = 2 * LPx + 1
    Ly = 2 * LPy + 1
    Mat = fill(NaN, Lx, Ly)
    El = PlaquetteList[begin]
    P = SpinConfig(Mat, El.S)
    filter!(x -> plaquetteIsInBounds(Mat, x...), path)

    emptyTilesList = getFreeTilesPath(P, path, PlaquetteList)

    # TilesDict = getFittingTilesDict(PlaquetteList)
    (; path, emptyTilesList)
end

function getStepDeleter(L, default = 5, tries = 5)
    function deleter(x, counter)
        counter[] += 1
        if counter[] > tries
            counter[] = 0
            return L
        end
        return default
    end
end

function getFreeTilesPath(P, path, PlaquetteList)
    # emptyTiles = Vector{Int}[]
    emptyTiles = Vector{CartesianIndex{2}}[]
    newP = copy(P)
    newP .= NaN
    testPlaq = first(PlaquetteList)
    for (i, j) in path
        Pij = getPlaquette(newP, i, j)
        addTile!(Pij, testPlaq)
        emptyTilesCurrent = findall(isnan, newP)
        push!(emptyTiles, emptyTilesCurrent)
    end
    return emptyTiles
end
