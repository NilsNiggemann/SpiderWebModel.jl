function tileMatrix!(Mat,UC)
    Lx,Ly = size(UC)
    for i in axes(Mat,1)
        for j in axes(Mat,2)
            Mat[i,j] = UC[i%Lx+1,j%Ly+1]
        end
    end

    S = SpinConfig(Mat,1/2)
end

function getTilings(Lx,Ly)
    UC = zeros(Float64,Lx,Ly)
    Configs = empty([SpinConfig(UC,1/2)])

    for i in eachindex(UC)
    end
    tileMatrix!(Mat,UC)
end

function constructAllConfigs(LPx,LPy,PlaquetteList)
    path = xdirecPath(LPx,LPy);
    
    Lx = 2*LPx+1
    Ly = 2*LPy+1
    Mat = fill(NaN,Lx,Ly)
    
    El = PlaquetteList[begin]
    P = SpinConfig(Mat,El.S)
    filter!(x -> plaquetteIsInBounds(P,x...),path)
    
    indices = shuffle(collect(eachindex(PlaquetteList)))
    
    AllConfigs_current = [P]
    AllConfigs_next = empty(AllConfigs_current)

    while iter < lastindex(path)-1
        
        iter = lastindex(tilingHistory)
        
        i,j = path[iter+1]

        for P in AllConfigs
            Pij = getPlaquette(P,i,j)
            Tiles = getFittingTiles(Pij,PlaquetteList)
            for tile in Tiles
                Pnew = copy(P)
                addTile!(Pnew,tile)
                push!(AllConfigs_next,Pnew)
            end
        end
        AllConfigs_current = AllConfigs_next
        AllConfigs_next = empty(AllConfigs_current)
    end
    return AllConfigs_current
end