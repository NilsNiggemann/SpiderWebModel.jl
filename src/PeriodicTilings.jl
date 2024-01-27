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


function addFullTile!(Pij::SpinConfig,Tile::SpinConfig,)
    addTile!(Pij,Tile)
    t22 = Tile[2,2]
    isnan(t22) && return Pij
    Pij[2,2] = Tile[2,2]
    return Pij
end

getXPath(Lx,Ly) = correctPath!(xdirecPath(cld(Lx,2),cld(Ly,2)),SpinConfig(fill(NaN,Lx,Ly),1/2))


function constructAllConfigs(Lx,Ly,PlaquetteList)
    LPx = cld(Lx,2)
    LPy = cld(Ly,2)
    path = xdirecPath(LPx,LPy)
    
    Mat = fill(NaN,Lx,Ly)
    
    El = PlaquetteList[begin]
    P = SpinConfig(Mat,El.S)

    correctPath!(path,P)
    
    AllConfigs_current = [zeros(Int,length(path))]

    AllConfigs_next = empty(AllConfigs_current)
    iter = 0
    TileListBuffer = collect(eachindex(PlaquetteList))

    while iter < lastindex(path)
        iter += 1
        i,j = path[iter]
        for history_buff in AllConfigs_current
            history = @view history_buff[1:iter-1]
            P = reconstructTiling!(P,history,PlaquetteList,path)
            Pij = getPlaquette(P,i,j)
            Tiles = getFittingTiles!(TileListBuffer,Pij,PlaquetteList)
            for Tile in Tiles
                newhistory = copy(history_buff)
                newhistory[iter] = Tile
                push!(AllConfigs_next,newhistory)
            end
        end
        display(stack(AllConfigs_next))
        AllConfigs_current = AllConfigs_next
        AllConfigs_next = empty(AllConfigs_current)
        println("progress: ", iter, "/",length(path), " = ", round(iter*100/length(path),digits = 1), "% \tnumber of Configs: ",length(AllConfigs_current) )
    end
    return AllConfigs_current
end

# function constructAllConfigs_fast(Lx,Ly,PlaquetteList)
#     LPx = cld(Lx,2)
#     LPy = cld(Ly,2)
#     path = xdirecPath(LPx,LPy)
    
#     Mat = fill(NaN,Lx,Ly)
    
#     El = PlaquetteList[begin]
#     P = SpinConfig(Mat,El.S)

#     correctPath!(path,P)
    
#     path_length = length(path)

#     CurrentConfigs = ElasticArray{Int}(undef, path_length, length(PlaquetteList))
#     fill!(CurrentConfigs,0)
#     NextConfigs = copy(CurrentConfigs)
#     CurrentConfigs[1,:] .= eachindex(PlaquetteList)

#     iter = 1
#     TileListBuffer = collect(eachindex(PlaquetteList))
    
#     numConfigs = length(PlaquetteList)

#     while iter < lastindex(path)
#         iter += 1
#         i,j = path[iter]

#         CurrentConfigs = @view CurrentConfigs[1:iter-1,1:numConfigs]

#         @info "" size(CurrentConfigs)
        
#         for (i,history) in enumerate(eachcol(CurrentConfigs))
#             println(history)
#             P = reconstructTiling!(P,history,PlaquetteList,path)
#             Pij = getPlaquette(P,i,j)
#             Tiles = getFittingTiles!(TileListBuffer,Pij,PlaquetteList)
            
#             numConfigs_old = size(CurrentConfigs,2)
#             numConfigs_new = numConfigs_old + length(Tiles)

#             resize!(NextConfigs,path_length,numConfigs_new)
#             fill!(NextConfigs,0)

#             # println("resizing from", (path_length,numConfigs_old)," to ",size(CurrentConfigs))
            
#             for Tile in Tiles
#                 NextConfigs[1:iter-1,i] .= history
#                 numConfigs_old += 1
#             end
#             CurrentConfigs[1:iter-1,numConfigs_old+1:end] .= history
            
#             nextConfig = @view CurrentConfigs[iter+1,numConfigs_old+1:end]
#             nextConfig .= Tiles
#             return CurrentConfigs
#         end
#         println("progress: ", iter, "/",path_length, " = ", round(iter*100/path_length,digits = 1), "% \tnumber of Configs: ",size(CurrentConfigs,2) )
#     end
#     return CurrentConfigs
# end

function applyStep!(P,tilingHistory,PlaquetteList,path,n)
    i,j = path[n]
    T = PlaquetteList[tilingHistory[n]]
    Pij = getPlaquette(P,i,j)
    addFullTile!(Pij,T)
    return P
end

function reconstructTiling!(P,tilingHistory,PlaquetteList,path)
    fill!(P,NaN)
    for i in eachindex(tilingHistory)
        P = applyStep!(P,tilingHistory,PlaquetteList,path,i)
    end
    return P
end

function reconstructTiling(tilingHistory,PlaquetteList,path)
    Lx = maximum(p[1] for p in path)+1
    Ly = maximum(p[2] for p in path)+1
    
    P = SpinConfig(fill(NaN,Lx,Ly),1/2)
    path = correctPath(path,P)
    reconstructTiling!(P,tilingHistory,PlaquetteList,path)
end

function reconstructTiling_xDirec(Lx::Int,Ly::Int,tilingHistory,PlaquetteList)
    P = SpinConfig(fill(NaN,Lx,Ly),1/2)
    LPx = cld(Lx,2)
    LPy = cld(Ly,2)
    path = correctPath!(xdirecPath(LPx,LPy),P)
    reconstructTiling!(P,tilingHistory,PlaquetteList,path)
end


