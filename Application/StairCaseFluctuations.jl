import SpiderWebModel as SW
using StaticArrays, CairoMakie
using MakieHelpers

function getStairCase(L) 
    UC = SA[
        1 1 1 0;
        0 0 1 0;
        1 0 1 1;
        1 0 0 0;
    ]
    Mat = zeros(Float64,L,L)
    mel(x) = x==1 ? 1/2 : -1/2
    for i in axes(Mat,1)
        for j in axes(Mat,2)
            Mat[i,j] = mel(UC[i%4+1,j%4+1])
        end
    end

    S = SW.SpinConfig(Mat,1/2)

end
##
Stair = getStairCase(8)
##
@time res = SW.getAllNeighborStates(Stair)


##
function plotPath()
    current = 1

    display(SW.plotApplPlaquettes(res.AllStates[current]))
    for i in 1:20
        current = rand(res.Neighbors[current])
        display(SW.plotApplPlaquettes(res.AllStates[current]))
    end
end
plotPath()
##

##
Stair = getStairCase(12)
# @profview res = SW.getAllNeighborStates(Stair)
##
@time res = SW.getAllNeighborStates(Stair)

# SW.swapStates!(res.AllStates,res.Nplus,res.Nminus,1,2)
# SW.swapStates!(res.AllStates,res.Nplus,res.Nminus,2,12)
GC.gc()
##
μ = 1
@time H = SW.H(res.AllStates,res.Neighbors,μ)
# SW.testNplusMinus(res.Nplus,res.Nminus)
# @info "" SW.ishermitian(H) length(filter(!=(0),H - H'))

# heatmap(-H.data |> Array ,axis = (;aspect=1))
##
GC.gc()
@time sol = SW.SolveH(H)
GC.gc()
mag = [SW.getMagnetization(res.AllStates,sol,CartesianIndex(i,j)) for i in axes(Stair,1),j in axes(Stair,2)]
@info "" mavg = sum(abs,mag)/length(mag) mavgBulk = sum(abs,mag[3:end-2,3:end-2])/length(mag[3:end-2,3:end-2])

with_theme(theme_SimpleTicks()) do 
    fig,ax,hm = heatmap(mag;axis = (;aspect=1,title = L"μ = %$μ",SW.getConfigAxis(Stair)...,xticks = [1,6,12],yticks = [1,6,12]),figure = (;size = 0.8 .*(400,300)))
    Colorbar(fig[1,2],hm,label = L"\langle S^z \rangle")
    save("exactFig/mag_mu=$μ.png",fig)
    fig
end
##

Sq = SW.getStructureFac(SW.spinConfig.(res.AllStates),sol)
##
# k = LinRange(0,2pi,40)
k = Sq.Sq[1].itp.ranges[1]
Sq_k = fetch.([Threads.@spawn Sq(kx,ky) for kx in k, ky in k])
with_theme(theme_PiTicks()) do 
    fig,ax,hm = heatmap(k,k,Sq_k;axis = (;aspect=1,title = L"μ = %$μ"),figure = (;size = 0.8 .*(400,300)))
    Colorbar(fig[1,2],hm,label = L"\mathcal{S}^{zz}(\mathbf{q})")
    save("exactFig/Sq_mu=$μ.png",fig)
    fig
end

##
function EnergyScaling(L)
    Stair = getStairCase(L)
    res = SW.getAllNeighborStates(Stair)
    μ = 0
    H = SW.H(res.AllStates,res.Nplus,res.Nminus,μ)
    sol = SW.SolveH(H)
    return sol.values[1]
end

Ls = 4:1:14
Es = ([EnergyScaling(L) for L in Ls])

##
let 
    fig = Figure()
    ax = Axis(fig[1,1],xlabel = L"L",ylabel = L"-E_0(L)",xscale = log10,yscale = log10)
    scatterlines!(ax,Ls,-Es,label = L"Staircase$$")
    lines!(ax,Ls,Ls.^2 ./16,label = L"E_0 =-L²/16",linestyle = :dash, color = :grey)
    axislegend(ax,position = :rb)
    fig
    save("exactFig/EnergyScaling.png",fig)
end