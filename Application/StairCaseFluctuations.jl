import SpiderWebModel as SW
using StaticArrays, CairoMakie
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
Flucs,_ = SW.generateAllFluctuations(Stair,SW.P1)
##
res = SW.getAllNeighborStates(Stair)



function adjFig(Flucs)
    dims1 = length(Flucs)
    dims2 = maximum(length.(Flucs))
    fig = Figure(size = 200 .*(dims2,dims1))
    S_ex = first(first(Flucs))
    suppresskwargs(i,j) = checkbounds(Bool,Flucs[i],j) ? (;) : (backgroundcolor = (:white,0.),xticklabelsvisible = false,yticklabelsvisible = false,xminorgridvisible = false,yminorgridvisible = false)
    axes = [Axis(fig[i,j];SW.getConfigAxis(S_ex)...,suppresskwargs(i,j)...) for i in 1:dims1,j in 1:dims2]
    for i in eachindex(Flucs)
        for j in eachindex(Flucs[i])
            SW.plotSpinConfig!(axes[i,j],Flucs[i][j])
            SW.plotApplPlaquettes!(axes[i,j],Flucs[i][j])
        end
    end
    fig
end
adjFig(Flucs)
##
SW.plotSpinConfig(res.AllConfigs[1])
plotApplPlaquettes!(current_axis(),res.AllConfigs[1])
display(current_figure())
for s in res.Nminus[1]
    SW.plotSpinConfig(res.AllConfigs[s]) 
    plotApplPlaquettes!(current_axis(),res.AllConfigs[s],SW.P1)
    plotApplPlaquettes!(current_axis(),res.AllConfigs[s],SW.P2,color = :blue)
    display(current_figure())
end
##
function plotPath()
    current = 1

    display(SW.plotApplPlaquettes(res.AllConfigs[current]))
    for i in 1:50
        neighbor = rand((res.Nminus,res.Nplus))
        isempty(neighbor[current]) && (i -=1;continue)
        current = rand(neighbor[current])
        display(SW.plotApplPlaquettes(res.AllConfigs[current]))
    end
end
plotPath()
##
Stair = getStairCase(12)
res = SW.getAllNeighborStates(Stair)

# SW.swapStates!(res.AllConfigs,res.Nplus,res.Nminus,1,2)
# SW.swapStates!(res.AllConfigs,res.Nplus,res.Nminus,2,12)
GC.gc()
##
μ = 0
H = SW.H(res.AllConfigs,res.Nplus,res.Nminus,μ)
# SW.testNplusMinus(res.Nplus,res.Nminus)
# @info "" SW.ishermitian(H) length(filter(!=(0),H - H'))

# heatmap(H,axis = (;aspect=1))
##
GC.gc()
@time sol = SW.SolveH(H)
GC.gc()
mag = [SW.getMagnetization(res.AllConfigs,sol,CartesianIndex(i,j)) for i in axes(Stair,1),j in axes(Stair,2)]
@info "" mavg = sum(abs,mag)/length(mag)
with_theme(theme_SimpleTicks()) do 
    fig,ax,hm = heatmap(mag;axis = (;aspect=1,title = L"μ = %$μ",SW.getConfigAxis(Stair)...,xticks = [1,6,12],yticks = [1,6,12]),figure = (;size = 0.8 .*(400,300)))
    Colorbar(fig[1,2],hm,label = L"\langle S^z \rangle")
    save("exactFig/mag_mu=$μ.png",fig)
    fig
end
##
using MakieHelpers

Sq = SW.getStructureFac(res.AllConfigs,sol)
##
# k = LinRange(0,2pi,40)
k = Sq.Sq[1].itp.ranges[1]
Sq_k = fetch.([Threads.@spawn Sq(kx,ky) for kx in k, ky in k])
with_theme(theme_PiTicks()) do 
    fig,ax,hm = heatmap(k,k,Sq_k;axis = (;aspect=1,title = L"μ = %$μ"),figure = (;size = 0.8 .*(400,300)))
    Colorbar(fig[1,2],hm,label = L"\mathcal{S}(\mathbf{q})")
    save("exactFig/Sq_mu=$μ.png",fig)
    fig
end