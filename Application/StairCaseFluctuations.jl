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

function plotApplPlaquettes!(ax,State,op;kwargs...)
    plaqs = SW.getApplicablePlaquettes_ns(State,op)
    points = Point2f.(plaqs)
    scatter!(ax,points,markersize = 13,color = :red;kwargs...)

end

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
            plotApplPlaquettes!(axes[i,j],Flucs[i][j],SW.P1)
            plotApplPlaquettes!(axes[i,j],Flucs[i][j],SW.P2,color = :blue)
            plotApplPlaquettes!(current_axis(),res.AllConfigs[s],SW.P2,color = :blue)
            # lines!(axes[i,j],1:10,color = :black)
        end
    end
    fig
end
adjFig(Flucs)
##
Stair = getStairCase(8)
res = SW.getAllNeighborStates(Stair)
##
SW.plotSpinConfig(res.AllConfigs[1])
plotApplPlaquettes!(current_axis(),res.AllConfigs[1],SW.P1)
plotApplPlaquettes!(current_axis(),res.AllConfigs[1],SW.P2,color = :blue)
display(current_figure())
for s in res.Nminus[1]
    SW.plotSpinConfig(res.AllConfigs[s]) 
    plotApplPlaquettes!(current_axis(),res.AllConfigs[s],SW.P1)
    plotApplPlaquettes!(current_axis(),res.AllConfigs[s],SW.P2,color = :blue)
    display(current_figure())
end
##
H = SW.H(res.AllConfigs,res.Nplus,res.Nminus,0)
heatmap(H,axis = (;aspect=1))

@info "" ishermitian(H) length(filter(!=(0),H - H'))