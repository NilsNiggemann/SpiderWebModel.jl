import SpiderWebModel as SW
using CairoMakie
using SpiderWebModel.StaticArrays
## plot 70 combinations

AllAllowedConfigs = SW.getAllGS(0.5)
##
let 
    sizex = Int(7)
    sizey = Int(10)

    layout = zeros(sizex,sizey)

    fig = Figure(resolution = 200 .* (sizex,sizey))
    allij = Tuple.(CartesianIndices(layout))
    xticks = yticks = [-1,0,1]
    axes = [Axis(fig[fj,fi]; xticklabelsvisible = false, yticklabelsvisible = false,xgridvisible = false,ygridvisible = false,aspect = 1,xticks ,yticks ) for (i,(fi,fj)) in enumerate(allij)]

    for (ax,Plaq) in zip(axes,AllAllowedConfigs)
        heatmap!(ax,xticks,yticks,Array(Plaq.Mat),colorrange = (-0.5,0.5))
    end
    # save("combs.pdf",fig,)
    fig
end

##
using CairoMakie.Makie.ColorSchemes
function plotSpinConfig!(ax,S;kwargs...)
    hm = heatmap!(ax,Array(S.Mat),colorrange = (-S.S,S.S),color = :RdBu_11)
    Color
end
function plotSpinConfig(S;kwargs...) 
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1,backgroundcolor = :grey)
    plotSpinConfig!(ax,S;kwargs...)
    return fig
end
function drawPlaquette!(ax,(i,j);kwargs...)
    points = Point2f.([(i-1,j-1),(i+1,j-1),(i+1,j+1),(i-1,j+1),(i-1,j-1)])
    lines!(ax,points;linestyle = :dash,color = :white,kwargs...)
end
drawPlaquette!((i,j);kwargs...) = drawPlaquette!(current_axis(),(i,j);kwargs...)
##

SpinConfig = let 
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1,backgroundcolor = :grey)
    S = SW.SpinConfig(-0.5ones(20,20),1/2)
    S[10,10] = 0.5
    S
end
SW.flipSpinsAlongLine!(SpinConfig,(10,14),1)
plotSpinConfig(SpinConfig)
##

AFM = let 
    S = SW.SpinConfig([0.5*(-1)^(i+j) for i in 1:40,j in 1:40],1/2)
end

SW.PlaquetteOperatorSave!(AFM,5,5)
plotSpinConfig(AFM)
##
SW.fulFillsConstraint(AFM)
##
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

function getApplicablePlaquettes(SpinConfig)
    plaqPos = [(i,j) for i in axes(SpinConfig.Mat,1) for j in axes(SpinConfig.Mat,2) if SW.CanApplyPlaquette(SpinConfig,i,j)]
    return plaqPos
end

let 
    StairCase = getStairCase(30)
    fig = plotSpinConfig(StairCase,markersize = 23)
    ax = current_axis()
    plaqs = getApplicablePlaquettes(StairCase)
    return plaqs
    points = Point2f.(plaqs)
    scatter!(ax,points,markersize = 13,color = :red)
    fig
end
##
let 
    StairCase = getStairCase(30)
    plaqs = getApplicablePlaquettes(StairCase)
    for R in plaqs
        SW.PlaquetteOperatorSave!(StairCase,R)
        break
    end
    plotSpinConfig(StairCase,markersize = 23)
end
##

function generateRandomGroundState(L,S=1/2;maxiter = 10000)
    AllowedPlaquettes = SW.getAllGS(S)
    Mat = rand(eachindex(AllowedPlaquettes),L,L)
    SC = SW.SpinConfig(Mat,S)
    sites_i = rand(1:L,maxiter)
    sites_j = rand(1:L,maxiter)

    for (i,j) in zip(sites_i,sites_j)
        P = SW.getPlaquette(SC,i,j)
        P .= -P
        SW.fulFillsConstraint(SC) && return SC
    end
    return nothing
    error("Could not find a ground state")
end
