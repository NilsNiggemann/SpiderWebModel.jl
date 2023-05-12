import SpiderWebModel as SW
using SpiderWebModel.SpinFRGLattices
using CairoMakie
using SpiderWebModel.StaticArrays
##
Lat = generateLUnitCells(4,SW.Basis)
R1 = Rvec(0,0,1)
let 
    P = SW.getPlaquette(R1)
    all(SW.isInPlaquetteOf.(Ref(R1),SW.getSitesFromPlaquette(P)))
    
end
##

let 
    P = SW.getPlaquette(R1)
    ri = Array([Point2f(SW.getCartesian(R)) for R in SW.getSitesFromPlaquette(P)] )
    # scatter(ri)
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1)
    scatter!(ax,ri,markersize = collect(5 .+ 3 .* eachindex(ri)))
    fig
end
##
function isInLBox(R,L)
    r = SW.getCartesian(R)
    -L <= r[1] <= L && -L <= r[2] <= L
end

function scatterRvecs!(ax,RV,Basis;kwargs...)
    NU = Basis.NUnique
    ri = [Point2f(getCartesian(R,Basis)) for R in RV]
    cols(R) = R.b == 1 ? :red : :black
    msize(R) = R.b == 1 ? 13 : 10
    
    scatter!(ax,ri,markersize = msize.(RV),color = cols.(RV);kwargs...)
end

function isOnBorder(R::Rvec_2D,L)
    r = getCartesian(R,SW.Basis)
    abs(r[1]) == L || abs(r[2]) == L
end

let 
    # sites = SW.generateLattice(10)
    sites = generateLUnitCells(4,SW.Basis)
    filter!(x-> isInLBox(x,2),sites )
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1)
    scatterRvecs!(ax,sites,SW.Basis)
    fig
end
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
    heatmap!(ax,Array(S.Mat),colorrange = (-S.S,S.S),color = :RdBu_11)
end
function plotSpinConfig(S;kwargs...) 
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1,backgroundcolor = :grey)
    plotSpinConfig!(ax,S;kwargs...)
    return fig
end
##
SpinConfig = let 
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1,backgroundcolor = :grey)
    S = SW.SpinConfig2(-0.5ones(20,20),1/2)
    S[10,10] = 0.5
    S
end
SW.flipSpinsAlongLine!(SpinConfig,(10,14),1)
plotSpinConfig(SpinConfig)
##

AFM = let 
    S = SW.SpinConfig2([0.5*(-1)^(i+j) for i in 1:40,j in 1:40],1/2)
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

    S = SW.SpinConfig(Mat)

end

function getApplicablePlaquettes(SpinConfig)
    plaqPos = [R for (R,S) in SpinConfig.Spins if SW.CanApplyPlaquette(SpinConfig,R)]
    return plaqPos
end

let 
    StairCase = getStairCase(30)
    fig = plotSpinConfig(StairCase,markersize = 23)
    ax = current_axis()
    plaqs = getApplicablePlaquettes(StairCase)
    points = Point2f.(getCartesian.(plaqs,Ref(SW.Basis)))
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