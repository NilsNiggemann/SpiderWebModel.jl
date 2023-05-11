import SpiderWebModel as SW
using SpiderWebModel.SpinFRGLattices
using CairoMakie
using SpiderWebModel.StaticArrays
##
Lat = SW.generateLattice(4)
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

    for (ax,Mat) in zip(axes,AllAllowedConfigs)
        heatmap!(ax,xticks,yticks,Array(Mat.x),colorrange = (-0.5,0.5))
    end
    save("combs.pdf",fig,)
    fig
end

## plot rot inv combinations

"""given a 3×3 matrix, return a 3×3 matrix with the same entries but rotated 90 degrees around the middle element"""
function rotl90(A)
    ind1, ind2 = axes(A)
    B = similar(A, (ind2,ind1))
    n = first(ind2)+last(ind2)
    for i=axes(A,1), j=ind2
        B[n-j,i] = A[i,j]
    end
    return B
end

function rotl90(A::SMatrix{3,3})
    A = A'
    return SMatrix{3,3}(
        A[3] , A[6] ,  A[9],
        A[2] , A[5] ,  A[8],
        A[1] , A[4] ,  A[7]
    )'
end

rotl90(S::SW.Plaquette) = SW.Plaquette(rotl90(S.x))
let 
    FilteredConfigs = empty(AllAllowedConfigs)

    
    r1 = rotl90
    r2 = rotl90 ∘ rotl90
    r3 = rotl90 ∘ rotl90 ∘ rotl90

    for c in AllAllowedConfigs
        confs = Set([c,r1(c),r2(c),r3(c)])
        if isempty(confs ∩ FilteredConfigs)
            push!(FilteredConfigs,c)
        end
    end
    
    sizex = Int(floor(sqrt(length(FilteredConfigs))))
    sizey = Int(ceil(sqrt(length(FilteredConfigs))))+1
    
    sizex = 3
    sizey = 7
    
    layout = zeros(sizex,sizey)

    fig = Figure(resolution = 200 .* (sizex,sizey))
    allij = Tuple.(CartesianIndices(layout))
    xticks = yticks = [-1,0,1]
    axes = [Axis(fig[fj,fi]; xticklabelsvisible = false, yticklabelsvisible = false,xgridvisible = false,ygridvisible = false,aspect = 1,xticks ,yticks ) for (i,(fi,fj)) in enumerate(allij)]
    # return length(axes),length(FilteredConfigs)
    # @assert length(axes) == length(FilteredConfigs) 
    for (ax,Mat) in zip(axes,FilteredConfigs)
        heatmap!(ax,xticks,yticks,Array(rotl90(Mat.x)),colorrange = (-0.5,0.5))
    end
    save("combsReduced.pdf",fig,)
    fig
end

##
using CairoMakie.Makie.ColorSchemes
function plotSpinConfig!(ax,S;kwargs...)
    rPoints = [Point2f(getCartesian(R,SW.Basis)) for R in keys(S.Spins)]
    vals = collect(values(S.Spins))
    scatter!(ax,rPoints,colormap = :greys,marker = :rect ,color = vals,markersize = 40;kwargs...) 
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
    S = SW.SpinConfig(Rvec_2D)
    L = 20
    Points = generateLUnitCells(L,SW.Basis)
    filter!(x-> x.b == 1 && isInLBox(x,L),Points )

    tileNums = rand(eachindex(AllAllowedConfigs),length(Points))
    @time for (R,tilenum) in zip(Points,tileNums)
        SW.setTile!(S,R,AllAllowedConfigs[tilenum])
    end
    S
end
plotSpinConfig(SpinConfig)
##

##
AFM = let 
    Rvecs = generateLUnitCells(10,SW.Basis)
    filter!(x->isInLBox(x,10),Rvecs )
    S = SW.SpinConfig(Rvec_2D)
    for R in Rvecs
        S[R] = R.b == 1 ? 1/2 : -1/2
    end
    S
end
plotSpinConfig(AFM)
##
using FRGLatticeEvaluation
let 
    Sij = SW.getCorrel(SpinConfig)
    Rij_vec = [SW.getCartesian(Ri)-SW.getCartesian(Rj) for Ri in keys(SpinConfig.Spins),Rj in keys(SpinConfig.Spins)]
    Sq = FourierStruct(Sij[:],Rij_vec[:],length(Rij_vec))
    k, sq = Fourier2D(Sq,xyplane,ext = 2pi,res = 100)
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1)
    heatmap!(ax,k ./pi,k ./pi,sq)
    fig
end