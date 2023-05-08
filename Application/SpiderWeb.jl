import SpiderWebModel as SW
using SpiderWebModel.SpinFRGLattices
using CairoMakie
using SpiderWebModel.StaticArrays
##
Lat = SW.generateLattice(4)
R1 = Rvec(0,0,1)
let 
    P = SW.getPlaquette(R1)
    SW.isInPlaquetteOf.(Ref(R1),P)
    
end
##

let 
    P = SW.getPlaquette(R1)
    ri = [Point2f(SW.getCartesian(R)) for R in P] 
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
    r = getCartesian(R,SW.BasisSW)
    abs(r[1]) == L || abs(r[2]) == L
end

let 
    # sites = SW.generateLattice(10)
    sites = generateLUnitCells(4,SW.Basis)
    filter!(x-> isInLBox(x,2),sites )
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1)
    scatterRvecs!(ax,sites,SW.BasisSW)
    fig
end
## plot 70 combinations

let 
    R0 = Rvec(0,0,1)
    sites = SW.getPlaquette(R0) ∪ [R0]
    
    combs = SW.getAllGS(0.5)
    # sizex = Int(ceil(√(length(combs))))
    # sizey = Int(floor(√(length(combs))))

    sizex = Int(7)
    sizey = Int(10)

    layout = zeros(sizex,sizey)

    fig = Figure(resolution = 200 .* (sizex,sizey))
    allij = Tuple.(CartesianIndices(layout))
    xticks = yticks = [-1,0,1]
    axes = [Axis(fig[fj,fi]; xticklabelsvisible = false, yticklabelsvisible = false,xgridvisible = false,ygridvisible = false,aspect = 1,xticks ,yticks ) for (i,(fi,fj)) in enumerate(allij)]

    for (ax,c) in zip(axes,combs)
        # scatterRvecs!(ax,sites,SW.BasisSW)
        S1,S2,S3,S4,S5,S6,S7,S8 = c
        Mat = [
            S4 S3 S2;
            S5 missing S1;
            S6 S7 S8
        ]
        heatmap!(ax,xticks,yticks,Mat,colorrange = (-0.5,0.5))
    end
    save("combs.pdf",fig,)
    fig
end