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
function plotSpinConfig!(ax,S::SW.SpinConfig;kwargs...)
    vals = filter!(!ismissing,unique(S.Mat))

    Amin = min(minimum(vals),-S.S)
    Amax = max(maximum(vals),S.S)
    us = Amin:Amax
    hm = heatmap!(ax,Array(S.Mat),colorrange = (Amin,Amax),colormap = cgrad(:linear_bgy_10_95_c74_n256, length(us), categorical = true))
    return hm
end

function plotSpinConfig(S;kwargs...) 
    fig = Figure()
    ax = Axis(fig[1,1],aspect = DataAspect(),backgroundcolor = :grey)
    hm = plotSpinConfig!(ax,S;kwargs...)
    us = filter!(!ismissing,unique(S))
    Colorbar(fig[1,2],hm,ticks = us)
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

##
function findOps(Conf)
    r = -1.:1.
    Allops = [SW.SpinConfig(SMatrix{3,3}(a,b,c,d,e,f,g,h,i),1/2) for a in r for b in r for c in r for d in r for e in r for f in r for g in r for h in r for i in r]
    filter!(x ->SW.CanApplyAnywhere(Conf,x),Allops)
    return Allops
end
##

let 

    StairCase = getStairCase(15)
    b = findOps(StairCase)
    plaqs = SW.getApplicablePlaquettes(StairCase,b[1])
    fig = plotSpinConfig(StairCase,markersize = 23)
    ax = current_axis()
    points = Point2f.(plaqs)
    scatter!(ax,points,markersize = 13,color = :red)
    fig
end
##
function generateRandomGroundState(L,S=1/2;maxiter = 1_000_000)
    AllowedPlaquettes = SW.getAllGS(S)

    PlaqetteMatrix = zeros(Int64,L,L)

    randBuffer = rand(1:70,maxiter)
    it = 0
    for i in 1:L
        for j in 1:L
            PlaqetteMatrix[i,j] = randBuffer[(i-1)*L+j]
            it +=1
        end
    end
    return nothing
    error("Could not find a ground state")
end
generateRandomGroundState(10)
##
S1 = SW.constructSpinConfigFromPlaquettes(4,4,AllAllowedConfigs)
##
using FRGLatticeEvaluation

##
using FFTW

"""given a matrix S_ij return the correlator C_ij = <S_ij S_00>"""
function getSij!(Sij,Conf::SW.SpinConfig,i,j)
    S = Conf.Mat
    si = S[i,j]
    for k in axes(S,1)
        for l in axes(S,2)
            Sij[k,l] = si*S[k,l]
        end
    end
    return Sij
end
getSij(Conf::SW.SpinConfig,i,j) = getSij!(similar(Conf.Mat),Conf,i,j)
function getfft(Conf::SW.SpinConfig)
    Sij = similar(Conf.Mat)
    Sq = zeros(ComplexF64,size(Conf.Mat))
    for i in axes(Conf.Mat,1)
        for j in axes(Conf.Mat,2)
            i == j || continue
            getSij!(Sij,Conf,i,j)
            FT = fft(Sij)
            Sq .+= FT
        end
    end
    return real.(Sq)
end

function plotfft(Conf::SW.SpinConfig)
    FT = getfft(Conf)

    k = fftfreq(size(Conf.Mat,1))
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1)
    heatmap!(ax,k,k,real(FT))
    fig
end
plotfft(AFM)

##
function test(LPx,LPy,PlaquetteList;randBuffSize=1_000_000,maxNumTries = 10_000,deleteRows = 4, randBuffer = rand(1:70,randBuffSize)
)
    Lx = 2*LPx+1
    Ly = 2*LPy+1

    Mat = Matrix{Union{Missing,Float64}}(undef,Lx,Ly)
    fill!(Mat,missing)

    El = PlaquetteList[begin]
    P = SW.SpinConfig(Mat,El.S)

    j = 2
    while j <= Ly-1
        i = 2
        while i <= Lx-1 
        
            Pij = SW.getPlaquette(P,i,j)
            
            tileNums = SW.getPossibleTiles(Pij,PlaquetteList)
            
            if isempty(tileNums)
                jmin = max(j-deleteRows,2)
                
                Mat[:,jmin:end] .= missing
                j = jmin
                i = 2
                display(plotSpinConfig(P))
                continue
            end
            T = PlaquetteList[SW.getRanTileNum!(randBuffer,tileNums)]
            Pij .= T
            display(plotSpinConfig(P))

            i += 2
        end
        j += 2

    end
    return P

end
##
randBuffer = rand(1:70,1_000_000)
##
test(8,5,AllAllowedConfigs;randBuffer=copy(randBuffer))