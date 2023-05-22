import SpiderWebModel as SW
using CairoMakie
using SpiderWebModel.StaticArrays
## plot 70 combinations

AllAllowedConfigs = SW.getAllGS(0.5)
AllConfigs = SW.getAllGS_noMissing(0.5)
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
    vals = filter!(x-> !ismissing(x) && !isnan(x),unique(S.Mat))
    isempty(vals) && (vals = [-1,1])
    Amin = min(minimum(vals),-S.S)
    Amax = max(maximum(vals),S.S)
    us = Amin:Amax
    hm = heatmap!(ax,Array(S.Mat),colorrange = (Amin,Amax),colormap = cgrad(:linear_bgy_10_95_c74_n256, length(us), categorical = true))

    translate!(hm, 0, 0, -100)
    return hm
end

function plotSpinConfig(S;kwargs...) 
    fig = Figure()
    ax = Axis(fig[1,1],
    aspect = DataAspect(),
    backgroundcolor = :grey,
    # xminorgridwidth = 2,
    # yminorgridwidth = 2,
    xminorgridcolor = :black,
    yminorgridcolor = :black,
    xminorgridvisible = true,
    yminorgridvisible = true,
    xgridvisible = false,
    ygridvisible = false,
    # xticks = (axes(S,1),string.(axes(S,1))) ,
    # yticks = (axes(S,2),string.(axes(S,2))) ,
    xminorticks = 0.5 .+ axes(S,1),
    yminorticks = 0.5 .+ axes(S,2),
    )
    hm = plotSpinConfig!(ax,S;kwargs...)
    us = filter!(x -> !ismissing(x) && !isnan(x) ,unique(S))
    isempty(us) && (us = [-1,1])
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
    S = SW.SpinConfig([0.5*(-1)^(i+j) for i in 1:100,j in 1:100],1/2)
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

function getConfigs(Lx,Ly,numConfigs = 10)
    S = fetch.([Threads.@spawn SW.constructSpinConfigFromPlaquettes(Lx,Ly,SW.ALLGS_S12_NOMISSING;maxNumTries = 4_000_000) for _ in 1:numConfigs])

    filter!(x -> SW.fulFillsConstraint(x,verbose = false) && !any(isnan,x),S)
    return S
end
Confs = getConfigs(22,22,150)
##
"""given a matrix, rotate it by 90 degrees"""
function rotate(Mat)
    Mat2 = zeros(size(Mat))
    for i in axes(Mat,1)
        row = Mat'[:,i]
        Mat2[:,i] .= reverse(row)
    end
    return Mat2

end
function rotate(Mat,n)
    if n == 0
        return Mat
    elseif n == 1
        return rotate(Mat)
    elseif n == 2
        return Mat |> rotate |> rotate
    elseif n == 3
        return Mat |> rotate |> rotate |> rotate
    end
    error("n must be 0,1,2,3")
end
rotate(S::SW.SpinConfig,n) = SW.SpinConfig(rotate(S.Mat,n),S.S)

function AllRots(Confs)
    Confs2 = [rotate(c,n) for c in Confs for n in 0:3]
    return Confs2
end

function randRots(Confs)
    Confs2 = [rotate(c,rand(0:3)) for c in Confs]
    return Confs2
end

##
using FRGLatticeEvaluation
using MakieHelpers

function getStructureFac(Confs)
    
    @time Sij = SW.getSij_fast(Confs)
    @time Rij = SW.getRij_vec(Confs[1])
    chik = FourierStruct(Sij,Rij,length(Sij))
end
##
function plotStructureFac(Confs)
    Confs2 = Confs
    # Confs2 = randRots(Confs)

    chik = getStructureFac(Confs2)

    # i = size(Confs2[1],1) ÷ 2
    # j = size(Confs2[1],2) ÷ 2

    # ij = LinearIndices(Confs2[1])[i,j]

    # Sij = SW.getSij(Confs2,ij)[:]
    # Sij2 = SW.getSij(Confs2,ij+1)[:]
    # append!(Sij,Sij2)
    # Rij = SW.getRij_vec(Confs2[1],ij)[:]
    # Rij2 = SW.getRij_vec(Confs2[1],ij+1)[:]
    # append!(Rij,Rij2)
    # chik = FourierStruct(Sij,Rij,1)

    kx = LinRange(0,2pi,80)
    ky = LinRange(0,2pi,80)
    chi = fetch.([Threads.@spawn chik(kx,ky) for kx in kx, ky in ky])
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1,xticks = PiTicks(),yticks = PiTicks())
    heatmap!(ax,kx,ky,chi)
    fig
end
plotStructureFac(Confs)
##

function generateFluctuations(StartConfig,maxiter = 500)

    Configs = Set([copy(StartConfig),])
    b1 = findOps(StartConfig)[1]
    b3 = findOps(StartConfig)[3]
    for i in 1:maxiter

        Conf = copy(rand(Configs))

        plaqs1 = SW.getApplicablePlaquettes(Conf,b1)
        plaqs3 = SW.getApplicablePlaquettes(Conf,b3)
        isempty(plaqs1) && isempty(plaqs3) && (@warn "no fluctuations possible"; break)

        if !isempty(plaqs1)
            rn2 = rand(eachindex(plaqs1))
            pij = SW.getPlaquette(Conf,plaqs1[rn2]...)
            pij .+= b1
            push!(Configs,Conf)
        end
        if !isempty(plaqs3)
            rn2 = rand(eachindex(plaqs3))
            pij = SW.getPlaquette(Conf,plaqs3[rn2]...)
            pij .+= b3
            push!(Configs,Conf)
        end
    end
    filter!(x -> SW.fulFillsConstraint(x,verbose = false),Configs)
    return Configs
end
a = generateFluctuations(Confs[10],500)
##
@profview plotStructureFac(collect(a))
# plotStructureFac([Confs[10]])

# @profview a = [test(18,18,SW.ALLGS_S12_NOMISSING;maxNumTries = 1_0000_000,triesDeleteRow = 18,deleteRows = 2) for i in 1:100]
@time a = [SW.constructSpinConfigFromPlaquettes(18,18,SW.ALLGS_S12_NOMISSING;maxNumTries = 1_0000_000,triesDeleteRow = 18,deleteRows = 2) for i in 1:100]
filter!(x -> !any(isnan.(x.Mat)),a)
##
@time a = SW.constructSpinConfigFromPlaquettes_old(14,14,SW.ALLGS_S12_NOMISSING;maxNumTries = 1_000_000,triesDeleteRow = 20)


