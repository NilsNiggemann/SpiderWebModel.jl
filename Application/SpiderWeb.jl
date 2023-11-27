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
# function findOps(Conf)
#     r = -1.:1.
#     ops(a,b,c,d,e,f,g,h,i) = SW.SpinConfig(SMatrix{3,3,Float64,9}(a,b,c,d,e,f,g,h,i),1/2)
#     it = Iterators.product(r,r,r,r,r,r,r,r,r)

#     Allops = [ops(a,b,c,d,e,f,g,h,i) for (a,b,c,d,e,f,g,h,i) in it]
    
#     # return reshape(Allops,length(Allops))
#     filter!(x ->SW.CanApplyAnywhere(Conf,x),reshape(Allops,length(Allops)))
#     return Allops
# end

##
function findOps(Conf)
    r = -1.:1.
    Allops = (SW.SpinConfig(SMatrix{3,3}(a,b,c,d,e,f,g,h,i),1/2) for a in r for b in r for c in r for d in r for e in r for f in r for g in r for h in r for i in r)
    # filter!(x ->SW.CanApplyAnywhere(Conf,x),Allops)
    Allops = [x for x in Allops if SW.CanApplyAnywhere(Conf,x)]
    return Allops
end
##
function plotApplPlaquettes(State)
    b = findOps(State)
    plaqs = SW.getApplicablePlaquettes_ns(State,b[1])
    plaqs2 = SW.getApplicablePlaquettes_ns(State,b[3])
    fig = plotSpinConfig(State,markersize = 23)
    ax = current_axis()
    points = Point2f.(plaqs)
    points2 = Point2f.(plaqs2)
    scatter!(ax,points,markersize = 13,color = :red)
    scatter!(ax,points2,markersize = 13,color = :lime)
    fig

end
##
let 

    StairCase = getStairCase(15)
    plotApplPlaquettes(StairCase)
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
generateRandomGroundState(5)
##
function getStepDeleter(L,default = 5,tries = 5)
    counter = 0
    function deleter(x)
        counter += 1
        if counter > tries
            return L
        end
        return default
    end
end
##
function getConfigs(L,numConfigs = 10;kwargs...)
    # paths = [SW.ydirecPathReverse(L),SW.xdirecPathReverse(L),SW.ydirecPath(L),SW.xdirecPath(L),SW.spiralPath(L)]

    Paths = (SW.xdirecPath(L),SW.ydirecPath(L),SW.xdirecPathReverse(L),SW.ydirecPathReverse(L))
    # Paths = (SW.xdirecPathReverse(L),SW.xdirecPathReverse(L))
    # Paths = (SW.spiralPath(L),)

    S = fetch.([Threads.@spawn SW.constructConfigPath(SW.DictAlgorithm(),L,L,SW.ALLGS_S12, path,maxiter = 400_000,deleteSteps = getStepDeleter(L+2,3,100),verbose = false;kwargs...) for _ in 1:numConfigs for path in Paths])

    filter!(x -> SW.fulFillsConstraint(x,verbose = false) && !any(isnan,x),S)
    return S
end
##
SW.Random.seed!(1234)
# @profview (@time Confs = getConfigs(18,20))
@time Confs = getConfigs(18,1000)
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
using LatticeFFTs
using MakieHelpers
using Statistics
using LatticeFFTs.Interpolations
##
function getStructureFac(Confs)
    plan = getLatticeFFTPlan(Confs[1].Mat,0)

    Sq = [getInterpolatedFFT(c.Mat,0,plan;Interpolation = BSpline(Constant())) for c in Confs]
end
##
function plotStructureFac(Confs;kwargs...)
    Confs2 = Confs
    # Confs2 = randRots(Confs)

    Sq = getStructureFac(Confs2)
    SSq(kx,ky) = mean(real(s(kx,ky)*s(-kx,-ky)) for s in Sq)
    # SSq(kx,ky) = mean(real(s(kx,ky)) for s in Sq)

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

    kx = LinRange(-0,2pi,200)
    ky = LinRange(-0,2pi,200)
    chi = fetch.([Threads.@spawn SSq(kx,ky) for kx in kx, ky in ky])
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1,xticks = PiTicks(),yticks = PiTicks())
    heatmap!(ax,kx,ky,chi;kwargs...)
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
a = generateFluctuations(getStairCase(16),500)
##
plotStructureFac(collect(a))
##
# plotStructureFac([Confs[10]])
using OrderedCollections
function generateAllFluctuations(StartConfig,operator = findOps(StartConfig)[1];maxiter = 500)

    Configs = [[copy(StartConfig)]]
    uniqeConfs = Set([copy(StartConfig)])
    for iter in 1:maxiter-1
        Confs = Configs[iter]
        newconf = empty(Confs)
        for Conf in Confs
            plaqs = SW.getApplicablePlaquettes(Conf,operator)
            # isempty(plaqs) && (@warn "no fluctuations possible"; break)
            
            for p in plaqs
                Conf2 = copy(Conf)
                pij = SW.getPlaquette(Conf2,p...)
                pij .+= operator
                if Conf2 ∉ uniqeConfs
                    push!(newconf,Conf2)
                    push!(uniqeConfs,Conf2)
                end
            end
            if isempty(newconf) 
                @info "" length(uniqeConfs)
                return Configs,uniqeConfs
            end
            # @assert i1 == i2 "i1 = $i1, i2 = $i2"
        end
        # @info "" length(newconf) length(uniqeConfs)
        push!(Configs,newconf)
    end
    # filter!(x -> SW.fulFillsConstraint(x,verbose = false),Configs)
    # @assert all(SW.fulFillsConstraint.(Configs,verbose = true))
    @warn "max iterations reached" length(uniqeConfs)
    return Configs,uniqeConfs
end
##
function generateFluctuations(StartConfig,b1,b3,maxiter = 500)

    Configs = Set([copy(StartConfig),])

    for i in 1:maxiter

        Conf = copy(rand(Configs))
        plaqs1 = SW.getApplicablePlaquettes_ns(Conf,b1)
        plaqs3 = SW.getApplicablePlaquettes_ns(Conf,b3)
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
##
a,ua = generateAllFluctuations(getStairCase(16),maxiter=3)
##
b = findOps(getStairCase(10))
##
# generateFluctuations(Confs[1],b[1],b[3],1)
ConfFluc = [collect(generateFluctuations(c,b[1],b[3],100)) for c in ua]
ConfFluc = append!(ConfFluc...)
##
a = SW.constructConfigPath(6,6,SW.ALLGS_S12,SW.ydirecPathReverse(6))
##
display.(plotSpinConfig.(a))
##
