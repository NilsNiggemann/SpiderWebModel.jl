import SpiderWebModel as SW
using CairoMakie
using SpiderWebModel.StaticArrays
## plot 70 combinations
using SpiderWebModel:plotSpinConfig

AllAllowedConfigs = SW.getAllGS(0.5)
AllConfigs = SW.getAllGS_noMissing(0.5)
##
let 
    sizex = Int(7)
    sizey = Int(10)

    layout = zeros(sizex,sizey)

    fig = Figure(size = 200 .* (sizex,sizey))
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
    function deleter(x,counter)
        counter[] += 1
        # println(counter[])
        if counter[] > tries
            counter[] = 0
            return L
        end
        return default
    end
end
##
function getLevels(path)
    origin = first(path)
    levels = [maximum(abs,x.-origin) for x in path]
end
function getPrevLevels(levels)
    mapnothing(x) = isnothing(x) ? 1 : x

    prevlevel = [mapnothing(findfirst(==(x-2),levels)) for x in levels]

    # prevlevel = [x + (i % x)  for (i,x) in enumerate(prevlevel)]
    # any(isnothing,prevlevel) && error("could not find all previous levels")
end

function getspiralDeleter(spiralpath,default = 3, tries = 5)
    # rsquares = [x[1]^2+x[2]^2 for x in spiralpath]
    levels = getLevels(spiralpath)
    # phis = [atan(x[2],x[1]) for x in spiralpath]
    
    prevlevel = getPrevLevels(levels)
    function deleter(x,counter)
        counter[] += 1
        if counter[] > tries
            counter[] = 0
            return max(1,x - prevlevel[x] - 3)
        end

        # if counter[] > tries ÷ 2
        #     return default + 1
        # end
        return default
    end
end
##

let 
    pt = SW.spiralPath(10)
    
    levels = getLevels(pt)
    lines(Point.(pt),color = levels,colormap = :flag,axis = (;aspect = 1))
    scatter!(Point.(pt),color = levels,colormap = :flag,markersize = 40)

    prevlevel = getPrevLevels(levels)
    # prevlevel = eachindex(levels)
    # return string.(prevlevel)
    text!(Point.(pt),text = string.(prevlevel),color = :white,align = (:center,:center))
    current_figure()
end
##
include("plotStructureFac.jl")
##
function getConfigs(L,numConfigs = 10;kwargs...)
    # paths = [SW.ydirecPathReverse(L),SW.xdirecPathReverse(L),SW.ydirecPath(L),SW.xdirecPath(L),SW.spiralPath(L)]

    Paths = (SW.xdirecPath(L),SW.ydirecPath(L),SW.xdirecPathReverse(L),SW.ydirecPathReverse(L))
    # Paths = (SW.xdirecPathReverse(L),)
    defaultDelete = 1
    tries = 60
    maxiter = 10_500_000
    
    S = fetch.([Threads.@spawn SW.constructConfigPath(SW.DictAlgorithm(),L,L,SW.ALLGS_S12, setup(path),maxiter = 200_000,deleteSteps = getStepDeleter(L+2,1,15),verbose = false;kwargs...) for _ in 1:numConfigs for path in Paths])


    filter!(x -> SW.fulFillsConstraint(x,verbose = false) && !any(isnan,x),S)
    @info "" L defaultDelete tries maxiter length(S) 
    return S
end
##
function getConfigsSpiral(L,numConfigs = 10;kwargs...)
    path = SW.spiralPath(L)
    defaultDelete = 1
    tries = 60
    maxiter = 18_500_000
    deleter = getspiralDeleter(path,defaultDelete,tries)
    setup = SW.setupCalc!(path,L,L,SW.ALLGS_S12)

    S = fetch.([Threads.@spawn SW.constructConfigPath(SW.DictAlgorithm(),L,L,SW.ALLGS_S12, setup,maxiter = maxiter,deleteSteps = deleter,verbose = false;kwargs...) for _ in 1:numConfigs])

    filter!(x -> SW.fulFillsConstraint(x,verbose = false) && !any(isnan,x),S)
    @info "" L defaultDelete tries maxiter length(S) 
    return S
end
##
SW.Random.seed!(345)
# @profview (@time Confs = getConfigs(25,100))
@time Confs = getConfigsSpiral(15,200)
##
using HDF5
let
    L = size(Confs[1],1)
    h5write("ConfsRaw/Confs$(L).h5","Confs",stack(Confs,dims = 3))
end
##
# Confs = [SW.SpinConfig(S,1/2) for S in eachslice(h5read("confs20.h5","Confs"),dims = 3)]
plotStructureFac(Confs,cbar = false)
##
save("Confs/SpiralPathSq_14.pdf",current_figure())
##
function visualizeConstruction(L)
    path = SW.spiralPath(L)
    setup = SW.setupCalc!(path,L,L,SW.ALLGS_S12)
    # setup(path) = SW.setupCalc!(path,L,L,SW.ALLGS_S12)
    S = SW.constructConfigPath(SW.DictAlgorithm(),L,L,SW.ALLGS_S12, setup,maxiter = 200,deleteSteps = getspiralDeleter(path,1,10),verbose = false,plotSteps = true)
    return S
end

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
a,ua = generateAllFluctuations(Confs[1],maxiter=2)
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
function takeFluctuations(Confs,b1 = findOps(Confs[1])[1])
    findNumFlucs = [length(SW.getApplicablePlaquettes_ns(c,b1)) for c in Confs]

    inds = sortperm(findNumFlucs)
    Confs = Confs[inds[end-50:end]]

    flucs = Set(empty(Confs))
    for Conf in Confs
        fl = generateAllFluctuations(Conf,b1,maxiter = 5)[2]
        flucs = flucs ∪ fl
    end

    return flucs
end

function takeFluctuations2(Conf,maxiter = 100,b = findOps(Confs[1]))
    b1,b3 = b[1],b[3]
    Conf2 = copy(Conf)
    for i in 1:maxiter
        plaqs1 = SW.getApplicablePlaquettes(Conf2,b1)
        plaqs3 = SW.getApplicablePlaquettes(Conf2,b3)
        isempty(plaqs1) && isempty(plaqs3) && (@warn "no fluctuations possible"; break)

        if !isempty(plaqs1)
            rn2 = rand(eachindex(plaqs1))
            pij = SW.getPlaquette(Conf2,plaqs1[rn2]...)
            pij .+= b1
        end

        if !isempty(plaqs3)
            rn2 = rand(eachindex(plaqs3))
            pij = SW.getPlaquette(Conf2,plaqs3[rn2]...)
            pij .+= b3
        end
    end
    return Conf2
end
function appendFluctuations(Confs,maxiter,numFlucs = 2,b= findOps(Confs[1]))
    newconfs = copy(Confs)
    for C in Confs
        for i in 1:numFlucs
            push!(newconfs,takeFluctuations2(C,maxiter,b))
            push!(newconfs,takeFluctuations2(C,maxiter,b))
        end
    end
    return newconfs
end
##
function readConfs(File,group="")
    Arr = h5read(File,group*"Confs")
    Confs = [SW.SpinConfig(Matrix(S),1/2) for S in eachslice(Arr,dims = 3)]
end
a = readConfs("ConfsRaw/Confs35.h5")
##
ConfFlucs = appendFluctuations(a,4,2,b)
plotStructureFac(ConfFlucs)