using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
using SpiderWebModel.Statistics
using SpiderWebModel.StaticArrays
##

function getBBCorrelator(BBVals,BVals,symReduc;index = lastindex(BBVals[1]))

    BBCorrelatorRaw = getBBCorrelator(BBVals,BVals,index)

    BBCorrelator = similar(BBCorrelatorRaw, length(symReduc.indicesMapping))
    for (i,k) in enumerate(symReduc.indicesMapping)
        BBCorrelator[i] = BBCorrelatorRaw[k]
    end
    return BBCorrelator
end

function getBBCorrelator(BBVals,BVals,index::Int)
    nth(x) = x[index]
    BBEnd = nth.(BBVals)
    BEnd =  nth(BVals)
    BBCorrelatorRaw = BBEnd .- BEnd^2
    return BBCorrelatorRaw
end
##
BB_Vec = Dict{Int,Vector{Matrix{Float64}}}()
nBras = Dict{Int,Vector{Int}}()

for i_arg in 1:24
    NumBins = 1
    binnumber,fileindex =  (i_arg-1) %NumBins +1  , (i_arg-1) ÷ NumBins +1

    infile = [joinpath(root,file) for (root,_,files) in walkdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_small2/") for file in files][fileindex]

    nBra = h5read(infile,"nBranch")
    L,_,NWalkers,NSteps = h5open(infile) do f
        size(f["SaveConfigs"])
    end
    res = SW.readResults(infile,NSteps÷ NumBins)[binnumber];

    outfile = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_ForwardWalking_2/L=$(L)/$i_arg/Spin1GFMC_L=$(L)_nBra=$(nBra)_NSteps=$(NSteps)_NW=$(NWalkers)_$(i_arg).h5"
    
    parentState = SW.stencilConfig(zeros(L,L),1;
        boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )

    refPlaq = SW.getCentralPlaquette(parentState)
    symReduc = SW.symmetryReducePlaquettes(parentState,refPlaq)
    GFMCPlaqs = collect(SW.plaquetteIterator(parentState))[symReduc.uniqueInds]

    mProj = 300 ÷ nBra

    res = SW.readResults(infile,NSteps ÷ NumBins)[binnumber];

    resB = h5read(outfile,"WeightsPlaquetteFlipOperator")

    resBB = h5read(outfile,"WeightsSpiderWebModel.BBOperator")

    maxIndex = findfirst(==(0),resBB[1,1,:])

    Gnp = SW.precomputeNormalizedAccWeight(res.TotalWeights,1,mProj)

    BBVals = [SW.get_observables_sfw(Gnp[1:maxIndex,:],resBB[:,j,1:maxIndex]',mean(res.TotalWeights)) for j in eachindex(GFMCPlaqs)]
    # BBVals = stack(stack([SW.get_observables_sfw(Gnp,resBB[:,j,:]',mean(res.TotalWeights)) for j in eachindex(GFMCPlaqs)]))

    BVals = SW.get_observables_sfw(Gnp,resB[:,begin,:]',mean(res.TotalWeights))
    # return BVals, BBVals
    # BBCorrelators = copy(BBVals)
    
    # for i in axes(BBCorrelators,2)
    #     BBCorrelators[:,i] .-= BVals.^2
    # end
    BBCorrelators = stack([
        getBBCorrelator(BBVals,BVals,symReduc,index = i)
        for i in 1:mProj
        ]
    )
    if haskey(BB_Vec,L)
        push!(BB_Vec[L],BBCorrelators)
        push!(nBras[L],nBra)
    else
        BB_Vec[L] = [BBCorrelators]
        nBras[L] = [nBra]
    end
end
##
for (L,BBCorrelators) in BB_Vec
    outfileFinal = "../Data/BBCorr2/BBCorr_L=$(L).h5"
    mkpath(dirname(outfileFinal))
    rm(outfileFinal,force=true)
    h5write(outfileFinal,"BBCorrelator",stack(BBCorrelators,dims=1))
    h5write(outfileFinal,"nBra",nBras[L])
    h5write(outfileFinal,"L",L)
end

