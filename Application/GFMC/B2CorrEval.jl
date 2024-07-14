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

##
parentState = SW.stencilConfig(zeros(30,30),1;
    boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)

B_Vec = Dict{Float64,Vector{Vector{Float64}}}()

for i_arg in 1:32
    μs = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    μ = μs[(i_arg-1)%length(μs)+1]
    L = 30
    τ = 0.2
    nBra = 1
    
    mProj = round(Int,150 ÷ τ)
    NSteps = 2_000
    equilibration_steps = 3_000
    NWalkers = 128*16
    scatter_fraction = 0.6

    outfile = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK/L=$(L)/$i_arg/Spin1GFMC_L=$(L)_tau=$(τ)_NSteps=$(NSteps)_NW=$(NWalkers)_mu=$(μ)_$(i_arg).h5"
    TotalWeights = h5read(outfile,"TotalWeights")

    Gnp = SW.precomputeNormalizedAccWeight(TotalWeights,1,mProj)
    energies = h5read(outfile,"energies")

    enproj = SW.getEnergies(TotalWeights,energies,1,100)
    
    outfileSFW = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_SFW/$i_arg/Spin1GFMC_L=$(L)_tau=$(τ)_NSteps=$(NSteps)_NW=$(NWalkers)_mu=$(μ)_$(i_arg).h5"

    resB = h5read(outfileSFW,"WeightsRandomPlaquetteFlipOperator")
    BVals = SW.get_observables_sfw(Gnp,resB[:,1,:]',SW.mean(TotalWeights)) ./length(collect(SW.plaquetteIterator(parentState)))

    # outfileSFW2 = "../Data/BBCorrRK/BBCorrRK.h5"
    if haskey(B_Vec,μ)
        push!(B_Vec[μ],BVals)
    else
        B_Vec[μ] = [BVals]
    end
    # h5write(outfileSFW2,"mu=$μ/B",BVals)

end
##
outfileFinal = "../Data/BBCorrRK/BBCorrRK.h5"
mkpath(dirname(outfileFinal))
rm(outfileFinal,force=true)
h5write(outfileFinal,"L",30)
h5write(outfileFinal,"τ",0.2)
h5write(outfileFinal,"nBra",1)

for (μ,BVals) in B_Vec

    h5write(outfileFinal,"BVals/mu=$μ/B",stack(BVals,dims=2))
    h5write(outfileFinal,"BVals/mu=$μ/mu",μ)

end