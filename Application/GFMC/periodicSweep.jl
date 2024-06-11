cd(@__DIR__)

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.precompile()
import SpiderWebModel as SW
using HDF5
nBra = 8
NSteps = 100_000
equilibration_steps = 30_000
NWalkers = 36*5
α = 0.13
λ = 1
##
for L in [16,18,20,22,24]
    Threads.@threads for i_arg in 1:32
        @info "starting run" L i_arg

        w_avg_estimate = L^2/4 #estimated average weight for each iteration, to reduce floating point errors

        outfile = "/storage/niggeni/Spiderweb/Data/Spin1GFMC_L=$(L)_nBra=$(nBra)_NSteps=$(NSteps)_NW=$(NWalkers)_alpha=$(α)_Lam=$(λ)_$(i_arg).h5"
        mkpath(dirname(outfile))
        @assert !isfile(outfile) "file already exists!"


        parentState = SW.stencilConfig(zeros(L,L),1;boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional())

        ψG = SW.PlaquetteNumberGuidingFunction(α)

        @info "starting run"  
        @time results = SW.startManyWalkerGFMC(parentState,NWalkers,NSteps,nBra,ψG,λ;outfile,equilibration_steps,w_avg_estimate)

    end
end
