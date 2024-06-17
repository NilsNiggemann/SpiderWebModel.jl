cd(@__DIR__)

using Pkg
Pkg.activate(@__DIR__)
#Pkg.instantiate()
#Pkg.precompile()

import SpiderWebModel as SW
using HDF5
nBra = 10
NSteps = 100_000
equilibration_steps = 100_000
NWalkers = 36*10
α = 0.15
λ = 0
ψG = SW.PlaquetteNumberGuidingFunction(α)

##
for L in [16,18,20]
    for i_arg in 1:8
        @info "starting run" L i_arg

        w_avg_estimate = L^2/6 #estimated average weight for each iteration, to reduce floating point errors

        outfile = "/storage/niggeni/Spiderweb/Data2/L=$(L)/Spin1GFMC_L=$(L)_nBra=$(nBra)_NSteps=$(NSteps)_NW=$(NWalkers)_alpha=$(α)_Lam=$(λ)_$(i_arg).h5"
        mkpath(dirname(outfile))
        @assert !isfile(outfile) "file already exists!"

        DT = SW.DiscreteTimeMethod(λ,nBra,w_avg_estimate)

        parentState = SW.stencilConfig(zeros(L,L),1;boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional())

        @info "starting run" 
        flush(stdout)
        @time results = SW.startManyWalkerGFMC(parentState,DT,NWalkers,NSteps,ψG;equilibration_steps=equilibration_steps,pre_equilibration_steps=equilibration_steps,outfile)

    end
end
