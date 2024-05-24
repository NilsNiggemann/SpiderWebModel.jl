import Pkg
Pkg.activate("Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel

##
S = SW.stencilConfig(zeros(16,16),1;
# boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
ψG = SW.PlaquetteNumberGuidingFunction(0.05)
nThermal = 1000
##

stochReconfRes = SW.repeatStochReconf(S,30000,ψG,40,1e-2;Nwalkers = 6,nbra = 10,error_threshold=1e-2,equilibration_steps=nThermal)
