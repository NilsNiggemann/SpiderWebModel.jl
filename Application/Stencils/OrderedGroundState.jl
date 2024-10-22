import Pkg
Pkg.activate(joinpath(@__DIR__,"../"))
cd(joinpath(@__DIR__,"../../"))
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel
include("plottingUtils.jl")

##
L = 30
S = SW.stencilConfig(zeros(L,L),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
S .= SW.periodicStateDenseLoops(L)
##
ψG = SW.CombinedVariationalFunction_1((SW.plaquetteCorrelationFunction(S),SW.localPlaquetteGuidingFunction(S,0.0)))
ψG2 = SW.fullVariationalFunction(S,0.0)

CT = SW.ContinuousTimeMethod(0.1,Hxx=SW.Hxx_RK(0.2))
##
@profview SW.startManyWalkerGFMC(S,CT,12,5000,ψG2;equilibration_steps=0)
@profview SW.startManyWalkerGFMC(S,CT,12,5000,ψG;equilibration_steps=0)
##
@time SW.startManyWalkerGFMC(S,CT,12,1000,ψG.Terms[2];equilibration_steps=0)

##
module MyArrayPartition

    struct MyPartition{T,N,F<:NamedTuple} <: AbstractArray{T,N}
        data::Vector{T}
        partitions::F
    end
    """ create a MyPartition object which wraps subarrays as a vector. Example: MyPartition(x = (5,), y = (3,3)) will wrap a 5 + 3*3 = 14 element vector with the partition arrays x and y begin views"""
    function MyPartition(T=Float64;kwargs...)
        totalLength = 0
        for (k,v) in kwargs
            totalLength += prod(v)
        end
        data = zeros(T,totalLength)
        for (k,v) in kwargs
            sz = size(v)
            data = reshape(data,v)
        end
    end

end