module SpiderWebModel
using StaticArrays, Random, Statistics
using OrderedCollections, Dictionaries, LinearAlgebra
using HDF5, H5Zblosc
import ChunkSplitters
import DataStructures
import CircularArrays
import DSP
import JuMP
import HiGHS

include("SpinConfig.jl")
include("Constraint.jl")
include("PeriodicTilings.jl")

include("Fluctuations/Fluctuations.jl")
include("Fluctuations/Fluctuations_BitVector.jl")
include("Fluctuations/RandomFluctuations.jl")

include("Plotting/plotSpinConfig.jl")
include("Plotting/PlaquetteFlips.jl")
include("Plotting/Fractons.jl")

include("GroundStateConstruction/TilingPaths.jl")
include("GroundStateConstruction/GroundStateConstruction.jl")
include("GroundStateConstruction/BigTilings.jl")

include("Observables/Observables.jl")
include("explicitStates/PeriodicStates.jl")

include("IO/saveHilbertSpace.jl")
end # module SpiderWebModel
