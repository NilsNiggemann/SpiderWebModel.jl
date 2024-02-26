module SpiderWebModel
using StaticArrays, Random, Statistics
using OrderedCollections, Dictionaries, LinearAlgebra
using Test
import ChunkSplitters
import DataStructures

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
end # module SpiderWebModel
