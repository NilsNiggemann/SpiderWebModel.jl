module SpiderWebModel
using StaticArrays, Random, Statistics
using OrderedCollections, Dictionaries, LinearAlgebra
using HDF5, H5Zblosc
using LinearAlgebra, SparseArrays, Arpack
using BitIntegers

import KrylovKit
import ChunkSplitters
import DataStructures
import CircularArrays
import DSP
import JuMP
import Gurobi
import ProgressMeter
import MathOptInterface as MOI
import Stencils


include("SpinConfig.jl")
include("Constraint.jl")
include("PeriodicTilings.jl")

include("Fluctuations/Fluctuations.jl")
include("Fluctuations/ConstructHilbertSpace.jl")
include("Fluctuations/RandomFluctuations.jl")
include("Fluctuations/ED.jl")

include("Plotting/plotSpinConfig.jl")
include("Plotting/PlaquetteFlips.jl")
include("Plotting/Fractons.jl")

include("GroundStateConstruction/TilingPaths.jl")
include("GroundStateConstruction/GroundStateConstruction.jl")
include("GroundStateConstruction/JuMPConstruction.jl")

include("Observables/Observables.jl")

include("explicitStates/PeriodicStates.jl")

include("IO/saveHilbertSpace.jl")

const GRB_ENV_REF = Ref{Gurobi.Env}()

function __init__()
    global GRB_ENV_REF
    GRB_ENV_REF[] = Gurobi.Env()
    return
end

end # module SpiderWebModel
