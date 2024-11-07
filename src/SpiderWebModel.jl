module SpiderWebModel
using StaticArrays, Random, Statistics
using OrderedCollections, Dictionaries, LinearAlgebra
using HDF5, H5Zblosc
using LinearAlgebra, SparseArrays, Arpack
using BitIntegers

import LoopVectorization
import KrylovKit
import ChunkSplitters
import DataStructures
import CircularArrays
import JuMP
import Gurobi
import ProgressMeter
import MathOptInterface as MOI
import Stencils
import StatsBase
import FFTW
import RecursiveArrayTools

include("SpinConfig.jl")
include("Constraint.jl")
include("PeriodicTilings.jl")
include("Observables/Observables.jl")

include("Fluctuations/StencilConfigs.jl")
include("Fluctuations/Fluctuations.jl")
include("Fluctuations/ConstructHilbertSpace.jl")
include("Fluctuations/ED.jl")
include("Fluctuations/GreensFunctionMonteCarlo.jl")
include("Fluctuations/Operators.jl")
include("Fluctuations/StraightForwardWalking.jl")
include("Fluctuations/RandomFluctuations.jl")

include("Variational/VariationalFunctions.jl")
# using .VariationalFunctions

include("Plotting/plotSpinConfig.jl")
include("Plotting/PlaquetteFlips.jl")
include("Plotting/Fractons.jl")

include("GroundStateConstruction/TilingPaths.jl")
include("GroundStateConstruction/GroundStateConstruction.jl")
include("GroundStateConstruction/JuMPConstruction.jl")


include("explicitStates/PeriodicStates.jl")

include("IO/saveHilbertSpace.jl")
include("IO/readGFMCResults.jl")

const GRB_ENV_REF = Ref{Gurobi.Env}()

function initGurobi()
    global GRB_ENV_REF
    GRB_ENV_REF[] = Gurobi.Env()
    return
end

end # module SpiderWebModel
