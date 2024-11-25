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
import Accessors

include("Types/Utils.jl")
include("Types/PeriodicMatrix.jl")
include("Types/SpinConfig.jl")
include("Types/StencilConfigs.jl")
include("Types/PlaquetteMoves.jl")
include("Types/EDTypes.jl")
include("Types/GFMCTypes.jl")
include("Types/Operators.jl")
include("Types/ObservableTypes.jl")
include("Types/VariationalWFTypes.jl")

include("Constraint.jl")
include("PeriodicTilings.jl")

include("ED/Fluctuations.jl")
include("ED/ConstructHilbertSpace.jl")
include("ED/ED.jl")
include("ED/Observables.jl")

include("GreenFunctionMonteCarlo/GreenFunctionMonteCarlo.jl")
include("GreenFunctionMonteCarlo/Observables.jl")
include("GreenFunctionMonteCarlo/StraightForwardWalking.jl")
include("GreenFunctionMonteCarlo/PlaquetteCorrOperator.jl")
include("GreenFunctionMonteCarlo/BqOperator.jl")
include("GreenFunctionMonteCarlo/RandomFluctuations.jl")

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
