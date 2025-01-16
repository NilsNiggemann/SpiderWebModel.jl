"""Abstract supertype for a walker in a GFMC simulation. A walker needs to have a method `get_config` which returns the configuration (for example the spin configuration), and a method `get_weights` which returns a vector of weights for each possible move."""
abstract type AbstractWalker end
struct SpiderWebWalker{C} <: AbstractWalker
    Config::C
    moves::Vector{Tuple{Int,Int,Int}}
    weights::Vector{Float64}
    Plaquette_positions::Vector{Tuple{Int,Int}}
    n_x::Vector{Float64}
    n_x´::Vector{Float64}
end
function spiderWebWalker(Config,Plaquette_positions)
    moves = Vector{Tuple{Int,Int,Int}}()
    weights = Vector{Float64}()
    # Plaquette_positions = collect(plaquetteIterator(Config))
    n_x = zeros(Float64,length(Plaquette_positions))
    n_x´ = zeros(Float64,length(Plaquette_positions))
    return SpiderWebWalker(copy(Config),moves,weights,Plaquette_positions,n_x,n_x´)
end

get_config(Walker::SpiderWebWalker) = Walker.Config
get_weights(Walker::SpiderWebWalker) = Walker.weights
plaquetteIterator(Walker::SpiderWebWalker) = Walker.Plaquette_positions

getOperatorRep(i,j,opNum) = i,j,opNum
# getOperatorRep(i,j,opNum) = CartesianIndex(i,j,opNum)
# getOperatorNumber(i,j,opNum,L) = LinearIndices((L,L,2))[i,j,opNum]

abstract type AbstractGFMCMethod end

struct DiscreteTimeMethod <: AbstractGFMCMethod 
    Λ::Float64
    nBranch::Int
    w_avg_estimate::Float64
end

DiscreteTimeMethod(;Λ=1.,nBranch,w_avg_estimate=1.) = DiscreteTimeMethod(Λ,nBranch,w_avg_estimate)

struct ContinuousTimeMethod{F2} <: AbstractGFMCMethod 
    τ::Float64
    nBranch::Int
    w_avg_estimate::Float64
    Hxx::F2
end
ContinuousTimeMethod(τ,nBranch::Integer,w_avg_estimate=1.,Hxx=Hxx_zero()) = ContinuousTimeMethod(float(τ),nBranch,float(w_avg_estimate),Hxx)
ContinuousTimeMethod(τ;nBranch=1,w_avg_estimate=1.,Hxx=Hxx_zero()) = ContinuousTimeMethod(τ,nBranch,float(w_avg_estimate),Hxx)

abstract type AbstractGFMCProblem end
struct SpiderwebGFMCProblem{MethodType<:AbstractGFMCMethod,T<:AbstractFloat,BufferType,C,F,W,O} <: AbstractGFMCProblem
    method::MethodType
    InitialState::C
    ψG::F
    Walkers::Vector{W}
    weights::Vector{T}
    Guiding_function_buffer::BufferType
    reconfiguration_buffer::Vector{T}
    Observables::O
end

function get_Guiding_function_buffer(problem::SpiderwebGFMCProblem)
    return problem.Guiding_function_buffer
end

abstract type AbstractGFMCObservables end
struct GFMCObservables{DT<:AbstractFloat,T,T2} <: AbstractGFMCObservables
    energies::Vector{DT}
    SaveConfigs::T
    TotalWeights::Vector{DT}
    reconfigurationTable::Matrix{Int}
    outfile::T2
end
get_reconfigurationTable(Observables::AbstractGFMCObservables) = Observables.reconfigurationTable
get_energies(Observables::AbstractGFMCObservables) = Observables.energies
get_TotalWeights(Observables::AbstractGFMCObservables) = Observables.TotalWeights
get_outfile(Observables::AbstractGFMCObservables) = Observables.outfile

abstract type AbstractOperator end
operatorname(X::T) where T <: AbstractOperator = string(T)