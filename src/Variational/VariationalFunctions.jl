# module VariationalFunctions

include("VariationalFunctions/GuidingFunctions.jl")
include("VariationalFunctions/LocalPlaquetteGuidingFunction.jl")
include("VariationalFunctions/OrderGuidingFunction.jl")
include("VariationalFunctions/PlaquetteCorrelationsGuidingFunctions.jl")
include("VariationalFunctions/PlaquetteNumberGuidingFunctions.jl")
include("VariationalFunctions/RBM.jl")
include("VariationalFunctions/RKFunction.jl")

include("StochasticReconfiguration.jl")

# export AbstractGuidingFunction
# export FullVariationalGuidingFunction, fullVariationalFunction
# export LocalPlaquetteGuidingFunction, localPlaquetteGuidingFunction
# export RKFunction
# export OrderGuidingFunction, orderGuidingFunction

function getOx_k(ψG,Walker::SpiderWebWalker,k)
    return getOx_k(ψG,get_config(Walker),k)
end


function nearbyInt(x1,x2,x_size)
    x_rsize = 1.0 / x_size

    dx = x1 - x2
    dx -= x_size * round(Int,dx * x_rsize)
end

isperiodic(S::Stencils.StencilArray) = Stencils.boundary(S) == Stencils.Wrap()
isperiodic(S::StencilSpinConfig) = isperiodic(parent(S))

function getDistReduction(S,ψG::AbstractGuidingFunction) 
    AllDists = Dict{SVector{2,Int},Int}()
    indicesMapping = collect(eachindex(ψG.params))
    uniqueInds = collect(indicesMapping)
    return (;AllDists,indicesMapping,uniqueInds)
end

function getReducedDist(ri,rj,Lx,Ly) 
    rpr = abs.(SVector(nearbyInt.(ri, rj,(Lx,Ly))))
    return sort(rpr)
end

centralPos(Lx,Ly) = ((Lx+1)//2,(Ly+1)//2)
centralPos(S::AbstractMatrix) = centralPos(size(S)...)
function getCentralPlaquette(S)
    allplaqs = collect(plaquetteIterator(S))
    central = centralPos(S)
    return allplaqs[findmin([norm(central .- r) for r in allplaqs])[2]]
end
function symmetryReducePlaquettes(S,R_ref)
    
    AllDists = Dict{SVector{2,Int},Int}()

    Allplaqs = collect(plaquetteIterator(S))

    indicesMapping = Int[]
    uniqueInds = Int[]
    Lx,Ly = size(S)
    ri = R_ref
    for (j,rj) in enumerate(Allplaqs)
        Rij = getReducedDist(ri,rj,Lx,Ly)
        if Rij ∉ keys(AllDists)
            uniqueInds = push!(uniqueInds,j)
            AllDists[Rij] = length(uniqueInds)
        end
        push!(indicesMapping,AllDists[Rij])
    end
    return (;AllDists,indicesMapping,uniqueInds)
end


function add_reconstructedFullParams!(ψG,indicesMapping,trimmedparams)
    for (i,k) in enumerate(indicesMapping)
        ψG.params[i] += trimmedparams[k]
    end
    return ψG
end

# end