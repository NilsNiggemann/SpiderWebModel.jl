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

function getOx_k(ψG::Union{FullVariationalGuidingFunction,LocalPlaquetteGuidingFunction},Walker::SpiderWebWalker,k)
    return getOx_k_plaqs(ψG,Walker.n_x,k)
end

function getOx_k(ψG::OrderGuidingFunction,Walker::SpiderWebWalker,k)
    α = get_alpha_i(ψG)

    if k in eachindex(α)
        return Walker.n_x[k]
    end
    # m = get_m_i(ψG)

    Ok = get_config(Walker)[k-lastindex(α)]

    return Ok
end

function getOx_k_plaqs(ψG::FullVariationalGuidingFunction,n::AbstractArray,k)
    par = ψG.params

    α = get_alpha_i(ψG)

    if k in eachindex(α)
        return n[k]
    end
    β = get_beta_ij(ψG)

    (i,j) = Tuple(CartesianIndices(β)[k-lastindex(α)])

    Ok = n[i]*n[j]

    return Ok
end
getOx_k_plaqs(::LocalPlaquetteGuidingFunction,n::AbstractArray,k) = n[k]

function getOx_k(ψG::RBM,x::AbstractMatrix,k)
    par = ψG.params

    paramtype = getParameterType(ψG,k)

    if paramtype == 1 
        return x[k]
    end
    bj = get_b_j(ψG)
    Wij = get_W_ij(ψG)
    
    if paramtype == 2
        j = k - ψG.N
        θj = get_theta_j(x,j,bj,Wij)
        return tanh(θj)
    elseif paramtype == 3
        i,j = Tuple(CartesianIndices(Wij)[k-ψG.N-length(bj)])
        θj = get_theta_j(x,j,bj,Wij)
        return x[i] * tanh(θj)
    end
    return Ok
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

function getDistReduction(S,ψG::FullVariationalGuidingFunction)
    
    AllDists = Dict{SVector{2,Int},Int}()
    if !isperiodic(S)
        indicesMapping = collect(eachindex(ψG.params))
        uniqueInds = collect(indicesMapping)
        return (;AllDists,indicesMapping,uniqueInds)
    end

    α = get_alpha_i(ψG)
    Allplaqs = collect(plaquetteIterator(S))

    # AllDists = Dict{Tuple{Int,Int},Int}()
    
    betaIndex = lastindex(α)
    indicesMapping = ones(Int,betaIndex)
    uniqueInds = [1]
    # indicesMapping = collect(eachindex(α))
    # uniqueInds = collect(eachindex(α))
    Lx,Ly = size(S)

    for (i,ri) in enumerate(Allplaqs)
        for (j,rj) in enumerate(Allplaqs)
            Rij = getReducedDist(ri,rj,Lx,Ly)
            # Rij = SVector(0,0)
            # Rij = (i,j)
            betaIndex += 1
            if Rij ∉ keys(AllDists)
                uniqueInds = push!(uniqueInds,betaIndex)
                AllDists[Rij] = length(uniqueInds)
            end
            push!(indicesMapping,AllDists[Rij])
        end
    end

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

function getDistReduction(S,ψG::LocalPlaquetteGuidingFunction)
    AllDists = Dict{SVector{2,Rational{Int}},Int}()
    if isperiodic(S)
        indicesMapping = ones(Int,length(ψG.params))
        uniqueInds = [1]
        return AllDists,indicesMapping,uniqueInds
    end
    
    α = get_alpha_i(ψG)
    Allplaqs = collect(plaquetteIterator(S))

    indicesMapping = Int[]
    uniqueInds = Int[]
    LxLy = size(S)
    r_Central = centralPos(S)
    for (i,ri) in enumerate(Allplaqs)
        x,y = ri .- r_Central
        if y < -x
            x,y = -y,-x
        end
        if y>x
            x,y = y,x
        end

        symMapped = SVector(x,y)
        if symMapped ∉ keys(AllDists)
            uniqueInds = push!(uniqueInds,i)
            AllDists[symMapped] = length(uniqueInds)
        end
        push!(indicesMapping,AllDists[symMapped])
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