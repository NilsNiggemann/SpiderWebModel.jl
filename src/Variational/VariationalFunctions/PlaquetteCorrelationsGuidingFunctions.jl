struct FullVariationalGuidingFunction{A<:AbstractArray} <: AbstractGuidingFunction
    params::A
end
guidingfunc_name(F::FullVariationalGuidingFunction) = "FullVariationalGuidingFunction"

function fullVariationalFunction(State,α::Real=0.1)
    plaqs = collect(plaquetteIterator(State))
    N = length(plaqs)
    params = zeros(N,N+1)
    ψ = FullVariationalGuidingFunction(params)
    get_alpha_i(ψ) .= α
    return ψ
end

function get_alpha_i(ψG::FullVariationalGuidingFunction)
    return @view ψG.params[:,begin]
end

function get_beta_ij(ψG::FullVariationalGuidingFunction)
    return @view ψG.params[:,begin+1:end]
end

function (ψG::FullVariationalGuidingFunction)(N□::AbstractArray) 
    return exp(guidingfunc_exponent(ψG,N□))
end

function guidingfunc_exponent(ψG::FullVariationalGuidingFunction,N□::AbstractArray) 
    α = get_alpha_i(ψG)
    β = get_beta_ij(ψG)
    exponent = α' * N□ + dot(N□,β,N□)
    return exponent
end

@inline guidingfuncRatio(ψG::FullVariationalGuidingFunction,Walker::SpiderWebWalker,move) = exp(guidingfuncRatio_exponent(ψG,n,n´,move))

function guidingfuncRatio_exponent(ψG::FullVariationalGuidingFunction,Walker::SpiderWebWalker,move,affectedPlaquettes)
    α = get_alpha_i(ψG)
    β = get_beta_ij(ψG)
    _guidingfuncRatio_exponent(α,β,Walker,affectedPlaquettes)
end

@inline updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::FullVariationalGuidingFunction,Λ=0) = updateWeightList_plaqs!(Walker,AffectedPlaquetteList,ψG,Λ)

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


function getOx_k(ψG::FullVariationalGuidingFunction,Walker::SpiderWebWalker,k)
    return getOx_k_plaqs(ψG,Walker.n_x,k)
end

function getNonSymmetric(ψG::AbstractGuidingFunction)
    indicesMapping = collect(eachindex(get_params(ψG)))
    uniqueInds = collect(indicesMapping)
    return SymmetryReducedWaveFunction(ψG,indicesMapping,uniqueInds)
end

function getDistReduction(S,ψG::FullVariationalGuidingFunction)
    
    AllDists = Dict{SVector{2,Int},Int}()
    if !isperiodic(S)
        indicesMapping = collect(eachindex(ψG.params))
        uniqueInds = collect(indicesMapping)
        return SymmetryReducedWaveFunction(ψG,indicesMapping,uniqueInds)
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

    return SymmetryReducedWaveFunction(ψG,indicesMapping,uniqueInds)

end
