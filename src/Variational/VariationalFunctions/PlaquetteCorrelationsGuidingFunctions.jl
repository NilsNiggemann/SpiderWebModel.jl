struct FullVariationalGuidingFunction{A<:AbstractArray} <: AbstractGuidingFunction
    params::A
end
guidingfunc_name(F::FullVariationalGuidingFunction) = "FullVariationalGuidingFunction"

function fullVariationalFunction(State,α::Real=0.1,type=Float32)
    plaqs = collect(plaquetteIterator(State))
    N = length(plaqs)
    params = zeros(type,N,N+1)
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

function (ψG::FullVariationalGuidingFunction)(N□::AbstractVector) 
    return exp(guidingfunc_exponent(ψG,N□))
end
function (ψG::FullVariationalGuidingFunction)(S::StencilSpinConfig) 
    N□ = getNPlaq(S)
    return exp(guidingfunc_exponent(ψG,N□))
end

function guidingfunc_exponent(ψG::FullVariationalGuidingFunction,N□::AbstractVector) 
    α = get_alpha_i(ψG)
    β = get_beta_ij(ψG)
    exponent = α' * N□ + dot(N□,β,N□)
    return exponent
end

function fill_GWF_buffer!(Buffer,ψG::FullVariationalGuidingFunction,Walker::SpiderWebWalker) 
    getNPlaq!(Walker)
    return Buffer
end
# function fill_GWF_buffer!(Buffer,ψG::FullVariationalGuidingFunction,Walker::SpiderWebWalker) 
#     ψG(Walker)
# end

function guidingfuncRatio(ψG::FullVariationalGuidingFunction,Walker::SpiderWebWalker,move,AffectedPlaquetteList::AbstractMatrix)
    α = get_alpha_i(ψG)
    β = get_beta_ij(ψG)

    i,j,opNum = move
    affectedPlaquettes = AffectedPlaquetteList[i,j]

    Config = get_config(Walker)

    applyPlaquette!(Config, i, j, opNum)
    getNPlaqfilled!(Walker,affectedPlaquettes)
    applyPlaquette!(Config, i, j, -opNum)

    n = Walker.n_x
    n´ = Walker.n_x´

    exponent = zero(eltype(α))
    exponent = dot(n,α) + dot(n,β,n) - dot(n´,α) - dot(n´,β,n´)
    return exp(-exponent)

    for i in affectedPlaquettes
        Δn = n´[i] - n[i]
        iszero(Δn) && continue
        
        exp_i = α[i]
        LoopVectorization.@turbo for j in eachindex(n,n´)
            exp_i += β[j,i]*(n´[j] + n[j])
        end
        exponent += exp_i*Δn
    end
    return exp(exponent)
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


function getOx_k(ψG::FullVariationalGuidingFunction,Walker::SpiderWebWalker,k)
    getNPlaq!(Walker)
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
