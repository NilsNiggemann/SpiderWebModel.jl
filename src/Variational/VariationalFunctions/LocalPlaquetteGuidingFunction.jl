struct LocalPlaquetteGuidingFunction{A<:AbstractVector} <: AbstractGuidingFunction
    params::A
end
guidingfunc_name(F::LocalPlaquetteGuidingFunction) = "LocalPlaquetteGuidingFunction"

function localPlaquetteGuidingFunction(State,α::Real=0.1,type=Float32)
    plaqs = collect(plaquetteIterator(State))
    N = length(plaqs)
    params = zeros(type,N)
    ψ = LocalPlaquetteGuidingFunction(params)

    get_alpha_i(ψ) .= α
    return ψ
end
function compute_GWF_buffer!(Buffer,ψG::LocalPlaquetteGuidingFunction,Walker::SpiderWebWalker) 
    getNPlaq!(Walker)
    return Buffer
end

# @inline updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::LocalPlaquetteGuidingFunction,Λ=0) = updateWeightList_plaqs!(Walker,AffectedPlaquetteList,ψG,Λ)

function get_alpha_i(ψG::LocalPlaquetteGuidingFunction)
    return ψG.params
end

function (ψG::LocalPlaquetteGuidingFunction)(N□::AbstractVector)
    α = get_alpha_i(ψG)
    exponent = α' * N□
    return exp(exponent)
end
(ψG::LocalPlaquetteGuidingFunction)(W::SpiderWebWalker) = ψG(W.n_x)
function (ψG::LocalPlaquetteGuidingFunction)(x::StencilSpinConfig) 
    α = get_alpha_i(ψG)
    exponent = zero(eltype(α))
    for (i,I) in enumerate(plaquetteIterator(x))
        applPlus, applMinus = P_applicable(x, I)
        n_i = applMinus + applPlus
        exponent += α[i] * n_i
    end
    return exp(exponent)
end

function guidingfuncRatio_log(ψG::LocalPlaquetteGuidingFunction,Walker::SpiderWebWalker,move,AffectedPlaquetteList)
    α = get_alpha_i(ψG)


    i,j,opNum = move

    affectedPlaquettes = AffectedPlaquetteList[i,j]

    Config = get_config(Walker)
    applyPlaquette!(Config,i,j,opNum)
    getNPlaqfilled!(Walker,affectedPlaquettes)
    applyPlaquette!(Config,i,j,-opNum)

    n = Walker.n_x
    n´ = Walker.n_x´

    exponent = zero(eltype(n))

    @inbounds @simd for ind in eachindex(affectedPlaquettes)
        i = affectedPlaquettes[ind]
        Δn = n´[i] - n[i]
        exponent += α[i]*Δn
    end

    return exponent
end
getOx_k_plaqs(::LocalPlaquetteGuidingFunction,n::AbstractArray,k) = n[k]

function getOx_k(ψG::LocalPlaquetteGuidingFunction,Walker::SpiderWebWalker,k)
    return getOx_k_plaqs(ψG,Walker.n_x,k)
end

function getDistReduction(S,ψG::LocalPlaquetteGuidingFunction)
    AllDists = Dict{SVector{2,Rational{Int}},Int}()
    if isperiodic(S)
        indicesMapping = ones(Int,length(ψG.params))
        uniqueInds = [1]
        return SymmetryReducedWaveFunction(ψG,indicesMapping,uniqueInds)
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
    
    return SymmetryReducedWaveFunction(ψG,indicesMapping,uniqueInds)

end
