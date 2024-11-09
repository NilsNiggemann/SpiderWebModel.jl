struct OrderGuidingFunction{T,T2<:AbstractVector{T}} <: AbstractGuidingFunction
    α_i::T2
    m_i::T2
    M_i::T2
    N::Int
    Nsites::Int
end
guidingfunc_name(F::OrderGuidingFunction) = "OrderGuidingFunction"
function orderGuidingFunction(State,α::Real=0.1,type=Float32)
    plaqs = collect(plaquetteIterator(State))
    N = length(plaqs)
    Nsites = length(State)
    α_i = zeros(type,N)
    m_i = zeros(type,Nsites)
    M_i = zeros(type,Nsites)

    ψ = OrderGuidingFunction(α_i,m_i,M_i,N,Nsites)
    get_alpha_i(ψ) .= α
    return ψ
end
get_alpha_i(ψ::OrderGuidingFunction) = ψ.α_i
get_m_i(ψ::OrderGuidingFunction) = ψ.m_i
get_M_i(ψ::OrderGuidingFunction) = ψ.M_i
get_params(ψ::OrderGuidingFunction) = RecursiveArrayTools.ArrayPartition(ψ.α_i,ψ.m_i,ψ.M_i)
# get_params(ψ::OrderGuidingFunction) = RecursiveArrayTools.ArrayPartition(ψ.α_i,ψ.m_i)
function (ψ::OrderGuidingFunction)(Walker::SpiderWebWalker)
    α = get_alpha_i(ψ)
    m = get_m_i(ψ)
    M = get_M_i(ψ)
    x = get_config(Walker)
    xlist = reshape(x, length(x))

    n = Walker.n_x
    exp_α = zero(eltype(get_params(ψ)))
    for i in eachindex(n)
        exp_α += α[i] * n[i]
    end

    logψ = zero(eltype(get_params(ψ)))

    for i in eachindex(xlist)
        logψ += m[i] * xlist[i] + 0*M[i] * xlist[i]^2
    end
    return exp(logψ + exp_α)
end

function guidingfuncRatio_exponent(ψG::OrderGuidingFunction,Walker::SpiderWebWalker,move::Tuple,affectedPlaquettes)
    α = get_alpha_i(ψG)
    m = get_m_i(ψG)
    M = get_M_i(ψG)

    Mat = parent(get_config(Walker))

    exp_α = zero(eltype(get_params(ψG)))
    n = Walker.n_x
    n´ = Walker.n_x´
    for i in affectedPlaquettes
        Δn = n´[i] - n[i]
        exp_α += α[i] * Δn
    end

    i,j,opSign = move
    
    sites = safe_parent_indices(Mat, (i, j))
    LI = LinearIndices(Mat)

    x = get_config(Walker)
    # for (ij,s) in zip(sites,P1_STENCIL)
    SpinNormalization = 1. /(2getSpin(x))
    exp_m = zero(exp_α)
    exp_M = zero(exp_α)
    
    @inbounds @simd for idx in eachindex(sites)
        i,j = sites[idx]
        s = P1_STENCIL[idx]*opSign *0.5
        I = LI[i,j]
        # exp_m += m[I]*s + M[I] * (2s*x[I] + s^2)
        exp_m += m[I]*s 
        # exp_M += M[I] * ((x[I] - s*x[I]))^2 - (x[I]^2)
        exp_M += M[I] * (-2x[I]*s +s^2)
    end

    return exp_α  + exp_m *SpinNormalization + exp_M * SpinNormalization^2
end
@inline updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::OrderGuidingFunction,Λ=0) = updateWeightList_plaqs!(Walker,AffectedPlaquetteList,ψG,Λ)


function getOx_k(ψG::OrderGuidingFunction,Walker::SpiderWebWalker,k)
    α = get_alpha_i(ψG)
    par = get_params(ψG)
    type,k = _getParamsTypeAndIndex(par,k)
    
    if type == 1
        return Walker.n_x[k]
    end
    x = get_config(Walker)
    if type == 2
        return x[k]
    end
    if type == 3
        return x[k]^2
    end
    throw(BoundsError(par,k))
end

function _symmetrize_sites(S,unitcell::Tuple{Int,Int})
    MappingDict = Dict{CartesianIndex{2},CartesianIndex{2}}()

    for (i,ri) in enumerate(CartesianIndices(S))
        x,y = Tuple(ri)
        x = x % unitcell[1]
        y = y % unitcell[2]
        symMapped = CartesianIndex(x,y)
        MappingDict[ri] = symMapped
    end
    return MappingDict

end

function symmetrize(S,ψG::OrderGuidingFunction,unitcell::Tuple{Int,Int})
    Map = _symmetrize_sites(S,unitcell)
    AllPlaquettes = collect(plaquetteIterator(S))
    uniqueSitesDict = Dict{CartesianIndex{2},Int}()
    indicesMapping = Int[]
    
    for (i,I) in enumerate(AllPlaquettes)
        I´ = Map[CartesianIndex(I)]
        if I´ ∉ keys(uniqueSitesDict)
            uniqueSitesDict[I´] = length(uniqueSitesDict) + 1
        end
        push!(indicesMapping,uniqueSitesDict[I´])
    end
    current_variational_par = length(uniqueSitesDict)
    empty!(uniqueSitesDict)
    for (i,I) in enumerate(CartesianIndices(S))
        I´ = Map[I]
        if I´ ∉ keys(uniqueSitesDict)
            uniqueSitesDict[I´] = length(uniqueSitesDict) + 1# + current_variational_par
        end
        push!(indicesMapping,uniqueSitesDict[I´] + current_variational_par)
    end
    current_variational_par = length(uniqueSitesDict)
    empty!(uniqueSitesDict)
    for (i,I) in enumerate(CartesianIndices(S))
        I´ = Map[I]
        if I´ ∉ keys(uniqueSitesDict)
            uniqueSitesDict[I´] = length(uniqueSitesDict) + 1# + current_variational_par
        end
        push!(indicesMapping,uniqueSitesDict[I´] + current_variational_par)
    end
    
    uniqueInds = findFirstUniqueIndices(indicesMapping)
    return SymmetryReducedWaveFunction(ψG,indicesMapping,uniqueInds)
    # uniqueInds = Map.uniqueInds

    # AllSites = collect(CartesianIndices(S))
    # for I in CartesianIndices(S)

end

"""given arr, return the indices of all the first unique elements"""
function findFirstUniqueIndices(arr::AbstractArray{T}) where T
    uniqueInds = Int[]
    uniqueVals = Set{T}()
    for (i,val) in enumerate(arr)
        if val ∉ uniqueVals
            push!(uniqueVals,val)
            push!(uniqueInds,i)
        end
    end
    return uniqueInds
end