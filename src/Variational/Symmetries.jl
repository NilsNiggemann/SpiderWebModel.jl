
"""Contains guiding wave function and information about symmetry equivalent parameters"""
struct SymmetryReducedWaveFunction{F<:AbstractGuidingFunction} <: AbstractGuidingFunction
    psi::F
    indicesMapping::Vector{Int}
    uniqueInds::Vector{Int}
end

@inline get_params(P::SymmetryReducedWaveFunction) = get_params(P.psi)
(psi::SymmetryReducedWaveFunction)(W::SpiderWebWalker) = psi.psi(W)

@inline guidingfunc_name(F::SymmetryReducedWaveFunction) = guidingfunc_name(F.psi)
@inline guidingfuncRatio(F::SymmetryReducedWaveFunction,Walker,move,Buffer) = guidingfuncRatio(F.psi,Walker,move,Buffer)

@inline getOx_k(ψG::SymmetryReducedWaveFunction,W,k) = getOx_k(ψG.psi,W,k)


abstract type AbstractSymop end
struct TranslationalSymmetry{T} <: AbstractSymop
    a1::T
    a2::T
end

function reduceParams!(ψSymm::SymmetryReducedWaveFunction,Symm::AbstractSymop,S::AbstractMatrix)
    indicesMapping = ψSymm.indicesMapping
    uniqueInds = ψSymm.uniqueInds
    params = get_params(ψSymm.psi)

    for par in indicesMapping
        indicesMapping[par] == par || continue

        type,k = _getParamsTypeAndIndex(params,par)
        equivalent_params = generate_equivalent(type,k,Symm,ψSymm.psi,S)
        equivalent_params .= remap_index.(type,equivalent_params,Ref(params))
        # indicesMapping[equivalent_params] .= par
        indicesMapping[equivalent_params] .= minimum(equivalent_params)
    end
    empty!(uniqueInds)
    append!(uniqueInds,findFirstUniqueIndices(indicesMapping))

    for i in eachindex(indicesMapping)
        ind = findfirst(==(indicesMapping[i]),uniqueInds)
        indicesMapping[i] = ind
    end
    return ψSymm
end

generate_equivalent(type,k,Symm::T1,ψSymm::T2,S) where {T1 <:AbstractSymop, T2 <: AbstractGuidingFunction} = error("Symmetry $T1 not implemented for wavefunction $T2")

function generate_equivalent_sites(site::siteType,T::TranslationalSymmetry,S::AbstractMatrix) where {siteType}
    (;a1,a2) = T
    Lx,Ly = size(S)
    newsites = Set(siteType[])
    wrap_offset = 0
    for n1 in -Lx:Lx
        for n2 in -Ly:Ly
            newsite = _translateAndWrap(site,n1*a1 + n2*a2,Lx,Ly,wrap_offset)
            push!(newsites,newsite)
        end
    end
    return collect(newsites)
end

function _translateAndWrap(I::siteType,v,Lx,Ly,wrap_offset=0) where {siteType}
    I´ = Tuple(I) .+ v
    I_wrap = wrap_indices_periodic(I´...,Lx,Ly,wrap_offset)
    return _convert(siteType,I_wrap)
end
_convert(x,y) = convert(x,y)
_convert(::Type{<:CartesianIndex},y) = CartesianIndex(y)

function generate_equivalent_site_pairs(I::siteType,J::siteType,T::TranslationalSymmetry,S::AbstractMatrix) where {siteType}
    (;a1,a2) = T
    Lx,Ly = size(S)
    newPairs = Set{Tuple{siteType,siteType}}()
    
    for n1 in -Lx:Lx
        for n2 in -Ly:Ly
            translation = n1*a1 + n2*a2
            I´ = _translateAndWrap(I,translation,Lx,Ly)
            J´ = _translateAndWrap(J,translation,Lx,Ly)
            push!(newPairs,(I´,J´))
            push!(newPairs,(J´,I´))
        end
    end
    return collect(newPairs)
end

index_to_site(k::Integer,S::AbstractMatrix) = CartesianIndices(S)[k]
site_to_index(site::siteType,S::AbstractMatrix) where {siteType} = LinearIndices(S)[site...]
site_to_index(site::CartesianIndex,S::AbstractMatrix) = LinearIndices(S)[Tuple(site)...]

function index_to_site_pair(k::Integer,S::AbstractMatrix,couplingMatrix::AbstractMatrix)
    I,J = Tuple(CartesianIndices(couplingMatrix)[k])

    return (index_to_site(I,S),index_to_site(J,S))
end

function site_pair_to_index(I::siteType,J::siteType,S::AbstractMatrix,couplingMatrix::AbstractMatrix) where siteType
    I´ = site_to_index(I,S)
    J´ = site_to_index(J,S)
    return LinearIndices(couplingMatrix)[I´,J´]
end

plaquette_to_index(plaquettesite,AllPlaqs) = findfirst(==(plaquettesite),AllPlaqs)

"""given the index k and the type of the parameter and the ArrayPartition containing all parameters, returns the linear index of the parameter. Is the inverse of _getParamsTypeAndIndex"""
Base.@propagate_inbounds function remap_index(partition,k,params::RecursiveArrayTools.ArrayPartition)
    partition == 1 && return k
    lens = length.(params.x[1:partition])
    Base.@boundscheck k > sum(lens) && error("Index out of bounds")
    return sum(lens[1:end-1]) + k
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