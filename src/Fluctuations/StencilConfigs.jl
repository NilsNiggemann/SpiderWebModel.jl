# struct StencilSpinConfig{T,MatType<:AbstractMatrix{T}} <: AbstractSpinConfig{T}
#     Mat::MatType
#     M::Int8
# end
struct StencilSpinConfig{MVal<:Val,T,MatType<:AbstractMatrix{T}} <: AbstractSpinConfig{T}
    Mat::MatType
    M::MVal
end

@inline Base.parent(S::StencilSpinConfig) = S.Mat
@inline getSpin(::StencilSpinConfig{Val{MVal}}) where {MVal} = MVal/2
@inline Base.copy(S::StencilSpinConfig) = StencilSpinConfig(copy(parent(S)), S.M)
function stencilConfig(A::AbstractMatrix{<:AbstractFloat}, S,paddingValue = Int8(typemax(Int8)); kwargs...)
    return stencilConfig(Int8.(2 .*A), S,paddingValue; kwargs...)
end

function stencilConfig(A::AbstractMatrix{Int8}, S,paddingValue = Int8(typemax(Int8)); kwargs...)
    St = Stencils.StencilArray(
        A,
        Stencils.Moore();
        boundary = Stencils.Remove(paddingValue),
        kwargs...,
    )
    M = Val(Int(2S)) # do not chose Int8 to avoid overflows
    return StencilSpinConfig(St, M)
end

stencilConfig(A::SpinConfig; kwargs...) = stencilConfig(parent(A), A.S; kwargs...)

function SpinConfig(S::StencilSpinConfig)
    return SpinConfig(parent(S)./2, getSpin(S))
end

allSpinsInBounds(S::StencilSpinConfig{Val{M}}; kwargs...) where M = allSpinsInBounds(S,M; kwargs...)

plotSpinConfig!(ax, S::StencilSpinConfig; kwargs...) = plotSpinConfig!(ax, SpinConfig(S); kwargs...)
plotSpinConfig(S::StencilSpinConfig; kwargs...) = plotSpinConfig(SpinConfig(S); kwargs...)

function plaquetteIterator(S::Stencils.StencilArray,shift=false)
    bound = Stencils.boundary(S)
    return stencil_plaquetteIterator(S,bound,shift)
end

function stencil_plaquetteIterator(S,::Stencils.Wrap,shift=false)
    inboundsInds = Base.Iterators.product(axes(S, 1), axes(S, 2))
    filterInds = Iterators.filter(ind -> isodd(sum(ind)+shift), inboundsInds)
end
stencil_plaquetteIterator(S,::Any,shift=false) = plaquetteIterator(parent(S),shift)    


@inline Base.@propagate_inbounds getPlaquetteStencil(
    S::Stencils.StencilArray,
    i::Int,
    j::Int,
) = Stencils.stencil(S, CartesianIndex(i, j))

@inline Base.@propagate_inbounds getPlaquetteStencil(S::StencilSpinConfig, i::Int, j::Int) =
    getPlaquetteStencil(parent(S), i, j)

@inline Base.@propagate_inbounds function getSitesFromStencil(S::Stencils.Moore)
    S6, S7, S8, S5, S1, S4, S3, S2 = S
    return SVector(S1, S2, S3, S4, S5, S6, S7, S8)
end

const P1_STENCIL = Int8.(SVector(-2, -2, 2, 2, 2, 2, -2, -2))

@inline Base.@propagate_inbounds function getPlaquetteSites(
    S::Stencils.StencilArray,
    i::Int,
    j::Int,
)
    S = getPlaquetteStencil(S, i, j)
    return getSitesFromStencil(S)
end

@inline getPlaquette(S::StencilSpinConfig, i, j) = getPlaquetteSites(S, i, j)

@inline Base.@propagate_inbounds getPlaquetteSites(S::StencilSpinConfig, i::Int, j::Int) =
    getPlaquetteSites(parent(S), i, j)

Base.@propagate_inbounds function applyPlaquette!(Config,i,j,sgn)
    applyPlaquette!(parent(Config), i, j, sgn)
    return Config
end

Base.@propagate_inbounds function applyPlaquette!(Mat::Stencils.AbstractStencilArray,i,j,sgn)
    # sites = Stencils.indices(Stencils.stencil(Mat), CartesianIndex(i, j))
    sites = safe_indices(Mat, (i, j))
    for (ij,s) in zip(sites,P1_STENCIL)
        i,j = ij
        Mat[i,j] += sgn *s
    end
    return Mat
end

@inline function safe_indices(A::Stencils.AbstractStencilArray,I)
    return safe_indices(A, Stencils.boundary(A), Stencils.padding(A), CartesianIndex(I))
end
@inline function safe_indices(A::Stencils.AbstractStencilArray{<:Any,R,<:Any,N}, ::Stencils.Wrap, ::Stencils.Conditional,I ::CartesianIndex
) where {R,N}
    inds = Stencils.indices(Stencils.stencil(A), I)
    radii = CartesianIndex(ntuple(_ -> -R, N))
    if checkbounds(Bool, A, I + radii) && checkbounds(Bool, A, I - radii)
        return inds
    else
        return get_wrappend_inds.(Ref(A), inds)
    end
end
@inline function safe_indices(A::Stencils.AbstractStencilArray,::Stencils.Remove,::Stencils.Conditional,I::CartesianIndex)
    inds = Stencils.indices(Stencils.stencil(A), I)
end
@inline function safe_indices(A::Stencils.AbstractStencilArray,::Any,::Any,I::CartesianIndex)
    error("setindex for generic Stencil not implemented yet")
end

# `Conditional` needs handling for specific boundary conditions.
# For Wrap we swap the side.
function get_wrappend_inds(A::Stencils.AbstractStencilArray{S,R}, I::Tuple) where {S,R}
    sz = Stencils.tuple_contents(S)
    wrapped_inds = map(I, sz) do i, s
        i < 1 ? i + s : (i > s ? i - s : i)
    end
    return wrapped_inds
end


"""returns a tuple of two booleans, corresponding to the applicability of P and P⁺. M=2S is an integer"""
@inline Base.@propagate_inbounds function P_applicable(Pij::Union{SVector{8},Stencils.Moore}, M::Integer)
    S1, S2, S3, S4, S5, S6, S7, S8 = Pij
    v = (S1, -S2, -S3, S4, S5, -S6, -S7, S8)

    p = all(<(M), v)
    pdagger = all(>(-M), v)
    return p, pdagger
end

"""checks whether F or F† can be applied. Pij are 2S_i in counterclockwise order. M = 2S is twice the spin."""
@inline Base.@propagate_inbounds function P_applicable(Pij::Union{SVector{8},Stencils.Moore}, ::Val{M}) where {M}
    # more performant than above since (S1 - 1) * (S2 + 1) ... is zero only if S1 == 1 and S2 == -1 ...
    # here we technically compute e.g. (2S1 - 2) * (2S2 + 2) ... with the same result
    S1, S2, S3, S4, S5, S6, S7, S8 = Pij
    p = !iszero((S1-M) * (S2+M) * (S3+M) * (S4-M) * (S5-M) * (S6+M) * (S7+M) * (S8-M))
    pdag = !iszero((S1+M) * (S2-M) * (S3-M) * (S4+M) * (S5+M) * (S6-M) * (S7-M) * (S8+M))
    return p, pdag
end

@inline Base.@propagate_inbounds function P_applicable(S::StencilSpinConfig, i::Int, j::Int)
    Pij = getPlaquetteSites(S, i, j)
    P_applicable(Pij, S.M)
end

@inline Base.@propagate_inbounds P_applicable(S::StencilSpinConfig, I::Union{<:NTuple{2},<:CartesianIndex{2}}) =
    P_applicable(S, I[1], I[2])

function getApplicablePlaquettes_separated!(
    plaqsPlus,
    plaqsMinus,
    Conf::StencilSpinConfig,
)
    empty!(plaqsPlus)
    empty!(plaqsMinus)
    for I in plaquetteIterator(Conf)
        applPlus, applMinus = P_applicable(Conf, I)
        if applPlus
            push!(plaqsPlus, Tuple(I))
        end
        if applMinus
            push!(plaqsMinus, Tuple(I))
        end
    end
    return plaqsPlus, plaqsMinus
end

function getApplicablePlaquettes!(
    plaqs,
    Conf::StencilSpinConfig,
)
    empty!(plaqs)
    for I in plaquetteIterator(Conf)
        applPlus, applMinus = P_applicable(Conf, I)
        if applPlus || applMinus
            push!(plaqs, Tuple(I))
        end
    end
    return plaqs
end

getApplicablePlaquettes(Conf::StencilSpinConfig) = getApplicablePlaquettes!( Tuple{Int,Int}[], Conf)

getApplicablePlaquettes(Conf::StencilSpinConfig,::Nothing) = getApplicablePlaquettes(Conf)