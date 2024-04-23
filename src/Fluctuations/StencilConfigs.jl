struct StencilSpinConfig{T,MatType<:AbstractMatrix{T}} <: AbstractSpinConfig{T}
    Mat::MatType
    M::Int8
end
@inline Base.parent(S::StencilSpinConfig) = S.Mat
@inline getSpin(S::StencilSpinConfig) = S.M/2
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
    M = Int8(2S)
    return StencilSpinConfig(St, M)
end

stencilConfig(A::SpinConfig; kwargs...) = stencilConfig(parent(A), A.S; kwargs...)

function SpinConfig(S::StencilSpinConfig)
    return SpinConfig(parent(S)./2, S.M/2)
end

allSpinsInBounds(S::StencilSpinConfig; kwargs...) = allSpinsInBounds(S,S.M; kwargs...)

plotSpinConfig!(ax, S::StencilSpinConfig; kwargs...) = plotSpinConfig!(ax, SpinConfig(S); kwargs...)
plotSpinConfig(S::StencilSpinConfig; kwargs...) = plotSpinConfig(SpinConfig(S); kwargs...)

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

@inline Base.@propagate_inbounds getPlaquetteSites(S::StencilSpinConfig, i::Int, j::Int) =
    getPlaquetteSites(parent(S), i, j)

Base.@propagate_inbounds function applyPlaquette!(Config,i,j,sgn)
    sites = Stencils.indices(Stencils.stencil(parent(Config)), CartesianIndex(i, j))
    for (ij,s) in zip(sites,P1_STENCIL)
        i,j = ij
        Config[i,j] += sgn *s
    end
    return Config
end


"""returns a tuple of two booleans, corresponding to the applicability of P and P⁺. M=2S is an integer"""
@inline Base.@propagate_inbounds function P_applicable(Pij::Union{SVector{8},Stencils.Moore}, M::Integer)
    S1, S2, S3, S4, S5, S6, S7, S8 = Pij
    v = (S1, -S2, -S3, S4, S5, -S6, -S7, S8)

    p = all(<(M), v)
    pdagger = all(>(-M), v)
    return p, pdagger
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