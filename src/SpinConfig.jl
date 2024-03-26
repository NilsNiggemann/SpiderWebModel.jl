import Base: size, getindex, setindex!, iterate, show, copy, hash

abstract type AbstractSpinConfig{T} <: AbstractMatrix{T} end
struct SpinConfig{T,MatType<:AbstractMatrix{T},T1<:Real} <: AbstractSpinConfig{T}
    Mat::MatType
    S::T1
end

@inline Base.parent(S::SpinConfig) = S.Mat
@inline getSpin(S::SpinConfig) = S.S

Base.@propagate_inbounds @inline Base.getindex(S::AbstractSpinConfig, i, j) = getindex(parent(S), i, j)
Base.@propagate_inbounds @inline Base.setindex!(S::AbstractSpinConfig, x, i, j) =
    setindex!(parent(S), x, i, j)
@inline Base.iterate(S::AbstractSpinConfig, i) = iterate(parent(S), i)
@inline Base.iterate(S::AbstractSpinConfig) = iterate(parent(S))

@inline Base.size(S::AbstractSpinConfig) = size(parent(S))
@inline Base.copy(S::AbstractSpinConfig) = SpinConfig(copy(parent(S)), getSpin(S))

@inline function plaquetteIterator(S::AbstractMatrix)
    inboundsInds = Base.Iterators.product(axes(S, 1)[begin+1:end-1], axes(S, 2)[begin+1:end-1])
    filterInds = Iterators.filter(ind -> isodd(sum(ind)), inboundsInds)
end

@inline function plaquetteIterator(S::AbstractSpinConfig)
    return plaquetteIterator(parent(S))
end

function booleanSpinConfig(Conf::AbstractMatrix, S::Real = 1 / 2)
    S == 1 / 2 || error("S must be 1/2")
    return SpinConfig(Conf .== 1 / 2, S)
end

booleanSpinConfig(Conf::AbstractSpinConfig) = booleanSpinConfig(Conf.Mat, Conf.S)

function floatSpinConfig(Conf::AbstractMatrix{Bool}, S::Real = 1 / 2)
    S == 1 / 2 || error("S must be 1/2")
    return SpinConfig(Conf .- 1 / 2, S)
end
floatSpinConfig(Conf::AbstractSpinConfig) = floatSpinConfig(parent(Conf), getSpin(Conf))
