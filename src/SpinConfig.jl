import Base: size, getindex, setindex!, iterate, show, copy, hash

struct SpinConfig{T,MatType<:AbstractMatrix{T},T1<:Real} <: AbstractMatrix{T}
    Mat::MatType
    S::T1
end

Base.@propagate_inbounds @inline Base.getindex(S::SpinConfig, i, j) = getindex(S.Mat, i, j)
Base.@propagate_inbounds @inline Base.setindex!(S::SpinConfig, x, i, j) = setindex!(S.Mat, x, i, j)
@inline Base.iterate(S::SpinConfig, i) = iterate(S.Mat, i)
@inline Base.iterate(S::SpinConfig) = iterate(S.Mat)
@inline Base.parent(S::SpinConfig) = S.Mat

@inline Base.size(S::SpinConfig) = size(S.Mat)
@inline Base.copy(S::SpinConfig) = SpinConfig(copy(S.Mat), S.S)

@inline function plaquetteIterator(S::AbstractMatrix)
    return Base.Iterators.product(axes(S,1)[begin+1:end-1], axes(S,2)[begin+1:end-1])
end

@inline function plaquetteIterator(S::SpinConfig)
    return plaquetteIterator(parent(S))
end

function booleanSpinConfig(Conf::AbstractMatrix, S::Real = 1 / 2)
    S == 1 / 2 || error("S must be 1/2")
    return SpinConfig(Conf .== 1 / 2, S)
end

booleanSpinConfig(Conf::SpinConfig) = booleanSpinConfig(Conf.Mat, Conf.S)

function floatSpinConfig(Conf::AbstractMatrix{Bool}, S::Real = 1 / 2)
    S == 1 / 2 || error("S must be 1/2")
    return SpinConfig(Conf .- 1 / 2, S)
end
floatSpinConfig(Conf::SpinConfig) = floatSpinConfig(Conf.Mat, Conf.S)
