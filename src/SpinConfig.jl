import Base: size, getindex, setindex!, iterate, show, copy, hash

struct SpinConfig{T,MatType<:AbstractMatrix{T},T1<:Real} <: AbstractMatrix{T}
    Mat::MatType
    S::T1
end

Base.getindex(S::SpinConfig, i, j) = getindex(S.Mat, i, j)
Base.setindex!(S::SpinConfig, x, i, j) = setindex!(S.Mat, x, i, j)
Base.iterate(S::SpinConfig, i) = iterate(S.Mat, i)
Base.iterate(S::SpinConfig) = iterate(S.Mat)

Base.size(S::SpinConfig) = size(S.Mat)
Base.copy(S::SpinConfig) = SpinConfig(copy(S.Mat), S.S)

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
