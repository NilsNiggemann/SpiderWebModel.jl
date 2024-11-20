"""
Abstract supertype for observables. An observable must be a function that takes a spin configuration and returns an array.
To define a subtype of O of `AbstractObservable`, one must define the following functions:

- obs(::O) returns the buffer array for the output of the observable.
- `O(out,Conf)`: Writes the observable for the given configuration to the preallocated array `out`. Returns out.
- Base.copy(O::O): Returns a copy of the observable.

If no preallocated array is given, the observable defaults to using the buffer array `obs(O)`.
"""
abstract type AbstractObservable end
obs_size(O::AbstractObservable) = size(obs(O))
(O::AbstractObservable)(Conf) = O(obs(O),Conf)
Base.show(io::IO,::MIME"text/plain",O::Obs) where {Obs <: AbstractObservable} = print(io, "$Obs,", " ∈ ", obs_size(O))
Base.display(io::IO,::MIME"text/plain",O::Obs) where {Obs <: AbstractObservable} = print(io, "$Obs,", " ∈ ", obs_size(O))

struct SqFFT{T<:FFTW.FFTWPlan} <: AbstractObservable
    Si::Matrix{ComplexF32}
    Sq::Matrix{ComplexF32}
    plan::T
end
"""
Given the dimension `dims`, allocates a functor that computes the FFT for a given spin configuration of size dims.
"""
function SqFFT(dims)
    Sq = zeros(ComplexF32,dims)
    Si = zeros(ComplexF32,dims)

    plan = FFTW.plan_fft(Si)

    return SqFFT(Si,Sq,plan)
end

obs(FFTSq::SqFFT) = FFTSq.Sq
copy(FFTSq::SqFFT) = SqFFT(copy(FFTSq.Si),copy(FFTSq.Sq),FFTSq.plan)

function (FFTSq::SqFFT)(out,Conf::Matrix{ComplexF32})
    mul!(out, FFTSq.plan, Conf)
    @inbounds for i in eachindex(out)
        out[i] = abs2(out[i])
    end
    out
end

function (FFTSq::SqFFT)(out,Conf)
    copyto!(FFTSq.Si,Conf)
    FFTSq(out,FFTSq.Si)
end