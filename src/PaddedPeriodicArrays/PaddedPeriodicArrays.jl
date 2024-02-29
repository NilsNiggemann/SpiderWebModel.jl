import CircularArrays

struct PeriodicPaddedMatrix{T, Mat <: AbstractMatrix{T},SubMat <: SubArray,PerMat} <: AbstractMatrix{T}
    A::SubMat
    A_period::PerMat
    A_pad::Mat
    pad::Int
end

"""implement a periodic matrix where the border is padded by ghost cells.
This makes accessing the matrix elements more efficient"""
function periodicPaddedMatrix(Ainit,pad)
    i1,i2 = extrema(axes(Ainit,1))
    j1,j2 = extrema(axes(Ainit,2))
    APeriod_temp = CircularArrays.CircularArray(Ainit)

    A_pad = Array(APeriod_temp[i1-pad:i2+pad,j1-pad:j2+pad])
    A = @view A_pad[i1+pad:i2+pad,j1+pad:j2+pad]

    # @info "" i1 i2 size(Ainit) size(A_pad) size(A)
    APeriod = CircularArrays.CircularArray(A)

    return PeriodicPaddedMatrix(A,APeriod,A_pad,pad)
end

Base.getindex(P::PeriodicPaddedMatrix, i::Integer) = getindex(P.A_pad,i)

Base.getindex(P::PeriodicPaddedMatrix, i::Integer, j::Integer) = getindex(P.A_pad,i+P.pad,j+P.pad)

Base.setindex!(P::PeriodicPaddedMatrix, x, I::Vararg{T,N}) where {T,N} = setindex!(P.A_period, x, I...)

Base.size(P::PeriodicPaddedMatrix) = size(P.A)

Base.copy(P::PeriodicPaddedMatrix) = periodicPaddedMatrix(P.A,P.pad)
Base.parent(P::PeriodicPaddedMatrix) = P.A_pad

function updateGhostCells!(P::PeriodicPaddedMatrix)
    for i in axes(P.A_pad,1), j in axes(P.A_pad,2)
        if i <= P.pad || j <= P.pad || i > size(P.A_pad,1)-P.pad || j > size(P.A_pad,2)-P.pad
            P.A_pad[i,j] = P.A_period[i-P.pad,j-P.pad]
        end
    end
end

@inline function Base.checkbounds(::Type{Bool},arr::PeriodicPaddedMatrix, I::Vararg{T,N}) where {T,N}
    # return checkbounds(Bool, arr.A, J...)
    return true
    # return all()
end

@inline function Base.checkbounds(arr::PeriodicPaddedMatrix, I::Vararg{T,N}) where {T,N}
    # return checkbounds(Bool, arr.A, J...)
    return true
    # return all()
end

function Base.checkbounds(A::PeriodicPaddedMatrix, I...)
    @inline
    nothing
end

# function spinConfig!(Conf::SpinConfig{<:PeriodicPaddedMatrix}, path::AbstractVector, InitialConf, plaqMap::PlaqMapping)
#     spinConfig!(Conf.Mat, path, InitialConf, plaqMap)
#     updateGhostCells!(Conf.Mat)
#     return Conf
# end

function spinConfig!(Conf::SpinConfig{T,<:PeriodicPaddedMatrix}, path::SBitVector, InitConf::AbstractMatrix, plaqMap::PlaqMapping) where T
    Conf.Mat.A_pad .= InitConf.Mat.A_pad

    for (i, op) in enumerate(path) 
        if op
            ij = plaqMap(i)
            flipPlaquette!(Conf, ij)
        end
    end
    updateGhostCells!(Conf.Mat)
    return Conf
end