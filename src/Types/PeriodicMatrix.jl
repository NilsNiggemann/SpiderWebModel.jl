struct PeriodicMatrix{T,Mat<:AbstractMatrix{T}} <: AbstractMatrix{T}
    UC::Mat
    Lx::Int
    Ly::Int
    offset::Int
end
@inline function plaquetteIterator(S::PeriodicMatrix,shift=false)
    inboundsInds = Base.Iterators.product(axes(S, 1), axes(S, 2))
    filterInds = Iterators.filter(ind -> isodd(sum(ind) +shift), inboundsInds)
end

PeriodicMatrix(Mat, Lx, Ly) = PeriodicMatrix(Mat, Lx, Ly, 0)
PeriodicMatrix(Mat) = PeriodicMatrix(Mat, size(Mat)..., 0)

function wrap_index_periodic(i, L)
    i = (i - 1) % L + 1
    if i <= 0
        i += L
    end
    return i
end

function wrap_indices_periodic(x, y, Lx, Ly, offset)
    xRegion = (x - 1) ÷ Lx

    y = y - offset * xRegion

    x = wrap_index_periodic(x, Lx)
    y = wrap_index_periodic(y, Ly)
    return x, y
end

function wrap_indices_periodic(x::UnitRange, y::UnitRange, Lx, Ly, offset)
    Inds = (wrap_indices_periodic(i, j, Lx, Ly, offset) for i in x for j in y)
end

function Base.getindex(P::PeriodicMatrix, i::Integer, j::Integer)
    Lx, Ly = size(P.UC)
    i, j = wrap_indices_periodic(i, j, Lx, Ly, P.offset)

    getindex(P.UC, i, j)
end

function Base.getindex(P::PeriodicMatrix, i, j)
    Lx, Ly = size(P.UC)
    inds = [CartesianIndex(wrap_indices_periodic(ii, jj, Lx, Ly, P.offset)) for ii in i, jj in j]

    getindex(P.UC, inds)
end

function Base.setindex!(P::PeriodicMatrix, x, i::Integer)
    setindex!(P.UC, x, CartesianIndices(P)[i])
end
function Base.setindex!(P::PeriodicMatrix, x, i, j)
    setindex!(P.UC, x, wrap_indices_periodic(i, j, size(P.UC)..., P.offset)...)
end

Base.size(P::PeriodicMatrix) = (P.Lx, P.Ly)
Base.copy(P::PeriodicMatrix) = PeriodicMatrix(copy(P.UC), P.Lx, P.Ly, P.offset)

@inline function Base.checkbounds(::Type{Bool}, arr::PeriodicMatrix, I...)
    i, j = wrap_indices_periodic(I..., size(arr)..., arr.offset)
    # checkbounds(Bool, arr, i,j) || throw_boundserror(arr, ij)
    true
end
