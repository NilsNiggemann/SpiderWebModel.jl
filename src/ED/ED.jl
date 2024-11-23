SolveH(H, range = 1:1) = eigen(H, range)
const SparseMat = Union{SparseMatrixCSC,Hermitian{<:Number,<:SparseMatrixCSC}}

function SolveH(H::SparseMat; kwargs...)
    if size(H) == (1, 1)
        return (; values = [float(real(only(H)))], vectors = [1.0])
    end
    values, vectors = eigs(H, nev = 1, which = :SR, explicittransform = :none; kwargs...)
    return (; values, vectors)
end

function SolveHKrylov(H; kwargs...)
    values, vectors, _ = KrylovKit.eigsolve(H, ones(size(H,1)), 1, :SR; kwargs...)
    return (; values, vectors)
end
