function setup_constraints(Lx,Ly)

    N = Lx * Ly
    C_ni = zeros(Int, N÷2, N)

    CI = CartesianIndices((Lx, Ly))
    LI = LinearIndices((Lx, Ly))

    idx_wrap(i,j) = LI[mod1(i, Lx), mod1(j, Ly)]
    idx_wrap(I) = idx_wrap(Tuple(I)...)

    signs = (1,1,-1,-1,1,1,-1,-1)
    
    neighborsites(i, j) = (
        (i, j+1),
        (i-1, j+1),
        (i-1, j),
        (i-1, j-1),
        (i, j-1),
        (i+1, j-1),
        (i+1, j),
        (i+1, j+1),
    )

    function neighborsites_linear(i,j)
        map(neighborsites(i,j)) do (ii,jj)
            idx_wrap(ii,jj)
        end
    end
    n = 1
    for x in 1:Lx, y in 1:Ly
        iseven(x+y) || continue
        neighbors = neighborsites_linear(x, y)
        for (k, neighbor) in enumerate(neighbors)
            C_ni[n, neighbor] = signs[k]
        end
        n += 1
    end
    return C_ni
end
##

function dfs(i, S, C, domain, residual, remaining_bound, subtree, solutions, explored, total, last_print, N, M)
    if i > N
        explored[] += one(explored[])
        print_progress(explored, total, last_print, length(solutions))

        if all(residual .== 0)
            push!(solutions, copy(S))
        end
        return
    end

    for s in domain
        S[i] = s

        for n in 1:M
            residual[n] += C[n,i]*s
        end

        feasible = true
        for n in 1:M
            if abs(residual[n]) > remaining_bound[n,i+1]
                feasible = false
                break
            end
        end

        if feasible
            dfs(i+1, S, C, domain, residual, remaining_bound, subtree, solutions, explored, total, last_print, N, M)
        else
            explored[] += Float64(subtree[i])
            print_progress(explored, total, last_print, length(solutions))
        end

        for n in 1:M
            residual[n] -= C[n,i]*s
        end
    end
end

#= function compressed_spinconfig(S::Vector{Int8})
    config = 0
    for (i, s) in enumerate(S)
        config += (s + 1) * Int64(3)^(i-1)
    end
    return config
end
function decompress_spinconfig(config::Int64, N::Int)
    S = Vector{Int8}(undef, N)
    for i in 1:N
        S[i] = (config % Int64(3)) - 1
        config ÷= Int64(3)
    end
    return S
end =#

function compress_spinconfig(S::AbstractVector{<:Integer}, ::Val{2})
    len = length(S)
    @assert iseven(len) "Length of S must be even for Val(2) compression"
    @assert len <= 78 "Length of S must be at most 78 for Val(2) compression"
    @assert all(s -> s in -1:1, S) "Spin values must be in the range -1, 0, 1"
    config1 = 0
    config2 = 0

    for i in eachindex(S)[1:end÷2]
        s = S[i]
        config1 += (s + 1) * Int64(3)^(i - 1)
        s2 = S[end÷2+i]
        config2 += (s2 + 1) * Int64(3)^(i - 1)
    end
    return (config1, config2)
end

function decompress_spinconfig(config::Tuple{Int64,Int64}, N::Int)
    @assert iseven(N) "N must be even for Val(2) decompression"
    S = Vector{Int8}(undef, N)
    config1 = config[1]
    config2 = config[2]
    for i in 1:N÷2
        S[i] = (config1 % Int64(3)) - 1
        config1 ÷= Int64(3)
        S[end÷2+i] = (config2 % Int64(3)) - 1
        config2 ÷= Int64(3)
    end
    return S
end

function constraintSolver(C::AbstractMatrix{<:Integer}, domain::AbstractVector{<:Integer})
    M, N = size(C)

    S = zeros(Int, N)
    solutions = Vector{Vector{Int8}}()

    maxspin = maximum(abs.(domain))
    absC = abs.(C)

    remaining_bound = zeros(Int, M, N+2)
    for n in 1:M
        for i in N:-1:1
            remaining_bound[n,i] = remaining_bound[n,i+1] + absC[n,i]*maxspin
        end
    end

    residual = zeros(Int, M)

    d = length(domain)
    subtree = [float(d)^(N-i) for i in 1:N]

    total = float(d)^N
    explored = Ref(0.0)
    last_print = Ref(0.0)

    dfs(1, S, C, domain, residual, remaining_bound, subtree, solutions, explored, total, last_print, N, M)
    print_progress(explored, total, last_print, length(solutions), true)
    println("\nDone.")
    return solutions
end

# first_nz_set: whether a nonzero spin has been placed by an ancestor.
# Symmetry breaking: when false, skip s=-1 (first nonzero must be +1).
function dfs_write_compressed!(i, first_nz_set::Bool, S, C, domain, residual, remaining_bound, subtree, io, found, explored, total, last_print, N, M)
    if i > N
        explored[] += one(explored[])
        if all(residual .== 0)
            found[] += 1
            config1, config2 = compress_spinconfig(S, Val(2))
            write(io, Int64(config1))
            write(io, Int64(config2))
        end
        print_progress(explored, total, last_print, found[])
        return
    end

    for s in domain
        # time-reversal symmetry breaking: once the first nonzero is fixed to +1,
        # the mirror image (-S) is excluded from the search.
        if !first_nz_set && s < 0
            explored[] += Float64(subtree[i])
            print_progress(explored, total, last_print, found[])
            continue
        end

        S[i] = s

        for n in 1:M
            residual[n] += C[n, i] * s
        end

        feasible = true
        for n in 1:M
            if abs(residual[n]) > remaining_bound[n, i + 1]
                feasible = false
                break
            end
        end

        if feasible
            dfs_write_compressed!(i + 1, first_nz_set || s != 0, S, C, domain, residual, remaining_bound, subtree, io, found, explored, total, last_print, N, M)
        else
            explored[] += Float64(subtree[i])
            print_progress(explored, total, last_print, found[])
        end

        for n in 1:M
            residual[n] -= C[n, i] * s
        end
    end
end

using Mmap

function constraintSolver_compressed_mmap(C::AbstractMatrix{<:Integer}, domain::AbstractVector{<:Integer};
        mmap_path::AbstractString = "compressed_solutions.bin", overwrite::Bool = false)
    M, N = size(C)
    @assert iseven(N) "N must be even for Val(2) compression"
    @assert N <= 78 "N must be at most 78 for Val(2) compression"
    @assert all(s -> s in -1:1, domain) "Domain values must be in the range -1, 0, 1"

    maxspin = maximum(abs.(domain))
    absC = abs.(C)
    remaining_bound = zeros(Int, M, N + 2)
    for n in 1:M
        for i in N:-1:1
            remaining_bound[n, i] = remaining_bound[n, i + 1] + absC[n, i] * maxspin
        end
    end

    d = length(domain)
    subtree = [float(d)^(N - i) for i in 1:N]
    total = float(d)^N

    path = abspath(mmap_path)
    if isfile(path)
        overwrite || error("Output file already exists: $(path). Set overwrite=true to replace it.")
        rm(path)
    end

    println("Streaming compressed solutions to disk...")
    io = open(path, "w+")
    S = zeros(Int, N)
    residual = zeros(Int, M)
    found = Ref(0)
    explored = Ref(0.0)
    last_print = Ref(0.0)
    dfs_write_compressed!(1, false, S, C, domain, residual, remaining_bound, subtree, io, found, explored, total, last_print, N, M)
    print_progress(explored, total, last_print, found[], true)
    println("\nDone writing ", found[], " compressed solutions.")
    flush(io)
    close(io)

    nbytes = filesize(path)
    recbytes = 2 * sizeof(Int64)
    @assert nbytes % recbytes == 0 "Corrupt compressed solution file size: $(nbytes) bytes"
    nsolutions = Int(nbytes ÷ recbytes)

    if nsolutions == 0
        return (path = path, nsolutions = 0, data = Matrix{Int64}(undef, 2, 0))
    end

    io_map = open(path, "r+")
    storage = Mmap.mmap(io_map, Matrix{Int64}, (2, nsolutions))
    close(io_map)

    return (path = path, nsolutions = nsolutions, data = storage)
end


function print_progress(explored, total, last_print, found_solutions, always_print = false)
    progress = Float64(explored[]) / Float64(total)
    if progress - last_print[] > 0.001 || progress == 1.0 || always_print
        barlen = 40
        filled = round(Int, progress * barlen)
        bar = repeat("█", filled) * repeat(" ", barlen - filled)

        print("\r[$bar] $(round(progress * 100, digits = 2))% ",
                "found solutions: $(found_solutions)", " explored: $(explored[]) / $(total)",
                )
        flush(stdout)
        last_print[] = progress
    end
end
using Printf
function print_progress(explored, total, last_print, found_solutions, always_print = false)
    current_time = time()
    progress = Float32(explored[]) / Float32(total)

    if current_time - last_print[] > 0.5 || progress == 1.0 || always_print
        barlen = 40
        filled = round(Int, progress * barlen)
        bar = repeat("█", filled) * repeat(" ", barlen - filled)

        # print("\r[$bar] $(round(progress * 100, digits = 2))% ",
        #         "found solutions: $(found_solutions) ",
        #         "explored: $(explored[]) / $(total)       ")
        print("\r[$bar] $(round(progress * 100, digits = 2))% ")
        @printf "found solutions: %i explored: %.3e / %.3e   " found_solutions  explored[] total
        flush(stdout)
        last_print[] = current_time
    end
end
##
##
C = setup_constraints(6,6)
domain = Int8[-1, 0, 1]
@time solutions = constraintSolver(C, domain)
##
@time sols = constraintSolver_compressed_mmap(C, domain, mmap_path = "compressed_solutions.bin", overwrite = true)
##
# newConstraints = add_constraints_to_matrix(C, 4, Dict(:STot => (-1, -1)))