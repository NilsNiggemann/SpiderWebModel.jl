# first_nz_set: whether a nonzero spin has been placed by an ancestor.
# Symmetry breaking: when false, skip s=-1 (first nonzero must be +1).
include("constraints.jl")
using StaticArrays, SparseArrays
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

using Printf
function print_progress(explored, total, last_print, found_solutions, always_print = false)
    current_time = time()
    progress = Float32(explored[]) / Float32(total)

    if current_time - last_print[] > 60 || always_print
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

function get_clusters(Lx,Ly)
    CI = CartesianIndices((Lx, Ly))
    LI = LinearIndices((Lx, Ly))

    idx_wrap(i,j) = LI[mod1(i, Lx), mod1(j, Ly)]
    idx_wrap(I) = idx_wrap(Tuple(I)...)

    signs = (1,1,-1,-1,1,1,-1,-1)
    
    neighborsites(i, j) = SA[
        (i, j+1),
        (i-1, j+1),
        (i-1, j),
        (i-1, j-1),
        (i, j-1),
        (i+1, j-1),
        (i+1, j),
        (i+1, j+1),
    ]

    function neighborsites_linear(i,j)
        map(neighborsites(i,j)) do (ii,jj)
            idx_wrap(ii,jj)
        end
    end

    clusters = SVector{8,Int}[]

    for x in 1:Lx, y in 1:Ly
        isodd(x+y) || continue
        neighbors = neighborsites_linear(x, y)
        push!(clusters, neighbors)
    end
    return clusters
end
# @time sols = constraintSolver_compressed_mmap(C, domain, mmap_path = "compressed_solutions.bin", overwrite = true)

function build_sparse_adjacency_matrix(Lx, Ly, sols_data, clusters, local_updates, domain)
    nconfigs = size(sols_data, 2)  # Number of configurations
    # adjacency = spzeros(Int, nconfigs, nconfigs)  # Sparse adjacency matrix
    adjacency_rows = Int[]
    adjacency_cols = Int[]
    Nsites = Lx * Ly
    # Create a hash set for quick lookup of compressed configurations
    # compressed_set = Set([(sols_data[1, i], sols_data[2, i]) for i in 1:nconfigs])
    compressed_dict = Dict([(sols_data[1, i], sols_data[2, i]) => i for i in 1:nconfigs])

    new_config = decompress_spinconfig((sols_data[1, 1], sols_data[2, 1]), Nsites) 
    for i in 1:nconfigs
        # Decompress the current configuration
        config = decompress_spinconfig((sols_data[1, i], sols_data[2, i]), Nsites) 

        for cluster in clusters
            # Apply the local update
            new_config .= config  # Start with the original configuration
            # S_i = config[cluster] 

            S_i_prime = @view new_config[cluster] 
            for move in local_updates
                S_i_prime .+= move  # Apply the local update to the cluster
                
                if new_config == config
                    print(stderr, S_i_prime, move)
                    error("Local update did not change the configuration. Check the local_updates and clusters definitions.")
                end
                is_valid_update = all(∈(domain),S_i_prime)
                if !is_valid_update
                    S_i_prime .-= move
                    continue  # Skip invalid configurations
                end

                # println(i)
                # Compress the new configuration
                compressed_new = compress_spinconfig(new_config, Val(2))
                S_i_prime .-= move  # Revert the change for the next iteration
                # Check if the new configuration exists in sols_data
                if haskey(compressed_dict, compressed_new)
                    # Find the index of the new configuration
                    j = compressed_dict[compressed_new]
                    # Add an edge to the adjacency matrix
                    # adjacency[i, j] = 1
                    push!(adjacency_rows, i)
                    push!(adjacency_cols, j)
                end
            end
        end
    end
    adjacency = sparse(adjacency_rows, adjacency_cols, ones(Int, length(adjacency_rows)), nconfigs, nconfigs)
    return adjacency
end

function find_sectors(H::SparseMatrixCSC{Int, Int})
    sectors = collect(1:size(H, 2))
    for (i,s) in enumerate(sectors)
        connected = findnz(H[:, i])[1]
        for c in connected
            sectors[c] = min(s,sectors[c])
        end
    end
    return sectors
end



function analyze_sectors_and_solutions(domain_list)
    Ls = [L for (L, S) in domain_list]
    S_values = [S for (L, S) in domain_list]
    total_solutions = [0 for _ in domain_list]
    num_sectors = [0 for _ in domain_list]

    for i in eachindex(domain_list)
        L = Ls[i]
        S = S_values[i]
        # Setup constraints and solve
        @info "Processing L=$L, S=$(join(S, ","))"
        C = setup_constraints(L, L)
        domain = S == 1 ? Int8[-1, 0, 1] : Int8[-1, 1]
        sols = constraintSolver_compressed_mmap(C, domain, mmap_path = "compressed_solutions.bin", overwrite = true)
        spin_fac = S == 1 ? 1 : 2  # correct local updates
        # Build adjacency matrix and find sectors
        clusters = get_clusters(L, L)
        local_updates = (
            spin_fac * SA[1, -1, -1, 1, 1, -1, -1, 1],
            spin_fac * SA[-1, 1, 1, -1, -1, 1, 1, -1],
        )
        H = build_sparse_adjacency_matrix(L, L, sols.data, clusters, local_updates, domain)
        sectors = find_sectors(H)
        unique_sectors = unique(sectors)
        # Determine spin type for formatting

        # Collect results

        total_solutions[i] = length(sectors)
        num_sectors[i] = length(unique_sectors)
        

    end

    return (;Ls, S_values, total_solutions, num_sectors)
end


function relative_sector_size(num_solutions,num_sectors)
    return num_sectors / num_solutions
end
##
using StaticArrays
Lx = 4
Ly = 4
C = setup_constraints(Lx, Ly)
# domain = Int8[-1, 1]
domain = Int8[-1,0, 1]
# @time solutions = constraintSolver(C, domain)
local_updates = (
    SA[1, -1, -1, 1, 1, -1, -1, 1],
    SA[-1, 1, 1, -1, -1, 1, 1, -1],
)

@time sols = constraintSolver_compressed_mmap(C, domain, mmap_path = "compressed_solutions.bin", overwrite = true)

clusters = get_clusters(Lx, Ly)
@profview H = build_sparse_adjacency_matrix(Lx, Ly, sols.data, clusters, local_updates, domain)
sectors = find_sectors(H)

unique_secs = unique(sectors)

##
# Example usage
domain_list = [
    (4, 0.5),
    (6, 0.5),
    (8, 0.5),
    # (10, 0.5),
    (4, 1),
    (6, 1),
    # (8, 0.5),
]
results = analyze_sectors_and_solutions(domain_list)

##    println("L   Domain       Unique Sectors   Total Solutions")
let 
    (;Ls, S_values, total_solutions, num_sectors) = results
    # Group results by spin type for relative analysis

    println("L   Spin   Unique Sectors   Total Solutions   Relative (%)")
    println("-------------------------------------------------------------------------------")

    for (L, S, total, sectors) in zip(Ls, S_values, total_solutions, num_sectors)
        relative_sectors = round(relative_sector_size(total, sectors) * 100, digits = 2)


        println(rpad(string(L), 4), rpad(S == 1 ? "Spin-1" : "Spin-1/2", 12), rpad(string(sectors), 17), rpad(string(total), 17), rpad(relative_sectors, 20))
    end
    
end