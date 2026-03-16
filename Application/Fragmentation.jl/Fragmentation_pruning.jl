# first_nz_set: whether a nonzero spin has been placed by an ancestor.
# Symmetry breaking: when false, skip s=-1 (first nonzero must be +1).
include("constraints.jl")
using StaticArrays, SparseArrays
function compress_spinconfig(S::AbstractVector{<:Integer}, ::Val{N_split}) where {N_split}
    @assert N_split >= 1 "N_split must be at least 1"
    len = length(S)
    @assert all(s -> s in -1:1, S) "Spin values must be in the range -1, 0, 1"

    chunk_len = cld(len, N_split)
    @assert chunk_len <= 39 "Per-split length must be at most 39 (Int64 base-3 limit); increase N_split"

    configs = zeros(Int64, N_split)
    for split in 1:N_split
        start_idx = (split - 1) * chunk_len + 1
        end_idx = min(split * chunk_len, len)
        start_idx > len && break

        acc = Int64(0)
        pow3 = Int64(1)
        for i in start_idx:end_idx
            acc += Int64(S[i] + 1) * pow3
            pow3 *= Int64(3)
        end
        configs[split] = acc
    end

    return Tuple(configs)
end

function decompress_spinconfig(config::NTuple{N_split,Int64}, N::Int) where {N_split}
    @assert N_split >= 1 "N_split must be at least 1"
    S = Vector{Int8}(undef, N)

    chunk_len = cld(N, N_split)
    @assert chunk_len <= 39 "Per-split length must be at most 39 (Int64 base-3 limit); increase N_split"

    for split in 1:N_split
        start_idx = (split - 1) * chunk_len + 1
        end_idx = min(split * chunk_len, N)
        start_idx > N && break

        acc = config[split]
        for i in start_idx:end_idx
            S[i] = Int8(acc % Int64(3) - 1)
            acc ÷= Int64(3)
        end
    end

    return S
end

# function compress_spinconfig(S::AbstractVector{<:Integer}, ::Val{2})
#     len = length(S)
#     @assert iseven(len) "Length of S must be even for Val(2) compression"
#     @assert len <= 78 "Length of S must be at most 78 for Val(2) compression"
#     return compress_spinconfig(S, Val{2}())
# end

# function decompress_spinconfig(config::Tuple{Int64,Int64}, N::Int)
#     @assert iseven(N) "N must be even for Val(2) decompression"
#     return decompress_spinconfig(config, N, Val{2}())
# end

function dfs_write_compressed!(i, first_nz_set::Bool, S, C, domain, residual, remaining_bound, subtree, io, found, explored, total, last_print, N, M, split_val::Val{N_split}) where {N_split}
    if i > N
        explored[] += one(explored[])
        if all(residual .== 0)
            found[] += 1
            compressed = compress_spinconfig(S, split_val)
            for part in compressed
                write(io, Int64(part))
            end
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
            dfs_write_compressed!(i + 1, first_nz_set || s != 0, S, C, domain, residual, remaining_bound, subtree, io, found, explored, total, last_print, N, M, split_val)
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
        mmap_path::AbstractString = "compressed_solutions.bin", overwrite::Bool = false,
        n_split::Union{Nothing,Int} = nothing)
    M, N = size(C)
    nsplit = isnothing(n_split) ? cld(N, 39) : n_split
    @assert nsplit >= 1 "n_split must be at least 1"
    @assert cld(N, nsplit) <= 39 "Per-split length exceeds 39; increase n_split"
    @assert all(s -> s in -1:1, domain) "Domain values must be in the range -1, 0, 1"
    split_val = Val(nsplit)

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
    dfs_write_compressed!(1, false, S, C, domain, residual, remaining_bound, subtree, io, found, explored, total, last_print, N, M, split_val)
    print_progress(explored, total, last_print, found[], true)
    println("\nDone writing ", found[], " compressed solutions.")
    flush(io)
    close(io)

    nbytes = filesize(path)
    recbytes = nsplit * sizeof(Int64)
    @assert nbytes % recbytes == 0 "Corrupt compressed solution file size: $(nbytes) bytes"
    nsolutions = Int(nbytes ÷ recbytes)

    if nsolutions == 0
        return (path = path, nsolutions = 0, data = Matrix{Int64}(undef, nsplit, 0), nsplit = nsplit)
    end

    io_map = open(path, "r+")
    storage = Mmap.mmap(io_map, Matrix{Int64}, (nsplit, nsolutions))
    close(io_map)

    return (path = path, nsolutions = nsolutions, data = storage, nsplit = nsplit)
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
    nsplit = size(sols_data, 1)
    split_val = Val(nsplit)
    # adjacency = spzeros(Int, nconfigs, nconfigs)  # Sparse adjacency matrix
    adjacency_rows = Int[]
    adjacency_cols = Int[]
    Nsites = Lx * Ly
    # Create a hash set for quick lookup of compressed configurations
    # compressed_set = Set([(sols_data[1, i], sols_data[2, i]) for i in 1:nconfigs])
    compressed_dict = Dict{Tuple,Int}()
    for i in 1:nconfigs
        key = ntuple(k -> sols_data[k, i], nsplit)
        compressed_dict[key] = i
    end

    first_key = ntuple(k -> sols_data[k, 1], nsplit)
    new_config = decompress_spinconfig(first_key, Nsites)
    for i in 1:nconfigs
        # Decompress the current configuration
        config_key = ntuple(k -> sols_data[k, i], nsplit)
        config = decompress_spinconfig(config_key, Nsites)

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
                compressed_new = compress_spinconfig(new_config, split_val)
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

function find_sectors(H::SparseMatrixCSC)
    sectors = collect(1:size(H, 2))
    for (i,s) in enumerate(sectors)
        connected = findnz(H[:, i])[1]
        for c in connected
            sectors[c] = min(s,sectors[c])
        end
    end
    return sectors
end

function get_connectivity(H::SparseMatrixCSC)
    connected = zeros(Int, size(H, 2))
    for i in axes(H, 2)
        connected[i] = length(findnz(H[:, i])[1])
    end
    return connected
end

function analyze_sectors_and_solutions(domain_list)
    Ls = [L for (L, S) in domain_list]
    S_values = [S for (L, S) in domain_list]
    total_solutions = [0 for _ in domain_list]
    num_sectors = [0 for _ in domain_list]
    mean_connectivity = [0.0 for _ in domain_list]
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
        mean_connectivity[i] = mean(get_connectivity(H))

    end

    return (;Ls, S_values, total_solutions, num_sectors, mean_connectivity)
end


function relative_sector_size(num_solutions,num_sectors)
    return num_sectors / num_solutions
end

function print_latex_table(results)
    (;Ls, S_values, total_solutions, num_sectors, mean_connectivity) = results
    total_solutions = correct_time_reversal.(total_solutions, S_values)
    num_sectors     = correct_time_reversal.(num_sectors, S_values)

    data = Dict{Tuple{Int,Float64}, NamedTuple}()
    for (L, S, total, sectors, mean_conn) in zip(Ls, S_values, total_solutions, num_sectors, mean_connectivity)
        data[(L, S)] = (total=total, sectors=sectors,
                        relative=round(relative_sector_size(total, sectors)*100, digits=2),
                        mean_conn=round(mean_conn, digits=2))
    end

    all_L    = sort(unique(Ls))
    all_spin = sort(unique(S_values))
    nspin    = length(all_spin)
    spin_label(s) = s == 1.0 ? raw"Spin-1" : raw"Spin-\nicefrac{1}{2}"
    get_val(L, s, field) = haskey(data, (L, s)) ? string(getfield(data[(L, s)], field)) : raw"\textemdash"

    categories = [("Sectors", :sectors), ("Total", :total), (raw"Rel.\ (\%)", :relative), (raw"\(\bar{z}\)", :mean_conn)]
    ncats = length(categories)

    # Total data columns = nspin * ncats
    ncols = nspin * ncats

    # Column spec: l | (c's grouped by category, separated by ||)
    col_groups = join([join(fill("c", nspin), "") for _ in 1:ncats], "||")
    println(raw"\begin{table}[htbp]")
    println(raw"  \centering")
    println("  \\begin{tabular}{l|$(col_groups)}")
    println(raw"  \toprule")

    # Row 1: category multicolumns
    cat_headers = join(["\\multicolumn{$(nspin)}{c}{$(cat_label)}" for (cat_label, _) in categories], " & ")
    println("  \$L\$ & $(cat_headers) \\\\")

    # Cmidrules under each category group
    cmidrules = join(["\\cmidrule(lr){$((k-1)*nspin+2)-$((k-1)*nspin+nspin+1)}" for k in 1:ncats], " ")
    println("  $(cmidrules)")

    # Row 2: spin sub-headers repeated for each category
    spin_headers = join(repeat([join([spin_label(s) for s in all_spin], " & ")], ncats), " & ")
    println("  & $(spin_headers) \\\\")
    println(raw"  \midrule")

    # Data rows
    for L in all_L
        vals = String[]
        for (_, field) in categories
            for s in all_spin
                push!(vals, get_val(L, s, field))
            end
        end
        println("  \$$(L)\$ & $(join(vals, " & ")) \\\\")
    end

    println(raw"  \bottomrule")
    println("  \\end{tabular}")
    println("  \\caption{Sector analysis for spin-1/2 and spin-1 on \$L\\times L\$ lattices.}")
    println(raw"  \label{tab:sectors}")
    println(raw"\end{table}")
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
H = build_sparse_adjacency_matrix(Lx, Ly, sols.data, clusters, local_updates, domain)
sectors = find_sectors(H)

unique_secs = unique(sectors)

##
# Example usage
domain_list = [
    (4, 0.5),
    (6, 0.5),
    (8, 0.5),
    (10, 0.5),
    (4, 1.),
    (6, 1.),
    # (8, 1),
]
results = analyze_sectors_and_solutions(domain_list)

##    println("L   Domain       Unique Sectors   Total Solutions")
function correct_time_reversal(num_solutions, Spin)
    if Spin == 0.5
        return 2num_solutions  # Account for time-reversal pairs
    else
        return 2num_solutions + 1  # Account for the zero-spin configuration which is its own pair
    end
end

let 
    (;Ls, S_values, total_solutions, num_sectors, mean_connectivity) = results
    total_solutions = correct_time_reversal.(total_solutions, S_values)
    num_sectors = correct_time_reversal.(num_sectors, S_values)

    # Index entries by (L, S) for side-by-side lookup
    data = Dict{Tuple{Int,Float64}, NamedTuple}()
    for (L, S, total, sectors, mean_conn) in zip(Ls, S_values, total_solutions, num_sectors, mean_connectivity)
        data[(L, S)] = (total=total, sectors=sectors, relative=round(relative_sector_size(total, sectors)*100, digits=2), mean_conn=round(mean_conn, digits=2))
    end

    all_L    = sort(unique(Ls))
    all_spin = sort(unique(S_values))  # e.g. [0.5, 1.0]
    spin_label(s) = s == 1.0 ? "Spin-1" : "Spin-1/2"
    get_val(L, s, field) = haskey(data, (L, s)) ? string(getfield(data[(L, s)], field)) : "—"

    categories = [("Sectors", :sectors), ("Total", :total), ("Rel(%)", :relative), ("MeanConn", :mean_conn)]

    # Column widths: L col + per-spin col within each category
    lw = 4   # L column
    cw = 14  # per-spin sub-column

    nspin = length(all_spin)

    # Row 1: category group headers
    row1 = rpad("L", lw)
    for (cat_label, _) in categories
        row1 *= rpad(cat_label, cw * nspin)
    end
    println(row1)

    # Row 2: spin sub-headers
    row2 = rpad("", lw)
    for _ in categories
        for s in all_spin
            row2 *= rpad(spin_label(s), cw)
        end
    end
    println(row2)

    println(repeat("-", lw + cw * nspin * length(categories)))

    for L in all_L
        row = rpad(string(L), lw)
        for (_, field) in categories
            for s in all_spin
                row *= rpad(get_val(L, s, field), cw)
            end
        end
        println(row)
    end
end

##
print_latex_table(results)