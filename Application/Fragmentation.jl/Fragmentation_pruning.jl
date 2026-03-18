# first_nz_set: whether a nonzero spin has been placed by an ancestor.
# Symmetry breaking: when false, skip s=-1 (first nonzero must be +1).
include("constraints.jl")
using StaticArrays, SparseArrays, Statistics
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

function plaquette_updates(S)
    S == 1 || S == 0.5 || error("Spin S must be either 1 or 0.5")
    s_fac = S == 1 ? 1 : 2
    P = s_fac*SA[1, -1, -1, 1, 1, -1, -1, 1]

    return (P, -P)
end

function get_domain(S)
    S == 1 || S == 0.5 || error("Spin S must be either 1 or 0.5")
    return S == 1 ? Int8[-1, 0, 1] : Int8[-1, 1]
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
        domain = get_domain(S)
        sols = constraintSolver_compressed_mmap(C, domain, mmap_path = "compressed_solutions.bin", overwrite = true)
        # Build adjacency matrix and find sectors
        clusters = get_clusters(L, L)
        local_updates = plaquette_updates(S)

        H = build_sparse_adjacency_matrix(L, L, sols.data, clusters, local_updates, domain)
        sectors = find_sectors(H)
        unique_sectors = unique(sectors)
        # Determine spin type for formatting

        # Collect results

        total_solutions[i] = correct_time_reversal(length(sectors), S_values[i])
        num_sectors[i] = correct_time_reversal(length(unique_sectors), S_values[i])
        mean_connectivity[i] = mean(get_connectivity(H))

    end

    return (;Ls, S_values, total_solutions, num_sectors, mean_connectivity)
end


function relative_sector_size(num_solutions,num_sectors,L)
    num_sectors<0 && return NaN
    num_solutions<0 && return NaN

    b_rel = (num_solutions/num_sectors)^(1/L^2)
    return b_rel
end

function print_latex_table(results)
    (;Ls, S_values, total_solutions, num_sectors, mean_connectivity) = results

    data = Dict{Tuple{Int,Float64}, NamedTuple}()
    for (L, S, total, sectors, mean_conn) in zip(Ls, S_values, total_solutions, num_sectors, mean_connectivity)
        data[(L, S)] = (total=total, sectors=sectors,
                        relative=round(relative_sector_size(total, sectors, L), digits=3),
                        mean_conn=round(mean_conn, digits=2))
    end
    # delete all values where data[(L, S)]
    spin_order = [0.5, 1.0]
    spin_label(s) = s == 1.0 ? raw"Spin-1" : raw"Spin-\nicefrac{1}{2}"
    function get_val(L, s, field)
        if !haskey(data, (L, s))
            return raw"\textemdash"
        end
        val = getfield(data[(L, s)], field)
        if val == -10 || isnan(val)
            return raw"\textemdash"
        end
        return string(val)
    end

    categories = [
        ("Total", :total),
        ("Sectors", :sectors),
        ("b/b_sec", :relative),
        ("avg. connectivity", :mean_conn),
    ]
    allL = sort(unique(Ls))
    println(raw"\begin{tabular}{l|cc||cc||cc||cc}")
    println(raw"  \toprule")
    println("  \$L\$ & \\multicolumn{2}{c}{Total} & \\multicolumn{2}{c}{Sectors} & \\multicolumn{2}{c}{b/b_{sec}} & \\multicolumn{2}{c}{avg. connectivity} " * "\\\\")
    println(raw"  \cmidrule{2-3} \cmidrule{4-5} \cmidrule{6-7} \cmidrule{8-9}")
    spin_headers = [spin_label(0.5), spin_label(1.0), spin_label(0.5), spin_label(1.0), spin_label(0.5), spin_label(1.0), spin_label(0.5), spin_label(1.0)]
    println("  & $(join(spin_headers, " & ")) " * "\\\\")
    println(raw"  \midrule")

    for L in allL
        vals = String[]
        for (_, field) in categories
            for s in spin_order
                val = get_val(L, s, field)
                push!(vals, val)
            end
        end
        println("  \$$(L)\$ & $(join(vals, " & ")) \\\\")
    end

    println(raw"  \bottomrule")
    println(raw"\end{tabular}")
end
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

##
appended_results = deepcopy(results)
push!(appended_results.Ls,12)
push!(appended_results.S_values,0.5)
push!(appended_results.total_solutions, 10103561614)  # obtained from counting
push!(appended_results.num_sectors, -10)  # unknown without adjacency analysis
push!(appended_results.mean_connectivity, NaN)  # 
push!(appended_results.Ls,8)
push!(appended_results.S_values,1.0)
push!(appended_results.total_solutions, 47067992003)
push!(appended_results.num_sectors, -10)  # unknown without adjacency analysis
push!(appended_results.mean_connectivity, NaN)  #


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

    # Index entries by (L, S) for side-by-side lookup
    data = Dict{Tuple{Int,Float64}, NamedTuple}()
    for (L, S, total, sectors, mean_conn) in zip(Ls, S_values, total_solutions, num_sectors, mean_connectivity)
        data[(L, S)] = (total=total, sectors=sectors, relative=round(relative_sector_size(total, sectors, L), digits=2), mean_conn=round(mean_conn, digits=2))
    end

    all_L    = sort(unique(Ls))
    all_spin = sort(unique(S_values))  # e.g. [0.5, 1.0]
    spin_label(s) = s == 1.0 ? "Spin-1" : "Spin-1/2"
    get_val(L, s, field) = haskey(data, (L, s)) ? string(getfield(data[(L, s)], field)) : "—"

    categories = [("Sectors", :sectors), ("Total", :total), ("b/b_sec", :relative), ("MeanConn", :mean_conn)]

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
print_latex_table(appended_results)

##
function countmap(v::AbstractVector)
    m = Dict{eltype(v), Int}()
    for x in v
        m[x] = get(m, x, 0) + 1
    end
    return m
end

function sector_analysis(Lx,Ly,S)
    C = setup_constraints(Lx, Ly)
    domain = get_domain(S)
    local_updates = plaquette_updates(S)
    file = "Solutions_$(Lx)x$(Ly)_$(S)"*tempname()
    mkpath(dirname(file))
    sols = constraintSolver_compressed_mmap(C, domain, mmap_path = file, overwrite = true)
    clusters = get_clusters(Lx, Ly)
    H = build_sparse_adjacency_matrix(Lx, Ly, sols.data, clusters, local_updates, domain)
    sectors = find_sectors(H)

    sector_sizes = countmap(sectors)
    sector_sizes_all = [sector_sizes[s] for s in sectors]

    return (;sols, sectors, H, sector_sizes, sector_sizes_all)
end

Lx = 6
Ly = 6

res_05 = sector_analysis(Lx, Ly, 0.5)
res_1 = sector_analysis(Lx, Ly, 1.0)

##


##
using CairoMakie, MakieHelpers
with_theme(theme_SimpleTicks()) do
    # sector_counts = countmap(sectors)
    # unique_sectors = sort(collect(keys(sector_counts)))
    # counts = [sector_counts[s] for s in unique_sectors]


    fig = Figure(size = 0.7 .*(800, 450))
    ax05 = Axis(fig[1, 1], xlabel="Sector size", ylabel="Count",
    #  title=L"S=1/2",
    yscale = log10,xscale = log10,
    xlabelvisible = false, 
    xticklabelsvisible = false, 
    )
    ax1 = Axis(fig[2, 1], xlabel="Sector size", ylabel="Count", 
    # title=L"S=1",
    yscale = log10,xscale = log10
    )
    linkxaxes!(ax05, ax1)
    # unique_sectors_05 =  sort(collect(keys(res_05.sector_sizes)))

    # counts05 = [res_05.sector_sizes[s] for s in unique_sectors_05]
    # counts1 = [res_1.sector_sizes[s] for s in unique_sectors_05]
    sector_counts = countmap(res_05.sectors)
    unique_sectors = sort(collect(keys(sector_counts)))
    counts = [sector_counts[s] for s in unique_sectors]


    hist!(ax05,counts, color=:black,bins=10)

    sector_counts = countmap(res_1.sectors)
    unique_sectors = sort(collect(keys(sector_counts)))
    counts = [sector_counts[s] for s in unique_sectors]

    hist!(ax1,counts, color=:black,bins=600)

    text!(ax05, (0.98,0.95), text=L"S=1/2", space = :relative, fontsize = 20, align = (:right, :top) )
    text!(ax1, (0.98,0.95), text=L"S=1", space = :relative, fontsize = 20, align = (:right, :top) )
    Label(fig[1, 1,TopLeft()], L"(a)$$", padding = (-40,0,-25, -20),fontsize = 16)
    Label(fig[2, 1,TopLeft()], L"(b)$$", padding = (-40,0,-25, -20),fontsize = 16)
    save("Application/figs/PaperFigs/sector_size_histograms.pdf", fig)
    fig
end
##
# decompress_spinconfig.([(x,) for (x,) in zip(sols.data[1, 1:10])], 36)
import SpiderWebModel as SW
S = SW.stencilConfig(zeros(Lx, Ly), 1, boundaryCondition = :periodic)

max_size = maximum(values(sector_sizes))
max_sec_num = 0
for (k,v) in sector_sizes
    if v == max_size
        max_sec_num = k
        break
    end
end

S[:].= decompress_spinconfig((sols.data[:, 5713]...,), Lx*Ly)
# S .= 0
# S[1,1] = 1
# S[4,6] = -1
SW.plotApplPlaquettes(S)

##
S_prime = copy(S)
circshift!(S_prime,S,(1,1))
SW.plotApplPlaquettes(S_prime)
Spr_flat = S_prime[:]
##
decompressed = decompress_spinconfig.([(sols.data[:, i]...,) for i in 1:size(sols.data, 2)], Lx*Ly)
shifted_Spr = findfirst(==(Spr_flat),decompressed)
##
sectors[shifted_Spr]
sector_sizes[sectors[shifted_Spr]]

