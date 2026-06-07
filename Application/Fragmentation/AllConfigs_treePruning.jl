include("constraints.jl")
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

##

# Thread-safe, pure-return DFS. No shared mutable state — safe to call from Threads.@threads.
function dfs_count_local!(i, first_nz_set::Bool, S, C, domain, residual, remaining_bound, N, M)
    if i > N
        return all(residual .== 0) ? 1 : 0
    end
    count = 0
    for s in domain
        if !first_nz_set && s < 0
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
            count += dfs_count_local!(i + 1, first_nz_set || s != 0, S, C, domain, residual, remaining_bound, N, M)
        end
        for n in 1:M
            residual[n] -= C[n, i] * s
        end
    end
    return count
end

# Splits the search tree at `split_depth` and processes subtrees in parallel threads.
# Increase split_depth for more tasks (better load-balancing with many threads).
function countSolutions_threaded(C::AbstractMatrix{<:Integer}, domain::AbstractVector{<:Integer};
        use_time_reversal::Bool = true, split_depth::Int = 5)
    M, N = size(C)
    maxspin = maximum(abs.(domain))
    absC = abs.(C)
    remaining_bound = zeros(Int, M, N + 2)
    for n in 1:M
        for i in N:-1:1
            remaining_bound[n, i] = remaining_bound[n, i + 1] + absC[n, i] * maxspin
        end
    end

    d = min(split_depth, N)
    task_list = Vector{Tuple{Vector{Int}, Vector{Int}, Bool}}()

    function collect_tasks!(i, first_nz_set, S, residual)
        if i > d
            push!(task_list, (copy(S), copy(residual), first_nz_set))
            return
        end
        for s in domain
            if !first_nz_set && s < 0
                continue
            end
            S[i] = s
            for n in 1:M; residual[n] += C[n, i] * s; end
            feasible = true
            for n in 1:M
                if abs(residual[n]) > remaining_bound[n, i + 1]
                    feasible = false; break
                end
            end
            if feasible
                collect_tasks!(i + 1, first_nz_set || s != 0, S, residual)
            end
            for n in 1:M; residual[n] -= C[n, i] * s; end
        end
    end

    S0 = zeros(Int, N)
    residual0 = zeros(Int, M)
    collect_tasks!(1, !use_time_reversal, S0, residual0)

    ntasks = length(task_list)
    println("Dispatching $ntasks subtree tasks across $(Threads.nthreads()) threads...")

    counts = zeros(Int, ntasks)
    progress_lock = ReentrantLock()
    completed = Ref(0)
    report_every = max(1, ntasks ÷ 200)

    Threads.@threads for t in eachindex(task_list)
        S_t, res_t, fnz_t = task_list[t]
        counts[t] = dfs_count_local!(d + 1, fnz_t, S_t, C, domain, res_t, remaining_bound, N, M)
        lock(progress_lock) do
            completed[] += 1
            if completed[] % report_every == 0 || completed[] == ntasks
                @printf "\r%.1f%% tasks done (%d/%d)   " 100.0 * completed[] / ntasks completed[] ntasks
                flush(stdout)
            end
        end
    end

    total_count = sum(counts)
    S_zero_state_possible = (0 in domain)
    nsolutions = use_time_reversal ? 2 * total_count - S_zero_state_possible : total_count
    println("\nTotal solutions: $nsolutions", use_time_reversal ? " (reconstructed from time reversal of half-space count ($total_count ))" : "")
    return nsolutions
end

##
C = setup_constraints(4,4)
domain = Int8[-1, 1]
# @time solutions = constraintSolver(C, domain)
##
# @time sols = constraintSolver_compressed_mmap(C, domain, mmap_path = "compressed_solutions.bin", overwrite = true)
@time sols = countSolutions_threaded(C, domain; split_depth = 8)
##
# newConstraints = add_constraints_to_matrix(C, 4, Dict(:STot => (-1, -1)))

##
# Ls = [4, 6, 8, 10]
# counts = [216, 5912, 350872, 37403668]
# bs = counts .^ (1 ./(Ls.^2))
# lines(1 ./Ls, bs)
# xlims!(0, 1/minimum(Ls) + 0.01)
# ylims!(1,1.5)
# current_figure()