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

function print_progress(explored, total, last_print, found_solutions, always_print = false)
    progress = Float64(explored[]) / Float64(total)
    if progress - last_print[] > 0.02 || progress == 1.0 || always_print
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
    solutions = Vector{Vector{Int}}()

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
    subtree = [big(d)^(N-i) for i in 1:N]

    total = big(d)^N
    explored = Ref(0.0)
    last_print = Ref(0.0)

    dfs(1, S, C, domain, residual, remaining_bound, subtree, solutions, explored, total, last_print, N, M)
    print_progress(explored, total, last_print, length(solutions), true)
    println("\nDone.")
    return solutions
end


##
C = setup_constraints(6,6)
domain = Int8[-1, 0, 1]
# domain = Int8[-1, 1]
@time solutions = constraintSolver(C, domain);