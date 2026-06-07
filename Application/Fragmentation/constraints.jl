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