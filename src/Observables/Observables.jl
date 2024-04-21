using LatticeFFTs
using LatticeFFTs
using LatticeFFTs.Interpolations

function getStructureFac(AllStates::AbstractVector{<:SpinConfig}, weights, tol = 0)
    plan = getLatticeFFTPlan(AllStates[1].Mat, 0)
    # weights = abs2.(Psi)
    # state_weight = collect(zip( AllStates,weights))[inds]

    Sq =
        fetch.([
            Threads.@spawn getInterpolatedFFT(
                c.Mat,
                0,
                plan;
                Interpolation = BSpline(Constant()),
            ) for c in AllStates
        ])
    # Sq = [getInterpolatedFFT(weight* c.Mat,0,plan;Interpolation = BSpline(Constant())) for (c,weight) in state_weight]
    SSq(kx, ky) = sum(w^2 * s(kx, ky) * s(-kx, -ky) for (w, s) in zip(weights, Sq))

    k = Sq[1].itp.ranges[1]
    Sq_k = fetch.([Threads.@spawn SSq(kx, ky) for kx in k, ky in k])

    return (; k, Sq = SSq, Sq_k)
end

function getStructureFac(
    AllStates::AbstractVector{<:SpinConfig},
    eigen::Union{Eigen,NamedTuple},
    tol = 0,
)
    Psi = eigen.vectors[:, 1]
    NumSites = length(AllStates[1])
    weights = abs2.(Psi) ./ NumSites
    return getStructureFac(AllStates, weights, tol)
end

function getEqualWeightStructureFac(AllStates)
    plan = getLatticeFFTPlan(AllStates[1].Mat, 0)
    Sq =
        fetch.([
            Threads.@spawn getInterpolatedFFT(
                c.Mat,
                0,
                plan;
                Interpolation = BSpline(Constant()),
            ) for c in AllStates
        ])

    NumSites = length(AllStates[1])
    weight(Nstates) = 1 / (Nstates * NumSites)

    SSq(kx, ky) = sum(s(kx, ky) * s(-kx, -ky) for s in Sq) * weight(length(Sq))

    function SSq(kx, ky, maxindex)
        sum(Sq[i](kx, ky) * Sq[i](-kx, -ky) for i = 1:maxindex) * weight(maxindex)
    end

    k = Sq[1].itp.ranges[1]
    Sq_k = fetch.([Threads.@spawn SSq(kx, ky) for kx in k, ky in k])

    return (; k, Sq = SSq, Sq_k)
end

getR(ij::CartesianIndex{2}) = float(SA[ij[1], ij[2]])

function getRij_vec(Config::SpinConfig, i)
    CI = CartesianIndices(Config)
    ri = getR(CI[i])
    rij = [ri - getR(j) for j in CI]
end

function getRij_vec(Config::SpinConfig)
    Ri = reshape(
        [float(SVector(Tuple(ij))) for ij in CartesianIndices(Config.Mat)],
        length(Config),
    )
    return [Ri[i] - Ri[j] for i in eachindex(Ri) for j = 1:i]
end

function getSij(Configs::AbstractVector{<:SpinConfig}, i, j)
    return mean(c[i] * c[j] for c in Configs)
end

function getSij(Configs::AbstractVector{<:SpinConfig}, i)
    return fetch.([Threads.@spawn getSij(Configs, i, j) for j in eachindex(Configs[1])])
end

function getSij(Configs::AbstractVector{<:SpinConfig})
    fac(i, j) = ifelse(i == j, 1, 2)
    return fetch.([
        Threads.@spawn fac(i, j) * getSij(Configs, i, j) for i in LinearIndices(Configs[1])
        for j = 1:i
    ])
end

function getMagnetization(AllStates, eigen::AbstractMatrix)
    ψ0 = eigen.vectors[:, 1]
    return getMagnetization(AllStates, ψ0)
end

function getMagnetization(AllStates, ψ0::AbstractVector)
    mag = zeros(AllStates |> first |> size)
    for n in eachindex(ψ0)
        Si = AllStates[n]
        mag .+= abs2(ψ0[n]) .* Si
    end
    return mag
end
