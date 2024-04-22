using LatticeFFTs
using LatticeFFTs
using LatticeFFTs.Interpolations

function getStructureFacWeights(AllStates::AbstractVector{<:SpinConfig}, weights, tol = 0)
    # weights = abs2.(Psi)
    # state_weight = collect(zip( AllStates,weights))[inds]
    Conf = parent(AllStates[1])
    Sq = similar(Conf, Complex{eltype(Conf)})
    Si = similar(Conf, Complex{eltype(Conf)})

    plan = LatticeFFTs.FFTW.plan_fft(Conf)
    # Sq(c) = mul!(outBuffer,plan,c)
    
    result = similar(Sq)
    avgweight = mean(weights)
    avgweight⁻¹ = 1 / avgweight
    for (c, weight) in zip(AllStates, weights)
        Si .= c
        mul!(Sq, plan, Si)
        result .+= abs2.(Sq) .* (avgweight⁻¹*weight)^2
    end
    # Sq = [getInterpolatedFFT(weight* c.Mat,0,plan;Interpolation = BSpline(Constant())) for (c,weight) in state_weight]
    # SSq(kx, ky) = sum(w^2 * s(kx, ky) * s(-kx, -ky) for (w, s) in zip(weights, Sq))
    # SSq(kx, ky) = sum(w^2 * abs2(Sq[kx, ky]) for (w, s) in zip(weights, Sq))

    # k = Sq[1].itp.ranges[1]
    # Sq_k = fetch.([Threads.@spawn SSq(Tuple(I)...) for I in CartesianIndices(Sq[1])])

    # return (; k, Sq = SSq, Sq_k)
    return result .* (avgweight^2)
end

function getStructureFac(
    AllStates::AbstractVector{<:SpinConfig},
    ψ0::AbstractVector,
    tol = 0,
)
    NumSites = length(AllStates[1])
    weights = abs2.(ψ0) ./ NumSites
    return getStructureFacWeights(AllStates, weights, tol)
end

function getEqualWeightStructureFac(AllStates)
    plan = getLatticeFFTPlan(Conf, 0)
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

function getSij(AllStates, ψ0::AbstractVector,i,j)
    Sij = 0.
    for n in eachindex(ψ0)
        Si = AllStates[n][i]
        Sj = AllStates[n][j]
        Sij += abs2(ψ0[n]) * Si*Sj
    end
    return Sij
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
