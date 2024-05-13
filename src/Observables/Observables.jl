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
    
    Lx,Ly = size(Conf)

    resultFull = similar(Sq,Lx+1,Ly+1)
    resultFull .= 0
    result = @view resultFull[1:end-1,1:end-1]

    avgweight = mean(weights)
    avgweight = 1
    avgweight⁻¹ = 1 / avgweight
    for (c, weight) in zip(AllStates, weights)
        Si .= c
        mul!(Sq, plan, Si)
        result .+= abs2.(Sq) .* (avgweight⁻¹*weight)
    end

    kx = (0:Lx) .*(2pi/Lx)
    ky = (0:Ly) .*(2pi/Ly)
    # Sq = [getInterpolatedFFT(weight* c.Mat,0,plan;Interpolation = BSpline(Constant())) for (c,weight) in state_weight]
    # SSq(kx, ky) = sum(w^2 * s(kx, ky) * s(-kx, -ky) for (w, s) in zip(weights, Sq))
    # SSq(kx, ky) = sum(w^2 * abs2(Sq[kx, ky]) for (w, s) in zip(weights, Sq))

    # k = Sq[1].itp.ranges[1]
    # Sq_k = fetch.([Threads.@spawn SSq(Tuple(I)...) for I in CartesianIndices(Sq[1])])
    result = result .* (avgweight)
    resultFull[end,1:end] .= resultFull[1,1:end]
    resultFull[1:end,end] .= resultFull[1:end,1]

    return (; kx,ky, Sq = resultFull)
end

function getSqCont(SqMat)
    Lx,Ly = size(SqMat)

    function Sq(kx,ky)
        # kx´ = mod2pi(kx)
        # ky´ = mod2pi(ky)
        ix = round(Int,kx/2pi*(Lx))
        iy = round(Int,ky/2pi*(Ly))
        ix = rem(ix,Lx)+1
        iy = rem(iy,Ly)+1
        # ix = rem2pi(kx´)*(Lx-1)+1
        return SqMat[ix,iy]
    end

    Sq(k) = Sq(k[1],k[2])
    return Sq
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
    NumSites = length(AllStates[1])
    weights = abs2.(normalize!(ones(length(AllStates)))) ./ NumSites
    println(weights[1])
    return getStructureFacWeights(AllStates, weights)
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


function getBi_square(AllStates,plaqMapping,ψ0::AbstractVector,I)
    res = 0.
    AllStatesDict = Dict(AllStates[i] => i for i in eachindex(AllStates))
    @inline function flipPlaquette(StateRep,ij)
        plaqstate = StateRep[ij]
        setindex(StateRep, ij, !plaqstate)
    end

    for n in eachindex(ψ0)
        x = AllStates[n]

        x´ = flipPlaquette(x, plaqMapping(I))
        if x´ in keys(AllStatesDict)
            m´ = AllStatesDict[x´]
            res += ψ0[n] * ψ0[m´]
        end
    end
    return res
end

function getBij_square(AllStates,plaqMapping,ψ0::AbstractVector,I,J)
    res = 0.
    AllStatesDict = Dict(AllStates[i] => i for i in eachindex(AllStates))
    @inline function flipPlaquette(StateRep,ij)
        plaqstate = StateRep[ij]
        setindex(StateRep, ij, !plaqstate)
    end

    for n in eachindex(ψ0)
        x = AllStates[n]

        x´ = flipPlaquette(x, plaqMapping(I))
        x´ in keys(AllStatesDict) || continue # do not allow virtual tunneling out of Hilbert space
        x´ = flipPlaquette(x´,plaqMapping(J))
        
        if x´ in keys(AllStatesDict)
            m´ = AllStatesDict[x´]
            res += ψ0[n] * ψ0[m´]
        end
    end
    return res
end