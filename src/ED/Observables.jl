function getStructureFac(
    AllStates::AbstractVector{<:SpinConfig},
    ψ0::AbstractVector,
    tol = 0,
)
    NumSites = length(AllStates[1])
    weights = abs2.(ψ0) ./ NumSites
    return getStructureFacWeights(AllStates, weights, tol)
end


function getStructureFacWeights(AllStates::AbstractVector{<:AbstractMatrix}, weights, tol = 0)
    # weights = abs2.(Psi)
    # state_weight = collect(zip( AllStates,weights))[inds]
    Conf = parent(AllStates[1])
    Sq = zeros(ComplexF32,size(Conf))
    Si = zeros(ComplexF32,size(Conf))
    plan = FFTW.plan_fft(Sq)
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

function getEqualWeightStructureFac(AllStates)
    NumSites = length(AllStates[1])
    weights = abs2.(normalize!(ones(length(AllStates)))) ./ NumSites
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

function getTauCorr(AllStates,eigen,tau,O::T) where T
    (;values,vectors) = eigen
    v0 = vectors[:,1]
    # S_I = getindex.(AllStates,I_Cart)

    # res = zeros(length(tau))
    res = [zero(O(AllStates[1])) for _ in tau]
    e0 = values[1]

    for n in eachindex(values)
        vn = @view vectors[:,n]
        S_0n = complex(zero(res[1]))
        for (i,Conf) in enumerate(AllStates)
            S_0n += O(Conf) * vn[i]*v0[i]
        end
        for (ti,t) in enumerate(tau)
            res[ti] = _add_res_ED!(res[ti],S_0n, exp(-(values[n]-e0)*t))
            # res[ti] += abs2(S_0n) * exp(-(values[n]-e0)*t)
        end
    end
    return res
end
_add_res_ED!(res_i::AbstractArray,S_0n::AbstractArray,expB) = (@. res_i += abs2(S_0n) * expB)
_add_res_ED!(res_i::Number,S_0n::Number,expB) = (res_i += abs2(S_0n) * expB)

function getGSObsED(AllStates,v0,O::T) where T
    res = zero(O(first(AllStates)))
    for (i,Conf) in enumerate(AllStates)
        res += O(Conf) * abs2(v0[i])
    end
    return res
end