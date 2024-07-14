function getStructureFacWeights(AllStates::AbstractVector{<:SpinConfig}, weights, tol = 0)
    # weights = abs2.(Psi)
    # state_weight = collect(zip( AllStates,weights))[inds]
    Conf = parent(AllStates[1])
    Sq = similar(Conf, Complex{eltype(Conf)})
    Si = similar(Conf, Complex{eltype(Conf)})
    plan = FFTW.plan_fft(Conf)
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

function _convertToInds(k, L)
    i =  k*L/(2pi)
    i = (i + L) % L # Ensure positive indices before modulo
    i = round(Int,i) +1
    i == L+1 && return 1
    return i
end

function getSqCont(SqMat;cutoffEnd = 1)
    SqPlot = @view SqMat[1:end-cutoffEnd,1:end-cutoffEnd]
    Lx,Ly = size(SqPlot)

    function Sq(kx,ky)
        ix = _convertToInds(kx,Lx)
        iy = _convertToInds(ky,Ly)
        return SqPlot[ix,iy]
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

copytont!(B, A) = LoopVectorization.vmapnt!(identity, B, A)
@views function getSqGFMC(res,p)
    Gnp = precomputeNormalizedAccWeight(res.TotalWeights,1,p)    # Gnp = ones(length(res.TotalWeights[nThermal:end]),p)

    Conf = res.SaveConfigs[:,:,begin,begin]
    NSites = length(Conf)
    Sq = similar(Conf, ComplexF32)
    
    Si = similar(Conf, ComplexF32)
    plan = FFTW.plan_fft(Si)

    function SqFunc(Conf)
        # Si .= Conf
        # copytont!(Si,Conf)
        copyto!(Si,Conf)
        mul!(Sq, plan, Si)
        for i in eachindex(Sq)
            Sq[i] = abs2(Sq[i])
        end
        Sq
        # Sq .= abs2.(Sq)
    end
    SaveConfs = res.SaveConfigs
    reconfTable = res.reconfigurationTable
    res = getObs(Gnp,SaveConfs,reconfTable,SqFunc,p÷2)
    newRes = similar(res,size(res).+1)
    newRes[begin:end-1,begin:end-1] .= res

    @views newRes[end,begin:end] .= newRes[begin,:]
    @views newRes[begin:end,end] .= newRes[:,begin]
    return real(newRes ./NSites)
    # obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
end
_getNbra(res,::Nothing) = res.nBra
_getNbra(res,nBra) = nBra

function getSqsGFMC(Results,p,nBra=nothing)

    Sqs = Vector{Matrix{Float64}}(undef,length(Results))
    Threads.@threads for i in eachindex(Results,Sqs)
        res = Results[i]
        _nBra = _getNbra(res,nBra)
        Sq = getSqGFMC(res,p÷_nBra)
        Sqs[i] = Sq
    end
    return Sqs
end

# struct SzOperator <: AbstractOperator
#     I::CartesianIndex{2}
# end
# SzOperator((i,j)) = SzOperator(CartesianIndex(i,j))

# (S::SzOperator)(conf::StencilSpinConfig) = conf[S.I] /2
# (S::SzOperator)(conf::AbstractMatrix{<:AbstractFloat}) = conf[S.I]