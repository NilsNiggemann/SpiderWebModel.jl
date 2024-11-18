"""
Abstract supertype for observables. An observable must be a function that takes a spin configuration and returns an array.
To define a subtype of O of `AbstractObservable`, one must define the following functions:

- obs(::O) returns the buffer array for the output of the observable.
- `O(out,Conf)`: Writes the observable for the given configuration to the preallocated array `out`. Returns out.
- Base.copy(O::O): Returns a copy of the observable.

If no preallocated array is given, the observable defaults to using the buffer array `obs(O)`.
"""
abstract type AbstractObservable end
obs_size(O::AbstractObservable) = size(obs(O))
(O::AbstractObservable)(Conf) = O(obs(O),Conf)
Base.show(io::IO,::MIME"text/plain",O::Obs) where {Obs <: AbstractObservable} = print(io, "$Obs,", " ∈ ", obs_size(O))
Base.display(io::IO,::MIME"text/plain",O::Obs) where {Obs <: AbstractObservable} = print(io, "$Obs,", " ∈ ", obs_size(O))

struct SqFFT{T<:FFTW.FFTWPlan} <: AbstractObservable
    Si::Matrix{ComplexF32}
    Sq::Matrix{ComplexF32}
    plan::T
end
"""
Given the dimension `dims`, allocates a functor that computes the FFT for a given spin configuration of size dims.
"""
function SqFFT(dims)
    Sq = zeros(ComplexF32,dims)
    Si = zeros(ComplexF32,dims)

    plan = FFTW.plan_fft(Si)

    return SqFFT(Si,Sq,plan)
end

obs(FFTSq::SqFFT) = FFTSq.Sq
copy(FFTSq::SqFFT) = SqFFT(copy(FFTSq.Si),copy(FFTSq.Sq),FFTSq.plan)

function (FFTSq::SqFFT)(out,Conf::Matrix{ComplexF32})
    mul!(out, FFTSq.plan, Conf)
    @inbounds for i in eachindex(out)
        out[i] = abs2(out[i])
    end
    out
end

function (FFTSq::SqFFT)(out,Conf)
    copyto!(FFTSq.Si,Conf)
    FFTSq(out,FFTSq.Si)
end

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
    SqFunc = SqFFT(size(Conf))
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

const ABSTRACTCOLLECTION = Union{AbstractRange,AbstractVector,Tuple}

@views function getSqGFMC(res,m_values::ABSTRACTCOLLECTION)
    pMax = maximum(m_values)
    Gnp = precomputeNormalizedAccWeight(res.TotalWeights,1,2pMax)

    Conf = res.SaveConfigs[:,:,begin,begin]
    NSites = length(Conf)

    SqFunc = SqFFT(size(Conf))
    SaveConfs = res.SaveConfigs
    reconfTable = res.reconfigurationTable
    res_m = getObs(Gnp,SaveConfs,reconfTable,SqFunc,m_values)
    newRes_m = [similar(real(res),size(res).+1) for res in res_m]

    for (res,newRes) in zip(res_m,newRes_m)
        newRes[begin:end-1,begin:end-1] .= res ./NSites

        @views newRes[end,begin:end] .= newRes[begin,:]
        @views newRes[begin:end,end] .= newRes[:,begin]
    end

    return newRes_m
    # obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
end
getSqGFMC(res,p,::Nothing) = getSqGFMC(res,p)
# @views function getSqGFMC(res,p,discardborder::Integer)
#     Gnp = precomputeNormalizedAccWeight(res.TotalWeights,1,p)    # Gnp = ones(length(res.TotalWeights[nThermal:end]),p)

#     Conf = res.SaveConfigs[begin+discardborder:end-discardborder,begin+discardborder:end-discardborder,begin,begin]

#     NSites = length(Conf)
#     Sq = similar(Conf, ComplexF32)
    
#     Si = similar(Conf, ComplexF32)
#     plan = FFTW.plan_fft(Si)

#     function SqFunc(Conf)
#         # Si .= Conf
#         # copytont!(Si,Conf)
#         copyto!(Si,Conf[begin+discardborder:end-discardborder,begin+discardborder:end-discardborder])
#         mul!(Sq, plan, Si)
#         for i in eachindex(Sq)
#             Sq[i] = abs2(Sq[i])
#         end
#         Sq
#         # Sq .= abs2.(Sq)
#     end
#     SaveConfs = res.SaveConfigs
#     reconfTable = res.reconfigurationTable
#     res = getObs(Gnp,SaveConfs,reconfTable,SqFunc,p÷2)
#     newRes = similar(res,size(res).+1)
#     newRes[begin:end-1,begin:end-1] .= res

#     @views newRes[end,begin:end] .= newRes[begin,:]
#     @views newRes[begin:end,end] .= newRes[:,begin]
#     return real(newRes ./NSites)
#     # obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
# end

function getSqsGFMC(Results,p)
    Sqs = Vector{Matrix{Float64}}(undef,length(Results))
    Threads.@threads for i in eachindex(Results,Sqs)
        res = Results[i]
        Sq = getSqGFMC(res,p)
        Sqs[i] = Sq
    end
    return Sqs
end

function getSqsGFMC(Results,p::ABSTRACTCOLLECTION)
    Sqs = Vector{Vector{Matrix{Float64}}}(undef,length(Results))
    Threads.@threads for i in eachindex(Results,Sqs)
        res = Results[i]
        Sq = getSqGFMC(res,p)
        Sqs[i] = Sq
    end
    return stack(stack(Sqs))
end
# struct SzOperator <: AbstractOperator
#     I::CartesianIndex{2}
# end
# SzOperator((i,j)) = SzOperator(CartesianIndex(i,j))

# (S::SzOperator)(conf::StencilSpinConfig) = conf[S.I] /2
# (S::SzOperator)(conf::AbstractMatrix{<:AbstractFloat}) = conf[S.I]



function getImagTimeCorr(Gnp,reconfigurationTable,ObsFunc::T,mtau=size(Gnp,2)÷4, m=size(Gnp,2)÷2) where {T}
    N = lastindex(reconfigurationTable,2)
    Obs = [float(zero(ObsFunc(1,2m))) for i in 1:mtau]
    O0 = float(zero(ObsFunc(1,2m)))
    # num = zero(Obs)
    denom = 0.

    Nw = size(reconfigurationTable,1)
    p = size(Gnp,2)

    BranchingMatrix = zeros(Int,Nw,m+1)

    WalkerMultiplicities = zeros(Int,Nw,m+1)
    for n in m+1:N
        Gn = Gnp[n,p]
        denom += Gn*Nw
        
        getBranchingMatrix!(BranchingMatrix,WalkerMultiplicities,reconfigurationTable,n,m)
        for α in 1:Nw
            # O0 = ObsFunc(BranchingMatrix[α,m],n-m)
            O0 = _set_to!(O0,ObsFunc(BranchingMatrix[α,m],n-m))
            for ntau in 0:mtau-1
                mult = WalkerMultiplicities[α,m-ntau]
                mult == 0 && continue
                Otau = ObsFunc(BranchingMatrix[α,m-ntau],n-m+ntau)
                Obs[ntau+1] = _add_numerator!(Obs[ntau+1],Gn*mult,O0,Otau)
                # Obs[ntau+1] += Gn*mult*O0*Otau
            end
        end
    end
    
    for i in eachindex(Obs)
        Obs[i] /= denom
    end
    return (Obs)
end
function getImagTimeCorr(result::GFMCObservables,ObsFunc::T,m,mtau = m÷2) where T
    SWFObsFunc = constructSWF_operator(result.SaveConfigs,ObsFunc)
    Gnp = precomputeNormalizedAccWeight(result.TotalWeights,1,m)
    getImagTimeCorr(Gnp,result.reconfigurationTable,SWFObsFunc,mtau,m)
end
function getImagTimeCorr(results::AbstractVector{<:GFMCObservables},ObsFunc::T,m,mtau = m÷2) where T
    stack(getImagTimeCorr.(results,ObsFunc,m,mtau))
end

_add_numerator!(Obsn::AbstractArray,Gnmult,O0::AbstractArray,Otau::AbstractArray) = (@. Obsn += Gnmult*O0*Otau)
_add_numerator!(Obsn::Number,Gnmult,O0::Number,Otau::Number) = (Obsn += Gnmult*O0*Otau)
_set_to!(a,b) = (a=b)
_set_to!(a::AbstractArray,b::AbstractArray) = (a.=b)

function getBranchingMatrix!(BranchingMatrix::AbstractMatrix,PopulationMatrix,reconfigurationTable::AbstractMatrix,n,projectionLength)
    # BranchingMatrix[:,begin] .= @view reconfigurationTable[:,begin]
    PopulationMatrix .= 0 
    for α in axes(reconfigurationTable,1)
        α´ = α
        for i_m in 0:projectionLength
            α´ = reconfigurationTable[α´,n-i_m]
            # println((; α,α´,i_m))
            BranchingMatrix[α,i_m+1] = α´
            PopulationMatrix[α´,i_m+1] += 1
        end
    end
    return (;BranchingMatrix,PopulationMatrix)
end

function getBranchingMatrix(reconfigurationTable,n,projectionLength) 
    BranchingMatrix = zeros(Int,size(reconfigurationTable,1),projectionLength+1)
    PopulationMatrix = zeros(Int,size(reconfigurationTable,1),projectionLength+1)
    getBranchingMatrix!(BranchingMatrix,PopulationMatrix,reconfigurationTable,n,projectionLength)
end

function constructSWF_operator(AllConfigs,OpFunc::T) where T
    function Obsfunc(α,n)
        conf = @view AllConfigs[:,:,α,n]
        return OpFunc(conf)
    end
end

# function getImagTimeCorr(results::GFMCObservables,ObsFunc::T,mtau=size(Gnp,2)÷4, m=size(Gnp,2)÷2) where {T}
#     Gnps = precomputeNormalizedAccWeight(results,1,m)

#     stack(getImagTimeCorr)
# end

function autocorrelation(SaveConfigs,O,lags,args...;kwargs...)

    dims = size(SaveConfigs)
    # configs = reshape(SaveConfigs,prod(dims[1:3]),dims[4])
    # configs = reshape(SaveConfigs[:,:,begin,:],dims[1]*dims[2],dims[4])
    # Otype = zero(O(SaveConfigs[:,:,begin,begin]))

    # E1 = [zero(Otype) for _ in 1: length(lags)]
    # E2 = [zero(Otype) for _ in 1: length(lags)]
    E1 = zeros(length(lags))
    E2 = zeros(length(lags))
        # var = 0.
    nsum = axes(SaveConfigs,4)[begin+maximum(lags):end]

    N = length(nsum)
    for i in nsum
        for (i_lag,lag) in enumerate(lags)
            O1 = 0.
            O2 = 0.
            Nw = length(axes(SaveConfigs,3))
            for j in axes(SaveConfigs,3)
                conf1 = @view SaveConfigs[:,:,j,i-lag]
                conf2 = @view SaveConfigs[:,:,j,i]
                O1 += O(conf1)
                O2 += O(conf2)
            end

            E1[i_lag] += O1*O2 / Nw
            E2[i_lag] += O1 /Nw
            # var += O1^2

        end
    end
    return (E1 .- E2.^2) ./ N
    # return StatsBase.autocor(configs,args...;kwargs...)
end