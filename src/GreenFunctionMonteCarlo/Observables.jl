add_elementwise!(x::AbstractArray,y) = (x .+= y)
add_elementwise!(x::Number,y::Number) = x + y

mult_elementwise!(x::Number,y::Number) = x * y
function mult_elementwise!(x::AbstractArray,y::Number)
    
    for i in eachindex(x)
        x[i] *= y
    end
    return x
end

divide_elementwise!(x::Number,y::Number) = x / y
divide_elementwise!(x::AbstractArray,y) = (x ./= y)

function getObs(Gnp,AllConfigs,reconfigurationTable,ObsFunc,m::Integer=size(Gnp,2)÷2)
    N = lastindex(AllConfigs,4)
    exampleConf = @view AllConfigs[:,:,begin,begin]
    Obs = ObsFunc(exampleConf)
    num = float(zero(Obs))
    denom = 0.
    GnO = float(zero(Obs))
    # fill!(Obs,zero(eltype(Obs)))

    Nw = size(reconfigurationTable,1)
    p = size(Gnp,2)
    # surviving_walker_mapping_list = zeros(Int,Nw)
    WalkerMultiplicities = zeros(Int,Nw)

    for n in m+1:N
        Gn = Gnp[n,p]
        # denom += Gn*Nw
        denom += Gn*Nw
        # reconfigList = @view reconfigurationTable[:,n-m]
        # surviving_walker_mapping!(surviving_walker_mapping_list,reconfigList)
        WalkerMultiplicities .= 0
        for α in 1:Nw
            α´ = α
            for i_m in 1:m
                α´ = reconfigurationTable[α´,n-i_m]
            end
            WalkerMultiplicities[α´] += 1
        end

        for α in 1:Nw
            mult = WalkerMultiplicities[α]
            mult == 0 && continue
            conf = @view AllConfigs[:,:,α,n-m]
            O = ObsFunc(conf)
            # surviving_index = surviving_walker_mapping_list[α´]
            # O = ObsFunc(AllConfigs[n-m][surviving_index])
            GnO = _set_to!(GnO,O)
            GnO = mult_elementwise!(GnO,Gn*mult)
            # @. num += Gn*O
            num = add_elementwise!(num,GnO)
        end
    end
    return divide_elementwise!(num,denom)
end
function getObs(result,ObsFunc,p::Integer)
    Gnp = precomputeNormalizedAccWeight(result.TotalWeights,1,p)
    return getObs(Gnp,result.SaveConfigs,result.reconfigurationTable,ObsFunc,p÷2)
end

function getObs_smallBuffer(Gnp,AllConfigs,reconfigurationTable,ObsFunc!::AbstractObservable,m_values::ABSTRACTCOLLECTION)
    N = lastindex(AllConfigs,4)

    pMax = maximum(m_values)

    Obs = obs(ObsFunc!)
    num_m = [zeros(size(Obs)) for _ in m_values]
    denom = 0.

    Nw = size(reconfigurationTable,1)
    p = size(Gnp,2)
    WalkerMultiplicities = zeros(Int,Nw)
    ObsBuffer = similar(Obs)
    
    for n in pMax+1:N
        Gn = Gnp[n,p]
        denom += Gn*Nw
        for (i_m,m) in enumerate(m_values)
            WalkerMultiplicities .= 0
            for α in 1:Nw
                α´ = α
                for i_m in 1:m
                    α´ = reconfigurationTable[α´,n-i_m]
                end
                WalkerMultiplicities[α´] += 1
            end


            for α in 1:Nw
                mult = WalkerMultiplicities[α]
                mult == 0 && continue

                conf = @view AllConfigs[:,:,α,n-m]
                ObsFunc!(ObsBuffer,conf)

                @. num_m[i_m] += ObsBuffer*Gn*mult
            end
        end
    end
    for i in eachindex(num_m)
        divide_elementwise!(num_m[i],denom)
    end
    return num_m
end

function getObs(Gnp,AllConfigs,reconfigurationTable,ObsFunc::AbstractObservable,m_values::ABSTRACTCOLLECTION)
    N = lastindex(AllConfigs,4)

    pMax = maximum(m_values)

    Obs = obs(ObsFunc)
    num_m = [zeros(size(Obs)) for _ in m_values]
    denom = 0.

    Nw = size(reconfigurationTable,1)
    p = size(Gnp,2)
    WalkerMultiplicities = zeros(Int,Nw)
    ObsBuffer = [similar(Obs) for α in 1:Nw, n in 1:(pMax)]
    
    wrap_idx(n) = (n-1) % (pMax) + 1
    obsBuffer(α,n) = ObsBuffer[α,wrap_idx(n)]

    _fill_obs_buffer!(ObsBuffer,1:pMax,ObsFunc,AllConfigs,pMax)

    for n in pMax+1:N
        Gn = Gnp[n,p]
        denom += Gn*Nw
        _fill_obs_buffer!(ObsBuffer,n-1,ObsFunc,AllConfigs,pMax)
        for (i_m,m) in enumerate(m_values)
            WalkerMultiplicities .= 0
            for α in 1:Nw
                α´ = α
                for i_m in 1:m
                    α´ = reconfigurationTable[α´,n-i_m]
                end
                WalkerMultiplicities[α´] += 1
            end


            for α in 1:Nw
                mult = WalkerMultiplicities[α]
                mult == 0 && continue

                O = obsBuffer(α,n-m)
                @. num_m[i_m] += O*Gn*mult
            end
        end
    end
    for i in eachindex(num_m)
        divide_elementwise!(num_m[i],denom)
    end
    return num_m
end

function _fill_obs_buffer!(ObsBuffer,nRange,ObsFunc!,AllConfigs,pMax)
    wrap_idx(n) = (n-1) % (pMax) + 1
    obsBuffer(α,n) = ObsBuffer[α,wrap_idx(n)]

    for n in nRange, α in axes(AllConfigs,3)
        conf = @view AllConfigs[:,:,α,n]
        ObsFunc!(obsBuffer(α,n),conf)
    end
    return
end

function surviving_walker_mapping!(mappingarr,reconfigList)
    fill!(mappingarr,0)
    surviving_walkers = 0
    for (α,α´) in enumerate(reconfigList)
        if α == α´
            surviving_walkers += 1
            mappingarr[α´] = surviving_walkers
        end
    end
    # println(reconfigList)
    # println(mappingarr)
    for α in eachindex(reconfigList,mappingarr)
        if mappingarr[α] == 0
            α´ = reconfigList[α]
            # println((α,α´,mappingarr[α´]))
            mappingarr[α] = mappingarr[α´]
        end
    end
    return mappingarr
end

function splitIntoBins(array,binsize)
    Iterators.partition(array,binsize)
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


@views function getSqGFMC(res,m_values::ABSTRACTCOLLECTION;UseBuffer=true)
    pMax = maximum(m_values)
    Gnp = precomputeNormalizedAccWeight(res.TotalWeights,1,2pMax)

    Conf = res.SaveConfigs[:,:,begin,begin]
    NSites = length(Conf)

    SqFunc = SqFFT(size(Conf))
    SaveConfs = res.SaveConfigs
    reconfTable = res.reconfigurationTable
    if UseBuffer
        res_m = getObs(Gnp,SaveConfs,reconfTable,SqFunc,m_values)
    else
        res_m = getObs_smallBuffer(Gnp,SaveConfs,reconfTable,SqFunc,m_values)
    end
    newRes_m = [similar(real(res),size(res).+1) for res in res_m]

    for (res,newRes) in zip(res_m,newRes_m)
        newRes[begin:end-1,begin:end-1] .= res ./NSites

        @views newRes[end,begin:end] .= newRes[begin,:]
        @views newRes[begin:end,end] .= newRes[:,begin]
    end

    return newRes_m
    # obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
end

function getSqsGFMC(Results,p;kwargs...)
    Sqs = Vector{Matrix{Float64}}(undef,length(Results))
    Threads.@threads for i in eachindex(Results,Sqs)
        res = Results[i]
        Sq = getSqGFMC(res,p;kwargs...)
        Sqs[i] = Sq
    end
    return Sqs
end

function getSqsGFMC(Results,p::ABSTRACTCOLLECTION;kwargs...)
    Sqs = Vector{Vector{Matrix{Float64}}}(undef,length(Results))
    Threads.@threads for i in eachindex(Results,Sqs)
        res = Results[i]
        Sq = getSqGFMC(res,p;kwargs...)
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