struct RBMSpin1{T} <: AbstractGuidingFunction
    a_i::Vector{T}
    A_i::Vector{T}
    b_j::Vector{T}
    w_ij::Matrix{T}
    W_ij::Matrix{T}
    N::Int
    hidden_density::Int
end

function RBMSpin1(N,hidden_density::Int,Type = Float32)
    M = N*hidden_density
    a_i = zeros(Type,N)
    A_i = zeros(Type,N)
    b_j = zeros(Type,M)
    w_ij = zeros(Type,N,M)
    W_ij = zeros(Type,N,M)

    ψ = RBMSpin1(a_i,A_i,b_j,w_ij,W_ij,N,hidden_density)
    return ψ
end
RBMSpin1(S::StencilSpinConfig,hidden_density::Int,Type = Float32) = RBMSpin1(length(S),hidden_density,Type)

get_alpha_i(ψG::RBMSpin1) = ψG.a_i
get_A_i(ψG::RBMSpin1) = ψG.A_i
get_b_j(ψG::RBMSpin1) = ψG.b_j
get_w_ij(ψG::RBMSpin1) = ψG.w_ij
get_W_ij(ψG::RBMSpin1) = ψG.W_ij

function get_params(ψG::RBMSpin1)
    return RecursiveArrayTools.ArrayPartition(ψG.a_i,ψG.A_i,ψG.b_j,ψG.w_ij,ψG.W_ij)
end

struct RBMSpin1Buffer{T,D}
    Θ::Vector{T}
    oneplustanhΘ::Vector{T}
    oneminustanhΘ::Vector{T}
    ΔΘ::Vector{T}
    x::Vector{T}
    Wji::Matrix{T}
    wjiSi_plus_WjiS²::D
end

function _compute_Weight_Sum_RBM1(wij,Wij,move,Conf::StencilSpinConfig)
    i_plaq,j_plaq,opSign = move

    sites = safe_parent_indices(parent(Conf), (i_plaq,j_plaq))
    res = zeros(length(axes(Wij,2)))

    LI = LinearIndices(Conf)
    for j in axes(Wij,2)
        for idx in eachindex(sites)
            I_x,I_y = sites[idx]
            i = LI[I_x,I_y]
            s = P1_STENCIL[idx]*opSign
            res[j] += wij[i,j]*s + Wij[i,j]*s^2
        end
    end

    return res
end

function _precompute_PlaqFlipWeights(wij,Wij,Conf::StencilSpinConfig)
    AllPlaqs = collect(plaquetteIterator(Conf))
    AllWeights = Dict((i,j,s) => _compute_Weight_Sum_RBM1(wij,Wij,(i,j,s),Conf) for (i,j) in AllPlaqs for s in (-1,1)) 
end

function _precompute_PlaqFlipWeights!(AllWeights,wij,Wij,Conf::StencilSpinConfig)
    AllPlaqs = collect(plaquetteIterator(Conf))
    for (i,j) in AllPlaqs
        for s in (-1,1)
            AllWeights[(i,j,s)] = _compute_Weight_Sum_RBM1(wij,Wij,(i,j,s),Conf)
        end
    end
    return AllWeights
end

function allocate_GWF_buffer(ψG::RBMSpin1,S::AbstractMatrix) 
    type = eltype(get_params(ψG))

    Θ = zeros(type,size(get_b_j(ψG)))
    oneplustanhΘ = zeros(type,size(get_b_j(ψG)))
    oneminustanhΘ = zeros(type,size(get_b_j(ψG)))
    ΔΘ = zeros(type,size(get_b_j(ψG)))
    x = zeros(type,length(S))
    x² = zeros(type,length(S))

    wij = get_w_ij(ψG)
    Wij = get_W_ij(ψG)

    Wji = Array(Wij')

    wjiSi_plus_WjiS² = _precompute_PlaqFlipWeights(wij,Wij,S)

    return RBMSpin1Buffer(Θ,oneplustanhΘ,oneminustanhΘ,ΔΘ,x,Wji,wjiSi_plus_WjiS²)
end

"""resets the relevant parameters of the buffer to initial values"""
function reset_GWF_buffer!(Buffer::RBMSpin1Buffer,ψG::RBMSpin1,S::AbstractMatrix)
    Buffer.x .= reshape(S,length(S))
    Buffer.Wji .= get_W_ij(ψG)'
    Buffer.wjiSi_plus_WjiS² = _precompute_PlaqFlipWeights!(Buffer.wjiSi_plus_WjiS²,get_w_ij(ψG),get_W_ij(ψG),S)
    return Buffer
end
function getParameterType(ψG::RBMSpin1,k)
    par = get_params(ψG)
    return _getParamsTypeAndIndex(par,k)[1]
end

guidingfunc_name(F::RBMSpin1) = "RBMSpin1"

# _get_Spin_normalization(x::StencilSpinConfig) = 1/(2*getSpin(x))
_get_Spin_normalization(x::StencilSpinConfig) = 1

Base.@propagate_inbounds function (ψG::RBMSpin1)(x::StencilSpinConfig)
    a = get_alpha_i(ψG)
    A = get_A_i(ψG)
    b = get_b_j(ψG)
    w = get_w_ij(ψG)
    W = get_W_ij(ψG)
    N = ψG.N
    @boundscheck checkbounds(x,1:N)
    @boundscheck checkbounds(w,1:N,axes(W,2))
    @boundscheck checkbounds(W,1:N,axes(W,2))


    aiσi = 0.

    for i in axes(W,1)
        xi = x[i]
        aiσi += a[i] * xi + A[i] * xi^2
    end
    # exp_term = exp(aiσi)
    
    # coshprod = 1.
    # cosh1inv = 1/cosh(1.)
    # logcosh1inv = log(1/cosh(1.))
    logcoshprod = 0.
    @inbounds for j in axes(W,2)
        Θ = get_Θ_j(x,j,b,w,W)
        # coshprod *= 2cosh(Θ)
        logcoshprod += log(cosh(Θ)) # get rid of factor 2 since it is just a global rescaling
    end
    return exp(aiσi + logcoshprod)
    # return (aiσi + logcoshprod)
end

Base.@propagate_inbounds function get_Θ_j(Config::StencilSpinConfig,j,b,w,W)
    θj = b[j]
    x = parent(Config)
    N = length(x)
    Is = CartesianIndices(x)

    @boundscheck checkbounds(x,1:N)
    @boundscheck checkbounds(W,1:N,j)
    @boundscheck checkbounds(w,1:N,j)

    @inbounds @simd for i in 1:N

        I = Is[i]
        xi = x[I]
        θj += (w[i,j] + W[i,j]*xi) * xi
        # θj += muladd(W[i,j],xi,w[i,j]) * xi
    end
    return θj
end

function compute_GWF_buffer!(Buffer::RBMSpin1Buffer,ψG::RBMSpin1,Walker::SpiderWebWalker)
    Config = get_config(Walker)
    Θ = Buffer.Θ
    b = get_b_j(ψG)
    x = Buffer.x
    x .= reshape(Config,length(Config))

    Is = CartesianIndices(x)
    W_ij = get_W_ij(ψG)
    w_ij = get_w_ij(ψG)
    # wji = Buffer.wji

    # mul!(Θ, wji, x)

    for j in eachindex(Θ)
        Θj = b[j]
        LoopVectorization.@turbo for i in eachindex(x)
            # Θj += w_ij[i,j]*x[i] + W_ij[i,j] * x²[i]
            # Θj += (w_ij[i,j] + W_ij[i,j]*x[i])*x[i]
            Θj += muladd(W_ij[i,j],x[i],w_ij[i,j])*x[i]

        end
        Θ[j] = Θj
    end
    Buffer.oneminustanhΘ .= tanh.(Θ)
    Buffer.oneplustanhΘ .= 1 .+ Buffer.oneminustanhΘ
    Buffer.oneminustanhΘ .= 1 .- Buffer.oneminustanhΘ
    
    return Buffer
end

function guidingfuncRatio_log(ψG::RBMSpin1,Walker::SpiderWebWalker,move::Tuple,Buffer::RBMSpin1Buffer)
    a = get_alpha_i(ψG)
    A = get_A_i(ψG)
    b = get_b_j(ψG)

    (;Θ,ΔΘ,oneplustanhΘ,oneminustanhΘ,Wji,x) = Buffer
    Conf = get_config(Walker)

    i,j,opSign = move

    sites = safe_parent_indices(parent(Conf), (i, j))

    LI = LinearIndices(Conf)

    exp_a = zero(eltype(Θ))
    exp_A = zero(eltype(Θ))
    # exp_a = 0.
    # exp_A = 0.

    @inbounds @simd for idx in eachindex(sites)
        i,j = sites[idx]
        s = P1_STENCIL[idx]*opSign
        I = LI[i,j]
        xi = Conf[I]
        exp_a += a[I]*s
        exp_A += A[I] * (2s*xi + s^2)
    end
    
    ΔΘ .= Buffer.wjiSi_plus_WjiS²[move] # save some time by using precomputed weights
    
    # fill!(ΔΘ,zero(eltype(ΔΘ)))

    for idx in eachindex(sites)
        I_x,I_y = sites[idx]
        i = LI[I_x,I_y]
        s = P1_STENCIL[idx]*opSign
        xi = x[i]
        fac = 2*s*xi
        fac == 0 && continue
        LoopVectorization.@turbo for j in eachindex(Θ)
            # ΔΘ[j] += wji[j,i]*s + Wji[j,i]*fac
            ΔΘ[j] += Wji[j,i]*fac
        end
    end

    coshprod = 1.
    LoopVectorization.@turbo for j in eachindex(Θ)
        Θj = Θ[j]
        oneplustanhΘj = oneplustanhΘ[j]
        oneminustanhΘj = oneminustanhΘ[j]
        
        ΔΘj = ΔΘ[j]

        eᶿ = exp(ΔΘj)
        e⁻ᶿ = inv(eᶿ)

        # coshprod *= 0.5*(eᶿ*(tanhΘj + 1) - e⁻ᶿ*(tanhΘj - 1))
        coshprod *= (eᶿ*oneplustanhΘj + e⁻ᶿ*oneminustanhΘj)
    end
    prefac = 1 / (2. ^ length(Θ))
    return exp_a + exp_A +log(coshprod*prefac)
end

function getOx_k(ψG::RBMSpin1,x::AbstractMatrix,k)
    paramtype = getParameterType(ψG,k)

    N = ψG.N
    M = N * ψG.hidden_density

    if paramtype == 1 
        return x[k]
    end
    k -= N
    if paramtype == 2
        return x[k]^2
    end
    k -= N

    bj = get_b_j(ψG)
    wij = get_w_ij(ψG)
    Wij = get_W_ij(ψG)
    
    if paramtype == 3
        j = k
        θj = get_Θ_j(x,j,bj,wij,Wij)
        return tanh(θj)
    end
    k -= M
    if paramtype == 4
        i,j = Tuple(CartesianIndices(Wij)[k])
        θj = get_Θ_j(x,j,bj,wij,Wij)
        return x[i] * tanh(θj)
    end
    k -= N*M
    if paramtype == 5
        i,j = Tuple(CartesianIndices(Wij)[k])
        θj = get_Θ_j(x,j,bj,wij,Wij)
        return x[i]^2 * tanh(θj)
    end
    return Ok
end
