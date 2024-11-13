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

struct RBMSpin1Buffer{T}
    Θ::Vector{T}
    coshΘ::Vector{T}
    ΔΘ::Vector{T}
    x²::Vector{T}
    wji::Matrix{T}
    Wji::Matrix{T}
end

function allocate_GWF_buffer(ψG::RBMSpin1,S::AbstractMatrix) 
    # Θ = similar(get_b_j(ψG))
    # ΔΘ = similar(get_b_j(ψG))
    # x² = zeros(eltype(get_params(ψG)),length(S))

    Θ = zeros(size(get_b_j(ψG)))
    coshΘ = zeros(size(get_b_j(ψG)))
    ΔΘ = zeros(size(get_b_j(ψG)))
    x² = zeros(length(S))

    wji = Array(get_w_ij(ψG)')
    Wji = Array(get_W_ij(ψG)')

    return RBMSpin1Buffer(Θ,coshΘ,ΔΘ,x²,wji,Wji)
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
    SpinNormalization = _get_Spin_normalization(x)
    # SpinNormalization = 1/2

    for i in axes(W,1)
        xi = x[i]*SpinNormalization
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

    SpinNormalization = _get_Spin_normalization(Config)
    @inbounds @simd for i in 1:N

        I = Is[i]
        xi = x[I] * SpinNormalization
        θj += (w[i,j] + W[i,j]*xi) * xi
        # θj += muladd(W[i,j],xi,w[i,j]) * xi
    end
    return θj
end

function fill_GWF_buffer!(Buffer::RBMSpin1Buffer,ψG::RBMSpin1,Walker::SpiderWebWalker)
    Config = get_config(Walker)
    Θ = Buffer.Θ
    b = get_b_j(ψG)
    x = reshape(Config,length(Config))
    x² = Buffer.x²

    Is = CartesianIndices(x)
    W_ij = get_W_ij(ψG)
    w_ij = get_w_ij(ψG)
    x² = Buffer.x²
    mul!(Θ, w_ij', x)

    x² .= x.*x 
    for j in eachindex(Θ)
        Θj_2 = zero(Θ[j])
        LoopVectorization.@turbo for i in eachindex(x)
            Θj_2 += W_ij[i,j] * x²[i]
        end
        Θ[j] = Θ[j] + Θj_2 +b[j]
    end
    Buffer.coshΘ .= cosh.(Θ)

    return Buffer
end

function guidingfuncRatio(ψG::RBMSpin1,Walker::SpiderWebWalker,move::Tuple,Buffer::RBMSpin1Buffer)
    a = get_alpha_i(ψG)
    A = get_A_i(ψG)
    b = get_b_j(ψG)

    (;Θ,ΔΘ,coshΘ,wji,Wji) = Buffer
    x = get_config(Walker)
    Mat = parent(x)

    i,j,opSign = move
    

    sites = safe_parent_indices(Mat, (i, j))

    LI = LinearIndices(Mat)

    exp_a = zero(eltype(Θ))
    exp_A = zero(eltype(Θ))
    # exp_a = 0.
    # exp_A = 0.

    @inbounds @simd for idx in eachindex(sites)
        i,j = sites[idx]
        s = P1_STENCIL[idx]*opSign
        I = LI[i,j]
        xi = x[I]
        exp_a += a[I]*s
        exp_A += A[I] * (2s*xi + s^2)
    end
    fill!(ΔΘ,zero(eltype(ΔΘ)))
    for idx in eachindex(sites)
        I_x,I_y = sites[idx]
        I = LI[I_x,I_y]
        s = P1_STENCIL[idx]*opSign
        xi = x[I]
        fac = (2*s*xi + s^2)
        LoopVectorization.@turbo for j in eachindex(Θ)
            ΔΘ[j] += wji[j,I]*s + Wji[j,I]*fac
        end
    end

    # logcoshprod = 0.
    coshprod = 1.
    LoopVectorization.@turbo for j in eachindex(Θ)
        Θj = Θ[j]
        coshΘj = coshΘ[j]
        ΔΘj = ΔΘ[j]
        # logcoshprod += log(cosh(Θj + ΔΘj)/cosh(Θj))
        # logcoshprod += log(cosh(Θj + ΔΘj)/coshΘj)
        coshprod *= cosh(Θj + ΔΘj)/coshΘj
    end
    # return exp(exp_a + exp_A + logcoshprod)
    return exp(exp_a + exp_A) *coshprod
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
