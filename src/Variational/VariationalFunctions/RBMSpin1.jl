struct RBMSpin1_1{T} <: AbstractGuidingFunction
    a_i::Vector{T}
    A_i::Vector{T}
    b_j::Vector{T}
    w_ij::Matrix{T}
    W_ij::Matrix{T}
    N::Int
    hidden_density::Int
end

function RBMSpin1_1(N,hidden_density::Int,Type = Float32)
    M = N*hidden_density
    a_i = zeros(Type,N)
    A_i = zeros(Type,N)
    b_j = zeros(Type,M)
    w_ij = zeros(Type,N,M)
    W_ij = zeros(Type,N,M)

    ψ = RBMSpin1_1(a_i,A_i,b_j,w_ij,W_ij,N,hidden_density)
    return ψ
end
RBMSpin1_1(S::StencilSpinConfig,hidden_density::Int,Type = Float32) = RBMSpin1_1(length(S),hidden_density,Type)

get_alpha_i(ψG::RBMSpin1_1) = ψG.a_i
get_A_i(ψG::RBMSpin1_1) = ψG.A_i
get_b_j(ψG::RBMSpin1_1) = ψG.b_j
get_w_ij(ψG::RBMSpin1_1) = ψG.w_ij
get_W_ij(ψG::RBMSpin1_1) = ψG.W_ij

function get_params(ψG::RBMSpin1_1)
    return RecursiveArrayTools.ArrayPartition(ψG.a_i,ψG.A_i,ψG.b_j,ψG.w_ij,ψG.W_ij)
end

function getParameterType(ψG::RBMSpin1_1,k)
    N = ψG.N
    M = N * ψG.hidden_density
    par = get_params(ψG)

    for j in 1:length(par.x)
        k -= length(par.x[j])
        if k <= 0
            return j
        end
    end
end

guidingfunc_name(F::RBMSpin1_1) = "RBMSpin1"

Base.@propagate_inbounds function (ψG::RBMSpin1_1)(x::AbstractMatrix)
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
        aiσi += a[i] * x[i] + A[i] * x[i]^2
    end
    exp_term = exp(aiσi)

    coshprod = 1.
    cosh1inv = 1/cosh(1.)

    @inbounds for j in axes(W,2)
        Θ = get_Θ_j(x,j,b,w,W)
        coshprod *= 2cosh(Θ)
        # coshprod *= cosh(Θ)*cosh1inv # get rid of factor 2 since it is just a global rescaling
    end

    return exp_term * coshprod
end

Base.@propagate_inbounds function get_Θ_j(x::AbstractMatrix,j,b,w,W)
    θj = b[j]
    N = length(x)
    Is = CartesianIndices(x)

    @boundscheck checkbounds(x,1:N)
    @boundscheck checkbounds(W,1:N,j)
    @boundscheck checkbounds(w,1:N,j)

    @inbounds @simd for i in 1:N
        I = Is[i]
        xi = x[I]
        # θj += (w[i,j] + W[i,j]*xi) * xi
        θj += muladd(W[i,j],xi,w[i,j]) * xi
    end
    return θj
end

function getOx_k(ψG::RBMSpin1_1,x::AbstractMatrix,k)
    par = get_params(ψG)
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
