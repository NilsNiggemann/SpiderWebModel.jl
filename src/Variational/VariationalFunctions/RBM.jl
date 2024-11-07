
struct RBM{T<:Number} <: AbstractGuidingFunction
    α_i::Vector{T}
    b_j::Vector{T}
    W_ij::Matrix{T}
    N::Int
    hidden_density::Int
end
function RBM(N::Int,hidden_density::Int,Type = Float32)
    α_i = zeros(Type,N)
    M = N*hidden_density
    b_j = zeros(Type,M)
    W_ij = zeros(Type,N,M)

    # return (α_i,b_j,W_ij,N,hidden_density)
    return RBM(α_i,b_j,W_ij,N,hidden_density)
end
RBM(S::StencilSpinConfig,hidden_density::Int,Type = Float32) = RBM(length(S),hidden_density,Type)

get_alpha_i(ψG::RBM) = ψG.α_i
get_b_j(ψG::RBM) = ψG.b_j
get_W_ij(ψG::RBM) = ψG.W_ij
function get_params(ψG::RBM)
    return RecursiveArrayTools.ArrayPartition(ψG.α_i,ψG.b_j,ψG.W_ij)
end

function RBM(x::AbstractMatrix,hidden_density::Int)
    N = length(x)
    M = N*hidden_density
    NParams = N + M + N*M

    params = zeros(NParams)
    # params = zeros(Float32,NParams)
    ψ = RBM(params,N,hidden_density)
    return ψ
end

function getParameterType(ψG::RBM,k)
    N = ψG.N
    M = N * ψG.hidden_density
    if k <= N
        return 1 # α
    elseif k <= N+M
        return 2 # b
    else
        return 3 # W
    end
end

guidingfunc_name(F::RBM) = "RBM"
function (ψG::RBM)(x::AbstractMatrix)
    a = get_alpha_i(ψG)
    b = get_b_j(ψG)
    W = get_W_ij(ψG)

    xlist = reshape(x, length(x))
    aiσi = a' * xlist
    exp_term = exp(aiσi)

    coshprod = 1.
    # cosh1inv = 1/cosh(1.)

    for j in axes(W,2)
        θj = get_theta_j(x,j,b,W)
        # θj = 1.
        coshprod *= 2cosh(θj)
        # coshprod *= cosh(θj)*cosh1inv # get rid of factor 2 since it is just a global rescaling
    end

    return exp_term * coshprod
end

function get_theta_j(x::AbstractMatrix,j,b,W)
    θj = b[j]
    N = length(x)
    # @inbounds @simd for i in axes(W,1)
    # for i in eachindex(IndexLinear(),x)
    # LoopVectorization.@turbo for i in axes(W,1)
    # for i in 1:N
    @inbounds @simd for i in 1:N
        θj += W[i,j]*x[i]
    end
    return θj
end

function getOx_k(ψG::RBM,x::AbstractMatrix,k)
    par = get_params(ψG)

    paramtype = getParameterType(ψG,k)

    if paramtype == 1 
        return x[k]
    end
    bj = get_b_j(ψG)
    Wij = get_W_ij(ψG)
    
    if paramtype == 2
        j = k - ψG.N
        θj = get_theta_j(x,j,bj,Wij)
        return tanh(θj)
    elseif paramtype == 3
        i,j = Tuple(CartesianIndices(Wij)[k-ψG.N-length(bj)])
        θj = get_theta_j(x,j,bj,Wij)
        return x[i] * tanh(θj)
    end
    return Ok
end
