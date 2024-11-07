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

function allocate_GWF_buffer(ψG::RBMSpin1,S::AbstractMatrix) 
    Θ = similar(get_b_j(ψG))
    ΔΘ = similar(get_b_j(ψG))
    x² = zeros(eltype(get_params(ψG)),length(S))
    return (;Θ,ΔΘ,x²)
end

function getParameterType(ψG::RBMSpin1,k)
    par = get_params(ψG)
    return _getParamsTypeAndIndex(par,k)[1]
end

guidingfunc_name(F::RBMSpin1) = "RBMSpin1"

Base.@propagate_inbounds function (ψG::RBMSpin1)(x::AbstractMatrix)
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

function fill_GWF_buffer!(Buffer,ψG::RBMSpin1,Config::AbstractMatrix)
    Θ = Buffer.Θ
    x = reshape(Config,length(Config))
    x² = Buffer.x²

    Is = CartesianIndices(x)
    W_ij = get_W_ij(ψG)
    w_ij = get_w_ij(ψG)
    x² = Buffer.x²
    mul!(Θ, w_ij, x)
    for j in eachindex(Θ)
        Θj = Θ[j]
        for i in eachindex(x)
            Θj += W_ij[i,j] * x²[i]
        end
        Θ[j] = Θj
    end
    return Θ

    x² .= x.*x
end

function _guidingfuncRatioRBM(ψG::RBMSpin1,Walker::SpiderWebWalker,move::Tuple,Buffer)
    a = get_alpha_i(ψG)
    A = get_A_i(ψG)
    b = get_b_j(ψG)
    w = get_w_ij(ψG)
    W = get_W_ij(ψG)

    (;Θ,x²,ΔΘ) = Buffer
    
    x = get_config(Walker)
    Mat = parent(x)

    i,j,opSign = move
    
    sites = safe_parent_indices(Mat, (i, j))

    LI = LinearIndices(Mat)

    exp_m = zero(eltype(Θ))

    @inbounds @simd for idx in eachindex(sites)
        i,j = sites[idx]
        s = P1_STENCIL[idx]*opSign
        I = LI[i,j]
        exp_m += a[I]*s + A[I] * (2s*x[I] + s^2)
    end

    for j in eachindex(Θ)
        ΔΘj = zero(eltype(Θ))
        for idx in eachindex(sites)
        # @inbounds @simd for idx in eachindex(sites)
            I_x,I_y = sites[idx]
            I = LI[I_x,I_y]
            s = P1_STENCIL[idx]*opSign
            ΔΘj += w[I,j]*s + W[I,j]*(2*s*x[I] + s^2)
        end
        ΔΘ[j] = ΔΘj
    end
    coshprod = 1.

    for j in eachindex(Θ)
        Θj = Θ[j]
        ΔΘj = ΔΘ[j]
        coshprod *= cosh(Θj + ΔΘj)/cosh(Θj)
    end
    return exp(exp_m)* coshprod
end

guidingfuncRatio(ψG::RBMSpin1,Walker::SpiderWebWalker,move,Buffer) = _guidingfuncRatioRBM(ψG,Walker,move,Buffer) 

# function updateWeightList!(Walker::SpiderWebWalker,Buffer,ψG::RBMSpin1,Λ=0)
#     (;Config,weights) = Walker
#     empty!(weights)
#     moves = getMoves!(Walker)
    
#     x = parent(parent(Config))
#     # x = parent(parent(get_config(Walker)))
#     fill_GWF_buffer!(Buffer,ψG,x)

#     for operator in moves
#         i,j,opNum = operator
#         weight = guidingfuncRatio(ψG,Walker,operator,Buffer)
#         push!(weights,weight)
#     end
#     if Λ != 0
#         push!(moves, (0,0,0))
#         push!(weights,Λ)
#     end
#     return weights
# end


function getOx_k(ψG::RBMSpin1,x::AbstractMatrix,k)
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
