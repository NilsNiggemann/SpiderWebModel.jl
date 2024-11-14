struct PlaquetteRBM{T} <: AbstractGuidingFunction
    a_i::Vector{T}
    A_i::Vector{T}
    b_j::Vector{T}
    w_ij::Matrix{T}
    W_ij::Matrix{T}
    N::Int
    hidden_density::Int
end

function PlaquetteRBM(N,hidden_density::Int,Type = Float32)
    M = N*hidden_density
    a_i = zeros(Type,N)
    A_i = zeros(Type,N)
    b_j = zeros(Type,M)
    w_ij = zeros(Type,N,M)
    W_ij = zeros(Type,N,M)

    ψ = PlaquetteRBM(a_i,A_i,b_j,w_ij,W_ij,N,hidden_density)
    return ψ
end
PlaquetteRBM(S::StencilSpinConfig,hidden_density::Int,Type = Float32) = PlaquetteRBM(length(collect(plaquetteIterator(S))),hidden_density,Type)

get_alpha_i(ψG::PlaquetteRBM) = ψG.a_i
get_A_i(ψG::PlaquetteRBM) = ψG.A_i
get_b_j(ψG::PlaquetteRBM) = ψG.b_j
get_w_ij(ψG::PlaquetteRBM) = ψG.w_ij
get_W_ij(ψG::PlaquetteRBM) = ψG.W_ij

function get_params(ψG::PlaquetteRBM)
    return RecursiveArrayTools.ArrayPartition(ψG.a_i,ψG.A_i,ψG.b_j,ψG.w_ij,ψG.W_ij)
end

function allocate_GWF_buffer(ψG::PlaquetteRBM,S::AbstractMatrix) 
    AffectedPlaquetteList = precomputeAffectedPlaquettes(S)

    Θ = zeros(size(get_b_j(ψG)))
    ΔΘ = zeros(size(get_b_j(ψG)))
    n² = zeros(ψG.N)

    return (;Θ,ΔΘ,n²,AffectedPlaquetteList)
end

function getParameterType(ψG::PlaquetteRBM,k)
    par = get_params(ψG)
    return _getParamsTypeAndIndex(par,k)[1]
end

guidingfunc_name(F::PlaquetteRBM) = "PlaquetteRBM"

# _get_Spin_normalization(x::StencilSpinConfig) = 1/(2*getSpin(x))

Base.@propagate_inbounds function _evaluate_plaquetteRBM(ψG::PlaquetteRBM,n::AbstractVector)
    a = get_alpha_i(ψG)
    A = get_A_i(ψG)
    b = get_b_j(ψG)
    w = get_w_ij(ψG)
    W = get_W_ij(ψG)
    N = ψG.N

    @boundscheck checkbounds(n,1:N)
    @boundscheck checkbounds(w,1:N,axes(w,2))

    αini = 0.

    for i in axes(w,1)
        αini += a[i] * n[i] + A[i] * n[i]^2
    end

    logthetasum = 0.
    @inbounds for j in axes(w,2)
        Θ = get_θ_j_plaq(n,j,b,w,W)
        logthetasum += log(1 + exp(Θ))
    end
    return exp(αini + logthetasum)
end
Base.@propagate_inbounds (ψG::PlaquetteRBM)(W::SpiderWebWalker) = _evaluate_plaquetteRBM(ψG,W.n_x)
Base.@propagate_inbounds function (ψG::PlaquetteRBM)(S::StencilSpinConfig) 
    n = getNPlaq(S)
    _evaluate_plaquetteRBM(ψG,n)
end

Base.@propagate_inbounds function get_θ_j_plaq(n::AbstractVector,j::Integer,b::AbstractVector,w::AbstractMatrix,W::AbstractMatrix)
    θj = b[j]
    @boundscheck checkbounds(w,firstindex(w),j)

    @inbounds @simd for i in axes(w,1)
    # for i in axes(w,1)
        θj += w[i,j] * n[i] + W[i,j] * n[i]^2
    end
    return θj
end

# __getTheta(W,ψG) = [get_θ_j_plaq(W.n_x,j,get_b_j(ψG),get_w_ij(ψG)) for j in 1:length(get_b_j(ψG))]


function compute_GWF_buffer!(Buffer,ψG::PlaquetteRBM,Walker::SpiderWebWalker)
    Θ = Buffer.Θ
    b = get_b_j(ψG)
    w_ij = get_w_ij(ψG)
    W_ij = get_W_ij(ψG)
    getNPlaq!(Walker)
    n = Walker.n_x
    n² = Buffer.n²
    mul!(Θ, w_ij', n)
    # Θ .+= b
    n² .= n.*n

    @boundscheck checkbounds(W_ij,eachindex(n²),eachindex(Θ))
    for j in eachindex(Θ)
        Θj_2 = zero(Θ[j])
        LoopVectorization.@turbo for i in eachindex(n²)
            Θj_2 += W_ij[i,j] * n²[i]
        end
        Θ[j] = Θ[j] + Θj_2 +b[j]
    end
    return Buffer
end

function guidingfuncRatio(ψG::PlaquetteRBM,Walker::SpiderWebWalker,move,Buffer)

    a = get_alpha_i(ψG)
    A = get_A_i(ψG)
    w = get_w_ij(ψG)
    W = get_W_ij(ψG)

    n = Walker.n_x
    n´ = Walker.n_x´

    (;Θ,ΔΘ,AffectedPlaquetteList) = Buffer

    i,j,opSign = move

    affectedPlaquettes = AffectedPlaquetteList[i,j]

    Config = get_config(Walker)

    applyPlaquette!(Config, i, j, opSign)
    getNPlaqfilled!(Walker,affectedPlaquettes)
    applyPlaquette!(Config, i, j, -opSign)

    exp_a = zero(eltype(Θ))
    exp_A = zero(eltype(Θ))

    # @boundscheck checkbounds(n,eachindex(affectedPlaquettes))
    # @boundscheck checkbounds(w,eachindex(affectedPlaquettes),eachindex(Θ))
    # @boundscheck checkbounds(w,i for i in affectedPlaquettes,eachindex(Θ))

    @inbounds @simd for ind in eachindex(affectedPlaquettes)
        i = affectedPlaquettes[ind]
        Δn = n´[i] - n[i]
        Δn² = n´[i]^2 - n[i]^2

        exp_a += a[i]* Δn
        exp_A += A[i] *Δn²
    end

    @inbounds @simd for j in eachindex(Θ)
    # for j in eachindex(Θ)
        ΔΘj = zero(eltype(Θ))
        # for idx in eachindex(sites)
        for ind in eachindex(affectedPlaquettes)
            i = affectedPlaquettes[ind]
            Δn = n´[i] - n[i]
            Δn² = n´[i]^2 - n[i]^2

            ΔΘj += w[i,j]*Δn + W[i,j]*Δn²
        end
        ΔΘ[j] = ΔΘj
    end

    logthetasum = 0.
    # for j in eachindex(Θ)
    LoopVectorization.@turbo for j in eachindex(Θ)
        Θj = Θ[j]
        ΔΘj = ΔΘ[j]
        arg = (1 + exp(Θj + ΔΘj))/(1 + exp(Θj))
        logthetasum += log(arg)
    end

    return exp(exp_a + exp_A + logthetasum)
end

function getOx_k(ψG::PlaquetteRBM,Walker::SpiderWebWalker,k)
    par = get_params(ψG)
    kold = k
    # paramtype,k = _getParamsTypeAndIndex(par,k)
    paramtype = getParameterType(ψG,k)
    getNPlaq!(Walker)
    n = Walker.n_x

    N = ψG.N
    M = N * ψG.hidden_density

    if paramtype == 1 
        return n[k]
    end

    k -= N

    if paramtype == 2
        return n[k]^2
    end
    k -= N

    bj = get_b_j(ψG)
    wij = get_w_ij(ψG)
    Wij = get_W_ij(ψG)
    
    if paramtype == 3
        j = k
        θj = get_θ_j_plaq(n,j,bj,wij,Wij)
        return 1/(1+exp(-θj))
    end
    k -= M

    if paramtype == 4
        i,j = Tuple(CartesianIndices(wij)[k])
        θj = get_θ_j_plaq(n,j,bj,wij,Wij)
        return n[i]/(1+exp(-θj))
    end
    k -= N*M

    if paramtype == 5
        i,j = Tuple(CartesianIndices(wij)[k])
        θj = get_θ_j_plaq(n,j,bj,wij,Wij)
        return n[i]^2/(1+exp(-θj))
    end

end
