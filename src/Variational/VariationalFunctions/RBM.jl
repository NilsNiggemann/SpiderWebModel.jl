struct RBM{A<:AbstractArray} <: AbstractGuidingFunction
    params::A
    N::Int
    hidden_density::Int
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

function get_alpha_i(ψG::RBM)
    N = ψG.N
    return @view ψG.params[1:N]
end

function get_b_j(ψG::RBM)
    N = ψG.N
    M = N * ψG.hidden_density
    bparams = @view ψG.params[N+1:N+M,:]
    return reshape(bparams,M)
end

function get_W_ij(ψG::RBM)
    N = ψG.N
    M = N * ψG.hidden_density
    Wparams = @view ψG.params[N+M+1:end]
    return reshape(Wparams,N,M)
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
    for i in axes(W,1)
        θj += W[i,j]*x[i]
    end
    return θj
end

