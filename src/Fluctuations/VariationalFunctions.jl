abstract type AbstractGuidingFunction end

struct RKFunction <: AbstractGuidingFunction end

(ψG::AbstractGuidingFunction)(x::AbstractWalker) = ψG(get_config(x))
struct PlaquetteNumberGuidingFunction <: AbstractGuidingFunction
    α::Float64
end
(ψG::PlaquetteNumberGuidingFunction)(ΔNPlaq::Real) = exp(ψG.α*ΔNPlaq)

guidingfunc_name(F::Function) = string(typeof(F))
guidingfunc_name(F::AbstractGuidingFunction) = string(typeof(F))
guidingfunc_name(F::RKFunction) = "RKFunction"
guidingfunc_name(F::PlaquetteNumberGuidingFunction) = "PlaquetteNumberGuidingFunction"

variational_parameters(P::PlaquetteNumberGuidingFunction) = Dict([:alpha=>P.α])
variational_parameters(::RKFunction) = Dict{Symbol,Float64}()
variational_parameters(P::AbstractGuidingFunction) = Dict([:params=>P.params])
variational_parameters(P::Function) = Dict([x => getproperty(P,x) for x in propertynames(P)])

guidingfuncRatio_exponent(ψG::RKFunction,n::AbstractArray,n´::AbstractArray,affectedPlaquettes) = 0.

function guidingfuncRatio_exponent(ψG::PlaquetteNumberGuidingFunction,n::AbstractArray,n´::AbstractArray,affectedPlaquettes)
    α = ψG.α
    exponent = zero(α)

    for i in affectedPlaquettes
        Δn = n´[i] - n[i]
        exponent += Δn
    end

    return exponent * α
end

struct FullVariationalGuidingFunction{A<:AbstractArray} <: AbstractGuidingFunction
    params::A
end
guidingfunc_name(F::FullVariationalGuidingFunction) = "FullVariationalGuidingFunction"

function fullVariationalFunction(State,α::Real=0.1)
    plaqs = collect(plaquetteIterator(State))
    N = length(plaqs)
    params = zeros(N,N+1)
    ψ = FullVariationalGuidingFunction(params)
    get_alpha_i(ψ) .= α
    return ψ
end

function get_alpha_i(ψG::FullVariationalGuidingFunction)
    return @view ψG.params[:,begin]
end

function get_beta_ij(ψG::FullVariationalGuidingFunction)
    return @view ψG.params[:,begin+1:end]
end

function (ψG::FullVariationalGuidingFunction)(N□::AbstractArray) 
    return exp(guidingfunc_exponent(ψG,N□))
end

function guidingfunc_exponent(ψG::FullVariationalGuidingFunction,N□::AbstractArray) 
    α = get_alpha_i(ψG)
    β = get_beta_ij(ψG)
    exponent = α' * N□ + dot(N□,β,N□)
    return exponent
end

@inline guidingfuncRatio(ψG::AbstractGuidingFunction,n::AbstractArray,n´::AbstractArray,affectedPlaquettes) = exp(guidingfuncRatio_exponent(ψG,n,n´,affectedPlaquettes))
@inline guidingfuncRatio(ψG::FullVariationalGuidingFunction,n::AbstractArray,n´::AbstractArray) = exp(guidingfuncRatio_exponent(ψG,n,n´))

function guidingfuncRatio_exponent(ψG::FullVariationalGuidingFunction,n::AbstractArray,n´::AbstractArray,affectedPlaquettes)
    α = get_alpha_i(ψG)
    β = get_beta_ij(ψG)

    exponent = zero(eltype(α))

    for i in affectedPlaquettes
        Δn = n´[i] - n[i]
        iszero(Δn) && continue
        
        exp_i = α[i]
        LoopVectorization.@turbo for j in eachindex(n,n´)
            exp_i += β[j,i]*(n´[j] + n[j])
        end
        exponent += exp_i*Δn
    end
    # exponent1 = α' *  n + dot(n,β,n)
    # exponent2 = α' *  n´ + dot(n´,β,n´)
    # exponent = exponent2 - exponent1

    return exponent
end

guidingfuncRatio_exponent(ψG::AbstractGuidingFunction,n::AbstractArray,n´::AbstractArray) = guidingfuncRatio_exponent(ψG,n,n´,eachindex(n,n´))

function updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::AbstractGuidingFunction,Λ=0)
    (;Config,weights) = Walker
    empty!(weights)
    moves = getMoves!(Walker)
    ψx = ψG(Walker)
    for operator in moves
        i,j,opNum = operator
        applyPlaquette!(Config, i, j, opNum)
        
        ψx´ = ψG(Walker)
        
        weight = ψx´/ψx
        push!(weights,weight)
        applyPlaquette!(Config, i, j, -opNum)
    end
    if Λ != 0
        push!(moves, (0,0,0))
        push!(weights,Λ)
    end
    return weights
end

function updateWeightList_plaqs!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::AbstractGuidingFunction,Λ=0)
    (;Config,weights) = Walker
    empty!(weights)
    moves = getMoves!(Walker)
    n_x = getNPlaq!(Walker)
    
    for operator in moves
        i,j,opNum = operator
        indices = AffectedPlaquetteList[i,j]
        applyPlaquette!(Config, i, j, opNum)
        
        n_x´ = getNPlaqfilled!(Walker,indices)
        
        weight = guidingfuncRatio(ψG,n_x,n_x´,indices)
        push!(weights,weight)
        applyPlaquette!(Config, i, j, -opNum)
    end
    if Λ != 0
        push!(moves, (0,0,0))
        push!(weights,Λ)
    end
    return weights
end

@inline updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::FullVariationalGuidingFunction,Λ=0) = updateWeightList_plaqs!(Walker,AffectedPlaquetteList,ψG,Λ)
@inline updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::PlaquetteNumberGuidingFunction,Λ=0) = updateWeightList_plaqs!(Walker,AffectedPlaquetteList,ψG,Λ)

struct LocalPlaquetteGuidingFunction{A<:AbstractVector} <: AbstractGuidingFunction
    params::A
end
guidingfunc_name(F::LocalPlaquetteGuidingFunction) = "LocalPlaquetteGuidingFunction"

function localPlaquetteGuidingFunction(State,α::Real=0.1)
    plaqs = collect(plaquetteIterator(State))
    N = length(plaqs)
    params = zeros(Float32,N)
    ψ = LocalPlaquetteGuidingFunction(params)

    get_alpha_i(ψ) .= α
    return ψ
end
@inline updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::LocalPlaquetteGuidingFunction,Λ=0) = updateWeightList_plaqs!(Walker,AffectedPlaquetteList,ψG,Λ)

function get_alpha_i(ψG::LocalPlaquetteGuidingFunction)
    return ψG.params
end

(ψG::LocalPlaquetteGuidingFunction)(N□::AbstractArray) = exp(guidingfunc_exponent(ψG,N□))

function guidingfunc_exponent(ψG::LocalPlaquetteGuidingFunction,N□::AbstractArray) 
    α = get_alpha_i(ψG)
    exponent = α' * N□
    return exponent
end

function guidingfuncRatio_exponent(ψG::LocalPlaquetteGuidingFunction,n::AbstractArray,n´::AbstractArray,affectedPlaquettes)
    α = get_alpha_i(ψG)

    exponent = zero(eltype(α))

    for ind in eachindex(affectedPlaquettes)
        i = affectedPlaquettes[ind]
        Δn = n´[i] - n[i]
        exponent += α[i]*Δn
    end

    return exponent
end

function updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::RKFunction,Λ=0)
    (;Config,weights,moves,n_x,n_x´) = Walker
    isempty(weights) && isempty(moves) || return weights #return if weights are already computed
    
    getMoves!(Walker)
    
    resize!(weights,length(moves))
    weights .= 1.

    if Λ != 0
        push!(moves, DIAGONAL_MOVE_ID)
        push!(weights,Λ)
    end
    return weights
end

function updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::T,Λ=0) where T
    (;Config,weights,moves,n_x,n_x´) = Walker
    isempty(weights) && isempty(moves) || return weights #return if weights are already computed
    
    getMoves!(Walker)
    # getNPlaq!(Walker)
    
    psix_inv = 1/ψG(get_config(Walker))
    
    for operator in moves
        i,j,opNum = operator
        indices = AffectedPlaquetteList[i,j]
        applyPlaquette!(Config, i, j, opNum)
        
        # getNPlaq!(Walker,indices)
        
        # N□ = getNPlaq_difference(n_x,n_x´,indices) 
        weight = ψG(get_config(Walker)) * psix_inv

        push!(weights,weight)
        applyPlaquette!(Config, i, j, -opNum)
    end

    
    if Λ != 0
        push!(moves, DIAGONAL_MOVE_ID)
        push!(weights,Λ)
    end
    return weights
end

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