abstract type AbstractGuidingFunction end

struct RKFunction <: AbstractGuidingFunction end

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

    for ind in eachindex(affectedPlaquettes)
        i = affectedPlaquettes[ind]
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

guidingfuncRatio(ψG::AbstractGuidingFunction,n::AbstractArray,n´::AbstractArray,affectedPlaquettes) = exp(guidingfuncRatio_exponent(ψG,n,n´,affectedPlaquettes))
guidingfuncRatio(ψG::FullVariationalGuidingFunction,n::AbstractArray,n´::AbstractArray) = exp(guidingfuncRatio_exponent(ψG,n,n´))

function guidingfuncRatio_exponent(ψG::FullVariationalGuidingFunction,n::AbstractArray,n´::AbstractArray,affectedPlaquettes)
    α = get_alpha_i(ψG)
    β = get_beta_ij(ψG)

    exponent = zero(eltype(α))

    for i in affectedPlaquettes
        Δn = n´[i] - n[i]
        Δn == 0 && continue
        
        exp_i = α[i]
        LoopVectorization.@turbo for j in eachindex(n,n´)
            exp_i += β[j,i]*(n´[j] + n[j])
        end
        exponent += exp_i*Δn
    end

    return exponent
end

guidingfuncRatio_exponent(ψG::AbstractGuidingFunction,n::AbstractArray,n´::AbstractArray) = guidingfuncRatio_exponent(ψG,n,n´,eachindex(n,n´))

function updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,weightfunc::AbstractGuidingFunction,Λ=0)
    (;Config,weights) = Walker
    empty!(weights)
    moves = getMoves!(Walker)
    n_x = getNPlaq!(Walker)
    
    for operator in moves
        i,j,opNum = operator
        indices = AffectedPlaquetteList[i,j]
        applyPlaquette!(Config, i, j, opNum)
         
        n_x´ = getNPlaqfilled!(Walker,indices)
        
        weight = guidingfuncRatio(weightfunc,n_x,n_x´,indices)
        push!(weights,weight)
        applyPlaquette!(Config, i, j, -opNum)
    end
    if Λ != 0
        push!(moves, (0,0,0))
        push!(weights,Λ)
    end
    return weights
end

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
    
    for operator in moves
        i,j,opNum = operator
        indices = AffectedPlaquetteList[i,j]
        applyPlaquette!(Config, i, j, opNum)
        
        weight = length(moves)

        push!(weights,weight)
        applyPlaquette!(Config, i, j, -opNum)
    end
    if Λ != 0
        push!(moves, DIAGONAL_MOVE_ID)
        push!(weights,Λ)
    end
    return weights
end
