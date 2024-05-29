
struct VariationalGuidingFunction{A<:AbstractArray} <: AbstractGuidingFunction
    params::A
end

function constructVariationalFunction(State,α::Real=0.1)
    plaqs = collect(plaquetteIterator(State))
    N = length(plaqs)
    params = zeros(N,N+1)
    ψ = VariationalGuidingFunction(params)
    get_alpha_i(ψ) .= α
    return ψ
end

function get_alpha_i(ψG::VariationalGuidingFunction)
    return @view ψG.params[:,begin]
end

function get_beta_ij(ψG::VariationalGuidingFunction)
    return @view ψG.params[:,begin+1:end]
end

function (ψG::VariationalGuidingFunction)(N□::AbstractArray) 
    return exp(guidingfunc_exponent(ψG,N□))
end

function guidingfunc_exponent(ψG::VariationalGuidingFunction,N□::AbstractArray) 
    α = get_alpha_i(ψG)
    β = get_beta_ij(ψG)
    exponent = α' * N□ + dot(N□,β,N□)
    return exponent
end

guidingfuncRatio(ψG::VariationalGuidingFunction,n::AbstractArray,n´::AbstractArray,affectedPlaquettes) = exp(guidingfuncRatio_exponent(ψG,n,n´,affectedPlaquettes))
guidingfuncRatio(ψG::VariationalGuidingFunction,n::AbstractArray,n´::AbstractArray) = exp(guidingfuncRatio_exponent(ψG,n,n´))

function guidingfuncRatio_exponent(ψG::VariationalGuidingFunction,n::AbstractArray,n´::AbstractArray,affectedPlaquettes)
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

guidingfuncRatio_exponent(ψG::VariationalGuidingFunction,n::AbstractArray,n´::AbstractArray) = guidingfuncRatio_exponent(ψG,n,n´,eachindex(n,n´))


function getWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,weightfunc::VariationalGuidingFunction,Λ)
    (;Config,moves,weights) = Walker
    empty!(weights)

    n_x = getNPlaq!(Walker)
    
    for operator in moves
        i,j,opNum = operator
        indices = AffectedPlaquetteList[i,j]
        applyPlaquette!(Config, i, j, opNum)
         
        n_x´ = getNPlaqfilled!(Walker,indices)
        
        weight = guidingfuncRatio(weightfunc,n_x,n_x´,indices)
        # ndiff = sum(n_x´ .- n_x)
        # weight = weightfunc(ndiff)
        push!(weights,weight)
        applyPlaquette!(Config, i, j, -opNum)
    end
    if Λ != 0
        push!(moves, (0,0,0))
        push!(weights,Λ)
    end
    return weights
end
