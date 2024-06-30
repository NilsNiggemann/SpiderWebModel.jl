abstract type AbstractDiagonalOperator end
(Hxx::AbstractDiagonalOperator)(x::SpiderWebWalker) = Hxx(get_config(x))

struct Hxx_zero <: AbstractDiagonalOperator end
(Hxx::Hxx_zero)(x::StencilSpinConfig) = 0.

struct Hxx_RK <: AbstractDiagonalOperator
    μ::Float64
end
(Hxx::Hxx_RK)(Walker::SpiderWebWalker) = Hxx.μ * length(Walker.moves)
function (Hxx::Hxx_RK)(x::AbstractMatrix)
    n = 0
    for I in plaquetteIterator(x)
        @inbounds applPlus, applMinus = P_applicable(Config, I)
        n += applPlus + applMinus
    end
    return n
end

struct Hxx_SIA <: AbstractDiagonalOperator
    U::Float64
end
(Hxx::Hxx_SIA)(Walker::SpiderWebWalker) = Hxx.U * sum(abs2, get_config(Walker))

struct CombinedOperator{F1,F2} <: AbstractDiagonalOperator
    Hxx1::F1
    Hxx2::F2
end
Base.:+(Hxx1::AbstractDiagonalOperator,Hxx2::AbstractDiagonalOperator) = CombinedOperator(Hxx1,Hxx2)

(Hxx::CombinedOperator)(Walker::SpiderWebWalker) = Hxx.Hxx1(Walker) + Hxx.Hxx2(Walker)

struct PlaquetteNumberOperator_1 <: 
    AbstractDiagonalOperator 
    I::CartesianIndex{2}
end
PlaquetteNumberOperator_1((i,j)) = PlaquetteNumberOperator_1(CartesianIndex(i,j))
function (N::PlaquetteNumberOperator_1)(x::StencilSpinConfig)
    applPlus, applMinus = P_applicable(x, N.I)
    n = applPlus + applMinus
    return n
end
