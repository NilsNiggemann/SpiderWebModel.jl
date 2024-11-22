"""provides a naive wrapper for a guiding function which does not use any buffer. Useful for debugging and testing"""
struct NaiveFunction{T<: AbstractGuidingFunction} <: AbstractGuidingFunction
    psi::T
end
(psi::NaiveFunction)(x::Any) = psi.psi(x)
(psi::NaiveFunction)(W::SpiderWebWalker) = psi.psi(W)

guidingfunc_name(F::NaiveFunction) = "NaiveFunction"
get_params(ψG::NaiveFunction) = get_params(ψG.psi)
guidingfuncRatio(ψG::NaiveFunction,W::SpiderWebWalker,move::Tuple{Int,Int,Int},Guiding_function_buffer) = guidingfuncRatio_naive(ψG.psi,W,move)
getOx_k(ψG::NaiveFunction,W::SpiderWebWalker,k) = getOx_k(ψG.psi,W,k)