"""provides a naive wrapper for a guiding function which does not use any buffer. Useful for debugging and testing"""
struct NaiveFunction{T<: AbstractGuidingFunction} <: AbstractGuidingFunction
    psi::T
end
(psi::NaiveFunction)(W::SpiderWebWalker) = psi.psi(W)

guidingfunc_name(F::NaiveFunction) = "NaiveFunction"
get_params(ψG::NaiveFunction) = get_params(ψG.psi)

getOx_k(ψG::NaiveFunction,W::SpiderWebWalker,k) = getOx_k(ψG.psi,W,k)