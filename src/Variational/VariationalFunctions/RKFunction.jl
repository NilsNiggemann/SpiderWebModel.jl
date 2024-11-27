struct RKFunction <: AbstractGuidingFunction end
(psi::RKFunction)(W::SpiderWebWalker) = 1.
(psi::RKFunction)(x::Any) = 1.

guidingfunc_name(F::RKFunction) = "RKFunction"
get_params(ψG::RKFunction) = nothing
guidingfuncRatio_log(ψG::RKFunction,Walker::SpiderWebWalker,move,Buffer) = 0.
allocate_GWF_buffer(ψG::RKFunction,S::StencilSpinConfig) = nothing
# function updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::RKFunction,Λ=0)
#     (;Config,weights,moves,n_x,n_x´) = Walker
#     isempty(weights) && isempty(moves) || return weights #return if weights are already computed
    
#     getMoves!(Walker)
    
#     resize!(weights,length(moves))
#     weights .= 1.

#     if Λ != 0
#         push!(moves, DIAGONAL_MOVE_ID)
#         push!(weights,Λ)
#     end
#     return weights
# end