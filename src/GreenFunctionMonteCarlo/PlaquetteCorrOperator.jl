"""Operator Δ + Σᵢⱼ cos(q rᵢⱼ) <PᵢPⱼ + PᵢPⱼ† + Pᵢ†Pⱼ + Pᵢ†Pⱼ†>.
The delta is to make sure the the operator is positive definite.
"""
struct BBqOperator_4 <: AbstractOperator 
end

operatorname(X::BBqOperator_4) = "BBqOperator_4"

function getAllTwoMoves(Walker::SpiderWebWalker)
    nx = Walker.n_x

    moves = copy(Walker.moves)
    # secondmoves = [empty(moves) for _ in eachindex(moves)]
    secondmoves = empty!([moves])

    Config = get_config(Walker)
    getMoves!(moves,Config)

    for move in moves
        i,j,s = move
        applyPlaquette!(Config, i, j, s)

        newmoves = getMoves!(copy(moves),Config)
        push!(secondmoves,newmoves)
        applyPlaquette!(Config, i, j, -s)
    end

    # return (;moves,secondmoves)
    return [(move1,move2) for (i,move1) in enumerate(moves) for move2 in secondmoves[i]]
end

function apply_operator!(Walker::SpiderWebWalker,O::BBqOperator_4,Guiding_function_buffer,ψG::AbstractGuidingFunction,q::SVector)
    AllMoves = getAllTwoMoves(Walker)
    if isempty(AllMoves)
        return 0.
    end
    # weights = getWeightListAll2Moves!(Walker,ψG,O,q)
    weights = getWeightListAll2Moves!(Walker,AllMoves,ψG,Guiding_function_buffer,O,q)

    w = sum(weights)
    if iszero(w)
        return w
    end
    moveidx = StatsBase.sample(StatsBase.Weights(weights))

    move1,move2 = AllMoves[moveidx]

    applyPlaquette!(Walker.Config, move1...)
    applyPlaquette!(Walker.Config, move2...)
    return w
end

function getWeightListAll2Moves!(Walker::SpiderWebWalker,AllMoves,ψG::AbstractGuidingFunction,Guiding_function_buffer,O::BBqOperator_4,q)
    nx = getNPlaq!(Walker)
    # (;weights) = Walker
    # empty!(weights)

    weights = [getWeight2Moves!(Walker,ψG,move1,move2,O,q) for (move1,move2) in AllMoves]
    # for I in plaquetteIterator(Config)
    #     appl = P_applicable(Config,I)

    #     for (operatorSign,valid) in zip((+1,-1),appl)
    #         if !valid 
    #             push!(weights,0.)
    #             continue
    #         end
    #         applyPlaquette!(Config, I[1], I[2], operatorSign)

    #         getMoves!(Walker)
    #         weights = updateWeightList!(Walker,Guiding_function_buffer,ψG)

    #         applyPlaquette!(Config, I[1], I[2], -operatorSign)
    #     end
    # end
    # push!(weights, (O.Delta + sum(nx))) # value of diagonal move
    # push!(weights, O.Delta) # value of diagonal move

    return weights
end

function getWeight2Moves!(Walker::SpiderWebWalker,ψG::AbstractGuidingFunction,move1,move2,O::BBqOperator_4,q)
    (;Config) = Walker

    ψx = ψG(Config)
    
    i_x, i_y, move_I = move1
    j_x, j_y, move_J = move2

    applyPlaquette!(Config, i_x, i_y, move_I)
    applyPlaquette!(Config, j_x, j_y, move_J)
    
    # n_x´ = getNPlaqfilled!(Walker,indices)

    # weight = guidingfuncRatio(ψG,n_x,n_x´,indices)
    ψx´ = ψG(Config)
    weight = ψx´/ψx
    # weight = 1

    qr_ij = q ⋅ SA[i_x-j_x, i_y-j_y]

    OperatorWeight = 2*(cos(qr_ij*0.5)^2)
    # OperatorWeight = 1

    applyPlaquette!(Config, j_x, j_y, -move_J)
    applyPlaquette!(Config, i_x, i_y, -move_I)

    return weight*OperatorWeight
end