
function initialize_GFMC(L,Spin = 1)
    Config = stencilConfig(zeros(Int8, L, L), Int8(Spin)) 
end


getOperatorRep(i,j,opNum) = i,j,opNum
# getOperatorRep(i,j,opNum) = CartesianIndex(i,j,opNum)
# getOperatorNumber(i,j,opNum,L) = LinearIndices((L,L,2))[i,j,opNum]

# function getOperatorFromNumber(Config, opNumber, operators)

# end

function getMoves!(
    moves,
    Conf
)
    empty!(moves)

    for I in plaquetteIterator(Conf)
        applPlus, applMinus = P_applicable(Conf, I)
        i,j = I
        if applMinus
            push!(moves, getOperatorRep(i,j,-1))
        end
        if applPlus
            push!(moves, getOperatorRep(i,j,1))
        end
    end
    return moves
end

import StatsBase
function performMarkovStep!(weights,moves,Config,weightfunc::T) where T
    empty!(weights)
    moves = getMoves!(moves,Config)
    ψₓ = weightfunc(Config)
    
    for operator in moves
        i,j,opNum = operator
        applyPlaquette!(Config, i, j, opNum)
        ψₓ´ = weightfunc(Config)
        push!(weights,ψₓ´/ψₓ)
        applyPlaquette!(Config, i, j, -opNum)
    end
    # TotalWeight = sum(weights)
    if isempty(weights)
        @info "No moves available" 
    end
    moveidx = StatsBase.sample(StatsBase.Weights(weights))
    move = moves[moveidx]
    applyPlaquette!(Config, move[1], move[2], move[3]) 
    return move,weights
end

function getLocalEnergy(weights)
    return -sum(1/weight for weight in weights)
end

function startSingleWalkerGFMC(InitialState,NSteps,weightfunc::T) where T
    # Config = initialize_GFMC(InitialState)
    # Config = copy(InitialState)
    Config = (InitialState)
    weights = Float64[]
    Allmoves = fill((Int8(0),Int8(0),Int8(0)),NSteps)
    LocalMoves = empty(Allmoves)
    energies = zeros(NSteps)
    TotalWeights = zeros(NSteps)
    
    for i in 1:NSteps
        move,weights = performMarkovStep!(weights,LocalMoves,Config,weightfunc)
        bx_TotalWeight = sum(weights)

        TotalWeights[i] = bx_TotalWeight
        
        Allmoves[i] = move
        energies[i] = getLocalEnergy(weights)
        
    end
    return (;Allmoves,energies,TotalWeights)
end

function normalizedAccWeight(weights,n,p)
    meanweight = mean(weights)
    # meanweight = 1
    prod(weights[n-j+1]/meanweight for j in 1:min(p,n))
end

function precomputeNormalizedAccWeight(weights,nThermal)
    meanweight = mean(weights)
    # meanweight = 1

    return [
        (cumprod( 
            weights[nThermal-i]/meanweight for i in lastindex(weights)-j:lastindex(weights)
        ))
        # for j in 1:length(weights)-1
        for j in 1:20
    ]
    # prod(weights[n-j]/meanweight for j in 0:p)
end
getAccumulatedWeight(weights,n) = cumprod(weights[n-j] for j in 1:n-1)

function getEnergy(weights,localEnergies,p,nthermalization)
    Gn(n) = normalizedAccWeight(weights,n,p)
    # EL_thermalized = @view localEnergies[nthermalization+1:end]
    N = lastindex(localEnergies)
    num = sum(Gn(n)*localEnergies[n] for n in nthermalization+1:N)
    denom = sum(Gn(n) for n in nthermalization+1:N)

    return num/denom
end