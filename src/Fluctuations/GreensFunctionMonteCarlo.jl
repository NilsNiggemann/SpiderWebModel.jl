
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
function performMarkovStep!(weights,moves,Config,weightfunc::T,Λ) where T
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
    push!(moves, (0,0,0))
    push!(weights,Λ)
    # TotalWeight = sum(weights)
    if isempty(weights)
        @info "No moves available" 
    end
    moveidx = StatsBase.sample(StatsBase.Weights(weights))
    move = moves[moveidx]
    if move != (0,0,0)
        applyPlaquette!(Config, move[1], move[2], move[3]) 
    end
    return move,weights
end

function getLocalEnergy(weights,Λ)
    return -sum(weights) + Λ
end


function NPlaquettes(Conf)
    moves = 0
    for I in plaquetteIterator(Conf)
        applPlus, applMinus = P_applicable(Conf, I)
        moves += applPlus + applMinus
    end
    return moves
end

function VaritationalFunc(α,NPlaq::Integer)
    return exp(α*NPlaq-20)
end

function VaritationalFunc(α)
    return Conf -> VaritationalFunc(α,NPlaquettes(Conf))
end

function startSingleWalkerGFMC(InitialState,NSteps,weightfunc::T,Λ) where T
    # Config = initialize_GFMC(InitialState)
    # Config = copy(InitialState)
    Config = deepcopy(InitialState)
    weights = Float64[]
    Allmoves = fill((Int8(0),Int8(0),Int8(0)),NSteps)
    LocalMoves = empty(Allmoves)
    energies = zeros(NSteps)
    TotalWeights = zeros(NSteps)
    
    for i in 1:NSteps
        move,weights = performMarkovStep!(weights,LocalMoves,Config,weightfunc,Λ)
        bx_TotalWeight = sum(weights)

        TotalWeights[i] = bx_TotalWeight
        
        Allmoves[i] = move
        energies[i] = getLocalEnergy(weights,Λ)
        
    end
    return (;Allmoves,energies,TotalWeights)
end

function normalizedAccWeight(weights,n,p)
    meanweight = mean(weights)
    # meanweight = 1
    prod(weights[n-j]/meanweight for j in 1:p)
end

function precomputeNormalizedAccWeight(weights,nThermal,PMax)
    # meanweight = 1
    N = lastindex(weights)
    bn = @view weights[nThermal:end]
    # meanweight = mean(bn)
    meanweight = mean(weights)
    
    Gnp = zeros(length(bn)-1,PMax)
    Gnp[:,1] .= bn[begin+1:end] ./meanweight 
    
    for p in 2:PMax
        # println("p = $p")
        Threads.@threads for n in p+1:lastindex(bn)-1
            # Gnp[n,p] = Gnp[n,p-1]*bn[n-p]/meanweight
            Gnp[n,p] = Gnp[n,p-1]*Gnp[n-p,1]
        end
    end
    return Gnp
end

function getEnergy(weights,localEnergies,p,nthermalization)
    Gn(n) = normalizedAccWeight(weights,n,p)
    # EL_thermalized = @view localEnergies[nthermalization+1:end]
    N = lastindex(localEnergies)
    num = sum(Gn(n)*localEnergies[n] for n in nthermalization+1:N)
    denom = sum(Gn(n) for n in nthermalization+1:N)

    return num/denom
end

function getEnergies(weights,localEnergies,nthermalization,PMax)
    Gnp = precomputeNormalizedAccWeight(weights,nthermalization,PMax)
    EL_thermalized = @view localEnergies[nthermalization+1:end]
    N = lastindex(localEnergies)
    num = fetch.([Threads.@spawn sum(Gnp[n,p]*EL_thermalized[n] for n in eachindex(EL_thermalized)) for p in 1:PMax])
    denom = fetch.([Threads.@spawn sum(Gnp[n,p] for n in eachindex(EL_thermalized)) for p in 1:PMax])

    return num ./denom
end