
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

function varitationalFunc(α,NPlaq::Integer,NPlaqEstimate)
    return exp(α*(NPlaq-NPlaqEstimate))
end

function ConstructVaritationalFunc(α,ConfEx=nothing)
    NPlaqEst = 0
    if ConfEx !== nothing
        NPlaqEst = NPlaquettes(ConfEx)
    end
    return Conf -> varitationalFunc(α,NPlaquettes(Conf),NPlaqEst)
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
    bn = @view weights[nThermal:end]
    # meanweight = mean(bn)
    # meanweight = mean(weights)
    meanweight = 1
    
    Gnp = zeros(length(bn)-PMax,PMax)
    # @views Gnp[:,1] .= bn[PMax+1:end] ./meanweight 
    for n in axes(Gnp,1)
        Gnp[n,1] = bn[n]/meanweight 
    end
    for p in 2:PMax
        # println("p = $p")
        # Threads.@threads 
        for n in axes(Gnp,1)[p:end]
            # Gnp[n,p] = Gnp[n,p-1]*bn[n-p]/meanweight
            Gnp[n,p] = Gnp[n,p-1]*Gnp[n-p+1,1]
        end
    end
    return Gnp
end

function getEnergy(weights,localEnergies,p,nthermalization)
    meanweight = mean(weights)
    num = 0.
    denom = 0.
    
    for n in nthermalization:length(localEnergies)
        Gnp = 1.
        for j in 1:p
            Gnp *= weights[n-j]/meanweight
        end
        num += Gnp*localEnergies[n]
        denom += Gnp
    end
    return num/denom
end

function getEnergies(weights,localEnergies,nthermalization,PMax)
    Gnp = precomputeNormalizedAccWeight(weights,nthermalization,PMax)
    EL_thermalized = @view localEnergies[nthermalization:end]
    N = lastindex(localEnergies)
    num = fetch.([Threads.@spawn sum(Gnp[n,p]*EL_thermalized[n-p] for n in axes(Gnp,1)[p+1:end]) for p in 1:PMax])
    denom = fetch.([Threads.@spawn sum(Gnp[n,p] for n in axes(Gnp,1)[p+1:end]) for p in 1:PMax])

    return num ./denom
end