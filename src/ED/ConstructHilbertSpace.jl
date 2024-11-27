function spinConfig!(Conf, path::AbstractVector, InitialConf, plaqMap::PlaqMapping)
    Conf .= InitialConf
    for (i, op) in enumerate(path) #this can probably be made faster
        if op
            ij = plaqMap(i)
            flipPlaquette!(Conf, ij)
        end
    end
    return Conf
end
function spinConfig(path, InitialConf, plaqMap)
    spinConfig!(copy(InitialConf), path, InitialConf, plaqMap)
end

function getNewStates!(states, Conf, StateRep::SBitVector, plaqMap::PlaqMapping)
    empty!(states)
    @inline function convertToStateRep(ij)
        # ij = plaqMap[CartesianIndex(plaqNum)]
        plaqstate = StateRep[ij]
        setindex(StateRep, ij, !plaqstate)
    end
    
    for (i,j) in plaquetteIterator(Conf)
        if canFlipPlaquette(Conf, i, j)
            val = convertToStateRep(plaqMap(i, j))
            push!(states, val)
        end
    end

    states
end
const STAIRCASE_GROWTH_FACTOR = 4e-7

function _generateHilbertSpace(
    InitialState,
    ::SBitVector{UIntType},
    growthfactor = STAIRCASE_GROWTH_FACTOR,
) where {UIntType}
    plaqMapping = PlaqMapping()
    # InitialState = booleanSpinConfig(InitialState) # is actually slower

    nThreads = 1 # multithreading not working yet as plaquette mapping is not thread safe
    # nThreads = Threads.nthreads()

    InitialState_rep = SBitVector{UIntType}(0, 0)

    CurrentStates = ([InitialState_rep])

    NewStates = empty!([InitialState_rep])

    AllStates = DataStructures.SwissDict(InitialState_rep => 1)

    Conf = copy(InitialState)

    neighbors_i = SBitVector{UIntType}[]

    Hrows = SBitVector{UIntType}[]
    Hcols = SBitVector{UIntType}[]

    sizehint!(Hrows, round(Int, growthfactor * exp10(size(Conf, 1))))
    sizehint!(Hcols, round(Int, growthfactor * exp10(size(Conf, 1))))

    while !isempty(CurrentStates)
        for i in eachindex(CurrentStates)
            StateRep = CurrentStates[i]
            Conf = spinConfig!(Conf, StateRep, InitialState, plaqMapping)

            neighbors_i = getNewStates!(neighbors_i, Conf, StateRep, plaqMapping)
            addVertex!(Hrows, Hcols, StateRep, neighbors_i)
            for s in neighbors_i
                if s ∉ keys(AllStates)
                    push!(NewStates, s)
                    AllStates[s] = length(AllStates) + 1
                end
            end
        end

        empty!(CurrentStates)
        append!(CurrentStates, NewStates)
        empty!(NewStates)
    end

    return (; Hrows, Hcols, AllStates, plaqMapping)
end


function generateHilbertSpace(InitialState, type = SBitVector{UInt64}(0, 0))
    (; Hrows, Hcols, AllStates, plaqMapping) = _generateHilbertSpace(InitialState, type)
    H = constructSparseMatrix(Hrows, Hcols, AllStates)
    AllStates = sortByValueOrder(AllStates)
    return HilbertSpace(AllStates, H, plaqMapping)
end

function addRKPotential!(Hilbert::HilbertSpace,μ)
    H = Hilbert.H.data
    H[diagind(H)] .= 0
    H[diagind(H)] .= -μ .* H*ones(size(H, 1))
    return Hilbert
end

function addVertex!(rows::AbstractVector, cols::AbstractVector, i, Neighbors)
    for neigh in Neighbors
        push!(rows, neigh)
        push!(cols, i)
    end
end

function constructSparseMatrix(rows, cols, AllStates)
    nThreads = Threads.nthreads()

    lenRows = length(rows)
    lenCols = length(cols)

    newrows = zeros(Int, lenRows)
    newcols = zeros(Int, lenCols)

    batches = ChunkSplitters.chunks(eachindex(newrows), n = nThreads)
    Threads.@threads for (iChunk, inds) in enumerate(batches)
        for i in inds
            row = rows[i]
            col = cols[i]
            row2 = AllStates[row]
            col2 = AllStates[col]
            newrows[i] = row2
            newcols[i] = col2
        end
    end
    return Symmetric(sparse(newrows, newcols, -1.0))
    # return SparseMatrixCSC(lenRows,lenCols,newrows,newcols,nzval)
end
