struct SBitVector{T} <: AbstractVector{Bool}
    x::T
    len::Int
end

function getlen(x::T, numBits = _numberOfBits(T)) where {T}
    return numBits - leading_zeros(x)
end

SBitVector(x) = SBitVector(x, getlen(x))

Base.size(v::SBitVector) = (v.len,)
# Base.display(v::SBitVector) = println(bitstring(v.x))
# Base.show(io::IO,v::SBitVector) = println(io,bitstring(v.x))
Base.typemax(x::SBitVector) = typemax(x.x)

Base.@propagate_inbounds function Base.getindex(v::SBitVector, i::Int)
    Base.@boundscheck 1 <= i <= length(v)
    return (v.x >> (i - 1)) % Bool
end

Base.@propagate_inbounds function Base.setindex(
    v::SBitVector{UIntType},
    i::Int,
    x::Bool,
) where {UIntType}
    numBits = _numberOfBits(UIntType)
    Base.@boundscheck 1 <= i <= numBits || throw(BoundsError(v, i))
    mask = UIntType(1) << (i - 1)
    res = ifelse(x, v.x | mask, v.x & ~mask)

    newlen = getlen(res, numBits)
    return SBitVector(res, newlen)
end

@inline _numberOfBits(::Type{UInt8}) = 8
@inline _numberOfBits(::Type{UInt16}) = 16
@inline _numberOfBits(::Type{UInt32}) = 32
@inline _numberOfBits(::Type{UInt64}) = 64
@inline _numberOfBits(::Type{UInt128}) = 128
@inline _numberOfBits(x) = sizeof(x) * 8

Base.:(==)(x::SBitVector, y::SBitVector) = x.x == y.x
Base.hash(x::SBitVector, h::UInt) = hash(x.x, h)

struct PlaqMapping{I<:Integer}
    d::OrderedDict{Tuple{I,I},I}
end

function (P::PlaqMapping)(ij::Tuple{I,I}) where {I<:Integer}
    index = get(P.d, ij, 0)
    if index == 0
        index = updatePlaqMapping!(P, ij)
    end
    return index
end

(P::PlaqMapping)(i::I, j::I) where {I<:Integer} = P((i, j))
(P::PlaqMapping)(i::Integer) = P.d.keys[i]

Base.length(P::PlaqMapping) = length(P.d)
Base.:(==)(P1::PlaqMapping, P2::PlaqMapping) = P1.d == P2.d

"""constructs arrays for mapping between plaquette position (i,j) and integers"""
function PlaqMapping()
    d = OrderedDict{Tuple{Int,Int},Int}()
    return PlaqMapping(d)
end

function updatePlaqMapping!(plaqMapping, plaqNum::Tuple{I,I}) where {I<:Integer}
    if plaqNum ∉ keys(plaqMapping.d)
        plaqMapping.d[plaqNum] = length(plaqMapping) + 1
    end
    return length(plaqMapping)
end

getPlaquettes(P::PlaqMapping) = P.d.keys

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

    for i in axes(Conf.Mat, 1), j in axes(Conf.Mat, 2)
        iseven(i + j) && continue
        if canFlipPlaquette(Conf, i, j)
            push!(states, convertToStateRep(plaqMap(i, j)))
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

    NewStates = [InitialState_rep]

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

struct HilbertSpace{StatesType,HType<:AbstractMatrix,PlaqMapType}
    AllStates::StatesType
    H::HType
    plaqMapping::PlaqMapType
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

    batches = ChunkSplitters.chunks(newrows, n = nThreads)

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
