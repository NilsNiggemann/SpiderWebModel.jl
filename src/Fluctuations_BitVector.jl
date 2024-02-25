using BitIntegers

struct SBitVector{T} <: AbstractVector{Bool}
    x::T
    len::Int
end

Base.size(v::SBitVector) = (v.len,)
# Base.display(v::SBitVector) = println(bitstring(v.x))
# Base.show(io::IO,v::SBitVector) = println(io,bitstring(v.x))
Base.typemax(x::SBitVector) = typemax(x.x)

Base.@propagate_inbounds function Base.getindex(v::SBitVector, i::Int)
    Base.@boundscheck 1 <= i <= length(v)
    return (v.x >> (i-1)) % Bool
end

Base.@propagate_inbounds function Base.setindex(v::SBitVector{UIntType}, i::Int, x::Bool) where {UIntType}
    Base.@boundscheck 1 <= i <= typemax(v) || throw(BoundsError(v, i))
    mask = UIntType(1) << (i-1)
    res = ifelse(x, v.x | mask, v.x & ~mask)

    newlen = _numberOfBits(UIntType) - leading_zeros(res) 
    return SBitVector(res, newlen)
end

@inline _numberOfBits(::Type{UInt8}) = 8
@inline _numberOfBits(::Type{UInt16}) = 16
@inline _numberOfBits(::Type{UInt32}) = 32
@inline _numberOfBits(::Type{UInt64}) = 64
@inline _numberOfBits(::Type{UInt128}) = 128
@inline _numberOfBits(x) = sizeof(x) * 8

Base.:(==)(x::SBitVector,y::SBitVector) = x.x == y.x
Base.hash(x::SBitVector, h::UInt) = hash(x.x, h) 

struct PlaqMapping{I<:Integer}
    d::OrderedDict{Tuple{I,I},I}
end

function (P::PlaqMapping)(ij::Tuple{I,I}) where {I <: Integer} 
    index = get(P.d,ij,0)
    if index == 0
        index = updatePlaqMapping!(P,ij)
    end
    return index
end

(P::PlaqMapping)(i::I,j::I) where {I <: Integer} = P((i,j))
(P::PlaqMapping)(i::Integer) = P.d.keys[i]

Base.length(P::PlaqMapping) = length(P.d)

"""constructs arrays for mapping between plaquette position (i,j) and integers"""
function PlaqMapping()
    d = OrderedDict{Tuple{Int,Int},Int}()
    return PlaqMapping(d)
end

function updatePlaqMapping!(plaqMapping,plaqNum::Tuple{I,I}) where {I<:Integer}
    if plaqNum ∉ keys(plaqMapping.d)
        plaqMapping.d[plaqNum] = length(plaqMapping) + 1
    end
    return length(plaqMapping)
end

getPlaquettes(P::PlaqMapping) = P.d.keys

function spinConfig!(Conf,path::AbstractVector,InitialConf,plaqMap::PlaqMapping)
    Conf .= InitialConf
    for (i,op) in enumerate(path) #this can probably be made faster
        if op
            ij = plaqMap(i)
            flipPlaquette!(Conf,ij)
        end
    end
    return Conf
end
spinConfig(path,InitialConf,plaqMap) = spinConfig!(copy(InitialConf),path,InitialConf,plaqMap)


function getNewStates!(states,Conf,StateRep::SBitVector,plaqMap::PlaqMapping)
    empty!(states)

    @inline function convertToStateRep(ij)
        # ij = plaqMap[CartesianIndex(plaqNum)]
        plaqstate = StateRep[ij]
        setindex(StateRep,ij,!plaqstate) 
    end

    for i in axes(Conf.Mat,1),j in axes(Conf.Mat,2)
        iseven(i+j) && continue
        if canFlipPlaquette(Conf,i,j)
            push!(states,convertToStateRep(plaqMap(i,j)))
        end
    end

    states
end

function getNewStates!(states,Conf,StateRep::BitVector,plaqMap::PlaqMapping)

    @inline function convertToStateRep(ij)
        # ij = plaqMap[CartesianIndex(plaqNum)]
        plaqstate = StateRep[ij]
        newState = copy(StateRep)
        newState[ij] = !plaqstate
        return newState
    end

    for i in axes(Conf.Mat,1),j in axes(Conf.Mat,2)
        iseven(i+j) && continue
        if canFlipPlaquette(Conf,i,j)
            push!(states,convertToStateRep(plaqMap(i,j)))
        end
    end

    states
end


function _generateHamiltonian(InitialState,::SBitVector{UIntType}) where {UIntType}
    plaqMapping = PlaqMapping()

    
    nThreads = 1 # multithreading not working yet as plaquette mapping is not thread safe
    # nThreads = Threads.nthreads()
    
    InitalState_rep = SBitVector{UIntType}(0,0)

    CurrentStates = ([InitalState_rep])

    NewStates_arr = [ [InitalState_rep] for _ in 1:nThreads]

    AllStates = Dict(InitalState_rep => 1)

    Conf_arr = [copy(InitialState) for _ in 1:nThreads]
    
    neighbors_i_arr = [SBitVector{UIntType}[] for _ in 1:nThreads]


    Hrows_arr = [SBitVector{UIntType}[] for _ in 1:nThreads]
    Hcols_arr = [SBitVector{UIntType}[] for _ in 1:nThreads]
    
    Conf = Conf_arr[1]
    neighbors_i = neighbors_i_arr[1]
    NewStates = NewStates_arr[1]
    Hrows = Hrows_arr[1]
    Hcols = Hcols_arr[1]
    
    while !isempty(CurrentStates)
        for i in eachindex(CurrentStates)
            StateRep = CurrentStates[i]
            Conf = spinConfig!(Conf,StateRep,InitialState,plaqMapping)

            neighbors_i = getNewStates!(neighbors_i,Conf,StateRep,plaqMapping)
            addVertex!(Hrows,Hcols,StateRep,neighbors_i)
            append!(NewStates,neighbors_i)

        end

        empty!(CurrentStates)
        
        for NewStates in NewStates_arr
            for s in NewStates
                if s ∉ keys(AllStates)
                    push!(CurrentStates,s)
                    AllStates[s] = length(AllStates) + 1
                end
            end
            empty!(NewStates)
        end

    end

    return (;Hrows_arr,Hcols_arr,AllStates,plaqMapping)

end

function generateHamiltonian(InitialState,type=SBitVector{UInt128}(0,0))
    (;Hrows_arr,Hcols_arr,AllStates,plaqMapping) = _generateHamiltonian(InitialState,type)
    rows = vvcat(Hrows_arr)
    cols = vvcat(Hcols_arr)
    H = constructSparseMatrix(rows,cols,AllStates)
    AllStates = sortByValueOrder(AllStates)
    return (;H,AllStates,plaqMapping)
end

function vvcat(vv::Vector{Vector{T}}) where {T}
    out = Vector{T}(undef, sum(length, vv))
    i = 0
    for v in vv, x in v
       @inbounds out[i+=1] = x
    end
    return out
end

function addVertex!(rows::AbstractVector,cols::AbstractVector,i,Neighbors)
    for neigh in Neighbors
        push!(rows,neigh)
        push!(cols,i)
    end
end

function constructSparseMatrix(rows,cols,AllStates)
    
    nThreads = Threads.nthreads()
    
    lenRows = length(rows)
    lenCols = length(cols)

    newrows = zeros(Int,lenRows)
    newcols = zeros(Int,lenCols)

    batches = ChunkSplitters.chunks(newrows,n=nThreads,split= :batch)
    
    Threads.@threads for (iChunk,inds) in enumerate(batches)
        for i in inds
            row = rows[i]
            col = cols[i]
            row2 = AllStates[row]
            col2 = AllStates[col]
            newrows[i] = row2
            newcols[i] = col2
        end
    end
    return (sparse(newrows,newcols,-1.))
    # return SparseMatrixCSC(lenRows,lenCols,newrows,newcols,nzval)
end