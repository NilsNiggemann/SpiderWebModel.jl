using BitIntegers
struct SBitVector{T} <: AbstractVector{Bool}
    x::T
    len::Int
end

Base.size(v::SBitVector) = (v.len,)
Base.display(v::SBitVector) = println(bitstring(v.x))
Base.show(io::IO,v::SBitVector) = println(io,bitstring(v.x))
Base.typemax(x::SBitVector) = typemax(x.x)

Base.@propagate_inbounds function Base.getindex(v::SBitVector, i::Int)
    Base.@boundscheck 1 <= i <= length(v)
    return (v.x >> (i-1)) % Bool
end

Base.@propagate_inbounds function Base.setindex(v::SBitVector{UIntType}, i::Int, x::Bool) where {UIntType}
    Base.@boundscheck 1 <= i <= length(v) || throw(BoundsError(v, i))
    mask = UIntType(1) << (i-1)
    res = ifelse(x, v.x | mask, v.x & ~mask)
    return SBitVector(res, v.len)
end

Base.:(==)(x::SBitVector,y::SBitVector) = x.x == y.x
Base.hash(x::SBitVector, h::UInt) = hash(x.x, h) #is it possible to use 128 bit hash?
# Base.hash(x::SBitVector, h::UInt128) = x.x #is it possible to use 128 bit hash?

function spinConfig!(Conf,path::AbstractVector,InitialConf,plaqMap::AbstractVector)
    Conf .= InitialConf
    for (i,op) in enumerate(path) #this can probably be made faster
        if op
            ij = plaqMap[i]
            flipPlaquette!(Conf,ij)
        end
    end
    return Conf
end

spinConfig(path,InitialConf,plaqMap) = spinConfig!(copy(InitialConf),path,InitialConf,plaqMap)

"""constructs arrays for mapping between plaquette position (i,j) and integers"""
function ConstructPlaqMapping(Lx,Ly)
    i = 0
    plaqMapping = zeros(Int,Lx,Ly)
    for x in 1:Lx, y in 1:Ly
        iseven(x+y) && continue
        i += 1
        plaqMapping[x,y] = i
    end

    inverseMapping = [findfirst(==(x),plaqMapping) for x in 1:i]

    return (;plaqMapping,inverseMapping)

end

function ConstructPlaqMapping(Conf::SpinConfig)
    i = 0
    plaqMapping = zeros(Int,size(Conf))
    for x in axes(plaqMapping,1), y in axes(plaqMapping,2)
        iseven(x+y) && continue
        plaquetteIsInBounds(Conf,x,y) || continue
        i += 1
        plaqMapping[x,y] = i
    end

    inverseMapping = [findfirst(==(x),plaqMapping) for x in 1:i]

    return (;plaqMapping,inverseMapping)

end


function getNewStates!(states,Conf,StateRep::SBitVector,plaqMap::AbstractMatrix)
    empty!(states)

    @inline function convertToStateRep(ij)
        # ij = plaqMap[CartesianIndex(plaqNum)]
        plaqstate = StateRep[ij]
        setindex(StateRep,ij,!plaqstate) 
    end

    for i in axes(Conf.Mat,1),j in axes(Conf.Mat,2)
        iseven(i+j) && continue
        if canFlipPlaquette(Conf,i,j)
            push!(states,convertToStateRep(plaqMap[i,j]))
        end
    end

    states
end

function getNewStates!(states,Conf,StateRep::BitVector,plaqMap::AbstractMatrix)
    empty!(states)

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
            push!(states,convertToStateRep(plaqMap[i,j]))
        end
    end

    states
end


function _generateHamiltonian(InitialState,::SBitVector{UIntType}) where {UIntType}
    Lx,Ly = size(InitialState)
    (;plaqMapping,inverseMapping) = ConstructPlaqMapping(InitialState)

    len = length(inverseMapping)

    
    nThreads = Threads.nthreads()
    
    InitalState_rep = SBitVector{UIntType}(0,len)

    @assert 2. ^len < typemax(InitalState_rep) "Type $Type is too small to represent all states. Try a larger type, such as UInt128, UInt256."

    CurrentStates = ([InitalState_rep])

    NewStates_arr = [ [InitalState_rep] for _ in 1:nThreads]

    AllStates = Dict(InitalState_rep => 1)

    Conf_arr = [copy(InitialState) for _ in 1:nThreads]
    
    neighbors_i_arr = [SBitVector{UIntType}[] for _ in 1:nThreads]


    Hrows_arr = [SBitVector{UIntType}[] for _ in 1:nThreads]
    Hcols_arr = [SBitVector{UIntType}[] for _ in 1:nThreads]
    
    
    while !isempty(CurrentStates)
        batches = ChunkSplitters.chunks(CurrentStates,n=nThreads,split= :batch)
        
        Threads.@threads for (iChunk,inds) in enumerate(batches)
            Conf = Conf_arr[iChunk]
            neighbors_i = neighbors_i_arr[iChunk]
            NewStates = NewStates_arr[iChunk]
            Hrows = Hrows_arr[iChunk]
            Hcols = Hcols_arr[iChunk]

            for i in inds
                StateRep = CurrentStates[i]
                Conf = spinConfig!(Conf,StateRep,InitialState,inverseMapping)

                neighbors_i = getNewStates!(neighbors_i,Conf,StateRep,plaqMapping)
                addVertex!(Hrows,Hcols,StateRep,neighbors_i)
                append!(NewStates,neighbors_i)

            end
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

    return (;Hrows_arr,Hcols_arr,AllStates)

end

function generateHamiltonian(InitialState,type=SBitVector{UInt128}(0,0))
    rowsArr,colsArr,AllStates = _generateHamiltonian(InitialState,type)
    rows = vvcat(rowsArr)
    cols = vvcat(colsArr)
    H = constructSparseMatrix(rows,cols,AllStates)
    return (;H,AllStates)
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