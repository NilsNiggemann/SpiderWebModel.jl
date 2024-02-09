struct SBitVector <: AbstractVector{Bool}
    x::UInt128
    len::Int
end
Base.size(v::SBitVector) = (v.len,)
Base.display(v::SBitVector) = println(bitstring(v.x))
Base.show(io::IO,v::SBitVector) = println(io,bitstring(v.x))

Base.@propagate_inbounds function Base.getindex(v::SBitVector, i::Int)
    Base.@boundscheck 1 <= i <= length(v)
    return (v.x >> (i-1)) % Bool
end

Base.@propagate_inbounds function Base.setindex(v::SBitVector, i::Int, x::Bool)
    Base.@boundscheck 1 <= i <= length(v) || throw(BoundsError(v, i))
    mask = UInt128(1) << (i-1)
    res = ifelse(x, v.x | mask, v.x & ~mask)
    return SBitVector(res, v.len)
end

Base.:(==)(x::SBitVector,y::SBitVector) = x.x == y.x
Base.hash(x::SBitVector, h::UInt) = hash(x.x, h) #is it possible to use 128 bit hash?

function spinConfig!(Conf,path::SBitVector,InitialConf,plaqMap::AbstractVector)
    Conf .= InitialConf
    for (i,op) in enumerate(path) #this can probably be made faster
        if op
            ij = plaqMap[i]
            flipPlaquette!(Conf,ij)
        end
    end
    return Conf
end

spinConfig(path::SBitVector,InitialConf,plaqMap) = spinConfig!(copy(InitialConf),path,InitialConf,plaqMap)

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

# function getNewStates!(Conf,State::UInt128,plaqMap::AbstractMatrix,inverseMap::AbstractVector)
#     StateRep = SBitVector(State,length(inverseMap))
#     getNewStates!(Conf,StateRep,plaqMap,inverseMap)
# end

function getNewStates(Conf,StateRep::SBitVector,plaqMap::AbstractMatrix)
    newStates = getApplicablePlaquettes(Conf)
    
    @inline function convertToStateRep(plaqNum)
        ij = plaqMap[CartesianIndex(plaqNum)]
        plaqstate = StateRep[ij]
        setindex(StateRep,ij,!plaqstate) 
    end

    return convertToStateRep.(newStates)
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

function generateHamiltonian(InitialState)
    Lx,Ly = size(InitialState)
    (;plaqMapping,inverseMapping) = ConstructPlaqMapping(Lx,Ly)

    len = length(inverseMapping)

    # InitalState_rep = UInt128(0)
    InitalState_rep = SBitVector(0,len)

    CurrentStates = ([InitalState_rep])
    NewStates = ([InitalState_rep])
    # AllStates = DataStructures.SwissDict(InitalState_rep => 1)
    AllStates = Dict(InitalState_rep => 1)

    Conf = copy(InitialState)
    
    neighbors_i = SBitVector[]
    sizehint!(neighbors_i,Lx*Ly÷2)

    Hrows = SBitVector[]
    Hcols = SBitVector[]
    
    sizehint!(Hrows,2^(2Lx-2))
    sizehint!(Hcols,2^(2Lx-2))

    while !isempty(CurrentStates)
        for State in CurrentStates
            Conf = spinConfig!(Conf,State,InitialState,inverseMapping)

            neighbors_i = getNewStates!(neighbors_i,Conf,State,plaqMapping)
            addVertex!(Hrows,Hcols,State,neighbors_i)

            append!(NewStates,neighbors_i)
        end
        # unique!(NewStates)
        empty!(CurrentStates)
        
        for s in NewStates
            if s ∉ keys(AllStates)
                push!(CurrentStates,s)
                AllStates[s] = length(AllStates) + 1
            end
        end
        empty!(NewStates)

        # @info "" length(AllStates) length(CurrentStates) length(NewStates)
        # CurrentStates,NewStates = NewStates,CurrentStates
    end
    # return Hrows,Hcols,AllStates
    # return sparse(Hrows,Hcols)
    # return Hrows,Hcols
    # return AllStates
    # H = Hermitian(constructSparseMatrix(Hrows,Hcols,AllStates))
    @time H = Hermitian(constructSparseMatrix(Hrows,Hcols,AllStates))
end

function addVertex!(rows::AbstractVector,cols::AbstractVector,i,Neighbors)
    for neigh in Neighbors
        push!(rows,neigh)
        push!(cols,i)
    end
end

function constructSparseMatrix(rows,cols,AllStates)
    
    nThreads = Threads.nthreads()
    newrows = zeros(Int,length(rows))
    newcols = zeros(Int,length(cols))

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
    # return newrows,newcols
    return (sparse(newrows,newcols,-1))
end