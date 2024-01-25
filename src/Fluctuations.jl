using LinearAlgebra, SparseArrays, Arpack
const P1 =  SA[
    -1.0  1.0   1.0;
    -1.0  0.0  -1.0;
    1.0  1.0  -1.0]
const P2 = -P1

const P1_SITES = -SVector(getSitesFromPlaquette(P1))/2
const P2_SITES = -SVector(getSitesFromPlaquette(P2))/2

struct PlaquetteFlip end

function plaquetteFlippable(plaqSites::SVector{8})
    return plaqSites == P1_SITES || plaqSites == P2_SITES    
end

function flipPlaquette!(Conf::AbstractMatrix,i,j)
    P = getPlaquette(Conf,i,j)

    P[1,1] = -P[1,1]
    P[1,2] = -P[1,2]
    P[1,3] = -P[1,3]
    P[2,1] = -P[2,1]
    # P[2,2] = -P[2,2]
    P[2,3] = -P[2,3]
    P[3,1] = -P[3,1]
    P[3,2] = -P[3,2]
    P[3,3] = -P[3,3]

    return Conf
end


function flipPlaquette!(Conf::AbstractMatrix,pos::Tuple)
    i,j = pos
    flipPlaquette!(Conf,i,j)
end

function flipPlaquette!(Conf::AbstractMatrix,pos::CartesianIndex{2})
    i,j = Tuple(pos)
    flipPlaquette!(Conf,i,j)
end

function flipPlaquette!(Conf::AbstractMatrix,pos::Int)
    ij = CartesianIndices(Conf)[pos]
    flipPlaquette!(Conf,ij)
end

struct LazyPath{P}
    path::Set{P}
end
Base.isequal(L1::LazyPath,L2::LazyPath) = L1.path == L2.path
Base.:(==)(L1::LazyPath,L2::LazyPath) = L1.path == L2.path
Base.hash(L::LazyPath) = hash(L.path)
Base.copy(L::LazyPath) = LazyPath(copy(L.path))
Base.empty!(L::LazyPath) = LazyPath(empty!(L.path))

function appendToPath!(path,pos)
    if pos ∈ path
        delete!(path,pos)
    else
        push!(path,pos)
    end
    return path
end

function Base.push!(L::LazyPath,pos)
    appendToPath!(L.path,pos)
    return L
end

function spinConfig!(Conf::AbstractMatrix,path,InitConf::AbstractMatrix)
    Conf .= InitConf
    for op in path
        flipPlaquette!(Conf,op)
    end
    return Conf
end
spinConfig!(Conf::AbstractMatrix,L::LazyPath,InitConf::AbstractMatrix) = spinConfig!(Conf,L.path,InitConf)

function spinConfig(InitConf::AbstractMatrix,path)
    NewConf = copy(InitConf)
    spinConfig!(NewConf,path,InitConf)
end

function generateAllConfigs(InitialState)
    alreadyDone = Set(Int[])
    startpath = empty(Set(1))

    AllPaths = Dict(
        startpath => 1
    )
    Conf = copy(InitialState)
    while length(alreadyDone) < length(AllPaths)
        for (path,num) in AllPaths
            # println(length(alreadyDone),"/",length(AllPaths))
            num in alreadyDone && continue
            
            Conf = spinConfig!(Conf,path,InitialState)
            appendToPaths!(AllPaths,Conf,path)

            push!(alreadyDone,num)

        end
    end
    return Dict(spinConfig(InitialState,path) => num for (path,num) in AllPaths)
    # return AllPaths
end

function appendToPaths!(AllPaths,Conf,path)

    plaqs = getApplicablePlaquettes(Conf)
    LI = LinearIndices(Conf)
    for p in plaqs
        newpath = copy(path)
        pInt = LI[CartesianIndex(p)]
        push!(newpath,pInt)
        
        ind = get(AllPaths,newpath,0)
        if ind == 0 
            AllPaths[newpath] = length(AllPaths)+1
        end
    end
    
    return AllPaths
end

function getNeighborStates!(AllConfigs,State,operator)

    NeighborStates = Int[]

    plaqs = getApplicablePlaquettes(State,operator)
    
    for p in plaqs
        NewState = copy(State)
        pij = getPlaquette(NewState,p...)
        pij .+= operator
        
        num = getStateNum!(AllConfigs,NewState)
        push!(NeighborStates,num)
    end
    
    return NeighborStates
end


function getStateNum!(AllConfigs::AbstractDict{T,Int},Conf::T) where T
    num = get(AllConfigs,Conf,0)
    if num == 0
        push!(AllConfigs,Conf => length(AllConfigs)+1)
        return length(AllConfigs)+1
    end
    return num
end

function getAllNeighborStates(StartConfig)
    AllConfigs = generateAllConfigs(StartConfig)
    
    Nplus = empty(Dict(1 => Int[]))
    Nminus = empty(Dict(1 => Int[]))

    len = length(AllConfigs)
    for (Conf,num) in AllConfigs
        neighbors = getNeighborStates!(AllConfigs,Conf,P1)
        Nplus[num] = neighbors
        
        neighbors = getNeighborStates!(AllConfigs,Conf,P2)
        Nminus[num] = neighbors
    end
    # @assert length(AllConfigs) == len "new Configs were generated!"
    return (;
    AllConfigs = sortByValueOrder(AllConfigs),
    Nplus = sortbyKeyOrder(Nplus),
    Nminus = sortbyKeyOrder(Nminus),
    )
end

function sortByValueOrder(D::Dict)
    ks = collect(keys(D))
    vals = collect(values(D))
    return ks[sortperm(vals)]
end

function sortbyKeyOrder(D::Dict)
    ks = collect(keys(D))
    vals = collect(values(D))
    return vals[sortperm(ks)]
end

function invertDict(D)
    D2 = Dict(v => k for (k,v) in D)
end

function H(AllStates,neighborplus,neighborminus,mu::T=0.) where {T<:Number}
    dim = length(AllStates)
    rows = Int[]
    cols = Int[]
    vals = T[]
    function addTerm!(n,m,val)
        push!(rows,n)
        push!(cols,m)
        push!(vals,val)
    end

    for n in 1:dim
        for m ∈ neighborplus[n]
            addTerm!(n,m,-one(T))
        end
        for m ∈ neighborminus[n]
            addTerm!(n,m,-one(T))
        end

        val = mu*(length(neighborplus[n] ∪ neighborminus[n]))
        addTerm!(n,n,val)
    end

    return Hermitian(sparse(rows,cols,vals))
end


function Hnonflip(AllStates,neighborplus,neighborminus)
    dim = length(AllStates)
    H = zeros(dim,dim)
    for n in axes(H,1), m in axes(H,2)
        if m ∈ neighborplus[n] || m ∈ neighborminus[n] 
            H[n,m] = 1
        end
    end
    return Hermitian(H)
end

function testNplusMinus(nplus,nminus)
    @testset "test nplus" failfast = true begin
        for (n,ns) in enumerate(nplus)
            for m in ns
                @test n in nminus[m]
            end
        end
    end
    @testset "test nminus" failfast = true begin
        for (n,ns) in enumerate(nminus)
            for m in ns
                @test n in nplus[m]
            end
        end
    end
end

SolveH(H,range=1:1) = eigen(H,range)
const SparseMat = Union{SparseMatrixCSC,Hermitian{<:Number,<:SparseMatrixCSC}}

function SolveH(H::SparseMat;kwargs...)
    values,vectors = eigs(H,nev=1,which=:SR;kwargs...)
    return (;values,vectors)
end

function SolveH(AllConfigs,Nplus,Nminus,mu,range = 1:1)
    H = H(AllConfigs,Nplus,Nminus,mu)
    return SolveH(H,range)
end

function getMagnetization(AllConfigs,eigen,i)
    ψ0 = eigen.vectors[:,1]
    mag = zero(eltype(eigen.vectors))
    for n in eachindex(ψ0)
        Si = AllConfigs[n][i]
        mag += abs2(ψ0[n])*Si
    end
    return mag
end

function plotApplPlaquettes!(ax,State;kwargs...)
    plaqs = getApplicablePlaquettes(State)
    points = Point2f.(plaqs)
    scatter!(ax,points,markersize = 13,color = :red;kwargs...)

end

function plotApplPlaquettes(State;heatmapkwargs = (;),kwargs...)
    fig = plotSpinConfig(State;heatmapkwargs...)
    plotApplPlaquettes!(current_axis(),State;kwargs...)
    fig
end

function flipSpinsAlongLine!(Conf,org,slope)
    slope ∈ (-Inf, Inf) && return flipSpinsAlongRow!(Conf,org[2])
    for i in axes(Conf.Mat,1), j in axes(Conf.Mat,2)
        if slope*(i-org[1]) == j-org[2]
            Conf[i,j] *= -1
        end
    end
    return Conf
end

function flipSpinsAlongRow!(Conf,i)
    Conf[i,:] .*= -1
    return Conf
end

function CanApply(Conf::SpinConfig,Op::AbstractMatrix,i,j)
    plaquetteIsInBounds(Conf,i,j) || return false
    P = getPlaquette(Conf,i,j)
    
    P .+= Op

    applicable = fulFillsConstraint(Conf)

    P .-= Op

    return applicable
end

function canFlipPlaquette(Conf::SpinConfig,i,j)
    plaquetteIsInBounds(Conf,i,j) || return false
    isodd(i+j) || return false

    P = getPlaquette(Conf,i,j)
    sites = SVector(getSitesFromPlaquette(P))
    return plaquetteFlippable(sites)
end

function CanApplyNonStrict(Conf::SpinConfig,Op,i,j)
    plaquetteIsInBounds(Conf,i,j) || return false
    isodd(i+j) || return false
    P = getPlaquette(Conf,i,j)
    
    P .+= Op

    applicable = all(x->abs(x)<=Conf.S,P)

    P .-= Op

    return applicable
end

function CanApplyAnywhere(Conf::SpinConfig,Op::AbstractMatrix)
    a1 = axes(Conf.Mat,1)
    a2 = axes(Conf.Mat,2)

    Opx,Opy = size(Op)
    for i in a1, j in a2
        firstindex(a1)+Opx <= i <= lastindex(a1)-Opx || continue
        firstindex(a2)+Opy <= j <= lastindex(a2)-Opy || continue

        if CanApply(Conf,Op,i,j)
            return true
        end
    end
    return false
end

"""assumes that Op is already an allowed operator"""
function getApplicablePlaquettes(Conf::SpinConfig,Op)
    plaqPos = [(i,j) for i in axes(Conf.Mat,1) for j in axes(Conf.Mat,2) if CanApplyNonStrict(Conf,Op,i,j)]
    return plaqPos
end

function getApplicablePlaquettes(Conf::SpinConfig)
    plaqPos = [(i,j) for i in axes(Conf.Mat,1) for j in axes(Conf.Mat,2) if canFlipPlaquette(Conf,i,j)]
    return plaqPos
end
