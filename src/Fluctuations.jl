using LinearAlgebra, SparseArrays, Arpack
const P1 =  SA[
    -1.0  1.0   1.0;
    -1.0  0.0  -1.0;
    1.0  1.0  -1.0]
const P2 = -P1

function getStateNum!(AllConfigs::AbstractDict{T,Int},Conf::T) where T
    num = get(AllConfigs,Conf,0)
    if num == 0
        push!(AllConfigs,Conf => length(AllConfigs)+1)
        return length(AllConfigs)+1
    end
    return num
end

function getNeighborStates!(AllConfigs,StartConfig,operator)

    NeighborStates = Int[]

    plaqs = getApplicablePlaquettes(StartConfig,operator)
    
    for p in plaqs
        Conf2 = copy(StartConfig)
        pij = getPlaquette(Conf2,p...)
        pij .+= operator
        
        num = getStateNum!(AllConfigs,Conf2)
        push!(NeighborStates,num)
    end
    
    return NeighborStates
end

function getAllNeighborStates(StartConfig)
    AllConfigs = Dict(StartConfig => 1)
    
    Nplus = empty(Dict(1 => Int[]))
    Nminus = empty(Dict(1 => Int[]))

    numConfigs_last = 0
    convergenceCounter = 0
    while true
        for (Conf,num) in AllConfigs
            # num in keys(Nplus) && continue
            neighbors = getNeighborStates!(AllConfigs,Conf,P1)
            Nplus[num] = neighbors
            
            neighbors = getNeighborStates!(AllConfigs,Conf,P2)
            Nminus[num] = neighbors
        end
        if length(AllConfigs) == numConfigs_last
            convergenceCounter +=1 
            convergenceCounter >= 2 && break
        end
        numConfigs_last = length(AllConfigs)
    end
    # @info "" length(AllConfigs) length(Nplus) length(Nminus) 
    
    # return (;AllConfigs,Nplus,Nminus)

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

function plotApplPlaquettes!(ax,State,op;kwargs...)
    plaqs = getApplicablePlaquettes(State,op)
    points = Point2f.(plaqs)
    scatter!(ax,points,markersize = 13,color = :red;kwargs...)

end

function plotApplPlaquettes!(ax,State;kwargs...)
    plotApplPlaquettes!(ax,State,P1;kwargs...)
    plotApplPlaquettes!(ax,State,P2;color = :blue,kwargs...)
end

function plotApplPlaquettes(State;kwargs...)
    fig = plotSpinConfig(State)
    plotApplPlaquettes!(current_axis(),State;kwargs...)
    fig
end

function CanApplyPlaquette(Conf::SpinConfig,i::Integer,j::Integer)
    plaquetteIsInBounds(Conf,i,j) || return false
    P = getPlaquette(Conf,i,j)
    □ = plaquetteOperator()
    
    P .+= □

    applicable = fulFillsConstraint(Conf)

    P .-= □

    return applicable
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

function CanApplyNonStrict(Conf::SpinConfig,Op::AbstractMatrix,i,j)
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
function getApplicablePlaquettes(Conf::SpinConfig,Op::AbstractMatrix)
    plaqPos = [(i,j) for i in axes(Conf.Mat,1) for j in axes(Conf.Mat,2) if CanApplyNonStrict(Conf,Op,i,j)]
    return plaqPos
end