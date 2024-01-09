using LinearAlgebra, SparseArrays, Arpack
const P1 =  SA[
    -1.0  1.0   1.0;
    -1.0  0.0  -1.0;
    1.0  1.0  -1.0]
const P2 = -P1

function generateAllFluctuations(StartConfig,operator = findOps(StartConfig)[1];maxiter = 500)

    Configs = [[copy(StartConfig)]]
    uniqeConfs = Set([copy(StartConfig)])
    for iter in 1:maxiter-1
        Confs = Configs[iter]
        newconf = empty(Confs)
        for Conf in Confs
            plaqs = getApplicablePlaquettes(Conf,operator)
            # isempty(plaqs) && (@warn "no fluctuations possible"; break)
            
            for p in plaqs
                Conf2 = copy(Conf)
                pij = getPlaquette(Conf2,p...)
                pij .+= operator
                if Conf2 ∉ uniqeConfs
                    push!(newconf,Conf2)
                    push!(uniqeConfs,Conf2)
                end
            end
            if isempty(newconf) 
                @info "" length(uniqeConfs)
                return Configs,uniqeConfs
            end
            # @assert i1 == i2 "i1 = $i1, i2 = $i2"
        end
        # @info "" length(newconf) length(uniqeConfs)
        push!(Configs,newconf)
    end
    # filter!(x -> fulFillsConstraint(x,verbose = false),Configs)
    # @assert all(fulFillsConstraint.(Configs,verbose = true))
    @warn "max iterations reached" length(uniqeConfs)
    return Configs,uniqeConfs
end

function getStateNum!(AllConfigs::Dict{T,Int},Conf::T) where T
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

    while length(Nplus) < length(AllConfigs)
        for (Conf,num) in AllConfigs
            # num in keys(Nplus) && continue
            neighbors = getNeighborStates!(AllConfigs,Conf,P1)
            Nplus[num] = neighbors
            
            neighbors = getNeighborStates!(AllConfigs,Conf,P2)
            Nminus[num] = neighbors
        end
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

function H(AllStates,neighborplus,neighborminus,mu=0)
    dim = length(AllStates)
    H = zeros(dim,dim)
    for n in axes(H,1), m in axes(H,2)
        if m ∈ neighborplus[n] || m ∈ neighborminus[n] 
            H[n,m] = -1
        end

        if n == m
            H[n,m] = mu*(length(neighborplus[n] ∪ neighborminus[n]))
        end
    end

    return sparse(Hermitian(H))
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

# function swapStates!(AllStates,Nplus,Nminus,i,j)
#     # temp = copy(AllStates[i])
#     # AllStates[i] .= AllStates[j]
#     # AllStates[j] .= temp
#     AllStates[i],AllStates[j] = AllStates[j],AllStates[i]
#     for neighborlist in (Nplus,Nminus)
#         for (i_n,n) in enumerate(neighborlist)
#             if n == i
#                 neighborlist[i_n] = j
#             elseif n == j
#                 neighborlist[i_n] = i
#             end
#         end
#     end

# end
SolveH(H,range=1:1) = eigen(H,range)
function SolveH(H::SparseMatrixCSC;kwargs...)
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
    plaqs = getApplicablePlaquettes_ns(State,op)
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