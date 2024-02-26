
import SpiderWebModel as SW

##
Stair = SW.getStairCase(14)
@time H = SW.generateHamiltonian(Stair, SW.SBitVector{UInt64}(0, 0))
# length(H.H)
##
ls = let
    ls = Int[]
    for i in 1:15
        Stair = SW.getStairCase(i)
        @time (; Hrows, AllStates) = SW._generateHamiltonian(Stair,
            SW.SBitVector{UInt64}(0, 0))
        push!(ls, length(Hrows))
    end
    ls
end
##
@profview H = SW.generateHamiltonian(Stair)
# length(H)
##
@time SW.SolveHKrylov(H.H)
##
@time begin
    a = SW.getAllNeighborStates(Stair)
    H1 = SW.H(a.AllStates, a.Neighbors, 0)
    # SW.SolveH(H1)
    length(H1)
end
