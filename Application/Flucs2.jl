
import SpiderWebModel as SW
##
Stair = SW.getStairCase(15)
x = SW.ConstructPlaqMapping(Stair)

plaq = x.plaqMapping[8,7]


State_rep0 = SW.SBitVector(0,length(x.inverseMapping))
State_rep = Base.setindex(State_rep0,plaq,true)
# ##

newSt = SW.getNewStates!(empty([State_rep]),copy(Stair),State_rep0,x.plaqMapping)
# ##
newc = SW.spinConfig(newSt[25],Stair,x.inverseMapping)
# # SW.plotApplPlaquettes(Stair)
SW.plotApplPlaquettes(newc)
##

Stair = SW.getStairCase(14)
x = SW.ConstructPlaqMapping(Stair)
##
@time H = SW.generateHamiltonian(Stair,SW.SBitVector{UInt128}(0,0))
length(H.H)
##
let 
    ls = Int[]
    for i in 1:14
        Stair = SW.getStairCase(14)
        H1 = SW.generateHamiltonian(Stair,SW.SBitVector{UInt128})
        push!(ls,length(H1.H)) 
    end
    ls
end
##
@profview H = SW.generateHamiltonian(Stair)
length(H)
# @time SW.SolveH(H)
##
# @time begin 
#     a = SW.getAllNeighborStates(Stair)
#     H1 = SW.H(a.AllStates,a.Neighbors,0)
#     SW.SolveH(H1)
# end