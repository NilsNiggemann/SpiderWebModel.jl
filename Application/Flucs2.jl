
Stair = SW.getStairCase(13)
x = SW.ConstructPlaqMapping(size(Stair)...)

plaq = x.plaqMapping[8,7]


State_rep0 = SW.SBitVector(0,length(x.inverseMapping))
State_rep = setindex(State_rep,plaq,true)
##

newSt = SW.getNewStates!(empty([State_rep]),copy(Stair),State_rep0,x.plaqMapping)
##
newc = SW.spinConfig(newSt[25],Stair,x.inverseMapping)
# SW.plotApplPlaquettes(Stair)
SW.plotApplPlaquettes(newc)
##

Stair = SW.getStairCase(14)
x = SW.ConstructPlaqMapping(size(Stair)...)
##
@time H = SW.generateHamiltonian(Stair)
##
@profview H = SW.generateHamiltonian(Stair)
length(H)
# @time SW.SolveH(H)
##
@time begin 
    a = SW.getAllNeighborStates(Stair)
    H1 = SW.H(a.AllStates,a.Neighbors,0)
    SW.SolveH(H1)
end