
Stair = SW.getStairCase(13)
x = SW.ConstructPlaqMapping(size(Stair)...)

plaq = x.plaqMapping[8,7]


State_rep0 = SW.SBitVector(0,length(x.inverseMapping))
State_rep = setindex(State_rep,plaq,true)
##

newSt = SW.getNewStates!(empty([State_rep])copy(Stair),State_rep0,x.plaqMapping)
##
newc = SW.spinConfig(newSt[25],Stair,x.inverseMapping)
# SW.plotApplPlaquettes(Stair)
SW.plotApplPlaquettes(newc)
##

Stair = SW.getStairCase(13)
x = SW.ConstructPlaqMapping(size(Stair)...)

@time H = SW.generateHamiltonian(Stair)
length(H)
##
@time SW.generateAllPaths(Stair) |> length