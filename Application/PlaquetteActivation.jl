import SpiderWebModel as SW
using StaticArrays, CairoMakie, Random
##
L = 9
function getBaseConf(L)
    plaqFlip = SW.SpinConfig(fill(NaN, L, L), 1 / 2)
    plaqFlip[(end ÷ 2):(end ÷ 2 + 3), (end ÷ 2 - 1):(end ÷ 2 + 2)] .= [1 0 0 0;
        1 1 0 0;
        0 1 1 1;
        1 0 0 1] .- 1 / 2

    @assert SW.fulFillsConstraint(plaqFlip, verbose = true)
    return plaqFlip
end
plaqFlip = getBaseConf(L)
SW.plotApplPlaquettes(plaqFlip)
##
let
    i, j = SW.getApplicablePlaquettes(plaqFlip, SW.P1)[1]
    Pij = SW.getPlaquette(plaqFlip, i, j)
    Pij .+= SW.P1
    SW.plotApplPlaquettes(plaqFlip)
end
##
let
    i, j = SW.getApplicablePlaquettes(plaqFlip, SW.P1)[1]

    Pij = SW.getPlaquette(plaqFlip, i, j)
    Pij .+= SW.P1
    SW.plotApplPlaquettes(plaqFlip)
end

##
function generateStates(L, numConfigs = 10)
    plaqFlip = getBaseConf(2L + 1)
    Path = SW.xdirecPathReverse(L)
    # Paths = (SW.xdirecPathReverse(L),)
    defaultDelete = 0
    tries = Inf
    maxiter = 10000
    set = SW.setupCalc!(Path, L, L, SW.ALLGS_S12)

    S = [SW.constructConfigPath(SW.DictAlgorithm(),
        plaqFlip,
        SW.ALLGS_S12,
        set;
        maxiter,
        deleteSteps = SW.getStepDeleter(L + 2, defaultDelete, tries),
        verbose = false,
        plotSteps = false) for _ in 1:numConfigs]

    filter!(x -> SW.fulFillsConstraint(x, verbose = false) && !any(isnan, x), S)
    @info "" L defaultDelete tries maxiter length(S)
    return S
end
Random.seed!(4034)
st = generateStates(4, 10)

##
res = SW.getAllNeighborStates(st[1])
##
display.(SW.plotApplPlaquettes.(res.AllConfigs))
