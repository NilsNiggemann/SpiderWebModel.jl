import SpiderWebModel as SW
using StaticArrays, CairoMakie
using MakieHelpers
import SpiderWebModel: getStairCase
##

res = SW.getAllNeighborStates(SW.getStairCase(8))
H = SW.H(res.AllStates, res.Neighbors, 1.0)
sol = SW.SolveH(H)
Sq = SW.getStructureFac(SW.spinConfig.(res.AllStates), sol)
##
paths = SW.generateRandomPaths(SW.getStairCase(10), 1200, 30; acceptanceRate = 0.9)
confs = SW.spinConfig.(paths)

Sq = SW.getEqualWeightStructureFac(confs);
# Sq = SW.getStructureFac(confs,sol.vectors[:,1]);
##
with_theme(theme_PiTicks()) do
    fig = Figure(size = 1.2 .* (400, 600))

    axConv = Axis(fig[1, 1], xticks = SimpleTicks(), yticks = SimpleTicks())
    ks = [(0.0, pi), (pi, 2pi)]
    IndRange = eachindex(confs)[100:end]
    for k in ks
        Sq1 = [real(Sq.Sq(k[1], k[2])) for maxIndex in IndRange]
        # Sq1 = [real(Sq.Sq(k[1],k[2],maxIndex)) for maxIndex in IndRange]
        lines!(axConv, IndRange, Sq1)
    end
    axSq = Axis(
        fig[2, 1],
        xlabel = L"q_x",
        ylabel = L"q_y",
        aspect = 1,
        xticks = PiTicks(0:(0.5pi):(2pi)),
        yticks = PiTicks(0:(0.5pi):(2pi)),
    )
    heatmap!(axSq, Sq.k, Sq.k, real(Sq.Sq_k))
    rowsize!(fig.layout, 1, Relative(0.35))
    fig
end
