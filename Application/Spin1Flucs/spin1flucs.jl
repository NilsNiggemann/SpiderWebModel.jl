import SpiderWebModel as SW
using CairoMakie
using MakieHelpers
using HDF5
##
StartState = SW.SpinConfig(zeros(100, 100), 1)
# StartState = SW.SpinConfig(Array(SW.getStairCase(24)),1)
SW.plotApplPlaquettes(StartState, SW.P2)
##
function getRandomConfs(StartState, N, acceptanceRate; prethermalizationSteps = 1000)
    confs = [[copy(StartState)] for _ = 1:Threads.nthreads()]

    nthreads = Threads.nthreads()
    Threads.@threads for ithread = 1:nthreads
        confs_i = confs[ithread]
        Conf = confs_i[1]
        FlipPlaqs = empty!([(0, 0)])
        iter = 0
        while true
            iter += 1
            PCurr = rand((SW.P1, SW.P2))

            FlipPlaqs = SW.getApplicablePlaquettes!(FlipPlaqs, Conf, PCurr)
            # FlipPlaqs = SW.getRandomSeparatedPlaquettes!(FlipPlaqs, sepPlaqs, Conf,PCurr)
            if isempty(FlipPlaqs)
                continue
            end
            i, j = rand(FlipPlaqs)
            P = SW.getPlaquette(Conf, i, j)
            P .+= PCurr
            if iter > prethermalizationSteps && rand() < acceptanceRate
                push!(confs_i, copy(Conf))
            end
            if length(confs_i) >= N / nthreads
                break
            end
        end
    end
    allconfs = confs[begin]
    for i = 2:nthreads
        append!(allconfs, confs[i])
    end
    return allconfs
end
##
@time confs = getRandomConfs(StartState, 2000, 0.0000001; prethermalizationSteps = 8000)
##
h5write("Application/ConfsRaw/Spin1Fluc_6.h5","Confs",stack(confs,dims=3))

# Sq = SW.getEqualWeightStructureFac(confs)
# ##
# with_theme(theme_PiTicks()) do
#     fig = Figure()
#     ax = Axis(fig[1, 1], aspect = 1)

#     hm = heatmap!(ax, Sq.k, Sq.k, real(Sq.Sq_k))
#     # @info "" sumrule = sum(real(Sq.Sq_k[1:end-1, 1:end-1])) sumRuleExpected = L^2 / 4
#     Colorbar(fig[1, 2], hm)
#     fig
#     save("Application/figs/Sq_S1_onlyS0Fluctuations.pdf", fig)
#     fig
# end
