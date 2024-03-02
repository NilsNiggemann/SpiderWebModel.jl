import SpiderWebModel as SW
using StaticArrays, CairoMakie
using MakieHelpers
##
function getConfigs(L, numConfigs = 10; kwargs...)
    # paths = [SW.ydirecPathReverse(L),SW.xdirecPathReverse(L),SW.ydirecPath(L),SW.xdirecPath(L),SW.spiralPath(L)]

    path = (SW.xdirecPath(L))
    # Paths = (SW.xdirecPathReverse(L),)
    defaultDelete = 1
    tries = 60
    maxiter = 10_500_000

    setup = SW.setupCalc!(path, L, L, SW.ALLGS_S12)

    S = fetch.([Threads.@spawn SW.constructConfigPath(SW.DictAlgorithm(),
                    L,
                    L,
                    SW.ALLGS_S12,
                    setup;
                    maxiter,
                    deleteSteps = SW.getStepDeleter(L + 2, 1, 15),
                    verbose = false,
                    kwargs...) for _ in 1:numConfigs])

    filter!(x -> SW.fulFillsConstraint(x, verbose = false) && !any(isnan, x), S)
    @info "" L defaultDelete tries maxiter length(S)
    return S
end
##
confs = getConfigs(12, 10)
##
function DoEDs(confs, μ = 1)
    En = Float64[]
    Sq_ks = Matrix{Float64}[]
    for c in confs
        res = SW.getAllNeighborStates(c)
        GC.gc()
        H = SW.H(res.AllConfigs, res.Nplus, res.Nminus, μ)
        println(size(H))
        if size(H)[1] == 1
            continue
        end
        GC.gc()
        sol = SW.SolveH(H)
        Sq = SW.getStructureFac(res.AllConfigs, sol)

        k = Sq.Sq[1].itp.ranges[1]
        Sq_k = fetch.([Threads.@spawn Sq(kx, ky) for kx in k, ky in k])
        push!(Sq_ks, Sq_k)
        push!(En, sol.values[1])
    end
    Sq_k = sum(Sq_ks)
    return (; En, Sq_ks, Sq_k)
end
μ = 0

edSol = DoEDs(confs, μ)
##

with_theme(theme_PiTicks()) do
    fig, ax, hm = heatmap(sum(edSol.Sq_ks[i] for i in eachindex(edSol.En));
        axis = (; aspect = 1, title = L"μ = %$μ"),
        figure = (; size = 0.8 .* (400, 300)))
    Colorbar(fig[1, 2], hm, label = L"\mathcal{S}^{zz}(\mathbf{q})")
    # save("exactFig/Sq_mu=$μ.png",fig)
    fig
end
