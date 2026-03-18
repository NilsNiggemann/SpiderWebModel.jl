import Pkg; Pkg.activate(@__DIR__)
using CairoMakie, MakieHelpers, HDF5, StatsBase

sectors_05 = h5read("Application/Data/Connectivity/sector_analysis_6x6.h5", "sectors_S05")
sectors_1 = h5read("Application/Data/Connectivity/sector_analysis_6x6.h5", "sectors_S10")
connectivities_05 = h5read("Application/Data/Connectivity/sector_analysis_28x28_staircase.h5", "connectivities_spinhalf")
connectivities_1 = h5read("Application/Data/Connectivity/sector_analysis_28x28_staircase.h5", "connectivities_spin1")
##
with_theme(theme_SimpleTicks()) do
    fig = Figure(size = 0.7 .*(800, 350))
    # ax05 = Axis(fig[1, 1], xlabel="Sector size", ylabel="Count",
    #  title=L"S=1/2",
    # yscale = log10,xscale = log10,
    # xlabelvisible = false, 
    # xticklabelsvisible = false, 
    # )
    ax1 = Axis(fig[1, 1], xlabel="Sector size", ylabel="Count", 
    # title=L"S=1",
    yscale = log10,xscale = log10
    )

    ax_prob = Axis(fig[1, 2], xlabel="Connectivity", ylabel="Probability")

    # linkxaxes!(ax05, ax1)

    
    sector_counts = countmap(sectors_1)
    unique_sectors = sort(collect(keys(sector_counts)))
    counts = [sector_counts[s] for s in unique_sectors]
    println(maximum(counts))

    hist!(ax1,counts, color=:black,bins=sort!(unique(counts)), label=L"S=1")

    sector_counts = countmap(sectors_05)
    unique_sectors = sort(collect(keys(sector_counts)))
    counts = [sector_counts[s] for s in unique_sectors]


    hist!(ax1,counts, color=(:red,0.8),bins=sort!(unique(counts)), label=L"S=1/2")
    println(maximum(counts))

    # text!(ax05, (0.98,0.95), text=L"S=1/2", space = :relative, fontsize = 20, align = (:right, :top) )
    # text!(ax1, (0.98,0.95), text=L"S=1", space = :relative, fontsize = 20, align = (:right, :top) )

    hist!(ax_prob,connectivities_05, color=(:red,0.9), bins=sort!(unique(connectivities_05)), label=L"S=1/2",normalization = :probability)
    hist!(ax_prob,connectivities_1, color=:black, bins=sort!(unique(connectivities_1)), label=L"S=1",normalization = :probability)


    axislegend(ax_prob, position = :rt, framevisible = false)
    Label(fig[1, 1,TopLeft()], L"(a)$$", padding = (-40,0,-25, -20),fontsize = 16)
    # Label(fig[2, 1,TopLeft()], L"(b)$$", padding = (-40,0,-25, -20),fontsize = 16)
    Label(fig[1, 2,TopLeft()], L"(b)$$", padding = (-40,0,-25, -20),fontsize = 16)


    save("Application/figs/PaperFigs/sector_size_histograms.pdf", fig)
    fig
end
##
##
