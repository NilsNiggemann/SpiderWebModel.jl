import SpiderWebModel as SW
using SpiderWebModel.HDF5
using MakieHelpers
using CairoMakie
using Statistics
include("../plottingUtils.jl")
cd(@__DIR__)
outfile_EnergyScaling = "../../Data/energy_mu_S1_3.h5"
# outfiles_Sq = [joinpath(root,file) for (root,_,files) in walkdir("../../Data/open_L24/") for file in files]
outfiles_Sq = [joinpath(root,file) for (root,_,files) in walkdir("../../Data/S0periodic_L28/mu=0.4/") for file in files]

EnGFMC = [h5read(file,"Energy") for file in outfiles_Sq]
SqsGFMC = [h5read(file,"StructureFactor") for file in outfiles_Sq]
inds = filter_outliers(EnGFMC)
EnGFMC = EnGFMC[inds]
SqsGFMC = SqsGFMC[inds]
errlines(mean(EnGFMC),std(EnGFMC))
##
using DataFrames

function collect_results(output_dir)
    sectorsDirs = readdir(output_dir)
    Sectors = Int[]
    mus = Float64[]
    Energy = Vector{Float64}[]
    StructureFactor = Array{Float64, 3}[]

    for sector in sectorsDirs
        sector_dir = joinpath(output_dir, sector)
        mu_files = readdir(sector_dir)

        for mu_file in mu_files
            file_path = joinpath(sector_dir, mu_file)

            h5open(file_path, "r") do f
                if all(k -> k in keys(f), ("mu", "energies", "structure_factors"))
                    mu = read(f["mu"])
                    energies = read(f["energies"])
                    structure_factors = read(f["structure_factors"])
                    sector = read(f["sector"])
                    push!(mus, mu)
                    push!(Sectors, sector)
                    push!(Energy, energies)
                    push!(StructureFactor, structure_factors)
                end
            end
        end
    end

    DF = DataFrame(Sector=Sectors, mu=mus, Energy=Energy, StructureFactor=StructureFactor)
    sort!(DF, [:Sector, :mu])
    return DF
end

##

results_df = collect_results("../../Data/EnergyScaling_S1_L20")
function plot_overview(results_df)
    fig = Figure(size =  1 .* (700, 500), fontsize = 22)
    
    toprow = fig[1, 1:4] = GridLayout()
    botrow = fig[2 ,1:4]= GridLayout()
    EnergyFig = botrow[1, 2:4] = GridLayout()

    SqFig = botrow[1, 1] = GridLayout()

    sectors = unique(results_df.Sector)
    colors = Makie.ColorSchemes.distinguishable_colors(length(sectors), lchoices = LinRange(0, 80, 100))

    for (i, sector) in enumerate(sectors)
        S = SW.get4x4PeriodicSpinConf(8,sector)
        if sector == 6
            S .= 2SW.getStairCase(size(S, 1)) # plot the staircase state, which is in the same sector
        end
        # S = SW.stencilConfig(conf[1:8, 1:8], 1)

        toprow[1, i] = ax_conf = Axis(fig, title=L"%$i";
            SW.getConfigAxis(S)...,
            yticklabelsvisible = false,
            xticklabelsvisible = false,
            spinecolors(colors[i])...,
            spinewidth = 4,
        )
        SW.plotSpinConfig!(ax_conf, S)
        colgap!(toprow, 0)
    end

    EnergyFig[1, 1:8] = ax = with_theme(theme_SimpleTicks()) do
        Axis(fig, xlabel=L"\mu", 
            ylabel=L"E_0/N_{\text{sites}}"
        )
    end

    for (i, sector) in enumerate(sectors)
        sector_data = results_df[results_df.Sector .== sector, :]
        L = size(sector_data.StructureFactor[1],1)

        mus = sector_data.mu
        energies = sector_data.Energy
        mean_energies = mean.(energies) ./ L^2 #./(1 .-mus)
        std_energies = std.(energies) ./ L^2 #./(1 .-mus)

        tline = errlines!(ax, mus, mean_energies, std_energies, label=L"%$i", markersize=0.1, color = colors[i], linewidth = 1)

        # ann_pos = max(2,round(Int,(i-1)/length(sectors) * length(mus)))

        # txt_pos = Point(mus[ann_pos], mean_energies[ann_pos])
        # txt_str = Makie.latexstring(sector)
        # t1 = scatter!(ax,txt_pos, color = :white,markersize=round(Int, length(txt_str)*20/3)*1px)
        # t2 = text!(ax,txt_pos, text=txt_str,color=colors[i],align=(:center,:center))
        # translate!(t1,0,0,100)
        # translate!(t2,0,0,2i+1)
        # translate!(tline,0,0,2i)
    end
    axislegend(ax, position = :rb, merge = true,nbanks=2)

    rowsize!(fig.layout, 1, Relative(0.25))
    colsize!(botrow, 1, Relative(0.4))

    SqFig[1, 1] = axSQ = Axis(fig, xlabel=L"q_x", ylabel=L"q_y", aspect = 1, xticks = PiTicks([0, pi]), yticks = PiTicks([0, pi]))
    
    let mu_plot = 0.1
        muSim = results_df.mu[findfirst(>=(mu_plot), results_df.mu)]
        resPl = only(filter(x->x.mu==muSim&&x.Sector == 1,results_df))
        # SqMat = dropmean(resPl.StructureFactor,dims=3)
        SqMat = mean(SqsGFMC)[:,:,1]
        
        Sq = SW.getSqCont(SqMat, cutoffEnd = 0)
        qx = qy = trueMomenta(-0.5pi, 1.5pi, size(SqMat, 1))
        hm = heatmap!(axSQ, qx, qy, Sq)
        scatter!(axSQ, Point(pi/2,pi/2), color = :black, markersize = 5)
        Colorbar(SqFig[0,1], hm, label = L"\mathcal{S}(\mathbf{q})", vertical = false,ticks = SimpleTicks())
    end
    Label(fig[1, 1,TopLeft()], L"(a)$$", fontsize = 22)
    Label(fig[2, 1,TopLeft()], L"(b)$$", fontsize = 22)
    Label(EnergyFig[1, 1,TopLeft()], L"(c)$$", fontsize = 22)
    fig
end

plot_overview(results_df)
