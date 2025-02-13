import SpiderWebModel as SW
using SpiderWebModel.HDF5
using MakieHelpers
using CairoMakie
using Statistics
include("../plottingUtils.jl")
cd(@__DIR__)
outfile_EnergyScaling = "../../Data/energy_mu_S1_2.h5"
outfiles_Sq = [joinpath(root,file) for (root,_,files) in walkdir("../../Data/open_L24/") for file in files]

EnGFMC = [h5read(file,"Energy") for file in outfiles_Sq]
SqsGFMC = [h5read(file,"StructureFactor") for file in outfiles_Sq]
errlines(mean(EnGFMC),std(EnGFMC))
##

function plot_overview(outfile_en,Sqs)
    fig = Figure(size = (700, 500),fontsize = 22)
    
    EnergyFig = fig[2,1:3] = GridLayout()
    SqFig = fig[2,4] = GridLayout()
    toprow = fig[1,1:4] = GridLayout()
    # Top row for configurations

    h5open(outfile_en, "r") do f
        
        sectors = sort(filter(!=("muRange"),keys(f)),by=x->parse(Int,x))
        # e_start = [haskeymean(f["$sector/energy"][1,:]) for sector in sectors]

        # perm = sortperm(e_start)
        # sectors = sectors[perm]

        colors = Makie.ColorSchemes.distinguishable_colors(length(sectors), lchoices = LinRange(0,80,100))

        for (i, sector) in enumerate(sectors)
            conf = read(f["$sector/conf"])
            if sector == "6"
                conf .= 4SW.getStairCase(size(conf,1)) #plot the staircase state, which is in the same sector
            end
            S = SW.stencilConfig(conf[1:8,1:8],1)
            toprow[1,i] = ax_conf = Axis(fig, title=L"%$i";SW.getConfigAxis(S)...,
            yticklabelsvisible=false,
            xticklabelsvisible=false,
            spinecolors(colors[i])...,
            spinewidth = 4
            )
            SW.plotSpinConfig!(ax_conf, S)
            colgap!(toprow, 0)
        end

        # Bottom row for energies
        EnergyFig[1,1:8] = ax = with_theme(theme_SimpleTicks()) do
            Axis(fig, xlabel=L"\mu", 
            ylabel=L"E_0/N_{\text{sites}}"
            # ylabel=L"E_0/(N_{\text{sites}}(1-\mu))"
            )
        end
        L = 20
        muRange = f["muRange"][:]
        for (i, sector) in enumerate(sectors)
            "energy" in keys(f[sector]) || continue
            Nsites = length(f[sector*"/conf"])
            energies = read(f["$sector/energy"])
            mean_energies = dropmean(energies, dims=2) ./Nsites# ./ (1 .-muRange)
            std_energies = dropstd(energies, dims=2) ./Nsites# ./ (1 .-muRange)

            errlines!(ax, muRange, mean_energies, std_energies, label=L"%$i", markersize=0.1,color = colors[i],linewidth = 2)
        end
        axislegend(ax, position=:lt, merge=true)
    end

    rowsize!(fig.layout,1,Relative(0.2))
    # colsize!(fig.layout,2,Relative(0.3))

    SqFig[1,1] = axSQ = Axis(fig, xlabel=L"q_x", ylabel=L"q_y",aspect=1,xticks = PiTicks([0,pi]),yticks = PiTicks([0,pi]))
    
    SqMat = mean(Sqs)[:,:,1]
    SqErr = std(Sqs)[:,:,1]
    
    Sq = SW.getSqCont(SqMat,cutoffEnd=0)
    qx = qy = trueMomenta(-0.5pi,1.5pi,size(SqMat,1))
    qs = Iterators.product(qx,qy)
    heatmap!(axSQ,qx,qy, Sq.(qs))

    fig
end

plot_overview(outfile_EnergyScaling,SqsGFMC)
