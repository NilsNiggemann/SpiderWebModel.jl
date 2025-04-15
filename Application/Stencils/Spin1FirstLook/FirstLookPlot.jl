import SpiderWebModel as SW
using SpiderWebModel.HDF5
using MakieHelpers
using CairoMakie
using Statistics
using DataFrames

include("../plottingUtils.jl")
cd(@__DIR__)
outfile_EnergyScaling = "../../Data/energy_mu_S1_3.h5"
# outfiles_Sq = [joinpath(root,file) for (root,_,files) in walkdir("../../Data/open_L24/") for file in files]
outfiles_Sq = [joinpath(root,file) for (root,_,files) in walkdir("../../Data/S0periodic_L28/mu=0.5/") for file in files]

EnGFMC = [h5read(file,"Energy") for file in outfiles_Sq]
SqsGFMC = [h5read(file,"StructureFactor") for file in outfiles_Sq]
inds = filter_outliers(EnGFMC)
# EnGFMC = EnGFMC[inds]
# SqsGFMC = SqsGFMC[inds]
errlines(mean(EnGFMC),std(EnGFMC))

##
function findFilesRecursive(directory)
    files =  [joinpath(root,file) for (root,_,files) in walkdir(directory) for file in files]
    filter!(isfile,files)
    return files
end

function read_data(directory)
    files = findFilesRecursive(directory)

    filekeys = Set{String}()
    function read_maybe(file,key)
        h5open(file, "r") do f
            if key in keys(f)
                return read(f[key])
            else
                return nothing
            end
        end
    end
    for file in files
        h5open(file, "r") do f
            union!(filekeys, keys(f))
        end
    end
    # return read_maybe.(files,Ref(first(filekeys)))
    
    data = DataFrame([k=>read_maybe.(files, Ref(k)) for k in filekeys])
    return data
end

function getEn(res;tau=10)
    dtau = res.τ
    tauIndex = ceil.(Int,tau/dtau)
    En = res.Energy
    En = [En[i][tauIndex[i]] for i in eachindex(En)]
    return En
end
##

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


@time res6x6 = read_data(ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/6x6Condensate_longprop/L=24")

res6x6.En = getEn(res6x6)

res_random_2 = read_data(ENV["MYSCRATCH"]*"Spiderweb/DataS1_CT_RK_equil/RandConf_2*_longprop")
res_random_2.En = getEn(res_random_2)
res_random_3 = read_data(ENV["MYSCRATCH"]*"Spiderweb/DataS1_CT_RK_equil/RandConf_3*_longprop")
res_random_3.En = getEn(res_random_3)
res_random_4 = read_data(ENV["MYSCRATCH"]*"Spiderweb/DataS1_CT_RK_equil/RandConf_4*_longprop")
res_random_4.En = getEn(res_random_4)
res_random_5 = read_data(ENV["MYSCRATCH"]*"Spiderweb/DataS1_CT_RK_equil/RandConf_5*_longprop")
res_random_5.En = getEn(res_random_5)
##
function convertDF(Sector,resDF)
    uniquemu = unique(resDF.mu)
    Energy = [filter(x->x.mu==mu,resDF).En for mu in uniquemu]
    StructureFactor = [mean(filter(x->x.mu==mu,resDF).StructureFactor) for mu in uniquemu]
    # return StructureFactor
    DataFrame(;Sector, mu = uniquemu, Energy,StructureFactor )
end
##
results_df = collect_results("../../Data/EnergyScaling_S1_L20")
results_df_others = convertDF("6x6",res6x6)
results_df_others = vcat(results_df_others, convertDF("rand_2",res_random_2))
results_df_others = vcat(results_df_others, convertDF("rand_3",res_random_3))
results_df_others = vcat(results_df_others, convertDF("rand_4",res_random_4))
results_df_others = vcat(results_df_others, convertDF("rand_5",res_random_5))

labelMap(Sector) = Dict(
    "6x6" => L"6×6",
    "rand_2" => L"\textrm{rand}_1",
    "rand_3" => L"\textrm{rand}_2",
    "rand_4" => L"\textrm{rand}_3",
    "rand_5" => L"\textrm{rand}_4",
)[Sector]
# results_df = vcat(results_df, res6x6_reduced)

function getSectorConfig(sector)
    
    if sector == "6x6"
        L = 12
        S = SW.stencilConfig(zeros(L,L),1,boundaryCondition = :periodic)
        return S .= 2SW.periodicState6x6Condensate(L)
    end
    S = SW.stencilConfig(zeros(36,36),1,boundaryCondition = :periodic)
    data = h5read(ENV["MYSCRATCH"]*"Spiderweb/confs/confs_fixed_conserved_2.h5","confs")
    if sector == "rand_2"
        S .= data[:,:,2]
    elseif sector == "rand_3"
        S .= data[:,:,3]
    elseif sector == "rand_4"
        S .= data[:,:,4]
    elseif sector == "rand_5"
        S .= data[:,:,5]
    else
        error("Sector $sector not found")
    end
    return S
end
# results_df.Config = getSectorConfig.(results_df.Sector, 8)
function plot_overview(results_df,results_df_others)
    fig = Figure(size =  1.1 .* (800, 500), fontsize = 22)
    toprow = fig[1, 1:4] = GridLayout()
    botrow = fig[2 ,1:4]= GridLayout()
    EnergyFig = botrow[1, 1:4] = GridLayout()

    # SqFig = botrow[1, 1] = GridLayout()
    SideRow = fig[1:2, 5] = GridLayout()

    sectors = unique(results_df.Sector)
    sectors_others = unique(results_df_others.Sector)
    
    colors = Makie.ColorSchemes.distinguishable_colors(length(sectors)+length(sectors_others), lchoices = LinRange(0, 80, 100))

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
        SW.plotSpinConfig!(ax_conf, S,constraintkwargs = (;markersize = 10))
        colgap!(toprow, 0)
    end
    Colorbar(toprow[1, 10], colorrange = [-1, 1],
    # ticks = [-1,0,1],
    ticks = ([-1,0,1], [L"|↓\rangle",L"|0\rangle",L"|↑\rangle"]), colormap = cgrad(:greys, 3, categorical = true),tellwidth = false,width = Relative(1))
    
     colsize!(toprow, 10, Relative(0.01))
     colgap!(toprow, 9, -70)

    EnergyFig[1, 1:9] = ax = with_theme(theme_SimpleTicks()) do
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

        tline = errlines!(ax, mus, mean_energies, std_energies, label=L"%$i", markersize=0.0, color = colors[i], linewidth = 1.5)
    end

    colors = colors[length(sectors)+1:end]
    
    
    for (i, sector) in enumerate(sectors_others)
        S = getSectorConfig(sector)
        # S = SW.stencilConfig(conf[1:8, 1:8], 1)

        SideRow[i, 1] = ax_conf = Axis(fig;
            SW.getConfigAxis(S)...,
            yticklabelsvisible = false,
            xminorgridwidth=i==1 ? 0.7 : 0.0,
            yminorgridwidth=i==1 ? 0.7 : 0.0,
            xticklabelsvisible = false,
            spinecolors(colors[i])...,
            spinewidth = 3,
        )
        SW.plotSpinConfig!(ax_conf, S,plotConstraints = false)
        # heatmap!(ax_conf, S, colormap = :greys)
        Label(SideRow[i, 1, Right()], labelMap(sector), padding = (0, 0, 0, 0), rotation = pi/2)
        rowgap!(SideRow, 0)
    end

    for (i, sector) in enumerate(sectors_others)
        sector_data = filter(x->x.Sector==sector&&x.mu>=-0.1,results_df_others)
        L = size(sector_data.StructureFactor[1],1)
        
        mus = sector_data.mu
        energies = sector_data.Energy
        mean_energies = mean.(energies) ./ L^2 #./(1 .-mus)
        std_energies = std.(energies) ./ L^2 #./(1 .-mus)
        
        # println(i," ", sector, " ", mean_energies)
        # tline = errlines!(ax, mus, mean_energies, std_energies, label=L"%$i", markersize=10, color = (colors[i],0.9), linewidth = 1.5)
        tline = scatter!(ax, mus, mean_energies, std_energies, label=L"%$i", markersize=20, color = (colors[i],0.9),marker =  '×')
    end

    # axislegend(ax, position = :rb, merge = true,nbanks=2)

    rowsize!(fig.layout, 1, Relative(0.25))
    colsize!(botrow, 1, Relative(0.4))
    colsize!(fig.layout, 5, Relative(0.2))
    colgap!(fig.layout, 4, 10)
    Label(fig[1, 1,TopLeft()], L"(a)$$", fontsize = 22)
    # Label(fig[2, 1,TopLeft()], L"(b)$$", fontsize = 22)
    Label(SideRow[1, 1,TopLeft()], L"(b)$$", fontsize = 22)
    Label(EnergyFig[1, 1,TopLeft()], L"(c)$$", fontsize = 22)
    save("../../figs/PaperFigs/EnergyScaling_S1_L20.pdf", fig)
    fig
end

plot_overview(results_df,results_df_others)
