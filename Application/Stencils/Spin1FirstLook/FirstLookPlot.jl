import SpiderWebModel as SW
using SpiderWebModel.HDF5
using MakieHelpers
using CairoMakie
using Statistics
using DataFrames
using MakieExtra

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
# results_df_others = vcat(results_df_others, convertDF("rand_3",res_random_3))
results_df_others = vcat(results_df_others, convertDF("rand_4",res_random_4))
results_df_others = vcat(results_df_others, convertDF("rand_5",res_random_5))

labelMap(Sector) = Dict(
    "6x6" => L"6×6",
    "rand_2" => L"\textrm{rand}_1",
    # "rand_3" => L"\textrm{rand}_2",
    "rand_4" => L"\textrm{rand}_2",
    "rand_5" => L"\textrm{rand}_3",
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
    fig = Figure(size =  1.0 .* (900, 550), fontsize = 22)
    toprow = fig[1, 1:6] = GridLayout()
    botrow = fig[2 ,1:4]= GridLayout()
    EnergyFig = botrow[1, 1:4] = GridLayout()

    # SqFig = botrow[1, 1] = GridLayout()
    SideRow = fig[2, 5:6] = GridLayout()

    sectors = unique(results_df.Sector)
    sectors_others = unique(results_df_others.Sector)
    
    colors = Makie.ColorSchemes.distinguishable_colors(length(sectors)+length(sectors_others), lchoices = LinRange(0, 80, 100))

        let
        S = SW.get4x4PeriodicSpinConf(8,6)

        ax_stair = toprow[1, 6] = Axis(fig;
            SW.getConfigAxis(S)...,
            yticklabelsvisible = false,
            valign = 0.5,
            halign = 0.5,
            height = Relative(0.9),
            width = Relative(0.9),
            xticklabelsvisible = false,
            spinecolors(colors[6])...,
            xticksvisible = false,
            yticksvisible = false,
            spinewidth = 2,
            alignmode = Mixed(left = 0,right = 0, top = -20,bottom = 20)
        )
        SW.plotSpinConfig!(ax_stair, S .= 2SW.getStairCase(size(S, 1)),constraintkwargs = (;markersize = 10))
        translate!(ax_stair.blockscene,0,0,-1000)
    end

    for (i, sector) in enumerate(sectors)
        S = SW.get4x4PeriodicSpinConf(8,sector)
        # if sector == 6
            # S .= 2SW.getStairCase(size(S, 1)) # plot the staircase state, which is in the same sector
            # Snew = copy(S)
            # for I in CartesianIndices(S)
            #     Ii,Ij = Tuple(I)
            #     Snew[end - Ii+1, Ij] = S[I]
            # end
            # S .= Snew
        # end
        # S = SW.stencilConfig(conf[1:8, 1:8], 1)

        shift_align = i <= 6 ? 50 : 0
        # shift_align = i == 6 ? 50 : 0
        toprow[1, i] = ax_conf = Axis(fig, title=L"%$i";
            SW.getConfigAxis(S)...,
            yticklabelsvisible = false,
            xticklabelsvisible = false,
            spinecolors(colors[i])...,
            xticksvisible = false,
            yticksvisible = false,
            spinewidth = 4,
            alignmode = Mixed(left = -shift_align,right = shift_align)

        )
        SW.plotSpinConfig!(ax_conf, S,constraintkwargs = (;markersize = 10))
        colgap!(toprow, 8)
    end


    # Colorbar(toprow[1, 10], colorrange = [-1, 1],
    # ticks = [-1,0,1],
    # ticks = ([-1,0,1], [L"|↓\rangle",L"|0\rangle",L"|↑\rangle"]), colormap = cgrad(:greys, 3, categorical = true),tellwidth = false,width = Relative(1))
    
    #  colsize!(toprow, 10, Relative(0.01))
    #  colgap!(toprow, 9, -70)

    EnergyFig[1, 1:9] = ax = with_theme(theme_SimpleTicks()) do
        Axis(fig, xlabel=L"\mu", 
            ylabel=L"E_0/(J'\times N_{\text{sites}})",
            # xgridvisible = false,
            # ygridvisible = false,
        )
    end

    EnergyFig[1, 1:9] = inset_En = Axis(fig,width = Relative(0.3),height = Relative(0.3),halign=0.17, valign=0.95,xticklabelsize= 16,yticklabelsize= 16,backgroundcolor = (:white,1.2))
    EnergyFig[1, 1:9] = inset_En2 = Axis(fig,width = Relative(0.3),height = Relative(0.3),halign=0.8, valign=0.2,xticklabelsize= 16,yticklabelsize= 16,backgroundcolor = (:white,1.2))

    translate!(inset_En.blockscene,0,0,1000)
    translate!(inset_En2.blockscene,0,0,1000)
    zoom_lines!(ax,inset_En)
    zoom_lines!(ax,inset_En2)

    leftzoom = (-0.1,0.15)
    rightzoom = (0.8,0.92)

    # for (i, sector) in enumerate(reverse(sectors))
    for i in reverse(eachindex(sectors))
        sector = sectors[i]
        sector_data = results_df[results_df.Sector .== sector, :]
        L = size(sector_data.StructureFactor[1],1)

        mus = sector_data.mu
        energies = sector_data.Energy
        mean_energies = mean.(energies) ./ L^2 #./(1 .-mus)
        std_energies = std.(energies) ./ L^2 #./(1 .-mus)

        tline = errlines!(ax, mus, mean_energies, std_energies; label=L"%$i", markersize=0.0, color = colors[i], linewidth = 1.5)
        # if i == 1
        #     translate!(tline,(0,0,10))
        # end
        zoominds = findall(x->x>leftzoom[1]&&x<leftzoom[2], mus)

        tline = errlines!(inset_En, mus[zoominds], mean_energies[zoominds], std_energies[zoominds]; label=L"%$i", markersize=0.0, color = colors[i], linewidth = 3)
        zoominds = findall(x->x>rightzoom[1]&&x<rightzoom[2], mus)

        tline = errlines!(inset_En2, mus[zoominds], mean_energies[zoominds], std_energies[zoominds]; label=L"%$i", markersize=0.0, color = colors[i], linewidth = 2)
    end

    colors = colors[length(sectors)+1:end]
    
    
    SideRowInds = Tuple.(CartesianIndices((2,2)))
    for (i, sector) in enumerate(sectors_others)
        S = getSectorConfig(sector)
        # S = SW.stencilConfig(conf[1:8, 1:8], 1)
        
        jax,iax = SideRowInds[i]
        SideRow[iax,jax] = ax_conf = Axis(fig;
            SW.getConfigAxis(S)...,
            yticklabelsvisible = false,
            xminorgridwidth=i==1 ? 0.7 : 0.0,
            yminorgridwidth=i==1 ? 0.7 : 0.0,
            xticklabelsvisible = false,
            xticksvisible = false,
            yticksvisible = false,
            spinecolors(colors[i])...,
            spinewidth = 3,
        )
        SW.plotSpinConfig!(ax_conf, S,plotConstraints = sector == "6x6",constraintkwargs = (;markersize = 10))
        # heatmap!(ax_conf, S, colormap = :greys)
        Label(SideRow[iax, jax, Top()], labelMap(sector), padding = (0, 0, 0, 0), rotation = 0)
        rowgap!(SideRow, 1)
    end

    for (i, sector) in enumerate(sectors_others)
        sector_data = filter(x->x.Sector==sector&&x.mu>=-0.1,results_df_others)
        L = size(sector_data.StructureFactor[1],1)
        
        mus = sector_data.mu
        println(mus)
        energies = sector_data.Energy
        mean_energies = mean.(energies) ./ L^2 #./(1 .-mus)
        std_energies = std.(energies) ./ L^2 #./(1 .-mus)
        
        # println(i," ", sector, " ", mean_energies)
        # tline = errlines!(ax, mus, mean_energies, std_energies, label=L"%$i", markersize=10, color = (colors[i],0.9), linewidth = 1.5)
        tline = scatter!(ax, mus, mean_energies, std_energies, label=L"%$i", markersize=30, color = (colors[i],0.9),marker =  '×')
        zoominds = findall(x->x>leftzoom[1]&&x<leftzoom[2], mus)
        tline = scatter!(inset_En, mus[zoominds], mean_energies[zoominds], std_energies[zoominds], label=L"%$i", markersize=40, color = (colors[i],0.9),marker =  '×')
        zoominds = findall(x->x>rightzoom[1]&&x<rightzoom[2], mus)

        tline = scatter!(inset_En2, mus[zoominds], mean_energies[zoominds], std_energies[zoominds], label=L"%$i", markersize=40, color = (colors[i],0.9),marker =  '×')
    end

    # axislegend(ax, position = :rb, merge = true,nbanks=2)
    elem_1 = [MarkerElement(color = :black, marker = '◼', markersize = 30,strokecolor = :black,strokewidth=2)]
    elem_2 = [MarkerElement(color = :grey, marker = '◼', markersize = 30,strokecolor = :black,strokewidth=2)]
    elem_3 = [MarkerElement(color = :white, marker = '◼', markersize = 30,strokecolor = :black,strokewidth=2)]

    Legend(fig[2, 5:6, Bottom()],
    [elem_1,elem_2,elem_3, ],
    [L"|-1⟩",L"|0⟩",L"|1⟩"],
    patchsize = (35, 35), rowgap = 10,nbanks=3,backgroundcolor = (:black,0.05),colgap=20,padding = (2,2,-2,-2),patchlabelgap=1)

    rowsize!(fig.layout, 1, Relative(0.25))
    colsize!(botrow, 1, Relative(0.4))
    colsize!(fig.layout, 5, Relative(0.2))
    colgap!(fig.layout, 4, 10)
    colgap!(SideRow, 1, 10)
    rowgap!(fig.layout, 1, 15)
    Label(fig[1, 1,TopLeft()], L"(a)$$", fontsize = 22, padding = (0, 50, 0, 0))
    Label(fig[1, 5,BottomRight()], "⤵", fontsize = 22,padding = (0, 210, 30 -5, 0),tellheight = false,tellwidth = false,rotation = 1.5pi)
    Label(fig[1, 5,BottomRight()], "⤴", fontsize = 22,padding = (0, 203, 40 -5, 0),tellheight = false,tellwidth = false,rotation = .0pi)
    
    # Label(fig[1, 5,BottomRight()], L"\mathcal{F}_{\square}", fontsize = 18,padding = (0, 173, 40 -25, 0),tellheight = false,tellwidth = false,rotation = .0pi)
    # Label(fig[2, 1,TopLeft()], L"(b)$$", fontsize = 22)
    Label(fig[2, 1,TopLeft()], L"(b)$$", fontsize = 22,padding = (0, 50, 0, 0))
    Label(fig[2, 5,TopLeft()], L"(c)$$", fontsize = 22, padding = (0, -12, 10, 0))
    save("../../figs/PaperFigs/EnergyScaling_S1_L20.pdf", fig)
    fig
end

# plot_overview(results_df,results_df_others)
plot_overview(results_df,results_df_others)
##
wL_large = 28
results_df_large = collect_results(ENV["MYSCRATCH"]*"/Spiderweb/4x4Comp/EnergyScaling_S1_L$L_large")


let 
    # meanenergies = mean(results_df_large.Energy)

    # fig = scatter(meanenergies)
    fig = Figure(size =  1.0 .* (900, 550), fontsize = 22)
    ax = with_theme(theme_SimpleTicks()) do
        Axis(fig[1,1], xlabel=L"# Sector",xticks = SimpleTicks(1:length(results_df_large.Energy)), 
            ylabel=L"E_0/N_{\text{sites}}",
        )
    end

    sectors = unique(results_df.Sector)
    sectors_others = unique(results_df_others.Sector)
    
    colors = Makie.ColorSchemes.distinguishable_colors(length(sectors)+length(sectors_others), lchoices = LinRange(0, 80, 100))[1:length(sectors)]

    # colors = Makie.ColorSchemes.distinguishable_colors(length(unique(results_df_large.Sector)))
    Snum = eachindex(results_df_large.Energy)

    enmean = mean.(results_df_large.Energy) / L_large^2
    enstd = std.(results_df_large.Energy) / L_large^2

    scatter!(ax,Snum, enmean , label = "GFMC", markersize = 10, marker = '◼', color = colors)
    errorbars!(ax, Snum, enmean, enstd,whiskerwidth = 10, color = colors)

    for (i,en) in enumerate(results_df_large.Energy)
        # scatter!(ax,i*ones(length(en)), en/L_large^2,color = colors[i], )
    end
    # errorbars!(eachindex(meanenergies),meanenergies,std(results_df_large.Energy))
    fig 
end
##
# resultsOld = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S2,CT,10,500,SW.PlaquetteNumberGuidingFunction(0.12);equilibration_steps=0,pre_equilibration_steps=1000) for _ in 1:10])
# resultsOld = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S2,CT,10,500,SW.PlaquetteNumberGuidingFunction(0.12);equilibration_steps=1000,pre_equilibration_steps=100) for _ in 1:10])
a = SW.findMaxFlipConf(copy(S2) .*= 1; numRuns=100, tau=0.1, Nwalkers=40, NSteps=100, ψG = SW.PlaquetteNumberGuidingFunction(0.2),mu = 0.0, pre_equilibration_steps = 10000, scatter_fraction = 1.)