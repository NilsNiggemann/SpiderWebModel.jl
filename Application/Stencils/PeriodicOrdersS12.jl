
##

a = generatePeriodic(6,1/2)
##
a_st = sort!(collect(filterConfs(a,1/2)),by=x->sum(abs,x))
##

a_configs = let
    L = 5
    confs = [makeConf(UC,L,1/2) for UC in a_st]

    # filter!(SW.fulFillsConstraint, confs)
    filter!(x->length(SW.getApplicablePlaquettes(x)) > 0,confs)
    @assert all(SW.fulFillsConstraint,confs)
    confs
end
##


reducedConfigs = getMaxFlipConfs(a_configs;Nwalkers = 28,NSteps = 1) #first reduction
##
reducedConfigs = getMaxFlipConfs(reducedConfigs;Nwalkers = 28,NSteps = 10)
##
reducedConfigs = getMaxFlipConfs(reducedConfigs;Nwalkers = 28,NSteps = 100)
##
##
function findFirstMiIndex(arr)
    minval = arr[begin]
    i_min = firstindex(arr)
    for (i,x) in enumerate(arr)
        if x < minval
            i_min = i
            minval = x
        end
    end
    return i_min
end


μ = 0.95

CT = SW.ContinuousTimeMethod(0.1,1,-0.266length(a_configs[1]),SW.Hxx_RK(μ))
ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))

##
@time res = findEnergies(upscale.(reducedConfigs,12),CT,ψG;Nwalkers = 28*2,NSteps = 1000)
##
with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"# Config $$",ylabel = L"E_0/N_\text{sites}")
    Nsites = 4*4
    errorbars!(eachindex(res.en),res.en ./Nsites,res.Δen ./Nsites,whiskerwidth=4,color = :black)
    scatter!(res.en ./Nsites,marker = '×',color = :black)
    fig
end
##
Allres = empty!([res])
mus_sectors = -0.1:0.05:0.99
for mu in mus_sectors
    CT = SW.ContinuousTimeMethod(0.1,1,-0.266length(a_configs[1]),SW.Hxx_RK(mu))
    ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-mu))
    res = findEnergies(upscale.(reducedConfigs,8),CT,ψG;Nwalkers = 28*2,NSteps = 800)
    push!(Allres,res)
end

##
with_theme(theme_SimpleTicks()) do 
    fig = Figure(fontsize = 22)
    E_scal = 0.2
    ens = getproperty.(Allres,:en)
    Δens = getproperty.(Allres,:Δen)
    mus = mus_sectors
    Elin = E_scal .* (mus .-1)

    ax = Axis(fig[2,1:3],xlabel = L"\mu",ylabel = L"E_0/(N_\text{sites}(1 - \mu ))")
    Nsites = 12*12

 
    minsectors = argmin.(ens)

    for i_sector in eachindex(reducedConfigs)
        en = getindex.(ens,i_sector) ./ Nsites ./(1 .-mus)
        Δen = getindex.(Δens,i_sector)./(1 .-mus) ./ Nsites
        scatterlines!(ax,mus,en ,marker = '×',linewidth = 2)
        errorbars!(ax,mus,en ,Δen,whiskerwidth=4)
    end

    for i in eachindex(mus)[3:end]
        if minsectors[i-1] != minsectors[i]
            vlines!(ax,mus[i-1],color = :black,linestyle = :dash,linewidth = 0.8)
        end
    end
    
    ax2 = Axis(fig[1,1],aspect=1,xticks = 1:2:4,yticks = 1:2:4)
    ax3 = Axis(fig[1,2],aspect=1,xticks = 1:2:4,yticks = 1:2:4)
    ax4 = Axis(fig[1,3],aspect=1,xticks = 1:2:4,yticks = 1:2:4)
    rowsize!(fig.layout,1,Relative(0.3))
    SW.plotApplPlaquettes!(ax2,reducedConfigs[1])
    SW.plotApplPlaquettes!(ax3,reducedConfigs[5])
    SW.plotApplPlaquettes!(ax4,reducedConfigs[6])
    fig
end
##
# perm = sortperm(res.en)
# # SW.plotApplPlaquettes(a_configs[perm][3])

# emin = minimum(res.en)
# numGS = findfirst(>(emin+3e-1),res.en[perm])

##

##
μ = -0.1

CT = SW.ContinuousTimeMethod(0.1,1,-0.266length(a_configs[1]),SW.Hxx_RK(μ))
ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))


@time results1 = [SW.startManyWalkerGFMC(upscale(reducedConfigs[end-5],8),CT,28*6,5000,ψG) for _ in 1:10]
##
@time results2 = [SW.startManyWalkerGFMC(upscale(reducedConfigs[end-6],8),CT,28*6,5000,ψG) for _ in 1:10]
##
plotEnergies(results1,CT;normalize=true,dense=true,τ = 10,color = :black)
plotEnergies!(results2,CT;normalize=true,dense=true,τ = 10,color = :red)
current_figure()