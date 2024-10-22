import Pkg
Pkg.activate(joinpath(@__DIR__,"../"))
cd(joinpath(@__DIR__,"../../"))
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel

include("plottingUtils.jl")
meanstd(x) = (mean(x),std(x))
function upscale(Conf,L)
    S = similar(Conf,L,L)
    per = SW.PeriodicMatrix(Conf,L,L)
    for I in CartesianIndices(S)
        S[I] = per[I]
    end
    S

end
##
function shuffleSector!(S,N)

    flips = (:diag,:anti)
    # flips = (:row,:col,:diag,:anti)
    spin = SW.getSpin(S)
    function transformSpins!(vec,sgn)
        any(==(sgn*spin),vec) && return
        vec .+= sgn
        return
    end

    function applyMove!(S)
        whichflip = rand(flips)
        sgn = -sign(sum(S))
        if sgn == 0 
            sgn = rand((-1,1))
        end
        if whichflip == :row
            i = rand(axes(S,1))
            row = @view S[i,begin+iseven(i):2:end]
            transformSpins!(row,sgn)
        elseif whichflip == :col
            j = rand(axes(S,2))
            col = @view S[begin+iseven(j):2:end,j]
            transformSpins!(col,sgn)
        elseif whichflip == :diag
            i = rand(axes(S,2)[1:2:end])
            diag = SW.getDiagonal(S,i,1,true)
            transformSpins!(diag,sgn)
        elseif whichflip == :anti
            i = rand(axes(S,2)[1:2:end])
            diag = SW.getDiagonal(S,i,-1,true)
            transformSpins!(diag,sgn)
        end
        return S
    end

    for iteration in 1:N
        applyMove!(S)
    end
    while sum(S) != 0
        applyMove!(S)
    end
    S .= 2S
    SW.fulFillsConstraint(S) || error("Constraint not fulfilled")
    return S
end
##
function generatePeriodic(L)
    S = SW.stencilConfig(zeros(L,L),1;
    boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )
    # SW.rand!(S)
    UCSIt = collect(
        Iterators.product((-1:1 for i in 1:L,j in 1:L)...)
    )
    function isGS(UC)
        S[:] .= 2 .*UC
        SW.fulFillsConstraint(S)
    end
    UCS = Matrix{Int8}[]

    for UC in UCSIt
        isGS(UC) || continue
        push!(UCS,copy(parent(parent(S))))
    end
    return UCS
end
a = generatePeriodic(4)

##
function filterConfs(UCs)
    Lx,Ly = size(UCs[1])
    S = SW.stencilConfig(zeros(Lx,Ly),1;
    boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )
    S_minus = copy(S)
    S_flip = copy(S)
    S_xflip = copy(S)
    S_yflip = copy(S)
    SUC = copy(parent(parent(S)))

    aSet = Set(empty(a))
    moves = empty!([(1,1,1)])
    function isUnique(UC)
        # S_periodic = SW.PeriodicMatrix(UC,Lx,Ly)
        S .= UC
        S_minus .= -UC
        S_flip .= UC'
        S_xflip .= @view UC[end:-1:1,:]
        S_yflip .= @view UC[:,end:-1:1]

        for S′ in (S,S_minus,S_flip,S_xflip,S_yflip)
            for T_x in 0:Lx-1, T_y in 0:Ly-1

                # SUC = parent(parent(S′))
                SUC_per = SW.PeriodicMatrix(parent(parent(S′)))
                # @info "" size(SUC_per[T_x:T_x+Lx-1,T_y:T_y+Ly-1])
                SUC .= @view SUC_per[T_x:T_x+Lx-1,T_y:T_y+Ly-1]

                if SUC in aSet
                    return false
                end
                if SUC' in aSet
                    return false
                end
                if -SUC in aSet
                    return false
                end
                if -SUC' in aSet
                    return false
                end
                SW.getMoves!(moves,S′)
                for (i,j,sgn) in moves
                    SW.applyPlaquette!(S′,i,j,sgn)
                    if SUC in aSet
                        return false
                    end
                    SW.applyPlaquette!(S′,i,j,-sgn)
                end
            end
        end
        return true
    end
    for UC in UCs
        isUnique(UC) || continue
        push!(aSet,UC)
    end
    return aSet
    
end
a_st = sort!(collect(filterConfs(a)),by=x->sum(abs,x))
##
function makeConf(UC,L)
    S = SW.stencilConfig(zeros(L,L),1;
    boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )
    S .= SW.getPeriodicState(UC,L,L)
    return S
end  

a_configs = let
    L = 4
    confs = [makeConf(UC,L) for UC in a_st]

    # filter!(SW.fulFillsConstraint, confs)
    filter!(x->length(SW.getApplicablePlaquettes(x)) > 0,confs)
    @assert all(SW.fulFillsConstraint,confs)
    confs
end
##

function getMaxFlipConfs(Configs;Nwalkers = 28*2,NSteps = 1000)
    filterConfs = Set(empty(Configs))
    mu = -0.8
    CT = SW.ContinuousTimeMethod(0.5,Hxx = SW.Hxx_RK(mu))
    ψG = SW.PlaquetteNumberGuidingFunction(0.5)
    Sbuff = copy(parent(parent(Configs[1])))
    ConfBuff = copy(Configs[1])
    moves = empty!([(1,1,1)])

    for S in Configs
        res = SW.startManyWalkerGFMC(S,CT,Nwalkers,NSteps,ψG)
        isNew = true
        maxMoves = 0
        maxConf = copy(ConfBuff)
        for c in eachslice(res.SaveConfigs,dims=(3,4))
            Sbuff .= c
            if Sbuff in filterConfs || -Sbuff in filterConfs
                isNew = false
                break
            end
            ConfBuff .= Sbuff
            SW.getMoves!(moves,ConfBuff)
            if length(moves) > maxMoves
                maxConf .= ConfBuff
                maxMoves = length(moves)
            end
        end
        if isNew
            push!(filterConfs,maxConf)
        end
    end
    filterConfs_arr = collect(filterConfs)

    sort!(filterConfs_arr,by = x->length(SW.getMoves!(moves,x)),rev = true)
end

reducedConfigs = getMaxFlipConfs(a_configs;Nwalkers = 28,NSteps = 1) #first reduction
##
reducedConfigs = getMaxFlipConfs(reducedConfigs;Nwalkers = 28,NSteps = 10)
##
reducedConfigs = getMaxFlipConfs(reducedConfigs;Nwalkers = 28,NSteps = 100)
##
reducedConfigs = makeConf.(filterConfs(parent.(parent.(reducedConfigs))),4)
reducedConfigs = getMaxFlipConfs(reducedConfigs;Nwalkers = 28*1,NSteps = 1000)
##
function plotConfs(Confs)
    numGS = length(Confs)
    with_theme(theme_SimpleTicks()) do
        fig = Figure(fontsize = 22,size = length(Confs) .*(12,12))
        axs = [Axis(fig[i,j],xlabelvisible=false,ylabelvisible=false,xticklabelsvisible=false,yticklabelsvisible=false,aspect=1)
        for i in 1:round(Int,sqrt(numGS),RoundUp),j in 1:round(Int,sqrt(numGS))]
        for i  in eachindex(Confs)
            ax = axs[i]
            # newConf = similar(Confs[i],4,4)
            # newConf[1:4,1:4] .= @view Confs[i][1:4,1:4]
            SW.plotApplPlaquettes!(ax,Confs[i])
            try
                rowgap!(fig.layout,i,0.)
                colgap!(fig.layout,i,0.)
            catch
            end
        end
        return fig
    end 
end
plotConfs(reducedConfigs)


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

function findEnergies(Configs,CT,ψG;Nwalkers = 28*2,NSteps = 2000,equilibration_steps = NSteps ÷8)
    en = zeros(length(Configs))
    Δen = zeros(length(Configs))

    Threads.@threads for i in eachindex(Configs,en)
        S = Configs[i]
        if length(SW.getApplicablePlaquettes(S)) == 0
            en[i] = 0
            Δen[i] = 0
            continue
        end
        results = [SW.startManyWalkerGFMC(S,CT,Nwalkers,NSteps,ψG;equilibration_steps,pre_equilibration_steps=NSteps) for _ in 1:6]

        energies = SW.getEnergies.(results,1,min(NSteps÷3, 300))
        energiesMean = mean.(energies)
        energiesStd = std.(energies)
        e0_index = findFirstMiIndex(energiesMean)
        en[i] = energiesMean[e0_index]
        Δen[i] = energiesStd[e0_index]
    end
    return (;en,Δen)
end
##

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