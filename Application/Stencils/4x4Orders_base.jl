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
function generatePeriodic(L,Spin)
    S = SW.stencilConfig(zeros(L,L),Spin;
    boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )
    # SW.rand!(S)
    UCSIt = 
        Iterators.product((-Spin:Spin for i in 1:L,j in 1:L)...)
    
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

##
function filterConfs(UCs,Spin)
    Lx,Ly = size(UCs[1])
    S = SW.stencilConfig(zeros(Lx,Ly),Spin;
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
function makeConf(UC,L,Spin)
    S = SW.stencilConfig(zeros(L,L),Spin;
    boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )
    S .= SW.getPeriodicState(UC,L,L)
    return S
end  
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
function plotConfs(Confs)
    numGS = length(Confs)
    with_theme(theme_SimpleTicks()) do
        fig = Figure(fontsize = 22,size = length(Confs) .*(12,12))
        axs = [Axis(fig[i,j]; SW.getConfigAxis(Confs[1])...,xlabelvisible=false,ylabelvisible=false,xticklabelsvisible=false,yticklabelsvisible=false,aspect=1,xticks = 1:4,yticks = 1:4)
        for i in 1:round(Int,sqrt(numGS),RoundUp),j in 1:round(Int,sqrt(numGS))]

        for i  in eachindex(Confs)
            ax = axs[i]
            # newConf = similar(Confs[i],4,4)
            # newConf[1:4,1:4] .= @view Confs[i][1:4,1:4]
            SW.plotApplPlaquettes!(ax,Confs[i])
            text!(Point(2,2),"$i",fontsize = 100,color = :green,)
            try
                rowgap!(fig.layout,i,0.)
                colgap!(fig.layout,i,0.)
            catch
            end
        end
        return fig
    end 
end


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

