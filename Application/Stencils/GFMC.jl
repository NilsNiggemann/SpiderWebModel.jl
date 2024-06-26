import Pkg
Pkg.activate(joinpath(@__DIR__ ,"../"))
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel
##
S = SW.stencilConfig(parent(SW.getStairCase(12)),1/2)
# S = SW.stencilConfig(SW.constructConfigPath(15,15,SW.ALLGS_S12),1/2)
HStair = SW.generateHilbertSpace(SW.SpinConfig(S))
##
# HTest = Array(HStair.H)
ExSol = SW.SolveHKrylov(HStair.H)
E0 = ExSol.values[1]

getRKWavefunction(ψ) = SW.normalize!(one.(ψ))

v0 = ExSol.vectors[1]
# v0 = getRKWavefunction(ExSol.vectors[1])
HConfs = SW.spinConfig.(HStair.AllStates,Ref(SW.SpinConfig(S)),Ref(HStair.plaqMapping))
##
function constructExactGuidingFunc(v0,AllStates)
    AllSTDict = Dict(SW.stencilConfig(parent(s),1/2)=>i for (i,s) in enumerate(AllStates))
    function psiG(Conf)
        ind = AllSTDict[Conf]
        return v0[ind] + 0.01
    end
end

function makeBlocks(arr;blocksize = nothing, numBlocks = 32)
    if blocksize === nothing
        blocksize = length(arr) ÷ numBlocks
    end
    blockedArr = collect(SW.splitIntoBins(arr,blocksize))
    
    minsize,maxsize = extrema(length.(blockedArr))
    
    if maxsize!=minsize
        pop!(blockedArr)
    end
    return blockedArr
end


##
magEx = SW.getMagnetization(HConfs, v0)
@time SqEx = SW.getStructureFac(HConfs,v0)
# @time SqEx = SW.getEqualWeightStructureFac(HConfs)
# @time SqEx = SW.getStructureFac(HConfs,SW.normalize!(ones(length(v0))))
#___________ManyWalkers_______________________
##
S = SW.stencilConfig(parent(SW.getStairCase(12)),1/2)
# ψG = SW.PlaquetteNumberGuidingFunction(0.197)
ψG = SW.PlaquetteNumberGuidingFunction(0.197)
# ψG = SW.fullVariationalFunction(S,0.197)

nThermal = 5000
# results = [SW.startManyWalkerGFMC(S,2,55_000,3,nThermal,SW.ConstructVaritationalFunc(0.197,S),0) for _ in 1:35]
nBra = 1
##
SW.Random.seed!(1234)
DT = SW.DiscreteTimeMethod(1.,3,1.)
@time resultsDT = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,DT,30,5000,ψG,equilibration_steps=nThermal) for i in 1:12])
##
SW.Random.seed!(1234)
CT = SW.ContinuousTimeMethod(0.01,10,14.,SW.Hxx_zero())

@time resultsCT = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,80,8000,ψG,equilibration_steps=nThermal) for i in 1:12])
# @time SW.startManyWalkerGFMC(S,2,3000,nBra,ψG,1;method = SW.ContinuousTimeMethod(),equilibration_steps=nThermal) 

##
# plotEnergies(resultsCT,nBra,E0,p=200,Emin = 1.05*E0,Emax = 0.9*E0)
# plotEnergies(resultsCT,1)

plotEnergies(resultsDT,DT.nBranch,color =:blue,E0,Emin=NaN,Emax=NaN,p=150,normalize=true)
plotEnergies!(resultsCT,CT.nBranch,E0,Emin=NaN,Emax=NaN,p=150,)
current_figure()

##

# outfile = "Data/temp/S12/"
# rm(outfile;recursive=true,force=true)
# mkpath(outfile)
# @time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,16,50000÷nBra,nBra,ψG,1;equilibration_steps=nThermal,outfile =string(outfile,i,".h5")) for i in 1:8])
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,1,100000,nBra,x->1,1;equilibration_steps=nThermal) for i in 1:42])
# @time SW.startManyWalkerGFMC(S,16,nThermal+500_00÷nBra,nBra,ψG,1;outfile = string(outfile,1,".h5"))
##
plotEnergies(results,nBra,E0,p=50,Emin = 1.02*E0,Emax = 0.93*E0)
##
@time results2 = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,1,100000,nBra,ψG,1;equilibration_steps=nThermal) for i in 1:42])
##
plotEnergies!(results2,nBra,E0,p=50,Emin = 1.02*E0,Emax = 0.93*E0,color = :blue)
current_figure()
##
plotEnergies(results2,1,E0,p=200,Emin = 1.002*E0,Emax = 0.99*E0,color = :blue)
##
nBra = 3
@time results3 = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,30,100000÷30÷nBra,nBra,ψG,1;equilibration_steps=nThermal) for i in 1:42])
##
plotEnergies!(results3,nBra,E0,p=200,Emin = 1.002*E0,Emax = 0.99*E0,color = :darkred)
current_figure()

##
#___________Observables_______________________
##
# getAvgMag(Conf) = sum(Conf) ./ (2*length(Conf))
function getMag(results,pmax,I)
    # I = CartesianIndex(I)
    Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights,1,pmax) for res in results]
    mz_func(Conf) = Conf[I] /2

    @views getObs(p) = [SW.getObs(Gnp[:,1:p],res.SaveConfigs,res.reconfigurationTable,mz_func,p÷2) for (res,Gnp) in zip(results,Gnps)]
    obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
end

function getMagRK(results,p,I)
    # I = CartesianIndex(I)
    TotWeights = [SW.mean(res.TotalWeights) for res in results]
    weights = [abs2.(res.TotalWeights ./ W) for (res,W) in zip(results,TotWeights)]
    # SW.normalize!.(weights,1)
    i,j = Tuple(I)
    @views getObs(ws,res) = SW.mean(w * mag/2 for (w,mag) in zip(ws,res.SaveConfigs[i,j,1,:]))

    obs = fetch.([Threads.@spawn getObs(ws,res) for (ws,res) in zip(weights,results)])
    return [obs for i in 1:p]
end

function getSiSj(results,pmax,I,J)
    # I = CartesianIndex(I)
    Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights,1,pmax) for res in results]

    getCorr(Conf) = (Conf[I]*Conf[J]) *0.25

    @views getObs(p) = [SW.getObs(Gnp[:,1:p],res.SaveConfigs,res.reconfigurationTable,getCorr,p÷2) for (res,Gnp) in zip(results,Gnps)]
    obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
end
##
results = resultsCT
nBra = CT.nBranch
mags1 = getMag(results,100÷nBra,CartesianIndex(3,3)) 
mags2 = getMag(results,100÷nBra,CartesianIndex(4,1)) 
mags3 = getMag(results,100÷nBra,CartesianIndex(2,3)) 
##
IJ_SS = (CartesianIndex(3,4),CartesianIndex(4,2))
SiSj = getSiSj(results,100÷nBra,IJ_SS[1],IJ_SS[2]) 
SiSjex = SW.getSij(HConfs,v0,IJ_SS[1],IJ_SS[2])
##

function errorBarLegend(size = 0.5;linekwargs = (;),markerkwargs = (;))
    center = 0.5
    ymin = center - size/2
    ymax = center + size/2

    errorbar = [Point2f(center, ymin), Point2f(center, ymax)]
    [
        LineElement(linepoints = errorbar),
        MarkerElement(points = errorbar, marker = :hline, markersize = 10;linekwargs...),
        MarkerElement(points = [Point2f(center, center),], marker = '●', markersize = 7;markerkwargs...)
    ]
end

with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"projection order $$",xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5))
    proj = nBra .*eachindex(mags1)
    Cols = (:blue,:green,:orange)
    Is = ((3,3),(4,1),(2,3))
    for (i,mags) in enumerate((mags1,mags2,mags3))
        mExact = magEx[Is[i]...]

        sgn = sign(mExact)
        m =  sgn .* mean.(mags)
        scatter!(ax,proj,m,label = L"GFMC$$",color = Cols[i], marker = '●',markersize = 5)
        err =sqrt.(var.(mags))
        errorbars!(ax,proj,m,err,whiskerwidth = 3.5,color = Cols[i])
        
        hlines!([sgn * mExact],color = Cols[i],label = L"exact $$")
        sgstring = sgn > 0 ? "" : "-"
        text!(ax,Point2(proj[end÷2],mExact-0.015),color = Cols[i],text = L"%$sgstring \langle S^z_{%$(Is[i])}\rangle")
    end
    chi = mean.(SiSj)

    chiScale = -4
    scatter!(ax,proj,chiScale .*chi,label = L"GFMC$$",color = :black, marker = '●',markersize = 5)
    err =sqrt.(var.(SiSj))
    errorbars!(ax,proj,chiScale .*chi,chiScale .*err,whiskerwidth = 3.5,color =  :black)
    hlines!([chiScale*SiSjex],color = :black,label = L"exact $$")

    i,j = Tuple.(IJ_SS)
    text!(ax,Point2(proj[end÷2],chiScale*SiSjex-0.015),color = :black,text = L"%$chiScale \times \langle S^z_{%$i} S^z_{%$j}\rangle",align = (:center,:center))

    return fig

    Legend(fig[1, 1], [[
        errorBarLegend(0.6)
    ],
    [LineElement(linepoints = [Point2f(0., 0.5), Point2f(1, 0.5)])]], [L"GFMC$$",L"exact $$"], tellheight = false, tellwidth = false,halign = :right, valign = :center,margin = (10,10,120,10))


    # axislegend(ax,merge=true,position = :rc,unique=true)
    # xlims!(ax,0.5,last(proj))
    # ylims!(ax,E0-1e-2,E0+3e-2)
    # save("Application/exactFig/GFMCEnergy.png",fig)
    fig
end
##
@views function getFullMag(res,p)
    # I = CartesianIndex(I)
    Gnp = SW.precomputeNormalizedAccWeight(res.TotalWeights,1,p)
    buffer = zeros(size(res.SaveConfigs[:,:,begin,begin]))
    function m(Conf)
        buffer .= Conf
    end
    res = SW.getObs(Gnp,res.SaveConfigs,res.reconfigurationTable,m,p÷2)
end
magFull = fetch.([Threads.@spawn getFullMag(res,20÷nBra) for res in results])

##
with_theme(theme_SimpleTicks()) do 
    m = mean(magFull) ./2
    fig,ax,hm = heatmap(m,colormap = :greys,axis=(;aspect=1,title = L"GFMC$$"),figure = (;size = 0.8 .*(360,500)))
    ax2 = Axis(fig[2,1],aspect=1,title = L"exact $$")
    heatmap!(ax2,magEx,colormap = :greys,colorrange = extrema(m))
    Colorbar(fig[1:2,2],hm,label = L"\langle S^z_i \rangle")
    fig
end
##

function getSq(res,p)
    @views Gnp = SW.precomputeNormalizedAccWeight(res.TotalWeights,1,p)
    # @views Gnp = ones(length(res.TotalWeights[nThermal:end]),p)

    Conf = @view res.SaveConfigs[:,:,begin,begin]
    NSites = length(Conf)
    Sq = similar(Conf, ComplexF64)
    
    Si = similar(Conf, ComplexF64)

    plan = SW.FFTW.plan_fft(Si)

    function SqFunc(Conf)
        Si.= Conf
        SW.mul!(Sq, plan, Si)

        Sq .= abs2.(Sq)
    end

    @views resSq = SW.getObs(Gnp,res.SaveConfigs,res.reconfigurationTable,SqFunc,p÷2)
    
    newRes = similar(resSq,Float64,size(resSq).+1)
    newRes[begin:end-1,begin:end-1] .= real.(resSq)

    @views newRes[end,begin:end] .= newRes[begin,:]
    @views newRes[begin:end,end] .= newRes[:,begin]
    newRes ./NSites
    
    # obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
end

##
SqsGFMC = fetch.([Threads.@spawn getSq(res,100÷nBra) for res in results])
##
with_theme(theme_PiTicks()) do 
    Sq = mean(real(SqsGFMC)) ./4
    kx = ky = 2pi .* LinRange(0,1,size(Sq,1))
    fig,ax,hm = heatmap(kx,ky,Sq,colormap = :viridis,axis=(;aspect=1,title = L"GFMC$$"),figure = (;size = (360,500)))
    ax2 = Axis(fig[2,1],aspect=1,title = L"exact $$")
    heatmap!(ax2,kx,ky,real(SqEx.Sq),colormap = :viridis,colorrange = extrema(Sq))
    Colorbar(fig[1:2,2],hm,label = L"\langle \mathcal{S}^{zz}(\textbf{q})\rangle")
    fig
end
##
#___________Spin-1_______________________

ψG = SW.PlaquetteNumberGuidingFunction(0.13)
nThermal = 1
SW.Random.seed!(1234)
S = SW.stencilConfig(zeros(20,20),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
# CT = SW.ContinuousTimeMethod(0.003,4,10.,SW.Hxx_SIA(0.0))
DT = SW.DiscreteTimeMethod(0.,4,67)
@time resultsDT = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,DT,30,2_000,ψG,equilibration_steps=0,pre_equilibration_steps=10_000,scatter_fraction=0.92) for i in 1:6])
##
results = resultsDT
nBra = DT.nBranch
plotEnergies(resultsDT,round(Int,DT.nBranch),p=200,nThermal=1000,label = L"discrete time$$")

##
# SqsGFMC = fetch.([Threads.@spawn (
# getSq(res,150÷nBra,1,1) .+ getSq(res,150÷nBra,2,2) .+ getSq(res,150÷nBra,1,2) .+ getSq(res,150÷nBra,2,1)) for res in results])
SqsGFMC = fetch.([Threads.@spawn getSq(res,500÷nBra) for res in results])
##
function SqFieldTheory(x,y)
    num = cos(x) - cos(y) +2sin(x)sin(y) 
    denom = (cos(x) - cos(y))^2 + (2sin(x)sin(y))^2
    return num^2/(sqrt(denom)+1e-30)
end

function SqFieldTheory2(kx,ky,k,b2)
    (sqrt(2)*sqrt(-((-4 + 4*cos(kx)*cos(ky) + cos(2*kx)*(1 - 2*cos(2*ky)) + cos(2*ky))*(40*b2 + k + 8*b2*(cos(2*kx) + 2*cos(kx)*(-4*cos(ky) + cos(kx)*cos(2*ky))))))*   (cos(kx) - cos(ky) + 2*sin(kx)*sin(ky))^2)/((cos(kx) - cos(ky))^2 + 4*sin(kx)^2*sin(ky)^2+1e-30)
end

function makeSqFTPlots(SqsGFMC,k=1,b2=0)
    SqFT_func(x,y)  = SqFieldTheory2(x,y,k,b2)
    with_theme(theme_PiTicks()) do 
        # Sq = sqrt.(var(real(SqsGFMC))) ./4
        Sq = mean(SqsGFMC) ./4
        kx = ky = 2pi .* LinRange(0,1,size(Sq,1))
        fig = Figure(fontsize = 22,size = (800,400))
        axMC = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"k_y",title = L"GFMC$$",aspect = 1)
        axerr = Axis(fig[1,2],xlabel = L"k_x",ylabel = L"k_y",title = L"std error$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

        axFT = Axis(fig[1,3],xlabel = L"k_x",ylabel = L"k_y",title = L"U(1) theory$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

        # axDiff = Axis(fig[1,4],xlabel = L"k_x",ylabel = L"k_y",title = L"Difference$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

        # err = abs.(Sq .- SqFieldTheory.(kx,kx'))
        err = sqrt.(var(SqsGFMC)) ./4
        hmMC = heatmap!(axMC,kx,ky,Sq,colormap = :viridis)
        SqFT = [SqFT_func(x,y) for x in kx, y in ky]
        hmFT = heatmap!(axFT,kx,ky,SqFT,colormap = :viridis)
        # heatmap!(axerr,kx,ky,err,colormap = :viridis,colorrange = extrema(!isnan,Sq))
        hmerr = heatmap!(axerr,kx,ky,err,colormap = :viridis)
        # heatmap!(axDiff,kx,ky,(Sq ./maximum(Sq)) .- (SqFT ./maximum(SqFT)),colormap = :viridis)

        Colorbar(fig[2,1],hmMC,label = L" \mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks())
        Colorbar(fig[2,2],hmerr,label = L"\sigma( \mathcal{S}^{zz}(\textbf{q}))",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks())
        Colorbar(fig[2,3],hmFT,label = L" \mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8), width = Relative(0.8),vertical=false,ticks = SimpleTicks())

        rowsize!(fig.layout,2,Relative(0.1))
        fig
    end
end
makeSqFTPlots(SqsGFMC,1,0.1)
##

projection_orders = [800,500]
SqsGFMC_p = [ fetch.([Threads.@spawn getSq(res,p÷nBra) for res in results]) for p in projection_orders]
##
# with_theme(merge(theme_PiTicks(),theme_dark())) do 
function plotCut(SqsGFMC_p,k=1,b2=0)
    SqFT_func(x,y)  = SqFieldTheory2(x,y,k,b2)

    with_theme(theme_PiTicks()) do 
        fig = Figure(fontsize = 18,size = (500,400))
        ax = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"\mathcal{S}(\mathbf{q})",title = L"GFMC$$",yticks = SimpleTicks())
        axquad = Axis(fig[2,1],xlabel = L"k_x^4",ylabel = L"\mathcal{S}(\mathbf{q})",yticks = SimpleTicks(),xticks = SimpleTicks())

        kx = 2pi .* LinRange(0,1,size(S,1).+1)
        
        kx1, zoomindex = findmin(x->abs(x-0.2pi),kx)
        kzoom = kx[begin:zoomindex]
        for (p,Sqs) in zip(projection_orders,SqsGFMC_p)
            Sqs ./= maximum(mean(Sqs))

            Sq = mean(Sqs)[1,:] #./4
            Sqerr = sqrt.(var(Sqs))[1,:] #./4
            # Sqerr ./= maximum(mean(Sqs))
            lines!(ax,kx,Sq,
            # color = p,colorrange = extrema(projection_orders),colormap = :viridis,linewidth = 3,
            label = L"$p=%$(p)$"
            )
            errorbars!(ax,kx,Sq,Sqerr,
            whiskerwidth = 3.5,
            # color = p,colorrange = extrema(projection_orders),colormap = :viridis
            )
            scatterlines!(axquad,kzoom.^4,Sq[begin:zoomindex],marker = '□',markersize = 8)
            errorbars!(axquad,kzoom.^4,Sq[begin:zoomindex],Sqerr[begin:zoomindex],whiskerwidth = 3.5)
        end
        
        kxFT = 2pi .* LinRange(0,1,300)
        SqFT = SqFT_func.(kxFT,0)
        SqFTFull = [SqFT_func(x,y) for x in kxFT, y in kxFT]
        SqFT ./= maximum(SqFTFull)
        lines!(ax,kxFT,SqFT,color = :black,linewidth = 2,linestyle = :dash,label = L"U(1)$$")

        _, zoomindex = findmin(x->abs(x-kx[zoomindex]),kxFT)
        kzoom = kxFT[begin:zoomindex]

        lines!(axquad,kzoom.^4,SqFT[begin:zoomindex],color = :black,linewidth = 2,linestyle = :dash,label = L"U(1)$$")

        lines!(axquad,kzoom.^4,kzoom .^4*1.15SqFT[zoomindex],color = :red,linestyle = :solid,label = L"$\simk^4$",linewidth = 0.8)
        axislegend(ax,merge=true)
        fig
    end
end
plotCut(SqsGFMC_p,1,1)
##

##
SW.Random.seed!(1234)
S = SW.stencilConfig(zeros(24,24),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
DT = SW.DiscreteTimeMethod(0.,6,size(S,1)^2/4)
ψG = SW.fullVariationalFunction(S,0.168)
stochReconfRes = SW.stochastic_reconfiguration(S,DT,i->round(Int,1000+ 30*i),ψG,10,i->max(2. - 1*log(i),0.1) ,SW.IterativeSRSolver();Nwalkers = 6*30,rel_tolerance=1e-8,equilibration_steps=100,pre_equilibration_steps=50_000,scatter_fraction=0.6)
ψG = SW.FullVariationalGuidingFunction(stochReconfRes.params_steps[:,:,end-2])
##
DT = SW.DiscreteTimeMethod(0.,12,size(S,1)^2/3.8)
scatter_fraction = 0.98
@time resultsDT = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,DT,40,10_000,ψG;equilibration_steps=0,pre_equilibration_steps=50_000,scatter_fraction) for i in 1:6])
##
# results = resultsDT
nBra = DT.nBranch
plotEnergies(resultsDT,DT.nBranch,p=500,nThermal=100,label = L"Continuous time$$")

##
equilib_plots(resultsDT;scatter_fraction,averageSteps=10,Ntrack=30,p = 300)
##
# U field repulsion

SW.Random.seed!(1234)
S = SW.stencilConfig(zeros(16,16),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)

CT = SW.ContinuousTimeMethod(0.01,3,10.,SW.Hxx_SIA(0.3))
ψG = SW.fullVariationalFunction(S,.2)
stochReconfRes = SW.stochastic_reconfiguration(S,CT,i->round(Int,2000+ 200*i),ψG,10,i->max(1. - 0.1*log(i),0.3) ,SW.IterativeSRSolver();Nwalkers = 6*50,rel_tolerance=1e-8,equilibration_steps=100,pre_equilibration_steps=50_000,scatter_fraction=0.6)
ψG = SW.FullVariationalGuidingFunction(stochReconfRes.params)
##
# ψG = SW.fullVariationalFunction(S,0.15)
# CT = SW.ContinuousTimeMethod(0.03,1,0.2655*prod(size(S)),SW.Hxx_zero())
CT = SW.ContinuousTimeMethod(0.03,1,23.2,SW.Hxx_SIA(0.3))
scatter_fraction = 0.8
@time resultsCT = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,500,10000,ψG;equilibration_steps=0,pre_equilibration_steps=50_000,scatter_fraction) for i in 1:6])
##
# results = resultsCT
plotEnergies(resultsCT,CT,nThermal=1,τ=2,normalize=true)
##
equilib_plots(resultsCT;scatter_fraction,averageSteps=10,Ntrack=80,p = 100)

##
SqsGFMC = fetch.([Threads.@spawn getSq(res,100÷nBra) for res in resultsCT])
makeSqFTPlots(SqsGFMC)
##
projection_orders = [800,300,100]
SqsGFMC_p = [ fetch.([Threads.@spawn getSq(res,p÷nBra) for res in resultsCT]) for p in projection_orders]
plotCut(SqsGFMC_p)