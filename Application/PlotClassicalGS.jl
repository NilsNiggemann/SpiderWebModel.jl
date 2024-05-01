import SpiderWebModel as SW
using CairoMakie
using MakieHelpers
using HDF5

##
files = readdir("Application/ConfsRaw/GurobiConfs",join=true)

Configs = 
    let 
    Confs = Vector{Matrix{Bool}}()

    for file in files
        h5open(file,"r") do f
            for k in keys(f["Confs"])
                confs = read(f["Confs/$k"])
                for conf in eachslice(confs,dims=3)
                    push!(Confs, collect(conf))
                end
            end
        end
    end
    Spinconfs = SW.floatSpinConfig.(Confs,0.5)
    @assert all(SW.fulFillsConstraint,Spinconfs)
    Spinconfs
end
##

SW.plotFractons(Configs[10])
##
function SqLargeN(qx,qy)
    cx = cos(qx)
    cy = cos(qy)
    sx = sin(qx)
    sy = sin(qy)
    c2x = cos(2*qx)
    c2y = cos(2*qy)

    num = 2*(cx -cy +2sx*sy)^2
    denom = 4 - 4*cx*cy - c2x*(1-2c2y) - c2y
    return num/denom
end
##
Spin1Configs = SW.SpinConfig.(eachslice(h5read("Application/ConfsRaw/Spin1Fluc_5.h5","Confs"),dims=3),1)

Sq = SW.getEqualWeightStructureFac(Configs)
SqSpin1 = SW.getEqualWeightStructureFac(Spin1Configs)

with_theme(theme_PiTicks()) do
    fig = Figure(size = 0.8 .*(750,300))
    ticks = ([pi,2pi],[L"$\pi$",L"$2\pi$"])

    ax = Axis(fig[1, 1], aspect = 1,xticks = ticks,yticks = ticks,xminorticksvisible = true ,xlabel = L"q_x",ylabel = L"q_y",yminorticksvisible = true,title = L"Large-N$$")
    ax2 = Axis(fig[1, 2], aspect = 1,yticklabelsvisible=false,xminorticksvisible = true ,xlabel = L"q_x",yminorticksvisible = true,xticks = ticks,yticks = ticks,title = L"class. spin-$1/2$")
    ax3 = Axis(fig[1, 3], aspect = 1,yticklabelsvisible=false,xminorticksvisible = true ,xlabel = L"q_x",yminorticksvisible = true,xticks = ticks,yticks = ticks,title = L"class. spin-$1$")

    SqFunc = SW.getSqCont(real(Sq.Sq)[1:end-1,1:end-1])

    Sq1Func = SW.getSqCont(real(SqSpin1.Sq)[1:end-1,1:end-1] ./2) # normalize by S(S+1)
    
    kx = collect(LinRange(pi/2,2.5pi,700))
    # kx = kx
    ky = kx
    kxSpin1 = SqSpin1.kx .+0.5pi# .+pi/size(SqSpin1.Sq,1)
    # kxSpin1 = 2pi .* (0.5:1)
    kySpin1 = kxSpin1

    SqLN = [SqLargeN(kx,ky)/4 for kx in kx , ky in ky]
    SqPlot = [SqFunc(kx,ky) for kx in kx , ky in kx]
    SqSpin1Plot = [Sq1Func(kx,ky) for kx in ky , ky in ky]
    @info "" sum(SqLN)/length(SqLN) sum(SqPlot)/length(SqPlot) sum(SqSpin1.Sq)/length(SqSpin1.Sq) length(kxSpin1)
    cRange = extrema(SqPlot)

    hm = heatmap!(ax,kx,ky,SqLN,colorrange = cRange)
    hm2 = heatmap!(ax2,kx,ky,SqPlot,colorrange = cRange)
    hm3 = heatmap!(ax3,kx,ky,SqSpin1Plot)
    # hm = halfhalfheatmap!(ax,kx,ky,SqFunc,SqLargeN,x->-x+3pi,normalize = true)
    Colorbar(fig[1, 4], hm2,ticks = [0,0.2,0.4,0.6],height = Relative(0.9),label = L"\mathcal{S}^{zz}(q)")
    save("Application/figs/Sq_comparison.png", fig,px_per_unit=3)
    # Colorbar(fig[1, 2], hm)
    fig
    
end