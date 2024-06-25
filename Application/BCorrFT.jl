using Cubature, CairoMakie, StaticArrays,MakieHelpers
function dispClass(kx,ky)
    4 - cos(2*kx) + cos(2*kx - 2*ky) - 2*cos(kx - ky) - cos(2*ky) - 2*cos(kx + ky) + cos(2*kx + 2*ky)
end
##
let 
    k = LinRange(-pi,pi,100)

    BB = fetch.([Threads.@spawn getBBCorrFieldTheory(ki,kj) for ki in k, kj in k])

    heatmap(k,k,BB)
    
end

##
# using GLMakie
# GLMakie.activate!()
# GLMakie.closeall() # close any open screen
##
with_theme(theme_PiTicks()) do 
# let
    k = LinRange(-pi,pi,200)

    z = BB = fetch.([Threads.@spawn dispClass(ki,kj) for ki in k, kj in k])
    x = y = k
    zmin, zmax = extrema(BB)
    cmap = :viridis
    # zlabel = L"\langle B(\mathbf{q})B(-\mathbf{q})\rangle"
    zlabel = L"$ε_1/J,$  $ε_2/J$"

    fig = Figure(size = (600, 450), fontsize = 22)
    # fig = Figure()
    ticklabels = (-pi:pi/2:pi,Makie.latexstring.(["-π","-π/2","0","π/2","π"]))
    ax = Axis3(fig[1, 1]; aspect = (1,1,0.4), perspectiveness = 0.0, elevation = π / 9,azimuth = 1.6π,
        xzpanelcolor = (:black, 0.15), yzpanelcolor = (:black, 0.15),
        zgridcolor = :grey, ygridcolor = :grey, xgridcolor = :grey,xlabel = L"q_x",ylabel = L"q_y",zlabel,zticks = SimpleTicks(),xticks = ticklabels,yticks = ticklabels,protrusions=(30,60,-0,-0),backgroundcolor = :gray97,ylabeloffset=70)
    # Box(fig[1, 1], strokewidth = 0)
    surface!(ax, x, y, 0 .*z; colormap = cmap, colorrange = (zmin, zmax),)
    sm = surface!(ax, x, y, z; colormap = cmap, colorrange = (zmin, zmax),
    transparency = true,alpha = 1,shininess = 30f0,backlight=10f0)
    
    # xm, ym, zm = minimum(ax.finallimits[])
    # contour!(ax, x, y, z; levels = 20, colormap = cmap, linewidth = 2,
    #     colorrange = (zmin, zmax), transformation = (:xy, zmin),
    #     transparency = true)
    wireframe!(ax, x[1:6:end], y[1:6:end], z[1:6:end,1:6:end]; overdraw = true, transparency = true,color = (:black, 0.4),linewidth = 0.5)
    
    # return fig
    # Colorbar(fig[1, 2], sm, height = Relative(0.5))
    zlims!(ax, zmin, zmax + 0.5)
    xlims!(ax, extrema(k)...)
    ylims!(ax, extrema(k)...)
    # return fig
    # colsize!(fig.layout, 1, Aspect(1, 1.0))
    save("exactFig/dispersion_classical.png",fig,px_per_unit = 2)
    fig
end

##
function ω(kx,ky)
    sx,cx = sincos(kx)
    sy,cy = sincos(ky)
    w2 = (cx - cy)^2 + 4*(sx*sy)^2
    return sqrt(w2)
end
ω(k) = ω(k[1],k[2])
function getBBCorrFieldTheory(kx,ky)
    f(k) = ω(k)*ω(SA[kx,ky]-k)
    return hcubature(f, SA[-pi/2,-pi], SA[pi/2,pi], reltol = 1e-6, abstol = 1e-6)[1]
end

with_theme(theme_PiTicks()) do 
    # let
    k = LinRange(-pi,pi,200)

    z = BB = fetch.([Threads.@spawn getBBCorrFieldTheory(ki,kj) for ki in k, kj in k])
        x = y = k
        zmin, zmax = extrema(BB)
        cmap = :viridis
        # zlabel = L"\langle B(\mathbf{q})B(-\mathbf{q})\rangle"
        zlabel = L"$ε_1/J,$  $ε_2/J$"
    
        fig = Figure(size = (600, 450), fontsize = 22)
        # fig = Figure()
        ticklabels = (-pi:pi/2:pi,Makie.latexstring.(["-π","-π/2","0","π/2","π"]))
        ax = Axis3(fig[1, 1]; aspect = (1,1,0.4), perspectiveness = 0.0, elevation = π / 9,azimuth = 1.6π,
            xzpanelcolor = (:black, 0.15), yzpanelcolor = (:black, 0.15),
            zgridcolor = :grey, ygridcolor = :grey, xgridcolor = :grey,xlabel = L"q_x",ylabel = L"q_y",zlabel,zticks = SimpleTicks(),xticks = ticklabels,yticks = ticklabels,protrusions=(30,60,-0,-0),backgroundcolor = :gray97,ylabeloffset=70)
        # Box(fig[1, 1], strokewidth = 0)
        surface!(ax, x, y, 0 .*z; colormap = cmap, colorrange = (zmin, zmax),)
        sm = surface!(ax, x, y, z; colormap = cmap, colorrange = (zmin, zmax),
        transparency = true,alpha = 1,shininess = 30f0,backlight=10f0)
        
        # xm, ym, zm = minimum(ax.finallimits[])
        # contour!(ax, x, y, z; levels = 20, colormap = cmap, linewidth = 2,
        #     colorrange = (zmin, zmax), transformation = (:xy, zmin),
        #     transparency = true)
        wireframe!(ax, x[1:6:end], y[1:6:end], z[1:6:end,1:6:end]; overdraw = true, transparency = true,color = (:black, 0.4),linewidth = 0.5)
        
        # return fig
        # Colorbar(fig[1, 2], sm, height = Relative(0.5))
        zlims!(ax, zmin, zmax + 0.5)
        xlims!(ax, extrema(k)...)
        ylims!(ax, extrema(k)...)
        # return fig
        # colsize!(fig.layout, 1, Aspect(1, 1.0))
        # save("exactFig/dispersion_classical.png",fig,px_per_unit = 2)
        fig
    end
end