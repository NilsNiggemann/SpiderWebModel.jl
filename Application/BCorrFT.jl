using Cubature, CairoMakie, StaticArrays,MakieHelpers
function dispClass(kx,ky)
    4 - cos(2*kx) + cos(2*kx - 2*ky) - 2*cos(kx - ky) - cos(2*ky) - 2*cos(kx + ky) + cos(2*kx + 2*ky)
end
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
SqLargeN(q) = SqLargeN(q[1],q[2])

##
using CairoMakie
CairoMakie.activate!()
with_theme(theme_PiTicks()) do
    k = LinRange(-0.5pi,1.5pi,1000)
    fig = Figure(size = 0.9 .*(380, 350), fontsize = 22,backgroundcolor = (:white,0.))
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect = 1;xticks = PiTicks([0,pi]), yticks = PiTicks([0,pi]),title = L"\mathcal{S}(\mathbf{q})")
    
    hm = heatmap!(ax,k,k,SqLargeN)
    Colorbar(fig[1,2],hm,height = Relative(1))
    colgap!(fig.layout,1,10)
    
    save("exactFig/SqLargeN.png",fig,px_per_unit = 3.)
    fig
end

##
with_theme(theme_PiTicks()) do
    k = LinRange(-pi,pi,100)

    BB = fetch.([Threads.@spawn getBBCorrFieldTheory(ki,kj) for ki in k, kj in k])

    heatmap(k,k,BB)
    
end

##
using GLMakie
GLMakie.activate!()
# GLMakie.closeall() # close any open screen

##
with_theme(theme_PiTicks()) do 
# let
    k = LinRange(-pi,pi,200)

    z = BB = fetch.([Threads.@spawn dispClass(ki,kj) for ki in k, kj in k])
    x = y = k

        
    k_WF = -pi:pi/10:pi

    z_WF = BB = fetch.([Threads.@spawn dispClass(ki,kj) for ki in k_WF, kj in k_WF])
    x = y = k

    zmin, zmax = extrema(BB)
    cmap = :viridis
    # zlabel = L"\langle B(\mathbf{q})B(-\mathbf{q})\rangle"
    zlabel = L"$ε_1/J,$  $ε_2/J$"
    # zlabel = ""

    fig = Figure(size = 0.8 .* (640, 550), fontsize = 22)
    # fig = Figure()
    ticklabels = (-pi:pi/2:pi,Makie.latexstring.(["-π","-π/2","0","π/2","π"]))
    ax = Axis3(fig[1, 1]; aspect = (1,1,0.6), perspectiveness = 0.0, elevation = π / 10,azimuth = 1.55π,
        xzpanelcolor = (:black, 0.15), yzpanelcolor = (:black, 0.15),
        zgridcolor = :grey, ygridcolor = :grey, xgridcolor = :grey,xlabel = L"q_x",ylabel = L"q_y",zlabel,zticks = SimpleTicks(),xticks = ticklabels,yticks = ticklabels,protrusions=(30,70,50,-10),ylabeloffset=70)
    # Box(fig[1, 1], strokewidth = 0)
    # surface!(ax, x, y, 0 .*z; colormap = :coolwarm, colorrange = (zmin, zmax),)
    surface!(ax, x, y, 0 .*z; colormap = :viridis,shininess = 30f0,backlight=10f0,specular = 0.8,label = L"\frac{ε_1}{J}")

    # meshpoints = [Point(x,y,) for (x,y) in Iterators.product(k,k)][:]
    # mesh!(ax, meshpoints,color = Pattern('x'), shading = NoShading,label = L"\frac{ε_1}{J}",depth_shift = -0.01)
    sm = surface!(ax, x, y, z; colormap = cmap, colorrange = (zmin, zmax),
    transparency = false,alpha = 1,shininess = 30f0,backlight=10f0,specular = 0.8,label = L"\frac{ε_2}{J}")
    # wireframe!(ax, k_WF, k_WF, 0*z_WF; transparency = true,color = (:black, 0.4),linewidth = 0.5,depth_shift=-0.04)
    
    # xm, ym, zm = minimum(ax.finallimits[])
    # contour!(ax, x, y, z; levels = 20, colormap = cmap, linewidth = 2,
    #     colorrange = (zmin, zmax), transformation = (:xy, zmin),
    #     transparency = true)

    
    wireframe!(ax, k_WF, k_WF, z_WF ; overdraw = false, transparency = true,color = (:black, 1),linewidth = 0.4,depth_shift = -4.0e-2)
    wireframe!(ax, k_WF, k_WF, 0z_WF ; overdraw = false, transparency = true,color = (:black, 1),linewidth = 0.4,depth_shift = -4.0e-2)
    
    # return fig
    # Colorbar(fig[1, 2], sm, height = Relative(0.5))
    zlims!(ax, zmin, zmax + 0.5)
    xlims!(ax, extrema(k)...)
    ylims!(ax, extrema(k)...)
    # return fig
    # Legend(fig[1,2], ax, margin = (-50, 0, 50,    0))
    # colsize!(fig.layout, 1, Aspect(1, 1.0))
    save("exactFig/dispersion_classical.png",fig,px_per_unit = 2)
    # save("exactFig/dispersion_classical.png",alpha_colorbuffer(fig))
    fig
end

##

with_theme(theme_PiTicks()) do 
    # let
    k = LinRange(-pi,pi,30)

    z = BB = fetch.([Threads.@spawn ω(ki,kj) for ki in k, kj in k])
    x = y = k
    zmin, zmax = extrema(BB)
    cmap = :viridis
    # zlabel = L"\langle B(\mathbf{q})B(-\mathbf{q})\rangle"
    zlabel = L"$\langle B^2(q) B^2(-q) \rangle$"

    fig = Figure(size = (600, 450), fontsize = 22)
    # fig = Figure()
    ticklabels = (-pi:pi/2:pi,Makie.latexstring.(["-π","-π/2","0","π/2","π"]))
    ax = Axis3(fig[1, 1]; aspect = (1,1,0.4), perspectiveness = 0.0, elevation = π / 9,azimuth = 1.6π,
        xzpanelcolor = (:black, 0.15), yzpanelcolor = (:black, 0.15),
        zgridcolor = :grey, ygridcolor = :grey, xgridcolor = :grey,xlabel = L"q_x",ylabel = L"q_y",zlabel,zticks = SimpleTicks(),xticks = ticklabels,yticks = ticklabels,protrusions=(30,60,-0,-0),backgroundcolor = :gray97,ylabeloffset=70)
    # Box(fig[1, 1], strokewidth = 0)
    sm = surface!(ax, x, y, z; colormap = cmap, colorrange = (zmin, zmax),
    transparency = true,alpha = 1,shininess = 30f0,backlight=10f0)
    
    # xm, ym, zm = minimum(ax.finallimits[])
    # contour!(ax, x, y, z; levels = 20, colormap = cmap, linewidth = 2,
    #     colorrange = (zmin, zmax), transformation = (:xy, zmin),
    #     transparency = true)
    wireframe!(ax, x[1:1:end], y[1:1:end], z[1:1:end,1:1:end]; overdraw = true, transparency = true,color = (:black, 0.4),linewidth = 0.5)
    
    # return fig
    # Colorbar(fig[1, 2], sm, height = Relative(0.5))
    zlims!(ax, zmin, zmax + 0.5)
    xlims!(ax, extrema(k)...)
    ylims!(ax, extrema(k)...)
    # return fig
    # colsize!(fig.layout, 1, Aspect(1, 1.0))
    # save("exactFig/B_corr_FT.png",current_figure())
    fig
end