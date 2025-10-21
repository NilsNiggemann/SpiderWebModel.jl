import SpiderWebModel as SW
# using CairoMakie
using MakieHelpers
include("../Stencils/plottingUtils.jl")
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

with_theme(theme_PiTicks()) do 
    # Sq = sqrt.(var(real(SqsGFMC))) ./4
    kx =ky= collect(trueMomenta(-0.5pi,1.5pi,100))
    
    fig = Figure(size = 0.8 .*(350,300))
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect = 1,xticks = PiTicks([0,pi]),yticks = PiTicks([0,pi]),)

    hmFT = heatmap!(ax,kx,ky,SqLargeN,colormap = :viridis)
    
    Colorbar(fig[1,2],hmFT,
    # label = L"\mathcal{S}(\textbf{q})",
    height = Relative(0.95),vertical=true,ticks = SimpleTicks(),flipaxis=true)
    save("../figs/PaperFigs/SqLargeN.pdf",fig)
    fig
end

##
let 
    using CairoMakie  # or GLMakie for GPU rendering
using CairoMakie.Makie.ColorSchemes

# Define a Value Suppressing Colormap (e.g., suppress mid-range values)
function value_suppressing_colormap()
    return cgrad([
        (:blue, 0.0),
        (:white, 0.5),
        (:red, 1.0)
    ]; scale = true)
end

# Generate grid data
nx, ny = 50, 50
x = range(-1, 1, length=nx)
y = range(-1, 1, length=ny)
z = [sin(3π * xi) * cos(3π * yi) for xi in x, yi in y]

# Flatten data for meshscatter
xf = repeat(x, inner=ny)
yf = repeat(y', outer=nx)
zf = vec(z)

# Plot using meshscatter to simulate heatmap
fig = Figure(resolution=(800, 600))
ax = Axis(fig[1, 1]; title="Heatmap using meshscatter with VCS")

# Plot each point colored by z-value
meshscatter!(ax, xf, yf, color=zf,
    colormap=value_suppressing_colormap(),
    markersize=5,
    shading=false
)

Colorbar(fig[1, 2], limits=extrema(zf), colormap=value_suppressing_colormap())

fig

    
end