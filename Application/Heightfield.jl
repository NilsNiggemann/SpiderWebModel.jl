import SpiderWebModel as SW
using CairoMakie
using MakieHelpers
##
θ(x) = x>=0 ? 1 : 0

function spin_HF(x,y,heightfield)
    if iseven(x+y)
        return heightfield(x+1,y) - heightfield(x,y+1) + heightfield(x-1,y) - heightfield(x,y-1)
    else
        return -heightfield(x+1,y+1) + heightfield(x-1,y+1) - heightfield(x-1,y-1) + heightfield(x+1,y-1)
    end
end

function linFunc(x) 

    a = -0.5
    b =0.5
         
    if x <b
        f = max(0.5x-a,0)
        f = min(f,1)
    else
        f = min(b-0.5x+1,1)
        f = max(f,0)
    end
    return f
end

# Fluc_Heightfield(x,y) = linFunc(x) * linFunc(y)

function Fluc_Heightfield(x,y)
    iseven(x+y) && return NaN
    x==y-1==0 && return 1
    return 0
end
# function Fluc_Heightfield(x,y)
#     base=3
#     top=1
#     height=1

#     # Slope of the pyramid sides
#     slope = (base - top) / (2 * height)
    
#     # if x==y==0
#     #     return height
#     # end
#     # Compute z based on x and y
#     if abs(x) <= base / 2 && abs(y) <= base / 2  # Inside base square
#         if abs(x) <= top / 2 && abs(y) <= top / 2  # Inside top square
#             return height  # Top surface
#         else  # Slanted sides
#             return height - slope * (max(abs(x) - top / 2, abs(y) - top / 2))
#         end
#     else
#         return 0  # Outside the base
#     end
# end

##
with_theme(theme_SimpleTicks()) do
    fig = Figure(size=(900, 700))
    SpinConstruction = GridLayout()
    Fluctuator = GridLayout()
    LineMove = GridLayout()
    Configs = GridLayout()
    fig.layout[1,1:4] = SpinConstruction
    fig.layout[2:3,1:2] = Fluctuator
    fig.layout[2:3,3:4] = LineMove
    fig.layout[4,1:4] = Configs

    SConf = SW.stencilConfig(zeros(3,3),1)

    # Panel (a): Schematic-like gridw
    SpinConstruction[1,1:2] = ax1 = Axis(fig,
    #  title=L"(a)$$";
    SW.getConfigAxis(SConf)...,
    aspect=1
    )
    # SpinConstruction[1,3:4] = axSpin2 = Axis(fig, title="(b)";SW.getConfigAxis(SConf)...)
    SpinConstruction[1,3:4] = axSpin2 = Axis(fig, 
    # title="(b)",
    aspect=1)
    text!(ax1, [1, 2, 3, 1, 3, 1.5, 2.5], [1, 1, 1, 3, 3, 2, 2], text=["h", "h", "h", "h", "h", "S", "h"], fontsize=20)
    # scatter!(ax1, [2, 2.5, 1.5, 2], [2.5, 2, 2, 1.5], color=:red, markersize=10)
    # Panel (b): 3D Surface Plots
    Fluctuator[1,1] = ax2 = Axis3(fig,aspect = (1,1,1))
    # Fluctuator[1,1] = axFluc = Axis(fig, title="(b)")
    
    x_discrete = -3:3
    y_discrete = -2:3

    x = LinRange(extrema(x_discrete)...,100)
    y = LinRange(extrema(y_discrete)...,100)

    x = x_discrete
    y = y_discrete
    # y = LinRange(extrema(y_discrete)...,100)

    z = [Fluc_Heightfield(x,y) for x in x, y in y]
    z_discrete = [Fluc_Heightfield(x,y) for x in x_discrete, y in y_discrete]

    # wireframe!(ax2, x_discrete, y_discrete, 0*z_discrete, color = :black)
    # surface!(ax2, x, y, z, colormap=:viridis,interpolate=false)

    box = Rect3(Point3f(-0.5), Vec3f(1))
    boxmarker(z) = Rect3(Point3f(-0.5*z), Vec3f(z))
    
    meshscatter!(ax2, x, y, z,
    # marker=vec(boxmarker.(z)),
    marker=box,
    markersize = 0.95,color = vec(z),colormap=:viridis)
    # x_discrete = x
    # y_discrete = y
    heatmap!(ax2,x_discrete,y_discrete, [spin_HF(x,y,Fluc_Heightfield) for x in x_discrete, y in y_discrete], transformation = (:xy, minimum(x_discrete)),colormap=:viridis)

    let 
        xgridlines = (minimum(x_discrete)-1:maximum(x_discrete)) .+0.5
        ygridlines = (minimum(y_discrete)-1:maximum(y_discrete)) .+0.5

        wireframe!(ax2, xgridlines,ygridlines, [minimum(x_discrete) for x in xgridlines, y in ygridlines], color = :black)
        # wireframe!(ax2, xgridlines,ygridlines, [0 for x in xgridlines, y in ygridlines], color = :black)
    xgridlines,ygridlines
    end

    # heatmap!(axSpin2,x_discrete,y_discrete, [Fluc_Heightfield(x,y) for x in x_discrete, y in y_discrete], transformation = (:xy, 0),colormap=:viridis)
    heatmap!(ax1,x_discrete,y_discrete, [Fluc_Heightfield(x,y) for x in x_discrete, y in y_discrete], transformation = (:xy, 0),colormap=:viridis)
    heatmap!(axSpin2,x_discrete,y_discrete, [spin_HF(x,y,Fluc_Heightfield) for x in x_discrete, y in y_discrete], transformation = (:xy, 0),colormap=:viridis)

    LineMove[1,1] = ax3 = Axis3(fig)
    z2 = [linFunc(y) for x in x, y in y]
    surface!(ax3, x, y, z2, colormap=:viridis,interpolate=false)
    

    function transformSpins!(vec,sgn)
        any(==(sgn*spin),vec) && return
        vec .+= sgn
        return
    end

    M⎸ = SW.stencilConfig(zeros(8,8),1)
    SW.transformSpinsAlongCol!(M⎸,3,x->x-2)
    
    M_ = SW.stencilConfig(zeros(8,8),1)
    SW.transformSpinsAlongRow!(M_,3,x->x-2)

    M╱ = SW.stencilConfig(zeros(8,8),1)

    # diag = SW.getDiagonal(S,3,1,true)
    # transformSpins!(diag,1)


    SW.transformSpinsAlongDiagonal!(M╱,3,1,x->x-2)

    M╲ = SW.stencilConfig(zeros(8,8),1)
    SW.transformSpinsAlongDiagonal!(M╲,5,-1,x->x-2)

    Confs = [M⎸,M_,M╱,M╲]
    # labels = ["M⎸","M_","M╱","M╲"]

    # $M_{\bm \vert}$, $M_{\bm -}$, $M_{\bm \diagup}$, and $M_{\bm \diagdown}$ 
    labels = [
        L"M_{\mathbf{|}}", 
        L"M_{\mathbf{-}}",
        L"M_{\mathbf{-}}",
        L"M_{\mathbf{\diagup}}", 
        L"M_{\mathbf{\diagdown}}"
        ]
    
    zlims!(ax2,extrema(x_discrete)...)

    for i in 1:4
        Configs[1, i] = ax = Axis(fig, title=labels[i];SW.getConfigAxis(Confs[i])...)
        # SW.plotSpinConfig!(ax, SLine)
        SW.plotFractons!(ax, Confs[i])
    end

    Label(SpinConstruction[1,1,TopLeft()],L"(a)$$")
    Label(SpinConstruction[1,3,TopLeft()],L"(b)$$")
    Label(Fluctuator[1,1,TopLeft()],L"(c)$$")
    Label(LineMove[1,1,TopLeft()],L"(d)$$")
    Label(Configs[1,1,TopLeft()],L"(e)$$")

    fig
end
