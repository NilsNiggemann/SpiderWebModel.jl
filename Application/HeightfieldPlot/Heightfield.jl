import SpiderWebModel as SW
# using CairoMakie
using MakieHelpers
using WGLMakie

function filledscatter!(ax,args...;marker,markersize = 30,color,circlekwargs=(;),kwargs...)
    scatter!(ax,args...;
    marker='⬤',
    markersize,
    color=:white,
    strokewidth=2,
    strokecolor=color,
    circlekwargs...
)

    scatter!(ax,args...;
        marker,
        markersize,
        color,
        kwargs...
    )

end

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

function Line_Heightfield(x,y)
    x,y = y,x
    iseven(x+y) && return NaN
    # return linFunc(x+1)
    if x==0 || abs(x)<=2
        # return 2-abs(x)
        return 1
    end
    return 0
end

##
with_theme(theme_SimpleTicks()) do
    fig = Figure(size= 0.8 .*(800, 900),fontsize = 19)
    SpinConstruction = GridLayout()
    Fluctuator = GridLayout()
    LineMove = GridLayout()
    Configs = GridLayout()
    fig.layout[1,1:4] = SpinConstruction
    fig.layout[2:3,1:2] = Fluctuator
    fig.layout[2:3,3:4] = LineMove
    fig.layout[4,1:4] = Configs

    SConf1 = SW.stencilConfig(zeros(3,3),1)
    SConf2 = SW.stencilConfig(zeros(4,4),1)
    # Panel (a): Schematic-like grid
    SpinConstruction[1,1:2] = axSpin1 = Axis(fig;
    #  title=L"(a)";
    xticklabelsvisible=false,
    yticklabelsvisible=false,
    SW.getConfigAxis(SConf1)...,
    aspect=1
    )


    # SpinConstruction[1,3:4] = axSpin2 = Axis(fig, title="(b)";SW.getConfigAxis(SConf)...)
    SpinConstruction[1,3:4] = axSpin2 = Axis(fig;
    # title="(b)",
    xticklabelsvisible=false,
    yticklabelsvisible=false,
    SW.getConfigAxis(SConf2)...,
    aspect=1)

    let # (a)
        
        SW.plotSpinConfig!(axSpin1, SConf1)
        
        scatter!(axSpin1, Point2[(2.5,2),(1.5,2)],marker = ['←', '→'],color = :black,markersize = 24)
        # scatter!(axSpin1, Point2[(2.75,2),(1.25,2)],marker = '⊕',color = :red,markersize = 30)
        filledscatter!(axSpin1, Point2[(2.75,2),(1.25,2)],marker = '+',color = :red,markersize = 24)
        
        
        scatter!(axSpin1, Point2[(2,1.5),(2,2.5)],marker = ['↑', '↓'],color = :black,markersize = 24)
        filledscatter!(axSpin1, Point2[(2,1.25),(2,2.75)],marker = '-',color = :blue,markersize = 24)

        text!(axSpin1, Point2[(2,2),(2,0.75),(0.7,2),(3.25,2),(2,3.25)], text=[L"S^z_⊠", L"h",L"h",L"h",L"h"], fontsize=24,align = (:center,:center))
    end
    
    let # (a)
        
        SW.plotSpinConfig!(axSpin2, SConf2)
    
        filledscatter!(axSpin2, Point2[(2.75, 3.75), (1.25, 2.25)], marker = '-', color = :blue, markersize = 24)
        
        filledscatter!(axSpin2, Point2[(1.25, 3.75), (2.75, 2.25)], marker = '+', color = :red, markersize = 24)
        
        scatter!(axSpin2, Point2[(2.5,2.5),(1.5,2.5),(1.5, 3.5), (2.5, 3.5)], marker = ['↖','↗', '↘','↙'], color = :black, markersize = 24)
        text!(axSpin2, Point2[(2, 3), (1, 4), (3, 4), (1, 2), (3, 2)], text = [L"S^z_□", L"h", L"h", L"h", L"h"], fontsize = 24, align = (:center, :center))
        ylims!(axSpin2,1.5,4.5)
        xlims!(axSpin2,0.5,3.5)
    end

    # Panel (b): 3D Surface Plots
    # Fluctuator[1,1] = axFluc = Axis(fig, title="(b)")
    min_height = 0.05
    
    x_discrete = -2:2
    y_discrete = -1:3

    z_bottom_Fluc = -2
    z_bottom_Line = -4

    x = LinRange(extrema(x_discrete)...,100)
    y = LinRange(extrema(y_discrete)...,100)

    x = x_discrete
    y = y_discrete

    let # (b)

        xgridlines = (minimum(x_discrete)-1:maximum(x_discrete)) .+0.5
        ygridlines = (minimum(y_discrete)-1:maximum(y_discrete)) .+0.5

        Fluctuator[1,1] = axFluc = Axis3(fig,aspect=:data,zlabel = L"h",xlabel = L"x",ylabel = L"y",xticklabelsvisible=false,yticklabelsvisible=false,xticks = xgridlines,yticks = ygridlines,zticks = [0,1],
        xlabeloffset = 5,
        ylabeloffset = 5,
        zlabeloffset = 20,
        )


        z = [Fluc_Heightfield(x,y) for x in x, y in y]
        z_discrete = [Fluc_Heightfield(x,y) for x in x_discrete, y in y_discrete]
    
        meshscatter!(axFluc, x, y, z,
        marker = Makie.Rect3D(Vec3f(-0.5, -0.5, -1.0 -min_height), Vec3f(1)),
        markersize = Vec3f.(1.0, 1.0, min_height .+vec(z)),
        # markersize = 0.95,
        color = vec(z),colormap=:viridis)

        heatmap!(axFluc,x_discrete,y_discrete, [spin_HF(x,y,Fluc_Heightfield) for x in x_discrete, y in y_discrete], transformation = (:xy, z_bottom_Fluc),colormap=:greys,colorrange = [-1,1])
    
        # x_markers = scatter!(axFluc,[Point(x,y,z_bottom) for x in x_discrete, y in y_discrete if iseven(x+y)],marker = '×',color = :lightgray,markersize = 0.5,markerspace = :data)

        # translate!(x_markers,0,0,0.1)

        let 

            xpoints = [Point(x,y,z_bottom_Fluc+0.01) for x in x_discrete, y in y_discrete if iseven(x+y)]
            text!(axFluc,xpoints;
            text = fill("×", length(xpoints)),
            color = :lightgray,
            # rotation = [i / 7 * 1.5pi for i in 1:7],
            align = (:center, :center),
            fontsize = 0.3,
            markerspace = :data
            )
            
        end

        
        # meshscatter!(axFluc,[Point(x,y,z_bottom) for x in x_discrete, y in y_discrete if iseven(x+y)],marker = '×',color = :lightgray,markersize = 50)

        wireframe!(axFluc, xgridlines,ygridlines, [z_bottom_Fluc for x in xgridlines, y in ygridlines], color = :black)
        rmnan(x) = isnan(x) ? 0 : x
        wireframe!(axFluc, xgridlines,ygridlines, [0 for x in xgridlines, y in ygridlines], color = :black)

    end

    let # (d)
    
        x_discrete = -5:5
        y_discrete = -4:4
                
        xgridlines = (minimum(x_discrete)-1:maximum(x_discrete)) .+0.5
        ygridlines = (minimum(y_discrete)-1:maximum(y_discrete)) .+0.5

        LineMove[1,1] = axLine = Axis3(fig,zlabel = L"h",xlabel = L"x",ylabel = L"y",xticklabelsvisible=false,yticklabelsvisible=false,aspect=:data,xticks = xgridlines,yticks = ygridlines,zticks = [0,1],
        xlabeloffset = 5,
        ylabeloffset = 5,
        
        )


        z2 = [Line_Heightfield(x,y) for x in x_discrete, y in y_discrete]
        # surface!(ax3, x, y, z2, colormap=:viridis,interpolate=false)
        meshscatter!(axLine, x_discrete, y_discrete, z2,
        marker = Makie.Rect3D(Vec3f(-0.5, -0.5, -1.0 -min_height), Vec3f(1)),
        markersize = Vec3f.(1.0, 1.0, min_height .+vec(z2)),
        color = vec(z2),colormap=:viridis)
    
        heatmap!(axLine,x_discrete,y_discrete, [spin_HF(x,y,Line_Heightfield) for x in x_discrete, y in y_discrete], transformation = (:xy, z_bottom_Line),colormap=:greys,colorrange = [-1,1])

        wireframe!(axLine, xgridlines,ygridlines, [z_bottom_Line for x in xgridlines, y in ygridlines], color = :black)
        rmnan(x) = isnan(x) ? 0 : x
        wireframe!(axLine, xgridlines,ygridlines, [0 for x in xgridlines, y in ygridlines], color = :black)


        let 

            xpoints = [Point(x,y,z_bottom_Line+0.01) for x in x_discrete, y in y_discrete if iseven(x+y)]
            text!(axLine,xpoints;
            text = fill("×", length(xpoints)),
            color = :lightgray,
            # rotation = [i / 7 * 1.5pi for i in 1:7],
            align = (:center, :center),
            fontsize = 0.2,
            markerspace = :data
            )
            
        end
    end

    function transformSpins!(vec,sgn)
        any(==(sgn),vec) && return
        vec .+= sgn
        return
    end

    M⎸ = SW.stencilConfig(zeros(8,8),1)
    SW.transformSpinsAlongCol!(M⎸,3,x->x-2)
    
    M_ = SW.stencilConfig(zeros(8,8),1)
    SW.transformSpinsAlongRow!(M_,3,x->x-2)

    M╱ = SW.stencilConfig(zeros(8,8),1)

    diag = SW.getDiagonal(M╱,3,1,true)
    transformSpins!(diag,-1)
    # SW.transformSpinsAlongDiagonal!(M╱,3,1,x->x-2)

    M╲ = SW.stencilConfig(zeros(8,8),1)
    diag = SW.getDiagonal(M╲,3,-1,true)
    transformSpins!(diag,-1)

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
    
    # zlims!(axLine,extrema(x_discrete)...)

    for i in 1:4
        Configs[1, i] = ax = Axis(fig, title=labels[i];
        SW.getConfigAxis(Confs[i])...,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        )
        # SW.plotSpinConfig!(ax, SLine)
        SW.plotFractons!(ax, Confs[i])
        # heatmap!(ax, Confs[i],colormap = :greys,colorrange = [-1,1],)
    end
    rowsize!(fig.layout,1,Relative(0.25))
    rowgap!(fig.layout,2,-100)
    rowgap!(fig.layout,3,-50)
    Label(SpinConstruction[1,1,TopLeft()],L"(a)")
    # Label(SpinConstruction[1,3,TopLeft()],L"(b)")
    Label(Fluctuator[1,1,TopLeft()],L"(b)")
    Label(LineMove[1,1,TopLeft()],L"(c)")
    Label(Configs[1,1,TopLeft()],L"(d)")
    save("../figs/PaperFigs/Heightfieldplot.png",fig,px_per_unit=3.)
    
    fig
end
