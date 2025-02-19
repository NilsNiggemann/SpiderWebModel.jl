import SpiderWebModel as SW
# using CairoMakie
using MakieHelpers
using WGLMakie
using ForwardDiff
using SpiderWebModel.StaticArrays
using LinearAlgebra
using Interpolations
using ImageFiltering
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

function spin_HF_continuum(x,y,heightfield)
    Hess = ForwardDiff.hessian(heightfield,SA[x,y])

    Eˣʸ = 2*(Hess[1,1] - Hess[2,2]) # 2*(∂²h/∂x² - ∂²h/∂y²)
    Eˣˣ = -4Hess[1,2] # -4∂²h/∂x∂y

    E = SA[
        Eˣˣ Eˣʸ;
        Eˣʸ -Eˣˣ
    ]
    return norm(E )
    # return abs2(Eˣˣ) + abs2(Eˣʸ)
    # return (Eˣˣ)

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

function pyramid(x, y; base=2.5, height=1.0)
    slope = height / (base / 2)
    x = x - 1
    if abs(x) <= base / 2 && abs(y) <= base / 2  # Inside base square
        return max(0, height - slope * max(abs(x), abs(y)))
    else
        return 0.0  # Outside the base
    end
end

function Fluc_Continuum_sharp(x, y; base=2.5, top=0.8, height=1)
    # Slope of the pyramid sides
    # Compute z based on x and y
    return min(top, pyramid(x, y; base=base, height=height))
end

# function Fluc_Continuum_1(x, y; smoothing=0.1,kwargs...)
#     neigh(phi) = SA[x,y] + smoothing*SA[cos(phi),sin(phi)]
#     return sum(Fluc_Continuum_sharp(neigh(phi)...;kwargs...) for phi in LinRange(0,2π,100))
# end
Fluc_Continuum_1 = let 
    
    x = -5:0.01:5
    y = -5:0.01:5
    z = ImageFiltering.imfilter([2.5Fluc_Continuum_sharp(x,y,top = 0.6) for x in x, y in y],Kernel.gaussian(3))

    itp = Interpolations.cubic_spline_interpolation((x,y),z)
    
    Fluc_Continuum_func(x,y) = itp(x,y)
    Fluc_Continuum_func(r::AbstractVector) = itp(r[1], r[2])
end


Line_Continuum_1 = let 
    
    x = -5:0.01:5
    y = -5:0.01:5
    z = ImageFiltering.imfilter([8Fluc_Continuum_sharp(0,y,top = 0.75) for x in x, y in y],Kernel.gaussian(5))

    itp = Interpolations.cubic_spline_interpolation((x,y),z)
    
    Line_Continuum_func(x,y) = itp(x,y)
    Line_Continuum_func(r::AbstractVector) = itp(r[1], r[2])
end

function manhattan_distance(x, y)
    return abs(x) + abs(y)
end

function Fluc_Heightfield_discrete(x,y;height=3)
    # Rotate (x, y) by 45 degrees
    iseven(x+y) && return NaN
    # dist = SW.norm(SW.SA[x,y] -SW.SA[1,0])
    dist = manhattan_distance(x,y-1)
    return floor(Int,max(0,height-dist/2))
end

function Line_Heightfield(x,y)
    x,y = y,x
    iseven(x+y) && return NaN
    # return linFunc(x+1)
    if abs(x)<=2
        # return 2-abs(x)
        return 1
    end
    return 0
end

##
with_theme(theme_SimpleTicks()) do
    fig = Figure(size= 0.8 .*(800, 1000),fontsize = 19)
    SpinConstruction = GridLayout()
    Fluctuator = GridLayout()
    LineMove = GridLayout()
    Configs = GridLayout()
    fig.layout[1,1:4] = SpinConstruction
    fig.layout[2:4,1:2] = Fluctuator
    fig.layout[2:4,3:4] = LineMove
    fig.layout[5,1:4] = Configs

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

    z_bottom_Fluc_cont = -2.5
    z_bottom_Fluc = -2
    z_bottom_bigFluc = -5
    z_bottom_line = -3


    let # (b)

        xgridlines = (minimum(x_discrete)-1:maximum(x_discrete)) .+0.5
        ygridlines = (minimum(y_discrete)-1:maximum(y_discrete)) .+0.5

        Fluctuator[1,1] = axFluc = Axis3(fig,aspect=:data,zlabel = L"h",xlabel = L"x",ylabel = L"y",xticklabelsvisible=false,yticklabelsvisible=false,xticks = xgridlines,yticks = ygridlines,zticks = [0,1],
        xlabeloffset = 5,
        ylabeloffset = 5,
        viewmode = :fitzoom,
        zlabeloffset = 20,
        protrusions=(0,0,0,0),

        )


        z_discrete = [Fluc_Heightfield(x,y) for x in x_discrete, y in y_discrete]
    
        meshscatter!(axFluc, x_discrete, y_discrete, z_discrete,
        marker = Makie.Rect3D(Vec3f(-0.5, -0.5, -1.0 -min_height), Vec3f(1)),
        markersize = Vec3f.(1.0, 1.0, min_height .+vec(z_discrete)),
        # markersize = 0.95,
        color = vec(z_discrete),colormap=:viridis)

        heatmap!(axFluc,x_discrete,y_discrete, [spin_HF(x,y,Fluc_Heightfield) for x in x_discrete, y in y_discrete], transformation = (:xy, z_bottom_Fluc),colormap=:greys,colorrange = [-1,1])
    
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
        # (x) = isnan(x) ? 0 : x
        wireframe!(axFluc, xgridlines,ygridlines, [0 for x in xgridlines, y in ygridlines], color = :black)

    end

    let # (c)
    
        x_discrete = -5:5
        y_discrete = -4:6
                
        xgridlines = (minimum(x_discrete)-1:maximum(x_discrete)) .+0.5
        ygridlines = (minimum(y_discrete)-1:maximum(y_discrete)) .+0.5

        LineMove[1,1] = axLine = Axis3(fig,zlabel = L"h",xlabel = L"x",ylabel = L"y",xticklabelsvisible=false,yticklabelsvisible=false,aspect=:data,xticks = xgridlines,yticks = ygridlines,zticks = [0,1,2,3],
        xlabeloffset = 5,
        ylabeloffset = 5,
        zlabeloffset = 25,
        viewmode = :fitzoom,
        protrusions=(0,0,0,0),
        )

        BigFluc(x,y) = Fluc_Heightfield_discrete(x,y,height = 3)

        z2 = [BigFluc(x,y) for x in x_discrete, y in y_discrete]
        # surface!(ax3, x, y, z2, colormap=:viridis,interpolate=false)
        meshscatter!(axLine, x_discrete, y_discrete, z2,
        marker = Makie.Rect3D(Vec3f(-0.5, -0.5, -1.0 -min_height), Vec3f(1)),
        markersize = Vec3f.(1.0, 1.0, min_height .+vec(z2)),
        color = vec(z2),colormap=:viridis)
    
        heatmap!(axLine,x_discrete,y_discrete, [spin_HF(x,y,BigFluc) for x in x_discrete, y in y_discrete], transformation = (:xy, z_bottom_bigFluc),colormap=:greys,colorrange = [-1,1])

        wireframe!(axLine, xgridlines,ygridlines, [z_bottom_bigFluc for x in xgridlines, y in ygridlines], color = :black)
        wireframe!(axLine, xgridlines,ygridlines, [0 for x in xgridlines, y in ygridlines], color = :black)


        let 

            xpoints = [Point(x,y,z_bottom_bigFluc+0.01) for x in x_discrete, y in y_discrete if iseven(x+y)]
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
    let # (d)
        
        x = LinRange(-1,4,300) .-0.5
        y = LinRange(-2,3,300) .-0.5

        Fluctuator[2,1] = axFluc_cont = Axis3(fig,aspect = :data,zlabel = L"h",xlabel = L"x",ylabel = L"y",xticklabelsvisible=false,yticklabelsvisible=false,
        xlabeloffset = 5,
        ylabeloffset = 5,
        zlabeloffset = 20,
        zticklabelsvisible = false,
        zticksvisible = false,
        viewmode = :fitzoom,
        protrusions=(0,40,0,0),
        xlabelcolor=:white,
        xlabelvisible=false,
        ylabelcolor=:white,
        ylabelvisible=false,
        )


        z = [Fluc_Continuum_1(x,y) for x in x, y in y]
    
        surface!(axFluc_cont, x, y, z,
        color = vec(z),colormap=:viridis,alpha = 1,shading=FastShading,shininess=10000f0,specular=1.0,fxaa=true,ssao=true)
        # return [spin_HF_continuum(x,y,Fluc_Continuum_1) for x in x, y in y]
        heatmap!(axFluc_cont,x,y, [spin_HF_continuum(x,y,Fluc_Continuum_1) for x in x, y in y], transformation = (:xy, z_bottom_Fluc_cont),colormap=:hot)
    
    end

    let # (e)
        
        x = LinRange(-3,3,300)
        y = LinRange(-3,3,300)

        LineMove[2,1] = axLineCont = Axis3(fig,aspect = :data,zlabel = L"h",xlabel = L"x",ylabel = L"y",xticklabelsvisible=false,yticklabelsvisible=false,
        xlabeloffset = 5,
        ylabeloffset = 5,
        zlabeloffset = 20,
        viewmode = :fitzoom,
        zticklabelsvisible = false,
        zticksvisible = false,
        protrusions=(0,-10,-10,0),
        xlabelcolor=:white,
        xlabelvisible=false,
        ylabelcolor=:white,
        ylabelvisible=false,

        )


        z = [Line_Continuum_1(x,y) for x in x, y in y]
    
        surface!(axLineCont, x, y, z,
        color = vec(z),colormap=:viridis,alpha = 1,shading=FastShading,shininess=10000f0,specular=1.0,fxaa=true,ssao=true)
        heatmap!(axLineCont,x,y, [spin_HF_continuum(x,y,Line_Continuum_1) for x in x, y in y], transformation = (:xy, z_bottom_line),colormap=:hot)
    
    end

    function transformSpins!(vec,sgn)
        any(==(sgn),vec) && return
        vec .+= sgn
        return
    end

    M = SW.stencilConfig(zeros(8,8),1)
    col = CartesianIndices(M)[3,1:2:end]

    M_ = SW.stencilConfig(zeros(8,8),1)
    row = CartesianIndices(M)[1:2:end,3]

    diag = SW.getDiagonal(M,3,1,true)

    antidiag = SW.getDiagonal(M,3,-1,true)

    changesites = [col,row,CartesianIndices(M)[diag.indices[1]],CartesianIndices(M)[antidiag.indices[1]]]
    # labels = ["M⎸","M_","M╱","M╲"]
    # $M_{\bm \vert}$, $M_{\bm -}$, $M_{\bm \diagup}$, and $M_{\bm \diagdown}$ 
    labels = [
        L"M_{\mathbf{|}}", 
        L"M_{\mathbf{-}}",
        L"M_{\mathbf{\diagup}}", 
        L"M_{\mathbf{\diagdown}}"
        ]
    
    # zlims!(axLine,extrema(x_discrete)...)
    let 
        bound_points = 0.5*[
            Point(-1,-1),
            Point(1,-1),
            Point(1,1),
            Point(-1,1),
            Point(-1,-1),
        ]
        function outline_site!(ax,site,args...;kwargs...)
            Ps = [Point(Tuple(site))-bp for bp in bound_points]
            lines!(Ps,args...;kwargs...)
        end

        for i in 1:4
            Configs[1, i] = ax = Axis(fig, title=labels[i];
            SW.getConfigAxis(M)...,
            xticklabelsvisible=false,
            yticklabelsvisible=false,
            )
            # SW.plotSpinConfig!(ax, SLine)
            SW.plotFractons!(ax, M)
            for s in changesites[i]
                outline_site!(ax,s;color=:red,linewidth = 2)
                SW.plotSiteHighlight!(ax,Tuple(s);color = (:red,0.3))
            end
            # heatmap!(ax, Confs[i],colormap = :greys,colorrange = [-1,1],)
        end
    end
 
    # rowsize!(SpinConstruction,1,Relative(1))
    rowsize!(fig.layout,1,Relative(0.18))
    # rowsize!(Fluctuator,2,Relative(0.4))

    rowsize!(Fluctuator,1,Relative(0.5))
    rowsize!(Fluctuator,2,Relative(0.5))

    rowsize!(LineMove,1,Relative(0.5))
    rowsize!(LineMove,2,Relative(0.5))

    rowgap!(fig.layout,1,-10)
    rowgap!(LineMove,1,-100)
    rowgap!(Fluctuator,1,-100)
    rowgap!(fig.layout,4,-10)
    # rowgap!(fig.layout,3,-50)
    Label(SpinConstruction[1,1,TopLeft()],L"(a)")
    # Label(SpinConstruction[1,3,TopLeft()],L"(b)")
    Label(Fluctuator[1,1,TopLeft()],L"(b)",padding = (0,0,-80,0))
    Label(LineMove[1,1,TopLeft()],L"(c)",padding = (0,0,-80,0))
    Label(Fluctuator[2,1,TopLeft()],L"(d)",padding = (0,0,-120,0))
    Label(LineMove[2,1,TopLeft()],L"(e)",padding = (0,0,-120,0))
    Label(Configs[1,1,TopLeft()],L"(f)")
    save("../figs/PaperFigs/Heightfieldplot.png",fig,px_per_unit=3.)
    
    fig
end
