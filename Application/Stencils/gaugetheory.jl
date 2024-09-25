using CairoMakie, MakieHelpers
import SpiderWebModel as SW
include("plottingUtils.jl")
##
a = []
push!(a,1)
push!(a,"a string")

##
b = a +2
##
with_theme(theme_PiTicks()) do 
    # k = LinRange(-0.5pi, 1.50pi, 32)
    k = trueMomenta(-0.5pi, 1.50pi, 400)

    # Define three different values for r
    r_values = [1e2,1,1e-2]
    # Create a 3x1 figure with shared y-axis
    
    fig = Figure(size = 100 .* (4, 1.8))
    for (i, r) in enumerate(r_values)
        # A = 2.3329237880493565
        # r = 25.664665666632807
        ax = Axis(fig[1, i], xlabel = L"q_x", aspect = 1,ylabelvisible = i == 1, yticklabelsvisible = i == 1,
        xticks = PiTicks([0,pi]),yticks = PiTicks([0,pi]),
        xminorticks = IntervalsBetween(2),yminorticks = IntervalsBetween(2),
        xminorticksvisible = true, yminorticksvisible = true,
        )
        heatmap!(ax, k, k, [SqFieldTheory(x, y, 1,r) for x in k, y in k])
        
        # Add label to the plot
        Label(fig[1,i,TopLeft()],string(('a':'z')[i],")"), color = :black, fontsize = 14,padding = (-20,0,-10,0))
    end
    
    fig
end