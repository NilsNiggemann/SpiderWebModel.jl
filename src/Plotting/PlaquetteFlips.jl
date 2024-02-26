function plotApplPlaquettes!(ax,State;square = (:green,0.0),heatmapkwargs=(;),kwargs...)
    plaqs = getApplicablePlaquettes(State)
    points = Point2f.(plaqs)
    plotSpinConfig!(ax,State;heatmapkwargs...)
    
    if square === true
        square = (:green,0.4)
    end
    
    if square[2] != 0.0
        for p in points
            px,py = p
            # band!(ax,[px-1.5,px+1.5],[py-1.5,py-1.5],[py+1.5,py+1.5],color = square)
            band!(ax,[px-1.5,px+1.5],[py+0.5,py+0.5],[py+1.5,py+1.5],color = square)
            band!(ax,[px-1.5,px+1.5],[py-1.5,py-1.5],[py-0.5,py-0.5],color = square)
            band!(ax,[px-1.5,px-0.5],[py-0.5,py-0.5],[py+0.5,py+0.5],color = square)
            band!(ax,[px+0.5,px+1.5],[py-0.5,py-0.5],[py+0.5,py+0.5],color = square)
            # plapoints = [Point2(-1,-1),Point2(-1,0),Point2(-1,1),Point2(0,1),Point2(0,-1),Point2(1,1),Point2(1,0),Point2(1,-1)] .+ p

            # scatter!(ax,plapoints,markersize = 30,marker = :rect,color = square;kwargs...)
        end
    end
    scatter!(ax,points,markersize = 13,color = :red;kwargs...)

end

function plotApplPlaquettes(State;heatmapkwargs = (;),kwargs...)
    fig = plotSpinConfig(State;heatmapkwargs...)
    plotApplPlaquettes!(current_axis(),State;kwargs...)
    fig
end
