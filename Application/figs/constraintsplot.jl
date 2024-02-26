using CairoMakie
import SpiderWebModel as SW
##
let
    S = SW.SpinConfig(fill(0, 10, 10), 1)
    # fig = Figure()
    # ax = Axis(fig[1,1],
    # aspect = DataAspect(),
    # backgroundcolor = :white,
    # xminorgridcolor = :black,
    # yminorgridcolor = :black,
    # xminorgridvisible = true,
    # yminorgridvisible = true,
    # xgridvisible = false,
    # ygridvisible = false,
    # xminorticks = 0.5 .+ axes(S,1)[1:end-1],
    # yminorticks = 0.5 .+ axes(S,2)[1:end-1],
    # limits = (0.5,10.5,0.5,10.5),
    # )

    # S = SW.SpinConfig(fill(0,20,20),1)
    points = Point2f[]
    for i in axes(S, 1)
        for j in axes(S, 2)
            if iseven(i + j) && SW.plaquetteIsInBounds(S, i, j)
                border = Point.([
                    (i - 1.5, j - 1.5),
                    (i + 1.5, j - 1.5),
                    (i + 1.5, j + 1.5),
                    (i - 1.5, j + 1.5),
                    (i - 1.5, j - 1.5)
                ])
                # lines!(border,color = :red)
                # band!([i-1.5,i+1.5],[j-1.5,j-1.5],[j+1.5,j+1.5],color = (:black,0.2))
                # scatter!([Point(i,j),],color = :red,markersize = 10)
                push!(points, Point2f(i, j))
                Pij = SW.getPlaquette(S, i, j)
                Pij .+= 1
                Pij[2, 2] -= 1
            end
        end
    end
    fig = SW.plotSpinConfig(S, colorrange = (minimum(S), maximum(S)))
    scatter!(points, color = :red, markersize = 10)
    fig
    # fig = SW.plotSpinConfig(S)

end
