function getCharges(Conf::SpinConfig, flipParity = false)
    Charges = zeros(size(Conf.Mat))

    for i in axes(Conf.Mat, 1), j in axes(Conf.Mat, 2)
        iseven(i + j + flipParity) || continue
        plaquetteIsInBounds(Conf, i, j) || continue
        P = getPlaquette(Conf, i, j)
        any(isnan, P) && continue

        c = constraint(P)

        if c ≠ 0
            Charges[i, j] = c
        end
    end
    return Charges
end

function getPlotPoints(State)
    ChargeMat = getCharges(State)

    col(i, j) = ChargeMat[i, j] > 0 ? :red : :blue
    size(i, j) = abs(ChargeMat[i, j]) * 18

    cols = [
        col(i, j) for
        i in axes(ChargeMat, 1), j in axes(ChargeMat, 2) if ChargeMat[i, j] ≠ 0
    ]
    sizes = [
        size(i, j) for
        i in axes(ChargeMat, 1), j in axes(ChargeMat, 2) if ChargeMat[i, j] ≠ 0
    ]

    points = [
        Point2(i, j) for
        i in axes(ChargeMat, 1), j in axes(ChargeMat, 2) if ChargeMat[i, j] ≠ 0
    ]

    return (; cols, sizes, points)
end

function getPlotPoints_all(State)
    ChargeMat = getCharges(State)

    col(i, j) = ChargeMat[i, j] > 0 ? (:red, 1) : (:blue, 1)
    size(i, j) = abs(ChargeMat[i, j]) * 15

    cols = [col(i, j) for i in axes(ChargeMat, 1), j in axes(ChargeMat, 2)][:]
    sizes = [size(i, j) for i in axes(ChargeMat, 1), j in axes(ChargeMat, 2)][:]

    points = [Point2(i, j) for i in axes(ChargeMat, 1), j in axes(ChargeMat, 2)][:]

    return (; cols, sizes, points)
end

function plotFractons!(ax, State; heatmapkwargs = (;), kwargs...)
    cols, sizes, points = getPlotPoints(State)
    plotSpinConfig!(ax, State; heatmapkwargs...)

    scatter!(ax, points, markersize = sizes, color = cols, marker = :rect; kwargs...)
    # scatter!(ax,points;color = cols,kwargs...)

end

function plotFractons(State; heatmapkwargs = (;), kwargs...)
    fig = plotSpinConfig(State; heatmapkwargs...)
    plotFractons!(current_axis(), State; kwargs...)
    fig
end
