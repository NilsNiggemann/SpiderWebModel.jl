using CairoMakie
using CairoMakie.Makie.ColorSchemes
function plotSpinConfig!(ax, S::AbstractSpinConfig; plotConstraints = true, kwargs...)
    vals = filter!(x -> !ismissing(x) && !isnan(x), unique(S.Mat))
    isempty(vals) && (vals = [-1, 1])
    Amin = min(minimum(vals), -getSpin(S))
    Amax = max(maximum(vals), getSpin(S))
    us = Amin:Amax
    # hm = heatmap!(ax,Array(S.Mat),colorrange = (Amin,Amax),colormap = cgrad(:linear_bgy_10_95_c74_n256, length(us), categorical = true);kwargs...)
    hm = heatmap!(
        ax,
        Array(parent(S)),
        colorrange = (Amin, Amax),
        colormap = cgrad(:grays, length(us), categorical = true);
        kwargs...,
    )
    if plotConstraints
        points =
            [Point(Tuple(I)...) for I in CartesianIndices(S.Mat) if iseven(sum(Tuple(I)))]
        scatter!(ax, points, marker = '×', color = :gray, markersize = 20)
    end
    translate!(hm, 0, 0, -100)
    return hm
end

function plotSpinConfig!(ax, S::AbstractSpinConfig{Bool}; kwargs...)
    plotSpinConfig!(ax, floatSpinConfig(S); kwargs...)
end

function getConfigAxis(S; kwargs...)
    (;
        aspect = DataAspect(),
        backgroundcolor = :grey,
        # xminorgridwidth = 2,
        # yminorgridwidth = 2,
        xminorgridcolor = :black,
        yminorgridcolor = :black,
        xminorgridvisible = true,
        yminorgridvisible = true,
        xgridvisible = false,
        ygridvisible = false,
        # xticks = (axes(S,1),string.(axes(S,1))) ,
        # yticks = (axes(S,2),string.(axes(S,2))) ,
        xminorticks = 0.5 .+ axes(S, 1),
        yminorticks = 0.5 .+ axes(S, 2),
        limits = (0.5, size(S, 1) + 0.5, 0.5, size(S, 2) + 0.5),
    )
end

function plotSpinConfig(S; kwargs...)
    fig = Figure(size = 1.2 .* (450, 400))
    ax = Axis(fig[1, 1]; getConfigAxis(S)...)
    hm = plotSpinConfig!(ax, S; kwargs...)
    # us = filter!(x -> !ismissing(x) && !isnan(x), unique(S))
    us = -S.S:S.S
    # isempty(us) && (us = [-1, 1])
    Colorbar(fig[1, 2], hm, ticks = us)
    return fig
end
