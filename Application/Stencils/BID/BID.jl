function HammingDistance(a::AbstractArray{T},b::AbstractArray{T}) where T
    return sum(a[i] != b[i] for i in eachindex(a,b))
end

HammingDistance([0,0,1,0,1],[0,0,0,3,1])

using Combinatorics, HDF5

function P(d::Function, r,C=1.)
    return P(d(r), r,C)
end

function P(d_r::Real, r,C=1.)
    return C / (2. ^d_r) * binomial(d_r, r)
end

function hamming_distance_distribution(X::AbstractMatrix)
    n = size(X, 1)
    dist_counts = Dict{Int, Int}()
    
    @views for i in 1:n-1
        for j in i+1:n
            dist = HammingDistance(X[i, :], X[j, :])
            if haskey(dist_counts, dist)
                dist_counts[dist] += 1
            else
                dist_counts[dist] = 1
            end
        end
    end
    
    return dist_counts
end
Confs = h5read("/home/nniggema/Documents/Projects/SpiderWebModel.jl/Application/ConfsRaw/GurobiConfs/Confs_50_2.h5","Confs")
X = collect(values(Confs))
X = [reshape(x,prod(size(x)[1:2]),size(x,3)) for x in X]
X  = hcat(X...)
Xrand  = rand(Bool,size(X))
##
# Example usage:
dist_distribution = hamming_distance_distribution(X)
dist_ref = hamming_distance_distribution(Xrand)

using CairoMakie, LsqFit

function get_hamming_distance_distribution(dist_counts::Dict{Int, Int})
    dists = collect(keys(dist_counts))
    counts = collect(values(dist_counts))
    
    sorted_indices = sortperm(dists)
    dists = dists[sorted_indices]
    counts = counts[sorted_indices]

    total_pairs = sum(counts)
    probabilities = counts #./ total_pairs
    
    return dists, probabilities
end

with_theme(theme_SimpleTicks()) do
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=L"r", ylabel=L"P(r)")

    r1, p1 = get_hamming_distance_distribution(dist_distribution)
    r2, p2 = get_hamming_distance_distribution(dist_ref)

    model_d(r, d) = d[1] + d[2] * r
    model_P(r, d) = P(model_d(r,d), r,d[3])
    model_P(rv::AbstractVector, d) = [P(model_d(r,d), r) for r in rv]

    # model_P(rv::AbstractVector, d) = exp.(d[1] .* (rv .- d[2]))

    p0 = [500., 0.1,0.9]
    
    fit1 = curve_fit(model_P, r1, p1, p0)
    fit2 = curve_fit(model_P, r2, p2, p0)

    scatter!(ax, r1, p1, color = :blue, label = "Data")
    scatter!(ax, r2, p2, color = :red, label = "Reference")

    strd(x;kwargs...) = string(round(x;digits=4))

    C1,d0_1,d0_2 = fit1.param
    C2,d0_1,d0_2 = fit2.param

    # lines!(ax, r1, model_P(r1, p0), color = :blue, label = "Fit")
    lines!(ax, r1, model_P(r1, fit1.param), color = :blue, label = L"\mathrm{BID}= %$(strd(d0_1))")
    lines!(ax, r2, model_P(r2, fit2.param), color = :red, label = L"\mathrm{BID}_{\mathrm{rand}} $= %$(strd(d0_2))$")
    axislegend(ax)

    fig
end
