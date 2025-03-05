using HDF5, CairoMakie, LsqFit, MakieHelpers
using SpecialFunctions

logbinomial(n,k) = logfactorial(n)-logfactorial(k)-logfactorial(n-k)

function HammingDistance(a::AbstractArray{T},b::AbstractArray{T}) where T
    return sum(a[i] != b[i] for i in eachindex(a,b))
end

HammingDistance([0,0,1,0,1],[0,0,0,3,1])


function P(d::Function, r,C=1.)
    return P(d(r), r,C)
end

# function P(d_r::Integer, r,C=1.)
#     return C / (2. ^d_r) * binomial(d_r, r)
# end
# function P(d_r::Real, r,C=1.)
#     return C / (2. ^d_r) * binomial(d_r, r)
# end

function P(d_r::Integer, r,C=1.)
    log_result = log(C) - d_r * log(2) + logbinomial(d_r, r)
    return exp(log_result)
end

# P(d_r::Real, r,C=1.) = P(d_r, r, C)
P(d_r::Real, r,C=1.) = P(round(Int, d_r), r, C)

strd(x;kwargs...) = string(round(x;digits=4))

function hamming_distance_distribution(X::AbstractMatrix)
    dist_counts = Dict{Int, Int}()
    
    @views for i in axes(X,2)
        for j in axes(X,2)
            j >= i && continue
            dist = HammingDistance(X[:,i], X[:,j])
            if haskey(dist_counts, dist)
                dist_counts[dist] += 1
            else
                dist_counts[dist] = 1
            end
        end
    end
    
    return dist_counts
end

function collect_confs_from_folder(folder_path::String)
    conf_files = readdir(folder_path,join=true)
    all_confs = []

    for file in conf_files
        if endswith(file, ".h5")
            try
                confs = h5read(file, "Confs")
                push!(all_confs, collect(values(confs)))
                
            catch
            end
        end
    end

    return all_confs
end

Confs = collect_confs_from_folder("../ConfsRaw/GurobiConfs")
X = [reshape(xx, prod(size(xx)[1:2]), size(xx, 3)) for x in Confs for xx in x]
X = Array(hcat(X...))
Xrand = rand(Bool, size(X))
# Example usage:
dist_distribution = hamming_distance_distribution(X)
dist_ref = hamming_distance_distribution(Xrand)
##
function get_hamming_distance_distribution(dist_counts::Dict{Int, Int})
    dists = collect(keys(dist_counts))
    counts = collect(values(dist_counts))
    
    sorted_indices = sortperm(dists)
    dists = dists[sorted_indices]
    counts = counts[sorted_indices]

    total_pairs = sum(counts)
    probabilities = counts ./ total_pairs
    
    return dists, probabilities
end

with_theme(theme_SimpleTicks()) do
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=L"r", ylabel=L"P(r)")

    r1, p1 = get_hamming_distance_distribution(dist_distribution)
    
    r2, p2 = get_hamming_distance_distribution(dist_ref)
    function getParams(d)
        BID = d[1]
        # BID = 1800
        d2 = d[2]
        C = d[3]
        return (;BID, d2, C)
    end

    function model_P(r, d)
        (;BID, d2, C) = getParams(d)
        return P(BID +d2*r, r,C)
    end

    model_P(rv::AbstractVector, d) = [model_P(r,d) for r in rv]
    # return model_P
    NSpins = size(X,1)
    # model_P(rv::AbstractVector, d) = exp.(d[1] .* (rv .- d[2]))

    params_start_Rand = [2500, 0.0,1.]
    params_start_Corr = [1800.5, 0.56,0.1]
    # p_0 = [NSpins, 0.01,1]
    
    
    scatter!(ax, r1, p1, color = :blue, label = L"data$$",markersize=3)
    scatter!(ax, r2, p2, color = :red, label = L"Random $N = %$(size(X,1))$",markersize=3)
    # lines!(ax, r2, model_P(r2, p_01), color = :red)
    # lines!(ax, r1, model_P(r1, params_start_Corr), color = :blue, label = "Fit")
    # return fig
    fit1 = curve_fit(model_P, r1, p1, params_start_Corr)
    # fit2 = curve_fit(model_P, r2, p2, params_start_Rand)

    # return fit1
    (;BID, d2, C) = getParams(fit1.param)
    println(fit1.param)
    # d0_2,d1_2,C_2 = fit2.param

    lines!(ax, r1, model_P(r1, fit1.param), color = :blue, label = L"\mathrm{BID}= %$(strd(BID))")
    lines!(ax, r2, model_P(r2, params_start_Rand), color = :red)
    axislegend(ax, position = :lt)

    fig
end
