
getR(ij::CartesianIndex{2}) = float(SA[ij[1],ij[2]])

function getRij_vec(Config::SpinConfig,i)
    CI = CartesianIndices(Config)
    ri = getR(CI[i])
    rij = [ri - getR(j) for j in CI]
end    

function getRij_vec(Config::SpinConfig)
    Ri = reshape([float(SVector(Tuple(ij))) for ij in CartesianIndices(Config.Mat)],length(Config))
    return [Ri[i] - Ri[j] for i in eachindex(Ri) for j in 1:i]
end

function getSij(Configs::AbstractVector{<:SpinConfig},i,j)
    return mean(c[i]*c[j] for c in Configs)
end

function getSij(Configs::AbstractVector{<:SpinConfig},i)
    return fetch.([Threads.@spawn getSij(Configs,i,j) for j in eachindex(Configs[1])])
end

function getSij(Configs::AbstractVector{<:SpinConfig})
    fac(i,j) = ifelse(i==j,1,2)
    return fetch.([Threads.@spawn fac(i,j)* getSij(Configs,i,j) for i in LinearIndices(Configs[1]) for j in 1:i])
end
