using LatticeFFTs
using LatticeFFTs
using LatticeFFTs.Interpolations

function getStructureFac(AllStates,eigen,tol=0)
    plan = getLatticeFFTPlan(AllStates[1].Mat,0)
    Psi = eigen.vectors[:,1]
    Nsites = length(AllStates[1])
    # weights = abs2.(Psi)
    # inds = findall(x->abs(x)>tol,weights)

    # state_weight = collect(zip( AllStates,weights))[inds]
    Sq = [getInterpolatedFFT(abs2(psin)* c.Mat,0,plan;Interpolation = BSpline(Constant())) for (c,psin) in zip( AllStates,Psi)]
    # Sq = [getInterpolatedFFT(weight* c.Mat,0,plan;Interpolation = BSpline(Constant())) for (c,weight) in state_weight]
    SSq(kx,ky) = real(sum(s(kx,ky)*s(-kx,-ky) for s in Sq))
end

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
