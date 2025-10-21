using Interpolations
using Optim

function getxi(Sq,I::CartesianIndex,dI::CartesianIndex)
    L = size(Sq,1)-1
    if Sq[I] < Sq[I+dI]
        return 0
    end
    xi_L = sqrt(Sq[I]/Sq[I+dI] -1 )
    return xi_L * L
end

function getxi(Sq,I::CartesianIndex)
    # neighbors = [(-1,-1)]
    neighbors = [(i,j) for i in -1:1,j in -1:1 if i != 0 || j != 0]
    # neighbors = [(1,0),(0,1),(-1,0),(0,-1)]
    return maximum(getxi(Sq,I,CartesianIndex(dI)) for dI in neighbors)
end

function getxi(Sq)
    # L = size(Sq,1)-1
    I = argmax(Sq)
    # return Sq[I] / L^2
    return getxi(Sq,I)
end

function getSq_tau(Sq::Array{T,3},tau::Real,Dtau::Real) where T
    if isinf(tau)
        for i in reverse(axes(Sq,3))
            S = @view Sq[:,:,i]
            any(isnan,S) && continue
            return S
        end
    else
        tauInd = ceil(Int,tau/Dtau)
    end
    return @view Sq[:,:,tauInd]
end

##
function getXis(Sqs)
    xis = zeros(size(Sqs,4))
    for (i,ii) in enumerate(axes(Sqs,4))
        Sq = @views dropmean(Sqs[:,:,:,ii],dims=3)

        xis[i] = getxi(Sq./4)
    end
    return xis
end
function getXis_err(Sqs,I)
    xis = zeros(size(Sqs)[3:4])
    for (i,ii) in enumerate(axes(Sqs,4))
        for (j,jj) in enumerate(axes(Sqs,3))
            Sq = @views Sqs[:,:,jj,ii]

            xis[j,i] = getxi(Sq./4,I)
        end
    end
    ximean = dropmean(xis,dims=1)
    xistd = dropstd(xis,dims=1)
    return (;ximean,xistd)
end

function getXis(Sqs::AbstractVector{<:AbstractMatrix},I)
    xis = zeros(length(Sqs))
    for (i,Sq) in enumerate(Sqs)
        xis[i] = getxi(Sq./4,I)
    end
    return xis
end


function detect_crossings(mus::Dict,xis::Dict)

    xis_intPol = Dict(L=> interpolate((mus[L],), xis[L], Gridded(Linear())) for L in keys(xis))

    Ls = sort(collect(keys(xis_intPol)))
    crossings = zeros(length(Ls)-1)
    yvals = zeros(length(Ls)-1)
    for i in eachindex(Ls)[1:end-1]
        L1 = Ls[i]
        L2 = Ls[i+1]
        xi_L1 = xis_intPol[L1]
        xi_L2 = xis_intPol[L2]

        mu_c = Optim.optimize(mu -> abs(xi_L1(mu) - xi_L2(mu)),0.0,1.0).minimizer
        crossings[i] = mu_c
        yvals[i] = xi_L1(mu_c)
    end

    return (;L = Ls[1:end-1],crossings,yvals)
end


function getXiLs(res::DataFrame,Q = pi/2)
    
    xis = Dict{Int,Vector{Float64}}()
    xis_err = Dict{Int,Vector{Float64}}()
    mus = Dict{Int,Vector{Float64}}()

    for L in unique(res.L)
        unique_mus = unique(res.mu)
        k = trueMomenta(0,2pi,L)
        i_k = findfirst(==(Q),k)
        Sq = [getSq(res,mu = mu,L = L,tau=10) for mu in unique_mus]
        xi_L = [getXis(Sq,CartesianIndex(i_k,i_k)) for Sq in Sq]
        xis[L] = mean.(xi_L) ./ L
        xis_err[L] = std.(xi_L) ./ L
        mus[L] = unique_mus
    end
end

function getXiCol(res::DataFrame,Q = (pi/2,pi/2);tau = 10)
    # return 
    
    
    xiL = map(zip(res.L,res.Sq,res.tau)) do (L,Sq,dtau)
        # println(typeof(Sq))
        k = trueMomenta(0,2pi,L)
        
        i_k = findfirst(==(Q[1]),k)
        j_k = findfirst(==(Q[2]),k)

        Sqtau = getSq_tau(Sq,tau,dtau)
        xi = getxi(Sqtau,CartesianIndex(i_k,j_k))
    end
end

function getSqMaxCol(res::DataFrame,Q = nothing;tau = 10)
    # return 
    
    Sqmax = map(zip(res.L,res.Sq,res.tau)) do (L,Sq,dtau)
        # println(typeof(Sq))
        k = trueMomenta(0,2pi,L)
        
        
        Sqtau = getSq_tau(Sq,tau,dtau)
        if isnothing(Q)
            return maximum(Sqtau)
        else
            i_k = findfirst(==(Q[1]),k)
            j_k = findfirst(==(Q[2]),k)
            return Sqtau[i_k,j_k]
        end
    end
end