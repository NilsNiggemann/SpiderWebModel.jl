import SpiderWebModel as SW
using CairoMakie
function getPeriodicTilings(states, offset,Lx,Ly)
    isValidTiling = zeros(Bool, length(states))
    Threads.@threads for i in eachindex(states)
        isValidTiling[i] = SW.isPeriodicTiling(3Lx, 3Ly, states[i], SW.ALLGS_S12, offset)
    end
    return findall(isValidTiling)
end

function getPeriodicTilings(states::AbstractVector{<:SW.SpinConfig}, Lx,Ly,offset)
    isValidTiling = zeros(Bool, length(states))
    Threads.@threads for i in eachindex(states)
        isValidTiling[i] = SW.isPeriodicTiling(states[i], 5Lx, 5Ly, offset)
    end
    return findall(isValidTiling)
end
function getDensity(state_lazy,Lx,Ly,Offset)
    state = SW.getPeriodicState_history(4Lx, 4Ly, Lx, Ly, state_lazy, SW.ALLGS_S12, Offset)
    density = SW.getApplicablePlaquettes(state) |> length
end

function getDensity(state::SW.SpinConfig,Lx,Ly,Offset)
    state = SW.getPeriodicState(state, 70Lx, 70Ly, Offset)
    density = SW.getApplicablePlaquettes(state) |> length
end

##
function get_all_periodic_UCs(Lx,Ly)
    uniqueUCS = empty!(Set([SW.PeriodicMatrix(zeros(Float64,Lx,Ly),2,2,0)]))
    @time "constructing Configs" UCpth = SW.constructAllConfigs(Lx, Ly, SW.ALLGS_S12)
    @time "fill empty sites" UC_Full = SW.fillEmptyStates(UCpth, Lx, Ly, SW.ALLGS_S12)

    for Offset in 0:max(Lx,Ly)-1
        @time "Offset = $Offset" periodicTilings = getPeriodicTilings(UC_Full,Lx,Ly, Offset)
        
        # UCs = UC_Full[periodicTilings]
        UCs_periodic = [SW.PeriodicMatrix(parent(UC_Full[ind]), Lx, Ly, Offset) for ind in periodicTilings]
        # states = [SW.getPeriodicState(UC, Lx, Ly, Offset) for UC in UCs_periodic]
        
        # @time "filtering Confs" unique_inds = find_unique(uniqueUCS,UCs_periodic, 0.5)
        # @info "" eltype(uniqueUCS) eltype(filteredUCs)
        union!(uniqueUCS,UCs_periodic)
    end
    return uniqueUCS

end

function findEnlargedUC(permat)
    Lx,Ly = size(permat.UC)
    
    # newUC = Lx*(permat.offset+1),Ly
    newUC = Lx * (Ly ÷ gcd(permat.offset, Ly)),Ly
    # newUC = permat[1:Lx*(permat.offset),1:Ly]
end

function filter_nontrivial(all_UCs)
    nontrivial_UCs = Set(empty(all_UCs))
    for UC in all_UCs
        Lx,Ly = findEnlargedUC(UC)
        S = SW.stencilConfig( UC[1:Lx,1:Ly], 1 / 2,boundaryCondition = :periodic)
        if SW.getApplicablePlaquettes(S) |> length > 0
            push!(nontrivial_UCs, UC)
        end
    end
    return nontrivial_UCs
end

function test_enlargedUC(conf)
    c1 = conf[23:102,12:103]
    c2 = enlarge_auto(conf)[23:102,12:103]
    return c1 == c2
end

function enlarge_single(Conf,n1,n2,S=1/2)
    Lx,Ly = size(Conf)
    UC = SW.PeriodicMatrix(parent(Conf),n1*Lx,n2*Ly)
    return SC(Matrix(UC),S)
end
function squarify_target(Conf,Ltarget)
    Lx,Ly = size(Conf)
    # @info "" Lx Ly Ltarget Lx*ceil(Int,Ltarget / Lx) Ly*ceil(Int,Ltarget / Ly)

    n = ceil(Int,Ltarget / Lx) 
    m = ceil(Int,Ltarget / Ly)
    # return n,m
    return enlarge_single(Conf,n,m)
end
function squarify(Conf)
    Lx,Ly = size(Conf)
    
    Lmax = max(Lx,Ly)
    return squarify_target(Conf,Lmax)
end

function enlarge_auto(UC;Lxfac=1,Lyfac=1)
    Lx,Ly = findEnlargedUC(UC)
    Ly = Lyfac .* Ly
    Lx = Lxfac .* Lx
    SW.PeriodicMatrix(UC[1:Lx,1:Ly],Lx,Ly)
end

function filter_time_reverse!(all_UCs)
    for UC in all_UCs
        if -UC in all_UCs
            delete!(all_UCs,-UC)
        end
    end
    return all_UCs
end

function filter_translations!(all_UCs)
    for UC in all_UCs
        Lx,Ly = size(UC)
        for Tx in 2:2:Lx
            for Ty in 2:2:Ly
                UC_per = UC[begin+Tx:end , begin+Ty:end]
                if UC_per in all_UCs
                    delete!(all_UCs,UC_per)
                end
            end
        end
    end
    return all_UCs
end
function rot90(A::AbstractMatrix, k::Integer=1)
    k = mod(k, 4)
    k == 0 && return copy(A)
    k == 1 && return reverse(permutedims(A), dims=1)
    k == 2 && return reverse(reverse(A, dims=1), dims=2)
    k == 3 && return reverse(permutedims(A), dims=2)
end

rotations(X) = (
    X,
    rot90(X, 1),
    rot90(X, 2),
    rot90(X, 3)
)

function ymirror(X)
    return reverse(X, dims=1)
end
function xmirror(X)
    return reverse(X, dims=2)
end

mirrors(X) = (
    X,
    ymirror(X),
    xmirror(X),
    ymirror(xmirror(X))
)

rotations(X::SW.PeriodicMatrix) = rotations(X.UC)
function translate(X, dx, dy)
    circshift(X, (dx, dy))
end

function canonical(X::AbstractArray{T}) where T
    # reps = typeof(parent(X))[]
    # reps = Vector{T}[]
    mats = Matrix{T}[]
    # Stot = sum(X)
    if sum(X) < 0
        X = -X
    end
    
    function sortFunc(mat)
        # return vec(mat)
        return vec(mat)
        # Lx,Ly = size(mat)
        # sg = Lx >= Ly
        # return sg ? vec(mat) : vec(rot90(mat,1))
    end
    numSites = length(X)

    for R in rotations(X)
        for Rprime in mirrors(R)
            for dx in 0:size(X,1)-1, dy in 0:size(X,2)-1
                Xprime = translate(R, dx, dy)
                Xprime = reduce_to_minimal(Xprime)
                # SCPrime = SW.stencilConfig(Xprime, 1 / 2; boundaryCondition = :periodic)
                SCPrime = SC(Xprime)
                Lx,Ly = size(Xprime)
                # length(SCPrime) > numSites && continue
                Lx < Ly && continue # enforce Lx >= Ly for minimal representation
                if SW.fulFillsConstraint(SCPrime)

                    # if length(SCPrime) < numSites
                    #     numSites = length(SCPrime)
                    #     empty!(mats)
                    # end
                
                    # push!(reps, vec(Xprime))
                    push!(mats, Xprime)
                end
            end
        end
    end
    # return mats[argmin(vec,mats)]
    return argmin(sortFunc,mats)
    # return mats[argmin(reps)]
    # return reshape(minimum(reps),size(X))
end
##
function SC(PM,S=1/2)
    SW.stencilConfig(PM, S; boundaryCondition = :periodic)
end
function SC(PM::SW.PeriodicMatrix,S=1/2)
    SW.stencilConfig(PM.UC, S; boundaryCondition = :periodic)
end

using SpiderWebModel.HDF5
function get_all_UCs(maxLx,maxLy)
    allConfs = Dict{Matrix{Int8},Tuple{Int,Int}}()
    for Lx in 2:maxLx
        for Ly in 2:maxLy
            println("Lx = $Lx,Ly = $Ly")
            all_UCs = filter_nontrivial(get_all_periodic_UCs(Lx,Ly))
            isempty(all_UCs) && continue
            enlarged_UCS = Set(enlarge_auto.(all_UCs))
            println(all(x->SW.fulFillsConstraint(SC(x)), enlarged_UCS))
            @info "" unique(size.(enlarged_UCS))
            filter_time_reverse!(enlarged_UCS)

            canonical_UCS = fetch.([Threads.@spawn canonical(UC.UC) for UC in enlarged_UCS])
            @time canonical_UCS_unique = unique(canonical_UCS)
            @info "" length(canonical_UCS_unique)
            canonical_UCS_unique = SC.(canonical_UCS_unique)
            sort!(canonical_UCS_unique, by = x -> length(SW.getApplicablePlaquettes(x))/length(x),rev=true)
            
            for S in canonical_UCS_unique
                if SW.fulFillsConstraint(S)
                    allConfs[S] = (Lx,Ly)
                else
                    @warn "invalid UC found"
                end
            end
        end
    end
    return allConfs
end
function has_subcell(x, lx, ly; atol=1e-12)
    Lx, Ly = size(x)
    (Lx % lx != 0 || Ly % ly != 0) && return false

    @inbounds for j in 1:Ly, i in 1:Lx
        ii = (i - 1) % lx + 1
        jj = (j - 1) % ly + 1
        if abs(x[i,j] - x[ii,jj]) > atol
            return false
        end
    end
    return true
end
function divisors(n::Int)
    ds = Int[]
    for d in 1:floor(Int, sqrt(n))
        if n % d == 0
            push!(ds, d)
            d != n ÷ d && push!(ds, n ÷ d)
        end
    end
    sort!(ds)
    return ds
end
function find_subcell(x; atol=1e-12)
    Lx, Ly = size(x)

    for lx in divisors(Lx)
        for ly in divisors(Ly)
            (lx == Lx && ly == Ly) && continue
            has_subcell(x, lx, ly; atol=atol) && return (lx, ly)
        end
    end
    return Lx,Ly
end

function reduce_to_minimal(x)
    Lx,Ly = find_subcell(x)
    
    return x[1:Lx,1:Ly]
end

a = get_all_UCs(6,6)
##
allConfs = unique(squarify_target.(reduce_to_minimal.(collect(keys(a))),12))
maxFlipConfs = [SW.findMaxFlipConf(x;numRuns=3,tau=0.1,Nwalkers=10,NSteps = 1000) for x in allConfs]
# allConfs = unique(canonical.(squarify_target.(collect(keys(a)),20)))
maxFlipConfs = SC.(squarify_target.(unique(canonical.(maxFlipConfs)),18))
sort!(maxFlipConfs, by = x -> length(SW.getApplicablePlaquettes(x))/length(x),rev=true)
##
# trim_confs = fft_trim_equivalents(allConfs)
# sort!(trim_confs, by = x -> length(SW.getApplicablePlaquettes(SC(x)))/length(x),rev=true)
##
display.(SW.plotApplPlaquettes.(maxFlipConfs[1:20]))
##
function no_overlapping_plaquettes(Conf)
    AllPlaqs = SW.getApplicablePlaquettes(Conf)
    length(AllPlaqs) == 0 && return true

    for p_i in AllPlaqs
        for p_j in AllPlaqs
            if p_i != p_j && !SW.plaquettesAreSeparated(p_i, p_j)
                return false
            end
        end
    end
    return true
end

function is_quantum_trivial(Conf;verbose=false)
    S = SW.SpinConfig(Conf)
    S_p = copy(S)
    AllPlaqs = SW.getApplicablePlaquettes(S_p)
    length(AllPlaqs) == 0 && return true


    for p_i in AllPlaqs
        S_p .= S
        SW.flipPlaquette!(S_p, p_i)
        AllPlaqs_p = SW.getApplicablePlaquettes(S_p)
        if AllPlaqs_p != AllPlaqs
            if verbose 
                display(SW.plotApplPlaquettes(S_p))
                # display(S_p - Conf)
                println("Found changing plaquette at $p_i")
                println("Original plaquettes: $AllPlaqs")
                println("New plaquettes: $AllPlaqs_p")

            end
            return false
        end
    end
    return true
end

function get_en(conf, mu,  i = nothing)
    @assert mu==0. "finite mu not supported"
    if !isnothing(i)
        println("Config $i")
    end
    triv = is_quantum_trivial(conf)
    Nsites = length(conf)
    if triv
        NPlaqs = SW.getApplicablePlaquettes(conf) |> length
        return -(1. -mu)/Nsites * NPlaqs
    end
    try
        H = SW.generateHilbertSpace(SW.SpinConfig(conf))
        # SW.addRKPotential!(H, mu)
        ExSol = SW.SolveHKrylov(H.H)
        return ExSol.values[1]/Nsites
    catch e
        @warn "Error computing energy" exception = e
        return NaN
    end
end


function get_ens_GFMC(S,mu;Nwalkers = 20 * 3,NSteps = 2000,NRuns=20)
    println("mu = $mu")
    CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = 0.1*length(S),Hxx = SW.Hxx_RK(mu))
    ψG = SW.RKFunction()
    results_en = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S,CT,Nwalkers,NSteps,50,ψG,estimate_w_avg=false,equilibration_steps=NSteps÷5,pre_equilibration_steps=10_000, scatter_fraction = 1.0).Energy for _ in 1:NRuns])
    # @time results_en = SW.measure_Sq_GFMC(S,CT,Nwalkers,NSteps,50,ψG.psi,equilibration_steps=NSteps÷5, estimate_w_avg=true)
    # return results_en
    # display(plotEnergies(results,CT))
    ens = stack(results_en)[end,:]

    return ens
end

import Statistics: mean, std

function e_trivial(conf, mu)
    Nsites = length(conf)
    NPlaqs = SW.getApplicablePlaquettes(conf) |> length
    return -(1. -mu) * NPlaqs/Nsites
end

function get_e_GFMC(conf, mu,  i = nothing;kwargs...)
    if !isnothing(i)
        println("Config $i")
    end
    # triv = is_quantum_trivial(conf)
    Nsites = length(conf)
    # if triv
    #     NPlaqs = SW.getApplicablePlaquettes(conf) |> length
    #     return -(1. -mu) * NPlaqs/Nsites,0
    # end
    try
        ens = get_ens_GFMC(conf,mu;kwargs...) 
        return mean(ens)/Nsites, std(ens)/Nsites
    catch e
        @warn "Error computing energy" exception = e
        return NaN
    end
end

##
ens = get_en.(maxFlipConfs,0.0,eachindex(maxFlipConfs))
##
ens_GFMC = get_e_GFMC.(maxFlipConfs,0.8,eachindex(maxFlipConfs),Nwalkers = 50,NRuns = 6)
# ens_GFMC = get_e_GFMC.(allConfs,0.0,eachindex(allConfs);Nwalkers = 80)
##
function errscatter!(ax,x,ytup;markersize=6,kwargs...)
    y, err = first.(ytup), last.(ytup)
    errorbars!(ax, x,y, err;whiskerwidth = 5,kwargs...,label = nothing)
    scatter!(ax, x,y, err;markersize,kwargs...)
end
errscatter!(x,ytup;kwargs...) = errscatter!(current_axis(), x, ytup;kwargs...)

using MakieHelpers
with_theme(theme_SimpleTicks()) do 
    fig = Figure()
    ax = Axis(fig[1,1],xlabel = L"Sector$$", ylabel = L"E/L^2")
    trivs = is_quantum_trivial.(maxFlipConfs) 
    inds = eachindex(maxFlipConfs)
    errscatter!(inds[.!trivs], ens_GFMC[.!trivs]; label = L"nontrivial $$", color = :red,markersize = 8)
    errscatter!(inds[trivs], ens_GFMC[trivs]; label = L"trivial $$", color = :black,markersize = 8)

    e_trivs = e_trivial.(maxFlipConfs[trivs],0.8)
    scatter!(inds[trivs], e_trivs; label = L"trivial exact $$", color = (:black,0.4),markersize = 30,marker = '×' )
    # scatter!(inds[trivs], e_trivs; label = L"trivial exact $$", color = :blue,markersize = 6,marker = '◯')

    # scatter!(ens, label = "Exact diagonalization", color = :blue,markersize = 6)
    hlines!(-(1-0.8)*1/10,color = :black, linestyle = :dash)
    axislegend(position = :rb)
    save("figs/PaperFigs/spin_half_sectors_ens_mu08.pdf", fig)
    fig
end
##
mins = findall(<(-0.016),first.(ens_GFMC))

