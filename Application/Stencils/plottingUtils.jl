using StaticArrays
using Optim
using MakieHelpers
using CairoMakie
import SpiderWebModel as SW
using SpiderWebModel.HDF5
dropmean(A; dims=:) = dropdims(mean(A; dims=dims); dims=dims)
dropstd(A; dims=:) = dropdims(std(A; dims=dims); dims=dims)

function errlines!(ax::MakieHelpers.Makie.AbstractAxis,x,y,err;bandkwargs = (;),kwargs...)
    l = lines!(ax,x,y;kwargs...)
    band!(ax,x,y .- err,y .+ err;color = (l.color[],0.3),bandkwargs...)
end
errlines!(x,y,err;bandkwargs = (;),kwargs...) = errlines!(current_axis(),x,y,err;bandkwargs,kwargs...)
errlines!(y,err;bandkwargs = (;),kwargs...) = errlines!(current_axis(),eachindex(y),y,err;bandkwargs,kwargs...)
errlines!(ax::MakieHelpers.Makie.AbstractAxis,y,err;bandkwargs = (;),kwargs...) = errlines!(ax,eachindex(y),y,err;bandkwargs,kwargs...)
function errlines(args...;axis = (;),kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1];axis...)
    errlines!(ax,args...;kwargs...)

    fig
end

function numberheatmap!(ax::MakieHelpers.Makie.AbstractAxis,x,y,z;plotheatmap=true,kwargs...)
    if plotheatmap
        heatmap!(ax,x,y,z;kwargs...)
    end
    for (i,xi) in enumerate(x), (j,yj) in enumerate(y)
        color = z[i, j] < mean(z) ? :white : :black
        text!(ax,xi,yj,text = string(z[i,j]),align = (:center, :center);color)
    end
end
numberheatmap!(args...;kwargs...) = numberheatmap!(current_axis(),args...;kwargs...)
function numberheatmap(args...;axis = (;),kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1];axis...)
    numberheatmap!(ax,args...;kwargs...)
    fig
end
numberheatmap!(z::AbstractMatrix;kwargs...) = numberheatmap!(axes(z,1),axes(z,2),z;kwargs...)
numberheatmap(z::AbstractMatrix;kwargs...) = numberheatmap(axes(z,1),axes(z,2),z;kwargs...)

_getkwargs(::Any) = (;xlabel = L"projection order $$")
_getkwargs(m::SW.ContinuousTimeMethod) = (;xlabel = L"\tau")
_getscaling(m::SW.DiscreteTimeMethod) = m.nBranch
_getscaling(i::Integer) = i
_getscaling(m::SW.ContinuousTimeMethod) = m.τ
function trueMomenta(kmin,kmax,L)
    nmin = floor(Int,L*kmin/(2pi))
    nmax = ceil(Int,L*kmax/(2pi))
    # return 1/(2pi*L*100) .* nmin:nmax
    return (nmin : nmax) .* 2pi/L
end
function roundToTrueMomentum(k,L)
    n = round(Int,k*L/(2pi))
    return n*2pi/L
end
function roundToTrueMomenta((kx,ky),L)
    nx = roundToTrueMomentum(kx,L)
    ny = roundToTrueMomentum(ky,L)
    return nx,ny
end

function plotEnergies(results,method,E0=NaN;normalize=true,axis=(;),kwargs...)
    ylabel = normalize ? L"E_0/L^2" : L"E_0"
    with_theme(theme_SimpleTicks()) do
        fig = Figure(fontsize = 22)
        ax = Axis(fig[1,1];ylabel ,xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5),_getkwargs(method)...,axis...)
        plotEnergies!(ax,results,method,E0;normalize,kwargs...)
        # ens = getfield.(obs,:E0)
        return fig
    end
end
_iscontinuous(m::Any) = false
_iscontinuous(m::SW.ContinuousTimeMethod) = true
function plotEnergies!(ax::Makie.Axis,results,method,E0=NaN;Emin=E0-1e-2,Emax=E0+2e-2,p=250,τ=nothing, nThermal=1,normalize=true,dense = _iscontinuous(method),legend = true,marker = '●',markersize = 5,label = L"GFMC$$",kwargs...)
    
    getnBra(i::Integer) = i
    getnBra(m::SW.AbstractGFMCMethod) = 1
    getnBra(m::SW.DiscreteTimeMethod) = m.nBranch
    nBra = getnBra(method)

    projSteps = p÷nBra

    if τ isa Real && method isa SW.ContinuousTimeMethod
        projSteps = round(Int,τ/method.τ)
    end

    ens = [SW.getEnergies(SW.get_TotalWeights(res),SW.get_energies(res),nThermal,projSteps) for res in results]
    en = mean(ens)
    # e0avg = minimum(en) -1e-5
    # en = log.(en .- e0avg)
    err = sqrt.(var(ens))
    if normalize
        NSpins = prod(size(results[1].SaveConfigs)[1:2])
        en = en ./NSpins
        err = err ./NSpins
        !isnan(E0) && (E0 = E0 /NSpins)
        !isnan(Emin) && (Emin = Emin /NSpins)
        !isnan(Emax) && (Emax = Emax /NSpins)
    end
    
    proj = _getscaling(method) .*(eachindex(en) .-1)
    if dense
        l = lines!(ax,proj,en,;label,color = :black,kwargs...)
        band!(ax,proj,en .- err,en .+ err;kwargs...,color = (l.color[],0.3))
    else
        scatterlines!(ax,proj,en;label,color = :black, marker,markersize ,kwargs...)
        errorbars!(ax,proj,en,err,whiskerwidth = 3.5,color = :black;kwargs...)
    end
    !isnan(E0)&& hlines!(ax,[E0],color = :red,label = L"exact $$")
    if legend
        axislegend(ax,merge=true,unique=true)
    end
    xlims!(ax,-0.5,last(proj))
    !isnan(Emin)&& !isnan(Emax) && ylims!(ax,Emin,Emax)
    return ax
end
plotEnergies!(results,method,E0=NaN;Emin=E0-1e-2,Emax=E0+2e-2,kwargs...) = plotEnergies!(current_axis(),results,method,E0;Emin=Emin,Emax=Emax,kwargs...)


function movingaverage(X::AbstractVector,numofele::Int)
    BackDelta = div(numofele,2) 
    ForwardDelta = isodd(numofele) ? div(numofele,2) : div(numofele,2) - 1
    len = lastindex(X)
    Y = similar(X)
    for n = eachindex(X)
        lo = max(firstindex(X),n - BackDelta)
        hi = min(len,n + ForwardDelta)
        @views Y[n] = mean(X[lo:hi])
    end
    return Y
end
function trackWalkerPath(reconfigTable,InitialWalkers,NSteps)
    walkerPopulation = zeros(Bool,size(reconfigTable,1),NSteps)
    currentWalkers = Set(InitialWalkers)
    nextWalkers = empty(currentWalkers)
    for step in 1:NSteps
        for walker in axes(reconfigTable,1)
            if reconfigTable[walker,step] in currentWalkers
                walkerPopulation[walker,step] = true
                push!(nextWalkers,walker)
            end
        end
        currentWalkers = copy(nextWalkers)
        empty!(nextWalkers)
    end
    return walkerPopulation
end

function equilib_plots(results;scatter_fraction,averageSteps = 100,Ntrack=50,p = 50,plotPopulation=false)
    fig = Figure(size = 100 .*(10, length(results)),theme = theme_SimpleTicks())
    initWalkers = collect(round(Int,scatter_fraction*size(results[1].reconfigurationTable,1)):size(results[1].reconfigurationTable,1))
    function getkw(i,title)
        if i == firstindex(results)
            return (;title,xlabelvisible = false,xticklabelsvisible = false)
        elseif i != lastindex(results)
            return (;xlabelvisible = false,xticklabelsvisible = false)
        end
        return (;)
    end

    axen = [Axis(fig[i,1],xlabel = "step";getkw(i,L"E_L/L^2")...) for i in eachindex(results)]
    axws = [Axis(fig[i,2],xlabel = "step";getkw(i,"weight")...) for i in eachindex(results)]
    axreconf = [Axis(fig[i,3],xlabel = "step";getkw(i,"walker heritage")...) for i in eachindex(results)]
    axSq = [Axis(fig[i,4],xlabel = L"q_x",aspect = 1,ylabel = L"q_y";xticks = PiTicks([0,pi,2pi]),yticks = PiTicks([0,pi,2pi]),getkw(i,L"\mathcal{S}(\mathbf{q})")...) for i in eachindex(results)]

    minweight,maxweight = (Inf,-Inf)
    en_min,en_max = (Inf,-Inf)
    
    markenergy = Inf
    markweight = -Inf
    confex = results[1].SaveConfigs[:,:,begin,begin]
    kx = 2pi .* LinRange(0,1,size(confex,1).+1)
    ky = 2pi .* LinRange(0,1,size(confex,2).+1)

    for (i,res) in enumerate(results)
        ws = movingaverage(res.TotalWeights,averageSteps)
        minweight = min(minimum(ws),minweight)
        maxweight = max(maximum(ws),maxweight)
        en = movingaverage(res.energies,averageSteps)
        en ./= length(confex)
        enAvg = mean(en[end÷4:end])
        markenergy = min(markenergy,enAvg)
        
        weightAvg = mean(ws[end÷4:end])
        markweight = max(markweight,weightAvg)

        en_min = min(minimum(en),en_min)
        en_max = max(maximum(en),en_max)
        lines!(axen[i],eachindex(en),en)
        lines!(axws[i],eachindex(ws),ws)

        hlines!(axen[i],[enAvg,],color = :black)

        hlines!(axen[i],[weightAvg,],color = :black)


        # WP = trackWalkerPath(res.reconfigurationTable,initWalkers,Ntrack)'
        if plotPopulation
            WP = SW.getBranchingMatrix(res.reconfigurationTable,Ntrack,Ntrack-1)
            heatmap!(axreconf[i],WP.PopulationMatrix',colormap = :jet)
        else
            WP = getBranchingHistory(res.reconfigurationTable,Ntrack)
            heatmap!(axreconf[i],WP',colormap = :jet)
        end
        hlines!(axreconf[i],minimum(initWalkers)-0.5,color = :black,linewidth = 1,linestyle = :dash)
        Sq = SW.getSqGFMC(res,p)
        heatmap!(axSq[i],kx,ky,Sq)
    end

    for (i,res) in enumerate(results)
        ylims!(axws[i],minweight,maxweight+1e-10)
        ylims!(axen[i],en_min,en_max+1e-10)
        hlines!(axen[i],markenergy,color = :red,linewidth = 1,linestyle = :dash)
        hlines!(axws[i],markweight,color = :red,linewidth = 1,linestyle = :dash)
    end

    fig
end

function getBranchingHistory(reconfigurationTable,NSteps)
    BranchingMatrix = fill(0,size(reconfigurationTable,1),NSteps)
    recView = @view reconfigurationTable[:,begin:NSteps]
    getBranchingHistory!(BranchingMatrix,recView)
end

function getBranchingHistory!(BranchingMatrix::AbstractMatrix,reconfigurationTable::AbstractMatrix)
    BranchingMatrix[:,begin] .= @view reconfigurationTable[:,begin]

    for n in axes(reconfigurationTable,2)[begin+1:end]
        for α in axes(reconfigurationTable,1)
            α´ = reconfigurationTable[α,n]
            BranchingMatrix[α,n] = BranchingMatrix[α´,n-1]
        end
    end
    return BranchingMatrix
end

function plotVarEn(stochReconfRes;normalization=1,movavg = 1,E_exact = NaN,alpha_index = 1,plotDiff = true)
    fig = Figure(theme = theme_SimpleTicks())

    ax = Axis(fig[1, 1], xlabel = "Iteration", ylabel = "Energy",xlabelvisible=false,xticklabelsvisible=false)
    axVar = Axis(fig[1,1], ylabel = "σE", yaxisposition=:right,yticklabelcolor=:red,ylabelcolor = :red,ygridstyle = :dash,xgridvisible = false,xticksvisible = false,xticklabelsvisible = false)

    ax2 = Axis(fig[2,1], xlabel = "Iteration", ylabel = "α")
    ax2norm = Axis(fig[2,1], ylabel = "||Δα||", yaxisposition=:right,yticklabelcolor=:red,ylabelcolor = :red,ygridstyle = :dash,xgridvisible = false,xticksvisible = false,xticklabelsvisible = false)
    linkxaxes!(ax,axVar, ax2,ax2norm)

    Epl =stochReconfRes.E0 ./ normalization


    x = eachindex(Epl)
    # errorbars!(ax,x,Epl,stochReconfRes.ΔE./ normalization,whiskerwidth=5,color = :black)
    Delta = stochReconfRes.ΔE./ normalization
    


    parsteps = stochReconfRes.params_steps
    parslice = parsteps[alpha_index,:]
    lastdim = length(size(stochReconfRes.params_steps))

    norms = let 
        if plotDiff
            pars = stochReconfRes.params_steps[1:100:end,:]
            diffs = diff(pars,dims = lastdim)
            SW.norm.(eachslice(diffs,dims = lastdim))
        else
            fill(NaN,length(x)-1)
        end
    end

    if movavg > 1
        Epl = movingaverage(Epl,movavg)
        Delta = movingaverage(Delta,movavg)
        norms = movingaverage(norms,movavg)
        parslice = movingaverage(parslice,movavg)
    end

    hlines!(ax,E_exact,color = :black,linestyle = :dash,linewidth = 2)

    band!(ax,x,Epl .-Delta, Epl .+ Delta,color = (:black,0.3))
    lines!(ax,x,Epl,color = :black)
    lines!(axVar,x,Delta,color = :red,linewidth = 1)
    lines!(ax2,x,parslice,color = :black)
    lines!(ax2norm,x[2:end],norms,color = :red)
    fig
end
function plotVarEn(filename::String;detectZero=true,kwargs...)
    E0 = h5read(filename,"E0")
    ΔE = h5read(filename,"ΔE")
    params_steps = h5read(filename,"params_steps")
    @views if detectZero
        idx = findfirst(iszero,E0)
        if !isnothing(idx)
            E0 = E0[1:idx-1]
            ΔE = ΔE[1:idx-1]
            params_steps = params_steps[:,1:idx-1]
        end
    end
    res = (;E0,ΔE,params_steps)
    plotVarEn(res;kwargs...)
end
function trueMomenta(kmin,kmax,L)
    nmin = floor(Int,L*kmin/(2pi))
    nmax = ceil(Int,L*kmax/(2pi))
    # return 1/(2pi*L*100) .* nmin:nmax
    return (nmin : nmax) .* 2pi/L
end
function reconstruct_momentumSpace(qs,arrvals;tol=1e-14)
    
    newarrvals = copy(arrvals)
    newqs = copy(qs)
    roundfunc(x) = round(x+1e-300,digits = -round(Int,log10(tol)))
    roundfunc(x::SVector) = roundfunc.(x)

    AllqsSet = Set(roundfunc.(qs))

    idx = 0
    
    for (q,val) in zip(newqs,newarrvals)
        idx +=1
        if idx > 5000
            @warn "reconstruct_momentumSpace: Could not reconstruct momentum space"
            return newqs,arrvals
        end

        qSymm = (SA[q[2],q[1]],SA[-q[1],q[2]],SA[q[1],-q[2]])
        for q´ in qSymm
            q´ = rem2pi.(q´,RoundNearest)
            q´round = roundfunc(q´)
            if !(q´round in AllqsSet)
                push!(newqs,q´)
                push!(newarrvals,val)
                push!(AllqsSet,q´round)
                # push!(newqsSet,round.(q´,digits = -round(Int,log10(tol))))
            end
            
        end

    end

    return newqs,newarrvals
end

function makeMatrix(qx,qy,qs,vals)
    FTrec = fill(NaN,length(qx),length(qy))
    
    for (i,qx) in enumerate(qx)
        for (j,qy) in enumerate(qy)
            idx = findfirst(==(SA[qx,qy]),qs)
            if !isnothing(idx)
                FTrec[i,j] = vals[idx]
            end
        end
    end
    return qx,qy,FTrec
end

uniqueTol(x;digits=16)  = unique(x -> round(x+1e-300;digits),x)
function makeMatrix(qs,vals) 
    qx = sort!(uniqueTol([q[1] for q in qs]))
    qy = sort!(uniqueTol([q[2] for q in qs]))
    # return sort!(unique(x->round(x+1e-300,digits=14),getindex.(qs,1)))
    makeMatrix(qx,qy,qs,vals)
end

function getLastSlice(arr::AbstractArray{T,N}) where {T,N}
    slicedims = tuple(collect(1 for i in 1:N-1)...)
    return view(arr,slicedims...,:)
end

function SqLargeN(qx,qy)
    cx = cos(qx)
    cy = cos(qy)
    sx = sin(qx)
    sy = sin(qy)
    c2x = cos(2*qx)
    c2y = cos(2*qy)

    num = 2*(cx -cy +2sx*sy)^2
    denom = 4 - 4*cx*cy - c2x*(1-2c2y) - c2y
    return num/denom
end

function SqFieldTheory_full(qx::Real, qy::Real, K::Real, W::Real, U::Real)
    cx = cos(qx)
    cy = cos(qy)
    sx = sin(qx)
    sy = sin(qy)
    
    numerator = (cx - cy + 2 * sx * sy)^2
    denominator = sqrt((cx - cy)^2 + 4 * sx^2 * sy^2) * sqrt(U / W + 4 * ((cx - cy)^2 + 4 * sx^2 * sy^2)) +1e-30
    
    return sqrt(K / (4 * W)) * numerator / denominator
end
SqFieldTheory(qx, qy, A,r) = SqFieldTheory_full(qx, qy, 4*A^2,1, r)
SqFieldTheory(q::AbstractVector,A::Real,r::Real) = SqFieldTheory(q[1],q[2],A,r)
SqFieldTheory(q::AbstractVector,coefs::AbstractVector) = SqFieldTheory(q[1],q[2],coefs[1],coefs[2])
SqFieldTheory(q::Tuple,coefs::AbstractVector) = SqFieldTheory(q[1],q[2],coefs[1],coefs[2])

function AsymFieldTheory_full(qx::Real, qy::Real, K::Real, W::Real, U::Real, p::Real)
    qx,qy = SA[
        cos(pi/2) sin(pi/2); 
        -sin(pi/2) cos(pi/2)
    ] * SVector(qx,qy)
    cx = cos(qx)
    cy = cos(qy)
    sx = sin(qx)
    sy = sin(qy)
    cxy = cos(qx+qy)

    numerator = sqrt(K)*(cx - cy + 2 * sx * sy)^2
    denominator1 = sqrt(U/4*(1+2p*(1+cxy)) + W*( (cx-cy)^2+4*sx^2*sy^2))
    denominator2 = sqrt((cx-cy)^2+4*sx^2*sy^2)
    return numerator / (denominator1 *denominator2 + 1e-30)
end
AsymFieldTheory(qx, qy, A,r,p) = AsymFieldTheory_full(qx, qy, 4*A^2,1, r,p)
AsymFieldTheory(q::AbstractVector,A::Real,r::Real,p) = AsymFieldTheory(q[1],q[2],A,r,p)
AsymFieldTheory(q::AbstractVector,coefs::AbstractVector) = AsymFieldTheory(q[1],q[2],coefs[1],coefs[2],coefs[3])
AsymFieldTheory(q::Tuple,coefs::AbstractVector) = AsymFieldTheory(q[1],q[2],coefs[1],coefs[2],coefs[3])


function optimizeCoeffs(SqMat,weightfunc=x->one(first(x)))
    q = trueMomenta(0., 2pi, size(SqMat, 1))[1:size(SqMat, 1)]

    function loss(v, w)
        l = 0.0
        v = abs(v)
        w = abs(w)
        for (i, qx) in enumerate(q), (j, qy) in enumerate(q)
            l += abs2(SqMat[i, j] - SqFieldTheory(qx, qy, v, w))*weightfunc(SA[qx,qy])
        end
        return l
    end

    loss(v) = loss(v[1], v[2])

    x0 = [1., 1.]

    res = optimize(loss, x0)
    @info res
    coefs = abs.(Optim.minimizer(res))
    return coefs
end
function optimizeCoeffsAsym(SqMat,weightfunc=x->one(first(x)))
    q = trueMomenta(0., 2pi, size(SqMat, 1) - 1)[1:end-1]

    function loss(v, w,r)
        l = 0.0
        v = abs(v)
        w = abs(w)
        r = abs(r)
        for (i, qx) in enumerate(q), (j, qy) in enumerate(q)
            l += abs2(SqMat[i, j] - AsymFieldTheory(qx, qy, v, w,r))*weightfunc(SA[qx,qy])
        end
        return l
    end

    loss(v) = loss(v[1], v[2],v[3])

    x0 = [1., 1.,1.]

    res = optimize(loss, x0)
    @info res
    coefs = abs.(Optim.minimizer(res))
    return coefs
end

function rasterCurve(curvePoints,grid,t)

    getPos(point) = findmin(x->SW.norm(SW.SVector(x.-point)),grid)[2]
    positions = getPos.(curvePoints)
    tnew = empty(t)
    posnew = empty(positions)
    for i in eachindex(t)
        p = positions[i]
        if p ∉ posnew
            push!(tnew, t[i])
            push!(posnew,p)
        end 
    end
    return tnew,posnew
end

function pointPath(p1::StaticArray,p2::StaticArray,res)
    Path = Vector{typeof(p1)}(undef,res)
    for i in eachindex(Path)
        Path[i] = p1 + i/res*(p2 -p1)
    end
    return Path
end
"""res contains the number of points along -pi,pi"""
function fetchKPath(points,res = 100)
    Path = Vector{typeof(points[begin])}(undef,0)
    # Path = []
    PointIndices = [1]
    for i in eachindex(points[begin:end-1])
        p1 = points[i]
        p2 = points[i+1]
        append!(Path,pointPath(p1,p2,round(Int,SW.norm(p1-p2)/2pi * res)))
        append!(PointIndices,length(Path)) # get indices corresponding to points
    end
    return PointIndices,Path
end

KPoints = Dict([
    "Γ" => SVector(0,0),
    "X" => SVector(pi,0),
    "M" => SVector(pi,pi),
    "X'" => SVector(0,pi)
    ])