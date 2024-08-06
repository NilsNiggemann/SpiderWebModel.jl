_getkwargs(::Any) = (;xlabel = L"projection order $$")
_getkwargs(m::SW.ContinuousTimeMethod) = (;xlabel = L"\tau")
_getscaling(m::SW.DiscreteTimeMethod) = m.nBranch
_getscaling(i::Integer) = i
_getscaling(m::SW.ContinuousTimeMethod) = m.nBranch*m.τ
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

function plotEnergies(results,method,E0=NaN;normalize=false,axis=(;),kwargs...)
    ylabel = normalize ? L"E_0/L^2" : L"E_0"
    with_theme(theme_SimpleTicks()) do
        fig = Figure(fontsize = 22)
        ax = Axis(fig[1,1];ylabel ,xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5),_getkwargs(method)...,axis...)
        plotEnergies!(ax,results,method,E0;normalize,kwargs...)
        # ens = getfield.(obs,:E0)
        return fig
    end
end

function plotEnergies!(ax::Makie.Axis,results,method,E0=NaN;Emin=E0-1e-2,Emax=E0+2e-2,p=250,τ=nothing, nThermal=1,normalize=false,dense = false,legend = true,marker = '●',markersize = 5,label = L"GFMC$$",kwargs...)
    
    getnBra(i::Integer) = i
    getnBra(m::SW.AbstractGFMCMethod) = m.nBranch
    nBra = getnBra(method)

    projSteps = p÷nBra

    if τ isa Real && method isa SW.ContinuousTimeMethod
        projSteps = round(Int,τ/method.τ/method.nBranch)
    end

    ens = [SW.getEnergies(res.TotalWeights,res.energies,nThermal,projSteps) for res in results]
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
    fig = Figure(size = 0.8 .*(1000, 1200),theme = theme_SimpleTicks())
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

function plotVarEn(stochReconfRes;normalization=1)
    fig = Figure(theme = theme_SimpleTicks())

    ax = Axis(fig[1, 1], xlabel = "Iteration", ylabel = "Energy",xlabelvisible=false,xticklabelsvisible=false)
    ax2 = Axis(fig[2,1], xlabel = "Iteration", ylabel = "α")

    Epl =stochReconfRes.E0 ./ normalization
    x = eachindex(Epl)
    errorbars!(ax,x,Epl,stochReconfRes.ΔE./ normalization,whiskerwidth=5)
    lines!(ax,x,Epl)
    parsteps = stochReconfRes.params_steps
    parslice = getLastSlice(parsteps)

    lines!(ax2,x,parslice)
    fig
end
function trueMomenta(kmin,kmax,L)
    nmin = floor(Int,L*kmin/(2pi))
    nmax = ceil(Int,L*kmax/(2pi))
    # return 1/(2pi*L*100) .* nmin:nmax
    return (nmin : nmax) .* 2pi/L
end
function getLastSlice(arr::AbstractArray{T,N}) where {T,N}
    slicedims = tuple(collect(1 for i in 1:N-1)...)
    return view(arr,slicedims...,:)
end