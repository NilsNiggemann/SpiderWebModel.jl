module SpiderWebModel
    using StaticArrays,Random,Statistics,LoopVectorization
    using OrderedCollections, LinearAlgebra
    using ElasticArrays
    using Test
    import Base:size,getindex,setindex!,iterate,show,copy,hash

    struct SpinConfig{T,MatType <: AbstractMatrix{T},T1<:Real} <:AbstractMatrix{T}
        Mat::MatType
        S::T1
    end

    Base.getindex(S::SpinConfig,i,j) = getindex(S.Mat,i,j)
    Base.setindex!(S::SpinConfig,x,i,j) = setindex!(S.Mat,x,i,j)
    Base.iterate(S::SpinConfig,i) = iterate(S.Mat,i)
    Base.iterate(S::SpinConfig) = iterate(S.Mat)

    Base.size(S::SpinConfig) = size(S.Mat)
    Base.copy(S::SpinConfig) = SpinConfig(copy(S.Mat),S.S)
    # function Base.hash(S::SpinConfig)
    #     S.S != 1/2 && return hash((S.Mat,S.S))
    # end
    # Base.show(io::IO, ::MIME"text/plain", S::SpinConfig) = show(io,S.Mat)
    # Base.show(io::IO, S::SpinConfig) = show(io,S.Mat)
    
    constraintSigns() = (1,1,-1,-1,1,1,-1,-1)
    function constraint(sites)
        cons = constraintSigns()
        return sum(sgn*S for (sgn,S) in zip(cons,sites))
    end

    """Order in which the spins are stored in the plaquette corresponding to the numbering convention"""
    plaquetteOrder(S1,S2,S3,S4,S5,S6,S7,S8,S9) = SMatrix{3,3}(S4,S5,S6,S3,S9,S7,S2,S1,S8)

    plaquetteOrder(S1,S2,S3,S4,S5,S6,S7,S8) = SMatrix{3,3}(S4,S5,S6,S3,NaN,S7,S2,S1,S8)
    plaquetteOrder(x) = plaquetteOrder(x...)

    function plaquette(sites,S=1/2)
        Mat = SMatrix{3,3}(plaquetteOrder(sites))
        return SpinConfig(Mat,S)
    end

    function getSitesFromPlaquette(Mat::AbstractMatrix)
        return (Mat[2,3],Mat[1,3],Mat[1,2],Mat[1,1],Mat[2,1],Mat[3,1],Mat[3,2],Mat[3,3])
    end

    function getSitesFromPlaquette(P::SpinConfig) 
        return getSitesFromPlaquette(P.Mat)
    end

    function constraint(P::SpinConfig)
        return constraint(getSitesFromPlaquette(P))
    end

    function getPlaquette(S::SpinConfig,i,j)
        Mat = @view S.Mat[i-1:i+1,j-1:j+1]
        return SpinConfig(Mat,S.S)
    end

    """generates all possible combinations 0,1 of the 8 spins in a plaquette """

    function getAllGS(S)
        allCombs = [plaquette((i,j,k,l,m,n,o,p),S) for i in -S:S for j in -S:S for k in -S:S for l in -S:S for m in -S:S for n in -S:S for o in -S:S for p in -S:S if constraint((i,j,k,l,m,n,o,p)) == 0]
        return allCombs
    end

    function getAllGS_noMissing(S)
        allCombs = [plaquette((i,j,k,l,m,n,o,p,q),S) for i in -S:S for j in -S:S for k in -S:S for l in -S:S for m in -S:S for n in -S:S for o in -S:S for p in -S:S for q in -S:S if constraint((i,j,k,l,m,n,o,p)) == 0]
        return allCombs
    end


    const ALLGS_S12 = getAllGS(0.5)
    const ALLGS_S12_NOMISSING = getAllGS_noMissing(0.5)

    _isZeroOrNaN(x) = x == 0 || isnan(x)
    # _isZeroOrMissing(x) = x === missing || x == 0

    # function addTile!(Pij::SpinConfig,Tile::SpinConfig,)
    #     for i in eachindex(Tile,Pij)
    #         ti = Tile[i]
    #         isnan(ti) && continue
    #         Pij[i] = ti
    #     end
    #     return Pij
    # end

    function addTile!(Pij::SpinConfig,Tile::SpinConfig,)
        Pij[1,1] = Tile[1,1]
        Pij[1,2] = Tile[1,2]
        Pij[1,3] = Tile[1,3]
        
        Pij[2,1] = Tile[2,1]
        #middle of tile is empty
        Pij[2,3] = Tile[2,3]
        
        Pij[3,1] = Tile[3,1]
        Pij[3,2] = Tile[3,2]
        Pij[3,3] = Tile[3,3]
        return Pij
    end


    function plaquetteIsInBounds(Conf::AbstractMatrix,iCenter::Integer,jCenter::Integer)
        for i in iCenter-1:iCenter+1, j in jCenter-1:jCenter+1
            if !checkbounds(Bool,Conf,i,j)
                return false
            end
        end
        return true
    end
    plaquetteIsInBounds(Conf::SpinConfig,iCenter::Integer,jCenter::Integer) = plaquetteIsInBounds(Conf.Mat,iCenter,jCenter)

    function allSpinsInBounds(Conf::SpinConfig;verbose=false)
        for x in Conf.Mat
            isnan(x) && continue
            if abs(x) > Conf.S 
                verbose && println("Spin larger than S")
                return false
            end
        end
        return true
    end

    """Assumes that constraint are only defined on every alternating site, starting from the first index"""
    function fulFillsConstraint_nonStrict(Conf::SpinConfig,flipParity=false;verbose = false)

        for i in axes(Conf.Mat,1), j in axes(Conf.Mat,2)
            iseven(i+j+flipParity) || continue 
            plaquetteIsInBounds(Conf,i,j) || continue
            P = getPlaquette(Conf,i,j)
            any(isnan,P) && continue

            c = constraint(P)

            if c ≠ 0
                verbose && println("Constraint not fulfilled at i,j = $((i,j))" )
                return false
            end
        end
        
        return true
    end

    function fulFillsConstraint(Conf::SpinConfig,flipParity=false;verbose = false)
        allSpinsInBounds(Conf;verbose) || return false
        return fulFillsConstraint_nonStrict(Conf,flipParity;verbose)
    end

    # function canPlaceTile(P::SpinConfig,T1)
    #     return all(_isZeroOrNaN, p -t for (p,t) in zip(P,T1))
    # end

    function canPlaceTile(P::SpinConfig,T1)
        for (p,t) in zip(P,T1)
            x = p - t
            if !(_isZeroOrNaN(x))
                return false
            end
        end
        return true
    end

    function canPlaceTile(P::Tuple,T1)
        Tsites = getSitesFromPlaquette(T1)
        return all(_isZeroOrNaN, p -t for (p,t) in zip(P,Tsites))
    end

    
    function TileFulFillsConstraint(S::SpinConfig,P,T1)
        OldP = SMatrix{3,3}(P)
        P .= T1
        cons = fulFillsConstraint_nonStrict(S)
        P .= OldP
        return cons
    end
    getPossibleTiles(P,PlaquetteList,Config) = [i for (i,t) in enumerate(PlaquetteList) if canPlaceTile(P,t) && TileFulFillsConstraint(Config,P,t)]

    function t_delete(j)
        if j <= 5
            return 1
        else
            return j-2
        end
    end
    function constructSpinConfigFromPlaquettes_old(LPx,LPy,PlaquetteList;maxNumTries = 10_000,triesDeleteRow = 10)

        Lx = 2*LPx+1
        Ly = 2*LPy+1
        Mat = fill(NaN,Lx,Ly)
        
        El = PlaquetteList[begin]
        P = SpinConfig(Mat,El.S)
        j = 2
        it = 0
        while j <= Ly-1
            i = 2
    
            
            tries = 0
            while i <= Lx-1 
                
                Pij = getPlaquette(P,i,j)
                it+=1
                SubConf = SubConfig(P,i-3:i+3,j-3:j+3)
                tileNums = getPossibleTiles(Pij,PlaquetteList,SubConf)
                
                if it > maxNumTries
                    println((i,j))
                    @warn "max Iterations reached" 
                    if isempty(tileNums)
                        @warn "no possible tiles"
                    end
                    return P
                end
                
                if isempty(tileNums)
                    tries += 1
    
                    if tries > triesDeleteRow
                        jdelete = t_delete(j)
                        
                        
                        P[:,jdelete:end] .= NaN
                        
                        j = max(jdelete,2)
                        i = 2
                        tries = 0
                    else
    
                        P[:,j:end] .= NaN
                        i = 2
                    end
    
                    continue
                end
                T = PlaquetteList[rand(tileNums)]
                Pij .= T
                i += 2
            end
            j += 2
        end
        @assert fulFillsConstraint(P,verbose=true) "initial configuration does not fulfill constraint"
        return P
    end

    function getIdxInBounds(i,L)
        return (i,j)
    end

    function SubConfig(S::SpinConfig,irange,jrange)
        imin, imax = extrema(irange)

        imin = max(imin,firstindex(S.Mat,1))
        imax = min(imax,lastindex(S.Mat,1))

        jmin, jmax = extrema(jrange)
        jmin = max(jmin,firstindex(S.Mat,2))
        jmax = min(jmax,lastindex(S.Mat,2))

        Matview = @view S.Mat[imin:imax,jmin:jmax]
        return SpinConfig(Matview,S.S)
    end

    getFittingTiles(P,PlaquetteList) = [i for (i,t) in enumerate(PlaquetteList) if canPlaceTile(P,t)]

    function getFittingTiles!(TileList,P,PlaquetteList)
        len = 0
        for (i,t) in enumerate(PlaquetteList)
            if canPlaceTile(P,t)
                len += 1
                TileList[len] = i
            end
        end
        return @view TileList[1:len]
    end

    function getfirstFittingTile(P,shuffleList,indices) 
        for (t,i) in zip(shuffleList,indices)
            if canPlaceTile(P,t)
                return i
            end
        end
        return 0
    end
    
    mapRightToLeft(i,L) = L-i+1

    function constructSpinConfigFromPlaquettes(LPx,LPy,PlaquetteList;maxNumTries = 10_000,triesDeleteRow = 10,deleteRows = 2,rightToLeft = false)

        Lx = 2*LPx+1
        Ly = 2*LPy+1
        Mat = fill(NaN,Lx,Ly)
        
        El = PlaquetteList[begin]
        P = SpinConfig(Mat,El.S)
        j = 2
        it = 0
        i_start(j) = 2 + isodd(j)
    
        while j <= Ly -1
            i = i_start(j)
    
            tries = 0
            while i <= Lx - 1
                i´ = rightToLeft ? mapRightToLeft(i,Lx) : i

                Pij = getPlaquette(P,i´,j)
                it+=1
                tileNums = getFittingTiles(Pij,PlaquetteList)
    
                if it > maxNumTries
                    @warn "max Iterations reached" 
                    if isempty(tileNums)
                        @warn "no possible tiles"
                    end
                    return P
                end
                
                if isempty(tileNums)
                    tries += 1
    
                    if  tries > triesDeleteRow
                        jdel = max(j-deleteRows+1,1)
                        P[:,jdel:end] .= NaN
                        j = max(j-deleteRows,2)
                        tries = 0
                        i = i_start(j)
                    # elseif tries < 3
                    #     # P[:,j+1:end] .= NaN
                    #     idel = max(i-2,i_start(j))
                    #     P[idel:end,j+1:end] .= NaN
                    #     i = idel
                    else
                        P[:,j+1:end] .= NaN
                        i = i_start(j)
                    end
    
                    continue
                end
                T = PlaquetteList[rand(tileNums)]
                Pij .= T
                i += 2
            end
            j += 1
        end

        @assert fulFillsConstraint(P,verbose=true) "initial configuration does not fulfill constraint"
        return P
    end

    function spiralPath(L)
        num_points = 4L^2
        coords = [(0, 0)]
        x, y = 0, 0
        dx, dy = 1, 0
        side_length = 1
        steps_in_side = 0
    
        for i in 2:num_points
            x += dx
            y += dy
            push!(coords, (x, y))
            steps_in_side += 1
    
            if steps_in_side == side_length
                steps_in_side = 0
                if dx == 1 && dy == 0
                    dx, dy = 0, 1
                elseif dx == 0 && dy == 1
                    dx, dy = -1, 0
                    side_length += 1
                elseif dx == -1 && dy == 0
                    dx, dy = 0, -1
                elseif dx == 0 && dy == -1
                    dx, dy = 1, 0
                    side_length += 1
                end
            end
        end
        filter!(x -> iseven(x[1]+x[2]), coords)
        for i in eachindex(coords)
            coords[i] = coords[i] .+ (L,L)
        end
        return coords
    end

    xdirecPath(LPx,LPy=LPx) = [(i,j) for j in 1:2LPy+1 for i in 1:2LPx+1 if iseven(i+j) ]
    xdirecPathReverse(LPx,LPy=LPx) = [(i,j) for j in 1:2LPy+1 for i in 2LPx+1:-1:1 if iseven(i+j) ]
    ydirecPath(LPx,LPy=LPx) = [(i,j) for i in 1:2LPx+1 for j in 1:2LPy+1 if iseven(i+j)]
    ydirecPathReverse(LPx,LPy=LPx) = [(i,j) for i in 1:2LPx+1 for j in 2LPy+1:-1:1 if iseven(i+j)]

    function correctPath!(path,Config)
        filter!(x -> plaquetteIsInBounds(Config,x...),path)
    end

    function correctPath(path,Config)
        filter(x -> plaquetteIsInBounds(Config,x...),path)
    end

    function constructConfigPath(LPx,LPy,PlaquetteList,
        path = xdirecPath(LPx,LPy);
        maxiter= 10000,
        deleteSteps = i->LPx,
        verbose = true,
        )
        
        Lx = 2*LPx+1
        Ly = 2*LPy+1
        Mat = fill(NaN,Lx,Ly)
        
        El = PlaquetteList[begin]
        P = SpinConfig(Mat,El.S)
        filter!(x -> plaquetteIsInBounds(P,x...),path)
        
        indices = shuffle(collect(eachindex(PlaquetteList)))
        shuffleList = PlaquetteList[indices]
        # path = spiralPath(LPx)
        tilingHistory = Int[]
        iter = 1
        TotIter = 1
        function applyStep!(P,n)
            i,j = path[n]
            T = PlaquetteList[tilingHistory[n]]
            Pij = getPlaquette(P,i,j)
            addTile!(Pij,T)
            return P
        end

        function resetFromCheckpoint!(P,iterNum)
            P.Mat .= NaN
            for n in 1:iterNum
                applyStep!(P,n)
            end
        end


        while iter < lastindex(path)-1
            
            iter = lastindex(tilingHistory)
            
            i,j = path[iter+1]
            TotIter += 1
            if TotIter > maxiter 
                if verbose 
                    @warn "maxiter reached"
                end
                break 
            end
            Pij = getPlaquette(P,i,j)
            if TotIter %2 == 0
                shuffle!(indices)
                for (i,t) in enumerate(indices)
                    shuffleList[i] = PlaquetteList[t]
                end
            end
            iT = getfirstFittingTile(Pij,shuffleList,indices)
            if iT == 0
                deleteat!(tilingHistory,max(1,iter-deleteSteps(iter)):iter)
                iter = lastindex(tilingHistory)
                resetFromCheckpoint!(P,iter)
                continue
            end
            T = PlaquetteList[iT]
            addTile!(Pij,T)
            push!(tilingHistory,iT)
            
        end
        # P.Mat .= NaN

        # confs = [copy(applyStep!(P,i)) for i in eachindex(tilingHistory)]
        # @assert fulFillsConstraint(P,verbose=true) "initial configuration does not fulfill constraint"
        # return confs
        return P
    end

    abstract type GroundStateAlgorithm end
    
    struct DictAlgorithm <: GroundStateAlgorithm end

    function getFittingTilesDict!(fittingTilesDict,P,PlaquetteList)::Vector{Int}
        tilenums = get(fittingTilesDict,P,nothing)
        if tilenums !== nothing
            return tilenums
        else
            newtiles = getFittingTiles(P,PlaquetteList)
            fittingTilesDict[P] = newtiles
            return newtiles
        end
    end
    using CairoMakie
    
    function setupCalc!(path,LPx,LPy,PlaquetteList)
        Lx = 2*LPx+1
        Ly = 2*LPy+1
        Mat = fill(NaN,Lx,Ly)
        El = PlaquetteList[begin]
        P = SpinConfig(Mat,El.S)
        filter!(x -> plaquetteIsInBounds(Mat,x...),path)

        emptyTilesList = getFreeTilesPath(P,path,PlaquetteList)

        # TilesDict = getFittingTilesDict(PlaquetteList)
        (;path,emptyTilesList)
    end

    # function getFittingTilesDict(PlaquetteList)

    #     P_init = (NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN)
        
    #     vals = unique([i  for P in PlaquetteList for i in P])

    #     fittingTilesDict = Dict(P => getFittingTiles(P,PlaquetteList) for P in Iterators.product(vals,vals,vals,vals,vals,vals,vals,vals))

    #     return fittingTilesDict
    # end
    function getStepDeleter(L,default = 5,tries = 5)
        function deleter(x,counter)
            counter[] += 1
            # println(counter[])
            if counter[] > tries
                counter[] = 0
                return L
            end
            return default
        end
    end
    
    function constructConfigPath(::DictAlgorithm, P::SpinConfig,PlaquetteList,setup;
        maxiter= 10000,
        deleteSteps = i->LPx,
        verbose = true,
        plotSteps = false,
        )
        
        path = setup.path
        emptyTilesList = setup.emptyTilesList
        
        # path = spiralPath(LPx)
        tilingHistory = Int[]
        sizehint!(tilingHistory,length(path))
        iter = 1
        TotIter = 1

        function applyStep!(P,n)
            i,j = path[n]
            T = PlaquetteList[tilingHistory[n]]
            Pij = getPlaquette(P,i,j)
            addTile!(Pij,T)
            return P
        end

        function resetToStep!(P,iterNum)
            emptyTiles = emptyTilesList[iterNum+1]
            for i in eachindex(emptyTiles)
                I = emptyTiles[i]
                P[I] = NaN
            end
        end

        P_init = getSitesFromPlaquette(getPlaquette(P, path[lastindex(tilingHistory)+1]...))
        fittingTilesDict = Dict(P_init => getFittingTiles(P_init,PlaquetteList))

        counter = Ref(0)
        successfulplacements = 0
        
        while iter < lastindex(path)-1
            
            iter = lastindex(tilingHistory)
            
            i,j = path[iter+1]
            TotIter += 1
            if TotIter > maxiter 
                if verbose 
                    @warn "maxiter reached"
                end
                break 
            end
            Pij = getPlaquette(P,i,j)
            PijSites = getSitesFromPlaquette(Pij)

            tileList = getFittingTilesDict!(fittingTilesDict,PijSites,PlaquetteList)
            if isempty(tileList)
                deleteat!(tilingHistory,max(2,iter-deleteSteps(iter,counter)):iter)
                iter = lastindex(tilingHistory)
                resetToStep!(P,iter)
                successfulplacements = 0
                continue
            end
            iT = rand(tileList)
            T = PlaquetteList[iT]
            addTile!(Pij,T)
            push!(tilingHistory,iT)
            if successfulplacements > 1
                counter[] = 0
            end
            successfulplacements +=1 
            plotSteps && display(plotSpinConfig(P))

        end
        return P
    end
        
    function constructConfigPath(Algo::DictAlgorithm,LPx,LPy,PlaquetteList,
        setup = setupCalc!(xdirecPath(LPx,LPy),LPx,LPy,PlaquetteList);
        kwargs...
        )
        
        Lx = 2*LPx+1
        Ly = 2*LPy+1
        Mat = fill(NaN,Lx,Ly)
        
        El = PlaquetteList[begin]
        P = SpinConfig(Mat,El.S)
        return constructConfigPath(Algo,P,PlaquetteList,setup; kwargs...)
    end

    using CairoMakie.Makie.ColorSchemes
    function plotSpinConfig!(ax,S::SpinConfig;plotConstraints=true,kwargs...)
        vals = filter!(x-> !ismissing(x) && !isnan(x),unique(S.Mat))
        isempty(vals) && (vals = [-1,1])
        Amin = min(minimum(vals),-S.S)
        Amax = max(maximum(vals),S.S)
        us = Amin:Amax
        # hm = heatmap!(ax,Array(S.Mat),colorrange = (Amin,Amax),colormap = cgrad(:linear_bgy_10_95_c74_n256, length(us), categorical = true);kwargs...)
        hm = heatmap!(ax,Array(S.Mat),colorrange = (Amin,Amax),colormap = cgrad(:grays, length(us), categorical = true);kwargs...)
        if plotConstraints
            points = [Point(Tuple(I)...) for I in CartesianIndices(S.Mat) if iseven(sum(Tuple(I)))]
            scatter!(ax,points,marker = '×',  color = :gray, markersize = 20)
        end
        translate!(hm, 0, 0, -100)
        return hm
    end
    
    function getConfigAxis(S;kwargs...)
        (;aspect = DataAspect(),
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
        xminorticks = 0.5 .+ axes(S,1),
        yminorticks = 0.5 .+ axes(S,2),
        limits = (0.5,size(S,1)+0.5,0.5,size(S,2)+0.5),
        )

    end

    function plotSpinConfig(S;kwargs...) 
        fig = Figure()
        ax = Axis(fig[1,1];
            getConfigAxis(S)...
        )
        hm = plotSpinConfig!(ax,S;kwargs...)
        us = filter!(x -> !ismissing(x) && !isnan(x) ,unique(S))
        isempty(us) && (us = [-1,1])
        Colorbar(fig[1,2],hm,ticks = us)
        return fig
    end

    function getFreeTilesPath(P,path,PlaquetteList)
        # emptyTiles = Vector{Int}[]
        emptyTiles = Vector{CartesianIndex{2}}[]
        newP = copy(P)
        newP .= NaN
        testPlaq = first(PlaquetteList)
        for (i,j) in path
            Pij = getPlaquette(newP,i,j)
            addTile!(Pij,testPlaq)
            emptyTilesCurrent = findall(isnan,newP)
            push!(emptyTiles,emptyTilesCurrent)
        end
        return emptyTiles
    end
    include("PeriodicTilings.jl")
    include("StructureFactor.jl")
    include("Fluctuations.jl")

end # module SpiderWebModel
