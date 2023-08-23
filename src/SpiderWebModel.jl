module SpiderWebModel
    using StaticArrays,Random,Statistics

    import Base:size,getindex,setindex!,iterate,show,copy

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


    function getCorrel(Conf::SpinConfig)
        S_ij = [si*sj for si in Conf.Mat, sj in Conf.Mat]
        return S_ij
    end

    plaquetteOperator() = SpinConfig(
        float(SA[
            +1 +1 -1;
            -1  0 -1;
            -1 +1 +1
        ]),
        1/2
    )

    function PlaquetteOperator!(Conf::SpinConfig, i::Integer,j::Integer)
        P = getPlaquette(Conf,i,j)
        Op = plaquetteOperator()
        P .+= Op
        return Conf
    end

    function plaquetteIsInBounds(Conf::SpinConfig,iCenter::Integer,jCenter::Integer)
        for i in iCenter-1:iCenter+1, j in jCenter-1:jCenter+1
            if !checkbounds(Bool,Conf.Mat,i,j)
                return false
            end
        end
        return true
    end

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
    function fulFillsConstraint_nonStrict(Conf::SpinConfig;verbose = false)

        for i in axes(Conf.Mat,1), j in axes(Conf.Mat,2)
            iseven(i+j) || continue 
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

    function fulFillsConstraint(Conf::SpinConfig;verbose = false)
        allSpinsInBounds(Conf;verbose) || return false
        return fulFillsConstraint_nonStrict(Conf;verbose)
    end

    function CanApplyPlaquette(Conf::SpinConfig,i::Integer,j::Integer)
        plaquetteIsInBounds(Conf,i,j) || return false
        P = getPlaquette(Conf,i,j)
        □ = plaquetteOperator()
        
        P .+= □

        applicable = fulFillsConstraint(Conf)

        P .-= □

        return applicable
    end

    function PlaquetteOperatorSave!(Conf::SpinConfig, i::Integer,j::Integer)
        CanApplyPlaquette(Conf,i,j) || return Conf
        PlaquetteOperator!(Conf,i,j)
    end

    function flipSpinsAlongLine!(Conf,org,slope)
        slope ∈ (-Inf, Inf) && return flipSpinsAlongRow!(Conf,org[2])
        for i in axes(Conf.Mat,1), j in axes(Conf.Mat,2)
            if slope*(i-org[1]) == j-org[2]
                Conf[i,j] *= -1
            end
        end
        return Conf
    end

    function flipSpinsAlongRow!(Conf,i)
        Conf[i,:] .*= -1
        return Conf
    end

    function CanApply(Conf::SpinConfig,Op::SpinConfig,i,j)
        plaquetteIsInBounds(Conf,i,j) || return false
        P = getPlaquette(Conf,i,j)
        
        P .+= Op

        applicable = fulFillsConstraint(Conf)

        P .-= Op

        return applicable
    end

    function CanApplyAnywhere(Conf::SpinConfig,Op::SpinConfig)
        a1 = axes(Conf.Mat,1)
        a2 = axes(Conf.Mat,2)

        Opx,Opy = size(Op)
        for i in a1, j in a2
            firstindex(a1)+Opx <= i <= lastindex(a1)-Opx || continue
            firstindex(a2)+Opy <= j <= lastindex(a2)-Opy || continue

            if CanApply(Conf,Op,i,j)
                return true
            end
        end
        return false
    end
    
    function getApplicablePlaquettes(Conf::SpinConfig,Op::SpinConfig)
        plaqPos = [(i,j) for i in axes(Conf.Mat,1) for j in axes(Conf.Mat,2) if CanApply(Conf,Op,i,j)]
        return plaqPos
    end
    
    function canTileUpRight(T1,T2)
        T1´ = T1[1:2,2:3]
        T2´ = T2[2:3,1:2]
        return all(_isZeroOrNaN,T1´-T2´)
    end
    
    function canTileUpLeft(T1,T2)
        T1´ = T1[1:2,1:2]
        T2´ = T2[2:3,2:3]
        return all(_isZeroOrNaN,T1´-T2´)
    end
    
    function canTileDownRight(T1,T2)
        T1´ = T1[2:3,2:3]
        T2´ = T2[1:2,1:2]
        return all(_isZeroOrNaN,T1´-T2´)
    end

    function canTileDownLeft(T1,T2)
        T1´ = T1[2:3,1:2]
        T2´ = T2[1:2,2:3]
        return all(_isZeroOrNaN,T1´-T2´)
    end
    
    function getTilings(T1,AllTilings)
        UpRight = [i for (i,T2) in enumerate(AllTilings) if canTileUpRight(T1,T2)]
        UpLeft = [i for (i,T2) in enumerate(AllTilings) if canTileUpLeft(T1,T2)]
        DownRight = [i for (i,T2) in enumerate(AllTilings) if canTileDownRight(T1,T2)]
        DownLeft = [i for (i,T2) in enumerate(AllTilings) if canTileDownLeft(T1,T2)]
        return (;UpRight,UpLeft,DownRight,DownLeft)
    end

    function canPlaceTile(P::SpinConfig,T1)
        return all(_isZeroOrNaN, p -t for (p,t) in zip(P,T1))
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

    function constructSpinConfigFromPlaquettes(LPx,LPy,PlaquetteList;maxNumTries = 10_000,triesDeleteRow = 10,deleteRows = 2)

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
    
                Pij = getPlaquette(P,i,j)
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


    function fillMissing!(P::SpinConfig,PlaquetteList)
    
        OldPij = fill(NaN,3,3)
        Lx,Ly = size(P.Mat)
        for i in 3:2:Lx-2, j in 3:2:Ly-2

            Pij = getPlaquette(P,i,j)
            OldPij .= Pij
            @assert fulFillsConstraint(P,verbose = true) "constraint not fulfilled in fillMissing"
            
            for T in PlaquetteList
                
                if canPlaceTile(Pij,T)
                    Pij .= T
                    if fulFillsConstraint(P)
                        break
                    else 
                        Pij .= OldPij
                    end
                end
            end
            @assert fulFillsConstraint(P,verbose = true) "constraint not fulfilled in fillMissing2"
        end
        
        return P
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


end # module SpiderWebModel
