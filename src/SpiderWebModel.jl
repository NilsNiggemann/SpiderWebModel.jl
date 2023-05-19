module SpiderWebModel
    using StaticArrays,Random

    import Base:size,getindex,setindex!,iterate,show,copy

    struct SpinConfig{T,MatType <: AbstractMatrix{T},T1<:Real} <:AbstractMatrix{T}
        Mat::MatType
        S::T1
    end

    Base.getindex(S::SpinConfig,i,j) = getindex(S.Mat,i,j)
    Base.setindex!(S::SpinConfig,x,i,j) = setindex!(S.Mat,x,i,j)
    Base.setindex!(::SpinConfig,::Missing,i,j) = nothing
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
    plaquetteOrder(S1,S2,S3,S4,S5,S6,S7,S8) = SMatrix{3,3}(S4,S5,S6,S3,missing,S7,S2,S1,S8)
    plaquetteOrder(x) = plaquetteOrder(x...)

    function plaquette(sites::NTuple{8,T},S=1/2) where T
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
        allCombs = [plaquette((i,j,k,l,m,n,o,p)) for i in -S:S for j in -S:S for k in -S:S for l in -S:S for m in -S:S for n in -S:S for o in -S:S for p in -S:S if constraint((i,j,k,l,m,n,o,p)) == 0]
        return allCombs
    end

    const ALLGS_S12 = getAllGS(0.5)

    _isZeroOrMissing(x) = x === missing || x == 0


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

    """Assumes that constraint are only defined on every alternating site, starting from the first index"""
    function fulFillsConstraint(Conf::SpinConfig;verbose = false)
        for x in Conf.Mat
            ismissing(x) && continue
            if abs(x) > Conf.S 
                verbose && println("Spin larger than S")
                return false
            end
        end

        for i in axes(Conf.Mat,1), j in axes(Conf.Mat,2)
            iseven(i+j) || continue 
            plaquetteIsInBounds(Conf,i,j) || continue
            P = getPlaquette(Conf,i,j)
            c = constraint(P)
            ismissing(c) && continue

            if c ≠ 0
                verbose && println("Constraint not fulfilled at i,j = $((i,j))" )
                return false
            end
        end
        return true
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

        for i in a1, j in a2
            firstindex(a1)+2 <= i <= lastindex(a1)-2 || continue
            firstindex(a2)+2 <= j <= lastindex(a2)-2 || continue

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
    
    # function sumNonMissing(sites)
    #     cons = constraintSigns()
    #     return sum(sgn*S for (sgn,S) in zip(cons,sites) if !ismissing(S))
    # end

    # function fillMissingConstraint!(sites)
    #     numMiss = count(ismissing,Sites)
    #     g
    #     return Conf
    # end

    
    function canTileUpRight(T1,T2)
        T1´ = T1[1:2,2:3]
        T2´ = T2[2:3,1:2]
        return all(_isZeroOrMissing,T1´-T2´)
    end
    
    function canTileUpLeft(T1,T2)
        T1´ = T1[1:2,1:2]
        T2´ = T2[2:3,2:3]
        return all(_isZeroOrMissing,T1´-T2´)
    end
    
    function canTileDownRight(T1,T2)
        T1´ = T1[2:3,2:3]
        T2´ = T2[1:2,1:2]
        return all(_isZeroOrMissing,T1´-T2´)
    end

    function canTileDownLeft(T1,T2)
        T1´ = T1[2:3,1:2]
        T2´ = T2[1:2,2:3]
        return all(_isZeroOrMissing,T1´-T2´)
    end
    
    function getTilings(T1,AllTilings)
        UpRight = [i for (i,T2) in enumerate(AllTilings) if canTileUpRight(T1,T2)]
        UpLeft = [i for (i,T2) in enumerate(AllTilings) if canTileUpLeft(T1,T2)]
        DownRight = [i for (i,T2) in enumerate(AllTilings) if canTileDownRight(T1,T2)]
        DownLeft = [i for (i,T2) in enumerate(AllTilings) if canTileDownLeft(T1,T2)]
        return (;UpRight,UpLeft,DownRight,DownLeft)
    end

    function canPlaceTile(P::SpinConfig,T1)
        return all(_isZeroOrMissing, p -t for (p,t) in zip(P,T1))
    end

    function setTile!(S::SpinConfig,i::Integer,j::Integer,T1)
        P = getPlaquette(S,i,j)
        P .= T1
    end

    getPossibleTiles(P,PlaquetteList) = [i for (i,t) in enumerate(PlaquetteList) if canPlaceTile(P,t)]
            
    function getRanTileNum!(randBuffer,tiles) 
        isempty(randBuffer) && error("too many steps")
        maxIndex = length(tiles)

        ranNumber = pop!(randBuffer)
        idx = ( (ranNumber - 1) % maxIndex) +1
        if idx == 0 || idx > maxIndex
            error("idx == $idx for ranNumber = $ranNumber and maxIndex = $maxIndex")
        end

        return tiles[idx]
    end

    function constructSpinConfigFromPlaquettes(LPx,LPy,PlaquetteList;randBuffSize=1_000_000,maxNumTries = 10_000,deleteRows = 4, randBuffer = rand(1:70,randBuffSize)
        )
        Lx = 2*LPx+1
        Ly = 2*LPy+1
    
        Mat = Matrix{Union{Missing,Float64}}(undef,Lx,Ly)
        fill!(Mat,missing)
    
        El = PlaquetteList[begin]
        P = SpinConfig(Mat,El.S)
    
        j = 2
        while j <= Ly-1
            i = 2
            while i <= Lx-1 
            
                Pij = getPlaquette(P,i,j)
                
                tileNums = getPossibleTiles(Pij,PlaquetteList)
                
                if isempty(tileNums)
                    jmin = max(j-deleteRows,2)
                    
                    Mat[:,jmin:end] .= missing
                    j = jmin
                    i = 2
                    continue
                end
                T = PlaquetteList[getRanTileNum!(randBuffer,tileNums)]
                Pij .= T
    
                i += 2
            end
            j += 2
    
        end
        @assert fulFillsConstraint(P,verbose=true) "initial configuration does not fulfill constraint"

        OldState = copy(P)
        for i in 1:randBuffSize
            i % 100 == 0 && println("i = $i")
            randPlaqlist = shuffle(PlaquetteList)
            fillMissing!(P,randPlaqlist)
            if !any(ismissing,P) 
                @assert fulFillsConstraint(P,verbose = true) "constraint not fulfilled"
                return P
            end
            P.Mat .= OldState.Mat
        end

        fulFillsConstraint(P,verbose = true) || warn("constraint not fulfilled")
        return P
    end

    function fillMissing!(P::SpinConfig,PlaquetteList)
    
        OldPij = Matrix{Union{Missing,Float64}}(undef,3,3)
        fill!(OldPij,missing)
        Lx,Ly = size(P.Mat)
        for i in 3:2:Lx-2, j in 3:2:Ly-2

            Pij = getPlaquette(P,i,j)
            OldPij .= Pij
            @assert fulFillsConstraint(P,verbose = true) "constraint not fulfilled in fillMissing"
            
            for T in PlaquetteList
                
                if canPlaceTile(Pij,T)
                    Pij.Mat .= T
                    if fulFillsConstraint(P)
                        break
                    else 
                        Pij.Mat .= OldPij
                    end
                end
            end
            @assert fulFillsConstraint(P,verbose = true) "constraint not fulfilled in fillMissing2"
        end
        
        return P
    end

end # module SpiderWebModel
