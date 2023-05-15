module SpiderWebModel
    using StaticArrays

    import Base:size,getindex,setindex!,iterate,show

    struct SpinConfig{T,MatType <: AbstractMatrix{T}} <:AbstractMatrix{T}
        Mat::MatType
        S::T
    end

    Base.getindex(S::SpinConfig,i,j) = getindex(S.Mat,i,j)
    Base.setindex!(S::SpinConfig,x,i,j) = setindex!(S.Mat,x,i,j)
    Base.setindex!(::SpinConfig,::Missing,i,j) = nothing
    Base.iterate(S::SpinConfig,i) = iterate(S.Mat,i)
    Base.iterate(S::SpinConfig) = iterate(S.Mat)

    Base.size(S::SpinConfig) = size(S.Mat)

    Base.show(io::IO, ::MIME"text/plain", S::SpinConfig) = print(io,S.Mat)
    Base.show(io::IO, S::SpinConfig) = print(io,S.Mat)

    function constraint(S1,S2,S3,S4,S5,S6,S7,S8)
        S1+S2-S3-S4+S5+S6-S7-S8
    end

    constraint(S) = constraint(S...)

    """Order in which the spins are stored in the plaquette corresponding to the numbering convention"""
    plaquetteOrder(S1,S2,S3,S4,S5,S6,S7,S8) = SMatrix{3,3}(S4,S5,S6,S3,missing,S7,S2,S1,S8)
    plaquetteOrder(x) = plaquetteOrder(x...)

    function plaquette(sites::NTuple{8,T},S=1/2) where T
        Mat = SMatrix{3,3}(plaquetteOrder(sites))
        return SpinConfig(Mat,S)
    end

    plaquette(S1,S2,S3,S4,S5,S6,S7,S8) = plaquette((S1,S2,S3,S4,S5,S6,S7,S8))
    plaquette(x::AbstractVector) = plaquette(x...)

    function getSitesFromPlaquette(Mat::AbstractMatrix)
        S4,S5,S6,S3,_,S7,S2,S1,S8 = Mat
        return (S1,S2,S3,S4,S5,S6,S7,S8)
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
        allCombs = [plaquette(SA[i,j,k,l,m,n,o,p]) for i in -S:S for j in -S:S for k in -S:S for l in -S:S for m in -S:S for n in -S:S for o in -S:S for p in -S:S if constraint(i,j,k,l,m,n,o,p) == 0]
        return allCombs
    end

    const ALLGS_S12 = getAllGS(0.5)

    _isZeroOrMissing(x) = x === missing || x == 0

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


    function setTile!(S::SpinConfig,i::Integer,j::Integer,T1)
        P = getPlaquette(S,i,j)
        P .= T1
    end

    function getCorrel(Conf::SpinConfig)
        S_ij = [si*sj for si in Conf.Mat, sj in Conf.Mat]
        return S_ij
    end

    plaquetteOperator() = SpinConfig(
        float(SA[
            -1 -1 +1;
            +1  0 +1;
            +1 -1 -1
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
        if any(x->abs(x)>Conf.S,Conf.Mat) 
            verbose && println("Spin larger than S")
            return false
        end
        
        for i in axes(Conf.Mat,1), j in axes(Conf.Mat,2)
            iseven(i+j) || continue 
            plaquetteIsInBounds(Conf,i,j) || continue
            P = getPlaquette(Conf,i,j)

            if constraint(P) != 0
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

end # module SpiderWebModel
