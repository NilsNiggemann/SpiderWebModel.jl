module SpiderWebModel
    using StaticArrays,SpinFRGLattices

    import Base:size,getindex,setindex!,iterate,show


    function SpiderWebBasis()
        a1 = SA[1,-1]
        a2 = SA[1,1]
        b = [SA[0,0],SA[1,0]]

        return Basis_Struct_2D(;a1,a2,b,NNdist = 1,NUnique = 2,SiteType= [1,2],refSites= [Rvec(0,0,1),Rvec(0,0,2)])
    end

    const Basis = SpiderWebBasis()
    getCartesian(x) = SpinFRGLattices.getCartesian(x,Basis)
    
    function constraint(S1,S2,S3,S4,S5,S6,S7,S8)
        S1+S2-S3-S4+S5+S6-S7-S8
    end

    constraint(S) = constraint(S...)

    function isInPlaquetteOf(R::Rvec_2D,Rp::Rvec_2D)
        d = dist(R,Rp,Basis)
        return d <= √(2)
    end

    """Order in which the spins are stored in the plaquette corresponding to the numbering convention"""
    plaquetteOrder(S1,S2,S3,S4,S5,S6,S7,S8) = SMatrix{3,3}(S4,S5,S6,S3,missing,S7,S2,S1,S8)
    plaquetteOrder(x) = plaquetteOrder(x...)


    struct Plaquette{T,Mat <: AbstractMatrix{T}} <: AbstractMatrix{T}
        Mat::Mat
    end
    
    Base.show(io::IO, ::MIME"text/plain", x::Plaquette) = print(io,x.Mat)
    Base.show(io::IO, x::Plaquette) = print(io,x.Mat)
    
    function plaquette(sites::NTuple{8,T}) where T
        Mat = SMatrix{3,3}(plaquetteOrder(sites))
        return Plaquette(Mat)
    end

    plaquette(S1,S2,S3,S4,S5,S6,S7,S8) = plaquette((S1,S2,S3,S4,S5,S6,S7,S8))
    plaquette(x::AbstractVector) = plaquette(x...)

    function getSitesFromPlaquette(Mat::SMatrix{3,3,T,9}) where T
        S4,S5,S6,S3,_,S7,S2,S1,S8 = Mat
        return (S1,S2,S3,S4,S5,S6,S7,S8)
    end

    function getSitesFromPlaquette(P::Plaquette) 
        return getSitesFromPlaquette(P.Mat)
    end

    Base.size(::Plaquette) = (8,)
    Base.getindex(P::Plaquette,i) = getindex(P.sites,i)
    Base.getindex(P::Plaquette,i,j) = getindex(P.Mat,i,j)
    Base.iterate(P::Plaquette,i) = iterate(P.sites,i)
    Base.iterate(P::Plaquette) = iterate(P.sites)

    function constraint(P::Plaquette)
        return constraint(getSitesFromPlaquette(P))
    end


    """generates all possible combinations 0,1 of the 8 spins in a plaquette """

    function getAllGS(S)
        allCombs = [plaquette(SA[i,j,k,l,m,n,o,p]) for i in -S:S for j in -S:S for k in -S:S for l in -S:S for m in -S:S for n in -S:S for o in -S:S for p in -S:S if constraint(i,j,k,l,m,n,o,p) == 0]
        return allCombs
    end

    const ALLGS_S12 = getAllGS(0.5)

    function getPlaquette(R::Rvec_2D)
        r = getCartesian(R)
        RV(r) = SpinFRGLattices.getRvec(r,Basis)
        S1 = r + SA[1,0] |> RV
        S2 = r + SA[1,1] |> RV
        S3 = r + SA[0,1] |> RV
        S4 = r + SA[-1,1] |> RV
        S5 = r + SA[-1,0] |> RV
        S6 = r + SA[-1,-1] |> RV
        S7 = r + SA[0,-1] |> RV
        S8 = r + SA[1,-1] |> RV
        return plaquette(S1,S2,S3,S4,S5,S6,S7,S8)
    end

    function getPlaquette(i,j)
        S1 = i+1,j |> CartesianIndex
        S2 = i+1,j+1 |> CartesianIndex
        S3 = i,j+1 |> CartesianIndex
        S4 = i-1,j+1 |> CartesianIndex
        S5 = i-1,j |> CartesianIndex
        S6 = i-1,j-1 |> CartesianIndex
        S7 = i,j-1 |> CartesianIndex
        S8 = i+1,j-1 |> CartesianIndex
        return plaquette(S1,S2,S3,S4,S5,S6,S7,S8)
    end


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

    struct SpinConfig{T,MatType <: AbstractMatrix{T}} <:AbstractMatrix{T}
        Mat::MatType
    end

    Base.getindex(S::SpinConfig,i,j) = getindex(S.Spins,i,j)
    Base.setindex!(S::SpinConfig,x,i,j) = setindex!(S.Spins,x,i,j)
    Base.setindex!(::SpinConfig,::Missing,i,j) = nothing
    Base.iterate(S::SpinConfig,i) = iterate(S.Spins,i)
    Base.iterate(S::SpinConfig) = iterate(S.Spins)

    Base.show(io::IO, ::MIME"text/plain", S::SpinConfig) = print(io,S.Spins)
    Base.show(io::IO, S::SpinConfig) = print(io,S.Spins)

    function getPlaquette(S::AbstractMatrix,i,j)
        Mat = @view S[i-1:i+1,j-1:j+1]
        return plaquette(Mat)
    end

    function setTile!(S::SpinConfig,i::Integer,j::Integer,T1)
        P = getPlaquette(i,j)
        for (r,t) in zip(P,T1)
            S[r] = t
        end
    end

    function getCorrel(Conf::SpinConfig)
        S_ij = [si*sj for si in Conf.Spins, sj in Conf.Spins]
        return S_ij
    end

    plaquetteOperatorSigns() = plaquette(+1,+1,-1,-1,+1,+1,-1,-1)
    function PlaquetteOperator!(Conf::SpinConfig, i::Integer,j::Integer)
        P = getPlaquette(Conf,i,j)
        signs = plaquetteOperatorSigns()
        P .+= signs
        return Conf
    end

    function plaquetteIsInBounds(Conf::SpinConfig,iCenter::Integer,jCenter::Integer)
        for i in iCenter-1:iCenter+1, j in jCenter-1:jCenter+1
            if !checkbounds(Bool,Conf.Spins,i,j)
                return false
            end
        end
        return true
    end

    """Assumes that constraint are only defined on every alternating site, starting from the first index"""
    function fulFillsConstraint(Conf::SpinConfig;verbose = false)
        for i in axes(Conf.Spins,1), j in axes(Conf.Spins,2)
            iseven(i+j) && continue 
            plaquetteIsInBounds(Conf,i,j) || continue
            P = getPlaquette(i,j)

            if constraint(P) != 0
                verbose && println("Constraint not fulfilled at $R = $(getCartesian(R))" )
                return false
            end
        end
        return true
    end

    function CanApplyPlaquette(Conf::SpinConfig,i::Integer,j::Integer)
        plaquetteIsInBounds(Conf,i,j) || return false
        P = getPlaquette(Conf,i,j)
        operations = plaquetteOperatorSigns()
        OldP = getSitesFromPlaquette(P)
        
        Pnew = OldP .+ operations
        
        if any(x->abs(x)>Conf.S, Pnew)
            return false
        end
        
        P .= Pnew

        applicable = fulFillsConstraint(Conf)

        P .= OldP

        return applicable
    end

    function PlaquetteOperatorSave!(Conf::SpinConfig, i::Integer,j::Integer)
        CanApplyPlaquette(Conf,i,j) || return Conf
        PlaquetteOperator!(Conf,i,j)
    end

    function 

end # module SpiderWebModel
