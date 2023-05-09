module SpiderWebModel
    using StaticArrays,SpinFRGLattices


    function SpiderWebBasis()
        a1 = SA[1,-1]
        a2 = SA[1,1]
        b = [SA[0,0],SA[1,0]]

        return Basis_Struct_2D(;a1,a2,b,NNdist = 1,NUnique = 2,SiteType= [1,2],refSites= [Rvec(0,0,1),Rvec(0,0,2)])
    end

    const BasisSW = SpiderWebBasis()
    getCartesian(x) = SpinFRGLattices.getCartesian(x,BasisSW)
    
    function constraint(S1,S2,S3,S4,S5,S6,S7,S8)
        S1+S2-S3-S4+S5+S6-S7-S8
    end

    constraint(S::AbstractArray) = constraint(S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[8])

    const ALLGS_S12 = getAllGS(0.5)

    function generateLattice(L)
        sites = [Rvec(i,j,nb) for i in 0:L-1 for j in 0:L-1 for nb in eachindex(BasisSW.b)]
    end

    struct SpinConfig2{RV<:Rvec,T<:Real}
        Spins::Dict{RV,T}
    end

    SpinConfig(::Type{<:R}) where R <: Rvec = SpinConfig2(Dict{R,Int}())
    
    import Base:getindex,setindex!

    Base.getindex(S::SpinConfig2,i) = getindex(S.Spins,i)
    Base.setindex!(S::SpinConfig2,x,i) = setindex!(S.Spins,x,i)
    Base.setindex!(S::SpinConfig2,x::Missing,i) = nothing

    function isInPlaquetteOf(R::Rvec_2D,Rp::Rvec_2D)
        d = dist(R,Rp,BasisSW)
        return d <= √(2)
    end

    struct Plaquette2{T}
        x::SMatrix{3, 3, Union{Missing, T}, 9}
    end
    """Order in which the spins are stored in the plaquette corresponding to the numbering convention"""
    plaquetteOrder(S1,S2,S3,S4,S5,S6,S7,S8) = SMatrix{3,3}(S4,S5,S6,S3,missing,S7,S2,S1,S8)
    plaquetteOrder(x::AbstractArray) = plaquetteOrder(x...)

    function reversePlaquetteOrder(x) 
        S4,S5,S6,S3,_,S7,S2,S1,S8 = x
        return SA[S1,S2,S3,S4,S5,S6,S7,S8]
    end

    Plaquette2(S1,S2,S3,S4,S5,S6,S7,S8) = Plaquette2(plaquetteOrder(S1,S2,S3,S4,S5,S6,S7,S8))
    Plaquette2(x::AbstractArray) = Plaquette2(x...)

    function constraint(P::Plaquette2)
        return constraint(reversePlaquetteOrder(P.x))
    end


    """generates all possible combinations 0,1 of the 8 spins in a plaquette """

    function getAllGS(S)
        allCombs = [Plaquette2(SA[i,j,k,l,m,n,o,p]) for i in -S:S for j in -S:S for k in -S:S for l in -S:S for m in -S:S for n in -S:S for o in -S:S for p in -S:S if constraint(i,j,k,l,m,n,o,p) == 0]
        return allCombs
    end


    function getPlaquette(R::Rvec_2D)
        R.b ≠ 1 && error("R must be on sublattice A")
        (;n1,n2) = R
        S2,S4,S6,S8 = Rvec(n1,n2+1,1),Rvec(n1-1,n2,1),Rvec(n1,n2-1,1),Rvec(n1+1,n2,1) #red sites
        S1,S3,S5,S7 = Rvec(n1,n2,2),Rvec(n1-1,n2,2),Rvec(n1-1,n2-1,2),Rvec(n1,n2-1,2) #black sites

        return Plaquette2(S1,S2,S3,S4,S5,S6,S7,S8)
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

    function setTile!(S::SpinConfig2,R,T1)
        P = getPlaquette(R)
        TO = plaquetteOrder(T1)
        for (r,t) in zip(P,TO)
            S[r] = t
        end
    end

end # module SpiderWebModel
