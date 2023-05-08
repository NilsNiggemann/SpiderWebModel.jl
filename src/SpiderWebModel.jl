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
    
    function constr(S1,S2,S3,S4,S5,S6,S7,S8)
        S1+S2-S3-S4+S5+S6-S7-S8
    end

    constr(S::AbstractArray) = constr(S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[8])

    """generates all possible combinations 0,1 of the 8 spins in a plaquette """

    function getAllGS(S)
        allCombs = [SA[i,j,k,l,m,n,o,p] for i in -S:S for j in -S:S for k in -S:S for l in -S:S for m in -S:S for n in -S:S for o in -S:S for p in -S:S if constr(i,j,k,l,m,n,o,p) == 0]
        return allCombs
    end

    const ALLGS_S12 = getAllGS(0.5)

    function generateLattice(L)
        sites = [Rvec(i,j,nb) for i in 0:L-1 for j in 0:L-1 for nb in eachindex(BasisSW.b)]
    end

    struct Spin{RV<:Rvec}
        m::Int
        R::RV
    end
    getSpinValue(S::Spin) = S.m
    struct SpinConfig1{RV<:Rvec}
        Spins::Dict{Rvec_2D,Spin{RV}}
    end

    SpinConfig1() = SpinConfig1(Dict{Rvec_2D,Int}())

    import Base.getproperty

    Base.getproperty(S::SpinConfig1,s::Symbol) = getproperty(S.Spins,s)

    function isInPlaquetteOf(R::Rvec_2D,Rp::Rvec_2D)
        d = dist(R,Rp,BasisSW)
        return d <= √(2)
    end

    function getPlaquette(R::Rvec_2D)
        R.b ≠ 1 && error("R must be on sublattice A")
        (;n1,n2) = R
        S2,S4,S6,S8 = Rvec(n1,n2+1,1),Rvec(n1-1,n2,1),Rvec(n1,n2-1,1),Rvec(n1+1,n2,1) #red sites
        S1,S3,S5,S7 = Rvec(n1,n2,2),Rvec(n1-1,n2,2),Rvec(n1-1,n2-1,2),Rvec(n1,n2-1,2) #black sites
        return (S1,S2,S3,S4,S5,S6,S7,S8)
    end

    function fulfillsConstraint(R::Rvec_2D,S::SpinConfig1)
        S1,S2,S3,S4,S5,S6,S7,S8 = getPlaquette(R)
        constr(S(S1),S(S2),S(S3),S(S4),S(S5),S(S6),S(S7),S(S8)) == 0
    end

end # module SpiderWebModel
