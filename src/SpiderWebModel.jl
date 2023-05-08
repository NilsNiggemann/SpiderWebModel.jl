module SpiderWebModel
    using StaticArrays,SpinFRGLattices

    function constr(S1,S2,S3,S4,S5,S6,S7,S8)
        S1+S2-S3-S4+S5+S6-S7-S8
    end

    constr(S::AbstractArray) = constr(S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[8])

    """generates all possible combinations 0,1 of the 8 spins in a plaquette """

    function getAllGS(S)
        allCombs = [SA[i,j,k,l,m,n,o,p] for i in -S:S for j in -S:S for k in -S:S for l in -S:S for m in -S:S for n in -S:S for o in -S:S for p in -S:S if constr(i,j,k,l,m,n,o,p) == 0]
        return allCombs
    end

    getAllGS(1/2)
    ##

    function SpiderWebBasis()
        a1 = SA[1,1]
        a2 = SA[-1,1]
        b = [SA[0,0],SA[0.5,0.5]]

        return Basis_Struct_2D(;a1,a2,b,NNdist = √(0.5),NUnique = 2,SiteType= [1,2],refSites= [Rvec(0,0,1),Rvec(0,0,2)])
    end

    function generateLattice(L)
        B = SpiderWebBasis()
        sites = [Rvec(i,j,nb) for i in 0:L-1 for j in 0:L-1 for nb in eachindex(B.b)]
    end

    function generatePlaquettes(L,S=1/2)
        B = SpiderWebBasis()
        sites = []
    end

end # module SpiderWebModel
