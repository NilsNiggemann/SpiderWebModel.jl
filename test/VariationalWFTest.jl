import SpiderWebModel as SW
using Test
using SpiderWebModel.StaticArrays

@testset "utils" begin
    ψG = SW.JastrowFunction(SW.stencilConfig(zeros(4,4),1),Float64)

    params = SW.get_params(ψG)
    @test length(params) == length(ψG.m_i) + length(ψG.v_ij) + length(ψG.α)

    @testset "indexMapping" failfast = true begin
        for i in eachindex(params)
            type,k = SW._getParamsTypeAndIndex(params,i )
            @test SW.remap_index(type,k,params) == i
        end
    end
end

##
function TestWFRatio(ψG,S;tol=1e-10)
    Buff = SW.allocate_GWF_buffer(ψG,S)
    W = SW.spiderWebWalker(copy(S),collect(SW.plaquetteIterator(S)))
    SW.getMoves!(W)
    psiname = SW.guidingfunc_name(ψG)

    SW.compute_GWF_buffer!(Buff,ψG,W)
    @testset "$psiname ψ(x´) / ψ(x)" begin
        for m in W.moves
            (i,j,op) = m

            Spr = copy(S)
            SW.applyPlaquette!(Spr,i,j,op)
            
            @test SW.guidingfuncRatio(ψG,W,m,Buff) - ψG(Spr)/ψG(S) ≈ 0 atol=tol
            
        end
    
    end
end


@testset "wavefunctions" begin
    
    
    S = SW.stencilConfig(zeros(8,8),1;
    boundaryCondition = :open
    )
    S .= SW.periodicStateDenseLoops(size(S,1))

    ψRBM = SW.RBMSpin1(S,1,Float64)
    SW.rand!(SW.get_params(ψRBM))  .*= 1e-2
    TestWFRatio(ψRBM,S)
    
    ψLocalPlaqs = SW.localPlaquetteGuidingFunction(S,0.,Float64)
    SW.rand!(SW.get_params(ψLocalPlaqs))
    TestWFRatio(ψLocalPlaqs,S)

    ψFull = SW.fullVariationalFunction(S,0.,Float64)
    SW.rand!(SW.get_params(ψFull))  .*= 1e-2
    TestWFRatio(ψFull,S)

    ψMag = SW.orderGuidingFunction(S,0.,Float64)
    SW.rand!(SW.get_params(ψMag))  .*= 1e-2
    TestWFRatio(ψMag,S)

    ψPlaqRBM = SW.PlaquetteRBM(S,1,Float64)
    SW.rand!(SW.get_params(ψPlaqRBM))  .*= 1e-2
    TestWFRatio(ψPlaqRBM,S)

    ψJastrow = SW.JastrowFunction(S,Float64)
    SW.rand!(SW.get_params(ψJastrow))  .*= 1e-2
    SW.get_v_ij(ψJastrow) .= SW.Symmetric(SW.get_v_ij(ψJastrow))
    TestWFRatio(ψJastrow,S)


end
##
@testset "Symmetries" begin
    S = SW.stencilConfig(zeros(12,12),1;
    boundaryCondition = :periodic
    )
    ψG = SW.JastrowFunction(S)

    params = SW.get_params(ψG)
    Symm = SW.TranslationalSymmetry(SA[2, 0], SA[0, 2])
    ψGSymm = SW.symmetrize(ψG,Symm,S)
    newparas = 1:length(ψGSymm.uniqueInds)
    
    SW.add_reconstructedFullParams!(ψG,ψGSymm.indicesMapping,newparas)

    @test ψG.m_i[10] - ψG.m_i[12] ≈ 0
    @test ψG.m_i[15] - ψG.m_i[17] ≈ 0

    vij_1 = SW.PeriodicMatrix(reshape(ψG.v_ij[1,:],size(S)))
    vij_2 = SW.PeriodicMatrix(reshape(ψG.v_ij[3,:],size(S)))[3:end+2,1:end]
    @test vij_1 == vij_2

    vij_2 = SW.PeriodicMatrix(reshape(ψG.v_ij[1+2*size(S,1),:],size(S)))[1:end,3:end+2]
    @test vij_1 == vij_2


    
    Symm = SW.TranslationalSymmetry(SA[2, 2], SA[-2, 2])

    ψGSymm = SW.symmetrize(ψG,Symm,S)
    newparas = 1:length(ψGSymm.uniqueInds)
    
    SW.add_reconstructedFullParams!(ψG,ψGSymm.indicesMapping,newparas)

    mxy = reshape(ψG.m_i,size(S))
    vIJ = reshape(ψG.v_ij,(size(S)...,size(S)...))

    @test mxy[3,3] - mxy[3-2,3+2] ≈ 0
    @test mxy[3,3] - mxy[3+2,3+2] ≈ 0

    @test vIJ[3,3,3,3] - vIJ[3-2,3+2,3-2,3+2] ≈ 0
    @test vIJ[3,3,3,3] - vIJ[3+2,3+2,3+2,3+2] ≈ 0

end