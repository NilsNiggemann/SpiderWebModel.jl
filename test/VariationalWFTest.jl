import SpiderWebModel as SW
using Test

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