import SpiderWebModel as SW
using Test

function TestWFRatio(ψG,S;tol=1e-10)
    Buff = SW.allocate_GWF_buffer(ψG,S)
    W = SW.spiderWebWalker(copy(S),collect(SW.plaquetteIterator(S)))
    SW.getMoves!(W)
    psiname = SW.guidingfunc_name(ψG)
    
    SW.fill_GWF_buffer!(Buff,ψG,W)
    @testset "$psiname ψ(x´) / ψ(x)" begin
        for m in W.moves
            (i,j,op) = m

            Spr = copy(S)
            SW.applyPlaquette!(Spr,i,j,op)
            
            @test SW.guidingfuncRatio(ψG,W,m,Buff) - ψG(Spr)/ψG(S) ≈ 0 atol=tol
            
        end
    
    end
end

S = SW.stencilConfig(zeros(8,8),1;
boundaryCondition = :open
)
S .= SW.periodicStateDenseLoops(size(S,1))

@testset "wavefunctions" begin
    
    ψRBM = SW.RBMSpin1(S,1,Float64)
    SW.rand!(SW.get_params(ψRBM))
    SW.get_params(ψRBM) .*= 1e-2
    TestWFRatio(ψRBM,S)
    
    ψLocalPlaqs = SW.localPlaquetteGuidingFunction(S,0.,Float64)
    SW.rand!(SW.get_params(ψLocalPlaqs))
    TestWFRatio(ψLocalPlaqs,S)

    ψFull = SW.fullVariationalFunction(S,0.,Float64)
    SW.rand!(SW.get_params(ψFull))
    SW.get_params(ψFull) .*= 1e-2

    ψMag = SW.orderGuidingFunction(S,0.,Float64)
    SW.rand!(SW.get_params(ψMag))
    SW.get_params(ψMag) .*= 1e-2


    TestWFRatio(ψMag,S)
end