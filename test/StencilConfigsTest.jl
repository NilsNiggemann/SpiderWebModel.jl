using Test
import SpiderWebModel as SW
@testset "Stencil Spin Configs" begin
    Sopen = SW.stencilConfig(zeros(18,18),1)

    @test length(collect(SW.plaquetteIterator(Sopen))) == 16^2 /2
    S = SW.stencilConfig(zeros(18,18),1;
    boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )
    conf = parent(S) 
    @test conf isa SW.Stencils.StencilArray
    @test SW.NPlaquettes(S) == 18^2
    
    P = (1,2)
    SW.applyPlaquette!(S,P...,1)
    
    plaq = SW.getPlaquetteSites(S,P...)
    @test plaq == SW.SVector{8,Int8}(2, -2, -2, 2, 2, -2, -2, 2)
    @test SW.fulFillsConstraint(S)
    
    allplaqs = collect(SW.plaquetteIterator(S))
    @test length(allplaqs) == 18^2 / 2

    @test SW.NPlaquettes(S) == 18^2 - 4*2 - 8 - 1
end

