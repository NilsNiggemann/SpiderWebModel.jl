using Test
import SpiderWebModel as SW
##
function reconfig_is_arraged(list)
    uniqueElements = unique(list)
    is_arranged = true
    for (α,α´) in enumerate(list)
        if α != α´ && α ∈ uniqueElements
            is_arranged = false
            break
        end
    end
    return is_arranged
end
@testset "minimizeReconfiguration" begin
    @test SW.minimizeReconfiguration!([1, 2, 3, 4, 5, 6, 6, 10, 9, 11, 11, 12]) == [1, 2, 3, 4, 5, 6, 6, 11, 9, 10, 11, 12]
 
    l1 = sort!([1,4,2,5,2,1,1,3,2,6,12,2,1])
    @test reconfig_is_arraged(SW.minimizeReconfiguration!(l1))
end

@testset "SpiderWebWalker" begin
    S = SW.stencilConfig(zeros(18,18),1;
    boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )
    allplaqs = collect(SW.plaquetteIterator(S))
    W = SW.spiderWebWalker(S,allplaqs)

    moves = SW.getMoves!(W)
    
    @test length(moves) == length(allplaqs) * 2
    
    Conf = SW.get_config(W)
    SW.applyPlaquette!(Conf,1,2,1)
    
    SW.applyPlaquette!(Conf,2,3,1)
    
    SW.applyPlaquette!(Conf,18,17,-1)
    
    moves = SW.getMoves!(W)
    ni = SW.getNPlaq!(W)

    @test length(ni) == length(allplaqs)
    @test all(n == 0. || n == 1. || n == 2. for n in ni)
    @test sum(ni) == SW.NPlaquettes(Conf) != SW.NPlaquettes(S)

    @test SW.fulFillsConstraint(Conf)

end