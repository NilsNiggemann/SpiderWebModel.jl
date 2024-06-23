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
