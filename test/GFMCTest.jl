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

##

function testResults(S,method)
    HStair = SW.generateHilbertSpace(SW.SpinConfig(S))
    ExSol = SW.SolveHKrylov(HStair.H)
    E0 = ExSol.values[1]

    v0 = ExSol.vectors[1]
    HConfs = SW.spinConfig.(HStair.AllStates,Ref(SW.SpinConfig(S)),Ref(HStair.plaqMapping))
    ψG = SW.PlaquetteNumberGuidingFunction(0.197)
    results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,method,12,1200,ψG,equilibration_steps=10,pre_equilibration_steps=100,scatter_fraction=0.5) for i in 1:18])
    # display(plotEnergies(results,method,E0))
    @testset "Energy" begin
        energies = last.(SW.getEnergies.(results,1,20))
        emean = SW.mean(energies)
        estd = SW.std(energies)
        @test emean - estd < E0 < emean + estd
    end
    
    Sqs = SW.getSqsGFMC(results,20) ./ 2
    SqEx = real(SW.getStructureFac(HConfs,v0).Sq)

    @testset "Structure factor" begin
        Sqmean = SW.mean(Sqs)
        SqDiff = abs.(Sqmean .- SqEx)
        Sqstd = SW.std(Sqs)
        @test sum(SqDiff .- Sqstd .>= 0) == length(SqDiff)
    end
end
##
SW.Random.seed!(1234)

S = SW.stencilConfig(parent(SW.getStairCase(10)),1/2,
boundaryCondition = :periodic
)
testResults(S,SW.DiscreteTimeMethod(0.,2,67))

