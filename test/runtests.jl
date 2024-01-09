using SpiderWebModel,Test
import SpiderWebModel as SW
##
@testset "plaquettes" begin
    c = [1,1,1,1,1,1,1,1] 
    P = SW.Plaquette(c...)
    @test SW.constraint(P) == SW.constraint(c) == 0

    c = [1,1,1,-1,1,1,1,1] .* 0.5
    P = SW.Plaquette(c...)
    @test SW.constraint(P) == SW.constraint(c)
end
##
@testset "tiling test" begin
    AllAllowedConfigs = SW.getAllGS(0.5)
    T1 = first(AllAllowedConfigs).x
    T2 = AllAllowedConfigs[2].x
    Tend = AllAllowedConfigs[end].x
    @test SW.canTileUpRight(T1,T1)
    @test SW.canTileDownRight(T1,T1)
    @test SW.canTileDownLeft(T1,T1)
    @test SW.canTileUpRight(T1,T2)
    @test !(SW.canTileUpRight(T1,Tend))
end
##
@testset "Hashing arrays" begin
    state = rand(-0.5:0.5,20,20)
    StateSet = Set([state])
    StateArr = [state]
    for i in eachindex(state)
        newstate = copy(state)
        push!(StateArr,newstate)
        newstate[i] = -newstate[i]
        push!(StateSet,newstate)
    end
    @test length(StateSet) == length(StateArr)
    
    for i in eachindex(state)
        newstate = copy(state)
        newstate[i] = -newstate[i]
        push!(StateSet,newstate)
    end
    
    @test length(StateSet) == length(StateArr)
end