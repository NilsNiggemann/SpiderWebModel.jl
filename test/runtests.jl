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