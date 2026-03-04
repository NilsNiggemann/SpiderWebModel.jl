using Pkg
Pkg.activate(".")
import SpiderWebModel as SW
using CairoMakie
using SpiderWebModel.LoopVectorization
ALLGS_1 = SW.getAllGS(1)
# ALLGS_12_test = SW.getAllGS(1 / 2)

function checkTiling_constraint_open_!(curr_config,states,i,j,k,l)
    S_bottom_left = states[i]
    S_bottom_right = states[j]
    S_top_left = states[k]
    S_top_right = states[l]
    Lx, Ly = size(states[begin])
    curr_config[1:Lx, 1:Ly] .= S_bottom_left
    curr_config[Lx+1:2Lx, 1:Ly] .= S_bottom_right
    curr_config[1:Lx, Ly+1:2Ly] .= S_top_left
    curr_config[Lx+1:2Lx, Ly+1:2Ly] .= S_top_right

    return SW.fulFillsConstraint(curr_config)
end

function check_partial(Conf,plaquettes)
    for I in plaquettes
        i,j = Tuple(I)
        P = SW.getPlaquette(Conf, i, j)
        c = SW.constraint(P)

        if c ≠ 0
            return false
        end
    end

    return true
end
check_partial(Conf,::Nothing) = true 

function checkAllTilings_periodic(states,bufferConfig_template,check_1=nothing,check_2=nothing,check_3=nothing)
    Lx,Ly = size(states[1])
    @assert Lx == Ly
    @assert iseven(Lx)
    
    counters = zeros(Int, length(states))
    for i in eachindex(states)
    # for i in eachindex(states)
        bufferConfig = copy(bufferConfig_template)

        bufferConfig[1:Lx, 1:Ly] .= states[i]
        bufferConfigparent = parent(parent(bufferConfig))

        for j in eachindex(states)
            LoopVectorization.@turbo bufferConfigparent[Lx+1:2Lx, 1:Ly] .= states[j]
            check_partial(bufferConfig,check_1) || continue
            for k in eachindex(states)
                LoopVectorization.@turbo bufferConfigparent[1:Lx, Ly+1:2Ly] .= states[k]
                check_partial(bufferConfig,check_2) || continue
                for l in eachindex(states)
                    LoopVectorization.@turbo bufferConfigparent[Lx+1:2Lx, Ly+1:2Ly] .= states[l]
                    check_partial(bufferConfig,check_3) || continue

                    # if SW.fulFillsConstraint(bufferConfig)
                    counters[i] += 1
                    # end
                end
            end
        end
        if i % 1000 == 0
            println("i = $i, counter = $(counters[i])")
            flush(stdout)
        end
        # println("i = $i, counter = $(counters[i])")
        # flush(stdout)
    end
    return sum(counters)
    # for i in 1:length(states)
    #     if checkTiling_constraint_open_!(bufferConfig,states, i, i,i,i)
    #         counter += 1
    #     end
    # end

end
##
Lx = 6
Ly = 2
@time a = SW.constructAllConfigs(Lx, Ly, ALLGS_1)

@time a_rec = SW.stencilConfig.(SW.fillEmptyStates(a, Lx, Ly, ALLGS_1))
a_rec = parent.(parent.(a_rec)) # Convert to arrays
##
Spin = 1.
bufferConfig = SW.stencilConfig(fill(Spin,2Lx, 2Ly),Spin,boundaryCondition=:periodic)

check1 = [(4,2), (5,3), (1,3), (8,2)]
check2 = [(3,1),(2,8), (2,4), (3,5)]
check3 = collect(SW.plaquetteIterator(bufferConfig, true)) 
check3 = setdiff!(check3, check1)
check3 = setdiff!(check3, check2)
check3 = setdiff!(check3, [(2,2), (3,3), (2,6) , (3,7), (6,2) , (7,3), (6,6), (7,7), ])
# SW.plotSpinConfig(bufferConfig)
# vlines!([4.5],color = :red)
# hlines!([4.5],color = :red)
# scatter!(Point.(check1), color = :blue, markersize = 14)
# scatter!(Point.(check2), color = :green, markersize = 14)
# scatter!(Point.(check3), color = :red, markersize = 14)
# current_figure()
##
checkAllTilings_periodic(rand(a_rec,1000),bufferConfig,check1,check2,check3)
times = Float64[]
lens = [1000,2000,2500,3000,3500,3750,4000,5000]
for len in lens
    confs = rand(a_rec, len)
    time = @elapsed checkAllTilings_periodic(confs,bufferConfig,check1,check2,check3)
    println("Number of tilings: $numTil")   
    push!(times, time)
end

##
scatterlines(lens, times .^0.25,  color = :blue)
