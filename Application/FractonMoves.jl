import SpiderWebModel as SW
using CairoMakie, GLMakie
##
St = SW.SpinConfig(SW.getStairCase(15) |> Array,1/2)
display(SW.plotFractons(St))

for x in [7,9,11,13,15,1,3,5]
    St[x,7] = -St[x,7]
    display(SW.plotFractons(St))
end
##


##

function AnimatespinConfig(StartState,path)
    time = Observable(1)
    xsites = [0,7,9,11,13,15,1,3,5]

    framerate = 1
    timestamps = eachindex(xsites)

    St = copy(StartState)

    function flipConf!(time)
        x = xsites[time]
        x == 0 && return St
        St[x,7] = -St[x,7]
        St
    end
    
    AllConfs = [copy(flipConf!(t)) for t in timestamps]
    AllCharges = [SW.getPlotPoints_all(Conf) for Conf in AllConfs]

    cols_arr = getproperty.(AllCharges,:cols)
    sizes_arr = getproperty.(AllCharges,:sizes)
    points_arr = getproperty.(AllCharges,:points)


    spinConf = @lift( AllConfs[$time] )
    cols = @lift( cols_arr[$time] )
    sizes = @lift( sizes_arr[$time] )
    points = @lift( points_arr[$time] )
    # charges = SW.getCharges(spinConf)
    # fig = SW.plotFractons(St)
    fig = Figure(size = (1000, 1000))
    ax = Axis(fig[1, 1];SW.getConfigAxis(St)...)

    hm = heatmap!(ax,spinConf;colormap = :grays,colorrange = (-0.5,0.5))
    constrPoints = [Point2(i,j) for i in axes(St.Mat,1),j in axes(St.Mat,2) if (i+j) %2 == 0]

    scatter!(ax,constrpoints,markersize = sizes,color = cols,marker = :rect)
    scatter!(ax,points,markersize = sizes,color = cols,marker = :rect)
    # scatter!(ax,points)
    translate!(hm, 0, 0, -100)


    record(fig, path, timestamps;
            framerate = framerate) do t
        time[] = t
    end
end

St = SW.SpinConfig(SW.getStairCase(15) |> Array,1/2)

AnimatespinConfig(St,"Animation/StairCaseFracton.mp4")
##

St = SW.SpinConfig(SW.getStairCase(15) |> Array,1/2)
SW.flipPlaquette!(St,12,7)
AnimatespinConfig(St,"Animation/StairCaseFlipFracton.mp4")

##

St = SW.constructConfigPath(SW.DictAlgorithm(),7,7,SW.ALLGS_S12)

AnimatespinConfig(St,"Animation/RandStateFracton.mp4")

##


function AnimatePlaqFlip(StartState,path)
    time = Observable(1)
    St = copy(StartState)
    Plaq = SW.getApplicablePlaquettes(St)[begin + 30]
    Pij = SW.getPlaquette(St,Plaq...)

    sites = [(0,0),(3,1),(3,2),(3,3),(2,3),(1,3),(1,2),(1,1),(2,1),(3,1)]
    framerate = 1
    timestamps = 1:9


    function flipConf!(time)
        x = sites[time]
        x == (0,0) && return St
        x = CartesianIndex(x)
        Pij[x] = -Pij[x]
        St
    end
    
    AllConfs = [copy(flipConf!(t)) for t in timestamps]
    AllCharges = [SW.getPlotPoints_all(Conf) for Conf in AllConfs]

    cols_arr = getproperty.(AllCharges,:cols)
    sizes_arr = getproperty.(AllCharges,:sizes)
    points_arr = getproperty.(AllCharges,:points)
    sites_arr = [Point2.(sites[2:t-1]) .+ Ref(Point2(Plaq) - Point2(2,2)) for t in 2:length(sites)]
    # return sites_arr, timestamps
    spinConf = @lift( AllConfs[$time] )
    cols = @lift( cols_arr[$time] )
    sizes = @lift( sizes_arr[$time] )
    points = @lift( points_arr[$time] )
    sites = @lift (sites_arr[$time])

    # charges = SW.getCharges(spinConf)
    # fig = SW.plotFractons(St)
    fig = Figure(size = (1000, 1000))
    ax = Axis(fig[1, 1];SW.getConfigAxis(St)...)

    hm = heatmap!(ax,spinConf;colormap = :grays,colorrange = (-0.5,0.5))
    constrPoints = [Point2(i,j) for i in axes(St.Mat,1),j in axes(St.Mat,2) if (i+j) %2 == 0]

    scatter!(ax,constrPoints,markersize = 25,color = :grey,marker = '×')
    scatter!(ax,points,markersize = sizes,color = cols,marker = :rect)
    scatter!(ax,Point2(Plaq),markersize = 20,color = :grey,marker = :circle)
    scatterlines!(ax,sites,markersize = 0,color = :grey,linestyle = :dash,linewidth = 4)
    # # scatter!(ax,points)
    translate!(hm, 0, 0, -100)


    record(fig, path, timestamps;
            framerate = framerate) do t
        time[] = t
    end
end

St = SW.SpinConfig(SW.getStairCase(15) |> Array,1/2)

AnimatePlaqFlip(St,"Animation/StairCasePlaqFlip.mp4")