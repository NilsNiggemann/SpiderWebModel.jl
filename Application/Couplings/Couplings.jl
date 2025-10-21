using SpinFRGLattices
using SpinFRGLattices.Octochlore: spin
using SpinFRGLattices.StructArrays
using CairoMakie
using SpinFRGLattices.StaticArrays
using MakieHelpers, GraphMakie
cd(@__DIR__)
##
function getCouplingsToS1(RefSite = Rvec(0, 0, 1))
    L = 4
    S1Terms = spin{Int,Rvec_2D}[]
    for n1 = (-L+RefSite.n1):(L+RefSite.n1), n2 = (-L+RefSite.n2):(L+RefSite.n2)
        isodd(n1 + n2) && continue
        sites = [
            Rvec(n1 + 1, n2, 1),
            Rvec(n1 + 1, n2 + 1, 1),
            Rvec(n1, n2 + 1, 1),
            Rvec(n1 - 1, n2 + 1, 1),
            Rvec(n1 - 1, n2, 1),
            Rvec(n1 - 1, n2 - 1, 1),
            Rvec(n1, n2 - 1, 1),
            Rvec(n1 + 1, n2 - 1, 1),
        ]
        sgns = [1, 1, -1, -1, 1, 1, -1, -1]
        plaq = [spin(r, s) for (r, s) in zip(sgns, sites)]
        for s1 in plaq
            for s2 in plaq
                if s1.site == RefSite
                    push!(S1Terms, spin(s1.fac * s2.fac, s2.site))
                elseif s2.site == RefSite
                    push!(S1Terms, spin(s1.fac * s2.fac, s1.site))
                end
            end
        end
    end
    return StructArray(S1Terms)
end

##
L = 1
AllPoints = [(n1, n2) for n1 = (-L):L for n2 = (-L):L]
inbox(x,L=1) = abs(x[1]) <= L && abs(x[2]) <= L
inbox(x::Rvec_2D,L=1) = abs(x.n1) <= L && abs(x.n2) <= L
filter!(inbox, AllPoints)

PlotPoints = Set(copy(AllPoints))

AllTerms = Dict()

function inBoxCouplings(p0)
    a = SpinFRGLattices.Octochlore.reduceCouplings(getCouplingsToS1(Rvec(p0..., 1)))
    filter!(x -> x.fac != 0, a)
    filter!(x -> inbox(x.site), a)
    filter!(x-> x.site != Rvec(p0..., 1), a)
    aDF = StructArray(;
        n1 = StructArray(a.site).n1,
        n2 = StructArray(a.site).n2,
        fac = a.fac,
    )
end

function getBonds(AllPoints)
    AllBonds = Dict()
    for p0 in AllPoints
        couplings = inBoxCouplings(p0)
        for term in couplings
            site = (term.n1, term.n2)
            fac = term.fac
            sitepair = sort(SA[p0, site])
            AllBonds[sitepair] = fac
        end
    end
    return AllBonds
end
##
AllBonds = getBonds(AllPoints)
AllBonds[(0,0),(0,0)] = 0

function getAdjacencyMatrix(AllBonds)
    nodes = Set()
    for (sitepair, _) in AllBonds
        push!(nodes, sitepair[1])
        push!(nodes, sitepair[2])
    end
    nodes = collect(nodes)
    node_indices = Dict(node => i for (i, node) in enumerate(nodes))
    n = length(nodes)
    adjacency_matrix = zeros(Float64, n, n)

    for (sitepair, fac) in AllBonds
        i, j = node_indices[sitepair[1]], node_indices[sitepair[2]]
        adjacency_matrix[i, j] = fac
        adjacency_matrix[j, i] = fac  # Assuming the graph is undirected
    end

    return adjacency_matrix, nodes
end

using GraphMakie
using GraphMakie.Graphs

function plotBondsAsGraph(AllBonds)
    adj_matrix, nodes = getAdjacencyMatrix(AllBonds)

    g = SimpleGraph(adj_matrix)

    fig = Figure(size = (400, 400), 
    backgroundcolor = :transparent
    )
    ax = Axis(
        fig[1, 1],
        aspect = 1,
        xgridcolor = :black,
        ygridcolor = :black,
        backgroundcolor = (:white, 0.0),
        xticks = [-1, 0, 1],
        yticks = [-1, 0, 1],
        xticklabelsvisible = false,
        yticklabelsvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        xgridvisible = false,
        ygridvisible = false,
        bottomspinevisible = false,
        topspinevisible = false,
        leftspinevisible = false,
        rightspinevisible = false,
    )

    curvedists = zeros(length(edges(g)))

    # sides
    curvedists[1] = 0.2
    curvedists[18] = -0.2
    curvedists[9] = 0.2
    curvedists[4] = -0.2

    # diags
    curvedists[8] = -0.2
    curvedists[5] = -0.085

    # bluediags
    curvedists[11] = -0.1
    curvedists[13] = 0.1
    curvedists[14] = 0.1
    curvedists[16] = -0.1
    
    # cross
    curvedists[17] = 0.2
    curvedists[10] = 0.08
    
    edgecolor = Dict(-4 => :blue, -2 => :blue, 2 => :red)

    edgeweights = [adj_matrix[e.src,e.dst] for e in edges(g)]

    edge_width = [3abs(e)^1.2 for e in edgeweights]
    edge_color = [edgecolor[e] for e in edgeweights]
    
    # edgecolors = [:black for i in 1:ne(g)]
    # edgecolors[4] = edgecolors[7] = :red

    # edgecolors = [:black for i in 1:ne(g)]
    # edgecolors[4] = edgecolors[7] = :red
    # node_marker = [isodd(sum(x)) ? '×' : '•' for x in nodes]
    node_marker = ['•' for x in nodes]

    p = graphplot!(
        ax,
        g;
        layout = Point.(nodes),
        # edgelabels = weights,
        nodecolor = :black,
        node_marker,
        # nlabels=repr.(1:nv(g)),
        # elabels = repr.(1:ne(g)),
        edge_color,
        edge_width,
        edgelabelsize = 10,
        node_size = 60,
        curve_distance=curvedists,
        curve_distance_usage=true
    )
    p.elabels_rotation[] = Dict(i => i == 2 ? 0.0 : Makie.automatic for i in 1:ne(g))
    # p.elabels_side[] = Dict(i => :right for i in [6,7])
    # p.elabels_offset[] = [Point2f(0.1, -0.1) for i in 1:ne(g)]
    # p.elabels_offset[][5] = Point2f(-0.4,0)
    # p.elabels_offset[] = p.elabels_offset[]

    p.elabels_shift[] = [0.58 for i in 1:ne(g)]
    # p.elabels_shift[][1] = 0.6
    # p.elabels_shift[][7] = 0.4
    # p.elabels_shift[] = p.elabels_shift[]

    # p.elabels_distance[] = Dict(8 => 30)
    # save("bonds_graph.pdf", fig)
    save("../figs/PaperFigs/couplings_All_plaq.svg", fig)
    return fig
end

plotBondsAsGraph(AllBonds)