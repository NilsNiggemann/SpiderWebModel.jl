import SpiderWebModel as SW
import LinearAlgebra as LA
cd(@__DIR__)
include("../plottingUtils.jl")
function get_S_diag!(S)
    Sdiag = SW.PeriodicMatrix(Int8[0 2 0 -2; 2 0 -2 0; 0 -2 0 2; -2 0 2 0])
    S .= Sdiag[axes(S)...]
    return S
end

function get_angle(q)
    ang = atan(q[2], q[1])
    return ang < 0 ? ang + 2pi : ang
end
function get_neighbors(S)
    movelist = empty!([(1,1,1)])
    moves = SW.getMoves!(movelist,S)
    allConfs = empty!([S])
    for move in moves
        Spr = copy(S)
        SW.applyPlaquette!(Spr, move)
        push!(allConfs, Spr)
    end
    return allConfs
end

function get_perturbative_spin_corrs(S,alpha)
    neighbors = get_neighbors(S)
    neighbors_fl = [ComplexF32.(parent(n)) for n in neighbors]
    Sq = zeros(size(S))
    Sq_buffer = zeros(ComplexF32,size(S))
    plan = SW.FFTW.plan_fft(Sq_buffer)
    Nsites = length(S)
    for m in neighbors_fl
        # for n in neighbors_fl
            # m == n || continue
            # return Sq_buffer,plan,m
        LA.mul!(Sq_buffer, plan, m)
        Sq .+= abs2.(Sq_buffer)
        # end
    end
    Sq .*= alpha^2
    Sq .+= abs2.(plan*S)
    return Sq/((length(neighbors_fl)+1)*Nsites)
end

function S_q_FT2(qx, qy, K, U, W)
    """
    Evaluate the structure factor at momentum (qx, qy).
    
    S(q) = sqrt(K) * (cx - cy + 2*sx*sy)^2 / 
           (sqrt((cx-cy)^2 + 4*sx^2*sy^2) * sqrt(U/4 + W*((cx-cy)^2 + 4*sx^2*sy^2)))
    
    Arguments:
        qx, qy: momentum components
        K, U, W: coupling parameters
    """
    cx = cos(qx)
    sx = sin(qx)
    cy = cos(qy)
    sy = sin(qy)
    
    # Common intermediate term
    delta_c = cx - cy + 2 * sx * sy
    L_sq = (cx - cy)^2 + 4 * sx^2 * sy^2
    
    # Numerator: sqrt(K) * (cx - cy + 2*sx*sy)^2
    num = sqrt(K) * delta_c^2
    
    # Denominator: sqrt(L_sq) * sqrt(U/4 + W*L_sq)
    denom = sqrt(L_sq) * sqrt(U/4 + W * L_sq)
    
    return num / denom
end

##
S = SW.stencilConfig(zeros(220,220),1,boundaryCondition = :periodic)
get_S_diag!(S)
##
Sq = get_perturbative_spin_corrs(S,2)
##
# filterSq = [s < maximum(Sq) ? s : NaN for s in Sq]
filterSq = [s < maximum(Sq)/1.05 ? s : 0 for s in Sq]
# filterSq_inv = [s >= maximum(Sq)/1.1 ? s : NaN for s in Sq]
with_theme(theme_PiTicks()) do 
    W = 1
    fig = Figure(size = 0.85 .*(610,400))
    ax_pert = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect = 1,xticks = PiTicks([0,pi]),yticks = PiTicks([0,pi]),
    title = L"$$pert. theory",
    )
    ax_FT = Axis(fig[1,3],xlabel = L"q_x",ylabel = L"q_y",aspect = 1,xticks = PiTicks([0,pi]),yticks = PiTicks([0,pi]),
    ylabelvisible = false,yticklabelsvisible = false,
    title = L"field theory $W=%$W$"
    )
    ax_scale = Axis(fig[1,4],xlabel = L"q_x",ylabel = L"q_y",aspect = 1,xticks = PiTicks([-0.5pi,0,0.5pi]),yticks = PiTicks([-0.5pi,0,0.5pi]),
    ylabelvisible = false,#yticklabelsvisible = false,
    title = L"$$pert. theory $\times q^{-4}$ ",)
    ax_cuts = Axis(fig[2,1:4],xlabel = L"φ",ylabel = L"\mathcal{S}(q)/q^4",
    xticks = PiTicks([0,pi,2pi]),yticks = SimpleTicks()
    )
    norm(x,y) = sqrt(x^2 + y^2) + 1e-16
    norm(q) = norm(q[1],q[2])

    SqCont = SW.getSqCont(Sq,cutoffEnd=0)

    # qx = qy = trueMomenta(-0.7pi,0.7pi,size(S,1))
    qx = qy = trueMomenta(-0.5pi,1.5pi,size(S,1))
    qxzoom = qyzoom = trueMomenta(-0.5pi,0.5pi,size(S,1))
    
    Sq_pert = [SqCont(x,y) for x in qx, y in qy]
    # Sq_FT = [S_q_FT2(x,y,0.01,20000,0) for x in qx, y in qy]
    Sq_FT = [SqFieldTheory_full(x,y,1,W,1) for x in qx, y in qy]
    # Sq_FT = [SqFieldTheory(x,y,1,1e10) for x in qx, y in qy]
    
    ring(r,phi) = r .* Point2f(cos(phi), sin(phi))
    
    radii = LinRange(0.1pi,0.25pi,4)
    
    
    Sq_scale = [SqCont(x,y)/norm(x,y)^4 for x in qxzoom, y in qyzoom]
    # Sq_scale_filtered = [s < maximum(Sq_scale)/1.0005 ? s : 0 for s in Sq_scale]
    # return Sq_scale
    Sq_scale_filtered = [s == maximum(Sq_scale) ? 0 : s for s in Sq_scale]
    
    hm_pert = heatmap!(ax_pert,qx,qy, Sq_pert,colorrange = extrema(filterSq), highclip= :red,rasterize=true)
    hm_pert.rasterize=5
    q_grid = collect(Iterators.product(qx, qy))

    maxqs = findall(x -> x == maximum(Sq_pert), Sq_pert)
    qmax = [Point2(q_grid[q]...) for q in maxqs]

    scatter!(ax_pert, qmax, color = :red, marker = :rect,markersize = 3) 

    hm_FT = heatmap!(ax_FT,qx,qy, Sq_FT,rasterize=true)
    hm_FT.rasterize=5
    hm_scale = heatmap!(ax_scale,qxzoom,qyzoom, Sq_scale, colorrange = extrema(Sq_scale_filtered), highclip= :red,rasterize=true)
    # return fig
    hm_scale.rasterize=10
    
    maxqs = findall(x -> x == maximum(Sq_scale), Sq_scale)
    qmax = [Point2(q_grid[q]...) for q in maxqs]
    scatter!(ax_scale, qmax, color = :red, marker = :rect,markersize = 8)
    Colorbar(fig[1,2],hm_pert,
    # label = L"\mathcal{S}(\textbf{q})",
    height = Relative(0.95),vertical=true,flipaxis=true,
    # ticks = ([maximum(filterSq) ], [string(maximum(Sq))]), minorticksvisible = true, minorticksize = 10, minortickwidth = 2, minortickcolor = :red,
    )
    colors = cgrad(:brg, length(radii), categorical = true)
    phis = LinRange(0,2pi,900)
    line_fatness = 3*1/size(Sq,1)

    for (idx,q_abs) in enumerate(radii)
        # cut_pert = [SqCont(r*cos(ϕ),r*sin(ϕ)) for ϕ in phi]
        # cut_FT = [S_q_FT2(r*cos(ϕ),r*sin(ϕ),0.01,20000,0) for ϕ in phi]
        # qs = ring.(r,phi)

        path_raw = [ q_abs .*(cos(phi), sin(phi)) for phi in phis]
        path = filter( q -> q_abs - line_fatness <= norm(q) <= q_abs + line_fatness, q_grid)
        phis_corrected = [get_angle(q) for q in path]
        perm = sortperm(phis_corrected)
        phis_corrected = phis_corrected[perm]
        path = path[perm]
        push!(phis_corrected, phis_corrected[1] + 2pi)
        push!(path, path[1])

        cut_scale = [SqCont(q[1],q[2])/norm(q)^4 for q in path]
        # cut_scale = [SqCont(q[1],q[2])/q_abs^4 for q in path]
        # lines!(ax_scale,path, color = colors[idx])
        # lines!(ax_cuts, phi, cut_pert, color=:black, label = r == radii[1] ? L"$$" : "")
        # lines!(ax_cuts, phi, cut_FT, color=:red, label = r == radii[1] ? L"$$FT$$" : "")
        scatter!(ax_scale, path, color = colors[idx], markersize = 1.5, )
        # scatterlines!(ax_cuts, phis_corrected, cut_scale, color = colors[idx],markersize = 4)
        lines!(ax_cuts, phis_corrected, cut_scale, color = colors[idx])
    end
    # rowsize!(fig.layout,1, Relative(0.55))
    # colsize!(fig.layout,5, Relative(0.4))
    Label(fig[1, 1,TopLeft()], L"(a)$$", padding = (-40,0,-25, -20),fontsize = 16)
    Label(fig[1, 3,TopLeft()], L"(b)$$", padding = (-40,15,-25, -20),fontsize = 16)
    Label(fig[1, 4,TopLeft()], L"(c)$$", padding = (-20,0,-25, -20),fontsize = 16)
    Label(fig[2, 1,TopLeft()], L"(d)$$", padding = (-40,0,-25, -20),fontsize = 16)
    save("../../figs/PaperFigs/SqPerturbation.pdf",fig,px_per_unit =3)
    fig
end
# heatmap!(filterSq_inv,colormap = :red)
# current_figure()
# heatmap(filterSq,colorrange = (0,0.0001))