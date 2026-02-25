import SpiderWebModel as SW

using SpiderWebModel.CairoMakie
using SpiderWebModel.FFTW
using Statistics
function getSS_Corrs(confs)
    Lx,Ly = size(confs[1])
    Nsites = Lx*Ly
    SiSj = zeros(Nsites,Nsites)
    Sq = zeros(Lx,Ly)
    for conf in confs
        conf_vec = reshape(conf,:)
        SiSj .+= conf_vec * conf_vec'
        Sq .+= abs2.(FFTW.fft(conf))
    end
    denom = length(confs)
    return SiSj ./ denom, Sq ./ denom
end
function SC(conf)
    return  SW.stencilConfig(conf, 1, boundaryCondition = :periodic)
end
function sc(conf)
    return  SW.stencilConfig(conf, 0.5, boundaryCondition = :periodic)
end
function trueMomenta(kmin,kmax,L)
    nmin = floor(Int,L*kmin/(2pi))
    nmax = ceil(Int,L*kmax/(2pi))
    # return 1/(2pi*L*100) .* nmin:nmax
    return (nmin : nmax) .* 2pi/L
end
function getSiSj_reconstructed(Sq)
    # kx = trueMomenta(-0.5pi,1.5pi,size(Sq,1))
    # ky = trueMomenta(-0.5pi,1.5pi,size(Sq,2))
    # SiSj_recon = zeros(size(Sq))
    return real(FFTW.ifft(Sq)) / prod(size(Sq))
end
function getSS_Tot(SiSj)
    Lx,Ly = size(SiSj)[1:2]
    SS_tot = zeros(Lx,Ly)
    for x in axes(SiSj,1),y in axes(SiSj,2)
        Sj = circshift(SiSj[x,y,:,:],(-x+1,-y+1))
        SS_tot .+= Sj
    end
    return SS_tot/(Lx*Ly)
end
function getDists_periodic(Lx,Ly)
    dists = zeros(Lx,Ly)
    for x in 1:Lx, y in 1:Ly
        dx = min(x-1,Lx-x+1)
        dy = min(y-1,Ly-y+1)
        dists[x,y] = sqrt(dx^2 + dy^2)
    end
    return dists
end
##
# Load or generate an S=0 configuration (spin-1/2 ground state)
L = 36
# S0_conf = SW.stencilConfig(zeros(L, L), 1, boundaryCondition = :periodic)
# S0_conf = sc(SW.getStairCase(L))
S0_conf = SW.get4x4PeriodicSpinConf(L,6)

# Set up GFMC parameters
τ = 0.3
NSteps = 6000
NWalkers = 2000
ψG = SW.RKFunction()  # or SW.FullVariationalGuidingFunction(params) if available
w_avg_estimate = 0.25 * L^2
mu = 0.9
Hxx = SW.Hxx_RK(mu)
CT = SW.ContinuousTimeMethod(τ, w_avg_estimate, Hxx)

# Run GFMC to sample configurations
results = SW.startManyWalkerGFMC(S0_conf, CT, NWalkers, NSteps, ψG; equilibration_steps=NSteps÷5, pre_equilibration_steps=10000, scatter_fraction=1.0)

# Access sampled configurations
sampled_configs = results.SaveConfigs


# confs = reshape(sampled_configs, L, L, :) ./2
confs = sampled_configs[:,:,1,:] ./ 2  # Convert from {-2,0,2} to {-1,0,1} 
println("Sampled $(size(confs, 3)) spin-1 configurations.")


##
SiSj, Sq = getSS_Corrs(eachslice(confs, dims=3))
# SiSj_tot = reshape(SiSj[2,:], L,L)
SiSj_tot = getSS_Tot(reshape(SiSj,L,L,L,L))
SiSj_recon = getSiSj_reconstructed(Sq)
dists = getDists_periodic(L,L)

##
# SiSj_tot = circshift(SiSj_tot,(-L÷2+1,-L÷2+1))
fig = Figure(size = (900,300))
ax1 = Axis(fig[1,2],title = L"\frac{1}{N} \sum_i \langle S(r_i) S(r_j-r_i) \rangle",xlabel = "Site Index",ylabel = "Site Index",aspect=1)
hm1 = heatmap!(ax1,SiSj_tot)
ax2 = Axis(fig[1,1],title = L"S(\mathbf{q})",xlabel = L"q_x",ylabel = L"q_y",aspect=1)
heatmap!(ax2,Sq)
ax3 = Axis(fig[1,3],title = L"FT(S(\mathbf{q}))",aspect=1)
# SiSj_recon = circshift(SiSj_recon,(-L÷2+1,-L÷2+1))
hm3 = heatmap!(ax3,SiSj_recon,colorrange = extrema(SiSj_tot))
Colorbar(fig[1,4],hm3)
ax4 = Axis(fig[1,5],xlabel = "Distance",ylabel = L"\frac{1}{N} \sum_i \langle S(r_i) S(r_j-r_i) \rangle",aspect=1,xscale = log10,yscale = log10)
xlims!(ax4,1,maximum(dists))
# scatter!(ax4,[dists[1,i] for i in 1:L],abs.([SiSj_tot[1,i] for i in 1:L]),markersize=9,color=:red)
scatter!(ax4,dists[:],abs.(SiSj_tot[:]),markersize=5,color=(:black,0.2))
unique_dists = sort(unique((dists[:])))
lines!(ax4,unique_dists,6.5 ./unique_dists.^3,color=:red,linestyle=:dash,label = L"\frac{1}{r^3}")
axislegend(ax4,position=:lb)
fig