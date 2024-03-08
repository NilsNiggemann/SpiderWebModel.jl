import SpiderWebModel as SW
using SpiderWebModel
using SpiderWebModel.CairoMakie
using SpiderWebModel.DSP
using SpiderWebModel.StaticArrays
import SpiderWebModel.LatticeFFTs
using SpiderWebModel.LatticeFFTs.FFTW
using MakieHelpers
using CairoMakie.Makie.ColorSchemes

##
iceRule = SA[
    -1 -1 1
    1 0 1
    1 -1 -1
]

function textheatmap!(ax, data; kwargs...)
    for i in axes(data, 1), j in axes(data, 2)
        isapprox(data[i, j], 0, atol = 1e-2) && continue
        # txtcolor = data[i, j] < 0 ? :white : :black
        txtcolor = data[i, j] < 0 ? :white : :white
        text!(
            ax,
            "$(MakieHelpers._simplify(data[i,j]))",
            position = (i, j),
            color = txtcolor,
            align = (:center, :center),
            fontsize = 10;
            kwargs...,
        )
    end
end
function plotCharges(Charge; kwargs...)
    fig, ax, hm = heatmap(
        parent(Charge),
        colormap = :balance,
        colorrange = ±(Charge),
        axis = SW.getConfigAxis(Charge),
        figure = (; size = (800, 800));
        kwargs...,
    )
    textheatmap!(ax, parent(Charge))
    Colorbar(fig[1, 2], hm)
    fig
end

getk(FT) = LinRange(0, (1 - 1 / size(FT, 1)) * 2pi, size(FT, 1))
mean(x) = sum(x) / length(x)

function plotChargeFT(Charges::AbstractVector; kwargs...)
    ChargeFFT = mean(real(fft(Charge)) for Charge in Charges)
    plotChargeFFT(ChargeFFT; kwargs...)
end

function plotChargeFT(Charge::AbstractMatrix; kwargs...)
    ChargeFFT = (real(fft(Charge)))
    plotChargeFFT(ChargeFFT; kwargs...)
end

function plotChargeFFT(ChargeFFT; kwargs...)
    with_theme(theme_PiTicks()) do
        fig = Figure()
        ax = Axis(fig[1, 1], title = L"C(\mathbf{k})", aspect = 1)

        k = getk(ChargeFFT)

        hm = heatmap!(
            ax,
            k,
            k,
            ChargeFFT,
            colormap = :balance,
            colorrange = ±(ChargeFFT);
            kwargs...,
        )
        Colorbar(fig[1, 2], hm, vertical = true)
        fig
    end
end

Mag(Config) = sum(Config)

StaggeredMag(Conf) =
    sum((-1)^(x + y) * Conf[x, y] for x in axes(Conf, 1), y in axes(Conf, 2))
numPosCharges(charge) = sum(x for x in charge if x > 0)

function deconvolve(Charge, kernelFFT, S_tot = 0, SStag = 0; normalize = true)
    ChargeFFT = FFTW.fft(Charge)

    # count(iszero,(kernelFFT)) <= 2 || error("Kernel has more than 2 zeros")

    SpinFFT = ChargeFFT ./ kernelFFT

    isnan(SpinFFT[1, 1]) && (SpinFFT[1, 1] = S_tot)
    isnan(SpinFFT[end÷2+1, end÷2+1]) && (SpinFFT[end÷2+1, end÷2+1] = SStag)
    for i = 1:size(SpinFFT, 1), j = 1:size(SpinFFT, 2)
        isnan(SpinFFT[i, j]) && (SpinFFT[i, j] = 0)
    end
    @info "" 2 * sum(real, SpinFFT) sum(imag, SpinFFT) size(SpinFFT) 2Mag(SpinFFT) StaggeredMag(
        SpinFFT,
    ) S_tot SStag

    # fillVal = -4.5
    with_theme(theme_PiTicks()) do
        k = getk(SpinFFT)

        display(heatmap(k, k, real.(SpinFFT), axis = (; aspect = 1)))
    end

    Spin = real.(FFTW.ifft(SpinFFT))
    Lx, Ly = size(Spin)
    newSpin = SW.CircularArrays.CircularArray(Spin)[Lx÷2+3:Lx+Lx÷2, Ly÷2+3:Ly+Ly÷2]
    # newSpin = SW.CircularArrays.CircularArray(Spin)[Lx÷2+3:Lx+Lx÷2, Ly÷2+3:Ly +Ly ÷ 2 ]
    # Spin = SW.SpinConfig(round.(Bool,newSpin),1/2)
    if normalize
        return SW.SpinConfig(copysign.(0.5, newSpin), 1 / 2)
    else
        Spin = SW.SpinConfig(round.(newSpin, digits = 10), 1 / 2)
    end
    # SW.floatSpinConfig(Spin)
end

##
# SW.flipPlaquette!(Conf,10,9)
# Conf = SW.SpinConfig(-ones(20,20)/2,1/2)
# Conf[11,9] = -Conf[11,9]
# Conf = SW.SpinConfig(ones(20,20)/2,1/2)
# Conf = SW.SpinConfig([-0.5+iseven(i+j) for i in 1:20,j in 1:20],1/2)

# SW.flipPlaquette!(Conf,19,3)
# Conf = SW.constructConfigPath(17,17,SW.ALLGS_S12)
Conf = SW.periodicState5x5(25)
# Conf = SW.periodicState6x6(24)
# Conf = SW.periodicState6x6_3(18)
# Conf = SW.getStairCase(24)
# for i in 1:24
#     SW.flipPlaquette!(Conf,rand(SW.getApplicablePlaquettes(Conf))...)
# end
SW.flipPlaquette!(Conf, 14, 9)
# SW.flipPlaquette!(Conf,7,16)
# for i in 1:4:size(Conf,1)
#     SW.flipSpinsAlongRow!(Conf,i)
# end
# for i in 1:2:size(Conf,2)
#     SW.flipSpinsAlongCol!(Conf,i)
# end
# for i in -2size(Conf,1)+1:6:2size(Conf,1)
#     SW.flipSpinsAlongDiagonal!(Conf,i,-1)
# end
# SW.flipSpinsAlongRow!(Conf,12)
# SW.flipSpinsAlongCol!(Conf,11)
# Conf = SW.SpinConfig(rand(-0.5:0.5,20,20),1/2)
SW.plotApplPlaquettes(Conf)
SW.plotFractons!(current_axis(), Conf)
current_figure()
##
±(x) = minmax(x, -x)
±(x::AbstractArray) = ±(maximum(abs, x))
Charge = conv(Conf, iceRule)#[3:23,3:23]
@info "" numPosCharges(Charge)
# Charge = conv(SW.booleanSpinConfig(Conf),iceRule)

plotCharges(enforceConstraint!(Charge))

##
plotChargeFT(Charge)
##
# paddedKernel = LatticeFFTs.padSusc(iceRule, size(Charge,1))
paddedKernel = LatticeFFTs.padSusc(iceRule, size(Charge, 1))
kernelFFT = fft(paddedKernel)
Spin = deconvolve(Charge, kernelFFT, Mag(Conf), StaggeredMag(Conf))
SW.plotApplPlaquettes(Spin)
SW.plotFractons!(current_axis(), Spin)
current_figure()
##

paddedKernel_large = LatticeFFTs.padSusc(iceRule, size(Charge, 1))
kernelFFT_large = fft(paddedKernel_large)
##
with_theme(theme_PiTicks()) do
    fig = Figure(size = (1000, 400))
    axre = Axis(fig[1, 1], title = L"$\text{Re}$ $c(\mathbf{k})$", aspect = 1)
    axim = Axis(fig[1, 2], title = L"$\text{Im}$ $c(\mathbf{k})$", aspect = 1)
    axabs = Axis(fig[1, 3], title = L"|c(\mathbf{k})|", aspect = 1)
    k = getk(kernelFFT_large)

    realPart = real(kernelFFT_large)
    hm = heatmap!(axre, k, k, realPart, colormap = :balance, colorrange = ±(realPart))
    Colorbar(fig[2, 1], hm, vertical = false)

    imPart = imag(kernelFFT_large)
    hm = heatmap!(axim, k, k, imPart, colormap = :balance, colorrange = ±(imPart))
    Colorbar(fig[2, 2], hm, vertical = false)

    # hm = heatmap!(axabs,k,k,abs.(kernelFFT_large).+1e-30,colormap = :viridis,colorscale = identity)

    hm = heatmap!(
        axabs,
        k,
        k,
        abs.(kernelFFT_large) .+ 1e-20,
        colormap = :viridis,
        colorscale = log10,
    )
    contour!(axabs, k, k, abs.(kernelFFT_large), levels = [3e-2], color = :black)
    Colorbar(fig[2, 3], hm, vertical = false)
    fig
end

##
fig = Figure()
ax = Axis(fig[1, 1], title = L"$\text{Re}$ $c_i$", aspect = 1)
newkernelFFT = copy(kernelFFT)
newkernelFFT[findall(iszero, newkernelFFT)] .= 1e-10
invkernelFFT = 1 ./ newkernelFFT
InvKernel = FFTW.ifft(invkernelFFT)
hm = heatmap!(ax, abs.(InvKernel), colorscale = log10)
Colorbar(fig[1, 2], hm, vertical = true)
fig
##
Spin = conv(complex(Charge), InvKernel)
##

Charge = rand(-1:0.5:1, 20, 20)
for i in axes(Charge, 1), j in axes(Charge, 2)
    iseven(i + j) && (Charge[i, j] = 0)
    # Charge[i,j] = 0
end
# plotCharges(Charge)

paddedKernel = LatticeFFTs.padSusc(iceRule, size(Charge, 1))
kernelFFT = fft(paddedKernel)
##
Spin = deconvolve(Charge, kernelFFT, normalize = true)
SW.plotFractons(Spin)
##
kernelAnalytical(kx, ky) = 2(cos(kx) - cos(ky) + cos(kx + ky) - cos(kx - ky))
kernelAnalytical(k::SVector{2}) = kernelAnalytical(k[1], k[2])

rotMat(θ) = SA[
    cos(θ) -sin(θ)
    sin(θ) cos(θ)
]
kernelRotated(kx, ky) = kernelAnalytical(rotMat(pi / 4) * SA[kx, ky])

# rootAnalytical(kx) = 2atan((2sin(kx) + sqrt(5)* sqrt(sin(kx)^2))/(1 + cos(kx)))
# rootAnalytical1(kx) = 2atan((2 + sqrt(5))*tan(kx/2))
# rootAnalytical2(kx) = 2atan((2 - sqrt(5))*tan(kx/2))

rootAnalytical1(kx) = 2atan((2 + sqrt(5)) * cot(kx / 2)) + pi
rootAnalytical2(kx) = 2atan((2 - sqrt(5)) * cot(kx / 2)) + pi

with_theme(theme_PiTicks()) do
    k = LinRange(0, 2pi, 1000)
    mat = [kernelAnalytical(i, j) for i in k, j in k]
    melFilter(x) = max(abs(x), 1e-10)
    absmat = melFilter.(mat)
    fig = Figure(size = (650, 400))
    ax1 = Axis(fig[1, 1], title = L"$c(\mathbf{k})$", aspect = 1)
    ax2 = Axis(
        fig[1, 2],
        title = L"$|c(\mathbf{k})|$",
        aspect = 1,
        yticksvisible = false,
        yticklabelsvisible = false,
    )
    hm1 = heatmap!(ax1, k, k, mat, colormap = :balance)
    hm2 = heatmap!(ax2, k, k, absmat, colormap = :jet, colorscale = log10)

    lines!(ax2, k, rootAnalytical1.(k), color = :black, linestyle = :dash, linewidth = 1)
    lines!(ax2, k, rootAnalytical2.(k), color = :black, linestyle = :dash, linewidth = 1)
    Colorbar(fig[2, 1], hm1, vertical = false)
    Colorbar(fig[2, 2], hm2, vertical = false)
    fig

end
##
let
    zeroThing = zeros(1000, 1000)

    k = getk(zeroThing)
    fig = Figure()
    ax = Axis(fig[1, 1], title = L"$\text{Re}$ $c(\mathbf{k})$", aspect = 1)
    for (i, kx) in enumerate(k), (j, ky) in enumerate(k)
        if isapprox(rootAnalytical1(kx), ky, rtol = 1e-3, atol = 2e-2) ||
           isapprox(rootAnalytical2(kx), ky, rtol = 1e-3, atol = 2e-2)
            zeroThing[i, j] = 1
        end
    end
    hm = heatmap!(ax, k, k, zeroThing, colormap = :balance, colorrange = ±(zeroThing))
    fig

end
##
function sampleCharges(N)
    getCharge() = conv(
        SW.constructConfigPath(
            20,
            20,
            SW.ALLGS_S12,
            SW.spiralPath,
            maxiter = 100000,
            verbose = false,
        ),
        iceRule,
    )
    charges = fetch.([Threads.@spawn getCharge() for i = 1:N])
    filter!(x -> !any(isnan, x), charges)
    @info "" length(charges)
    chargeffts = fetch.([Threads.@spawn abs.(fft(charge)) for charge in charges])
end
function getAvgCharges(N)
    chargeffts = sampleCharges(N)
    avgCharge = mean(chargeffts)
    plotChargeFFT(avgCharge; colormap = :balance, colorrange = ±(avgCharge))
end
getAvgCharges(1000)
##
function enforceConstraint!(Charge)
    for i in axes(Charge, 1), j in axes(Charge, 2)
        iseven(i + j) && (Charge[i, j] = 0)
    end
    return Charge
end
function modifyFourier!(ChargeFFT)
    k = getk(ChargeFFT)
    for (i, kx) in enumerate(k), (j, ky) in enumerate(k)
        ChargeFFT[i, j] *= abs2(kernelAnalytical(kx, ky))
        # ChargeFFT[i,j] *= abs(kernelAnalytical(kx,ky))
    end
    ChargeFFT
end
function getStructureFac(ChargeFFT)
    function c(kx, ky)
        res = kernelAnalytical(kx, ky)
        if res ≈ 0
            return copysign(1e-10, res)
        else
            return res
        end
    end

    # CircFFTs = SW.CircularArrays.CircularArray.(ChargeFFT)
    k = getk(first(ChargeFFT))
    return [
        mean(C[i, j]^2.0 for C in ChargeFFT) / c(kx, ky)^2 for (i, kx) in enumerate(k),
        (j, ky) in enumerate(k)
    ]
    # return mean(abs.(C) for C in ChargeFFT)
end
##
Charges = mean(abs.(fft(enforceConstraint!(rand(-4:4, 50, 50)))) for i = 1:200)
##
plotChargeFFT(Charges, ; colormap = :viridis, colorrange = (0, maximum(Charges)))
##
Charges =
    real.(modifyFourier!.(fft(enforceConstraint!(rand(-4:4, 100, 100))) for i = 1:2000))

SF = getStructureFac(Charges)
plotChargeFFT(SF, colormap = :viridis, colorrange = extrema(SF))
##
ChargeFFTs = sampleCharges(10000)
plotChargeFFT(mean(ChargeFFTs))
##
SF = getStructureFac(ChargeFFTs)
plotChargeFFT(SF, colormap = :viridis, colorrange = extrema(SF))
##
Conf = SW.constructConfigPath(
    10,
    10,
    SW.ALLGS_S12,
    SW.xdirecPath,
    maxiter = 10000000,
    verbose = false,
)
# Conf = rand(-0.5:0.5,20,20)
##

Charge = conv(Conf, iceRule)

plotCharges(Charge; colorrange = (-4, 4))
