"""given a matrix, rotate it by 90 degrees"""
function rotate(Mat)
    Mat2 = zeros(size(Mat))
    for i in axes(Mat, 1)
        row = Mat'[:, i]
        Mat2[:, i] .= reverse(row)
    end
    return Mat2
end
function rotate(Mat, n)
    if n == 0
        return Mat
    elseif n == 1
        return rotate(Mat)
    elseif n == 2
        return Mat |> rotate |> rotate
    elseif n == 3
        return Mat |> rotate |> rotate |> rotate
    end
    error("n must be 0,1,2,3")
end
rotate(S::SW.SpinConfig, n) = SW.SpinConfig(rotate(S.Mat, n), S.S)

function AllRots(Confs)
    Confs2 = [rotate(c, n) for c in Confs for n = 0:3]
    return Confs2
end

function randRots(Confs)
    Confs2 = [rotate(c, rand(0:3)) for c in Confs]
    return Confs2
end

##
using LatticeFFTs
using MakieHelpers
using Statistics
using LatticeFFTs.Interpolations
##
function getStructureFac(Confs)
    plan = getLatticeFFTPlan(Confs[1].Mat, 0)

    Sq = [
        getInterpolatedFFT(c.Mat, 0, plan; Interpolation = BSpline(Constant())) for
        c in Confs
    ]
end
##
function plotStructureFac(Confs; cbar = false, kwargs...)
    Confs2 = Confs
    # Confs2 = randRots(Confs)

    Sq = getStructureFac(Confs2)
    SSq(kx, ky) = mean(real(s(kx, ky) * s(-kx, -ky)) for s in Sq)
    # SSq(kx,ky) = mean(real(s(kx,ky)) for s in Sq)

    # i = size(Confs2[1],1) ÷ 2
    # j = size(Confs2[1],2) ÷ 2

    # ij = LinearIndices(Confs2[1])[i,j]

    # Sij = SW.getSij(Confs2,ij)[:]
    # Sij2 = SW.getSij(Confs2,ij+1)[:]
    # append!(Sij,Sij2)
    # Rij = SW.getRij_vec(Confs2[1],ij)[:]
    # Rij2 = SW.getRij_vec(Confs2[1],ij+1)[:]
    # append!(Rij,Rij2)
    # chik = FourierStruct(Sij,Rij,1)

    kx = LinRange(-0, 2pi, 200)
    ky = LinRange(-0, 2pi, 200)
    chi = fetch.([Threads.@spawn SSq(kx, ky) for kx in kx, ky in ky])
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1, 1], aspect = 1, xticks = PiTicks(), yticks = PiTicks())
    hm = heatmap!(ax, kx, ky, chi; kwargs...)
    if cbar
        Colorbar(fig[1, 2], hm)
    end
    fig
end
