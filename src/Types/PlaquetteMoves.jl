const P1 = SA[
    -1.0 1.0 1.0
    -1.0 0.0 -1.0
    1.0 1.0 -1.0
]
const P2 = -P1
function getSitesFromPlaquette(Mat::AbstractMatrix)
    return (
        Mat[2, 3],
        Mat[1, 3],
        Mat[1, 2],
        Mat[1, 1],
        Mat[2, 1],
        Mat[3, 1],
        Mat[3, 2],
        Mat[3, 3],
    )
end
const P1_SITES = -SVector(getSitesFromPlaquette(P1)) / 2
const P2_SITES = -SVector(getSitesFromPlaquette(P2)) / 2

const P1_SITES_BOOL = P1_SITES .== 1 / 2
const P2_SITES_BOOL = P2_SITES .== 1 / 2

"""Order in which the spins are stored in the plaquette corresponding to the numbering convention"""
plaquetteOrder(S1, S2, S3, S4, S5, S6, S7, S8, S9) =
    SMatrix{3,3}(S4, S5, S6, S3, S9, S7, S2, S1, S8)

function plaquetteOrder(S1, S2, S3, S4, S5, S6, S7, S8)
    SMatrix{3,3}(S4, S5, S6, S3, NaN, S7, S2, S1, S8)
end
plaquetteOrder(x) = plaquetteOrder(x...)

function plaquette(sites, S = 1 / 2)
    Mat = SMatrix{3,3}(plaquetteOrder(sites))
    return SpinConfig(Mat, S)
end

