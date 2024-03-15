const P1 = SA[
    -1.0 1.0 1.0
    -1.0 0.0 -1.0
    1.0 1.0 -1.0
]
const P2 = -P1

const P1_SITES = -SVector(getSitesFromPlaquette(P1)) / 2
const P2_SITES = -SVector(getSitesFromPlaquette(P2)) / 2

const P1_SITES_BOOL = P1_SITES .== 1 / 2
const P2_SITES_BOOL = P2_SITES .== 1 / 2

function plaquetteFlippable(plaqSites::SVector{8})
    return plaqSites == P1_SITES || plaqSites == P2_SITES
end

function plaquetteFlippable(plaqSites::SVector{8,Union{Integer,Bool}})
    return plaqSites == P1_SITES_BOOL || plaqSites == P2_SITES_BOOL
end

function flipPlaquette!(Conf::AbstractMatrix, i, j)
    P = getPlaquette(Conf, i, j)

    P[1, 1] = -P[1, 1]
    P[1, 2] = -P[1, 2]
    P[1, 3] = -P[1, 3]
    P[2, 1] = -P[2, 1]
    # P[2,2] = -P[2,2]
    P[2, 3] = -P[2, 3]
    P[3, 1] = -P[3, 1]
    P[3, 2] = -P[3, 2]
    P[3, 3] = -P[3, 3]

    return Conf
end

function flipPlaquette!(Conf::AbstractMatrix{Bool}, i, j)
    P = getPlaquette(Conf, i, j)

    P[1, 1] = !P[1, 1]
    P[1, 2] = !P[1, 2]
    P[1, 3] = !P[1, 3]
    P[2, 1] = !P[2, 1]
    # P[2,2] = !P[2,2]
    P[2, 3] = !P[2, 3]
    P[3, 1] = !P[3, 1]
    P[3, 2] = !P[3, 2]
    P[3, 3] = !P[3, 3]

    return Conf
end

function flipPlaquette!(Conf::AbstractMatrix, pos::Tuple)
    i, j = pos
    flipPlaquette!(Conf, i, j)
end

function flipPlaquette!(Conf::AbstractMatrix, pos::CartesianIndex{2})
    i, j = Tuple(pos)
    flipPlaquette!(Conf, i, j)
end

function flipPlaquette!(Conf::AbstractMatrix, pos::Int)
    ij = CartesianIndices(Conf)[pos]
    flipPlaquette!(Conf, ij)
end

function sortByValueOrder(D)
    ks = collect(keys(D))
    vals = collect(values(D))
    return ks[sortperm(vals)]
end

function sortbyKeyOrder(D)
    ks = collect(keys(D))
    vals = collect(values(D))
    return vals[sortperm(ks)]
end

function invertDict(D)
    D2 = Dict(v => k for (k, v) in D)
end


function flipSpinsAlongLine!(Conf, org, slope)
    slope ∈ (-Inf, Inf) && return flipSpinsAlongRow!(Conf, org[2])
    for i in axes(Conf.Mat, 1), j in axes(Conf.Mat, 2)
        if slope * (i - org[1]) == j - org[2]
            Conf[i, j] *= -1
        end
    end
    return Conf
end

function flipSpinsAlongDiagonal!(Conf, org, slope)
    j = org
    for i in axes(Conf.Mat, 1)
        j += slope
        if checkbounds(Bool, Conf, i, j)
            Conf[i, j] *= -1
        end
    end
    return Conf
end

function flipSpinsAlongRow!(Conf, i, offset = 0)
    Conf[1+offset:2:end, i] .*= -1
    return Conf
end

function flipSpinsAlongCol!(Conf, i, offset = 0)
    Conf[i, 1+offset:2:end] .*= -1
    return Conf
end

function CanApply(Conf::SpinConfig, Op::AbstractMatrix, i, j)
    plaquetteIsInBounds(Conf, i, j) || return false
    P = getPlaquette(Conf, i, j)

    P .+= Op

    applicable = fulFillsConstraint(Conf)

    P .-= Op

    return applicable
end

function canFlipPlaquette(Conf::SpinConfig, i, j)
    isodd(i + j) || return false
    plaquetteIsInBounds(Conf, i, j) || return false

    # P = getPlaquette(Conf,i,j)
    # sites = SVector(Mat[2,3],Mat[1,3],Mat[1,2],Mat[1,1],Mat[2,1],Mat[3,1],Mat[3,2],Mat[3,3])
    sites = SVector(
        Conf[i, j+1],
        Conf[i-1, j+1],
        Conf[i-1, j],
        Conf[i-1, j-1],
        Conf[i, j-1],
        Conf[i+1, j-1],
        Conf[i+1, j],
        Conf[i+1, j+1],
    )
    return plaquetteFlippable(sites)
end

function CanApplyNonStrict(Conf::SpinConfig, Op, i, j)
    isodd(i + j) || return false
    plaquetteIsInBounds(Conf, i, j) || return false
    @inbounds P = getPlaquette(Conf, i, j)

    OpSites = getSitesFromPlaquette(Op)

    sites = SVector(
        Conf[i, j+1],
        Conf[i-1, j+1],
        Conf[i-1, j],
        Conf[i-1, j-1],
        Conf[i, j-1],
        Conf[i+1, j-1],
        Conf[i+1, j],
        Conf[i+1, j+1],
    )

    Psites = sites .+ OpSites

    applicable = all(x -> -Conf.S <= x <= Conf.S, Psites)

    return applicable
end

CanApplyNonStrict(Conf::SpinConfig, ::Nothing, i, j) = canFlipPlaquette(Conf, i, j)

function CanApplyAnywhere(Conf::SpinConfig, Op::AbstractMatrix)
    a1 = axes(Conf.Mat, 1)
    a2 = axes(Conf.Mat, 2)

    Opx, Opy = size(Op)
    for i in a1, j in a2
        firstindex(a1) + Opx <= i <= lastindex(a1) - Opx || continue
        firstindex(a2) + Opy <= j <= lastindex(a2) - Opy || continue

        if CanApply(Conf, Op, i, j)
            return true
        end
    end
    return false
end

"""assumes that Op is already an allowed operator"""
function getApplicablePlaquettes(Conf::SpinConfig, Op)
    plaqPos = [
        (i, j) for i in axes(Conf.Mat, 1) for
        j in axes(Conf.Mat, 2) if CanApplyNonStrict(Conf, Op, i, j)
    ]
    return plaqPos
end

function getApplicablePlaquettes(Conf::SpinConfig)
    plaqPos = [
        (i, j) for i in axes(Conf.Mat, 1) for
        j in axes(Conf.Mat, 2) if canFlipPlaquette(Conf, i, j)
    ]
    return plaqPos
end
getApplicablePlaquettes(Conf::SpinConfig,::Nothing) = getApplicablePlaquettes(Conf)