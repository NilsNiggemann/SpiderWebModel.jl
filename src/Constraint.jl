CONSTRAINT_SIGNS = (1, 1, -1, -1, 1, 1, -1, -1)
constraintSigns() = CONSTRAINT_SIGNS

function getSitesFromPlaquette(P::AbstractSpinConfig)
    return getSitesFromPlaquette(P.Mat)
end

function constraint(P::AbstractSpinConfig)
    return constraint(getSitesFromPlaquette(P))
end


function constraint(sites)
    cons = CONSTRAINT_SIGNS
    return sum(sgn * S for (sgn, S) in zip(cons, sites))
end

Base.@propagate_inbounds function getPlaquette(S::AbstractSpinConfig, i, j)
    Mat = @view parent(S)[(i-1):(i+1), (j-1):(j+1)]
    return SpinConfig(Mat, getSpin(S))
end

"""generates all possible combinations 0,1 of the 8 spins in a plaquette """

function getAllGS(S)
    allCombs = [
        plaquette((i, j, k, l, m, n, o, p), S) for i = (-S):S for j = (-S):S for
        k = (-S):S for l = (-S):S for m = (-S):S for n = (-S):S for o = (-S):S for
        p = (-S):S if constraint((i, j, k, l, m, n, o, p)) == 0
    ]
    return allCombs
end

function getAllGS_noMissing(S)
    allCombs = [
        plaquette((i, j, k, l, m, n, o, p, q), S) for i = (-S):S for j = (-S):S for
        k = (-S):S for l = (-S):S for m = (-S):S for n = (-S):S for o = (-S):S for
        p = (-S):S for q = (-S):S if constraint((i, j, k, l, m, n, o, p)) == 0
    ]
    return allCombs
end

const ALLGS_S12 = getAllGS(0.5)
const ALLGS_S12_NOMISSING = getAllGS_noMissing(0.5)

function plaquetteIsInBounds(Conf::AbstractMatrix, iCenter::Integer, jCenter::Integer)
    i = (iCenter-1):(iCenter+1)
    j = (jCenter-1):(jCenter+1)
    return checkbounds(Bool, Conf, i, j)
end

function plaquetteIsInBounds(Conf::AbstractSpinConfig, iCenter::Integer, jCenter::Integer)
    plaquetteIsInBounds(parent(Conf), iCenter, jCenter)
end

function allSpinsInBounds(Conf,Spin; verbose = false)
    for x in Conf
        isnan(x) && continue
        if abs(x) > Spin
            verbose && println("Spin larger than S")
            return false
        end
    end
    return true
end

function allSpinsInBounds(Conf::AbstractSpinConfig; verbose = false)
    allSpinsInBounds(parent(Conf.Mat), getSpin(Conf); verbose)
end

"""Assumes that constraint are only defined on every alternating site, starting from the first index"""
function fulFillsConstraint_nonStrict(Conf::AbstractSpinConfig, flipParity = false; verbose = false)
    for I in plaquetteIterator(Conf,!flipParity)
        i,j = Tuple(I)
        P = getPlaquette(Conf, i, j)
        any(isnan, P) && continue

        c = constraint(P)

        if c ≠ 0
            verbose && println("Constraint not fulfilled at i,j = $((i,j))")
            return false
        end
    end

    return true
end

function fulFillsConstraint(Conf::AbstractSpinConfig, flipParity = false; verbose = false)
    allSpinsInBounds(Conf; verbose) || return false
    return fulFillsConstraint_nonStrict(Conf, flipParity; verbose)
end
