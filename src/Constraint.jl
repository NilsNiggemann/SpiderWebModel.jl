CONSTRAINT_SIGNS = (1,1,-1,-1,1,1,-1,-1)
constraintSigns() = CONSTRAINT_SIGNS

function constraint(sites)
    cons = CONSTRAINT_SIGNS
    return sum(sgn*S for (sgn,S) in zip(cons,sites))
end

"""Order in which the spins are stored in the plaquette corresponding to the numbering convention"""
plaquetteOrder(S1,S2,S3,S4,S5,S6,S7,S8,S9) = SMatrix{3,3}(S4,S5,S6,S3,S9,S7,S2,S1,S8)

plaquetteOrder(S1,S2,S3,S4,S5,S6,S7,S8) = SMatrix{3,3}(S4,S5,S6,S3,NaN,S7,S2,S1,S8)
plaquetteOrder(x) = plaquetteOrder(x...)

function plaquette(sites,S=1/2)
    Mat = SMatrix{3,3}(plaquetteOrder(sites))
    return SpinConfig(Mat,S)
end

function getSitesFromPlaquette(Mat::AbstractMatrix)
    return (Mat[2,3],Mat[1,3],Mat[1,2],Mat[1,1],Mat[2,1],Mat[3,1],Mat[3,2],Mat[3,3])
end

function getSitesFromPlaquette(P::SpinConfig) 
    return getSitesFromPlaquette(P.Mat)
end

function constraint(P::SpinConfig)
    return constraint(getSitesFromPlaquette(P))
end

function getPlaquette(S::SpinConfig,i,j)
    Mat = @view S.Mat[i-1:i+1,j-1:j+1]
    return SpinConfig(Mat,S.S)
end

"""generates all possible combinations 0,1 of the 8 spins in a plaquette """

function getAllGS(S)
    allCombs = [plaquette((i,j,k,l,m,n,o,p),S) for i in -S:S for j in -S:S for k in -S:S for l in -S:S for m in -S:S for n in -S:S for o in -S:S for p in -S:S if constraint((i,j,k,l,m,n,o,p)) == 0]
    return allCombs
end

function getAllGS_noMissing(S)
    allCombs = [plaquette((i,j,k,l,m,n,o,p,q),S) for i in -S:S for j in -S:S for k in -S:S for l in -S:S for m in -S:S for n in -S:S for o in -S:S for p in -S:S for q in -S:S if constraint((i,j,k,l,m,n,o,p)) == 0]
    return allCombs
end


const ALLGS_S12 = getAllGS(0.5)
const ALLGS_S12_NOMISSING = getAllGS_noMissing(0.5)

function plaquetteIsInBounds(Conf::AbstractMatrix,iCenter::Integer,jCenter::Integer)
    i = iCenter-1:iCenter+1
    j = jCenter-1:jCenter+1
    return checkbounds(Bool,Conf,i,j)
end

plaquetteIsInBounds(Conf::SpinConfig,iCenter::Integer,jCenter::Integer) = plaquetteIsInBounds(Conf.Mat,iCenter,jCenter)

function allSpinsInBounds(Conf::SpinConfig;verbose=false)
    for x in Conf.Mat
        isnan(x) && continue
        if abs(x) > Conf.S 
            verbose && println("Spin larger than S")
            return false
        end
    end
    return true
end

"""Assumes that constraint are only defined on every alternating site, starting from the first index"""
function fulFillsConstraint_nonStrict(Conf::SpinConfig,flipParity=false;verbose = false)

    for i in axes(Conf.Mat,1), j in axes(Conf.Mat,2)
        iseven(i+j+flipParity) || continue 
        plaquetteIsInBounds(Conf,i,j) || continue
        P = getPlaquette(Conf,i,j)
        any(isnan,P) && continue

        c = constraint(P)

        if c ≠ 0
            verbose && println("Constraint not fulfilled at i,j = $((i,j))" )
            return false
        end
    end
    
    return true
end

function fulFillsConstraint(Conf::SpinConfig,flipParity=false;verbose = false)
    allSpinsInBounds(Conf;verbose) || return false
    return fulFillsConstraint_nonStrict(Conf,flipParity;verbose)
end

