module S

using StaticArrays
abstract type AbstractSymop <: Function end

struct SymmetryTransformation{Dim,T<:Real,Dim²} <: AbstractSymop
    WMatrix::SMatrix{Dim,Dim,T,Dim²}
end

getTranslation(S::SymmetryTransformation) = S.WMatrix[1:2,3]
getRotation(S::SymmetryTransformation) = SMatrix{2,2}(S.WMatrix[1:2,1:2])


SymmetryTransformation(t::SVector{2,T},W::SMatrix{2,2,T,4}) where T <: Real = SymmetryTransformation(SA[
    W[1,1] W[1,2] t[1];
    W[2,1] W[2,2] t[2];
    0 0 1
    ]
) 


import Base.*
*(S::SymmetryTransformation,S2::SymmetryTransformation) = SymmetryTransformation(S.WMatrix*S2.WMatrix)

import Base.^
^(S::SymmetryTransformation,i::Number) = SymmetryTransformation(S.WMatrix^i)

"""given an origin and a transformation matrix, returns a SymmetryTransformation"""
function Symmetrytransformation_origin(origin::AbstractVector{T},matrix::AbstractMatrix{T}) where T <: Real
    return SymmetryTransformation(-matrix*origin +origin,matrix)
end

# getOrigin(S::SymmetryTransformation) = inv(getRotation(S)+I)*getTranslation(S)

function (S::SymmetryTransformation{3})(vec::T) where T
    x,y = vec
    x,y,_ = S.WMatrix * SA[x,y,one(x)]
    return convert(T,(x,y))
end

function (S::SymmetryTransformation{3})(I::CartesianIndex{2})
    return CartesianIndex(S(Tuple(I)))
end

parseAsSMatrix(v::AbstractVector{<:AbstractString}) = SMatrix{2,2,Float64,4}([parse(Float64,s) for s in v])'

function Base.show(io::IO, x::SymmetryTransformation{D,T}) where {D,T}
    println(io,"$D dim SiteTransformation{",T,"}")
    Base.print_matrix(io,x.WMatrix)
end
Base.show(io::IO, ::MIME"text/plain", x::SymmetryTransformation) = show(io,x)

ispuretranslation(S::SymmetryTransformation) = abs(det(getRotation(S))) == 0

function Translation(x::Real,y::Real)
    return SymmetryTransformation(SA[
        1 0 x;
        0 1 y;
        0 0 1
        ]
    )
end

function Rotation(θ::Real)
    return SymmetryTransformation(SA[
        cos(θ) -sin(θ) 0;
        sin(θ) cos(θ) 0;
        0 0 1
        ]
    )
end


function generateSymms(irreps::AbstractVector{<:SymmetryTransformation};digits = 14,maxiter = 1000)
    syms = Set(irreps)
    roundsyms = Set(copy(irreps)) # rounded version to avoid duplicates

    iteration = 0
    for s1 in syms
        for s2 in syms
            if iteration > maxiter
                @warn "Too many iterations"
                return collect(syms)
            end
            Snew = s1*s2
            
            isInUnitCell(getTranslation(Snew)) || continue
            roundSnew = round(Snew; digits)
            if roundSnew ∉ roundsyms
                iteration += 1
                push!(syms,Snew)
                push!(roundsyms,roundSnew)
            end
        end
    end
    return collect(syms)
end

function isInUnitCell(v::AbstractVector{T}) where T
    return all(0 .<= v .< 1)
end

function Base.round(s::SymmetryTransformation;kwargs...)
    T = round.(s.WMatrix;kwargs...)
    SymmetryTransformation(T)
end

function translate_back_to_UC(s::SymmetryTransformation) 
    org = getTranslation(s)
    rot = getRotation(s)
    org = org - floor.(org)
    return SymmetryTransformation(org,rot)
end

end

##

a1 = S.Translation(1,1)
a2 = S.Translation(1,-1)

b = round(S.Rotation(π/2))

syms = S.generateSymms([a1,a2,b])

##
x = rand(12,12)

sites = collect(CartesianIndices(x))
##
syms[1](sites[1])