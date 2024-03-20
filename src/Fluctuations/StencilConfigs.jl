function stencilConfig(A::AbstractArray,S::T;kwargs...) where T
    St = Stencils.StencilArray(A,Stencils.Moore();boundary = Stencils.Remove(T(typemax(T))),kwargs...)
    return SpinConfig(St,S)
end

stencilConfig(A::SpinConfig;kwargs...) = stencilConfig(parent(A),A.S;kwargs...)


@inline Base.@propagate_inbounds getPlaquetteStencil(S::Stencils.StencilArray,i::Int,j::Int) = Stencils.stencil(S,CartesianIndex(i,j))

@inline Base.@propagate_inbounds getPlaquetteStencil(S::SpinConfig,i::Int,j::Int) = getPlaquetteStencil(parent(S),i,j)

@inline Base.@propagate_inbounds function getSitesFromStencil(S::Stencils.Moore)
    S6,S7,S8,S5,S1,S4,S3,S2 = S
    return SVector(S1,S2,S3,S4,S5,S6,S7,S8)
end

@inline Base.@propagate_inbounds function getPlaquetteSites(S::Stencils.StencilArray,i::Int,j::Int)
    S = getPlaquetteStencil(S,i,j)
    return getSitesFromStencil(S)
end

@inline Base.@propagate_inbounds getPlaquetteSites(S::SpinConfig,i::Int,j::Int) = getPlaquetteSites(parent(S),i,j)

"""returns a tuple of two booleans, corresponding to the applicability of P and P⁺. M=2S is an integer"""  
@inline function P_applicable(Pij::Union{SVector{8},Stencils.Moore},M::Integer)
    S1, S2, S3, S4, S5, S6, S7, S8 = Pij
    v = (S1,-S2,-S3,S4,S5,-S6,-S7,S8)

    p = all(<(M),v)
    pdagger = all(>(-M),v)
    return p ,pdagger
end

@inline function P_applicable(S::SpinConfig, i::Int, j::Int)
    Pij = getPlaquetteSites(S, i, j)
    P_applicable(Pij, S.S)
end

@inline P_applicable(S::SpinConfig, I::Union{<:NTuple{2},<:CartesianIndex{2}}) = P_applicable(S, I[1], I[2])

function getApplicablePlaquettes!(plaqsPlus,plaqsMinus,Conf::SpinConfig{T,<:Stencils.StencilArray}) where T
    empty!(plaqsPlus)
    empty!(plaqsMinus)
    for I in plaquetteIterator(Conf)
        applPlus,applMinus = P_applicable(Conf, I)
        if applPlus
            push!(plaqsPlus, Tuple(I))
        end
        if applMinus
            push!(plaqsMinus, Tuple(I))
        end
    end
    return plaqsPlus,plaqsMinus
end
