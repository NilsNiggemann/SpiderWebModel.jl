function stencilConfig(A::AbstractArray,S::Real;kwargs...)
    St = Stencils.StencilArray(A,Stencils.Moore();boundary = Stencils.Remove(NaN),kwargs...)
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


@inline function P_plus_applicable(Pij::Union{SVector{8},Stencils.Moore},Spin::Real)
    S1, S2, S3, S4, S5, S6, S7, S8 = Pij
    max(S1,-S2,-S3,S4,S5,-S6,-S7,S8) <= Spin-1
end

@inline function P_minus_applicable(Pij::Union{SVector{8},Stencils.Moore},Spin::Real)
    S1, S2, S3, S4, S5, S6, S7, S8 = Pij
    max(-S1,S2,S3,-S4,-S5,S6,S7,-S8) <= Spin-1
end

@inline function P_plus_applicable(S::SpinConfig, i::Int, j::Int)
    Pij = getPlaquetteSites(S, i, j)
    P_plus_applicable(Pij, S.S)
end
@inline function P_minus_applicable(S::SpinConfig, i::Int, j::Int)
    Pij = getPlaquetteSites(S, i, j)
    P_minus_applicable(Pij, S.S)
end
@inline P_plus_applicable(S::SpinConfig, I::Union{<:NTuple{2},<:CartesianIndex{2}}) = P_plus_applicable(S, I[1], I[2])
@inline P_minus_applicable(S::SpinConfig, I::Union{<:NTuple{2},<:CartesianIndex{2}}) = P_plus_applicable(S, I[1], I[2])

function getApplicablePlaquettes!(plaqsPos,Conf::SpinConfig{T,<:Stencils.StencilArray},SignFlip::Bool) where T
    empty!(plaqsPos)
    if SignFlip
        for I in plaquetteIterator(Conf)
            if P_plus_applicable(Conf, I)
                push!(plaqsPos, Tuple(I))
            end
        end
    else
        for I in plaquetteIterator(Conf)
            if P_minus_applicable(Conf,I)
                push!(plaqsPos, Tuple(I))
            end
        end
    end
    return plaqsPos
end

function getApplicablePlaquettes2!(plaqsPos,Conf::SpinConfig{T,<:Stencils.StencilArray},SignFlip::Bool) where T
    if SignFlip
        Stencils.mapstencil!(x->P_plus_applicable(x,Conf.S),plaqsPos,parent(Conf))
    else
        Stencils.mapstencil!(x->P_minus_applicable(x,Conf.S),plaqsPos,parent(Conf))
    end
    return plaqsPos
end

getApplicablePlaquettes(Conf::SpinConfig{T,<:Stencils.StencilArray},SignFlip::Bool) where T = getApplicablePlaquettes!(Vector{Tuple{Int,Int}}(),Conf,SignFlip)
getApplicablePlaquettes2(Conf::SpinConfig{T,<:Stencils.StencilArray},SignFlip::Bool) where T = getApplicablePlaquettes2!(deepcopy(parent(Conf)),Conf,SignFlip)