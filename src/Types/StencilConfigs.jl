struct StencilSpinConfig{MVal<:Val,T,MatType<:AbstractMatrix{T}} <: AbstractSpinConfig{T}
    Mat::MatType
    M::MVal
end

@inline Base.parent(S::StencilSpinConfig) = S.Mat
@inline getSpin(::StencilSpinConfig{Val{MVal}}) where {MVal} = MVal/2
@inline Base.similar(S::StencilSpinConfig,dims::Tuple) = StencilSpinConfig(similar(parent(S),dims), S.M)
@inline Base.similar(S::StencilSpinConfig,dims::Vararg{Int64, N}) where { N} = StencilSpinConfig(similar(parent(S),dims...), S.M)
# @inline Base.similar(S::StencilSpinConfig,dims...) = StencilSpinConfig(similar(parent(S),dims...), S.M)
# @inline Base.copy(S::StencilSpinConfig) = StencilSpinConfig(copy(parent(S)), S.M)
"""
create a StencilSpinConfig with Spin magnitude S, initialized by an array of Floats i.e.

    [
        0.5 -0.5 0.5;
        -0.5 0.5 -0.5;
        0.5 -0.5 0.5
    ] 

for efficient Monte Carlo simulations.
Supported boundary conditions are 

    :open # Open boundary conditions: No plaquette at the boundary may be flipped.
    :open_soft # Open boundary conditions: Ignore out-of bounds spins so that plaquettes at the boundary may be flipped. 
    :periodic # Periodic boundary conditions: the StencilArray is wrapped around the edges.
Boundary properties of the StencilArray may also be passed directly as kwargs, which will override defaults.
"""
function stencilConfig(A::AbstractMatrix{<:AbstractFloat}, S; kwargs...)
    @assert all(in(-S:S), A) "Spin configuration must be initialized within -S:S"
    return stencilConfig(Int8.(2 .*A), S;kwargs...)
end
function parseBoundaryCondition(boundaryCondition)
    paddingValue = Int8(typemax(Int8))
    if boundaryCondition == :open
        boundary = Stencils.Remove(paddingValue)
        padding = Stencils.Conditional()
    elseif boundaryCondition == :open_soft
        boundary = Stencils.Remove(Int8(0))
        padding = Stencils.Halo(:out)
    elseif boundaryCondition == :periodic
        boundary = Stencils.Wrap()
        padding = Stencils.Conditional()
    else
        error("boundaryCondition $boundaryCondition not implemented!")
    end
    return paddingValue,boundary,padding
end

function construct_stencilConfig(A::AbstractMatrix{Int8}, S,boundary = Stencils.Remove(paddingValue),paddingValue = Int8(typemax(Int8));kwargs...)
    St = Stencils.StencilArray(
        A,
        Stencils.Moore();
        boundary,
        kwargs...,
    )
    M = Val(Int(2S)) # do not chose Int8 to avoid overflows
    return StencilSpinConfig(St, M)
end
function stencilConfig(A::AbstractMatrix{Int8},S;boundaryCondition=:open,kwargs...)
    @assert all(in(-2S:2S), A) "Spin configuration must be initialized within -S:S"
    paddingValue,boundary,padding = parseBoundaryCondition(boundaryCondition)
    return construct_stencilConfig(A,S,boundary,paddingValue;padding,kwargs...)
end
stencilConfig(A::SpinConfig; kwargs...) = stencilConfig(parent(A), A.S; kwargs...)

function SpinConfig(S::StencilSpinConfig)
    boundary = Stencils.boundary(parent(S))
    return SpinConfig(S, boundary)
end
function SpinConfig(S::StencilSpinConfig,::Stencils.Wrap)
    return SpinConfig(PeriodicMatrix(parent(S)./2), getSpin(S))
end

function SpinConfig(S::StencilSpinConfig,::Stencils.Any)
    return SpinConfig(parent(S)./2, getSpin(S))
end

allSpinsInBounds(S::StencilSpinConfig{Val{M}}; kwargs...) where M = allSpinsInBounds(S,M; kwargs...)

plotSpinConfig!(ax, S::StencilSpinConfig; kwargs...) = plotSpinConfig!(ax, SpinConfig(S); kwargs...)
plotSpinConfig(S::StencilSpinConfig; kwargs...) = plotSpinConfig(SpinConfig(S); kwargs...)

function plaquetteIterator(S::Stencils.StencilArray,shift=false)
    bound = Stencils.boundary(S)
    padding = Stencils.padding(S)
    return stencil_plaquetteIterator(S,bound,padding,shift)
end

function stencil_plaquetteIterator(S,::Stencils.Wrap,::Stencils.Conditional,shift=false)
    inboundsInds = Base.Iterators.product(axes(S, 1), axes(S, 2))
    filterInds = Iterators.filter(ind -> isodd(sum(ind)+shift), inboundsInds)
end
function stencil_plaquetteIterator(S,::Stencils.Remove,::Stencils.Halo,shift=false)
    inboundsInds = Base.Iterators.product(axes(S, 1), axes(S, 2))
    filterInds = Iterators.filter(ind -> isodd(sum(ind)+shift), inboundsInds)
end
stencil_plaquetteIterator(S,::Any,::Any,shift=false) = plaquetteIterator(parent(S),shift)


@inline Base.@propagate_inbounds getPlaquetteStencil(
    S::Stencils.StencilArray,
    i::Int,
    j::Int,
) = Stencils.stencil(S, CartesianIndex(i, j))

@inline Base.@propagate_inbounds getPlaquetteStencil(S::StencilSpinConfig, i::Int, j::Int) =
    getPlaquetteStencil(parent(S), i, j)

@inline Base.@propagate_inbounds function getSitesFromStencil(S::Stencils.Moore)
    S6, S7, S8, S5, S1, S4, S3, S2 = S
    return SVector(S1, S2, S3, S4, S5, S6, S7, S8)
end

const P1_STENCIL = Int8.(SVector(-2, -2, 2, 2, 2, 2, -2, -2))

@inline Base.@propagate_inbounds function getPlaquetteSites(
    S::Stencils.StencilArray,
    i::Int,
    j::Int,
)
    S = getPlaquetteStencil(S, i, j)
    return getSitesFromStencil(S)
end

@inline getPlaquette(S::StencilSpinConfig, i, j) = getPlaquetteSites(S, i, j)

@inline Base.@propagate_inbounds getPlaquetteSites(S::StencilSpinConfig, i::Int, j::Int) =
    getPlaquetteSites(parent(S), i, j)

inverse_move(move) = (move[1], move[2], -move[3])

Base.@propagate_inbounds function applyPlaquette!(Config,i,j,sgn)
    applyPlaquette!(parent(Config), i, j, sgn)
    return Config
end
Base.@propagate_inbounds @inline applyPlaquette!(Config,move) = applyPlaquette!(Config,move[1], move[2], move[3])
Base.@propagate_inbounds function applyPlaquette!(Mat::Stencils.AbstractStencilArray,i,j,sgn)
    # sites = Stencils.indices(Stencils.stencil(Mat), CartesianIndex(i, j))
    sgntype = convert(eltype(Mat), sgn)
    sites = safe_parent_indices(Mat, (i, j))
    # sites = Stencils.indices(Mat, (i, j))
    # for (ij,s) in zip(sites,P1_STENCIL)
    for idx in safe_iterate_sites(Mat, (i, j))
        i,j = sites[idx]
        s = P1_STENCIL[idx]
        Mat[i,j] += sgntype *s
    end
    return Mat
end
@inline function safe_parent_indices_linear(A::Stencils.AbstractStencilArray,I,::Stencils.Conditional)
    inds = safe_parent_indices(A, CartesianIndex(I))

    LI = LinearIndices(A)
    linear_inds = map(inds) do (i, j)
        LI[i, j]
    end

    return linear_inds
end

@inline function safe_parent_indices_linear(A::Stencils.AbstractStencilArray,I,::Stencils.Halo)
    inds = safe_parent_indices(A, CartesianIndex(I))

    LI = LinearIndices(A)
    linear_inds = map(inds) do (i, j)
        checkbounds(Bool,LI,i,j) ? LI[i, j] : 0
    end

    return linear_inds
end
@inline safe_parent_indices_linear(A::Stencils.AbstractStencilArray,I) = safe_parent_indices_linear(A, I, Stencils.padding(A))
@inline safe_parent_indices_linear(A::StencilSpinConfig,I) = safe_parent_indices_linear(parent(A), I)

getStencilRadii(::Stencils.AbstractStencilArray{<:Any,R,<:Any,N}) where {R,N} = CartesianIndex(ntuple(_ -> -R, N))

@inline safe_parent_indices(S::StencilSpinConfig,I) = safe_parent_indices(parent(S), I)
@inline safe_parent_indices(A::Stencils.AbstractStencilArray,I) = Stencils.indices(A, I)

@inline safe_iterate_sites(S::StencilSpinConfig,I) = safe_iterate_sites(parent(S), I)
@inline function safe_iterate_sites(A::Stencils.AbstractStencilArray,I)
    return safe_iterate_sites(A, Stencils.boundary(A), Stencils.padding(A), CartesianIndex(I))
end
@inline function safe_iterate_sites(A::Stencils.AbstractStencilArray,::Any,::Stencils.Conditional,I)
    inds = eachindex(Stencils.stencil(A))
end
@inline function safe_iterate_sites(A::Stencils.AbstractStencilArray,::Any,::Stencils.Halo,I)
    site_inds = Stencils.indices(A, I)
    
    # inds = Iterators.filter(i->checkbounds(Bool,A,CartesianIndex(site_inds[i])),eachindex(site_inds))
    inds = SmallCollections.SmallVector{8,Int}(i for (i,I) in enumerate(site_inds) if checkbounds(Bool,A,CartesianIndex(I)))
end

# `Conditional` needs handling for specific boundary conditions.
# For Wrap we swap the side.
function get_wrappend_inds(A::Stencils.AbstractStencilArray{S,R}, I::Tuple) where {S,R}
    sz = Stencils.tuple_contents(S)
    wrapped_inds = map(I, sz) do i, s
        i < 1 ? i + s : (i > s ? i - s : i)
    end
    return wrapped_inds
end


"""returns a tuple of two booleans, corresponding to the applicability of P and P⁺. M=2S is an integer"""
@inline Base.@propagate_inbounds function P_applicable(Pij::Union{SVector{8},Stencils.Moore}, M::Integer)
    S1, S2, S3, S4, S5, S6, S7, S8 = Pij
    v = (S1, -S2, -S3, S4, S5, -S6, -S7, S8)

    p = all(<(M), v)
    pdagger = all(>(-M), v)
    return p, pdagger
end

"""checks whether F or F† can be applied. Pij are 2S_i in counterclockwise order. M = 2S is twice the spin."""
@inline Base.@propagate_inbounds function P_applicable(Pij::Union{SVector{8},Stencils.Moore}, ::Val{M}) where {M}
    # more performant than above since (S1 - 1) * (S2 + 1) ... is zero only if S1 == 1 and S2 == -1 ...
    # here we technically compute e.g. (2S1 - 2) * (2S2 + 2) ... with the same result
    S1, S2, S3, S4, S5, S6, S7, S8 = Pij
    p = !iszero((S1-M) * (S2+M) * (S3+M) * (S4-M) * (S5-M) * (S6+M) * (S7+M) * (S8-M))
    pdag = !iszero((S1+M) * (S2-M) * (S3-M) * (S4+M) * (S5+M) * (S6-M) * (S7-M) * (S8+M))
    return p, pdag
end

@inline Base.@propagate_inbounds function P_applicable(S::StencilSpinConfig, i::Int, j::Int)
    Pij = getPlaquetteSites(S, i, j)
    P_applicable(Pij, S.M)
end

@inline Base.@propagate_inbounds P_applicable(S::StencilSpinConfig, I::Union{<:NTuple{2},<:CartesianIndex{2}}) =
    P_applicable(S, I[1], I[2])

function getApplicablePlaquettes_separated!(
    plaqsPlus,
    plaqsMinus,
    Conf::StencilSpinConfig,
)
    empty!(plaqsPlus)
    empty!(plaqsMinus)
    for I in plaquetteIterator(Conf)
        applPlus, applMinus = P_applicable(Conf, I)
        if applPlus
            push!(plaqsPlus, Tuple(I))
        end
        if applMinus
            push!(plaqsMinus, Tuple(I))
        end
    end
    return plaqsPlus, plaqsMinus
end

function getApplicablePlaquettes!(
    plaqs,
    Conf::StencilSpinConfig,
)
    empty!(plaqs)
    for I in plaquetteIterator(Conf)
        applPlus, applMinus = P_applicable(Conf, I)
        if applPlus || applMinus
            push!(plaqs, Tuple(I))
        end
    end
    return plaqs
end

getApplicablePlaquettes(Conf::StencilSpinConfig) = getApplicablePlaquettes!( Tuple{Int,Int}[], Conf)

getApplicablePlaquettes(Conf::StencilSpinConfig,::Nothing) = getApplicablePlaquettes(Conf)