
struct SimpleJastrowFunction{T<:Number} <: AbstractGuidingFunction
    m_i::Vector{T}
    v_ij::Matrix{T}
end
function SimpleJastrowFunction(N::Int,Type = Float64)
    m_i = zeros(Type,N)
    v_ij = zeros(Type,N,N)
    return SimpleJastrowFunction(m_i,v_ij)
end
SimpleJastrowFunction(S::StencilSpinConfig,Type = Float64) = SimpleJastrowFunction(length(S),Type)

get_m_i(ψG::SimpleJastrowFunction) = ψG.m_i
get_v_ij(ψG::SimpleJastrowFunction) = ψG.v_ij

function get_params(ψG::SimpleJastrowFunction)
    return RecursiveArrayTools.ArrayPartition(ψG.m_i,get_v_ij(ψG))
end

function (ψG::SimpleJastrowFunction)(x::StencilSpinConfig)
    _evaluate_jastrow_simple(ψG,x)
end

guidingfunc_name(F::SimpleJastrowFunction) = "SimpleJastrowFunction"
(ψG::SimpleJastrowFunction)(W::SpiderWebWalker) = _evaluate_jastrow_simple(ψG,get_config(W))

function _evaluate_jastrow_simple(ψG,x::AbstractMatrix)
    m = get_m_i(ψG)
    v = get_v_ij(ψG)

    exp_m = zero(eltype(m))
    for i in eachindex(IndexLinear(),x)
        exp_m += m[i] * x[i]
    end
    x_lin = reshape(x, length(x))
    exp_v = 0.5*dot(x_lin,v,x_lin)

    return exp(exp_m + exp_v)
end

struct SimpleJastrow_GWF_Buffer{T<:Number}
    x_i::Vector{T}
    h_i::Vector{T}
    prefac_moves::Matrix{T}
end

function GWFBuffer_set_to!(A::SimpleJastrow_GWF_Buffer,B::SimpleJastrow_GWF_Buffer)
    A.x_i .= B.x_i
    A.h_i .= B.h_i
    # A.prefac_moves .= B.prefac_moves
    # A.safe_parent_indices .= B.safe_parent_indices
    return A
end


function allocate_GWF_buffer(ψG::SimpleJastrowFunction{T},S::AbstractMatrix) where T
    x_i = zeros(T,length(S))
    h_i = zeros(T,length(S))
    
    prefac_moves = _precompute_prefac_moves(ψG,S)
    return SimpleJastrow_GWF_Buffer(x_i,h_i,prefac_moves)
end

function _precompute_prefac_moves(ψG::SimpleJastrowFunction,Conf::StencilSpinConfig)
    AllPlaqs = collect(plaquetteIterator(Conf))
    vij = get_v_ij(ψG)
    AllWeights = fill(NaN,size(Conf))
    for (i,j) in AllPlaqs
        AllWeights[i,j] = _precompute_jastrow_weight(vij,(i,j),Conf)
    end
    return AllWeights
end

function premove_update_GWF_buffer!(Buffer::SimpleJastrow_GWF_Buffer,ψG::SimpleJastrowFunction,Walker::SpiderWebWalker) 
    return Buffer
end

Base.@propagate_inbounds function compute_GWF_buffer!(Buffer::SimpleJastrow_GWF_Buffer,ψG::SimpleJastrowFunction,Walker::SpiderWebWalker)

    x = Buffer.x_i

    x .= reshape(get_config(Walker),length(x))

    v = get_v_ij(ψG)
    h = Buffer.h_i
    # fill!(h,0.0)
    @boundscheck checkbounds(x,axes(v,1))
    @boundscheck checkbounds(h,axes(v,1))

    LoopVectorization.@turbo for i in axes(v,1)
        hi = 0.
        for j in axes(v,2)
            hi += x[j]*v[j,i]
        end
        h[i] = hi
    end
    # mul!(h,v,x) 
    return Buffer
end


function post_move_update_GWF_buffer!(Buffer::SimpleJastrow_GWF_Buffer,ψG::SimpleJastrowFunction,Walker::SpiderWebWalker,move::Tuple)
    Config = get_config(Walker)

    i,j,opSign = move

    affected_sites = safe_parent_indices(Config, (i, j))
    safe_iterator = safe_iterate_sites(Config,(i,j))

    v = get_v_ij(ψG)
    h = Buffer.h_i
    
    # LoopVectorization.@turbo for (index,i) in enumerate(affected_sites)

    # _jastrow_update_kernel!(safe_iterator,Config,h,v,affected_sites,opSign)
    LI = LinearIndices(Config)
    for index in safe_iterator
        
        i = LI[CartesianIndex(affected_sites[index])]
        s = P1_STENCIL[index]*opSign
        
        LoopVectorization.@turbo for j in eachindex(h)
            h[j] += v[j,i]*s
        end
    end

    return Buffer
end

function guidingfuncRatio_log(ψG::SimpleJastrowFunction,Walker::SpiderWebWalker,move::Tuple,Buffer::SimpleJastrow_GWF_Buffer)
    m = get_m_i(ψG)
    Config = get_config(Walker)

    (;h_i,prefac_moves) = Buffer
    i,j,opSign = move
    prefac = prefac_moves[i,j]

    # sites = Buffer.safe_parent_indices[i, j]
    sites = safe_parent_indices_linear(Config, (i, j))

    exp_h = zero(eltype(get_params(ψG)))
    exp_m = zero(eltype(get_params(ψG)))

    @inbounds @simd for idx in safe_iterate_sites(Config,(i,j))
        I = sites[idx]
        s = P1_STENCIL[idx]
        exp_h += h_i[I]*s
        exp_m += m[I]*s
    end

    return opSign*(exp_h + exp_m) + prefac
end

function getOx_k(ψG::SimpleJastrowFunction,W::SpiderWebWalker,k)
    par = get_params(ψG)
    x = get_config(W)
    paramtype,k = _getParamsTypeAndIndex(par,k)
    if paramtype == 1
        return x[k]
    end

    vij = get_v_ij(ψG)
    
    if paramtype == 2
        i,j = Tuple(CartesianIndices(vij)[k])
        return 0.5*x[i] * x[j]
    end
    error("Invalid parameter type")
end

function generate_equivalent(type,k,T::TranslationalSymmetry,ψ::SimpleJastrowFunction,S::AbstractMatrix) 
    
    if type == 1
        site = index_to_site(k,S)
        sites = generate_equivalent_sites(site,T,S)
        return site_to_index.(sites,Ref(S))
    elseif type == 2
        siteI,siteJ = index_to_site_pair(k,S,ψ.v_ij)
        sites = generate_equivalent_site_pairs(siteI,siteJ,T,S)
        i_sites = getindex.(sites,1)
        j_sites = getindex.(sites,2)
        sites2 = generate_equivalent_site_pairs(siteJ,siteI,T,S)
        append!(i_sites,getindex.(sites2,1))
        append!(j_sites,getindex.(sites2,2))
        return site_pair_to_index.(i_sites,j_sites,Ref(S),Ref(ψ.v_ij))
    else
        error("Invalid type")
    end

end

function generate_equivalent(type,k,T::ExchangeSymmetry,ψ::SimpleJastrowFunction,S::AbstractMatrix) 
    if type == 2
        i,j = index_to_site_pair(k,S,ψ.v_ij)
        return [k,site_pair_to_index(j,i,S,ψ.v_ij)]
    else
        return [k]
    end
end