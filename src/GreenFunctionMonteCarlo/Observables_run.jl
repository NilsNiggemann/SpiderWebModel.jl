struct SqObs_Buffers{T_high<:AbstractFloat,T_low<:AbstractFloat,FFTType<:SqFFT} <: AbstractGFMCObservables
    TotalWeights::CircularArrays.CircularVector{T_high}
    energies::CircularArrays.CircularVector{T_high}
    FFTBuffers::Vector{FFTType}
    SqBuffers::CircularArrays.CircularArray{T_low,4}
    Gnps::CircularArrays.CircularMatrix{T_high}
    reconfigurationTable::CircularArrays.CircularMatrix{Int}
    PopulationMatrix::CircularArrays.CircularMatrix{Int}
end

struct GFMCObservables_StructureFac{T<:AbstractFloat,BuffType<:SqObs_Buffers,T2} <: AbstractGFMCObservables
    Energy::Vector{T}
    Sq_numerator::Array{T,3}
    obs_denominator::Vector{T}
    en_denominator::Vector{T}
    Buffers::BuffType
    outfile::T2
end
get_reconfigurationTable(O::GFMCObservables_StructureFac) = O.Buffers.reconfigurationTable
get_energies(O::GFMCObservables_StructureFac) = O.Buffers.energies
get_TotalWeights(O::GFMCObservables_StructureFac) = O.Buffers.TotalWeights
get_outfile(O::GFMCObservables_StructureFac) = O.outfile

function create_Sq_Buffers(InitConfig,NWalkers,m_proj)
    p_proj = 2m_proj #projection of forward walking needs to be twice the projection of the wavefunction
    Lx,Ly = size(InitConfig)
    
    energies = CircularArrays.CircularArray(zeros(p_proj))
    TotalWeights = CircularArrays.CircularArray(zeros(p_proj))
    reconfigurationTable = CircularArrays.CircularArray(zeros(Int,NWalkers,p_proj))
    PopulationMatrix = CircularArrays.CircularArray(zeros(Int,NWalkers,m_proj))
    Gnps = CircularArrays.CircularArray(zeros(Float64,p_proj,p_proj))
    SqBuffers = CircularArrays.CircularArray(zeros(Float32,Lx,Ly,NWalkers,m_proj))

    FFTBuffers = fetch.([Threads.@spawn SqFFT((Lx,Ly)) for i in 1:NWalkers])
    
    return SqObs_Buffers(TotalWeights,energies,FFTBuffers,SqBuffers,Gnps,reconfigurationTable,PopulationMatrix)
end

function setup_Sq_Observables(InitConfig,NWalkers,NSteps,m_proj,outfile::Nothing)
    Lx,Ly = size(InitConfig)

    Buffers = create_Sq_Buffers(InitConfig,NWalkers,m_proj)

    Energy = zeros(Float64,m_proj)
    en_denominator = zeros(Float64,m_proj)

    Sq_numerator = zeros(Lx,Ly,m_proj)
    obs_denominator = zeros(Float64,m_proj)

    return GFMCObservables_StructureFac(Energy,Sq_numerator,obs_denominator,en_denominator,Buffers,outfile)
end
function setup_Sq_Observables(InitConfig,NWalkers,NSteps,m_proj,filename::String)
    p_proj = 2m_proj
    Lx,Ly = size(InitConfig)
    Buffers = create_Sq_Buffers(InitConfig,NWalkers,m_proj)

    h5open(filename,"cw") do file
        file["NWalkers"] = NWalkers
        file["NSteps"] = NSteps
        Energy = createMMapArray(file,"Energy",Float64,(m_proj,))
        Sq_numerator = createMMapArray(file,"Sq_numerator",Float64,(Lx,Ly,m_proj))
        obs_denominator = createMMapArray(file,"obs_denominator",Float64,(m_proj,))
        en_denominator = createMMapArray(file,"en_denominator",Float64,(m_proj,))
        return GFMCObservables_StructureFac(Energy,Sq_numerator,obs_denominator,en_denominator,Buffers,filename)
    end
end

function setup_Sq_problem(InitialState::StencilSpinConfig,method::AbstractGFMCMethod,Nwalkers::Integer,NSteps,m_proj,nThreads,ψG,outfile)
    setup = setup_many_walker_GFMC(InitialState,Nwalkers,nThreads)
    Guiding_function_buffer = allocate_GWF_buffers_threads(ψG,InitialState,Nwalkers)
    
    (;Walkers,weights,reconfiguration_buffer) = setup
    Observables = setup_Sq_Observables(InitialState,Nwalkers,NSteps,m_proj,outfile)

    return SpiderwebGFMCProblem(method,InitialState,ψG,Walkers,weights,Guiding_function_buffer,reconfiguration_buffer,Observables)
end


function compute_Sq_Walkers!(SqBuffers::AbstractArray{T,4},Walkers,n,FFTBuffers::AbstractVector{<:AbstractObservable}) where T
    @assert axes(SqBuffers,3) == eachindex(FFTBuffers) == eachindex(Walkers)
    pMax = size(SqBuffers,4)

    Threads.@threads for α in eachindex(FFTBuffers)
        conf = get_config(Walkers[α])
        Sq_view = @view SqBuffers[:,:,α,n]

        ObsFunc! = FFTBuffers[α]
        Sq = obs(ObsFunc!)
        ObsFunc!(Sq,conf)
        Sq_view .= Sq
    end
    return
end

function updateGnp!(Gnp,TotalWeights,n)
    nMax,pMax = size(Gnp)
    for p in 1:pMax
        if n-p < 1
            Gnp[n,p] = 0
            continue
        end
        
        Gnp[n,p] = prod(@view TotalWeights[n-p:n])
    end
    return
end
# function updateEnergies!(Observables::GFMCObservables_StructureFac,i,Walkers::AbstractVector{<:AbstractWalker},method)
#     energies = get_energies(Observables)
#     TotalWeights = get_TotalWeights(Observables)
#     nMax = length(energies)
#     energies[i] = getLocalEnergyWalkers_before(weights,Walkers,method)
#     TotalWeights[i] = mean(weights)
#     return nothing
# end

function updateEnergies!(Observables::GFMCObservables_StructureFac,i,Walkers::AbstractVector{<:AbstractWalker},weights,method)
    energies = get_energies(Observables)
    TotalWeights = get_TotalWeights(Observables)
    energies[i] = getLocalEnergyWalkers_before(weights,Walkers,method)
    TotalWeights[i] = mean(weights)
    Gnps = Observables.Buffers.Gnps
    Energy = Observables.Energy
    en_denominator = Observables.en_denominator
    updateGnp!(Gnps,TotalWeights,i)
    getEnergy_step!(Energy,en_denominator,Gnps,energies,i)

    return nothing
end

function saveObservables!(Observables::GFMCObservables_StructureFac,n,Walkers::AbstractVector{<:SpiderWebWalker})

    (;Sq_numerator,obs_denominator,en_denominator,Energy) = Observables
    (;reconfigurationTable,energies,Gnps,TotalWeights,SqBuffers,PopulationMatrix,FFTBuffers) = Observables.Buffers

    compute_Sq_Walkers!(SqBuffers,Walkers,n,FFTBuffers)

    Nw = length(Walkers)
    
    m_max = size(Sq_numerator,3)
    getPopulationMatrix!(PopulationMatrix,reconfigurationTable,n,m_max-1)
    
    Nw⁻¹ = 1/Nw

    m_values = 0:m_max-1
    Threads.@threads for m_index in eachindex(m_values)
        m = m_values[m_index]
        Gnp = Gnps[n,1+2m]
        obs_denominator[m_index] += Gnp
        # obs_denominator[m_index] += Gnp*Nw
        @views for α in 1:Nw
            mult = PopulationMatrix[α,m_index]
            mult == 0 && continue
            mult *= Nw⁻¹
            O = SqBuffers[:,:,α,n-m]
            @. Sq_numerator[:,:,m_index] += O*Gnp*mult
        end
    end

end

function getEnergy_step!(Energy,en_denominator,Gnp,localEnergies,n)
    for p in eachindex(Energy)
        en_denominator[p] += Gnp[n,p]
        Energy[p] += Gnp[n,p]*localEnergies[n]
    end
    return Energy
end

function getPopulationMatrix!(PopulationMatrix,reconfigurationTable::AbstractMatrix,n,projectionLength)
    PopulationMatrix .= 0 
    nMax = size(reconfigurationTable,2)
    for α in axes(reconfigurationTable,1)
        α´ = α
        for i_m in 0:projectionLength
            if n-i_m < 1
                break
            end
            α´ = reconfigurationTable[α´,n-i_m]
            PopulationMatrix[α´,i_m+1] += 1
        end
    end
    return PopulationMatrix
end

function measure_Sq_GFMC(InitialState::StencilSpinConfig,method::AbstractGFMCMethod,Nwalkers::Integer,nSteps::Integer,mProj,ψG; equilibration_steps = 0, pre_equilibration_steps = equilibration_steps ÷ 5, scatter_fraction = 0.8,initializer = UnguidedWalkInitializer(pre_equilibration_steps,scatter_fraction),nThreads=2*Threads.nthreads(),outfile = nothing,estimate_w_avg=true,kwargs...)
    prob = setup_Sq_problem(InitialState,method,Nwalkers,nSteps,mProj,nThreads,ψG,outfile)

    initializeGFMC!(prob,nThreads,initializer)
    runGFMC!(prob,equilibration_steps;nThreads,reconfigure=true,save_energies = true,saveObservables = false)

    if estimate_w_avg
        w_avg = get_w_avg_estimate(prob)
        prob = set_w_avg_estimate(prob,w_avg)
    end
    outfile = get_outfile(prob.Observables)
    saveParameters(outfile,equilibration_steps,method,ψG)
    
    fill_all_Buffers!(prob,nThreads)
    runGFMC!(prob,nSteps;nThreads)

    return prob.Observables    
end
function set_w_avg_estimate(prob::T,w_avg_estimate)::T where {T<:SpiderwebGFMCProblem}
    newmethod = set_w_avg_estimate(prob.method,w_avg_estimate)
    (;InitialState,ψG,Walkers,weights,Guiding_function_buffer,reconfiguration_buffer,Observables) = prob
    SpiderwebGFMCProblem(newmethod,InitialState,ψG,Walkers,weights,Guiding_function_buffer,reconfiguration_buffer,Observables)
end

function normalized_Sq(Observables::GFMCObservables_StructureFac)
    numerator = copy(Observables.Sq_numerator)
    denominator = Observables.obs_denominator
    for i in eachindex(denominator)
        @. numerator[:,:,i] ./= denominator[i]
    end
    return numerator
end

function normalized_En(Observables::GFMCObservables_StructureFac)
    numerator = copy(Observables.Energy)
    denominator = Observables.en_denominator
    @. numerator ./= denominator
    return numerator
end

function ensure_numeric_w_avg!(method)
    if !isfinite(get_w_avg_estimate(method))
        set_w_avg_estimate!(method,0.)
        return true
    end
    return false
end

function get_w_avg_estimate(prob)
    Observables = prob.Observables
    method = prob.method

    energies = get_energies(Observables)
    TotalWeights = get_TotalWeights(Observables)

    Energy = normalized_En(Observables)[begin]
    
    fill!(energies,zero(eltype(energies)))
    fill!(TotalWeights,zero(eltype(TotalWeights)))
    
    w_avg_estimate = -Energy
    return w_avg_estimate
end