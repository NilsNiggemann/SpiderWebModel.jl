struct UnguidedWalkInitializer <: AbstractGFMCInitializer 
    pre_equilibration_steps::Int
    scatter_fraction::Float64
end

function initialize!(Walkers::AbstractVector{<:SpiderWebWalker},I::UnguidedWalkInitializer)
    if I.pre_equilibration_steps > 0 && I.scatter_fraction > 0
        random_init_walkers!(Walkers,I.pre_equilibration_steps,I.scatter_fraction)
    end
    return
end

function random_init_walkers!(Walkers::AbstractVector{<:SpiderWebWalker},equilibration_steps,fraction=1.0)
    Nw = length(Walkers)

    Threads.@threads for α in eachindex(Walkers)[1:round(Int,fraction*Nw)]
        Walker = Walkers[α]
        
        for _ in 1:equilibration_steps
            movepos = Tuple(rand(Walker.Plaquette_positions))
            movesgn = rand(1:2)
            P_applicable(Walker.Config, movepos)[movesgn] || continue
            applyPlaquette!(Walker.Config, movepos[1], movepos[2], (1,-1)[movesgn])
        end
    end
end

"""draw random configurations with given pre-given weights to accelerate equilibration times"""
struct WeightedConfigsInitializers{T1<:AbstractVector{<:AbstractMatrix},T2<:AbstractVector} <: AbstractGFMCInitializer
    configs::T1
    weights::T2
end

function WeightedConfigsInitializers(SaveConfigs::AbstractArray{<:Number,4},TotalWeights::AbstractVector)
    configs = collect(eachslice(SaveConfigs,dims=(3,4)))
    configsVec = reshape(configs,length(configs))
    weights = [w for w in TotalWeights for _ in axes(SaveConfigs,3)]
    return WeightedConfigsInitializers(configsVec,weights)
end

function WeightedConfigsInitializers(resultsArr::AbstractVector,weight::Symbol=:TotalWeights)

    w1 = WeightedConfigsInitializers(resultsArr[begin].SaveConfigs,getproperty(resultsArr[begin],weight))
    for i in eachindex(resultsArr)[2:end]
        w2 = WeightedConfigsInitializers(resultsArr[i].SaveConfigs,getproperty(resultsArr[i],weight))
        append!(w1.configs,w2.configs)
        append!(w1.weights,w2.weights)
    end
    return w1
end

function initialize!(Walkers::AbstractVector{<:SpiderWebWalker},I::WeightedConfigsInitializers)
    for Walker in Walkers
        rand_conf = StatsBase.sample(I.configs,StatsBase.Weights(I.weights))
        get_config(Walker) .= rand_conf
    end
end

struct CombinedInitializer{I1 <: AbstractGFMCInitializer,I2 <: AbstractGFMCInitializer} <: AbstractGFMCInitializer
    I1::I1
    I2::I2
end

function initialize!(Walkers::AbstractVector{<:SpiderWebWalker},I::CombinedInitializer)
    initialize!(Walkers,I.I1)
    initialize!(Walkers,I.I2)
end

function initializeGFMC!(prob::AbstractGFMCProblem,initializer)
    
    (;Walkers,Observables,method,ψG) = prob
    initialize!(Walkers,initializer)

    return prob,Observables
end

struct VariableTimePropagation{T<:AbstractVector,Proptype<:ContinuousTimeMethod} <: AbstractGFMCInitializer 
    tauRange::T
    Propagator::Proptype
end

function initializeGFMC!(prob::AbstractGFMCProblem,initializer::VariableTimePropagation)
    
    StepNum = 1
    nThreads = length(prob.Walkers)
    accumulated_time = 50.
    for τ in initializer.tauRange
        prob = Accessors.@set prob.method.τ = τ
        
        if accumulated_time >= 20
            fill_all_Buffers!(prob,nThreads)
            accumulate_time = 0.
        end
        
        accumulated_time += τ
        runGFMC!(prob,StepNum,nThreads,true,false,false)
    end
    
    (;Observables) = prob
    return prob,Observables
end
