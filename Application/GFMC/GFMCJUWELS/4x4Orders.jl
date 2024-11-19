#!/bin/bash
#=
#!/bin/bash

#SBATCH --account=pmfrg

#SBATCH --job-name=2ECompGFMC                 # replace name
#SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de  # replace email address
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=90GB         # memory , more means less gc time
#SBATCH --time=0-24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=END
#SBATCH --output=/p/project/pmfrg/niggemann1/JobsOutput/Spiderweb/GFMC/EnergySweep/4x4Spin1_%a.out    # File to which standard Out- will be written

jutil env activate -p pmfrg
cd $PROJECT/niggemann1
module --force purge
module load Stages/2024  
module load GCCcore/.12.3.0

module load Julia/1.9.3
export JULIA_DEPOT_PATH=/p/scratch/pmfrg/niggemann1/.julia/

julia -O3 -t $SLURM_CPUS_PER_TASK /p/project/pmfrg/niggemann1/Jobs/SpiderWebModel.jl/Application/GFMC/GFMCJUWELS/4x4Orders.jl ${SLURM_ARRAY_TASK_ID}
exit
=#

import Pkg
Pkg.activate(@__DIR__)
cd(@__DIR__)
import SpiderWebModel as SW
using SpiderWebModel.HDF5
function upscale(Conf,L)
    S = similar(Conf,L,L)
    per = SW.PeriodicMatrix(Conf,L,L)
    for I in CartesianIndices(S)
        S[I] = per[I]
    end
    S
end

##
function generatePeriodic(L)
    S = SW.stencilConfig(zeros(L,L),1;
    boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )
    # SW.rand!(S)
    UCSIt = collect(
        Iterators.product((-1:1 for i in 1:L,j in 1:L)...)
    )
    function isGS(UC)
        S[:] .= 2 .*UC
        SW.fulFillsConstraint(S)
    end
    UCS = Matrix{Int8}[]

    for UC in UCSIt
        isGS(UC) || continue
        push!(UCS,copy(parent(parent(S))))
    end
    return UCS
end
a = generatePeriodic(4)

function filterConfs(UCs)
    Lx,Ly = size(UCs[1])
    S = SW.stencilConfig(zeros(Lx,Ly),1;
    boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )
    S_minus = copy(S)
    S_flip = copy(S)
    S_xflip = copy(S)
    S_yflip = copy(S)
    SUC = copy(parent(parent(S)))

    aSet = Set(empty(a))
    moves = empty!([(1,1,1)])
    function isUnique(UC)
        # S_periodic = SW.PeriodicMatrix(UC,Lx,Ly)
        S .= UC
        S_minus .= -UC
        S_flip .= UC'
        S_xflip .= @view UC[end:-1:1,:]
        S_yflip .= @view UC[:,end:-1:1]

        for S′ in (S,S_minus,S_flip,S_xflip,S_yflip)
            for T_x in 0:Lx-1, T_y in 0:Ly-1

                # SUC = parent(parent(S′))
                SUC_per = SW.PeriodicMatrix(parent(parent(S′)))
                # @info "" size(SUC_per[T_x:T_x+Lx-1,T_y:T_y+Ly-1])
                SUC .= @view SUC_per[T_x:T_x+Lx-1,T_y:T_y+Ly-1]

                if SUC in aSet
                    return false
                end
                if SUC' in aSet
                    return false
                end
                if -SUC in aSet
                    return false
                end
                if -SUC' in aSet
                    return false
                end
                SW.getMoves!(moves,S′)
                for (i,j,sgn) in moves
                    SW.applyPlaquette!(S′,i,j,sgn)
                    if SUC in aSet
                        return false
                    end
                    SW.applyPlaquette!(S′,i,j,-sgn)
                end
            end
        end
        return true
    end
    for UC in UCs
        isUnique(UC) || continue
        push!(aSet,UC)
    end
    return aSet
    
end
a_st = sort!(collect(filterConfs(a)),by=x->sum(abs,x))
##
function makeConf(UC,L)
    S = SW.stencilConfig(zeros(L,L),1;
    boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )
    S .= SW.getPeriodicState(UC,L,L)
    return S
end

a_configs = let
    L = 4
    confs = [makeConf(UC,L) for UC in a_st]

    # filter!(SW.fulFillsConstraint, confs)
    filter!(x->length(SW.getApplicablePlaquettes(x)) > 0,confs)
    @assert all(SW.fulFillsConstraint,confs)
    confs
end
##

function getMaxFlipConfs(Configs;Nwalkers = 28*2,NSteps = 1000)
    filterConfs = Set(empty(Configs))
    mu = -0.8
    CT = SW.ContinuousTimeMethod(0.5,Hxx = SW.Hxx_RK(mu))
    ψG = SW.PlaquetteNumberGuidingFunction(0.5)
    Sbuff = copy(parent(parent(Configs[1])))
    ConfBuff = copy(Configs[1])
    moves = empty!([(1,1,1)])

    for S in Configs
        res = SW.startManyWalkerGFMC(S,CT,Nwalkers,NSteps,ψG)
        isNew = true
        maxMoves = 0
        maxConf = copy(ConfBuff)
        for c in eachslice(res.SaveConfigs,dims=(3,4))
            Sbuff .= c
            if Sbuff in filterConfs || -Sbuff in filterConfs
                isNew = false
                break
            end
            ConfBuff .= Sbuff
            SW.getMoves!(moves,ConfBuff)
            if length(moves) > maxMoves
                maxConf .= ConfBuff
                maxMoves = length(moves)
            end
        end
        if isNew
            push!(filterConfs,maxConf)
        end
    end
    filterConfs_arr = collect(filterConfs)

    sort!(filterConfs_arr,by = x->length(SW.getMoves!(moves,x)),rev = true)
end

reducedConfigs = getMaxFlipConfs(a_configs;Nwalkers = 48,NSteps = 1) #first reduction
##
reducedConfigs = getMaxFlipConfs(reducedConfigs;Nwalkers = 48,NSteps = 10)
##
reducedConfigs = getMaxFlipConfs(reducedConfigs;Nwalkers = 48,NSteps = 100)
##
reducedConfigs = makeConf.(filterConfs(parent.(parent.(reducedConfigs))),4)
reducedConfigs = getMaxFlipConfs(reducedConfigs;Nwalkers = 48*1,NSteps = 1000)
##

function findEnergies(Configs,CT,ψG;Nwalkers = 28*2,NSteps = 2000,equilibration_steps = NSteps ÷8,outfileDIR = nothing)
    en = zeros(length(Configs))
    Δen = zeros(length(Configs))
    getOutfile(i) = isnothing(outfileDIR) ? nothing : joinpath(outfileDIR,"$(i).h5")

    Threads.@threads for i in eachindex(Configs,en)
        S = Configs[i]
        if length(SW.getApplicablePlaquettes(S)) == 0
            en[i] = 0
            Δen[i] = 0
            continue
        end
        results = [SW.startManyWalkerGFMC(S,CT,Nwalkers,NSteps,ψG;equilibration_steps,pre_equilibration_steps=NSteps,outfile = getOutfile(i_st)) for i_st in 1:6]

        energies = SW.getEnergies.(results,1,min(NSteps÷3, 300))
        energiesMean = mean.(energies)
        energiesStd = std.(energies)
        e0_index = findFirstMiIndex(energiesMean)
        en[i] = energiesMean[e0_index]
        Δen[i] = energiesStd[e0_index]
    end
    return (;en,Δen)
end
##
L = 16
Nwalkers = 48*20
NSteps = 2000
outfileDIR = "/p/scratch/pmfrg/niggemann1/Spiderweb/DataS1_CT_RK_equil/SectorComp/L=$(L)/runs/"
mkpath(outfileDIR)


mus_sectors = -0.1:0.05:0.99

AllresEn = zeros(length(mus_sectors),length(reducedConfigs))
AllresΔEn = zeros(length(mus_sectors),length(reducedConfigs))

Threads.@threads for (i) in eachindex(mus_sectors)
    mu = mus_sectors[i]
    CT = SW.ContinuousTimeMethod(0.1,1,-0.266length(a_configs[1]),SW.Hxx_RK(mu))
    ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-mu))
    res = findEnergies(upscale.(reducedConfigs,L),CT,ψG;Nwalkers,NSteps,outfileDIR)
    AllresEn[i,:] = res.en
    AllresΔEn[i,:] = res.Δen
end

outfileTotal = "/p/scratch/pmfrg/niggemann1/Spiderweb/DataS1_CT_RK_equil/SectorComp/L=$(L)/result.h5"

h5write(outfileTotal,"energies",AllresEn)
h5write(outfileTotal,"Δenergies",AllresΔEn)
h5write(outfileTotal,"mus",mus_sectors)
AllConfs = stack(parent.(parent.(reducedConfigs)),dims=3)
h5write(outfileTotal,"UCs",AllConfs)