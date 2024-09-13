#!/bin/bash
#=
#!/bin/bash
# SBATCH --dependency=afterok:8745821
#SBATCH --job-name=eval
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=3-10:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/GFMCCTRK/eval%a.out
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=220GB
# SBATCH --qos=cont
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
~/.bashrc
julia -O3 -t $SLURM_CPUS_PER_TASK /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/evalScript.jl $SLURM_ARRAY_TASK_ID
exit
=#
using Pkg
Pkg.activate(@__DIR__)

@time using HDF5
flush(stdout)
@time import SpiderWebModel as SW
flush(stdout)

SECTOR_NAME  = "Condensate"
L=32
τ = 0.10
outfileDIR = ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/L=$(L)/periodic_RK_Full_$(SECTOR_NAME)"
folders = readdir(outfileDIR,join=true)
function getEns(results)
    en = [SW.getEnergies(res.TotalWeights,res.energies,1,1000÷res.nBra) for res in results]
end

let 
    
    for folder in folders
        GC.gc()
        GC.gc()
        files = readdir(folder,join=true)

        @time AllResults = vcat(SW.readResults.(files)...);

        μmatch = match(r"mu_(-?\d+\.?\d*)", folder)
        μ = parse(Float64, μmatch.captures[1])

        println("μ = ",μ)
        flush(stdout)

        outfileTotal = ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/eval/L=$(L)/mu=$(μ)/Spin1GFMC_Eval_periodic_$(SECTOR_NAME)_mu$(μ)_L$(L).h5"

        mkpath(dirname(outfileTotal))
        # @assert !isfile(outfileTotal) "file already exists!"
        if isfile(outfileTotal)
            println("file already exists!")
            continue
        end


        en  = stack(getEns(AllResults))
        h5open(outfileTotal,"cw") do file
            file["energies"] = en
            file["nBra"] = AllResults[1].nBra
            file["L"] = L
            file["mu"] = μ
            file["tau"] = τ
        end
        
        Threads.@threads for projectionSteps in (10,20,50,100,200,250,500)
        # for projectionSteps in (20,40)
            SqsGFMC = stack(SW.getSqsGFMC(AllResults,projectionSteps),dims=3)
            h5open(outfileTotal,"cw") do file
                file["SqsGFMC/$projectionSteps"] = SqsGFMC
            end
            println(L, " ",projectionSteps)
            flush(stdout)
        
        end
        GC.gc()
        GC.gc()
    end

end
