#!/bin/bash
#=
#!/bin/bash

#SBATCH --account=pmfrg

#SBATCH --job-name=RKSpin1                 # replace name
#SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de  # replace email address
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=90GB         # memory , more means less gc time
#SBATCH --time=0-24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=END
#SBATCH --output=/p/project/pmfrg/niggemann1/JobsOutput/Spiderweb/GFMC/RKSpin1_%a.out    # File to which standard Out- will be written

jutil env activate -p pmfrg
cd $PROJECT/niggemann1
module --force purge
module load Stages/2024  
module load GCCcore/.12.3.0

module load Julia/1.9.3
export JULIA_DEPOT_PATH=/p/scratch/pmfrg/niggemann1/.julia/

julia -O3 -t $SLURM_CPUS_PER_TASK /p/project/pmfrg/niggemann1/Jobs/SpiderWebModel.jl/Application/Spin1Flucs/spin1flucs_large.jl ${SLURM_ARRAY_TASK_ID}
exit
=#
cd(@__DIR__)

i_arg = parse(Int, ARGS[1])
using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5

##
# L = 40
# equilibration_steps = 0
# samplingRate=1e-3
##
L = 80
equilibration_steps = 600_000_000
samplingRate=1e-9
##
outfile = "/p/scratch/pmfrg/niggemann1/Spiderweb/DataRK/Spin1RK_L=$(L)_equilibration_steps=$(equilibration_steps)_samplingRate=$(samplingRate)_$(i_arg).h5"
mkpath(dirname(outfile))
@assert !isfile(outfile) "file already exists!"
##
S = SW.stencilConfig(zeros(L,L),1,;boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional())
##
@time confs = SW.getRandConfs(S, 100,Threads.nthreads(); equilibration_steps,samplingRate)
confArr = reshape(confs,L,L,size(confs,3)*size(confs,4))
##
h5write(outfile,"confs",confArr)
h5write(outfile,"equilibration_steps",equilibration_steps)
h5write(outfile,"samplingRate",samplingRate)

