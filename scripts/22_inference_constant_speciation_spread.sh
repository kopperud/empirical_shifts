#!/usr/bin/env sh
#SBATCH --job-name=constant_speciation_spread
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b.kopperud@lmu.de
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=logs/constant_speciation_spread.log
#SBATCH --error=logs/constant_speciation_spread.err
#SBATCH --qos=low_prio_res
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=46
#SBATCH --partition=lemmium

#module load R/4.2.3 gnu openblas
#module load openblas
module load gnu/7
#module load R/4.3.2
module load R/4.2.3

export R_HOME="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R"
export LD_LIBRARY_PATH="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R/lib"

julia --threads ${SLURM_CPUS_PER_TASK} scripts/22_inference_constant_speciation_spread.jl > output/screen.txt
