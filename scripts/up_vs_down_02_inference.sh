#!/usr/bin/env sh
#SBATCH --job-name=up_vs_down
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=b.kopperud@lmu.de
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=logs/up_vs_down.log
#SBATCH --error=logs/up_vs_down.err
#SBATCH --qos=low_prio_res
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=80
#SBATCH --partition=krypton

#module load R/4.2.3 gnu openblas
module load R/4.3.2 gnu openblas

export R_HOME="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R"
export LD_LIBRARY_PATH="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R/lib"
#echo ${SLURM_CPUS_PER_TASK} > output/ntasks.txt

julia --threads ${SLURM_CPUS_PER_TASK} scripts/up_vs_down_inference.jl > logs/up_vs_down.txt
