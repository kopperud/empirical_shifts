#!/usr/bin/env sh
#SBATCH --job-name=backbone
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b.kopperud@lmu.de
#SBATCH --mem=200GB
#SBATCH --output=logs/backbone.log
#SBATCH --error=logs/backbone.err
#SBATCH --qos=low_prio_res
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=lemmium

#module load R/4.2.3 gnu openblas
#module load openblas
module load gnu/7
#module load R/4.3.2
module load R/4.2.3

export R_HOME="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R"
export LD_LIBRARY_PATH="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R/lib"
echo ${SLURM_CPUS_PER_TASK} > output/ntasks.txt

julia --threads ${SLURM_CPUS_PER_TASK} scripts/30_inference_backbone.jl > output/screen.txt
