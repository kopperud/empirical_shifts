#!/usr/bin/env sh
#SBATCH --job-name=downshift
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b.kopperud@lmu.de
#SBATCH --mem=160GB
#SBATCH --output=logs/downshift.log
#SBATCH --error=logs/downshift.err
#SBATCH --qos=normal_prio
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --partition=lemmium

#module load R/4.2.3 gnu openblas
#module load openblas
module load gnu/7
#module load R/4.3.2
module load R/4.2.3

export R_HOME="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R"
export LD_LIBRARY_PATH="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R/lib"
echo ${SLURM_CPUS_PER_TASK} > output/ntasks.txt

julia --threads ${SLURM_CPUS_PER_TASK} scripts/32_inference_downshift.jl > output/screen.txt
