#!/bin/sh
#slurm sbatch requesting one K20 GPU. CUDA_VISIBLE_DEVICES will be set.
#SBATCH --nodes 1
#SBATCH --partition gpu
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32G
#SBATCH --gres=gpu
#SBATCH --time=20:00:00

module load singularity

srun SINGULARITY_COMMAND