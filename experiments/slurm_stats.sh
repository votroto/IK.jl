#!/bin/bash
#SBATCH --time=23:59:00
#SBATCH --ntasks=13
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G

module load Gurobi
module load Julia

srun -l --multi-prog ./experiments/slurm_stats.conf