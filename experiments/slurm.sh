#!/bin/bash
#SBATCH --time=03-00:00:00
#SBATCH --ntasks=13
#SBATCH --partition=cpulong
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G

module load Gurobi
module load Julia

SRC="./experiments/run_stats_experiment.jl"

echo "start"

srun --ntasks=1 -c 4 --exclusive julia --project "$SRC" "feasible" "kuka" "warm" 1000 &

srun --ntasks=1 -c 4 --exclusive julia --project "$SRC" "feasible" "rand_6rad" "warm" 1000 &
srun --ntasks=1 -c 4 --exclusive julia --project "$SRC" "feasible" "rand_4rad" "warm" 1000 &
srun --ntasks=1 -c 4 --exclusive julia --project "$SRC" "feasible" "rand_orth" "warm" 1000 &

srun --ntasks=1 -c 4 --exclusive julia --project "$SRC" "feasible" "rand_6rad" "cold" 1000 &
srun --ntasks=1 -c 4 --exclusive julia --project "$SRC" "feasible" "rand_4rad" "cold" 1000 &
srun --ntasks=1 -c 4 --exclusive julia --project "$SRC" "feasible" "rand_orth" "cold" 1000 &

srun --ntasks=1 -c 4 --exclusive julia --project "$SRC" "feasible" "icub_v2_7" "warm" 1000 &
srun --ntasks=1 -c 4 --exclusive julia --project "$SRC" "feasible" "icub_v2_8" "warm" 1000 &
srun --ntasks=1 -c 4 --exclusive julia --project "$SRC" "feasible" "icub_v2_9" "warm" 1000 &
srun --ntasks=1 -c 4 --exclusive julia --project "$SRC" "feasible" "icub_v2_10" "warm" 1000 &

srun --ntasks=1 -c 4 --exclusive julia --project "$SRC" "uniform" "kuka" "warm" 1000 &
srun --ntasks=1 -c 4 --exclusive julia --project "$SRC" "uniform" "icub_v2_7" "warm" 1000 &

echo "waiting"
wait
echo "done"
