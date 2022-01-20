#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --array=1-11
#SBATCH --cpus-per-task=9
#SBATCH --mem=32G
#SBATCH --account=def-egreenbl
#SBATCH --job-name=align_genome
#SBATCH --output=output/%x-%j.out
#SBATCH --mail-user=keeganfl@student.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# list should contain all of the
# Commands you wish to run.

echo "Starting task $SLURM_ARRAY_TASK_ID"
commands=$(sed -n "${SLURM_ARRAY_TASK_ID}p" step_18_genome_align_list)

# Then execute all of the commands in parrallel.
eval $commands
