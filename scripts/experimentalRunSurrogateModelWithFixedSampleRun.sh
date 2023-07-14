#!/bin/bash
#SBATCH --job-name=smallerExperimentalRun
#SBATCH --time=72:00:00
#SBATCH --mem=5G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-5780
# Load the required modules
# module load gcccore/10.2.0
# module load cmake/3.18.4
# module load eigen/3.3.8
# # Move into folder and run

cd ..
./cpp_code/main experimentalRunSurrogateModelWithFixedSample ${SLURM_ARRAY_TASK_ID} 10





