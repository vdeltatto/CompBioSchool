#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gpus=1
#SBATCH --partition=gpu_h100
#SBATCH --time=03:00:00

module purge
module load 2024
module load foss/2024a
module load CMake/3.29.3-GCCcore-13.3.0
module load CUDA/12.6.0
module load GROMACS/2024.3-foss-2024a-CUDA-12.6.0-PLUMED-2.9.2

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#command to run energy minimization:
#gmx mdrun -deffnm en_min -s en_min.tpr -ntmpi 1 -ntomp 2 -v

#command to run NVT minimization:
#gmx mdrun -deffnm nvt -s nvt.tpr -ntmpi 1 -ntomp 2 -v

#command to run NPT minimization:
#gmx mdrun -deffnm npt -s npt.tpr  -ntmpi 1 -ntomp 2 -v

#command to run production:
#gmx mdrun -deffnm prod -s prod.tpr  -ntmpi 1 -ntomp 2 -v -nb gpu -pme gpu -pin on
