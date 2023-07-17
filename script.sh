#!/bin/bash

#SBATCH --partition=multiple_il
#SBATCH --job-name=sgCavity
#SBATCH --nodes=1           # number of nodes to use
#SBATCH --ntasks-per-node=40 # number of tasks to run on each node
#SBATCH --time=5:00:00

#SBATCH --error=job.%J.err 
#SBATCH --output=job.%J.out

#SBATCH --mail-type=ALL
#SBATCH --mail-user=mingliang.zhong@kit.edu

module load compiler/intel/2022.2.1
module load  mpi/openmpi/4.1

# edit executable

mpirun ./cavity2d
