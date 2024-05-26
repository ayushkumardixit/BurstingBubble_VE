# Viscoelastic bursting bubbles
 Here, we simulate viscoelastic bursting bubbles

Locally: 
qcc -O2 -Wall -disable-dimensions -fopenmp burstingBubbleVE_v4.c -o burstingBubbleVE_v4 -lm
export OMP_NUM_THREADS=4
./burstingBubbleVE_v4

Remotely: (like on Snellius)
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions burstingBubbleVE_v4_Snellius.c -o burstingBubbleVE_v4_Snellius -lm

sbatch script: 

#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=1:00:00
#SBATCH --partition=genoa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=v.sanjay@utwente.nl

srun --mpi=pmi2 -n 32 --gres=cpu:32 --mem-per-cpu=1750mb burstingBubbleVE_v4_Snellius