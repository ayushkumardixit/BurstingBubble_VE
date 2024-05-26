#!/bin/bash

MAXlevel="10"
De="0.0"
Ohs="1e-2"
Ohp="1e-2"
Oha="1e-4"
Bond="1e-3"
tmax="4"

# qcc -fopenmp -Wall -O2 burstingBubbleVE.c -o burstingBubbleVE -lm
# export OMP_NUM_THREADS=1
# ./burstingBubbleVE

qcc -Wall -O2 burstingBubbleVE.c -o burstingBubbleVE -lm
./burstingBubbleVE $MAXlevel $De $Ohs $Ohp $Oha $Bond $tmax