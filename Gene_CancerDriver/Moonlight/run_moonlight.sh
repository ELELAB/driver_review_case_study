#!/bin/bash 
NTHREADS=4
export OMP_NUM_THREADS=$NTHREADS
tsp -L jonathan184243 -N $NTHREADS Rscript moonlightstartscrript.R
 
