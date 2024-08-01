#!/bin/bash

#SBATCH -J TIC
#SBATCH -n 36

#SBATCH -t 02-23:00:00
#SBATCH -o test_mpi-%j.err
#SBATCH -e test_mpi-%j.out
#SBATCH -D .

echo "EXECUTING MPI!"
mpirun valgrind --leak-check=full --show-reachable=yes --log-file=nc.vg.%p ../FBL_INV -fname blah

