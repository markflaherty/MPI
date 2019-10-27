#!usr/bin/env bash

set -e 

numP = "2"
numN = $1
numN = "100"
if [[ "$OSTYPE" == "linux-gnu" ]]; then
	module load mpi/openmpi-x86_64
mpicc -g -Wall -o genprimes genprimes.c -std=c99
mpirun -n $numP ./genprimes $1 > outputs/$1

cat outputs/$numN