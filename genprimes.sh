#!/usr/bin/env bash

set -e 

numP="2"
if [[ "$OSTYPE" == "linux-gnu" ]]; then
	module load mpi/openmpi-x86_64
fi
mpicc -g -Wall -o genprimes genprimes.c -std=c99
mpirun -n $numP ./genprimes $1 > output/$1

#cat output/$numN