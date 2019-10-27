#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
int main (int argc, char ** argv) {
  int i;
  int n = 100;
  int index;
  int size;
  int prime;
  int count;
  int g;
  int first;
  long int high;
  long int low;
  int rank;
  int comm_size;
  char * hit;
  double t;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);
  t = -MPI_Wtime();
  low  = 2 + (long int)(rank) * (long int)(n - 1) / (long int)comm_size;
  high = 1 + (long int)(rank + 1) * (long int)(n - 1) / (long int)comm_size;
  size = high - low + 1;
  hit = (char *) calloc(size, sizeof(char));
  if (hit == NULL) {
   MPI_Finalize();
   exit(1);
  }
  if (rank == 0){
  	index = 0;
  } 
  prime = 2;
  do {
    if (prime*prime > low) {
      first = prime*prime-low;
    } 
    else {
      if ((low%prime) == 0){
      	first = 0;
      } 
      else{
      	first = prime - (low % prime);
      } 
    }
    
    for (i = first; i < size; i += prime){
    	hit[i] = 1;
    }
    if (rank == 0) {
      while (hit[++index]);
      prime = index + 2;
    }
    if(comm_size > 1){
    	MPI_Bcast(&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
    }
  } while (prime*prime<=n);
  count = 0;
  for (i = 0; i < size; i++){
  	if (hit[i] == 0){
  		count++;
  		printf("%d\n", i+low);
  	} 

  }
  if (comm_size > 1) {
    MPI_Reduce(&count, &g, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  } 
  else {
    g = count;
  }
  
  t += MPI_Wtime();
  if (rank == 0) {
    printf("In %f seconds we found %d primes less than or equal to %d.\n",
		t, g, n);
  }
  
  MPI_Finalize();
  return 0;
}