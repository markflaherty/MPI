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
      printf("%d\n",prime);
    }
    if(comm_size > 1){
    	MPI_Bcast(&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
    }
  } while (prime*prime<=n);
  
  count = 0;
  
  for (i = 0; i < size; i++){
  	if (hit[i] == 0){
  		count++;
  		if(rank != 0){
  			printf("%d\n",i);
  		}
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






/*
int main (int argc, char ** argv) {

  int i, p, n, index, rank, comm_size, size, prime, count, global, first;
  long int high, low;
  n = 100;
  char* hit;
  double time;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);
  time = -MPI_Wtime();
  low  = 2+(long int)(rank)*(long int)(n - 1)/(long int)comm_size;
  high = 1 + (long int)(rank + 1) * (long int)(n - 1)/(long int)comm_size;
  size = high-low+1;
  hit = (char*) calloc(size, sizeof(char));
  if (hit == NULL) {
   MPI_Finalize();
   exit(1);
  }
  if ((2 + (n - 1 / comm_size)) < (int) sqrt((double) n)) {
    if (rank == 0) printf("Too many processes.\n");
    MPI_Finalize();
    exit(1);
  }
  do {
    if (prime*prime > low) {
      first = prime*prime - low;
    } 
    else {
      if ((low%prime) == 0) 
      	first = 0;
      else 
      	first = prime - (low % prime);
    }
    for (i = first; i < size; i += prime) 
    	hit[i] = 1;
    if (rank == 0) {
      while (hit[++index]);
      prime = index + 2;
    }
    if (comm_size > 1) 
    	MPI_Bcast(&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
  } while (prime * prime <= n);
  count = 0;
  for(p = 0; p < size; p++){
  	if(hit[p] == 0){
  		count++;
  	}
  }
  if (comm_size > 1) {
    MPI_Reduce(&count, &global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  } 
  else {
    global = count;
  }
  time += MPI_Wtime();
  if (rank == 0) {
    printf("In %f seconds we found %d primes less than or equal to %d.\n",
		time, count, n);
  }
  MPI_Finalize();
  return 0;
}

*/

  /*
  if (rank == 0) 
  	index = 0;
  prime = 2;
  while(prime*prime <= n){
  	if (prime*prime > low) {
      	first = prime*prime - low;
    } 
    else{
    	if(low%prime == 0){
    		first = 0;
    	}
    	else{
    		first = prime-(low%prime);
    	}
    }
    for(i = first; i < size; i+= prime){
    	hit[i] = 1;
    }
    if(rank == 0){
    	while(hit[++index]);
    	prime = index+2;
    }
    if(comm_size > 1){
    	MPI_Bcast(&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
    }
  }
  */









  /*
  // Check for the command line argument.
  if (argc != 2) {
    if (comm_rank == 0) printf("Please supply a range.\n");
    MPI_Finalize();
    exit(1);
  }
  */
  //n = atoi(argv[1]);
  /*
  // Bail out if all the primes used for sieving are not all held by
  // process zero.
  if ((2 + (n - 1 / comm_size)) < (int) sqrt((double) n)) {
    if (comm_rank == 0) printf("Too many processes.\n");
    MPI_Finalize();
    exit(1);
  }
  */
  // Figure out this process's share of the array, as well as the integers
  // represented by the first and last array elements.

/*
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int* prime(int start, int end, int rank, int* arr, int* a){
	int i,j,k,m, p;
	int lim = (int)((end-start)+1)/2;
	int count = 0;
	int distance = end-start;
	printf("%d\t %d\t %s\t %d\t %d\n", rank, start, "yes",count,end);
	if(arr == NULL){
		for(j = 0; j < distance; j++){
			int curr = arr[i];
			if(curr == 0){
				continue;
			}
			for(k = i; k < end; k++){
				if(a[k] != 0){
					if(arr[k] != curr && arr[k]%curr == 0){
						arr[k] = 0;
					}
				}
			}
		}
	}
	else{
		for(p = 0; p < sizeof(a)/sizeof(int); p++){
			int curr = a[p];
			for(i = 0; i < distance; i++){
				if(arr[i]%curr == 0){
					arr[i] = 0;
				}
			}
		}

	}

	/*
	for(i = start; i <= lim; i++){
		int curr = a[i];
		if(curr == 0){
			continue;
		}
		for(k = i; k < end; k++){
			if(a[k] != 0){
				if(a[k] != curr && a[k]%curr == 0){
					a[k] = 0;
				}
			}
		}
	}
	for(m = start; m < end; m++){
		if(a[m] != 0){
			count++;
		}
	}
	*/
/*
	for(m = 0; m < end; m++){
		if(arr[m] == 0){
			count++;
		}
	}
	int* ret = malloc(count*sizeof(int));
	int iter = 0;
	for(p = 0; p < end; p++){
		if(arr[m] != 0){
			ret[iter] = arr[m];
			iter++;
		}
	} 

	return ret;
}

int limit = 100;
int sum = 0;
int main(int argc, char *argv[]){
	int j,p,m,tasks,rank,myStart,myEnd;
	int *nums;
	MPI_Status stat;
	MPI_Init(&argc,&argv); 
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&tasks);
	int jump = (limit/tasks);
	int offset = 0;
	myStart = rank*jump;
	myEnd = myStart+jump;
	int first_elem = rank*(limit-2)/tasks + 2;
    int last_elem = (rank+1)*(limit-2)/tasks - 1 + 2;
    int size = last_elem - first_elem + 1;
    int array[size];
    int next;
    int index_first_multiple;
    int global;
    int local;
    for(j = 0; j < size; j++){
        array[j] = 1; 
    }
	int i = 2;
	while(i*i <= limit){
		if(first_elem%i == 0){
			index_first_multiple = 0;
		}
		else{
			index_first_multiple = i - first_elem%i;
		}
		for(p = index_first_multiple; p < size; p+=i){
			array[i] = 0;
		}
		if(rank == 0){
			array[i-2] = 0;
		}
		if(rank == 0){
			next = i+1;
			while(array[next-2] == 0)
				next+1;
			i = next;
		}
		MPI_Bcast (&i, 1, MPI_INT, 0, MPI_COMM_WORLD);
		//array[i-2] = 0;
	}
	local = 0;
	global = 0;
	for(m = 0; m < size; m++){
		if(array[m]==1){
			local++;
		}
	}
	MPI_Reduce (&local, &global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	printf("%d\n", global);
	/*
	if(rank == 0){
		nums = malloc(limit*sizeof(int));
		for(i = 2; i <= limit; i++){
			nums[i] = i;
		}
		for(p = 1; p < tasks; p++){
			MPI_Send(&nums[offset], jump, MPI_INT,p, 0, MPI_COMM_WORLD);
			offset+=jump;
		}
		int* a = prime(myStart, myEnd, rank, NULL, NULL);
	}
	else{
		nums = malloc(jump*sizeof(int));
		MPI_Recv(nums,jump, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
		int* a = prime(myStart, myEnd, rank, nums, send);
		int w = 0;
		for(m = 0; m < limit; m++){
			if(send[m] != NULL){
				continue;
			}
			else{
				if(w == sizeof(a)/sizeof(int)){
					break;
				}
				else{
					send[m] = a[w];
					w++;
				}
			}
		}

		//MPI_Gather(a, jump, MPI_INT, send, jump, MPI_INT, MPI_COMM_WORLD);
	}
	/*
	if(rank == 0){
		int count = prime(myStart, myEnd);
		MPI_Reduce(&count, &sum, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
		printf(" %d\t %d\t %d\t %d\t %d\t %d\n",count,tasks,jump, rank, myStart, myEnd);
	}
	for(i = 1; i < tasks; i++){
		if(i == rank){
			int count = prime(myStart, myEnd);
			MPI_Reduce(&count, &sum, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
			printf(" %d\t %d\t %d\t %d\t %d\t %d\n",count,tasks,jump, rank, myStart, myEnd);
		}
	}
	*/
	//int count = prime(myStart, myEnd, rank);
	//int count = 0;
	//printf(" %d\t %d\t %d\t %d\t %d\t %d\n",count,tasks,jump, rank, myStart, myEnd);
	//MPI_Reduce(&count, &sum, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
	/*
	printf(" %d\t %d\t %d\t %d\t %d\t %d\n",count,tasks,jump, rank, myStart, myEnd);
	MPI_Reduce(&count, &sum, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
	*/
/*
	MPI_Finalize();
		

}
*/

