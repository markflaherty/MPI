#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "this_MPI.h"
int main(int argc, char** argv)		{
    int     count;                /* local prime count */
    double  elapsed_time;         /* parallel execution time */
    int     first;                /* index of first multiple */
    int     global_count;         /* global prime count */
    int     high_value;           /* highest value on this proc */
    int     i;
    int     id;                   /* process id number */
    int     index;                /* index of current prime */
    int     low_value;            /* lowest value on this proc */
    int     n;                    /* sieving from 2, ..., n */
    int     p;                    /* number of processes */
    int     proc0_size;           /* size of proc 0's subarray */
    int     prime;                /* current prime */
    int     size;                 /* elements in marked string */
    int     first_value_index;
    int     prime_step;
    int     prime_doubled;
    int     sqrt_n;
    int     prime_multiple;
    int     num_per_block;
    int     block_low_value;
    int     block_high_value;
    int     first_index_in_block; 
    char*   marked;               /* portion of 2, ..., n */
    char*   primes;

    MPI_Init(&argc, &argv);

    /* start the timer */
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (argc != 2)    {
        if (id == 0) /* parent process */
            printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    } /* if (argc != 2) */

    n = atoi(argv[1]);

    /* 
     * Figure out this process's share of the array, as well as the 
     * integers represented by the first and last array elements 
     */
    low_value  = BLOCK_FIRST + BLOCK_LOW(id, p, n - 1)  * BLOCK_STEP;
    high_value = BLOCK_FIRST + BLOCK_HIGH(id, p, n - 1) * BLOCK_STEP;
    size       = BLOCK_SIZE(id, p, n - 1);

    /*
     * bail out if all the primes used for sieving are not all 
     * help by process 0
     */
    proc0_size = (n - 1) / p;

    if ((2 + proc0_size) < (int)sqrt((double)n))    {
        if (id == 0) /* parent process */
            printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    } /* if */
    
    // compute primes from 2 to sqrt(n);
    sqrt_n = sqrt(n);
    primes = (char*)calloc(sqrt_n + 1, 1);
    for (prime_multiple = 2; 
         prime_multiple <= sqrt_n; 
         prime_multiple += 2)    {
        primes[prime_multiple] = 1;
    } /* for */

    for (prime = 3; prime <= sqrt_n; prime += 2)    {
        if (primes[prime] == 1)
            continue;

        for (prime_multiple = prime << 1;
             prime_multiple <= sqrt_n; 
             prime_multiple += prime)    {
            primes[prime_multiple] = 1;
        }
    } /* for */

    /* 
     * allocate this process' share of the array 
     */
    marked = (char*)calloc(size * sizeof(char), 1);
    if (marked == NULL)    {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    } /* if */

    num_per_block    = 1024 * 1024;
    block_low_value  = low_value;
    block_high_value = MIN(high_value, 
                           low_value + num_per_block * BLOCK_STEP);
    
    for (first_index_in_block = 0;
         first_index_in_block < size; 
         first_index_in_block += num_per_block)    {
        for (prime = 3; prime <= sqrt_n; prime++)       {
            if (primes[prime] == 1)
                continue;
            if (prime * prime > block_low_value)   {
                first = prime * prime;
            }
           else   {
                if (!(block_low_value % prime))    {
                    first = block_low_value;
                }
                else    {
                    first = prime - (block_low_value % prime) + 
                            block_low_value;
                }
           }
        
           /*
            * optimization - consider only odd multiples 
            *                of the prime number
            */
           if ((first + prime) & 1) // is odd 
              first += prime;

           first_value_index = (first - BLOCK_FIRST) / BLOCK_STEP - 
                               BLOCK_LOW(id, p, n - 1);
           prime_doubled     = prime << 1;
           prime_step        = prime_doubled / BLOCK_STEP;
           for (i = first; i <= high_value; i += prime_doubled)   {
               marked[first_value_index] = 1;
               first_value_index += prime_step;
           } /* for */
        }
        
        block_low_value += num_per_block * BLOCK_STEP;
        block_high_value = MIN(high_value, 
                          block_high_value + num_per_block * BLOCK_STEP); 
    } /* for first_index_in_block */


    /* 
     * count the number of prime numbers found on this process 
     */
    count = 0;
    for (i = 0; i < size; i++)
        if (!marked[i])
            count++;

    MPI_Reduce(&count, &global_count, 1, MPI_INT, 
               MPI_SUM, 0, MPI_COMM_WORLD);

    /*
     * stop the timer 
     */
    elapsed_time += MPI_Wtime();

    /* print the results */
    if (id == 0)   {
        global_count += 1; /* add first prime, 2 */
        printf("%d primes are less than or equal to %d\n", 
               global_count, n);
        printf("Total elapsed time: %10.6fs\n", 
               elapsed_time);
    } /* if */

    MPI_Finalize();

    return 0;
}

















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

