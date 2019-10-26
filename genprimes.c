
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#define UPPER_BOUND 10000000000
#define MODE 1  // 0 to forsake CL-ARG
#define PRINT 0 // 1 to print primes.
#define PAGE 2500


#define SetBit(A,k)     ( A[(k/32)] |= (1 << (k%32)) )
#define ClearBit(A,k)   ( A[(k/32)] &= ~(1 << (k%32)) )  
#define TestBit(A,k)    ( A[(k/32)] & (1 << (k%32)) ) 
#define MAX( a,b)  ( (a > b) ? a : b  )
#define MIN( a,b)  ( (a < b) ? a : b  )


int getRootPrimes(long limit , int** rootPrimes){
    int i =0, count =0;
    long llsqrt = ceil(sqrt(limit));
    int prime =2;
    int *arr = (int *)calloc(ceil(limit/32), sizeof(int));
    while(1){
      while (prime <= llsqrt && TestBit(arr,prime))
            prime ++;
      for (i = prime * prime; i <= limit; i += prime)
                SetBit(arr,i);
      prime ++;
      if(prime > llsqrt) break;
    }
    for(i=2;i<limit;i++){
        if(TestBit(arr,i))
        count+=1;
     }   
     
    (*rootPrimes) = (int *)calloc(count, sizeof(int));
    count =0;
    for(i=2;i<limit;i++){
        if(!TestBit(arr,i)){
        (*rootPrimes)[count] = i;
        count+=1;
        }
     }
     return count;
}

int small(long limit){
    int* rootPrimes ;
    long base;
    return getRootPrimes(limit , &rootPrimes);
}

int removeComposites(long baseIndex, int limit , int* arr, int* rootPrimes, int rootCount)
{       int i =0, j =0, count = 0;
        long start, dprime;
        for( i = 0; i< rootCount; i++)
        {   
            dprime = 2*rootPrimes[i];
            start  = MAX(   (baseIndex+ (rootPrimes[i] - baseIndex%rootPrimes[i]))   ,  rootPrimes[i]*rootPrimes[i]  ) -  baseIndex;
            
            
            if( start%2 == 0 )
                start+=rootPrimes[i];
            if(rootPrimes[i]==2)
                dprime = 2;
            
            for(j = start ; j < limit; j+=dprime )
                 SetBit(arr,j);
        }
        for( i = 0; i< limit; i++){
            if(!TestBit(arr,i)){
                count+=1;
               printf(" %ld \n", (i+baseIndex));  
             }
        }
        return count;
}


int main(int argc, char *argv[]){
    
    long limit, n_hi, n_lo;
    if (argc > 1 && MODE == 1)
        limit = (long)pow(10,atoi(argv[1]));
    else
        limit = UPPER_BOUND;

    int p, id , root  = 0;
    double wtime;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if (id == 0)
        wtime = MPI_Wtime();
    if (id != p - 1)
       n_hi = (long)(limit / p) * (id + 1);
    else
       n_hi = limit;
    n_lo = (long)(limit / p) * id;
    

    long lsqrt = (long) ceil(sqrt(limit));
    
    if (lsqrt < 10000){
     if(id == root){  
        int ans = small(limit);
        wtime = MPI_Wtime() - wtime;
        printf("         N        Pi          Time\n");
        printf("  %10ld    %10d  %16f\n", limit, ans, wtime);
        }
        MPI_Finalize();
        exit(0);
    }

    long part_size = n_hi - n_lo;
    int* rootPrimes ;
    long base;
    int rootCount = getRootPrimes(lsqrt , &rootPrimes);
    int d = 0;
    long i;
    int *arr = (int *)calloc((int)ceil((n_hi-n_lo)/32), sizeof(int));
    if(n_lo == 0){
        SetBit(arr,1);
        SetBit(arr,4);       
    }
    SetBit(arr,0);

   long j = 0, count = 0;
   long start, dprime;
    
    for( i = 0; i < rootCount; i++)
    {   
        dprime = 2 * rootPrimes[i];
        if(n_lo == 0) start = ((long)rootPrimes[i])*((long)rootPrimes[i]);
        else start  = MAX(   ( (long)ceil( ((double)n_lo) / rootPrimes[i]) * rootPrimes[i])   ,  ((long)rootPrimes[i])*((long)rootPrimes[i])  ) -  n_lo;  

        if( start % 2 == 0 && rootPrimes[i] > 2)
            start += rootPrimes[i];
        if(rootPrimes[i] == 2)
            dprime = 2;

        for(j = start ; j < part_size; j += dprime )
                SetBit(arr,j);
    }

    int* all_range;
    if(id == root){
        all_range = (int*) calloc((limit/32), sizeof(int));
    }

    MPI_Gather(arr, (part_size/32), MPI_INT, all_range, (part_size/32), MPI_INT, root, MPI_COMM_WORLD);

   long total_count = 1;
    if(id == root){
        for(i = 1; i < limit; i += 2){
            if (!TestBit(all_range,i)){
                    total_count += 1;
                    if(PRINT == 1)
                    {   printf("%ld \n", i);
                    }
            }
        }
   }
    if (id == root) {
        wtime = MPI_Wtime() - wtime;
        printf("         N        Pi          Time\n");
        printf("  %10ld %10ld  %16f\n", limit, total_count, wtime);
    }
       
    MPI_Finalize();
    exit(0);
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
