#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int prime(int start, int end){
	int x = end;
	int i,j,k,m;
	int lim = (int)(x+1)/2;
	int count = 0;
	int a[end-start];
	for(j = start; j <= end; j++){
		a[j] = j;
	}
	for(i = start; i <= lim; i++){
		int curr = a[i];
		if(curr == 0){
			continue;
		}
		for(k = i; k <= end; k++){
			if(a[k] != 0){
				if(a[k] != curr && a[k]%curr == 0){
					a[k] = 0;
				}
				else if(a[k] != curr){
					if(a[k]%2 == 0 || a[k]%3 == 0 || a[k]%5 == 0 || a[k]%7 == 0){
						a[k] = 0;
					}
				}
			}
		}
	}
	for(m = start; m <= end; m++){
		if(a[m] != 0){
			count++;
		}
	}
	return count;
}
int limit = 100;
int main(int argc, char *argv[]){
	int i;
	int tasks;
	int rank;
	int myStart;
	int myEnd;
	int sum = 0;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&tasks);
	int jump = (limit/tasks);
	myStart = rank*jump;
	myEnd = myStart+jump;
	int count = prime(myStart, myEnd);
	printf(" %d\t %d\t %d\t %d\t %d\t %d\n",count,tasks,jump, rank, myStart, myEnd);
	MPI_Reduce(&count, &sum, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
	if(rank == 0)
		printf("%d\n", sum);
	MPI_Finalize();
		

}
