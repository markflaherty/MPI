#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int prime(int start, int end, int rank){
	int i,j,k,m;
	int lim = (int)((end-start)+1)/2;
	int count = 0;
	int distance = end-start;
	int a[distance+1];
	for(j = start; j < end; j++){
		a[j] = j;
	}
	for(i = start; i <= lim; i++){
		int curr = a[i];
		printf("%d\t %d\t %d\t %d\n", rank, start, end, curr);
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
	return count;
}
int limit = 100;
int sum = 0;
int main(int argc, char *argv[]){
	int i,tasks,rank,myStart,myEnd;
	MPI_Init(&argc,&argv); 
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&tasks);
	int jump = (limit/tasks);
	myStart = rank*jump;
	myEnd = myStart+jump;
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
	int count = prime(myStart, myEnd, rank);
	//int count = 0;
	//printf(" %d\t %d\t %d\t %d\t %d\t %d\n",count,tasks,jump, rank, myStart, myEnd);
	//MPI_Reduce(&count, &sum, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
	/*
	printf(" %d\t %d\t %d\t %d\t %d\t %d\n",count,tasks,jump, rank, myStart, myEnd);
	MPI_Reduce(&count, &sum, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
	*/
	MPI_Finalize();
		

}
