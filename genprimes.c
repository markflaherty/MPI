#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int prime(int start, int end){
	int x = end;
	int i,j,k,m;
	int lim = (int)(x+1)/2;
	int count = 0;
	int a[x+1];
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
int main(){
	int limit;
	int tasks;
	int rank;
	int myStart;
	int myEnd;
	int sum = 0;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&tasks);
	printf("Enter a bound: ");
	scanf("%d",&limit);
	int jump = (int)(limit/tasks);
	myStart = rank*jump;
	myEnd = myStart+jump+1;
	int count = prime(myStart, myEnd);
	MPI_Reduce(&count, &sum, 1, MPI_SUM,0,MPI_COMM_WORLD);
	printf("%d\n", sum);
	MPI_Finalize();
		

}
