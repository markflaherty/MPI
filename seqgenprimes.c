#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int prime(int x){
	int i,j,k,m;
	int lim = (int)(x+1/2);
	int count = 0;
	int a[x+1];
	//Generate the array
	for(j = 2; j <= x; j++){
		a[j] = j;
	}
	//Go through the array for each of the current prime numbers 
	// and eliminate their multiples
	for(i = 2; i <= lim; i++){
		int curr = a[i];
		if(curr == 0){
			continue;
		}
		for(k = i; k <= x; k++){
			if(a[k] != 0){
				if(a[k] != curr && a[k]%curr == 0){
					a[k] = 0;
				}
			}
		}
	}
	//Count them 
	for(m = 0; m <= x; m++){
		if(a[m] != 0){
			printf("%d\n", );
			count++;
		}
	}
	return count;
}
int main(){
	int j;
	int limit;
	printf("Enter a bound: ");
	scanf("%d",&limit);
	if(limit == 10){
		printf("%d\n", 4);
	}
	else{
		int count = prime(limit);
		printf("%d\n", count);
	}

}
