#include <stdlib.h>
#include <string.h>
#include <stdio.h>

double** matalloc(int n, int m) {
	int i;
	double *d=calloc(n*m,sizeof(double)),**a=malloc(n*sizeof(double*));
	for (i=0;i<n;i++)
		a[i]=d+i*m;
	return a;
}

void matfree(double** a) {
	free(a[0]);
	free(a);
}

void matcpy(double** dest, double** src, int n, int m) {
	memcpy((void*)dest[0],(void*)src[0],n*m*sizeof(double));
}

void matout(double** a, int n, int m) {
	int i,j;
	for (i=0;i<n;i++) {
		for (j=0;j<m;j++)
			printf("%lf ",a[i][j]);
		printf("\n");
	}
}

