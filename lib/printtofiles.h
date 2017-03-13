#include <stdio.h>

void printtofiles(int a) {
	sprintf(allfilen,"data/all-%d.dat",k);
	sprintf(tfilen,"data/table-%d.dat",k);
	FILE *allfile=fopen(allfilen,"w"),*tfile=fopen(tfilen,"w");
	int i,j;
    for (i=0;i<alln-1;i++) {
        for (j=1;j<allm-1;j++)
            fprintf(allfile,"%lf ",allu[i][j]);
        fprintf(allfile,"\n");
        fprintf(tfile,"%lf %lf\n",AAs[i],Psis[i]);
    }
	fclose(allfile);
	fclose(tfile);
	if (a!=0)
		exit(-1);
}
