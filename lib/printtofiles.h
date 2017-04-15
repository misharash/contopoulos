#include <stdio.h>

void printtofiles(int a) {
    sprintf(infilen,"data/inside-%d.dat",k);
    sprintf(outfilen,"data/outside-%d.dat",k);
    sprintf(tfilen,"data/table-%d.dat",k);
    FILE *infile=fopen(infilen,"w"),*outfile=fopen(outfilen,"w"),*tfile=fopen(tfilen,"w");
    int i,j;
    for (i=0;i<=N;i++) {
	for (j=0;j<=N;j++) {
	    fprintf(infile,"%lf ",inu[i][j+1]);
	    fprintf(outfile,"%lf ",outu[i+1][j]);
	}
	fprintf(infile,"\n");
	fprintf(outfile,"\n");
	fprintf(tfile,"%lf %lf\n",Psis[i],AAs[i]);
    }
    fclose(infile);
    fclose(outfile);
    fclose(tfile);
    if (a!=0)
	exit(-1);
}
