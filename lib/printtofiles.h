#include <stdio.h>

void printtofiles(int a) {
	sprintf(infilen,"data/inside-%d.dat",k);
	sprintf(outfilen,"data/outside-%d.dat",k);
	sprintf(tfilen,"data/table-%d.dat",k);
	FILE *infile=fopen(infilen,"w"),*outfile=fopen(outfilen,"w"),*tfile=fopen(tfilen,"w");
	int i,j;
	for (i=1;i<inn-1;i++) {
        for (j=1;j<inm-1;j++)
            fprintf(infile,"%lf ",inu[i][j]);
        fprintf(infile,"\n");
    }
    for (i=1;i<outn-1;i++) {
        for (j=0;j<outm-1;j++) {
            fprintf(outfile,"%lf ",outu[i][j]);
            fprintf(tfile,"%lf %lf\n",AAs[j],Psis[j]);
        }
        fprintf(outfile,"\n");
    }
	fclose(infile);
	fclose(outfile);
	fclose(tfile);
	if (a!=0)
		exit(-1);
}
