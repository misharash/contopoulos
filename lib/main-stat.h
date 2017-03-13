#include <stdio.h>

int main_loop() {
	int i,j,kk;
	//fill f's
	for (i=1;i<alln-1;i++)
		for (j=1;j<allm-1;j++)
			allf[i][j]=AA(allu[i][j]);
    allf[alln-2][allm-2]=0; //top right angle is special
    //f's for LC
    for (j=1;j<allm-1;j++)
        allf[N][j]=(allf[N+1][j]-allf[N-1][j])/2/h;
	//main loop
	double diff=1e20,diff2;
	double **oldallu=matalloc(alln,allm);
	k=0;
    //what is diff? (not actual in monopole)
	while (diff>1e-3) {
		diff=0;
		k++;
		//solve many times
		diff2=1e20;
		kk=0;
		while (diff2>1e-2) {
			matcpy(oldallu,allu,alln,allm);
			solvePSR(alla,allb,allc,alld,alle,allf,allu,allmask,alln,allm,0.9999);
			kk++;
			diff2=0;
			//rewrite f's and count difference
			for (i=1;i<alln-1;i++)
				for (j=1;j<allm-1;j++) {
					allf[i][j]=AA(allu[i][j]);
					diff2+=sqr(allu[i][j]-oldallu[i][j]);
				}
            allf[alln-2][allm-2]=0; //top right angle is special
            //f's for LC
            for (j=1;j<allm-1;j++)
                allf[N][j]=(allf[N+1][j]-allf[N-1][j])/2/h;
			diff2=sqrt(diff2)/N;
			printf("Solved %d: diff %lf\n",kk,diff2);
            //printtofiles(-1); //temp
		}
		//fill array for sorting
		for (j=1;j<=N+1;j++) {
			dp[j-1].Psi=allu[N][j];
			dp[j-1].AA=(allu[N+1][j]-allu[N-1][j])/2/h;
		}
		//sort
		qsort(dp,N+1,sizeof(datapair),dpcmp);
		//fill AAs, Psis
		for (i=0;i<=N;i++) {
			AAs[i]=dp[i].AA;
			Psis[i]=dp[i].Psi;
			printf("%lf - %lf; ",Psis[i],AAs[i]);
		}
		//print to files
		printtofiles(0);
		break; //temporarily!
		//print number of loops and difference
		diff=sqr(diff/N);
		printf("\n%d %lf\n",k,diff);
	}
	return 0;
}
