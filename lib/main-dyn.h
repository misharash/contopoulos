#include <stdio.h>

int main_loop() {
	int i,j;
	//main loop
	double diff=1e20;
	k=0;
	int MAXITS=100;
	while ((diff>1e-3)&&(k<MAXITS)) {
		diff=0;
		k++;
		//solve inside
		solvePSR(ina,inb,inc,ind,ine,inf,inu,inmask,inn,inm,0.9999);
		printf("Solved inside\n");
		//set outside boundaries as Psiopen
		//bottom
		j=0;
		for (i=0;i<outn;i++)
			outu[i][j]=inu[N][1];
		//right
		i=N+1;
		for (j=0;j<outm;j++)
			outu[i][j]=inu[N][1];
		//solve outside
		solvePSR(outa,outb,outc,outd,oute,outf,outu,outmask,outn,outm,0.5);
		printf("Solved outside\n");
		//fill array for sorting
		for (j=1;j<=N+1;j++) {
			dp[j-1].Psi=inu[N][j]*l1+outu[1][j-1]*l2;
			diff+=sqr(inu[N][j]-outu[1][j-1]);
			dp[j-1].AA=mu1*AA(inu[N][j])+mu2*AA(outu[1][j-1])+mu3*(inu[N][j]-outu[1][j-1]);
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
		//print number of loops and difference
		diff=sqrt(diff/N);
		printf("\n%d %lf\n",k,diff);
		break; //temp!!
	}
	return 0;
}
