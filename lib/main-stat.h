#include <stdio.h>

int main_loop() {
	//file to write evolution
	FILE *evolfile=fopen("data/evol.dat","w");
	int i,j,kk;
	//fill f's
	for (i=1;i<inn-1;i++)
		for (j=1;j<inm-1;j++)
			inf[i][j]=-AA(inu[i][j]);
	for (i=1;i<outn-1;i++)
		for (j=1;j<outm-1;j++)
			outf[i][j]=-AA(outu[i][j]);
	//main loop
	double diff=1e20,diff2,diffl;
	double **oldinu=matalloc(inn,inm),**oldoutu=matalloc(outn,outm);
	k=0;
	while (diff>1e-3) {
		diff=0;
		k++;
		//solve inside many times
		diff2=1e20;
		kk=0;
		while (diff2>1e-3) {
			matcpy(oldinu,inu,inn,inm);
			solvePSR(ina,inb,inc,ind,ine,inf,inu,inmask,inn,inm,0.9999);
			kk++;
			diff2=0;
			//rewrite f's and count difference
			for (i=1;i<inn-1;i++)
				for (j=1;j<inm-1;j++) {
					diffl=inu[i][j]-oldinu[i][j];
					inu[i][j]=oldinu[i][j]+ewti*diffl;
					inf[i][j]=-AA(inu[i][j]);
					diff2+=sqr(diffl);
				}
			diff2=sqrt(diff2)/N;
			printf("Solved inside %d: diff %lf\n",kk,diff2);
		}
		//set outside boundaries as Psiopen
		//bottom
		j=0;
		for (i=0;i<outn;i++)
			outu[i][j]=inu[N][1];
		//right
		i=N+1;
		for (j=0;j<outm;j++)
			outu[i][j]=inu[N][1];
		//solve outside many times
		diff2=1e20;
		kk=0;
		while (diff2>1e-3) {
			matcpy(oldoutu,outu,outn,outm);
			solvePSR(outa,outb,outc,outd,oute,outf,outu,outmask,outn,outm,0.9999);
			kk++;
			diff2=0;
			//rewrite f's and count difference
			for (i=1;i<outn-1;i++)
				for (j=1;j<outm-1;j++) {
					diffl=outu[i][j]-oldoutu[i][j];
					outu[i][j]=oldoutu[i][j]+ewto*diffl;
					outf[i][j]=-AA(outu[i][j]);
					diff2+=sqr(diffl);
				}
			diff2=sqrt(diff2)/N;
			printf("Solved outside %d: diff %lf\n",kk,diff2);
		}
		//fill array for sorting
		for (j=1;j<=N+1;j++) {
			dp[j-1].Psi=(inu[N][j]+outu[1][j-1])/2;
			diff+=sqr(inu[N][j]-outu[1][j-1]);
			dp[j-1].AA=mu1*AA(inu[N][j])+mu2*AA(outu[1][j-1])+mu3*(inu[N][j]-outu[1][j-1]);
		//fprintf(stderr,"%d %lf\n",k,diff);
		}
		//sort
		qsort(dp,N+1,sizeof(datapair),dpcmp);
		//fill AAs, Psis
		for (i=0;i<=N;i++) {
			AAs[i]=dp[i].AA;
			Psis[i]=dp[i].Psi;
			printf("%lf - %lf; ",Psis[i],AAs[i]);
		}
		//to make AA'(Psi) continuous (maybe wrong)
		//AAs[0]=0;
		//integrate to get A0
		A0=0;
		for (i=0;i<N;i++)
			A0+=(AAs[i]+AAs[i+1])*(Psis[i]-Psis[i+1]);
		A0=sqrt(A0);
		//count derivatives
		A0s=AAs[0]/A0;
		A0ss=-(A0+A0s*tPsi)/sqr(tPsi);
		//print to files
		printtofiles(0);
		//break; //temporarily!
		//print number of loops and difference
		diff=sqrt(diff/N);
		printf("\n%d %lf\n",k,diff);
		fprintf(evolfile,"%d %lf\n",k,diff);
	}
	return 0;
}
