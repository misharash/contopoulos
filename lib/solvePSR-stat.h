#include <stdio.h>

//R_LC=1
void solvePSR(double **a, double **b, double **c, double **d, double **e,
	double **f, double **u, double** mask, int n, int m, const double rjac)
{
	const int MAXITS=10000;
	const double EPS=1e-3;
	double anormf=0.0,omega=1.0,anorm,resid;
	int j,l,k,ipass,jsw,lsw;
	for (j=1;j<n-1;j++)
		for (l=1;l<m-1;l++) {
			resid=a[j][l]*u[j+1][l]+b[j][l]*u[j-1][l]
					+c[j][l]*u[j][l+1]+d[j][l]*u[j][l-1]
					+e[j][l]*u[j][l]-f[j][l];
			resid*=mask[j][l]/e[j][l];
			//if (fabs(resid)>1e-4)
				//printf("j=%d l=%d resid=%lf\n",j,l,fabs(resid));
			anormf += fabs(resid);
		}
	//printf("%lf\n",anormf);
	for (k=0;k<MAXITS;k++) {
		anorm=0.0;
		jsw=1;
		for (ipass=0;ipass<2;ipass++) {
			lsw=jsw;
			for (j=1;j<n-1;j++) {
				for (l=lsw;l<m-1;l+=2) {
					resid=a[j][l]*u[j+1][l]+b[j][l]*u[j-1][l]
						+c[j][l]*u[j][l+1]+d[j][l]*u[j][l-1]
						+e[j][l]*u[j][l]-f[j][l];
					resid*=mask[j][l]/e[j][l];
					anorm += fabs(resid);
					u[j][l] -= omega*resid;
				}
				lsw=3-lsw;
			}
			jsw=3-jsw;
			omega=(k == 0 && ipass == 0 ? 1.0/(1.0-0.5*rjac*rjac) :
				1.0/(1.0-0.25*rjac*rjac*omega));
		}
		//printf("anorm=%lf\n",anorm);
		if (anorm < EPS*anormf+1e-15) {
			printf("%d steps\n",k);
			return;
		}
	}
	printf("MAXITS exceeded");
}

