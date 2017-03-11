#include <math.h>

double AA(double Psi) {
	double dPsi=Psi-Psiopen;
	if (Psi>Psiopen) return -Psiopen*Psiopen*dPsi/sigma2*exp(-dPsi*dPsi/sigma2);
	if (Psi<0) return 0;
	return 2*Psi*(1-Psi/Psiopen)*(2-Psi/Psiopen);
	//return 2*Psi*(1-Psi/Psiopen)*(2-Psi/Psiopen)/1.05; //distort this
}
/*
double AA(double Psi) {
	double dPsi=Psi-Psis[0];
	if (dPsi>0)
		return AAs[0]*exp(-dPsi*dPsi/sigma2);
	int i;
	for (i=1;(i<=N)&&(Psis[i]>=Psi);i++);
	if (fabs(Psis[i]-Psis[i-1])>1e-5)
		return AAs[i-1]+(AAs[i]-AAs[i-1])*(Psi-Psis[i-1])/(Psis[i]-Psis[i-1]);
	return AAs[i];
}
*/
