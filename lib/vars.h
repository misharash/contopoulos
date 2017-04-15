//global variables
int N=30; //base grid size (for 1R_LC)
double *AAs,*Psis; //table of (Psy,AA') on light cylinder
double sigma2=0.01; //dispersion of gaussian
double Psiopen=1.742; //Psiopen for monopole current distribution
double **allu; //data array
int alln,allm,k; //limits and counters
char *allfilen,*tfilen; //file names to write
double h; //just 1/N
datapair *dp; //array for sorting Psi and AA together
double **alla,**allb,**allc,**alld,**alle,**allf,**allmask; //coefficient arrays for solver
double ewt=0.25; //evolution weight
