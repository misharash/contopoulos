//global variables
int N=80; //grid size
double mu1=0.5,mu2=0.5,mu3=-.4; //weight coefficients for AA
double l1=0.5,l2=0.5; //weight coefficients for Psi
double *AAs,*Psis; //table of (Psy,AA') on light cylinder
double A0; //value of A on last open field line
double sigma2=0.01; //dispersion of gaussian
double Psiopen=1.742; //Psiopen for monopole current distribution
double **inu,**outu; //data arrays
int inm,outm,inn,outn,k; //limits and counters
char *infilen,*outfilen,*tfilen; //file names to write
double h; //just 1/N
datapair *dp; //array for sorting Psi and AA together
double **ina,**inb,**inc,**ind,**ine,**inf,**inmask; //coefficient arrays for inside solver
double **outa,**outb,**outc,**outd,**oute,**outf,**outmask; //coefficient arrays for outside solver
double ewti=0.3,ewto=0.04;
