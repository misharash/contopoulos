#include <math.h>
#include <stdlib.h>

int AA_init() {
    h=1.0/N;
    int i;
    //creating array for sorting AAs, Psis together
    dp=malloc((10*N+1)*sizeof(datapair));
    //Take AA' from fig. 1c
    //double Psiopen=1.742;
    AAs=calloc(10*N+1,sizeof(double));
    Psis=calloc(10*N+1,sizeof(double));
    for (i=0;i<=10*N;i++) {
	Psis[i]=(double)(10*N-i)/(10*N)*Psiopen;
	//AAs[i]=2*Psis[i]*(1-Psis[i]/Psiopen)*(2-Psis[i]/Psiopen);
	AAs[i]=2*Psis[i]*(1-Psis[i]/Psiopen)*(2-Psis[i]/Psiopen)*1.1; //distort this
    }
    return 0;
}

int all_init() {
    int i,j;
    double z,x;
    //solver init
    alln=10*N+2;
    allm=10*N+3;
    alla=matalloc(alln,allm);
    allb=matalloc(alln,allm);
    allc=matalloc(alln,allm);
    alld=matalloc(alln,allm);
    alle=matalloc(alln,allm);
    allf=matalloc(alln,allm);
    allu=matalloc(alln,allm);
    allmask=matalloc(alln,allm);
    //fill mask
    for (i=1;i<alln-1;i++)
	    for (j=1;j<allm-1;j++)
		    allmask[i][j]=1;
    //filling initial data
    for (i=1;i<alln-1;i++) {
	x=h*i;
	for (j=1;j<allm-1;j++) {
	    z=h*(j-1);
	    alla[i][j]=(sqr(x)-1)/sqr(h)+(sqr(x)+1)/x/(2*h);
	    allb[i][j]=(sqr(x)-1)/sqr(h)-(sqr(x)+1)/x/(2*h);
	    allc[i][j]=(sqr(x)-1)/sqr(h);
	    alld[i][j]=(sqr(x)-1)/sqr(h);
	    alle[i][j]=-4*(sqr(x)-1)/sqr(h);
	    //initial u is dipole
	    //allu[i][j]=sqr(x)/pow(sqr(x)+sqr(z),1.5);
	    //allu[i][j]=Psiopen*(1-z/hypot(x,z)); //not dipole but monopole
	    //allu[i][j]=0.8*Psiopen*(1-z/hypot(x,z)); //distorted monopole
	}
    }
    //create exact dipolar area here to begin, size of about N/5
    //assume m=1
    int da=N/5;
    for (i=1;i<=da;i++) {
	x=i*h;
	for (j=1;j<=da+1;j++) {
	    z=h*(j-1);
	    //inu[i][j]=sqr(x)/pow(sqr(x)+sqr(z),3.0/2.0);
	    allu[i][j]=Psiopen*(1-z/hypot(x,z)); //not dipole but monopole
	    allmask[i][j]=0;
	}
    }
    //Timokhin's form of equation on LC
    i=N;
    x=1;
    for (j=1;j<allm-1;j++) {
        alla[i][j]=allb[i][j]=4/sqr(h);
        allc[i][j]=alld[i][j]=2/sqr(h);
        alle[i][j]=-12/sqr(h);
    }
	//at the Z axis already zeros, it's okay
	//Neumann boundary condtions
    //top -- radial field
    j=10*N+1;
    z=h*(j-1);
    for (i=1;i<alln-1;i++) {
	x=h*i;
        alla[i][j]-=x/z*allc[i][j];
        allb[i][j]+=x/z*allc[i][j];
        alld[i][j]+=allc[i][j];
	allc[i][j]=0;
        //other don't change
    }
    //bottom inside -- dPsi/dZ=0
    j=1;
    z=h*(j-1);
    for (i=1;i<=N;i++) {
	x=h*i;
	allc[i][j]+=alld[i][j];
	alld[i][j]=0;
	allmask[i][j]=0; //temp for monopole check
	allu[i][j]=Psiopen; //temp
	//other coefficients don't change
    }
    //bottom outside -- now (monopole) just Psiopen, then will think
    for (i=N+1;i<alln-1;i++) {
        allu[i][j]=Psiopen;
        allmask[i][j]=0;
    }
    //right -- radial field
    i=alln-2;
    x=h*i;
    for (j=1;j<allm-1;j++) {
	z=h*j;
        allc[i][j]-=z/x*alla[i][j];
        alld[i][j]+=z/x*alla[i][j];
        allb[i][j]+=alla[i][j];
	alla[i][j]=0;
        //other don't change
	}
    //top right angle is special
    i=alln-2;
    j=allm-2;
    alla[i][j]=allc[i][j]=0;
    allb[i][j]=alld[i][j]=-1;
    alle[i][j]=2;
    return 0;
}
