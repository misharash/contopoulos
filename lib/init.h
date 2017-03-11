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

int inside_init() {
	int i,j;
	double z,x;
	//inside solver init
	inn=N+2;
	inm=10*N+3;
	ina=matalloc(inn,inm);
	inb=matalloc(inn,inm);
	inc=matalloc(inn,inm);
	ind=matalloc(inn,inm);
	ine=matalloc(inn,inm);
	inf=matalloc(inn,inm);
	inu=matalloc(inn,inm);
	inmask=matalloc(inn,inm);
	//fill mask
	for (i=0;i<inn;i++)
		for (j=0;j<inm;j++)
			inmask[i][j]=1;
	//filling initial data inside
	for (i=1;i<inn-1;i++) {
		x=h*i;
		for (j=1;j<inm-1;j++) {
			z=h*(j-1);
			ina[i][j]=(sqr(x)-1)/sqr(h)+(sqr(x)+1)/x/(2*h);
			inb[i][j]=(sqr(x)-1)/sqr(h)-(sqr(x)+1)/x/(2*h);
			inc[i][j]=(sqr(x)-1)/sqr(h);
			ind[i][j]=(sqr(x)-1)/sqr(h);
			ine[i][j]=-4*(sqr(x)-1)/sqr(h);
			//initial u is dipole
			//inu[i][j]=sqr(x)/pow(sqr(x)+sqr(z),1.5);
			inu[i][j]=Psiopen*(1-z/hypot(x,z)); //not dipole but monopole
		}
	}
	//create exact dipolar area here to begin, size of about N/20 (???)
	//assume m=1
	int da=N/5;
	for (i=1;i<=da;i++) {
		x=i*h;
		for (j=1;j<=da+1;j++) {
			z=h*(j-1);
			//inu[i][j]=sqr(x)/pow(sqr(x)+sqr(z),3.0/2.0);
			inu[i][j]=Psiopen*(1-z/hypot(x,z)); //not dipole but monopole
			inmask[i][j]=0;
		}
	}
	//at the axis already zeros, it's okay
	//Neumann boundary condtions
    //top -- radial field
    j=10*N+1;
    z=h*(j-1);
    for (i=1;i<inn-1;i++) {
		x=h*i;
        ina[i][j]-=x/z*inc[i][j];
        inb[i][j]+=x/z*inc[i][j];
        ind[i][j]+=inc[i][j];
		inc[i][j]=0;
        //other don't change
	}
	//bottom -- do it like in notes2.pdf
	j=1;
	z=h*(j-1);
	for (i=1;i<inn-1;i++) {
		x=h*i;
		inc[i][j]+=ind[i][j];
		ind[i][j]=0;
		inmask[i][j]=0; //temp for monopole check
		inu[i][j]=Psiopen; //temp
		//other coefficients don't change
	}
	//right: this boundary condition is just degradated equation, so we can do nothing but approximate derivative as left derivative
	i=N;
	for (j=1;j<inm-1;j++) {
		ina[i][j]=0;
		inb[i][j]=-2/h;
		inc[i][j]=0;
		ind[i][j]=0;
		ine[i][j]=2/h;
	}
	return 0;
}

int outside_init() {
	int i,j;
	double x,z;
	//outside solver init
	outn=9*N+3;
	outm=10*N+2;
	outa=matalloc(outn,outm);
	outb=matalloc(outn,outm);
	outc=matalloc(outn,outm);
	outd=matalloc(outn,outm);
	oute=matalloc(outn,outm);
	outf=matalloc(outn,outm);
	outu=matalloc(outn,outm);
	outmask=matalloc(outn,outm);
	//fill mask
	for (i=0;i<outn;i++)
		for (j=0;j<outm;j++)
			outmask[i][j]=1;
	for (i=1;i<outn-1;i++) {
		x=1.+h*(i-1);
		for (j=1;j<outm-1;j++) {
			z=h*j;
			outa[i][j]=(sqr(x)-1)/sqr(h)+(sqr(x)+1)/x/(2*h);
			outb[i][j]=(sqr(x)-1)/sqr(h)-(sqr(x)+1)/x/(2*h);
			outc[i][j]=(sqr(x)-1)/sqr(h);
			outd[i][j]=(sqr(x)-1)/sqr(h);
			oute[i][j]=-4*(sqr(x)-1)/sqr(h);
			//initial u is dipole
			//outu[i][j]=sqr(x)/pow(sqr(x)+sqr(z),1.5);
			outu[i][j]=Psiopen*(1-z/hypot(x,z)); //not dipole but monopole
		}
	}
	//lower boundary value will be rewritten after solving equation inside
	//Neumann boundary conditions
    //top -- radial field
    j=10*N;
    z=h*j;
    for (i=1;i<outn-1;i++) {
		x=h*(i-1);
        outa[i][j]-=x/z*outc[i][j];
        outb[i][j]+=x/z*outc[i][j];
        outd[i][j]+=outc[i][j];
		outc[i][j]=0;
        //other don't change
	}
	//right -- radial field
    i=outn-2;
    x=h*(i-1);
    for (j=1;j<outm-1;j++) {
		z=h*j;
        outd[i][j]-=z/x*outa[i][j];
        outc[i][j]+=z/x*outa[i][j];
        outb[i][j]+=outa[i][j];
		outa[i][j]=0;
        //other don't change
	}
	//left is also degradated equation
	i=1;
	for (j=1;j<outm-1;j++) {
		outa[i][j]=2/h;
		outb[i][j]=0;
		outc[i][j]=0;
		outd[i][j]=0;
		oute[i][j]=-2/h;
	}
    //top right angle is special
    i=outn-2;
    j=outm-2;
    outa[i][j]=outc[i][j]=0;
    outb[i][j]=outd[i][j]=-1;
    oute[i][j]=2;
	return 0;
}
