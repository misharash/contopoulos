#include <math.h>
#include <stdlib.h>

int AA_init() {
	h=1.0/N;
	int i;
	//creating array for sorting AAs, Psis together
	dp=malloc((N+1)*sizeof(datapair));
	//Take AA' from fig. 1c
	//double Psiopen=1.742;
	AAs=calloc(N+1,sizeof(double));
	Psis=calloc(N+1,sizeof(double));
	for (i=0;i<=N;i++) {
		Psis[i]=(double)(N-i)/N*Psiopen;
		//AAs[i]=2*Psis[i]*(1-Psis[i]/Psiopen)*(2-Psis[i]/Psiopen);
		AAs[i]=2*Psis[i]*(1-Psis[i]/Psiopen)*(2-Psis[i]/Psiopen)*1.1; //distort this
		//take A=-2*Psi => AA'=4*Psi
		//Psis[i]=1e6*(N-i);
		//AAs[i]=4*Psis[i];
	}
	return 0;
}

int inside_init() {
	int i,j;
	double zin,z,xin,x;
	//inside solver init
	inn=N+2;
	inm=N+2;
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
		xin=h*i;
		for (j=1;j<inm-1;j++) {
			zin=h*(j-1);
			ina[i][j]=(1-sqr(xin))/sqr(h)-(1-sqr(xin))/(2*xin*h)-2*xin/(2*h);
			inb[i][j]=(1-sqr(xin))/sqr(h)+(1-sqr(xin))/(2*xin*h)+2*xin/(2*h);
			inc[i][j]=(1-sqr(xin))*(sqr(sqr(1-zin))/sqr(h)-2*sqr(1-zin)*(1-zin)/(2*h));
			ind[i][j]=(1-sqr(xin))*(sqr(sqr(1-zin))/sqr(h)+2*sqr(1-zin)*(1-zin)/(2*h));
			ine[i][j]=-(1-sqr(xin))/sqr(h)*2*(1+sqr(sqr(1-zin)));
			//inf[i][j]=-AA(inu[i][j]); //is to be rewritten every time!
			//initial u is dipole
			x=xin;
			z=1/(1-zin)-1;
			//inu[i][j]=sqr(x)/pow(sqr(x)+sqr(z),1.5);
			inu[i][j]=Psiopen*(1-z/hypot(x,z)); //not dipole but monopole
		}
	}
	//create exact dipolar area here to begin, size of about N/20 (???)
	//assume m=1
	int da=N/5;
	//int da=25;
	for (i=1;i<=da;i++) {
		x=i*h;
		for (j=1;j<=da+1;j++) {
			zin=h*(j-1);
			z=1/(1-zin)-1;
			//inu[i][j]=sqr(x)/pow(sqr(x)+sqr(z),3.0/2.0);
			inu[i][j]=Psiopen*(1-z/hypot(x,z)); //not dipole but monopole
			//printf("%lf %lf %lf\n",x,z,inu[i][j]);
			inmask[i][j]=0;
		}
	}
	//boundary values already zeros, it's okay
	//Neumann boundary condtions
	//bottom -- do it like in notes2.pdf
	j=1;
	zin=h*(j-1);
	for (i=1;i<inn-1;i++) {
		xin=h*i;
		inc[i][j]=(1-sqr(xin))*(sqr(sqr(1-zin))*2/sqr(h));
		ind[i][j]=0;
		inmask[i][j]=0; //temp for monopole check
		inu[i][j]=Psiopen; //temp
		//other coefficients don't change
	}
	//right: this boundary condition is just degradated equation, so we can do nothing but approximate derivative as left derivative
	i=N;
	for (j=1;j<inm-1;j++) {
		ina[i][j]=0;
		inb[i][j]=2/h;
		inc[i][j]=0;
		ind[i][j]=0;
		ine[i][j]=-2/h;
	}
	return 0;
}

int outside_init() {
	int i,j;
	double xout,x,zout,z;
	//outside solver init
	outn=N+2;
	outm=N+1;
	outa=matalloc(outn,outm);
	outb=matalloc(outn,outm);
	outc=matalloc(outn,outm);
	outd=matalloc(outn,outm);
	oute=matalloc(outn,outm);
	outf=matalloc(outn,outm);
	outu=matalloc(outn,outm);
	outmask=matalloc(inn,inm);
	//fill mask
	for (i=0;i<outn;i++)
		for (j=0;j<outm;j++)
			outmask[i][j]=1;
	for (i=1;i<outn-1;i++) {
		xout=h*(i-1);
		for (j=1;j<outm-1;j++) {
			zout=h*j;
			outa[i][j]=(1-1/sqr(1-xout))*(sqr(sqr(1-xout))/sqr(h)-3*sqr(1-xout)*(1-xout)/(2*h))-2*(1-xout)/(2*h);
			outb[i][j]=(1-1/sqr(1-xout))*(sqr(sqr(1-xout))/sqr(h)+3*sqr(1-xout)*(1-xout)/(2*h))+2*(1-xout)/(2*h);
			outc[i][j]=(1-1/sqr(1-xout))*(sqr(sqr(1-zout))/sqr(h)-2*sqr(1-zout)*(1-zout)/(2*h));
			outd[i][j]=(1-1/sqr(1-xout))*(sqr(sqr(1-zout))/sqr(h)+2*sqr(1-zout)*(1-zout)/(2*h));
			oute[i][j]=-(1-1/sqr(1-xout))*2/sqr(h)*(sqr(sqr(1-xout))+sqr(sqr(1-zout)));
			//outf[i][j]=-AA(outu[i][j]); //is to be rewritten every time!
			//initial u is dipole
			x=1/(1-xout);
			z=1/(1-zout)-1;
			//outu[i][j]=sqr(x)/pow(sqr(x)+sqr(z),1.5);
			outu[i][j]=Psiopen*(1-z/hypot(x,z)); //not dipole but monopole
		}
	}
	//upper boundary value is already set as zero
	//lower boundary value will be rewritten after solving equation inside
	//Neumann boundary conditions
	//left is also degradated equation
	i=1;
	for (j=1;j<outm-1;j++) {
		outa[i][j]=-2/h;
		outb[i][j]=0;
		outc[i][j]=0;
		outd[i][j]=0;
		oute[i][j]=2/h;
	}
	//right -- changed to Dirichlet, also will be rewritten after solving equation inside
	return 0;
}
