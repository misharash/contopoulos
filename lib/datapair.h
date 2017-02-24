typedef struct _datapair {
	double Psi,AA;
} datapair;

int dpcmp(const void *a,const void *b) { //to qsort datapairs by _de_creasing Psi
	datapair *aa=(datapair*)a,*bb=(datapair*)b;
	if (aa->Psi==bb->Psi)
		return 0;
	if (aa->Psi>bb->Psi)
		return -1;
	return 1;
}
