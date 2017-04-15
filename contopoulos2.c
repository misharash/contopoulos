#include <stdlib.h> //for malloc
#include <signal.h> //SIGINT handling
#include <fenv.h> //for NaN tracking
//for mkdir
#include <sys/stat.h>
#include <sys/types.h>
double sqr(double a) {return a*a;}

#include "lib/datapair.h"
#include "lib/matops.h"
#include "lib/vars.h"
#include "lib/AA-interp.h"
#include "lib/solvePSR-stat.h"
#include "lib/printtofiles.h"
#include "lib/init.h"
#include "lib/main-stat.h"

int main() {
	feenableexcept(FE_INVALID | FE_OVERFLOW); //for NaN tracking
	//malloc
	allfilen=malloc(50);
	tfilen=malloc(50);
	//mkdir
	mkdir("data",0775);
	//call inits
	AA_init();
	all_init();
	//write arrays to file when SIGINT
	signal(SIGINT,printtofiles);
	//main loop
	main_loop();
	//output to file result arrays
	printtofiles(0);
	return 0;
}
