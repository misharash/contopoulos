# Based on The Astrophysical Journal 1999 Jan 20 'The axisymmetric pulsar magnetosphere' Ioannis Contopoulos, Demosthenes Kazanas, and Christian Fendt
Now changed to cartesian coordinates, boundary at 10R_LC and radial field condition on it (like in A.N.Timokhin's paper 'On the force-free magnetosphere of aligned rotator', MNRAS, 5 Feb 2008)
* contopoulos1.c -- variant that calculates AA'(Psi) inside the so-called 'elliptic solver' from Numerical Recipes, works not stable enough so is outdated and may be not compatible with some headers
* contopoulos2.c -- variant that is closer to the article: it calculates AA'(Psi) once, puts it into coefficient matrix and gives to elliptic solver, then repeats the same for stable function AA'(Psi) until the difference norm is small enough
* plot1.py -- plot an image for monopole test (like in article)
* plot1a.py -- plot an image for monopole test covering the whole computational domain
* plot1a-evol.py -- plot images for monopole test covering the whole computational domain showing evolution of the solution during recomputation of AA'(x,z)
* lib/ -- parts that are shared between different variants:
	* AA-interp.h -- provides function AA'(Psi) that gets value by linear interpolation from grid (that we get from the light cylinder), for monopole there is an exact function
	* datapair.h -- provides a data type for sorting pairs (Psi,AA) by decresing Psi
	* init.h -- initializes coefficient arrays
	* main-dyn.h -- main loop for variant 1, so is outdated and maybe not compatible with other headers
	* main-stat.h -- main loop for variant 2
	* matops.h -- matrix operations on pure C
	* printtofiles.h -- output results to files
	* solvePSR-dyn.h -- modified Numerical Recipes 'elliptic solver' for variant 1, so is outdated and maybe not compatible with other headers
	* solvePSR-stat.h -- modified Numerical Recipes 'elliptic solver' for variant 2
	* vars.h -- all global variables