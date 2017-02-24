#Based on The Astrophysical Journal 1999 Jan 20 'The axisymmetric pulsar magnetosphere' Ioannis Contopoulos, Demosthenes Kazanas, and Christian Fendt
* contopoulos1.c -- variant that calculates AA'(Psi) inside the so-called 'elliptic solver' from Numerical Recipes
* contopoulos2.c -- variant that is closer to the article: it calculates AA'(Psi) once, puts it into coefficient matrix and gives to elliptic solver, then repeats the same for stable function AA'(Psi) until the difference is small enough
* plot1.py -- plot an image for monopole test (like in article)
* lib/ -- parts that are shared between different variants:
	* AA-interp.h -- provides function AA'(Psi) that gets value by linear interpolation from grid (that we get from the light cylinder)
	* datapair.h -- provides a data type for sorting pairs (Psi,AA) by decresing Psi
	* init.h -- initializes coefficient arrays
	* main-dyn.h -- main loop for variant 1
	* main-stat.h -- main loop for variant 2
	* matops.h -- matrix operations on pure C
	* printtofiles.h -- output results to files
	* solvePSR-dyn.h -- modified Numerical Recipes 'elliptic solver' for variant 1
	* solvePSR-stat.h -- modified Numerical Recipes 'elliptic solver' for variant 2
	* vars.h -- all global variables