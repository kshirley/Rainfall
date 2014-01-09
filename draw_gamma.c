#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include "kshirley.h"
#include <Rmath.h>
#include <R.h>

void draw_gamma(int *params, double *Sinv_in, double *lhs_in, double *rate_in, double *alpha_in){
	
  int t, j, J=params[0], T=params[1];
	double **Sinv, **lhs, **rate; // These are input from R and will be read back out
	double **z1, **z2, **z3, **temp, alpha=alpha_in[0]; // These variables are for the C code, and stay here
	
	Sinv = dmatrix(J,J);
	fill_dmatrix(Sinv_in,J,J,Sinv);

	lhs = dmatrix(T,J);
	fill_dmatrix(lhs_in,T,J,lhs);
	
	rate = dmatrix(T,1);
	fill_dmatrix(rate_in,T,1,rate);

  // set up objects we'll need:
	z1 = dmatrix(1,J);
  z2 = dmatrix(J,1);
	z3 = dmatrix(1,1);
	temp = dmatrix(1,J);
	
	for (t=0; t<T; t++){		
		dsm(lhs,t,t,0,J-1,temp);
	  dmm(temp,1,J,Sinv,J,J,z1);
		dtp(temp,1,J,z2);
		dmm(z1,1,J,z2,J,1,z3);
		rate[t][0] = (z3[0][0] + alpha)/2;
	}
  
	// read back to R:
	unfill_dmatrix(rate,T,1,rate_in);
	
  // Free memory:
  free_dmatrix(Sinv,J);
  free_dmatrix(lhs,T);
  free_dmatrix(rate,T);
  free_dmatrix(z1,1);
	free_dmatrix(z2,J);
	free_dmatrix(z3,1);
	free_dmatrix(temp,1);
		
}

/* compiling instructions
 R CMD SHLIB draw_gamma.c
*/


