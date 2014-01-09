#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include "kshirley.h"
#include <Rmath.h>
#include <R.h>

void chol_lower();
void chol_inverse();

void zdraw(int *params, double *Smu_in, double *Pxb_in, double *Sig_in, double *Prior_in){
	
  int i, j, s, t, S=params[0], T=params[1];
	double **Smu, **Pxb, **Sig, **Prior; // These are input from R and will be read back out
	double **Sigma, *p, **B, **C, **z_right, **z_mean, **z, **x; // Just variables for the C code
	
	GetRNGstate();	

	Smu = dmatrix(S,T);
	fill_dmatrix(Smu_in,S,T,Smu);

	Pxb = dmatrix(S,T);
	fill_dmatrix(Pxb_in,S,T,Pxb);

	Sig = dmatrix(S,T);
	fill_dmatrix(Sig_in,S,T,Sig);

	Prior = dmatrix(S,S);
	fill_dmatrix(Prior_in,S,S,Prior);

	Sigma = dmatrix(S,S); // Used to create covariance matrix for each time point
	p = dvector(S); // empty vector of length S for cholesky decomposition diagonal
	B = dmatrix(S,S); // placeholder for transpose in cholesky inverse
	C = dmatrix(S,S); // placeholder for result of cholesky inverse
	z_right = dmatrix(S,1); // right-hand matrix for the mean vector for each time point
	z_mean = dmatrix(S,1); // mean vector for each time point
	z = dmatrix(S,1); // set up the vector of length S to generate independent random normals  
	x = dmatrix(S,1);	// place holder for matrix of correlated random normal draws

	
	for (t=0; t<T; t++){

		// Construct the covariance matrix for this value of t from Prior (S x S) and Sig (S x T)
		// Note: Each column of Sig will form a diagonal (S x S) matrix that is added to Prior
		for (i=0; i<S; i++){
			for (j=i; j<S; j++){
				if (i == j) Sigma[i][j] = Prior[i][j] + Sig[i][t];
				else {
				  Sigma[i][j] = Prior[i][j];
					Sigma[j][i] = Prior[j][i];
				}
			}
		}				
				
    // C becomes inverse of Sigma, via Cholesky decomposition:
	  chol_inverse(Sigma,S,p,B,C); // Now Sigma is changed to L^(-1), where L is the lower cholesky decomp of Sigma
	  
    // compute the right-hand (S x 1) matrix in the mean of z_t:
		for (s=0; s<S; s++) z_right[s][0] = Smu[s][t] + Pxb[s][t];
		
    // compute the mean of z_t
		dmm(C,S,S,z_right,S,1,z_mean);
		
		// Now C is the lower triangle Cholesky decomp of the original Sigma^(-1)
		chol_lower(C,S,p);
	
    // draw S independent normal RVs:
	  for (s=0; s<S; s++) z[s][0] = rnorm(0,1);
	
	  // Multiply the lower triangular cholesky square root by the independent normal vector:
	  dmm(C,S,S,z,S,1,x); // x is the (S x 1) matrix of correlated random normals
		
		// Write over one of the mean matrices to send output back to R:
		for (s=0; s<S; s++) Smu[s][t] = x[s][0] + z_mean[s][0]; // Smu will hold the (S x T) matrix of new draws
	  
	}

  
	// read the result back to input vector z_in:
	unfill_dmatrix(Smu,S,T,Smu_in);

  // Free memory:
	free_dmatrix(Smu,S);
	free_dmatrix(Pxb,S);
	free_dmatrix(Sig,S);
	free_dmatrix(Prior,S);
	free_dmatrix(Sigma,S);
	free_dvector(p);
  free_dmatrix(B,S);
	free_dmatrix(C,S);
  free_dmatrix(z_right,S);
	free_dmatrix(z_mean,S);
  free_dmatrix(z,S);
  free_dmatrix(x,S);
	
	PutRNGstate();
	
}

/* compiling instructions
 R CMD SHLIB zdraw.c
*/




void chol_lower(double **a, int n, double *p){
  
	int i, j, k;
	double sum=0.0;

	// Cholesky decomposition of Sigma:
	for (i=1; i<=n; i++){
		for (j=i; j<=n; j++){
			for (sum=a[i-1][j-1], k=i-1; k>=1; k--) sum -= a[i-1][k-1]*a[j-1][k-1];
			if (i == j){
				p[i-1]=sqrt(sum);
			} else a[j-1][i-1]=sum/p[i-1];
		}
	}
	
  // Get lower triangular matrix from Cholesky decomposition:
	for (i=0; i<n; i++){
		for (j=i; j<n; j++){
			if (i == j) a[i][j] = p[i];
			else a[i][j] = 0.0;
		}
	}
	
}


// inputs b and c are placeholders, a is inverted, and the output is c:
void chol_inverse(double **a, int n, double *p, double **b, double **c){
  
	int i, j, k;
	double sum=0.0;
	
	// Cholesky decomposition of Sigma:
	for (i=1; i<=n; i++){
		for (j=i; j<=n; j++){
			for (sum=a[i-1][j-1], k=i-1; k>=1; k--) sum -= a[i-1][k-1]*a[j-1][k-1];
			if (i == j){
				p[i-1]=sqrt(sum);
			} else a[j-1][i-1]=sum/p[i-1];
		}
	}
	
  // Invert the lower triangle:
  for (i=1; i<=n; i++){
		a[i-1][i-1] = 1.0/p[i-1];
		for (j=i+1; j<=n; j++){
			sum = 0.0;
			for (k=i; k<j; k++) sum -= a[j-1][k-1]*a[k-1][i-1];
			a[j-1][i-1] = sum/p[j-1];
		}
	}
  
  // set upper triangular elements to zero:
	for (i=0; i<n; i++){
		for (j=i+1; j<n; j++) a[i][j] = 0.0;
	}
	
	dtp(a,n,n,b); // b becomes transpose of a
	dmm(b,n,n,a,n,n,c);	// a-transpose %*% a = c, where c is the inverse of a
}


