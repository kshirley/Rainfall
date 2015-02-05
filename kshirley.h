/*
 *  kshirley.h
 *  
 *
 *  Created by Kenneth Shirley on 8/22/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>

double dmean();
double imean();
void dsum();
void isum();

int *ivector();
int **imatrix();
int ***i3array();
double *dvector();
double **dmatrix();
double ***d3array();
double ****d4array();

void free_ivector();
void free_imatrix();
void free_i3array();
void free_dvector();
void free_dmatrix();
void free_d3array();
void free_d4array();

void fill_ivector();
void fill_dvector();
void fill_imatrix();
void fill_dmatrix();
void fill_i3array();
void fill_d3array();
void fill_d4array();
void unfill_dmatrix();
void unfill_d4array();
void unfill_i3array();
void unfill_imatrix();

double unif_rand();
double norm_rand();
double exp_rand();
void dmm();
void imm();
void dsm();
void ism();
void dtp();
void itp();
void dtsq();

int sample_old();
void sample();

// function to calculate mean of doubles
double dmean(double *x, int first, int last){
  double s=0.0, m=0.0;
  int i, n;
  n = last - first + 1;
  for (i=first; i<last+1; i++) s += x[i];
  m = s/n;
  return m;
}

// function to calculate mean of integers
double imean(int *x, int first, int last){
  double s=0.0, m=0.0;
  int i, n;
  n = last - first + 1;
  for (i=first; i<last+1; i++) s += x[i];
  m = s/n;
  return m;
}



double *dvector(int n){
  double *v;
  v = (double *)malloc(n * sizeof(double));
  return v;
}

int *ivector(int n){
  int *v;
  v = (int *)malloc(n * sizeof(int));
  return v;
}

double **dmatrix(int nrow, int ncol){
  double **mat;
  int i;
  mat = (double **)malloc(nrow * sizeof(double *));
  for (i=0; i<nrow; i++){
    mat[i] = (double *)malloc(ncol * sizeof(double));
  }
  return mat;
}

int **imatrix(int nrow, int ncol){
  int **mat;
  int i;
  mat = (int **)malloc(nrow * sizeof(int *));
  for (i=0; i<nrow; i++){
    mat[i] = (int *)malloc(ncol * sizeof(int));
  }
  return mat;
}

int ***i3array(int nrow, int ncol, int ndep){
  int ***mat;
  int i, j;
  mat = (int ***)malloc(nrow * sizeof(int **));
  for (i=0; i<nrow; i++){
    mat[i] = (int **)malloc(ncol * sizeof(int *));
    for (j=0; j<ncol; j++){
      mat[i][j] = (int *)malloc(ndep * sizeof(int));
    }
  }
  return mat;
}

double ***d3array(int nrow, int ncol, int ndep){
  double ***mat;
  int i, j;
  mat = (double ***)malloc(nrow * sizeof(double **));
  for (i=0; i<nrow; i++){
    mat[i] = (double **)malloc(ncol * sizeof(double *));
    for (j=0; j<ncol; j++){
      mat[i][j] = (double *)malloc(ndep * sizeof(double));
    }
  }
  return mat;
}

double ****d4array(int nrow, int ncol, int ndep, int nup){
  double ****mat;
  int i, j, k;
  mat = (double ****)malloc(nrow * sizeof(double ***));
  for (i=0; i<nrow; i++){
    mat[i] = (double ***)malloc(ncol * sizeof(double **));
    for (j=0; j<ncol; j++){
      mat[i][j] = (double **)malloc(ndep * sizeof(double *));
      for (k=0; k<ndep; k++){
        mat[i][j][k] = (double *)malloc(nup * sizeof(double));
      }
    }
  }
  return mat;
}

// Free a dvector
void free_dvector(double *v){
  free(v);
}

// Free a dmatrix
void free_dmatrix(double **v, int nrow){
  int i;
  for (i=0; i<nrow; i++) free(v[i]);
  free(v);
}

// Free a d3array
void free_d3array(double ***v, int nrow, int ncol){
  int i, j;
  for (i=0; i<nrow; i++){
    for (j=0; j<ncol; j++){
      free(v[i][j]);
    }
    free(v[i]);
  }
  free(v);
}

// Free an i3array
void free_i3array(int ***v, int nrow, int ncol){
  int i, j;
  for (i=0; i<nrow; i++){
    for (j=0; j<ncol; j++){
      free(v[i][j]);
    }
    free(v[i]);
  }
  free(v);
}

// Free a d4array
void free_d4array(double ****v, int nrow, int ncol, int ndep){
  int i, j, k;
  for (i=0; i<nrow; i++){
    for (j=0; j<ncol; j++){
      for (k=0; k<ndep; k++){
        free(v[i][j][k]);
      }
      free(v[i][j]);
    }
    free(v[i]);
  }
  free(v);
}

// Free an ivector
void free_ivector(int *v){
  free(v);
}

// Free an imatrix
void free_imatrix(int **v, int nrow){
  int i;
  for (i=0; i<nrow; i++) free(v[i]);
  free(v);
}

// Random Uniform(0,1)
double unif_rand(){
  return drand48();
}

// Random Normal(0,1)
double norm_rand(){
  double u1, u2, z;
  u1 = drand48();
  u2 = drand48();
  z = sqrt(-2.0*log(u1)) * cos(2.0*M_PI*u2);  
  return(z);
}

// Random Exponential(1)
double exp_rand(){
  double u;
  u = drand48();
  return(-log(u));
}

// Double matrix multiplier
void dmm(double **A, int n1, int m1, double **B, int n2, int m2, double **C){
  int i, j, k;
  double s;
  
  for (i=0; i<n1; i++){
    for (j=0; j<m2; j++){
      s = 0.0;
      for (k=0; k<m1; k++){
        s += A[i][k]*B[k][j];}
      C[i][j] = s;}}
}

// Integer matrix multiplier
void imm(int **A, int n1, int m1, int **B, int n2, int m2, int **C){
  int i, j, k;
  int s;
  
  for (i=0; i<n1; i++){
    for (j=0; j<m2; j++){
      s = 0;
      for (k=0; k<m1; k++){
        s += A[i][k]*B[k][j];}
      C[i][j] = s;}}
}

// Sub double matrix creator:
void dsm(double **A, int n1, int n2, int m1, int m2, double **C){
  int i, j;
  int nrow=n2-n1+1, ncol=m2-m1+1;
  
  for (i=0; i<nrow; i++){
    for (j=0; j<ncol; j++){
      C[i][j] = A[n1+i][m1+j];}}
}

// Sub integer matrix creator:
void ism(int **A, int n1, int n2, int m1, int m2, int **C){
  int i, j;
  int nrow=n2-n1+1, ncol=m2-m1+1;
  
  for (i=0; i<nrow; i++){
    for (j=0; j<ncol; j++){
      C[i][j] = A[n1+i][m1+j];}}
}

// Transpose of a double matrix:
void dtp(double **X, int n, int m, double **Y){
  int i, j;
  
  for (i=0; i<n; i++){
    for (j=0; j<m; j++){
      Y[j][i] = X[i][j];
    }
  }
}

// Transpose of a square matrix of doubles, replacing the original matrix:
void dtsq(double **X, int n){
  int i, j;
  double temp;
	
  for (i=0; i<n; i++){
    for (j=i+1; j<n; j++){
      temp = X[i][j];
      X[i][j] = X[j][i];
      X[j][i] = temp;
    }
  }
}


// Transpose of an integer matrix:
void itp(int **X, int n, int m, int **Y){
  int i, j;
  
  for (i=0; i<n; i++){
    for (j=0; j<m; j++){
      Y[j][i] = X[i][j];}}
}

// Fill a double matrix with an input vector:
void fill_dmatrix(double *input, int n, int m, double **output){
  int i,j;
  for (i=0; i<n; i++){
    for (j=0; j<m; j++){
      output[i][j] = input[j*n + i];
    }
  }
}

// Unfill a double matrix into an output vector:
void unfill_dmatrix(double **input, int n, int m, double *output){
  int i,j;
  for (i=0; i<n; i++){
    for (j=0; j<m; j++){
      output[j*n + i] = input[i][j];
    }
  }
}

// Unfill an integer matrix into an output vector:
void unfill_imatrix(int **input, int n, int m, int *output){
  int i,j;
  for (i=0; i<n; i++){
    for (j=0; j<m; j++){
      output[j*n + i] = input[i][j];
    }
  }
}


// Fill an integer matrix with an input vector:
void fill_imatrix(int *input, int n, int m, int **output){
  int i,j;
  for (i=0; i<n; i++){
    for (j=0; j<m; j++){
      output[i][j] = input[j*n + i];
    }
  }
}

// Fill an integer i3array with an input vector:
void fill_i3array(int *input, int n, int m, int l, int ***output){
  int i,j,k;
  for (i=0; i<n; i++){
    for (j=0; j<m; j++){
      for (k=0; k<l; k++){
	output[i][j][k] = input[k*m*n + j*n + i];
      }
    }
  }
}

// Fill a double vector with an input vector:
void fill_dvector(double *input, int n, double *output){
  int i;
  for (i=0; i<n; i++) output[i] = input[i];
}

// Fill a double vector with an input vector:
void fill_ivector(int *input, int n, int *output){
  int i;
  for (i=0; i<n; i++) output[i] = input[i];
}

// Fill a double d3array with an input vector:
void fill_d3array(double *input, int n, int m, int l, double ***output){
  int i,j,k;
  for (i=0; i<n; i++){
    for (j=0; j<m; j++){
      for (k=0; k<l; k++){
	output[i][j][k] = input[k*m*n + j*n + i];
      }
    }
  }
}

// Fill a double d4array with an input vector:
void fill_d4array(double *input, int n, int m, int r, int s, double ****output){
  int i,j,k,l;
  for (i=0; i<n; i++){
    for (j=0; j<m; j++){
      for (k=0; k<r; k++){
	for (l=0; l<s; l++){
	  output[i][j][k][l] = input[l*r*m*n + k*m*n + j*n + i];
	}
      }
    }
  }
}

// Unfill a double d4array into an output vector:
void unfill_d4array(double ****input, int n, int m, int r, int s, double *output){
  int i,j,k,l;
  for (i=0; i<n; i++){
    for (j=0; j<m; j++){
      for (k=0; k<r; k++){
	for (l=0; l<s; l++){
	  output[l*r*m*n + k*m*n + j*n + i] = input[i][j][k][l];
	}
      }
    }
  }
}

// Unfill an integer i3array into an output vector:
void unfill_i3array(int ***input, int n, int m, int r, int *output){
  int i,j,k;
  for (i=0; i<n; i++){
    for (j=0; j<m; j++){
      for (k=0; k<r; k++){
	output[k*m*n + j*n + i] = input[i][j][k];
      }
    }
  }
}


// Sum of a double vector:
void dsum(double *x, int n, double *z){
  int i;
  for (i=0; i<n; i++) z[0] += x[i];
}

// Sum of an integer vector:
void isum(int *x, int n, int *z){
  int i;
  for (i=0; i<n; i++) z[0] += x[i];
}

// Replicate the Sample function from R, with replacement always:
/* int sample_old(int first, int last, double *prob){ */
/*   int j, l = last - first + 1; */
/*   int z; */
/*   double x, s;   */
  
/*   x = drand48(); */
/*   s = 0.0; */
/*   for (j = 0; j < l; j++){ */
/*     s += prob[j]; */
/*     if (x < s){ */
/*       z = j + first; */
/*       break; */
/*     } */
/*   } */
/*   return z; */
/* } */


// Replicate the Sample function from R, with replacement always:
void sample(int first, int last, double *prob, int *z){
  int j, l=last-first+1;
  double x, s;
  
  x = drand48();
  s = 0.0;
  for (j=0; j<l; j++){
    s += prob[j];
    if (x<s){
      z[0] = j+first;
      break;
    }
  }
}

/*
// Generate a random dirichlet draw:
void rdirichlet(int n, double *alpha, double *y){
int i;
double y_sum=0.0, *z;	
z = dvector(n);
		
// random sample from gamma distribution:
for (i=0; i<n; i++){
y[i] = rgamma(alpha[i],1);
y_sum += y[i];
}	
for (i=0; i<n; i++){
z[i] = y[i]/y_sum;
y[i] = z[i];
}
	
free_dvector(z);
}
*/





