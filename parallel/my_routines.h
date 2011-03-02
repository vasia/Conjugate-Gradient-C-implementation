#ifndef MY_ROUTINES_H
#define MY_ROUTINES_H

void matvec(double *Ax, double **Adata, double *x, int n);

void jacobi_precond(double *M, double **Adata, int n);

void psolve(double *Minvx, double **Mdata, double *x, int n);

void generate_randsys(double **A, double *x, double *b, int n);

void generate_randvec(double *vec, int n);

void generate_randmat(double **mat, int n);
	
#endif                          /* CSR_PROBLEM_H */
