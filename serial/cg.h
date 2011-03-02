#ifndef CG_H
#define CG_H

#define	MAX_NPROCS	100

void axpy(double *dest, double alpha, double *x, double *y, int n);

double ddot(double *x, double *y, int n);

int precond_cg(void (*matvec) (double *Ax, double **Adata, double *x, int n),
               void (*psolve) (double *Minvx, double *Mdata, double *x,
                               int n), double **Adata, double *Mdata, double *b,
               double *x, double rtol, int n, double *rhist, int maxiter);

void dummy_psolve(double *Minvx, void *Mdata, double *x, int n);

#endif                          /* CG_H */
