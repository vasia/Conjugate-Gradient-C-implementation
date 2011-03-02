#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cblas.h"

/* Υπολογισμός του dest = a*x + y
 *
 * Ορίσματα:
 *  dest  - Διάνυσμα όπου αποθηκεύεται το αποτέλεσμα.  Μπορεί να έιναι το x ή το y.
 *  a 	  - σταθερά
 *  x, y  - Διανύσματα εισόδου
 *  n     - Μέγεθος διανυσμάτων
 */
void axpy(int n, double a, double *x, int incX, double *y, int incY)
{
    int i;
    for (i = 0; i < n; ++i)
        y[i] = a * x[i] + y[i];
}


/* Υπολογισμός του εσωτερικού γινομένου δύο διανυσμάτων x'*y
 *
 * Ορίσματα:
 *  x, y - Διανύσματα εισόδου
 *  n    - Μέγεθος διανυσμάτων
 */
double ddot(int n, double *x, int incX, double *y, int incY)
{
    int i;
    double final_sum = 0;

    for (i = 0; i < n; ++i)
        final_sum += x[i] * y[i];

    return final_sum;
}

/* Επίλυση του Ax = b χρησιμοποιώντας preconditioned conjugate-gradient.
 * 
 * Ορίσματα:
 *  matvec(Ax, Adata, x, n) - Πολλ/σμός Adata*x στο Ax
 *  psolve(Minvx, Mdata, x, n) - Εφαρμογή του preconditioner
 *  Adata - Ο πίνακας του συστήματος
 *  Mdata - Ο preconditioner
 *  b     - Δεξί μέλος του συστήματος
 *  x     - Αποτέλεσμα
 *  rtol  - Tolerance στο υπόλοιπο (||r||/||b||)
 *  n     - Μέγεθος Συστήματος
 *  maxiter - Μέγιστος αριθμός επαναλήψεων
 *
 * Επιστρέφει:
 *  Αριθμό επαναλήψεων για επιτυχή εκτέλεση, -1 σε περίπτωση αποτυχίας
 */




int precond_cg_blas(void (*matvec_blas) (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const int m, const int n, const double alpha, 	       const double  *Adata, const int lda,const double  *x, const int incX, const double beta, double  *Ax, const int incY),
               void (*psolve) (double *Minvx, double *Mdata, double *x, int n), 
	       double *Adata, double *Mdata, double *b,
	       double *x, double rtol, int n, int maxiter)
{
    const int nbytes = n * sizeof(double);

    double bnorm2;              /* ||b||^2 */
    double rnorm2;              /* Νόρμα υπολοίπου στο τετράγωνο */
    double rz, rzold;           /* r'*z για 2 διαδοχικές επαναλήψεις */
    double alpha, beta;
    double rz_local,rnorm2_local,bnorm2_local;

    double *s;                  /* Κατεύθυνση αναζήτησης */
    double *r;                  /* Υπόλοιπο         */
    double *z;                  /* Προσωρινό διάνυσμα */

    int i = 0,j;                /* Τρέχουσα επανάληψη */

    s =  malloc(nbytes);
    r =  malloc(nbytes);
    z =  malloc(nbytes);

    bnorm2    = ddot(n, b, 1, b, 1);
  
    memset(x, 0, nbytes);	//αρχικοποίηση λύσης
    memcpy(r, b, nbytes);	//και υπολοίπου - r0=b-A*x0 (x0=0)

    psolve(z, Mdata, r, n);	//εφαρμογή του preconditioner - z0 = (M στην -1)*r0s
    memcpy(s, z, nbytes);	//αρχικοποίηση κατεύθυνσης αναζήτησης	- p0 = z0

    /* Αρχικοποίηση rz και rnorm2 */
    rz        = ddot(n, r, 1, z, 1);
    rnorm2    = ddot(n, r, 1, r, 1);

    for (i = 0; i < maxiter ; ++i) {

	matvec_blas(CblasRowMajor, CblasNoTrans, n, n, 1.0, Adata, n, s, 1, 0.0, z, 1);

        
        // Ddot
        alpha = rz / ddot(n, s, 1, z, 1);	//ak = rkT*zk/pkT*A*pk
        axpy(n, alpha, s, 1, x, 1);	//xk+1 = xk + ak*pk
        axpy(n, -alpha, z, 1, r, 1);	//rk+1 = rk - ak*A*pk
  
        psolve(z, Mdata, r, n);		//zk+1 = (M στην -1)*rk+1

        rzold = rz;
        
	rz = ddot(n, r, 1, z, 1);  		//rTk+1*zk+1
        beta = -rz / rzold;		//β = rTk+1*zk+1/rTk*zk 
	axpy(n, -beta-1, s, 1, s, 1);	//pk+1 = zk+1+βk*pk
	axpy(n, 1, z, 1, s, 1);	//pk+1 = zk+1+βk*pk
        
        
        // check error 
	rnorm2     = ddot(n, r, 1, r, 1);
        if(rnorm2 <= bnorm2 * rtol * rtol)
        	break;
        
    
    }
    
    free(z);
    free(r);
    free(s);

    if (i >= maxiter)
        return -1;


    return i;
}

