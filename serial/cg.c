#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Υπολογισμός του dest = a*x + y
 *
 * Ορίσματα:
 *  dest  - Διάνυσμα όπου αποθηκεύεται το αποτέλεσμα.  Μπορεί να έιναι το x ή το y.
 *  a 	  - σταθερά
 *  x, y  - Διανύσματα εισόδου
 *  n     - Μέγεθος διανυσμάτων
 */
void axpy(double *dest, double a, double *x, double *y, int n)
{
    int i;
    for (i = 0; i < n; ++i)
        dest[i] = a * x[i] + y[i];
}


/* Υπολογισμός του εσωτερικού γινομένου δύο διανυσμάτων x'*y
 *
 * Ορίσματα:
 *  x, y - Διανύσματα εισόδου
 *  n    - Μέγεθος διανυσμάτων
 */
double ddot(double *x, double *y, int n)
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
int precond_cg(void (*matvec) (double *Ax, double **Adata, double *x, int n),
               void (*psolve) (double *Minvx, double *Mdata, double *x, int n), 
	       double **Adata, double *Mdata, double *b,
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

    s = (double *) malloc(nbytes);
    r = (double *) malloc(nbytes);
    z = (double *) malloc(nbytes);

 
    bnorm2    = ddot(b, b, n);
    
    memset(x, 0, nbytes);	//αρχικοποίηση λύσης
    memcpy(r, b, nbytes);	//και υπολοίπου - r0=b-A*x0 (x0=0)

    psolve(z, Mdata, r, n);	//εφαρμογή του preconditioner - z0 = (M στην -1)*r0

    memcpy(s, z, nbytes);	//αρχικοποίηση κατεύθυνσης αναζήτησης	- p0 = z0

    /* Αρχικοποίηση rz και rnorm2 */
    rz        = ddot(r, z, n);
    rnorm2    = ddot(r, r, n);
    
        
    printf("rz=%2.15f,alpha=%2.15f,rnorm=%2.15f,bnorm=%f\n",rz,alpha,rnorm2,bnorm2);

    for (i = 0; i < maxiter ; ++i) {
        printf("Επανάληψη #%d\n", i);

        matvec(z, Adata, s, n);	//z:=A*pk 

        
        /* Ddot*/
        alpha = rz / ddot(s, z, n);	//ak = rkT*zk/pkT*A*pk
        axpy(x, alpha, s, x, n);	//xk+1 = xk + ak*pk
        axpy(r, -alpha, z, r, n);	//rk+1 = rk - ak*A*pk
  
        psolve(z, Mdata, r, n);		//zk+1 = (M στην -1)*rk+1

        rzold = rz;
        

        rz = ddot(r, z, n);  		//rTk+1*zk+1
        beta = -rz / rzold;		//β = rTk+1*zk+1/rTk*zk 
        axpy(s, -beta, s, z, n);	//pk+1 = zk+1+βk*pk

        printf("(%d)rz=%2.15f,alpha=%2.15f,rnorm=%2.15f\n",i,rz,alpha,rnorm2);
        
        
        /* check error */
        rnorm2     = ddot(r, r, n);
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

