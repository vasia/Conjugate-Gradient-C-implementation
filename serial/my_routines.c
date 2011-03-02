#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif

/* Πολλαπλασιασμός πίνακα επί διάνυσμα
 */
void matvec(double *Ax, double **Adata, double *xvect, int n)
{
    int i, j;		
    
    for (i = 0; i < n; ++i) {
        Ax[i] = 0;
        for (j = 0; j < n; ++j) {
            Ax[i] += Adata[i][j]*xvect[j];
        }
    }
}


/* Δημιουργία ενός Jacobi preconditioner - Τα μηδενικά δεν αποθηκεύονται 
 * Mij = Aij, i=j
 * Mij = 0, i<>j
 */
void jacobi_precond(double *M, double **Adata, int n)
{
    int i;		
    
    for (i = 0; i < n; ++i) {
           M[i] += Adata[i][i];
    }
}

/* Εφαρμογή του Jacobi preconditioner
 * 
 */
void psolve(double *Minvx, double *Mdata, double *x, int n)
{
    int i;		
    
    for(i=0; i<n; i++){
	Minvx[i] = 1/Mdata[i]*x[i];
    }
}

/* Παραγωγή ενός τυχαίου διανύσματος vec μεγέθους n - solution vector
 * 
 */
void generate_randvec(double *vec, int n)
{
    int i;	

    for(i=0; i<n; i++){
	vec[i] = (double)random()/(RAND_MAX+1.0);
    }
	
}

/* Παραγωγή ενός τυχαίου συμμετρικού θετικά ορισμένου πίνακα A μεγέθους n
 * 
 * Παράγεται ένας τυχαίος πίνακας C και πολ/ζεται με τον ανάστροφό του ώστε
 * να δώσει ως αποτέλεσμα συμμετρικό θετικά ορισμένο πίνακα
 */
void generate_randmat(double **mat, int n)
{
    int i, j, k;	
    double **randmat;
    double **mattrans;	

    randmat = malloc(n * sizeof(double*));
 		if( randmat == NULL)  	    printf("randmat--Out of memory");  
	  	  
		randmat[0] = malloc(n * n *sizeof(double));
		for(i = 1; i < n; i++)
			randmat[i] = randmat[0] + i * n;	
	
     mattrans = malloc(n * sizeof(double*));
 		if( mattrans == NULL)  	    printf("mattrans--Out of memory");  
	  	  
		mattrans[0] = malloc(n * n *sizeof(double));
		for(i = 1; i < n; i++)
			mattrans[i] = mattrans[0] + i * n;		

    for(i=0; i<n; i++){
	for(j=0; j<n; j++){
		randmat[i][j] = (double)random()/(RAND_MAX+1.0);
    	}
    }

    for(i=0; i<n; i++){
	for(j=0; j<n; j++){
		mattrans[i][j] = randmat[j][i];
    	}
    }

    for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			for(k=0;k<n;k++){
				mat[i][k] += mattrans[i][j]*randmat[j][k];
			}
		}
    }

    //free memory
    free(randmat[0]); free(randmat);
    free(mattrans[0]); free(mattrans);
	
}

/* Παραγωγή ενός τυχαίου συστήματος A*x=b μεγέθους n
 * 
 */
void generate_randsys(double **A, double *x, double *b, int n)
{
    int i, j;

    generate_randmat(A, n);
    generate_randvec(x, n);	

    for(i=0;i<n;i++)
	    for(j=0;j<n;j++)
    		b[i]+= A[i][j]*x[j];
}

/*int main(int argc, char **argv)
{
	int i, j;	
	double **Adata;	//Ο πίνακας του συστήματος
	double *x;	//Η λύση του συστήματος (που παράγεται για το τυχαίο σύστημα)
	double *cg_x;	//Η λύση του συστήματος που προκύπτει από τη cg
	double *b;	//Το δεξί μέλος του συστήματος
	int n = atoi(argv[1]);	//Το μέγεθος του συστήματος

	//Memory allocation
	Adata = calloc(n, sizeof(double*));
 		if( Adata == NULL)  	    printf("Adata--Out of memory");  
	  	  
		Adata[0] = calloc(n * n, sizeof(double));
		for(i = 1; i < n; i++)
			Adata[i] = Adata[0] + i * n;

	x = calloc(n, sizeof(double));
	cg_x = calloc(n, sizeof(double));
	b = calloc(n, sizeof(double));

	//Παραγωγή τυχαίου συστήματος
	generate_randsys(Adata, x, b, n);

	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			printf("Adata[%d][%d] = %lf\t", i, j, Adata[i][j]);  
		}
		printf("\n");
	}
}*/

