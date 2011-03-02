#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "cg.c"
#include "my_routines.c"
#include <sys/time.h>
#include <time.h>


int main(int argc, char **argv)
{
	struct timeval tvs, tve, tvs1, tve1;
  	double starttime, endtime, sysstart, sysend;
	int i, j;	
	double *Adata;	//Ο πίνακας του συστήματος
	double *Mdata;	//O preconditioner
	double *x;	//Η λύση του συστήματος (που παράγεται για το τυχαίο σύστημα)
	double *cg_x;	//Η λύση του συστήματος που προκύπτει από τη cg
	double *b;	//Το δεξί μέλος του συστήματος
	int n = atoi(argv[1]);	//Το μέγεθος του συστήματος
	int maxiter = 1000;	//max iterations before returning
        double rtol     = 1e-8;
	double abs_error = 0.0;
	FILE *f;

	//Memory allocation
	Adata = malloc(n*n*sizeof(double));
 	if( Adata == NULL)  	    printf("Adata--Out of memory");  

	Mdata = calloc(n, sizeof(double));
	x = calloc(n, sizeof(double));
	cg_x = calloc(n, sizeof(double));
	b = calloc(n, sizeof(double));

	
	gettimeofday(&tvs1, NULL);
  	sysstart = tvs1.tv_sec*1000.0 + tvs1.tv_usec/1000.0;

	//Άνοιγμα αρχείου
	if( (f=fopen(argv[2], "r"))== NULL ){
		perror("fopen");
		exit(EXIT_FAILURE);
	}

	//Παραγωγή τυχαίου συστήματος
	generate_randsys_blas(Adata, x, b, n, f);

	fclose(f);
	
	gettimeofday(&tve1, NULL);
  	sysend = tve1.tv_sec*1000.0 + tve1.tv_usec/1000.0;


	//Εύρεση preconditioner
	jacobi_precond_blas(Mdata, Adata, n);

 	gettimeofday(&tvs, NULL);
  	starttime = tvs.tv_sec*1000.0 + tvs.tv_usec/1000.0;

	//Επίλυση με cg
	int iter = precond_cg_blas(matvec_blas, psolve, Adata, Mdata, b, cg_x, rtol, n, maxiter);
               
	gettimeofday(&tve, NULL);
  	endtime = tve.tv_sec*1000.0 + tve.tv_usec/1000.0;

   	printf("Xronos ekteleshs = %lf msec\n", endtime - starttime); 
	printf("Επαναλήψεις = %d\n", iter);

	//free memory
	free(Adata); free(Mdata); free(x); free(cg_x); free(b); 
}
