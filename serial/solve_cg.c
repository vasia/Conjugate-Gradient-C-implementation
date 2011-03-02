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
	double **Adata;	//Ο πίνακας του συστήματος
	double *Mdata;	//O preconditioner
	double *x;	//Η λύση του συστήματος (που παράγεται για το τυχαίο σύστημα)
	double *cg_x;	//Η λύση του συστήματος που προκύπτει από τη cg
	double *b;	//Το δεξί μέλος του συστήματος
	int n = atoi(argv[1]);	//Το μέγεθος του συστήματος
	int maxiter = 1000;	//max iterations before returning
        double rtol     = 1e-8;
	double abs_error = 0.0;

	//Memory allocation
	Adata = calloc(n, sizeof(double*));
 		if( Adata == NULL)  	    printf("Adata--Out of memory");  
	  	  
		Adata[0] = calloc(n * n, sizeof(double));
		for(i = 1; i < n; i++)
			Adata[i] = Adata[0] + i * n;

	Mdata = calloc(n, sizeof(double));
	x = calloc(n, sizeof(double));
	cg_x = calloc(n, sizeof(double));
	b = calloc(n, sizeof(double));

	//Παραγωγή τυχαίου συστήματος
	gettimeofday(&tvs1, NULL);
  	sysstart = tvs1.tv_sec*1000.0 + tvs1.tv_usec/1000.0;

	generate_randsys(Adata, x, b, n);
	
	gettimeofday(&tve1, NULL);
  	sysend = tve1.tv_sec*1000.0 + tve1.tv_usec/1000.0;

	/*for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			printf("Adata[%d][%d] = %lf\t", i, j, Adata[i][j]);  
		}
		printf("\n");
	}*/

	//Εύρεση preconditioner
	jacobi_precond(Mdata, Adata, n);

 	gettimeofday(&tvs, NULL);
  	starttime = tvs.tv_sec*1000.0 + tvs.tv_usec/1000.0;

	//Επίλυση με cg
	int iter = precond_cg(matvec, psolve, Adata, Mdata,
  		b, cg_x, rtol, n, maxiter);
               
	gettimeofday(&tve, NULL);
  	endtime = tve.tv_sec*1000.0 + tve.tv_usec/1000.0;

   	printf("Xronos ekteleshs = %lf msec\n", endtime - starttime); 

	printf("Xronos dhmiourgias tyxaioy systhmatos = %lf msec\n", sysend - sysstart); 

	//for(i=0; i<n; i++)
		//printf("x(real) = %lf - xout = %lf, dif = %lf\n", x[i], cg_x[i], x[i]-cg_x[i]);

	//Εύρεση απόλυτου error
	for(i=0; i<n; i++){
		if(fabs(x[i]-cg_x[i]) > abs_error)
			abs_error = fabs(x[i]-cg_x[i]);
	}

	printf("Επαναλήψεις = %d, Abs_error = %lf\n", iter, abs_error);

	//free memory
	free(Adata[0]); free(Adata);
    	free(Mdata); free(x); free(cg_x); free(b); 
}
