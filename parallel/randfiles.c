#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "my_routines.c"
#include <string.h>

/*
*	Δημιουργία αρχείων με τυχαία συστήματα
*
*/
int main(int argc, char **argv)
{
		
	double **Adata;	//Ο πίνακας του συστήματος
	double *x;	//Η λύση του συστήματος (που παράγεται για το τυχαίο σύστημα)
	double *b;	//Το δεξί μέλος του συστήματος
	int n = atoi(argv[1]);	//Το μέγεθος του συστήματος
	int i,j;
	FILE *f1;
	

	//Memory allocation
	Adata = calloc(n, sizeof(double*));
 		if( Adata == NULL)  	    printf("Adata--Out of memory");  
	  	  
		Adata[0] = calloc(n * n, sizeof(double));
		for(i = 1; i < n; i++)
			Adata[i] = Adata[0] + i * n;

	
	x = calloc(n, sizeof(double));
	b = calloc(n, sizeof(double));

	//Παραγωγή τυχαίου συστήματος
	generate_randsys(Adata, x, b, n);
	
	//εγγραφη σε αρχειο
	if( (f1=fopen("sys8192", "w"))== NULL ){
		perror("fopen");
		exit(EXIT_FAILURE);
	}
	
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			fprintf(f1, "%lf", Adata[i][j]);
			fprintf(f1, " ");
		}
		fprintf(f1, "\n");
	}

	fclose(f1);


	
	//free memory
	free(Adata[0]); free(Adata);
    	free(x); free(b);
}
