#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "cg_MPI.c"
#include "my_routines_MPI.c"
//#include "timers.h"
//#include <sys/time.h>
//#include <time.h>
#include <mpi.h>

int main(int argc, char **argv)
{
	struct timeval tvs, tve, tvs1, tve1;
	int i, j;	
	double *Adata;	//Ο πίνακας του συστήματος
	double *AdataLocal;	//Ο πίνακας του συστήματος
	double *Mdata;	//O preconditioner
	double *x;	//Η λύση του συστήματος (που παράγεται για το τυχαίο σύστημα)
	double *cg_x;	//Η λύση του συστήματος που προκύπτει από τη cg
	double *b;	//Το δεξί μέλος του συστήματος
	double *MdataLocal;	//O preconditioner
	double *bLocal;	//Το δεξί μέλος του συστήματος
	int n = atoi(argv[1]);	//Το μέγεθος του συστήματος
	int nLocal, MnLocal;	//local μεγεθος συστήματος και preconditioner
	int maxiter = 1000;	//max iterations before returning
        double rtol     = 1e-8;
	double abs_error = 0.0;
	FILE *f;

	//MPI-related variables
	int size,rank, source, dest;
  	MPI_Status status;
	int *sendcnts; 
        int *displs;
	int *sendcntsb; 
        int *displsb;
	int myoffset;

	//timers
	double exec_timer = 0.0;
	double comm_timer = 0.0;
	
	//Arxikopoihsh MPI
 	MPI_Init(&argc,&argv);

	//Swzoume sto size ton ari8mo twn diergasiwn
  	MPI_Comm_size(MPI_COMM_WORLD,&size);

	//Ka8e diergasia pairnei to rank ths
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	//Υπολογισμός local size συστήματος
	nLocal = n*n/size;
	MnLocal = n/size;
	if(n*n%size>0){
		if(rank<n*n%size)
			nLocal++;
	}
	if(n%size>0){
		if(rank<n%size)
			MnLocal++;
	}
	//printf("rank = %d, nLocal = %d\n", rank, nLocal);
	//Memory allocation
	//AdataLocal = malloc(nLocal*sizeof(double));
	AdataLocal = malloc(MnLocal*n*sizeof(double));
 	if( AdataLocal == NULL)  	    printf("AdataLocal--Out of memory");  

	MdataLocal = calloc(MnLocal, sizeof(double));
	cg_x = calloc(MnLocal, sizeof(double));
	bLocal = calloc(MnLocal, sizeof(double));


	//Άνοιγμα αρχείου
	if(rank==0){
		Adata = malloc(n*n*sizeof(double));
 		if( Adata == NULL)  	    printf("Adata--Out of memory");  
		b = calloc(n, sizeof(double));

		if( (f=fopen(argv[2], "r"))== NULL ){
			perror("fopen");
			exit(EXIT_FAILURE);
		}

		//Παραγωγή τυχαίου συστήματος
		x = calloc(n, sizeof(double));
		generate_randsys_blas(Adata, x, b, n, f);

		fclose(f);

		//Προετοιμασία για την κλήση της MPI_Scatterv()
		//1. Scatter matrix A
		sendcnts = malloc(size*sizeof(int));
 		if( sendcnts == NULL)  	    printf("sendcnts--Out of memory");  
		displs = malloc(size*sizeof(int));
 		if( displs == NULL)  	    printf("displs--Out of memory");  
		
		int tmpn = n/size;
		int offset = 0;
		for(i=0; i<size; i++){
			if(n%size>0){
				if(i<n%size)
					sendcnts[i] = (tmpn+1)*n;
				else	
					sendcnts[i] = tmpn*n;
			}
			else
				sendcnts[i] = tmpn*n;
		displs[i] = offset;
		offset += sendcnts[i];
		}

		//2. Send row - offsets
		myoffset = 0;
		int stat_send;
		for(i=1; i<size; i++)
			stat_send = MPI_Send(&displs[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		
	}

	else{	//rank != 0
		//Receive row-offsets
		int stat_recv = MPI_Recv(&myoffset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		myoffset = myoffset/n;
	}
	
	//3. Scatter vector b
	// *** Οι ίδιοι πίνακες sendcntsb, displsb θα χρησιμοποιηθούν και στη κλήση της MPI_Allgatherv 
	// μέσα στη ρουτίνα precond_cg_blas *** //
	sendcntsb = malloc(size*sizeof(int));
	if( sendcntsb == NULL)  	    printf("sendcntsb--Out of memory");  
	displsb = malloc(size*sizeof(int));
	if( displsb == NULL)  	    printf("displsb--Out of memory");  

	int tmpnb = n/size;
	int offsetb = 0;
	for(i=0; i<size; i++){
		if(n%size>0){
			if(i<n%size)
				sendcntsb[i] = tmpnb+1;
			else	
				sendcntsb[i] = tmpnb;
		}
		else
			sendcntsb[i] = tmpnb;
	displsb[i] = offsetb;
	offsetb += sendcntsb[i];
	}

	//Μοίρασμα πίνακα συστήματος
	int stat = MPI_Scatterv(Adata, sendcnts, displs, MPI_DOUBLE, AdataLocal, n*MnLocal, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//Μοίρασμα διανύσματος b
	int statb = MPI_Scatterv(b, sendcntsb, displsb, MPI_DOUBLE, bLocal, MnLocal, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//Εύρεση preconditioner
	jacobi_precond_blas(MdataLocal, AdataLocal, MnLocal, n, myoffset);

	MPI_Barrier(MPI_COMM_WORLD);
	//Επίλυση με cg
	start_timer(exec_timer);
	int iter = precond_cg_blas(matvec_blas, psolve, AdataLocal, MdataLocal, bLocal, cg_x, rtol, MnLocal, n, maxiter, sendcntsb, displsb, &comm_timer);
	MPI_Barrier(MPI_COMM_WORLD);
	stop_timer(exec_timer);

   //	printf("Xronos ekteleshs = %lf msec\n", endtime - starttime); 

	//printf("Xronos dhmiourgias tyxaioy systhmatos = %lf msec\n", sysend - sysstart); 

	//for(i=0; i<n; i++)
		//printf("x(real) = %lf - xout = %lf, dif = %lf\n", x[i], cg_x[i], x[i]-cg_x[i]);

	if(rank ==0 ) {printf("Επαναλήψεις = %d\n", iter);
		 printf("Run time:%f, Communication time:%f\n",timer_duration(exec_timer), timer_duration(comm_timer));
	}
	//free memory
	free(AdataLocal); free(MdataLocal); free(cg_x); free(bLocal); free(sendcntsb); free(displsb);
	 MPI_Finalize();
}
