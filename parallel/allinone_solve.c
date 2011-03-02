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

	//solve parameters
	double bnorm2;              /* ||b||^2 */
	double rnorm2;              /* Νόρμα υπολοίπου στο τετράγωνο */
	double rz, rzold;           /* r'*z για 2 διαδοχικές επαναλήψεις */
	double alpha, beta;
	double *s;                  /* Κατεύθυνση αναζήτησης */
	double *sGlobal;                  /* Κατεύθυνση αναζήτησης */
	double *r;                  /* Υπόλοιπο         */
	double *z;                  /* Προσωρινό διάνυσμα */

	//timers
	double allgather_timer = 0.0;
	double ddot_timer1 = 0.0;
	double ddot_timer2 = 0.0;
	double ddot_timer3 = 0.0;
	int stat;



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
	stat = MPI_Scatterv(Adata, sendcnts, displs, MPI_DOUBLE, AdataLocal, n*MnLocal, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//Μοίρασμα διανύσματος b
	int statb = MPI_Scatterv(b, sendcntsb, displsb, MPI_DOUBLE, bLocal, MnLocal, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//Εύρεση preconditioner
	jacobi_precond_blas(MdataLocal, AdataLocal, MnLocal, n, myoffset);

	MPI_Barrier(MPI_COMM_WORLD);
	//Επίλυση με cg
	start_timer(exec_timer);
//	int iter = precond_cg_blas(matvec_blas, psolve, AdataLocal, MdataLocal, bLocal, cg_x, rtol, MnLocal, n, maxiter, sendcntsb, displsb, &comm_timer);

    int iter;
    const int nbytes = MnLocal * sizeof(double);
    sGlobal = malloc(n*sizeof(double));
    s =  malloc(nbytes);
    r =  malloc(nbytes);
    z =  malloc(nbytes);

    bnorm2    = ddot(MnLocal, bLocal, 1, bLocal, 1, &comm_timer);
    
    memset(cg_x, 0, nbytes);	//αρχικοποίηση λύσης
    memcpy(r, bLocal, nbytes);	//και υπολοίπου - r0=b-A*x0 (x0=0)

    psolve(z, MdataLocal, r, MnLocal);	//εφαρμογή του preconditioner - z0 = (M στην -1)*r0

    memcpy(s, z, nbytes);	//αρχικοποίηση κατεύθυνσης αναζήτησης	- p0 = z0

    /* Αρχικοποίηση rz και rnorm2 */
    rz        = ddot(MnLocal, r, 1, z, 1, &comm_timer);
    rnorm2    = ddot(MnLocal, r, 1, r, 1, &comm_timer);

   for (i = 0; i < maxiter ; ++i) {

//	start_timer(*comm_timer);
	start_timer(allgather_timer);
	stat = MPI_Allgatherv(s, MnLocal, MPI_DOUBLE, sGlobal, sendcntsb, displsb, MPI_DOUBLE, MPI_COMM_WORLD);
//	stop_timer(*comm_timer);
	stop_timer(allgather_timer);
	matvec_blas(MnLocal, n, AdataLocal, sGlobal, z);
        
        // Ddot
        alpha = rz / ddot(MnLocal, s, 1, z, 1, &ddot_timer1);	//ak = rkT*zk/pkT*A*pk
        axpy(MnLocal, alpha, s, 1, cg_x, 1);	//xk+1 = xk + ak*pk
        axpy(MnLocal, -alpha, z, 1, r, 1);	//rk+1 = rk - ak*A*pk
  
        psolve(z, MdataLocal, r, MnLocal);		//zk+1 = (M στην -1)*rk+1

        rzold = rz;
        
	rz = ddot(MnLocal, r, 1, z, 1, &ddot_timer2);  	//rTk+1*zk+1
        beta = -rz / rzold;		//β = rTk+1*zk+1/rTk*zk 
	axpy(MnLocal, -beta-1, s, 1, s, 1);	//pk+1 = zk+1+βk*pk
	axpy(MnLocal, 1, z, 1, s, 1);	//pk+1 = zk+1+βk*pk
        
        
        // check error
	rnorm2     = ddot(MnLocal, r, 1, r, 1, &ddot_timer3);
        if(rnorm2 <= bnorm2 * rtol * rtol)
        	break;
        
    
    }
    
    free(z);
    free(r);
    free(s);

    if (i >= maxiter)
        iter = -1;
    else
    	iter = i;


	MPI_Barrier(MPI_COMM_WORLD);
	stop_timer(exec_timer);

if(rank ==0 ){
/*	 printf("Allgather time:%f, ddot time1:%f, ddot time2:%f, ddot time3:%f\n",timer_duration(allgather_timer), timer_duration(ddot_timer1), timer_duration(ddot_timer2), timer_duration(ddot_timer3));
	 printf("Run time:%f\n",timer_duration(exec_timer));
	 printf("Επαναλήψεις = %d\n", iter);*/
       	 printf("%f&%f&%f&%f&%f\n",timer_duration(allgather_timer), timer_duration(ddot_timer1), timer_duration(ddot_timer2), timer_duration(ddot_timer3), timer_duration(exec_timer));
}
	//free memory
	free(AdataLocal); free(MdataLocal); free(cg_x); free(bLocal); free(sendcntsb); free(displsb);
	 MPI_Finalize();
}
