//#define _REENTRANT

#include <time.h>               /* clock */

/* A simple wrapper around the clock() timer */
typedef struct _Timer {
    clock_t clock_holder;
    clock_t duration_clocks;
} Timer;


extern double	comp_timer,init_timer,mat_mpi_timer,ddot_mpi_timer,mpi_timer,runtime_timer,matvec_timer,ddot_timer,ax_timer,psolve_timer;

#define	start_timer(timer)	((timer) -= MPI_Wtime())
#define	stop_timer(timer)	((timer) += MPI_Wtime())

#define	timer_duration(timer)	(timer)
