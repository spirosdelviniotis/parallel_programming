/*
 * Author: Spiros Delviniotis
 * File: main.c
 * Project: MPI 
 */

/* TO DO */
// Add timer for calculation of execution time
// Make code robust
// Remove prints and unusable comments
// Handle return values from MPI functions


#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#else
int omp_get_thread_num(void) { return 0; }
int omp_get_num_threads(void) { return 1; }
#endif

#define NXPROB      1000                 /* x dimension of problem grid */
#define NYPROB      1000                 /* y dimension of problem grid */
#define STEPS       100                /* number of time steps */


struct Parms {
	float cx;
	float cy;
} parms = {0.1, 0.1};


int main(int argc, char *argv[])
{
	void	inidat(),
		update();
	float *table_u; /* array for grid */
	int	offset,			/* for sending rows of data */
		ix, iy, iz, it;		/* loop variables */
	double 	start_time = 0,		/* start time */
		end_time = 0,		/* end time */
		process_clock = 0,	/* process's duration */
		master_clock  = 0;	/* master's duration */



	table_u = (float*) malloc((2 * NXPROB * NYPROB) * sizeof(float)); //free must be added!
	if (table_u == NULL) {
		printf("Main ERROR: Allocation memory.\n");
		return(EXIT_FAILURE);
	}

	for (iz = 0; iz < 2; iz++) {
		for (ix = 0; ix < NXPROB; ix++) {
			for (iy = 0; iy < NYPROB; iy++) {
				offset = iz * NXPROB * NYPROB + ix * NYPROB + iy;
				*(table_u + offset) = 0.0;
			}
		}
	}

	/* Initialize Table */
	inidat(NXPROB, NYPROB, table_u);
	iz = 0;
	
	for (it = 1; it <= STEPS; it++){
		update(1, NXPROB - 2 + iz*NXPROB*NYPROB, table_u + (1-iz)*NXPROB*NYPROB);  //?????
		
		iz = 1 - iz;
	}
	
	/* Free resources */
	free(table_u);
	
}


/*****************************************************************************
 *  subroutine update
 *****************************************************************************/
void update(int start, int end, int ny, float *u1, float *u2)
{
	int ix, iy;
	
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		int total = omp_get_num_threads();
		printf("	This is thread %d of %d\n", tid, total);
		#pragma omp for
			for (ix = start; ix <= end; ix++) {
				for (iy = 1; iy <= ny-2; iy++) { 
					*(u2+ix*ny+iy) = *(u1+ix*ny+iy)  + 
							parms.cx * (*(u1+(ix+1)*ny+iy) +
							*(u1+(ix-1)*ny+iy) - 
							2.0 * *(u1+ix*ny+iy)) +
							parms.cy * (*(u1+ix*ny+iy+1) +
							*(u1+ix*ny+iy-1) - 
							2.0 * *(u1+ix*ny+iy));
				}
			}
	}
}
	

/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, float *u) {
	int ix, iy;

	for (ix = 0; ix <= nx-1; ix++) {
		for (iy = 0; iy <= ny-1; iy++) {
			*(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
		}
	}
}
