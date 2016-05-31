/*
 * Author: Spiros Delviniotis
 * File: main.c
 * Project: MPI 
 */


#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NXPROB      1000                 /* x dimension of problem grid */
#define NYPROB      1000                 /* y dimension of problem grid */
#define STEPS       100                /* number of time steps */
#define MAXWORKER   8                  /* maximum number of worker tasks */
#define MINWORKER   3                  /* minimum number of worker tasks */
#define BEGIN       1                  /* message tag */
#define DTAG        0                  /* message tag */
#define UTAG        1                  /* message tag */
#define LTAG        2                  /* message tag */
#define RTAG        3                  /* message tag */
#define NONE        -1                 /* indicates no neighbor */
#define DONE        4                  /* message tag */
#define MASTER      0                  /* taskid of first process */
#define UP          0
#define DOWN        1
#define LEFT        2
#define RIGHT       3


struct Parms {
	float cx;
	float cy;
} parms = {0.1, 0.1};


int main(int argc, char *argv[])
{
	void	inidat(),
		prtdat(),
		update_outside_table(),
		update_inside_table();
	float *table_u; /* array for grid */
	int	taskid, /* this task's unique id */
		rc,
		numworkers,		/* number of worker processes */
		numtasks,		/* number of tasks */
		offset,			/* for sending rows of data */
		ix, iy, iz, it,		/* loop variables */
		sub_table_dimention,	/* Inner sub table dimention for task's grids */
		sub_x, sub_y,		/* Sub table dimention for task's grids with extra 2 levels perimetrical */
		nbrs[4],
		size,			/* Size of idim specifying the number of processes in each dimension */
		dims[2],		/* Array of size ndims */
		reorder = 1,		/* Ranking may be reordered (true) or not (false) (logical) */
		periods[2] = {0, 0};	/* Logical  array  of  size ndims specifying whether the grid is periodic */
	double 	start_time = 0,		/* start time */
		end_time = 0,		/* end time */
		process_clock = 0,	/* process's duration */
		master_clock  = 0;	/* master's duration */

	/* Logs for Debug */
	FILE * fp[9];
	int temp;

	for (temp = 0; temp < 9; temp++) {
		char str[10];
		sprintf(str, "%d", temp);
		fp[temp] = fopen(str, "w");
	}

	/* First, find out my taskid and how many tasks are running */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	/* Initialization */
	numworkers = numtasks;

	if (NXPROB * NYPROB % numworkers != 0) {
		printf("Main ERROR: Number of tasks must be: %dx%d %% %d != 0 .\n", NXPROB, NYPROB, numworkers);
		return(EXIT_FAILURE);
	}

	/* Create Cartesian Topology */
	MPI_Comm cartcomm;
	size = sqrt(numworkers);
	dims[0] = size;
	dims[1] = size;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm);

	/* Blocks  the  caller  until  all  group members have called it */
	MPI_Barrier(cartcomm);

	/* Returns  the  shifted source and destination ranks, given a  shift direction and amount */
	MPI_Cart_shift(cartcomm, 0, 1, &nbrs[UP], &nbrs[DOWN]);
	MPI_Cart_shift(cartcomm, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]);

	sub_table_dimention = sqrt(NXPROB * NYPROB / numtasks);
	sub_x = sub_table_dimention + 2;
	sub_y = sub_table_dimention + 2;

	table_u = (float*) malloc((2 * sub_x * sub_y) * sizeof(float)); //free must be added!
	if (table_u == NULL) {
		printf("Main ERROR: Allocation memory.\n");
		MPI_Abort(cartcomm, rc);
		return(EXIT_FAILURE);
	}

	for (iz = 0; iz < 2; iz++) {
		for (ix = 0; ix < sub_x; ix++) {
			for (iy = 0; iy < sub_y; iy++) {
				offset = iz * sub_x * sub_y + ix * sub_y + iy;
				*(table_u + offset) = 0.0;
				//printf("	[%d]\n",  offset);
			}
		}
	}

	/* Initialize Table */
	inidat(sub_table_dimention, sub_table_dimention, sub_y, (table_u + sub_y + 1));

	iz = 0;
	MPI_Request req[8];
	MPI_Status status[8];

	/* Datatype Definition */
	MPI_Datatype COL_INT;
	MPI_Type_vector(sub_table_dimention, 1, sub_y, MPI_FLOAT, &COL_INT);
	MPI_Type_commit(&COL_INT);
	
	start_time = MPI_Wtime();  /* Start Timer */
	
	for (it = 1; it <= STEPS; it++){
//		printf(" --> In Step[%d]: nbrs[UP] = %d.\n", it, nbrs[UP]);
//		printf(" --> In Step[%d]: nbrs[DOWN] = %d.\n", it, nbrs[DOWN]);
//		printf(" --> In Step[%d]: nbrs[LEFT] = %d.\n", it, nbrs[LEFT]);
//		printf(" --> In Step[%d]: nbrs[RIGHT] = %d.\n", it, nbrs[RIGHT]);
		if (nbrs[UP] >= 0){								//up
			MPI_Isend(table_u + iz*sub_x*sub_y + sub_y + 1, sub_table_dimention, MPI_FLOAT, nbrs[UP], DTAG, cartcomm, &req[0]);
			MPI_Irecv(table_u + iz*sub_x*sub_y + 1, sub_table_dimention, MPI_FLOAT, nbrs[UP], UTAG, cartcomm, &req[1]);
		}

		if (nbrs[DOWN] >= 0){							//down
			MPI_Isend(table_u + iz*sub_x*sub_y + sub_table_dimention*sub_y + 1, sub_table_dimention , MPI_FLOAT, nbrs[DOWN], UTAG, cartcomm, &req[2]);
			MPI_Irecv(table_u + iz*sub_x*sub_y + (sub_table_dimention+1)*sub_y + 1, sub_table_dimention , MPI_FLOAT, nbrs[DOWN], DTAG, cartcomm, &req[3]);
		}

		if (nbrs[LEFT] >= 0){								//left
			MPI_Isend(table_u + iz*sub_x*sub_y + sub_y + 1, 1, COL_INT, nbrs[LEFT], RTAG, cartcomm,&req[4]);
			MPI_Irecv(table_u + iz*sub_x*sub_y + sub_y, 1, COL_INT, nbrs[LEFT], LTAG, cartcomm, &req[5]);
		}

		if (nbrs[RIGHT] >= 0 ){							//right
			MPI_Isend(table_u + iz*sub_x*sub_y + sub_y + sub_table_dimention, 1, COL_INT, nbrs[RIGHT], LTAG, cartcomm,&req[6]);
			MPI_Irecv(table_u + iz*sub_x*sub_y + sub_y + sub_table_dimention + 1, 1, COL_INT , nbrs[RIGHT], RTAG, cartcomm,&req[7]);
		}

		update_inside_table(sub_table_dimention - 2, table_u + iz*sub_x*sub_y, table_u + (1-iz)*sub_x*sub_y);

		if(nbrs[UP] >= 0){
			MPI_Wait(&req[0],&status[0]);
			MPI_Wait(&req[1],&status[1]);
		}
		if(nbrs[DOWN] >= 0){
			MPI_Wait(&req[2],&status[2]);
			MPI_Wait(&req[3],&status[3]);
		}
		if(nbrs[LEFT] >= 0){
			MPI_Wait(&req[4],&status[4]);
			MPI_Wait(&req[5],&status[5]);
		}
		if(nbrs[RIGHT] >= 0){
			MPI_Wait(&req[6],&status[6]);
			MPI_Wait(&req[7],&status[7]);
		}

		update_outside_table(sub_table_dimention, table_u + iz*sub_x*sub_y, table_u + (1-iz)*sub_x*sub_y);
		iz = 1 - iz;
//		printf(" -> End of i loop.\n");
	}
	
	end_time = MPI_Wtime();  /* Stop Timer */
	
	prtdat(fp[taskid], sub_x, sub_y, table_u, 1, 1, sub_table_dimention, taskid);

	process_clock = end_time - start_time;
	MPI_Reduce(&process_clock, &master_clock, 1, MPI_DOUBLE, MPI_MAX, 0, cartcomm);
	MPI_Barrier(cartcomm);

	if (taskid == MASTER){
		printf("Total time elapsed = %lf \n", master_clock);
	}

	/* Free resources */
	free(table_u);
	MPI_Type_free(&COL_INT);
	MPI_Finalize();
	
}


/**************************************************************************
 *  subroutines update_calculation, update_inside_table, update_outside_table
 ****************************************************************************/
void update_calculation(int ix, int iy, int y, float *u1, float *u2)
{

	*(u2 + ix * y + iy) = *(u1 + ix * y + iy) +
		parms.cx * (*(u1 + (ix + 1) * y + iy) +
		*(u1 + (ix - 1) * y + iy) -
		2.0 * *(u1 + ix * y + iy)) +
		parms.cy * (*(u1 + ix * y + iy + 1) +
		*(u1 + ix * y + iy - 1) -
		2.0 * *(u1 + ix * y + iy));

}


void update_inside_table(int end, float *u1, float *u2)
{
	int i, j;
	for (i = 2; i <= end + 1; i++) { //end - 1 ???
		for (j = 2; j <= end + 1; j++) { //end - 1 ???
			update_calculation(i, j, end + 4, u1, u2); //end+4 ????
		}
	}
}


void update_outside_table(int end, float *u1, float *u2)
{
	int i;
	for (i = 1; i <= end; i++) {
		update_calculation(1, i, end + 2, u1, u2);
		update_calculation(end, i, end + 2, u1, u2);
		update_calculation(i, 1, end + 2, u1, u2);
		update_calculation(i, end, end + 2, u1, u2);
	}
}


/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, int y, float *u)
{
	int ix, iy;

	for (ix = 0; ix <= nx - 1; ix++) {
		for (iy = 0; iy <= ny - 1; iy++) {
			*(u + ix * (y) + iy) = (float) (ix * (nx - ix - 1) * iy * (ny - iy - 1) + 10);
			// here must be copied data from original initialized table.
			//printf("	[%d] = %6.1f", (int) (u + ix * (y) + iy), *(u + ix * (y) + iy));
		}
	}
}


/**************************************************************************
 * subroutines prtdat
 **************************************************************************/
void prtdat(FILE *fp, int x, int y, float *u, int start1, int start2, int end, int taskid)
{
	int j, i;
	for (i = 0; i < x; i++) {
		for (j = 0; j < y; j++) {
			fprintf(fp, "%6.1f ", *(u + i * y + j));
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n\n");
	for (i = start1; i <= start1 + end - 1; i++) {
		for (j = start2; j <= start2 + end - 1; j++) {
			fprintf(fp, "%6.1f ", *(u + i * y + j));
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
	printf("taskid in fprint = %d\n", taskid);
}

