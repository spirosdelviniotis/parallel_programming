/*
 * Author: Spiros Delviniotis
 * File: main.c
 * Project: MPI 
*/


#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NXPROB      6                 /* x dimension of problem grid */
#define NYPROB      6                 /* y dimension of problem grid */
#define STEPS       100                /* number of time steps */
#define MAXWORKER   8                  /* maximum number of worker tasks */
#define MINWORKER   3                  /* minimum number of worker tasks */
#define BEGIN       1                  /* message tag */
#define LTAG        2                  /* message tag */
#define RTAG        3                  /* message tag */
#define NONE        0                  /* indicates no neighbor */
#define DONE        4                  /* message tag */
#define MASTER      0                  /* taskid of first process */


struct Parms { 
	float cx;
	float cy;
} parms = {0.1, 0.1};


int main (int argc, char *argv[])
{
	void	inidat(),
		prtdat(), 
		update();
	float	u[2][NXPROB][NYPROB];       /* array for grid */
	float	*table_u;		    /* array for grid */
	int	taskid,                     /* this task's unique id */
		numworkers,                 /* number of worker processes */
		numtasks,                   /* number of tasks */
		averow,rows,offset,extra,   /* for sending rows of data */
		dest, source,               /* to - from for message send-receive */
		left,right,		    /* neighbor tasks */
		msgtype,                    /* for message types */
		rc,start,end,               /* misc */
		i,ix,iy,iz,it,              /* loop variables */
		sub_table_dimention,	    /* Inner sub table dimention for task's grids */
		sub_x, sub_y;		    /* Sub table dimention for task's grids with extra 2 levels perimetrical */
	MPI_Status status;
	
	/* Logs for Debug */
	FILE *fp[9];
	int temp;

	for(temp = 0; temp<9; temp++){
		char str[10];
		sprintf(str, "%d", temp);
		fp[temp] = fopen(str, "w");
	}
	
	
	/* First, find out my taskid and how many tasks are running */
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
	
	/* Initialization */
	numworkers = numtasks;
	
	if ( NXPROB*NYPROB% numworkers != 0 ){
		printf("Main ERROR: Number of tasks must be: %dx%d %% %d != 0 .\n",NXPROB,NYPROB,numworkers);
		return (EXIT_FAILURE);
	}
	
	sub_table_dimention = sqrt(NXPROB*NYPROB/numtasks);
	sub_x = sub_table_dimention + 2;
	sub_y = sub_table_dimention + 2;
	
	printf("---> sub_dimention = %d\n", sub_table_dimention);
	printf("---> numtasks = %d\n", numtasks);
	
	table_u = (float*)malloc((2*sub_x*sub_y)*sizeof(float));  //free must be added!
	if (table_u == NULL){
		printf("Main ERROR: Allocation memory.\n");
		return (EXIT_FAILURE);
	}
	
	for (iz=0; iz < 2; iz++){
		for (ix=0; ix < sub_x; ix++){
			for (iy=0; iy < sub_y; iy++){
				offset = iz*sub_x*sub_y+ix*sub_y+iy;
				*(table_u + offset) = 0.0;
				//printf("	[%d]\n",  offset);
			}
		}
	}
	
		
	
	/**********************************************************************/
	numworkers = numtasks - 1;

	if (taskid == MASTER) {
		/************************* master code *******************************/
		/* Check if numworkers is within range - quit if not */
		if ((numworkers > MAXWORKER) || (numworkers < MINWORKER)) {
			printf("ERROR: the number of tasks must be between %d and %d.\n",MINWORKER+1,MAXWORKER+1);
			printf("Quitting...\n");
			MPI_Abort(MPI_COMM_WORLD, rc);
			exit(1);
		}
		printf ("Starting mpi_heat2D with %d worker tasks.\n", numworkers);

		/* Initialize grid */
		printf("Grid size: X= %d  Y= %d  Time steps= %d\n",NXPROB,NYPROB,STEPS);
		printf("Initializing grid and writing initial.dat file...\n");
		inidat(NXPROB, NYPROB, u);
		//new_inidat(sub_table_dimention, sub_table_dimention, sub_y, (table_u + sub_y + 1));
		prtdat(NXPROB, NYPROB, u, "initial.dat");

		/* Distribute work to workers.  Must first figure out how many rows to */
		/* send and what to do with extra rows.  */
		averow = NXPROB/numworkers;
		extra = NXPROB%numworkers;
		printf("--> averow: %d, extra: %d\n",averow,extra);
		offset = 0;
		for (i=1; i<=numworkers; i++)
		{
			rows = (i <= extra) ? averow+1 : averow; 
			/* Tell each worker who its neighbors are, since they must exchange */
			/* data with each other. */  
			if (i == 1) 
			   left = NONE;
			else
			   left = i - 1;
			if (i == numworkers)
			   right = NONE;
			else
			   right = i + 1;

			/*  Now send startup information to each worker  */
			dest = i;
			MPI_Send(&offset, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&left, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&right, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&u[0][offset][0], rows*NYPROB, MPI_FLOAT, dest, BEGIN, MPI_COMM_WORLD);

			printf("Sent to task %d: rows = %d, offset = %d ",dest,rows,offset);
			printf("left = %d, right = %d\n",left,right);
			offset = offset + rows;
		}
		
		/* Now wait for results from all worker tasks */
		for (i=1; i<=numworkers; i++)
		{
			source = i;
			msgtype = DONE;
			MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&u[0][offset][0], rows*NYPROB, MPI_FLOAT, source, msgtype, MPI_COMM_WORLD, &status);
		}

		/* Write final output, call X graph and finalize MPI */		
		printf("Writing final.dat file and generating graph...\n");
		prtdat(NXPROB, NYPROB, &u[0][0][0], "final.dat");
		printf("Click on MORE button to view initial/final states.\n");
		printf("Click on EXIT button to quit program.\n");
      
		MPI_Finalize();
	}   /* End of master code */

	/************************* workers code **********************************/
	if (taskid != MASTER) 
	{
		new_inidat(sub_table_dimention, sub_table_dimention, sub_y, (table_u + sub_y + 1));
		
		/* Initialize everything - including the borders - to zero */
		for (iz=0; iz<2; iz++)
			for (ix=0; ix<NXPROB; ix++) 
				for (iy=0; iy<NYPROB; iy++) 
					u[iz][ix][iy] = 0.0;

		/* Receive my offset, rows, neighbors and grid partition from master */
		source = MASTER;
		msgtype = BEGIN;
		MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&left, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&right, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&u[0][offset][0], rows*NYPROB, MPI_FLOAT, source, msgtype, MPI_COMM_WORLD, &status);

		/* Determine border elements.  Need to consider first and last columns. */
		/* Obviously, row 0 can't exchange with row 0-1.  Likewise, the last */
		/* row can't exchange with last+1.  */
		if (offset == 0) 
		   start = 1;
		else 
		   start = offset;
		
		if ((offset + rows) == NXPROB) 
		   end = start + rows - 2;
		else 
		   end = start + rows - 1;

		/* Begin doing STEPS iterations.  Must communicate border rows with */
		/* neighbors.  If I have the first or last grid row, then I only need */
		/*  to  communicate with one neighbor  */
		printf("Task %d received work. Beginning time steps...\n",taskid);
		iz = 0;
		for (it = 1; it <= STEPS; it++)
		{
			if (left != NONE)
			{
				MPI_Send(&u[iz][offset][0], NYPROB, MPI_FLOAT, left, RTAG, MPI_COMM_WORLD);
				source = left;
				msgtype = LTAG;
				MPI_Recv(&u[iz][offset-1][0], NYPROB, MPI_FLOAT, source, msgtype, MPI_COMM_WORLD, &status);
			}
			if (right != NONE)
			{
				MPI_Send(&u[iz][offset+rows-1][0], NYPROB, MPI_FLOAT, right, LTAG, MPI_COMM_WORLD);
				source = right;
				msgtype = RTAG;
				MPI_Recv(&u[iz][offset+rows][0], NYPROB, MPI_FLOAT, source, msgtype, MPI_COMM_WORLD, &status);
			}
			/* Now call update to update the value of grid points */
			update(start,end,NYPROB,&u[iz][0][0],&u[1-iz][0][0]);
			iz = 1 - iz;
		}

		/* Finally, send my portion of final results back to master */
		MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
		MPI_Send(&rows, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
		MPI_Send(&u[iz][offset][0], rows*NYPROB, MPI_FLOAT, MASTER, DONE, MPI_COMM_WORLD);
		fprint(fp[taskid], sub_x, sub_y, table_u, 1, 1, sub_table_dimention, taskid);

		/* Free resources */
		free(table_u);
		MPI_Finalize();
	}
}


/**************************************************************************
 *  subroutine update
 ****************************************************************************/

void update_calculation(int ix, int iy, int y, float *u1, float *u2) {

	*(u2+ix*y+iy) = *(u1+ix*y+iy)  +
			parms.cx * (*(u1+(ix+1)*y+iy) +
			*(u1+(ix-1)*y+iy) -
			2.0 * *(u1+ix*y+iy)) +
			parms.cy * (*(u1+ix*y+iy+1) +
			*(u1+ix*y+iy-1) -
			2.0 * *(u1+ix*y+iy));

}


void update(int start, int end, int ny, float *u1, float *u2)
{
   int ix, iy;
   for (ix = start; ix <= end; ix++) 
      for (iy = 1; iy <= ny-2; iy++) 
         *(u2+ix*ny+iy) = *(u1+ix*ny+iy)  + 
                          parms.cx * (*(u1+(ix+1)*ny+iy) +
                          *(u1+(ix-1)*ny+iy) - 
                          2.0 * *(u1+ix*ny+iy)) +
                          parms.cy * (*(u1+ix*ny+iy+1) +
                         *(u1+ix*ny+iy-1) - 
                          2.0 * *(u1+ix*ny+iy));
}

/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, float *u) {
	int ix, iy;

	for (ix = 0; ix <= nx-1; ix++){
		for (iy = 0; iy <= ny-1; iy++){
			*(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
			printf("	[%d] = %6.1f", (int)(u+ix*ny+iy), *(u+ix*ny+iy));
		}
		printf("\n---- OLD ----\n");
	}
}

void new_inidat(int nx, int ny, int y , float *u) {
	int ix, iy;
	
	for (ix = 0; ix <= nx-1; ix++){
		for (iy = 0; iy <= ny-1; iy++){
			*(u+ix*(y)+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1) + 10); 
			// here must be copied data from original initialized table.
			printf("	[%d] = %6.1f", (int)(u+ix*(y)+iy), *(u+ix*(y)+iy));
		}
		printf("\n---- NEW ----\n");
	}
}


/**************************************************************************
 * subroutine prtdat
 **************************************************************************/

void fprint(FILE *fp,int x,int y,float *u, int start1, int start2, int end, int taskid)
{
	int j,i;
	for(i=0;i<x;i++){
		for(j=0;j<y;j++){
			fprintf(fp, "%6.1f ",*(u+i*y+j));
		}
		fprintf(fp,"\n");
	}

	fprintf(fp, "\n\n");
	for (i = start1; i <= start1+end-1; i++){
		for (j = start2; j <= start2+end-1; j++){
			fprintf(fp, "%6.1f ",*(u+i*y+j));
		}
		fprintf(fp,"\n");
	}

	fclose(fp);
	printf("taskid in fprint = %d\n", taskid);
}


void prtdat(int nx, int ny, float *u1, char *fnam) {
int ix, iy;
FILE *fp;

fp = fopen(fnam, "w");
for (iy = ny-1; iy >= 0; iy--) {
  for (ix = 0; ix <= nx-1; ix++) {
    fprintf(fp, "%6.1f", *(u1+ix*ny+iy));
    if (ix != nx-1) 
      fprintf(fp, " ");
    else
      fprintf(fp, "\n");
    }
  }
fclose(fp);
}


