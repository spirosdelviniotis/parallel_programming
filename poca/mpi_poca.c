#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#else
int omp_get_thread_num(void) { return 0; }
int omp_get_num_threads(void) { return 1; }
#endif

#define NXPROB      6                 /* x dimension of problem grid */
#define NYPROB      6                 /* y dimension of problem grid */
#define STEPS       100                /* number of time steps */
#define MAXWORKER   8                  /* maximum number of worker tasks */
#define MINWORKER   3                  /* minimum number of worker tasks */
#define BEGIN       1                  /* message tag */
#define LTAG        2                  /* message tag */
#define RTAG        3                  /* message tag */
#define UTAG		4
#define DTAG		5
#define NONE        -1                  /* indicates no neighbor */
#define DONE        4                  /* message tag */
#define MASTER      0                  /* taskid of first process */

struct Parms {
  float cx;
  float cy;
} parms = {0.1, 0.1};

int main (int argc, char *argv[])
{
void inidat(), prtdat(), inner_update(), outer_update(), fprint();
								/* array for grid */
int	taskid,                     /* this task's unique id */
	numworkers,                 /* number of worker processes */
	numtasks,                   /* number of tasks */
	up,down,left,right,        	/* neighbor tasks */
	rc,
	ix,iy,iz,it,            	/* loop variables */
	d,x,y;						/* block size variables */

int periods[2]={0,0}, 			/* cartesian topology variables */
	dims[2],
	reorder=1,
	size ; 						/* size of dimmensions for cartesian topology */
double 	t1 = 0,					/* MPI time variables */
	t2 = 0,
	global = 0,				/* time variables */
	local  = 0;	

	FILE *fp[9];
	int temp;

	for(temp = 0; temp<9; temp++){
		char str[10];
		sprintf(str, "%d", temp);
		fp[temp] = fopen(str, "w");
	}


	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
	
	numworkers = numtasks;
	size = sqrt(numworkers);
	
	dims[0] = size;
	dims[1] = size;
	
	MPI_Comm newcomm;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &newcomm);

	MPI_Barrier(newcomm);

	MPI_Cart_shift(newcomm , 0 , 1,&up,&down);
	MPI_Cart_shift(newcomm , 1 , 1,&left,&right);

	if(NXPROB*NYPROB% numworkers != 0 ){
	   MPI_Abort(newcomm, rc);
	   exit(1);

	}
	d = sqrt(NXPROB*NYPROB/numworkers);
	x = d+2;
	y = d+2;

	float *u;	

	u = malloc((2*x*y)*sizeof(float));

	for(iz=0;iz<2;iz++){
	   for(ix=0;ix<x;ix++){
		   for(iy=0;iy<y;iy++){
			*(u+iz*x*y+ix*y+iy) = 0.0;
		   }
	   }
	}

	inidat(d,d,y,(u+y+1));
     
	iz = 0;
	MPI_Request 	req[8];
	MPI_Status   status[8];

	MPI_Datatype COL_INT;
	MPI_Type_vector(d,1,y,MPI_FLOAT,&COL_INT);
	MPI_Type_commit(&COL_INT);
	t1 = MPI_Wtime();
      
	for (it = 1; it <= STEPS; it++){

		if (up >= 0){								//up
			MPI_Isend(u+iz*x*y+y+1, d , MPI_FLOAT, up,DTAG, newcomm ,&req[0]);
			MPI_Irecv(u+iz*x*y+1, d , MPI_FLOAT, up,UTAG, newcomm, &req[1]);
		}

		if (down >= 0){							//down
			MPI_Isend(u+iz*x*y+d*y+1, d , MPI_FLOAT, down, UTAG, newcomm,&req[2]);
			MPI_Irecv(u+iz*x*y+(d+1)*y+1, d , MPI_FLOAT, down,DTAG,newcomm,&req[3]);
		}

		if (left >= 0){								//left
			MPI_Isend(u+iz*x*y+y+1, 1, COL_INT, left,RTAG, newcomm,&req[4]);
			MPI_Irecv(u+iz*x*y+y, 1, COL_INT, left,LTAG, newcomm, &req[5]);
		}

		if (right >= 0 ){							//right
			MPI_Isend(u+iz*x*y+y+d, 1, COL_INT, right,LTAG, newcomm,&req[6]);
			MPI_Irecv(u+iz*x*y+y+d+1, 1, COL_INT , right,RTAG, newcomm,&req[7]);
		}

		inner_update(2, d-2, 2, u+iz*x*y, u+(1-iz)*x*y);

		if(up >= 0){
			MPI_Wait(&req[0],&status[0]);
			MPI_Wait(&req[1],&status[1]);
		}
		if(down >= 0){
			MPI_Wait(&req[2],&status[2]);
			MPI_Wait(&req[3],&status[3]);
		}
		if(left >= 0){
			MPI_Wait(&req[4],&status[4]);
			MPI_Wait(&req[5],&status[5]);
		}
		if(right >= 0){
			MPI_Wait(&req[6],&status[6]);
			MPI_Wait(&req[7],&status[7]);
		}

		outer_update(1,d,1, u+iz*x*y, u+(1-iz)*x*y);
		iz = 1 - iz;
	}
	
	t2 = MPI_Wtime();

	fprint(fp[taskid], x, y, u, 1, 1, d, taskid);

	local = t2 - t1;

	MPI_Reduce(&local,&global,1,MPI_DOUBLE,MPI_MAX,0,newcomm);
	MPI_Barrier(newcomm);

	if(taskid == 0){
	  printf("Total time elapsed = %lf \n",global);
	}

	free(u);

	MPI_Type_free(&COL_INT);
	MPI_Finalize();
}



/**************************************************************************
 *  subroutine update
 ****************************************************************************/


void update_calc(int ix, int iy, int y, float *u1, float *u2) {

	*(u2+ix*y+iy) = *(u1+ix*y+iy)  +
		parms.cx * (*(u1+(ix+1)*y+iy) +
        *(u1+(ix-1)*y+iy) -
        2.0 * *(u1+ix*y+iy)) +
        parms.cy * (*(u1+ix*y+iy+1) +
        *(u1+ix*y+iy-1) -
        2.0 * *(u1+ix*y+iy));

}

void inner_update(int start, int end, int ny, float *u1, float *u2)
{
	int i,j;
	printf("	>> START of inner_update with: end = %d.\n", end);
	//#pragma omp parallel
    {
		//int tid = omp_get_thread_num();
        //int total = omp_get_num_threads();
        //#pragma omp for
		for (i = 2; i <= 1+end; i++) {
			for (j = 2; j <= 1+end; j++){
				update_calc(i, j, end+4, u1, u2);
				//printf("	>>  i = %d, j = %d.\n",i,j);
			}
		}
	//printf("Inner Update END!\n");
	}
}

void outer_update(int start, int end, int ny, float *u1, float *u2)
{
	int i;
	for (i = 1; i <= end; i++) {
		update_calc(1, i, end+2, u1, u2);
		update_calc(end, i, end+2, u1, u2);
		update_calc(i, 1, end+2, u1, u2);
		update_calc(i, end, end+2, u1, u2);
	}
}

/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny,int y , float *u) {
int ix, iy;
for (ix = 0; ix <= nx-1; ix++)
  for (iy = 0; iy <= ny-1; iy++)
     *(u+ix*(y)+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1)) + 10;
}


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

