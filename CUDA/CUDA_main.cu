#ifndef CUDA_MAIN_H
#define CUDA_MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include "common/inc/lcutil.h"
#include "common/inc/timestamp.h"

#define NXPROB      3600                /* x dimension of problem grid */
#define NYPROB      3600                /* y dimension of problem grid */
#define STEPS       10000                /* number of time steps */


extern "C" double calculation_GPU();


int main(int argc,char *argv[])
{
	double time_lapse;
	
	time_lapse = calculation_GPU();
	
	printf("Total time elapsed for:\n\tSteps: %d\n\tTable [%d]x[%d] = %lf ms. \n", STEPS, NXPROB, NYPROB, time_lapse);
}


/* subroutine inidat */
void inidat(int nx, int ny, float *u) {
	int ix, iy;

	for (ix = 0; ix <= nx-1; ix++) 
		for (iy = 0; iy <= ny-1; iy++)
			*(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
}


/* This code runs on device */
__global__ void update(int y, float *u1, float *u2)
{	
	int ix = 0;
	int iy = 0;
	
	struct Parms { 
		float cx;
		float cy;
	} parms = {0.1, 0.1};
	
	/* Get coordinators */
	ix = blockIdx.x*blockDim.x+threadIdx.x;
	iy = blockIdx.y*blockDim.y+threadIdx.y;
  	
	*(u2 + ix * y + iy) = *(u1 + ix * y + iy) +
		parms.cx * (*(u1 + (ix + 1) * y + iy) +
		*(u1 + (ix - 1) * y + iy) -
		2.0 * *(u1 + ix * y + iy)) +
		parms.cy * (*(u1 + ix * y + iy + 1) +
		*(u1 + ix * y + iy - 1) -
		2.0 * *(u1 + ix * y + iy));
	
	__syncthreads();
}


/* This code runs on host */
extern "C" double calculation_GPU()
{
  	float	*table_1,
		*table_2,
		*table_host;
	int	it,
		size,
		iz = 0;
	timestamp start_time;
	
	/* Calculate size and allocate memory in device for the arrays */
	size = NXPROB*NYPROB*sizeof(float);
	CUDA_SAFE_CALL(cudaMalloc((void**)&table_1,(long) size));
	CUDA_SAFE_CALL(cudaMalloc((void**)&table_2,(long) size));

	/* Allocate memory in host for the array */
	table_host = (float*)malloc(size);
	if (table_host == NULL) {
		printf("Main ERROR: Allocation memory.\n");
		exit(-1);
	}
	
	/* Initialize table_host with zeros and then call inidat*/
	memset(table_host, 0, NXPROB*NYPROB*sizeof(float));
	inidat(NXPROB, NYPROB, table_host);
	
	/* Copy table_1 and table_2 to GPU */
	CUDA_SAFE_CALL(cudaMemcpy(table_1, table_host, size, cudaMemcpyHostToDevice));	
	CUDA_SAFE_CALL(cudaMemcpy(table_2, table_host, size, cudaMemcpyHostToDevice));
	 
	/* Create N blocks of N threads each */
	dim3 NumberOfThreads(NXPROB);			
	dim3 NumberOfBlocks(NYPROB);
	
	/* Start the Clock */
	start_time = getTimestamp();
	
	/* Go! */
	for (it = 1; it <= STEPS; it++)
	{       
		if ( iz==0 ){
			update<<<NumberOfBlocks,NumberOfThreads>>>(NYPROB, table_1, table_2);
		}
		else {
			update<<<NumberOfBlocks,NumberOfThreads>>>(NYPROB, table_2, table_1);
		}
		
		/* Swap table pointers for next loop */
		iz = 1 - iz;
        	
		/* Sync Cuda Threads */
		CUDA_SAFE_CALL(cudaThreadSynchronize());		
	}
	
	/* Copy table with results to table_host from GPU */
	CUDA_SAFE_CALL(cudaMemcpy(table_host, table_2, NXPROB*NYPROB*sizeof(float), cudaMemcpyDeviceToHost));
	
	/* Free Resources */
	CUDA_SAFE_CALL(cudaFree(table_1) );	
	CUDA_SAFE_CALL(cudaFree(table_2) );
	free(table_host);

	return getElapsedtime(start_time);
}


#endif	// CUDA_MAIN_H

