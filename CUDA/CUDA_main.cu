#ifndef CUDA_MAIN_H
#define CUDA_MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include "common/inc/lcutil.h"
#include "common/inc/timestamp.h"

#define NXPROB      1000                 /* x dimension of problem grid */
#define NYPROB      1000                 /* y dimension of problem grid */
#define STEPS       100                /* number of time steps */


int main(int argc, char* argv[])
{
  	float	*table_1,
		*table_2,
		*table_host;
	int	it,
		iz,
		size;
	timestamp start_time;

	size = NXPROB*NYPROB*sizeof(float);
	
	CUDA_SAFE_CALL(cudaMalloc((void**)&table_1,(long) size));
	CUDA_SAFE_CALL(cudaMalloc((void**)&table_2,(long) size));

	table_host = (float*)malloc(size);
	if (table_host == NULL) {
		printf("Main ERROR: Allocation memory.\n");
		return(EXIT_FAILURE);
	}
	
	/* Initialize table_host */
	inidat(NXPROB, NYPROB, table_host);

	CUDA_SAFE_CALL(cudaMemcpy(table_1, table_host, size, cudaMemcpyHostToDevice));	
	CUDA_SAFE_CALL(cudaMemcpy(table_2, table_host, size, cudaMemcpyHostToDevice));
	
	iz = 0; 
	start_time = getTimestamp();

	for (it = 1; it <= STEPS; it++)
	{       
		if(iz==0){
			update<<<NYPROB,NXPROB>>>(NYPROB, table_1, table_2);
		}else{
			update<<<NYPROB,NXPROB>>>(NYPROB, table_2, table_1);
		}		
		iz = 1 - iz;
        		
		CUDA_SAFE_CALL(cudaThreadSynchronize());		
	}
			
	printf("Time elapsed = %6.4lf in ms\n",getElapsedtime(start_time));

	CUDA_SAFE_CALL(cudaMemcpy(table_host, table_2, NXPROB*NYPROB*sizeof(float), cudaMemcpyDeviceToHost));
	
	CUDA_SAFE_CALL(cudaFree(table_1) );	
	CUDA_SAFE_CALL(cudaFree(table_2) );
	free(table_host);	
}


__global__ void update(int ny, float *u1, float *u2)
{	
	struct Parms { 
		float cx;
		float cy;
	} parms = {0.1, 0.1};
//	
//	??????????????????????????
//	int index = blockIdx.x*blockDim.x + threadIdx.x;
//	int ix = index / (ny-2) + 1;
//	int iy = index % (ny-2) + 1;
 
  	
	*(u2+ix*ny+iy) =  *(u1+ix*ny+iy)  + 
                  	parms.cx * (*(u1+(ix+1)*ny+iy) +
			*(u1+(ix-1)*ny+iy) - 
			2.0 * *(u1+ix*ny+iy)) +
			parms.cy * (*(u1+ix*ny+iy+1) +
			*(u1+ix*ny+iy-1) - 
			2.0 * *(u1+ix*ny+iy));	
	
	__syncthreads();
}


void inidat(int nx, int ny, float *u) {
	int ix, iy;

	for (ix = 0; ix <= nx-1; ix++) {
		for (iy = 0; iy <= ny-1; iy++) {
			*(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
		}
	}
}
#endif	// CUDA_MAIN_H

