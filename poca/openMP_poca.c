#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#else
int omp_get_thread_num(void) { return 0; }
int omp_get_num_threads(void) { return 1; }
#endif

#define	NXPROB		20
#define	NYPROB		20
#define STEPS		1000

struct Parms { 
  float cx;
  float cy;
} parms = {0.1, 0.1};

long long
timeval_diff(struct timeval *difference,
             struct timeval *end_time,
             struct timeval *start_time
            )
{
  struct timeval temp_diff;

  if(difference==NULL)
  {
    difference=&temp_diff;
  }

  difference->tv_sec =end_time->tv_sec -start_time->tv_sec ;
  difference->tv_usec=end_time->tv_usec-start_time->tv_usec;

  /* Using while instead of if below makes the code slightly more robust. */

  while(difference->tv_usec<0)
  {
    difference->tv_usec+=1000000;
    difference->tv_sec -=1;
  }

  return 1000000LL*difference->tv_sec+
                   difference->tv_usec;

}
/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, float *u) {
int ix, iy;

for (ix = 0; ix <= nx-1; ix++) 
  for (iy = 0; iy <= ny-1; iy++)
     *(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
}
/**************************************************************************
 *  subroutine update
 ****************************************************************************/
void update(int start, int end, int ny, float *u1, float *u2)
{
   int ix, iy;
   
   
   #pragma omp parallel
    {
        //int tid = omp_get_thread_num();
        //int total = omp_get_num_threads();
		
        #pragma omp for   
		//printf("This is thread %d of %d\n", tid, total);
		   for (ix = start; ix <= end; ix++){ 
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
int main(int arg , char *argv[])
{
	float *u;
	int ix,iy,iz,it;	
	
	u = malloc((2*NXPROB*NYPROB)*sizeof(float));

	for(iz=0;iz<2;iz++){
	   for(ix=0;ix<NXPROB;ix++){
		   for(iy=0;iy<NYPROB;iy++){
			*(u+iz*NXPROB*NYPROB+ix*NYPROB+iy) = 0.0;			
		   }
	   }
	}
	
	inidat(NXPROB, NYPROB, u);
		
	iz = 0;	
	struct timeval start, end;
	gettimeofday(&start, NULL);
	for (it = 1; it <= STEPS; it++){	
		
		update(1,NXPROB-2,NYPROB,u+iz*NXPROB*NYPROB, u+(1-iz)*NXPROB*NYPROB);
		
		iz = 1 - iz;
	}	
	gettimeofday(&end, NULL);
	
	printf("Iterations: %d , Time elapsed: %lf seconds\n",
			STEPS,((end.tv_sec  - start.tv_sec) * 1000000u + 
         end.tv_usec - start.tv_usec) / 1.e6);
	free(u);
}
