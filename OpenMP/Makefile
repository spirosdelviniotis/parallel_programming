## -*- Makefile -*-
##
## User: Spiros Delviniotis
##
##


#### Compiler and tool definitions shared by all build targets #####
OpenMP_main: OpenMP_main.c
	#mpicc -fopenmp -o OpenMP_main OpenMP_main.c -lm
	mpicc -o Serial_OpenMP_main OpenMP_main.c -lm
clean:
	rm -f OpenMP_main
