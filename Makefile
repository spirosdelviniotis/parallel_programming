## -*- Makefile -*-
##
## User: Spiros Delviniotis
##
## This file is generated automatically.
##


#### Compiler and tool definitions shared by all build targets #####
main: main.c
	mpicc -o main main.c

clean:
	rm -f main
