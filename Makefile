## -*- Makefile -*-
##
## User: Spiros Delviniotis
##
## This file is generated automatically.
##


#### Compiler and tool definitions shared by all build targets #####
exe: main.c
	mpicc -o ring main.c

clean:
	rm -f exe
