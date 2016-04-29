## -*- Makefile -*-
##
## User: Spiros Delviniotis
##
## This file is generated automatically.
##


#### Compiler and tool definitions shared by all build targets #####
ring: ring.c
        mpicc -o ring ring.c

clean:
        rm -f ring
