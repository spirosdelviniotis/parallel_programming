#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/MPI/MPI_main.o \
	${OBJECTDIR}/MPI_Heat_2D/mpi_heat2D.o \
	${OBJECTDIR}/OpenMP/OpenMP_main.o \
	${OBJECTDIR}/OpenMP/Pure_OpenMP_main.o \
	${OBJECTDIR}/Serial_heat2D.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/parallel_programming_di

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/parallel_programming_di: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.c} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/parallel_programming_di ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/MPI/MPI_main.o: MPI/MPI_main.c 
	${MKDIR} -p ${OBJECTDIR}/MPI
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MPI/MPI_main.o MPI/MPI_main.c

${OBJECTDIR}/MPI_Heat_2D/mpi_heat2D.o: MPI_Heat_2D/mpi_heat2D.c 
	${MKDIR} -p ${OBJECTDIR}/MPI_Heat_2D
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MPI_Heat_2D/mpi_heat2D.o MPI_Heat_2D/mpi_heat2D.c

${OBJECTDIR}/OpenMP/OpenMP_main.o: OpenMP/OpenMP_main.c 
	${MKDIR} -p ${OBJECTDIR}/OpenMP
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/OpenMP/OpenMP_main.o OpenMP/OpenMP_main.c

${OBJECTDIR}/OpenMP/Pure_OpenMP_main.o: OpenMP/Pure_OpenMP_main.c 
	${MKDIR} -p ${OBJECTDIR}/OpenMP
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/OpenMP/Pure_OpenMP_main.o OpenMP/Pure_OpenMP_main.c

${OBJECTDIR}/Serial_heat2D.o: Serial_heat2D.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Serial_heat2D.o Serial_heat2D.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/parallel_programming_di

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
