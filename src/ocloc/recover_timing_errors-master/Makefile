SHELL = /bin/bash	# Needs to be defined

# This makefile takes care of the compilation of the following program:
# ./SRC/TIMIING_ERR_INV.F90 		-> ./BIN/Timing_err_inv

# Some notes:
# -c flag makes sure the linking is not done until all the source scripts providing modules are compiled.
# Never only remove the .mod files in the MODULES directory! they are used...
# If you do remove them, make sure you also remove the corresponding object .o files, because the make utility checks if these oject files are updated (or not)

######################################################################################################################
####################################################### COMPILER #####################################################
######################################################################################################################

FOR=gfortran

######################################################################################################################
############################################## MAKE MAIN PROGRAMMES ##################################################
######################################################################################################################

# Directories
M = ./MOD
S = ./SRC
B = ./BIN
O = ./OBJ

##### Libraries and their locations (Not all libraries are used. Libraries, and their locations, depends on your cluster/machine)
SACLIBDIR=-L$(SACHOME)/lib
SACLIBS=-lsac -lsacio
SACMODS=-I$(SACHOME)/include

##### Flags (Not all flags are used. Flags are also compiler specific)
DEBUGFLAGS=-g -Wall -fbounds-check -fbacktrace -fmax-errors=1 #-fcheck=all ##-fno-range-check
OPTIONS=-ffree-line-length-none
OPTFLAGS=-funroll-all-loops -fpic -O3  #-fcheck=all
OPENMPFLAGS=-fopenmp
CPPFLAGS=-traditional-cpp

###################################### DEPENDENCIES OF MAIN PROGRAMS ON OBJECT FILES ###################################################

TIMING_ERR_INV_OBJECTS = $O/File_operations.o $O/Array_operations.o $O/Strings.o $O/Precision.o $O/Snr_routines.o $O/Return_taper.o $O/Cubspl.o $O/Ppvalu.o $O/Interv.o

###################################### DEPENDENCIES OF MAIN PROGRAMS ON OBJECT FILES ###################################################

$B/Timing_err_inv		: $S/TIMING_ERR_INV.F90 ${TIMING_ERR_INV_OBJECTS}
			${FOR} $< ${TIMING_ERR_INV_OBJECTS} ${OPTIONS} ${SACLIBDIR} ${SACLIBS} ${SACMODS} ${DEBUGFLAGS} -I$M -o $@

$O/Array_operations.o		: $S/Array_operations.F90
			${FOR} -c $< ${DEBUGFLAGS} -J$M -o $@

$O/File_operations.o		: $S/File_operations.F90 $O/Strings.o
			${FOR} -c $< ${DEBUGFLAGS} -J$M -o $@

$O/Strings.o			: $S/Strings.F90 $O/Precision.o
			${FOR} -c $< ${DEBUGFLAGS} -J$M -o $@

$O/Precision.o			: $S/Precision.F90
			${FOR} -c $< ${DEBUGFLAGS} -J$M -o $@

$O/Return_taper.o		: $S/Return_taper.F90
			${FOR} -c $< ${DEBUGFLAGS} -J$M -o $@

$O/Snr_routines.o		: $S/Snr_routines.F90
			${FOR} -c $< ${DEBUGFLAGS} -J$M -o $@

$O/Cubspl.o:			 $S/Cubspl.f90 $O/Ppvalu.o
			${FOR} -c $< ${OPTFLAGS} -o $@

$O/Ppvalu.o:			 $S/Ppvalu.f90 $O/Interv.o
			${FOR} -c $< ${OPTFLAGS} -o $@

$O/Interv.o:		 	 $S/Interv.f90
			${FOR} -c $< ${OPTFLAGS} -o $@

clean				:
			rm $B/* $O/* $M/*
