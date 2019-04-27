## Simple make file, It does NOT check for interdependencies of modules.
## It probably is best to do make clean; make if you change a module.


## List main program module last
## 

SRC_DIRS ?=   #./src
#BUILD_DIR ?= ./build
BUILD_DIR ?=  #./src
RUN_DIR ?= ../run


#SOURCES:= $(shell find src -name '*.f90')

SOURCES = precise.f90 \
	constants.f90 \
	analytic.f90 \
	input.f90 \
	initial_conditions.f90 \
	geo.f90 \
	inv.f90 \
	solver.f90


## Define name of main program
PROGRAM = test


# Compiler
FF = gfortran


# Delete program
# Linux
RM = rm -f
# DOSe
#RM = del


## Compiler options
#
# for the Intel Fortran 90 compiler
# CFLAGS = -c -fast -heap-arrays
#
# for the gfortran or g95 compiler 
# 	-g for gdb
# 	-02 is safer  # -03 is maybe faster
CFLAGS = -c -O3 


#Linker Options
#
# for the Intel Fortran 90 compiler
# LDFLAGS = -fast -heap-arrays
#
# for the gfortran or g95 compiler
LDFLAGS = -O2


## messing with directory structure below this line
#OBJECTS = $(SRC_DIRS)/$(SOURCES:.f90=.o)
OBJECTS := $(SOURCES:.f90=.o)

all: $(SOURCES) $(RUN_DIR)/$(PROGRAM)


$(RUN_DIR)/$(PROGRAM): $(OBJECTS)
	$(FF) $(LDFLAGS) $(OBJECTS) -o $@


$(OBJECTS) : %.o : %.f90
	$(FF)  $(CFLAGS) $< -o $@ 

#$(OBJECTS): %.o: %.f90
#	$(FF) $(CFLAGS) $< -o $(BUILD_DIR)/$@


clean:
	$(RM) *.o *.mod 

realclean:
	$(RM) *.o  *.mod $(RUN_DIR)/$(PROGRAM)