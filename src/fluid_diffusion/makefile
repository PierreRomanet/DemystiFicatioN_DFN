# GNU compilation

LIBS=-L/usr/local/lib -Llibhdf5_fortran.a

F90=/Users/pierre/Dropbox/bib_postdoc_Pierre_Romanet/Softwares/Fluid_diffusion/2D/library/hdf5-1.14.2/hdf5/bin/h5fc

F90FLAGS=$(INCPATH) -ffree-line-length-none -O3 -w -fopenmp -c
CFLAGS=$(INCPATH)


SRC=src
DST=build


OBJ = variables.o read_input.o initialisation.o FVM.o ode_bs.o solver.o

software: $(OBJ) master.f90
	$(F90) -ffree-line-length-none -O3 -w -fopenmp $(INCPATH) $(LIBS) $(OBJ) master.f90 -o ../$(DST)/DiffUnion 

read_input.o: variables.o read_input.f90
	$(F90) $(F90FLAGS) $(FLAG) read_input.f90

solver.o: variables.o FVM.o ode_bs.o solver.f90
	$(F90) $(F90FLAGS) solver.f90
	
ode_bs.o: variables.o FVM.o ode_bs.f90
	$(F90) $(F90FLAGS) ode_bs.f90

initialisation.o: variables.o initialisation.f90
	$(F90) $(F90FLAGS) initialisation.f90 

%.o: %.f90
	$(F90) $(F90FLAGS)   -c $< -o  $@

%.o: %.f
	$(F90) $(F90FLAGS)   -c $< -o  $@

variables.o: variables.f90
	$(F90) $(F90FLAGS) variables.f90

clean:
	rm *.o *.mod ../$(DST)/DiffUnion
