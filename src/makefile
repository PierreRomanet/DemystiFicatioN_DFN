# Make file for the software
F90 = h5fc

# Libraries
LIBS=-L../library -lhfmm2d -llapack -lrefblas

# Generation des objets
FFLAGS =

# Generation de l'executable
LDFLAGS = -fopenmp

FLAG = -c -O1 -Wall 
OBJ_HMAT = initialisation_hmat.o ACA_algorithm.o block_cluster_tree.o cluster_tree.o get_info_matrices.o H_matrix_construction.o H_matrix_recompression.o calculation_hmat.o eval_kernel_hmatrix.o

LIB = libhfmm2d.a
OBJ = variables.o configuration.o initialisation.o fluid_diffusion.o RateStateAgingModeIII.o RateStateAgingModeII.o hmatrixModeII.o hmatrixModeIII.o fmmModeIII.o directModeIII.o directModeII.o nrtype.o nrutil.o ode_bs.o solver.o


software: $(OBJ) master.f90
	$(F90) -O3 $(LDFLAGS) $(OBJ_HMAT) $(OBJ) $(LIBS) master.f90 -o ../build/DemystiFicatioN 

initialisation.o: initialisation_hmat.o RateStateAgingModeIII.o RateStateAgingModeII.o hmatrixModeII.o hmatrixModeIII.o fmmModeIII.o directModeIII.o directModeII.o eval_kernel_hmatrix.o variables.o libhfmm2d.a initialisation.f90
	$(F90) $(FLAG) initialisation.f90 

solver.o: ode_bs.o fluid_diffusion.o solver.f90 
	$(F90) $(FLAG) solver.f90 
	
#-----------------------------------------------------------------------------------------
# ODE solver
ode_bs.o: nrtype.o nrutil.o ode_solver/ode_bs.f90
	$(F90) $(FLAG) ode_solver/ode_bs.f90
	
nrtype.o: ode_solver/nrtype.f90
	$(F90) $(FLAG) ode_solver/nrtype.f90

nrutil.o: ode_solver/nrutil.f90
	$(F90) $(FLAG) ode_solver/nrutil.f90

	

#-----------------------------------------------------------------------------------------
# FRICTION LAWS
RateStateAgingModeIII.o: fluid_diffusion.o friction/RateStateAgingModeIII.f90
	$(F90) $(FLAG) friction/RateStateAgingModeIII.f90

RateStateAgingModeII.o: fluid_diffusion.o friction/RateStateAgingModeII.f90
	$(F90) $(FLAG) friction/RateStateAgingModeII.f90

# RateStateAgingModeII_R.o: friction/RateStateAgingModeII_R.f90
#	$(F90) $(FLAG) friction/RateStateAgingModeII_R.f90
# 
# RateStateSlipModeII.o: friction/RateStateSlipModeII.f90
#	$(F90) $(FLAG) friction/RateStateSlipModeII.f90
# 
# RateStateSlipModeIII.o: friction/RateStateSlipModeIII.f90
#	$(F90) $(FLAG) friction/RateStateSlipModeIII.f90

#-----------------------------------------------------------------------------------------
# FINITE VOLUME METHOD
fluid_diffusion.o: fluid_diffusion.f90
	$(F90) $(FLAG) fluid_diffusion.f90




#-----------------------------------------------------------------------------------------
# KERNEL CALCULATION
hmatrixModeII.o: calculation_hmat.o kernel_calculation/hmatrixModeII.f90
	$(F90) $(FLAG) kernel_calculation/hmatrixModeII.f90
	
hmatrixModeIII.o: calculation_hmat.o kernel_calculation/hmatrixModeIII.f90
	$(F90) $(FLAG) kernel_calculation/hmatrixModeIII.f90

fmmModeIII.o: kernel_calculation/fmmModeIII.f90
	$(F90) $(FLAG) kernel_calculation/fmmModeIII.f90
	
directModeII.o: kernel_calculation/directModeII.f90
	$(F90) $(FLAG) kernel_calculation/directModeII.f90
	
directModeIII.o: kernel_calculation/directModeIII.f
	$(F90) $(FLAG) kernel_calculation/directModeIII.f


#-----------------------------------------------------------------------------------------
# H-matrix specific
initialisation_hmat.o: cluster_tree.o block_cluster_tree.o H_matrix_construction.o H_matrix_recompression.o hmatrix/initialisation_hmat.f90
	$(F90) $(FLAG) hmatrix/initialisation_hmat.f90

cluster_tree.o: hmatrix/cluster_tree.f90
	$(F90) $(FLAG) hmatrix/cluster_tree.f90

eval_kernel_hmatrix.o: hmatrix/eval_kernel_hmatrix.f90
	$(F90) $(FLAG) hmatrix/eval_kernel_hmatrix.f90

ACA_algorithm.o: get_info_matrices.o hmatrix/ACA_algorithm.f90
	$(F90) $(FLAG) hmatrix/ACA_algorithm.f90

block_cluster_tree.o: hmatrix/block_cluster_tree.f90
	$(F90) $(FLAG) hmatrix/block_cluster_tree.f90

get_info_matrices.o: hmatrix/get_info_matrices.f90
	$(F90) $(FLAG) hmatrix/get_info_matrices.f90

H_matrix_construction.o: ACA_algorithm.o hmatrix/H_matrix_construction.f90
	$(F90) $(FLAG) $(LIBS) hmatrix/H_matrix_construction.f90

H_matrix_recompression.o: hmatrix/H_matrix_recompression.f90
	$(F90) $(FLAG) $(LIBS) hmatrix/H_matrix_recompression.f90

calculation_hmat.o: get_info_matrices.o hmatrix/calculation_hmat.f90
	$(F90) $(FLAG) hmatrix/calculation_hmat.f90

%.o: %.f90
	$(F90) $(FLAG) -c $< -o	 $@
	
friction/%.o: friction/%.f90
	$(F90) $(FLAG) -c $< -o	 $@
	
hmatrix/%.o: hmatrix/%.f90
	$(F90) $(FLAG) -c $< -o	 $@

%.o: %.f
	$(F90) $(FLAG) -c $< -o	 $@
	
	
	
	
	
#-----------------------------------------------------------------------------------------
# Libraries (to modify later)
liblapack.a: blas.f zlapack.f clapack.f dlapack.f ilapack.f slapack.f install.f
	$(F90) $(FLAG) blas.f zlapack.f clapack.f dlapack.f ilapack.f slapack.f install.f
	ar rvs liblapack.a blas.o zlapack.o clapack.o dlapack.o ilapack.o slapack.o install.o

libhfmm2d.a: 
	cd ../library/fmm_library; make ; mv libhfmm2d.a ../libhfmm2d.a


#-----------------------------------------------------------------------------------------
clean:
	rm *.o *.mod ../library/fmm_library/*.o ../library/libhfmm2d.a
