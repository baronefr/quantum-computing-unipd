lib_core :
	gfortran -c lib/mydebugger.f90

lib_build : lib_core
	gfortran -c lib/matrix_r8.f90
	gfortran -c ising.f90

rg : lib_build
	@echo " > exr: HW8 (RG)"
	@echo " just run, the parameters are defined inside exr_rg.f90"
	gfortran -c lib/rg.f90
	gfortran mydebugger.o matrix_r8.o rg.o ising.f90 exr_rg.f90 -o rg.x -lblas -llapack

dmrg : lib_build
	@echo " > exr: HW8 (DMRG)"
	@echo " just run, the parameters are defined inside exr_dmrg.f90"
	gfortran -c lib/idmrg.f90
	gfortran mydebugger.o matrix_r8.o idmrg.o ising.f90 exr_dmrg.f90 -o dmrg.x -lblas -llapack -O5

clean:
	@echo " > cleaning compile files..."
	rm -f *.o
	rm -f *.mod
	rm -f lib/*.mod

cleanx: clean
	@echo " > cleaning exe ..."
	rm -f *.x
