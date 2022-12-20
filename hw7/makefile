lib_core :
	gfortran -c lib/mydebugger.f90

lib_build : lib_core
	gfortran -c lib/matrix_r8.f90 -llapack
	gfortran -c lib/matrix_r8_packed.f90 -llapack

exr : lib_build
	@echo " > exr: HW7"
	@echo " just run, the parameters are defined inside exr.f90"
	gfortran mydebugger.o matrix_r8_packed.o ising.f90 exr.f90 -o exr.x -llapack

clean:
	@echo " > cleaning compile files..."
	rm -f *.o
	rm -f *.mod
	rm -f lib/*.mod

cleanx: clean
	@echo " > cleaning exe ..."
	rm -f *.x
