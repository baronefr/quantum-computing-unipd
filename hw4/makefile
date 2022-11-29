lib_core :
	gfortran -c lib/mydebugger.f90
	gfortran -c lib/mytools.f90

lib_build : lib_core
	gfortran -c lib/matrix_c16.f90 -llapack
	gfortran -c lib/matrix_r8.f90 -llapack

exr : lib_build
	@echo " > exr: time independent harmonic oscillator"
	@echo "          -> tridiagonal hamiltonian"
	gfortran mydebugger.o matrix_r8.o mytools.o harmonic.f90 exr-tri.f90  -o exr.x -llapack -O8

extra : lib_build
	@echo " > exr: time independent harmonic oscillator"
	@echo "          -> pentadiagonal hamiltonian"
	gfortran mydebugger.o matrix_r8.o mytools.o harmonic.f90 exr-penta.f90  -o exr.x -llapack -O8

clean:
	@echo " > cleaning compile files..."
	rm -f *.o
	rm -f *.mod
	rm -f lib/*.mod

cleanx: clean
	@echo " > cleaning exe ..."
	rm -f *.x