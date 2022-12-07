lib_core :
	gfortran -c lib/mydebugger.f90
	gfortran -c lib/mytools.f90

lib_build : lib_core
	gfortran -c lib/matrix_r8.f90 -llapack
	gfortran -c lib/harmonic.f90

init : lib_build
	@echo " > preparation of initial state"
	gfortran mydebugger.o matrix_r8.o mytools.o harmonic.o  get_psi0.f90 -o init.x -llapack
	@echo " > call with :     ./init.x  *Nx*  data/psi0_"

exr : lib_build
	@echo " > exr: time-dependent harmonic oscillator"
	gfortran mydebugger.o matrix_r8.o trotter_suzuki.f90 exr.f90 -o exr.x -llapack -lfftw3 -lm
	@echo " > call with :     ./exr.x   *Nt*  *tmax*  *tau*  data/psi_"

clean:
	@echo " > cleaning compile files..."
	rm -f *.o
	rm -f *.mod
	rm -f lib/*.mod

cleanx: clean
	@echo " > cleaning exe ..."
	rm -f *.x