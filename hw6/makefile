lib_core :

lib_build : lib_core
	gfortran -c lib/density_matrix.f90 -fdefault-integer-8
	gfortran -c lib/quantum_states.f90 -fdefault-integer-8

exr : lib_build
	@echo " > exr: density matrix-test"
	gfortran density_matrix.o test_partial_trace.f90 -o test.x -fdefault-integer-8
	@echo " > call with :     ./test.x   *filename*   *system*"

bench : lib_build
	@echo " > exr: efficiency test"
	gfortran quantum_states.o efficiency.f90 -o bench.x -llapack -fdefault-integer-8

clean:
	@echo " > cleaning compile files..."
	rm -f *.o
	rm -f *.mod
	rm -f lib/*.mod

cleanx: clean
	@echo " > cleaning exe ..."
	rm -f *.x