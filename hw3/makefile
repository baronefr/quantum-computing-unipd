lib_core :
	gfortran -c lib/mydebugger.f90
	gfortran -c lib/mytools.f90

lib_build : lib_core
	gfortran -c lib/matrix_c16.f90
	gfortran -c lib/matmul_loops.f90



exr1 : lib_build
	@echo " > exr 1"
	gfortran -c lib/matmul_loops.f90
	gfortran matmul_loops.o exr1_launcher.f90 -o exr1.x -llapack

exr2 : lib_build
	@echo " > exr 2"
	gfortran -c lib/matrix_c16.f90
	gfortran mydebugger.o mytools.o matrix_c16.o exr2.f90 -o exr2.x -llapack
	@echo " [info]  exr 2 args :     matrix_size(int) [ output_file(string) ]"

exr3 : lib_build
	@echo " > exr 3"
	gfortran -c lib/matrix_c16.f90
	gfortran mydebugger.o mytools.o matrix_c16.o exr3.f90 -o exr3.x -llapack
	@echo " [info]  exr 3 args :     matrix_size(int) iter(int) matrix_type(char*1) [ output_file(string) ]"



clean:
	@echo " > cleaning ..."
	rm -f *.o
	rm -f *.mod
	rm -f *.x