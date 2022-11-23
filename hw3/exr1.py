#!/bin/python3

# ========================================
#  QIC  >>  assignment 3 / exr 1
#  UniPD, AY 2022/23, Physics of Data
# ----------------------------------------
#   coder : Barone Francesco
#   dated : 09 November 2022
# ========================================


import os
import sys
import math
from subprocess import run



# I/O --------------------------------------------------------------
BENCHFILES_ORIGIN = "data/O{}.bench"

# TEST PARAMETERS --------------------------------------------------

# define the optimization levels to test
OPTI_LEVELS = [ 0, 1, 2, 5]

# define the class of iterables to benchmark
def benchmark_iterable(opti, bfile) -> list:
    Nmin, Nmax = 100, 10000

    arr = [];  dim = Nmin; 
    while dim < Nmax:
        if dim <= 2000: rep = 3
        else:          rep = 1
        arr.append( (dim,rep,bfile) )

        dim += int( 10**math.floor( math.log10(dim) ) )
        
    return arr


# initial checks before executing the code
def preliminaries() -> None:
    # create directory for benchmarks, if not exists
    directory = os.path.dirname(BENCHFILES_ORIGIN)
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        print("target path already exists")



# internals (do not edit) ------------------------------------------
EXE_FILENAME = "exr1.x"
OPTI_BENCHFILES = [ BENCHFILES_ORIGIN.format(str(i)) for i in OPTI_LEVELS ]


# define the actions needed to compile the executable
def make(opt_level = 0) -> None:
    print(" [make] making fortran program with args : ", locals())

    # instructions to be executed to compile the exe
    COMPILE_COMMANDS = [ "gfortran -c lib/matmul_loops.f90 -O{}",
                         "gfortran exr1_launcher.f90 matmul_loops.o -o {}" ]

    os.system( COMPILE_COMMANDS[0].format(str(opt_level)) )
    os.system( COMPILE_COMMANDS[1].format(EXE_FILENAME)   )

    print(" [make] checking exe ... ", end='')
    if os.path.exists(EXE_FILENAME):
        print("ok".format(EXE_FILENAME))
    else:
        print("failed")
        raise Exception("[FATAL] {} does not exist".format(EXE_FILENAME))


# define the housekeeping operations at the end of iterables
def clean() -> None:
    os.system("rm {}".format(EXE_FILENAME))
    os.system("rm matmul_launcher.o matmul_loops.o")





# ------------------------------------------------------------------
#  MAIN
# ------------------------------------------------------------------
if __name__ == "__main__":
    #  process flags
    argv = sys.argv[1:]
    if "--dry-run" in argv:  dry_run = True
    else:  dry_run = False

    #  DO PRELIMINARIES
    try:   preliminaries()
    except Exception as err:
        print(err)
        sys.exit(1) # exit with error code 1

    #  EXECUTE THE BENCHMARKS
    print("doing benchmarks")
    try: 
        for opti, bfile in zip(OPTI_LEVELS, OPTI_BENCHFILES):
            # compile with the optimization flag
            make(opt_level = opti)
            print("target file is", bfile)

            for eargs in benchmark_iterable(opti, bfile):
                action = "./" + EXE_FILENAME + " " + " ".join( str(ar) for ar in eargs)

                if dry_run: print('[sys] ', action )
                else:
                    p = run( action.split() )
                    if p.returncode != 0:
                        raise Exception( f'runtime error : { p.returncode }' )

            clean() # housekeeping

    except Exception as err:
        print(err)
        sys.exit(2) # exit with error code 2

    sys.exit(0)