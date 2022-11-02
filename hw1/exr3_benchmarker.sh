#!/bin/bash -e

# ========================================
#  QIC  >>  assignment 1 / exr 3
#  UniPD, AY 2022/23, Physics of Data
# ----------------------------------------
#   coder : Barone Francesco
#   dated : 21 October 2022
# ========================================

# This script iterates the execution of exr3, feeding as first argument
# the dimension of the square matrices to be multiplied.
# It is up to the .x to write the result of benchmark
# in file ./data/benchmarks.txt


optimize_level=(0 1 2 3 5)  # list of optimizations to test
for opt in "${optimize_level[@]}"
do

## --- repeat this block for various optimization flags
echo " [i] testing optimization $opt"
echo " -> gfortran -c matmul_loops.f90 exr3.f90 -O${opt}"

echo "#${opt}" > data/benchmarks.txt  # mark optimization in bench file

# compile
gfortran -c matmul_loops.f90 exr3.f90 -O${opt}
gfortran exr3.o matmul_loops.o -o exr3.x

# save benchmark file, if any
#mv data/benchmarks.txt data/benchmarks.bak

# echo current time, just to see how much it takes
starttime=$(date)
echo " [i] start at $starttime"

# call the program with different matrix sizes
for (( c=100; c<=2000; c+=100 )) 
do
   echo "-- dim = $c"
   ./exr3.x $c
done

for (( c=3000; c<=10000; c+=1000 ))
do
   echo "-- dim = $c"
   ./exr3.x $c
   #  Be careful when benchmarking higher values...
   #   at dim=3000 takes around 60s per function with O0 flag.
done

# echo time when all is done
echo "started at $starttime"
echo "finished at $(date)"

# housekeeping
rm exr3.x
rm exr3.o
rm matmul_loops.o
rm matmul_loops.mod


done # end of opt level iteration

exit 0
