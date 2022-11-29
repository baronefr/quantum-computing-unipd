#!/bin/bash -e

# ========================================
#  QIC  >>  assignment 4
#  UniPD, AY 2022/23, Physics of Data
# ----------------------------------------
#   coder : Barone Francesco
#   dated : 22 November 2022
# ========================================

# This script iterates the execution of exr for various parameter N
# values. The outputs (eigenvalues & vectors) are stored in ./data/* .

# preprocess args ------------------
do_extra=0
while getopts e opt; do
   case $opt in
     e ) do_extra=1 ;;
  esac
done
shift $((OPTIND-1))


arg_xmax='20.0'
if [[ "${do_extra}" -eq 1 ]] ; then
   echo "doing extra!"
   sleep 1s

   make extra
   arg_data='data/penta/'
else
   echo "doing std method!"
   sleep 1s

   make exr
   arg_data='data/tri/'
fi

make clean

# echo current time, just to see how much it takes
starttime=$(date)
echo " [i] start at $starttime"

# call the program with different matrix sizes
for (( c=101; c<=1001; c+=100 )) 
do
   echo "-- dim = $c"
   ./exr.x $c $arg_xmax $arg_data
done

for (( c=2001; c<=10001; c+=1000 ))
do
   echo "-- dim = $c"
   ./exr.x $c $arg_xmax $arg_data
done

# echo time when all is done
echo "started at $starttime"
echo "finished at $(date)"

# housekeeping
make cleanx

exit 0
