# compile the exr
make exr

# created the required directories
mkdir -p data/
mkdir -p data/test
mkdir -p data/answ

# run the script to generate the random states
python3 test_generate_states.py

# run the exr for different matrices and system values
./test.x 00.dat          1
./test.x 00.dat          2
./test.x bell0011.dat    1
./test.x bell0011.dat    2
./test.x rand5.dat       1
./test.x rand5.dat       3
./test.x rand5.dat       4
./test.x rand10.dat      1
./test.x rand10.dat      7
./test.x rand10.dat      10

# check the traced matrices with Qiskit
python3 test_check_answ.py

# clean the exe
make cleanx