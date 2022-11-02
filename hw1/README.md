# Homework 1

## EXR 1

> **Target file**: `exr1.f90`

In this program, I just create a vector `y` containing a discrete sine function sampled at step dx. Then I compute the derivative of $y$ with **4th order central difference**.

Compile the code with
```bash
gfortran exr1.f90 -o exr1.x
```
<br><br>

## EXR 2

> **Target files**: `exr2a.f90`, `exr2b.f90`

In this exercises, we test two (trivial) problems that occur when handling large value variables: overflow (for int) & infinities (for real & double).

Compile the code with
```bash
gfortran exr2a.f90 -o exr2a.x -fno-range-check
gfortran exr2b.f90 -o exr2b.x -fno-range-check
```

Notice that without the `no-range-check` flag **the compiler warns about the overflow**.  That's exactly what we are trying to test!


### Output

```
baronefr@cayde:~/uni/courses/quantum-computing/homework/hw1$ ./exr2a.x
 twomillion increment test -------
  sum with int*2 = -31615
  sum with int*4 =     2000001
 [i] direct overflow test ------
  int*2 :  32767 + 1  !=  -32768
  int*4 :  2147483647 + 1  !=  -2147483648
```

```
baronefr@cayde:~/uni/courses/quantum-computing/homework/hw1$ ./exr2b.x
 [ real ] multiplying   3.14159259E+32   1.41421360E+21
      ->  ans:          Infinity
 [double] multiplying   3.1415926535897933E+032   1.4142135381698609E+021
      ->  ans:    4.4428828621216639E+053
 [i] about real precision
  single precision max value   3.40282347E+38
  double precision max value   1.7976931348623157E+308
```

### What is happening?

The signed 2-complement representation of 2000000 is `00000000 00011110 10000100 10000000`.

However, if you try to store it as a `integer*2` only the least significant bits will be memorized, i.e. `10000100 10000000`, which is the 2 complement representation of -31616. Therefore, the sum to one will return -31615, as printed in the output. We need `integer*4` to store the correct representation of the number, avoiding this error.

<br>

Finally, the exr 2b shows just a precision loss in the least significant digits.

<br><br>

## EXR 3

> **Target files**: `exr3.f90`, `matmul_loops.f90`, `exr3_benchmarker.sh`

Compile with
```bash
gfortran -c matmul_loops.f90 exr3.f90
gfortran exr3.o matmul_loops.o -o exr3.x
```

*Remark*: You must provide as first argument the dimension $N$ of the (square) matrix to benchmark. Example:
```bash
./exr3.x 1000
```
Execute the benchmarks with `exr3_benchmarker.sh`. This script will store all the timing results to `data/benchmarks.txt`.

<br><br>

## Bibliography

- https://web.stanford.edu/class/me200c/tutorial_77/10_arrays.html
- https://www.adt.unipd.it/corsi/Bianco/www.pcc.qub.ac.uk/tec/courses/f90/stu-notes/F90_notesMIF_5.html
