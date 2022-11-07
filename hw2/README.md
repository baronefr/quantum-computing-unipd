# Homework 2

## EXR 1

> **Target files**: `exr1.f90`, `mydebugger.f90`

Compile with
```bash
gfortran -c mydebugger.f90 exr1.f90
gfortran exr1.o mydebugger.o -o exr1.x
```

Output:
```
DEBUG : T
 LEVEL :           0
 -- list of verbosity levels --
 DL_CHECK    :           0
 DL_FATAL    :           1
 DL_ERROR    :           2
 DL_WARN     :           3
 DL_INFO     :           4
 DL_INFOP    :           5
 DL_PEDANTIC :           6
 -----------
[checkpoint] your message
[ -!- ]    (real)    1.00000000    
[checkpoint] now we will test higher debug levels
 -> [pause] press enter to continue

 -> resume
 
[ error ] print me (2) error message
[warning] print me (3) warning
[ info  ] print me (4) info
[ info  ] print me (5) info but green
[  dbg  ] print me (6) pedantic
   (int4)            6
[checkpoint] last test, have a good day
```

<br><br>

## EXR 2

> **Target files**: `exr2.f90`, `matmul_loops.f90`, `mydebugger.f90`

Compile with
```bash
gfortran -c mydebugger.f90 matmul_loops.f90 exr2.f90
gfortran exr2.o matmul_loops.o mydebugger.o -o exr2.x
```

Output:
```
[checkpoint] testing matrices...
   (int4)          700
[checkpoint]  testing IKJ
[  dbg  ] start matrix product
[ info  ] matrix product completed
[checkpoint]  testing IJK
[  dbg  ] start matrix product
[ info  ] matrix product completed
[  dbg  ] You see that there are many messages... What about errors?
[checkpoint]  testing JIK (doing INTENTIONALLY an error)
[ fatal ] wrong input shapes for matrix product
STOP (fatal error)
```

<br><br>

## EXR 3

> **Target files**: `exr3.f90`, `type_complex_matrix`, `mydebugger.f90`

Compile with
```bash
gfortran -c mydebugger.f90 mod_matrix_c8.f90 exr3.f90
gfortran mydebugger.o mod_matrix_c8.o exr3.f90 -o exr3.x
```

Output (cropped on right margin):
```
[  dbg  ] allocating memory
[  dbg  ] allocating memory
 Trace of B is             (51.6117783,50.2991295)
       (7.629394531E-06,0.131537676)           (0.755605221,0.458650112)           (0.532767057,0.218959093)
           (0.688980818,0.702206612)           (0.987145424,0.954414845)           (0.851269484,0.289316177)
           (0.385769367,0.626861334)           (0.659053326,0.709319830)           (0.538661242,0.280288696)
[  dbg  ] allocating memory
[ info  ] writing CMat to data/matrix.txt
[ info  ] file written data/matrix.txt
 The matrix written to file
      (7.629394531E-06,-0.131537676)          (0.688980818,-0.702206612)          (0.385769367,-0.626861334)
          (0.755605221,-0.458650112)          (0.987145424,-0.954414845)          (0.659053326,-0.709319830)
          (0.532767057,-0.218959093)          (0.851269484,-0.289316177)          (0.538661242,-0.280288696)
[ info  ] reading CMat from data/matrix.txt
[ info  ] parsed matrix of size ->  100  100
[  dbg  ] allocating memory
[ info  ] reading values
[ info  ] CMat reading completed
 The matrix taken from file
      (7.629394531E-06,-0.131537676)          (0.688980818,-0.702206612)          (0.385769367,-0.626861334)
          (0.755605221,-0.458650112)          (0.987145424,-0.954414845)          (0.659053326,-0.709319830)
          (0.532767057,-0.218959093)          (0.851269484,-0.289316177)          (0.538661242,-0.280288696)
 The matrices are the same ...
  ... but now I change an element
[ error ] matrices are NOT equal
[  dbg  ] deallocating memory
[  dbg  ] deallocating memory
```

<br><br>

## Bibliography

- Here I found some **known issues in FORTRAN 90**. [Take a look](https://www.cs.rpi.edu/~szymansk/OOF90/bugs.html), I've encountered one of them so far! 
