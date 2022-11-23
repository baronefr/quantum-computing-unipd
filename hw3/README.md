# Homework 3



## EXR 1

> **Target files**: `exr1_launcher.f90`, `exr1.py`, `exr1_plotter.py`, `lib/matmul_loops.f90`

For standalone execution, compile with
```bash
make exr1
```
and execute with
```bash
./exr1 1000 1
```

To execute the automatic script launch (with python3)
```bash
./exr1.py
```

Finally, to make the plots (outputs in `img/`), execute (with python3)
```bash
./exr1_plotter.py
```




## EXR 2

> **Target files**: `exr2.f90`, `lib/matrix_c16.f90`

Compile with
```bash
make exr2
```

Execute:
```bash
./exr2.x 1000 data/exr2_herm.txt
```

Output:
```
baronefr@cayde:/mnt/beta/baronefr/uni/courses/quantum-computing/homework/hw3$ ./exr2.x 1000 data/exr2_herm.txt
[  dbg  ] handle dist selection
   (int4)            3
[  dbg  ] allocating memory
[ info  ] generate rand values with LAPACK
[ info  ] rand done
 preview of matrix:
             (-0.91879412032584240,0.0000000000000000)            (-1.3247792688980957,-0.29606247222751958)            (-0.54770114411521764,0.49141802627109871)
             (-1.3247792688980957,0.29606247222751958)              (0.98378091409430135,0.0000000000000000)             (1.0017683601099252,-0.41573792647671542)
           (-0.54770114411521764,-0.49141802627109871)              (1.0017683601099252,0.41573792647671542)              (0.90320784043070224,0.0000000000000000)
[ info  ] eigenvalues ok
 average spacing =   0.17811352298425667
 appending to data/exr2_herm.txt
[ info ] writing file data/exr2_herm.txt
[  dbg  ] deallocating memory
```




## EXR 3

> **Target files**: `exr3.f90`, `lib/matrix_c16.f90`

Compile with
```bash
make exr3
```

Execute:
```bash
./exr3.x 1000 50 h data/exr3_herm.txt
./exr3.x 1000 50 d data/exr3_diag.txt
```

Output:
```
baronefr@cayde:/mnt/beta/baronefr/uni/courses/quantum-computing/homework/hw3$ python3 exr3_plotter.py 
loading file data/exr3_diag.txt (diag)
 imported 1998000 values
loading file data/exr3_herm.txt (herm)
 imported 1998000 values
fit coeffs : [12.82202943  2.75436433  2.54867489  1.33552144]
  with err : [1.46859158 0.11571082 0.06397999 0.03113055]
/mnt/beta/baronefr/uni/courses/quantum-computing/homework/hw3/exr3_plotter.py:43: RuntimeWarning: overflow encountered in exp
  return a * (x**alpha) * np.exp(-b *(x**beta))
fit coeffs : [3.81700849 2.65254384 0.164551   0.68728475]
  with err : [0.22996449 0.06128966 0.01361065 0.01230137]
fit coeffs : [4.07067961 1.55269721]
  with err : [0.0260219  0.00671097]
fit coeffs : [1.79477786 1.93092409]
  with err : [0.00504232 0.00766586]
```