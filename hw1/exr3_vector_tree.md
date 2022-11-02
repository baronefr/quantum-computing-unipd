# EXR3 vector tree check

For general reference on the `dump-tree` options, check the [GCC](https://gcc.gnu.org/onlinedocs/gcc/Developer-Options.html) documentation.

For this test, I have removed all the loop matrix multiplication from exr3.f90, except for the **loop JKI**. The code will run only once this function.


## no compiler optimization

To check the optimized tree without any optmization, compile with
```
gfortran -fdump-tree-optimized exr3.f90 -O0
```

In the *.optimized* file you would see something like
```
  <bb 17> :
  _52 = (integer(kind=8)) jj_82;
  _53 = stride.15_126 * _52;
  _54 = offset.16_129 + _53;
  _55 = (integer(kind=8)) ii_81;
  _56 = stride.13_123 * _55;
  _57 = _54 + _56;
  _58 = (*m3.0_124)[_57];
  _59 = (integer(kind=8)) kk_83;
  _60 = stride.3_98 * _59;
  _61 = offset.4_101 + _60;
  _62 = (integer(kind=8)) ii_81;
  _63 = stride.1_95 * _62;
  _64 = _61 + _63;
  _65 = (*m1.0_96)[_64];
  _66 = (integer(kind=8)) jj_82;
  _67 = stride.9_112 * _66;
  _68 = offset.10_115 + _67;
  _69 = (integer(kind=8)) kk_83;
  _70 = stride.7_109 * _69;
  _71 = _68 + _70;
  _72 = (*m2.0_110)[_71];
  _73 = _65 * _72;
  _74 = (integer(kind=8)) jj_82;
  _75 = stride.15_126 * _74;
  _76 = offset.16_129 + _75;
  _77 = (integer(kind=8)) ii_81;
  _78 = stride.13_123 * _77;
  _79 = _76 + _78;
  _80 = _58 + _73;
  (*m3.0_124)[_79] = _80;
  ii_203 = ii_81 + 1;
  goto <bb 16>; [INV]
```
which is not optimized at all.


## with O5 optimization

To dump the optimized tree, run
```
gfortran -fdump-tree-optimized exr3.f90 -O5
```

The relevant section of the *.optimized* tree is
```
  <bb 22> [local count: 44957436]:
  # ivtmp.176_401 = PHI <0(21), ivtmp.176_413(22)>
  vect__100.106_386 = MEM <vector(4) real(kind=4)> [(real(kind=4) *)vectp_m3.105_377 + ivtmp.176_401 * 1];
  vect__97.109_397 = MEM <vector(4) real(kind=4)> [(real(kind=4) *)vectp_m1.108_387 + ivtmp.176_401 * 1];
  vect__96.110_399 = vect__97.109_397 * vect_cst__398;
  vect__95.111_400 = vect__100.106_386 + vect__96.110_399;
  vect_cst__411 = {m2.0__I_lsm0.149_197, m2.0__I_lsm0.149_197, m2.0__I_lsm0.149_197, m2.0__I_lsm0.149_197};
  vect__222.119_435 = MEM <vector(4) real(kind=4)> [(real(kind=4) *)_567 + ivtmp.176_401 * 1];
  vect__223.120_436 = vect_cst__411 * vect__222.119_435;
  vect__224.121_437 = vect__95.111_400 + vect__223.120_436;
  MEM <vector(4) real(kind=4)> [(real(kind=4) *)vectp_m3.105_377 + ivtmp.176_401 * 1] = vect__224.121_437;
  ivtmp.176_413 = ivtmp.176_401 + 16;
  if (ivtmp.176_413 == _417)
    goto <bb 23>; [16.67%]
  else
    goto <bb 22>; [83.33%]
```
which is much more vector... :)
