# Quantum Information and Computing

The weekly homeworks of QIC [class](https://en.didattica.unipd.it/off/2021/LM/SC/SC2443/000ZZ/SCP8082721/N0) (Prof. Simone Montangero).

<br>



### Topic 1: introduction to FORTRAN and LAPACK


#### **[HW 1](slides/hw1-Barone-SLIDES.pdf)** (due 1 nov)

Evaluation: ⭐⭐⭐⭐/4

1. **Warm-up Fortran program** | I choose to implement the 4th order central difference.
2. **Integer and real numbers** | Overflow issues of variable types.
3. **Matrix product CPU time benchmarking** | Benchmark the CPU time of matrix product functions: nested loops vs `matmul`.


#### **[HW 2](slides/hw2-Barone-SLIDES.pdf)** (due 8 nov)

Evaluation: ⭐⭐⭐⭐/4

1. **Checkpoint subroutine** | Write a checkpoint subroutine for debugging.
2. **Documentation** | Rewrite HW1/EXR3 with documentation.
3. **Complex matrix type** | Define a Complex Matrix type, with appropriate routines for initialization, deallocation and some math interfaces.


#### **[HW 3](slides/hw3-Barone-SLIDES.pdf)** (due 22 nov)

Evaluation: ⭐⭐⭐⭐/4

1. **Matrix product CPU time benchmarking** | Write a python script to automatize the benchmarking process (I already did this with a bash script in HW1).
2. **Diagonalize a complex matrix** | Use LAPACK routines to generate random hermitian matrices and diagonalize.
3. **Random matrix theory** | Study the eigenvalue spacing of random hermitian matrices (Wiger-Dyson statistics).

<br>


### Topic 2: solutions to Schrödinger equation


#### **[HW 4](slides/hw4-Barone-SLIDES.pdf)** (due 29 nov)

Evaluation: ⭐⭐⭐⭐/4

* Solution of the **time independent Schrödinger equation**. Finite difference methods have been used to compute the second order derivative.


#### **[HW 5](slides/hw5-Barone-SLIDES.pdf)** (due 6 dec)

Evaluation: ⭐⭐⭐⭐/4

* Solution of the **time dependent Schrödinger equation** with Trotter-Suzuki expansion.

<br>


### Topic 3: many-body systems


#### **[HW 6](slides/hw6-Barone-SLIDES.pdf)** (due 13 dec)

Evaluation: ⭐⭐⭐⭐/4

* Description of systems with **density matrices**: representation in FORTRAN and implementation of **partial trace** subroutine.


#### **[HW 7](slides/hw7-Barone-SLIDES.pdf)** (due 20 dec)

Evaluation: ⭐⭐⭐⭐/4

* **Transverse field Ising model**: exact diagonalization of $N$ spin-1/2 particles on a **one-dimensional lattice**.


#### **[HW 8](slides/hw8-Barone-SLIDES.pdf)** (due 6 jan)

Evaluation: ⭐⭐⭐⭐/4
<br>
My evaluation:⭐⭐⭐☆/4
> I could have done something better, to improve the convergence speed. Maybe I will play with this HW again later on...

* Transverse field Ising model with **Real Space Renormalization Group** and **Infinite Density Matrix Renormalization Group**.


***

<h5 align="center">Quantum Information and Computing<br>AY 2022/2023 - University of Padua</h5>

<p align="center">
  <img src="https://raw.githubusercontent.com/baronefr/baronefr/main/shared/2022_unipd.png" alt="" height="70"/>
  &emsp;
  <img src="https://raw.githubusercontent.com/baronefr/baronefr/main/shared/2022_pod.png" alt="" height="70"/>
</p>
