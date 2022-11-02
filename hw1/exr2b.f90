! ========================================
!  QIC  >>  assignment 1 / exr 2B
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 21 October 2022
! ========================================

program real_precision

    implicit none

    integer :: i, j
    real             :: x1_single, x2_single, rs
    double precision :: x1_double, x2_double, rd   ! aka real*8

    !!  Remark |  As Fortran does not have a builtin Pi constant,
    !  I will use the trick suggested
    !  here         https://community.intel.com/t5/Intel-Fortran-Compiler/Pi-value-in-double-precision/m-p/987489
    !  and here     https://stackoverflow.com/questions/2162461/fortran77-compiler-treatment-of-pi-4-d0datan1-d0
    !  which exploits the expression pi = 4*arctan(1) .
    !  It is up to the compiler to evaluate pi at the right precision,
    !  given the expression that we provide in the variable declaration.


    ! real  -------------------------------------------
    x1_single = 4.0*datan(1.0d0)*(10d0**32)
    x2_single = sqrt(real(2d0))*(10d0**21)
    rs = x1_single*x2_single

    print *, "[ real ] multiplying", x1_single, x2_single
    print *, "     ->  ans: ", rs

    ! double ------------------------------------------
    x1_double = 4.0d0*datan(1.0d0)*(10d0**32)
    x2_double = sqrt(real(2d0))*(10d0**21)
    rd = x1_double*x2_double

    print *, "[double] multiplying", x1_double, x2_double
    print *, "     ->  ans: ", rd

    ! useful things to know ---------------------------
    print *, "[i] about real precision"
    print *, " single precision max value", huge(x1_single)
    print *, " double precision max value", huge(x1_double)

end program real_precision
