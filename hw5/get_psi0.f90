! ========================================
!  QIC  >>  assignment 5
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 28 November 2022
! ========================================


! This file is meant as an attachment to HW5.
! In this program, we will compute the first 10 eigenstates of
! the time-independent harmonic oscillator, then save them to a
! binary file, which will be imported in exr.f90 .


program main
    use matrix_r8
    use mytools
    use harmonic
    implicit none

    type(real8_matrix) :: Hp
    double precision, dimension(:), allocatable :: eigval
    type(real8_matrix) :: eigvec
    
    integer :: N, max_k
    double precision :: omega = 1d0, xmax = 10d0
    character(len=256) :: outfile
    character(len=128) :: filename
    character(len=512) :: cmdarg

    parameter( max_k = 10 )

    ! set global debug verbosity
    DEBUG_LEVEL = DL_WARN

    ! -------------------------------------------
    print *, "script called with args:"

    !  ARG 1:  dimension of matrix (int)
    call get_command_argument(1, cmdarg)
    read(cmdarg,'(I32)') N
    print *, "   N =", N

    !  ARG 2:  output file prefix [string]
    call getarg(2, filename)
    print *, "   filename = ", filename

    ! -------------------------------------------

    Hp = harmonic_tridiagonal_packed(N, omega, xmax)
    call RMatPrintHead(Hp, Hp%size(1), 3)

    call eigenvv_first_banded_packed(max_k, Hp, eigval, eigvec)
    ! remark: H is destroyed, but eigenv* are placed in [eigval, eigvec]
    
    ! write eigenvectors to file
    write(outfile, "(AI0AF0.1AF0.1A)") trim(adjustl(filename)), N, '-', xmax, '-', omega, ".eigvec"
    print *, "target filename: ", outfile
    call RMatDump(eigvec, outfile) ! this is a binary file

    ! -------------------------------------------
    
    print *, 'preview of eigvector:'
    call RMatPrintHead(eigvec, 6, 1)

    ! -------------------------------------------

    call Delete(Hp)
end program
