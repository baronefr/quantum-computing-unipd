! ========================================
!  QIC  >>  assignment 4
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 19 November 2022
! ========================================



program main
    use matrix_r8
    use mytools
    use harmonic
    implicit none

    type(real8_matrix) :: Hp
    double precision, dimension(:), allocatable :: eigval
    type(real8_matrix) :: eigvec
    
    integer :: N, max_k
    double precision :: omega = 1d0, xmax = 20d0
    character(len=256) :: outfile
    character(len=128) :: outprefix
    character(len=512) :: cmdarg

    parameter( max_k = 200 )

    ! set global debug verbosity
    DEBUG_LEVEL = DL_WARN

    ! -------------------------------------------
    print *, "script called with args:"

    !  ARG 1:  dimension of matrix (int)
    call get_command_argument(1, cmdarg)
    read(cmdarg,'(I32)') N
    print *, "   N =", N

    !  ARG 2:  max value of x (real)
    call get_command_argument(2, cmdarg)
    read(cmdarg, *) xmax
    print *, "   xmax =", xmax

    !  ARG 3:  output file prefix [string]
    call getarg(3, outprefix)
    print *, "   outprefix = ", outprefix

    ! -------------------------------------------

    Hp = harmonic_pentadiagonal_packed(N, omega, xmax)
    call RMatPrintHead(Hp, Hp%size(1), 3)

    call eigenvv_first_banded_packed(max_k, Hp, eigval, eigvec)
    ! remark: H is destroyed, but eigenv* are placed in [eigval, eigvec]
    
    ! write eigenvectors to file
    write(outfile, "(AI0AF0.1AF0.1A)") trim(adjustl(outprefix)), N, '-', xmax, '-', omega, ".eigvec"
    print *, "target filename: ", outfile
    call RMatDump(eigvec, outfile) ! this is a binary file

    ! write eigenvalues to file
    write(outfile, "(AI0AF0.1AF0.1A)") trim(adjustl(outprefix)), N, '-', xmax, '-', omega, ".eigval"
    print *, "target filename: ", outfile
    call dump_real8_array_txt(outfile, eigval, 1, max_k, .false.) ! this is a csv file

    ! -------------------------------------------

    call Delete(Hp)
end program
