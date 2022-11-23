! ========================================
!  QIC  >>  assignment 3 / exr 3
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 10 November 2022
! ========================================


program main
    use matrix_c16
    use mytools
    implicit none

    type(complex16_matrix) :: A
    double precision, dimension(:), allocatable :: eigv
    double precision :: avg
    integer :: iter

    ! size of matrix
    integer :: N

    ! etc
    integer :: icount


    ! ---------------------------------------
    ! parse argument from cmd
    character(len=512) :: cmdarg
    character(len=256) :: outfile
    character(len=16)  :: cmtype

    !  ARG 1:  dimension of matrix [int]
    call get_command_argument(1, cmdarg)
    read(cmdarg,'(I32)') N
    !  ARG 2:  dimension of matrix [int]
    call get_command_argument(2, cmdarg)
    read(cmdarg,'(I32)') iter
    !  ARG 3:  type of matrix to generate [string]
    call getarg(3, cmtype)
    !  ARG 4:  output file (optional) [string]
    call getarg(4, outfile)
    ! ---------------------------------------


    ! set debug verbosity -------------------
    DEBUG_LEVEL = DL_ERROR

    ! set manually a random seed ------------
    !   remark:  each element must be between 0 and 4095,
    !            and ISEED(4) must be odd
    call iseed_manual(iseed, 22, 7, 1999, 117)
    !   remark:  I made a function in mytools to do this :)



    ! ---------------------------------------
    allocate( eigv(N) )

    do icount = 1, iter

        if (trim(cmtype) == 'h') then
            ! init hermitian matrix
            A = InitHermitian(N, 's')
        else if (trim(cmtype) == 'd') then
            ! init diagonal matrix
            A = InitDiagonalReal(N, 's')
        else
            stop 'wrong arg flag (cmtype)'
        end if

        !print *, "preview of matrix:"
        !call CMatPrintHead(A, 3, 3)
        
        ! compute the eigenvalues
        eigv = eigenvalues(A)

        ! compute spacing // done in-place, discard first value
        eigv(2:N) = eigv(2:N) - eigv(1:(N-1))
        eigv(1) = 0

        ! normalize spacing
        avg = sum(eigv)/(N-1)
        print *, "average spacing = ", avg
        eigv = eigv/avg


        call write_array_txt(outfile, eigv, 2, N, .not.icount==1)

        call Delete(A)
    end do

    
end program ! -------------------------------