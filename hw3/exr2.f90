! ========================================
!  QIC  >>  assignment 3 / exr 2
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 09 November 2022
! ========================================

!  useful Refs: LAPACK documentation
! https://netlib.org/lapack/explore-html/d0/d50/group__complex16_o_t_h_e_rauxiliary_ga81089fc5e8f0586d35f60128d35bef39.html
! https://netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_gaf23fb5b3ae38072ef4890ba43d5cfea2.html


program main
    use matrix_c16
    use mytools
    implicit none

    type(complex16_matrix) :: A
    double precision, dimension(:), allocatable :: eigv
    double precision :: avg

    ! size of matrix
    integer :: N

    ! etc
    integer :: ii
    integer, parameter :: fileunit = 22
    logical :: file_isopened


    ! ---------------------------------------
    ! parse argument from cmd
    character(len=512) :: cmdarg
    character(len=256) :: outfile
    !  ARG 1:  dimension of matrix [int]
    call get_command_argument(1, cmdarg)
    read(cmdarg,'(I32)') N
    !  ARG 2:  output file (optional) [string]
    call getarg(2, outfile)
    ! ---------------------------------------


    ! set debug verbosity -------------------
    DEBUG_LEVEL = DL_PEDANTIC

    ! set manually a random seed ------------
    !   remark:  each element must be between 0 and 4095,
    !            and ISEED(4) must be odd
    call iseed_manual(iseed, 22, 7, 1999, 117)
    !   remark:  I made a function in mytools to do this :)

    ! ---------------------------------------
    allocate( eigv(N) )

    ! init hermitian matrix
    A = InitHermitian(N, 'n')

    ! Remark: for a diagonal matrix, use
    !  A = InitDiagonalReal(N, 'n') 

    ! just check the heaad of the matrix
    print *, "preview of matrix:"
    call CMatPrintHead(A, 3, 3)

    ! compute the eigenvalues
    eigv = eigenvalues(A)

    ! compute spacing // done in-place, discard first value
    eigv(2:N) = eigv(2:N) - eigv(1:(N-1))
    eigv(1) = 0

    ! normalize spacing
    avg = sum(eigv)/(N-1)
    print *, "average spacing = ", avg
    eigv = eigv/avg


    ! ---------------------------------------
    ! write spacing to file
    if (trim(outfile) == '') then
        print *, "skipping file operation, arg is empty"
    else
        print *, "appending to ", outfile
        open(fileunit, file=outfile, action='write')
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then
            call checkpoint(level=DL_INFO, msg = 'writing file '//outfile)
            do ii = 2,N ! recall to discard the first zero
                write(fileunit, *) eigv(ii)
            end do
        else
            call checkpoint(level=DL_ERROR, msg = 'file '//outfile//' cannot be opened')
        end if
    end if
    ! ---------------------------------------

    call Delete(A)
end program ! -------------------------------