! ========================================
!  QIC  >>  assignment 6
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 9 December 2022
! ========================================



program main
    use density_matrix_module

    type(density_matrix) :: dm
    type(density_matrix) :: rho
    integer :: system = 1

    character(len=512) :: cmdfname
    character(len=512) :: cmdarg
    character(len=512) :: outfile

    character(*), parameter :: input_prefix = 'data/test/'
    character(*), parameter :: output_prefix = 'data/answ/'



    !  ARG 1: density matrix filename
    print *, '[i] cmd args ----------------'
    call getarg(1, cmdfname)
    print *, "  fname = ", trim(cmdfname)

    !  ARG 2:  system to be traced on
    call get_command_argument(2, cmdarg)
    read(cmdarg, *) system
    print *, "  system =", system

    ! --------------------------------------



    ! read test density matrix from file
    dm = dmread( trim(adjustl(input_prefix))//trim(adjustl(cmdfname)) )
    if( dm%size.le.9) then
        ! print this matrix only if it is small ...
        call dmprint(dm, dm%size, dm%size)
    end if


    ! compute the partial trace on the selected system
    rho = partial_trace(dm, system)
    call dmprint(rho, rho%size, rho%size)

    write(outfile, "(AI0AA)") trim(adjustl(output_prefix)), system, '+', trim(adjustl(cmdfname))
    print *, 'writing to', outfile

    call dmwrite(rho, outfile)

end program