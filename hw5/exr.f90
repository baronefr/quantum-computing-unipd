! ========================================
!  QIC  >>  assignment 5
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 28 November 2022
! ========================================



program main
    
    use trotter_suzuki
    use matrix_r8

    type(trotter_suzuki_sim) :: tss
    type(real8_matrix) :: psi_zero
    integer :: Nt
    real*8 :: xmax, tmax, tau
    integer, parameter :: chosen_eigenstate = 1 ! use the first eigenstate stored in init dataset

    integer :: ii
    character(len=512) :: outfile
    character(len=512) :: cmdarg

    integer, parameter :: fileunit = 22

    DEBUG_LEVEL = DL_WARN

    ! -------------------------------------------
    print *, "exe called with args : "
    !  ARG 1:  number of time steps (int)
    call get_command_argument(1, cmdarg)
    read(cmdarg,'(I32)') Nt
    print *, "   Nt =", Nt

    !  ARG 2:  tmax (real)
    call get_command_argument(2, cmdarg)
    read(cmdarg, *) tmax
    print *, "   tmax =", tmax

    !  ARG 3:  tau (real)
    call get_command_argument(3, cmdarg)
    read(cmdarg, *) tau
    print *, "   tau =", tau

    !  ARG 4:  output file prefix [string]
    call getarg(4, outfile)
    print *, "   prefix = ", outfile

    ! -------------------------------------------

    ! read psi_0 from file
    !psi_zero = RMatRead('data/psi0_omega3_10001-10.0-3.0.eigvec')
    psi_zero = RMatRead('data/psi0_10001-10.0-1.0.eigvec')
    xmax = 10d0 ! depends on chosen dataset

    ! init tss object
    tss = init_trotter_suzuki(xmax, tmax, tau, &
                              psi_zero%size(1), Nt, &
                              psi_zero%val(:,chosen_eigenstate) )
    call normalize(tss%psi) ! normalize psi, to be sure...

    ! print the details of this simulation
    call tellme_ts(tss)


    ! open file to dump tss details

    ! LOG FILE STRUCTURE ---------------------
    ! 
    !  (integer) Nx
    !  (integer) Nt
    !  (real)    xmax
    !  (real)    tmax
    !  (real)    omega
    !  (real)    tau
    !  {  (real) time    (double complex)  psi  }  x Nt
    ! ---------------------------------------
    write(cmdarg, "(AI0AI0AF0.1A)") trim(adjustl(outfile)), psi_zero%size(1), '-', Nt, '-', tau, ".dat"
    print *, ""
    print *, "target file : ", cmdarg
    open(fileunit, file = cmdarg, action='write', form="unformatted")



    ! write metadata
    write(fileunit) tss%Nx
    write(fileunit) tss%Nt
    write(fileunit) tss%xmax
    write(fileunit) tss%tmax
    write(fileunit) tss%omega
    write(fileunit) tss%tau

    ! write the initial state
    write(fileunit) tss%time
    write(fileunit) tss%psi




    print *, ' +++ computing +++'

    ! compute (& save) the time evolution of psi
    do ii = 1, Nt
        call time_evolution(tss)
        call normalize(tss%psi)

        write(fileunit) tss%time
        write(fileunit) tss%psi
    end do

    close(fileunit)


    print *, ''
    print *, ' [info] final psi (head):', tss%psi(1:6)

    call destroy_ts(tss)
end program
