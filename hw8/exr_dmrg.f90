! ========================================
!  QIC  >>  assignment 8
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 24 December 2022
! ========================================




program main
    use matrix_r8
    use ising
    use idmrg_algo

    implicit none

    type(idmrga) :: dmobj
    real*8 :: lambda_max, lambda, thres, prev_energy
    integer :: ii

    character(len=256) :: filename
    character(len=256) :: cmdfname

    ! parameters of simulation
    integer :: max_m = 3
    integer :: lambda_steps = 30
    integer :: max_iter = 1000   ! max iterations per DMRG run

    ! convergence thres of DMRG
    thres = 1e-10

    ! lambda parameter range
    lambda_max = 3d0

    ! system info
    DEBUG_LEVEL = DL_INFO


    !  ARG 1: suffix filename
    print *, '[i] cmd args ----------------'
    call getarg(1, cmdfname)
    print *, "  label = ", trim(cmdfname)

    
    ! open output file
    write(filename, "(AI0AI0AA)") 'data/dmrg', max_m, '-', lambda_steps, trim(cmdfname), ".hw8"
    print *, 'output file : ', filename




    ! execute DMRG for all lambdas
    do ii = 0, lambda_steps

        lambda = ii*(lambda_max/lambda_steps)
        print *, '[!] lambda =', lambda
        
        
        ! init infinite dmrg
        dmobj = idmrg_init(lambda, max_m = max_m, max_iter = max_iter )


        if( ii.eq.0 ) then
            ! write to file the parameters
            call write_header_dmrg(dmobj, filename)
        end if
    

        prev_energy = huge(prev_energy)
        do while ( .not.dmobj%converged )
            
            ! execute idmrg step
            call idmrg_exe(dmobj)
            !print *, '  gs = ', dmobj%energy  ! print the energy @ every iteration

            if( abs(dmobj%energy - prev_energy) .lt. thres ) then
                print *, 'convergence reached at iteration', dmobj%counter
                dmobj%converged = .true.
                exit
            end if

            prev_energy = dmobj%energy
        end do

        print *, 'converged energy = ', dmobj%energy

        print *, 'write entry to file'
        call write_entry_dmrg(dmobj, filename)

        call idmrg_delete(dmobj)

    end do    ! repeat with next lambda



end program
