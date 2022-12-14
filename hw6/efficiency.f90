! ========================================
!  QIC  >>  assignment 6
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 9 December 2022
! ========================================



program main
    use quantum_states_module

    type(qstate) :: psi
    integer :: iid, jjn, kkr, rep

    integer :: dmax, nmax
    double precision :: time_tic, time_toc, time   ! for cpu timing

    ! files to use for benchmark storage
    character(len=256) :: fname_sep = 'data/separable.csv'
    character(len=256) :: fname_gen = 'data/general.csv'

    ! rng seed
    iseed(1) = 22
    iseed(2) = 22
    iseed(3) = 22
    iseed(4) = 25

    ! parameters for the benchmark  [ WARNING: 32GB ram required ]
    dmax = 14
    nmax = 8
    rep = 100 ! how many times to repeat the benchmark for separable state

    ! ------------------------------------------
    ! benchmark separable states init
    ! ------------------------------------------
    open  (22, file=fname_sep, action='write', position='append')

    ! benchmark separable states init (more times, cause it takes 0 time lol)
    do iid = 2, dmax
        do jjn = 2, nmax
            print *, 'allocate ', iid, '*', jjn

            call CPU_TIME(time_tic)
            do kkr = 1, rep

                call allocate_separable(iid, jjn, psi)

                ! free memory
                call qdeall(psi)
        
            end do
            call CPU_TIME(time_toc)
            time = time_toc-time_tic

            write (22,'(*(G0.8,:,","))') iid, jjn, time/rep

        end do
    end do

    close (22)

    

    ! ------------------------------------------
    ! benchmark general states init
    ! ------------------------------------------
    open  (22, file=fname_gen, action='write', position='append')

    do iid = 2, dmax
        do jjn = 2, nmax
            print *, 'allocate ', iid, '^', jjn

            call CPU_TIME(time_tic)
            call allocate_general(iid, jjn, psi)
            call CPU_TIME(time_toc)
            time = time_toc-time_tic

            ! free memory
            call qdeall(psi)

            allocate( psi%data(400) )
            psi%data(10) = 4
            deallocate( psi%data )

            ! write benchmarks to file
            write (22,'(*(G0.8,:,","))') iid, jjn, time
            
        end do
    end do

    close (22)
    

end program