! ========================================
!  QIC  >>  assignment 3 / exr 1
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 09 November 2022
! ========================================


!  to execute this script, provide the following cmd args:
!
!   -  dim [int]          dimension of the square matrices to be multiplied
!   -  rep [int]          how many times to repeat the same benchmark
!                           (takes the average for rep > 1)
!   -  outfile [string]   the name of the output file to store the benchmarks
!

program matrix_test
    use matmul_loops  ! link from other file
    implicit none
    
    ! matrix definitions
    real, dimension(:,:), allocatable :: m1, m2, m3
    integer :: dim   ! for matrix dimension

    ! for benchmarking
    integer :: rep, rep_idx                     ! to repeat N times the operations
    double precision :: time_tic, time_toc      ! for cpu timing
    double precision, dimension(4) :: dts = 0   ! to store time intervals
    
    ! to parse argument from cmd (useful with scripts)
    character(len=256) :: outfile
    character(len=256) :: cmdarg

    !  ARG 1:  dimension of matrix (int)
    call get_command_argument(1, cmdarg)
    read(cmdarg,'(I32)') dim

    !  ARG 2:  number of benchmarks to execute (int)
    call get_command_argument(2, cmdarg)
    read(cmdarg,'(I32)') rep

    !  ARG 3:  output file (string)
    call getarg(3,outfile)

    ! printing ARGS
    print *, "exe called with:"
    print *, " dim     : ", dim
    print *, " rep     : ", rep
    print *, " outfile : ", outfile

    ! -------------------------------------------

    ! allocate the matrices
    allocate(m1(dim,dim))
    allocate(m2(dim,dim))
    allocate(m3(dim,dim))

    ! init randomly
    call random_number(m1)
    call random_number(m2)
    ! m3 is reduced to 0 at every computation

    ! repeat the measure
    do rep_idx = 1, rep

        ! [test 1] intuitive loop -------------------------
        call CPU_TIME(time_tic)
            call matmul_loop_ijk(m1,m2,m3)
        call CPU_TIME(time_toc)
        dts(1) = dts(1) + (time_toc-time_tic)

        ! [test 2] intuitive swap - -----------------------
        call CPU_TIME(time_tic)
            call matmul_loop_jik(m1,m2,m3)
        call CPU_TIME(time_toc)
        dts(2) = dts(2) + (time_toc-time_tic)

        ! [test 3] more clever -----------------------
        call CPU_TIME(time_tic)
            call matmul_loop_jki(m1,m2,m3)
        call CPU_TIME(time_toc)
        dts(3) = dts(3) + (time_toc-time_tic)

        ! [test 4] matmul ----------------------------
        call CPU_TIME(time_tic)
            m3 = matmul(m1,m2)
        call CPU_TIME(time_toc)
        dts(4) = dts(4) + (time_toc-time_tic)
    end do

    deallocate(m1)
    deallocate(m2)
    deallocate(m3)

    ! take the average time
    dts = dts / rep

    ! print to screen (just to take a look while executing)
    print *, "!n,   ijk,   jik,   jki,   matmul"
    print '(*(G0.6,:,"  "))', dim, dts
    
    ! print to file, if arg is not empty (store data to be plotted)
    if (trim(outfile) == '') then
        print *, "skipping file operation, arg is empty"
    else
        print *, "appending to ", outfile
        open  (22, file=outfile, action='write', position='append')
        write (22,'(*(G0.8,:,","))') dim, dts   ! print as comma separated values
        close (22)
    end if
end program matrix_test