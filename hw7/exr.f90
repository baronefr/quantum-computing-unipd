! ========================================
!  QIC  >>  assignment 7
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 14 December 2022
! ========================================



program main
    use matrix_r8_pack
    use ising

    implicit none

    type(real8_pack_matrix) :: H
    double precision, dimension(:), allocatable :: eigval
    double precision, dimension(:,:), allocatable :: eigvec
    real*8 :: lambda_max, lambda
    integer :: ii

    character(len=256) :: filename
    integer, parameter :: fileunit = 22
    logical :: file_isopened
    double precision :: time_tic, time_toc, time   ! for cpu timing

    ! parameters of simulation
    integer*8 :: N = 2_8
    integer :: lambda_steps = 30
    integer :: neig = 10
    lambda_max = 3d0


    ! system info
    DEBUG_LEVEL = DL_WARN




    
    ! open output file
    write(filename, "(AI0AI0AI0A)") 'data/', N, '-', lambda_steps, '-',  neig, ".hw7"
    print *, 'output file : ', filename

    open(fileunit, file=filename,  action='write', form='unformatted') ! status='new',
    inquire(unit=fileunit, opened=file_isopened)

    if(file_isopened) then
        call checkpoint(level=DL_INFO, msg = 'writing (binary)'//filename)

        write(fileunit) N 
        write(fileunit) lambda_steps + 1
        write(fileunit) neig

        close(fileunit)

    else
        call checkpoint(level=DL_ERROR, msg = 'file '//filename//' cannot be opened')
        call exit(1)
    end if



    ! compute
    print *, '+++ EXE SIMULATION (packed mode) +++'
    call CPU_TIME(time_tic)
    do ii = 0, lambda_steps

        lambda = ii*(lambda_max/lambda_steps)
        print *, '[!] lambda =', lambda

        H = generate_packed_hamiltonian(N, lambda)
        print *, '  hamiltonian ok, diagonalizing ...'
        !call PMatDump(H, 'ciccio.dat') ! for debug

        call eigenvv_pack(neig, H, eigval, eigvec, do_eigvec = .false.)
        ! remark: eigvec is not referenced

        ! dump to file
        open(fileunit, file=filename, status="old", action="write", &
                       form='unformatted', position='append')
        write(fileunit) eigval(1:neig)
        close(fileunit)

        deallocate( eigval )
    end do

    call CPU_TIME(time_toc)
    time = time_toc-time_tic
    print *, 'total CPU time (min)', time/60d0

    print *, 'data written in', filename
end program

