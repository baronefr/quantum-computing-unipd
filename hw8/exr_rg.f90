! ========================================
!  QIC  >>  assignment 8
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 16 December 2022
! ========================================



program main
    use matrix_r8
    use ising
    use rg_algo

    implicit none

    type(rga) :: rgobj
    type(real8_matrix) ::  sigmax
    real*8 :: lambda_max, lambda
    integer :: ii

    character(len=256) :: filename

    ! parameters of simulation
    integer*8 :: N = 6_8
    integer :: lambda_steps = 30
    

    lambda_max = 3d0


    ! other requirements
    sigmax = InitZero( (/2_8, 2_8/) )
    sigmax%data(2,1) = 1d0;       sigmax%data(1,2) = 1d0;


    ! system info
    DEBUG_LEVEL = DL_WARN

    
    ! open output file
    write(filename, "(AI0AI0A)") 'data/rg', N, '-', lambda_steps,  ".hw8"
    print *, 'output file : ', filename


    do ii = 0, lambda_steps

        lambda = ii*(lambda_max/lambda_steps)
        print *, '[!] lambda =', lambda
    
        rgobj = rg_init(N)

        if( ii.eq.0) then
            call write_header(rgobj, filename)
        end if
    
        ! init system of size N
        rgobj%H = ising_hamiltonian(N, lambda)
        !     H_int = A .tensor. B
        rgobj%A = kron_l1( 2_8**(N-1), sigmax )
        rgobj%B = kron_r1( sigmax, 2_8**(N-1) )

        call rg_exe(rgobj)

        call write_entry(rgobj, lambda, filename)

    end do

end program
