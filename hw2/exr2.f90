! ========================================
!  QIC  >>  assignment 2 / exr 2
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 28 October 2022
! ========================================


!  This program tests the matrix product with loops.
! In this version, the matmul_loops has a DOCUMENTATION.
! Eventually, the module matmul_loops uses explicitly 
! the CHECKPOINTS, as requested.

program main
    use mydebugger
    use matmul_loops
    implicit none
    
    ! matrix definitions
    real, dimension(:,:), allocatable :: m1, m2, m3
    integer :: dim   ! for matrix size

    dim = 700
    call checkpoint(msg='testing matrices...', var = dim)

    ! -------------------------------------------

    ! allocate the two matrices
    !   ... I am making the shapes different this time
    allocate(m1(dim,dim-5))
    allocate(m2(dim-5,dim-10))
    allocate(m3(dim,dim-10))

    ! init randomly the matrices
    call random_number(m1)
    call random_number(m2)


    ! HIGH VERBOSITY tests ----------------------
    DEBUG_LEVEL = 6
    call checkpoint(msg = ' testing IKJ')
    call matmul_loop_ikj(m1,m2,m3)
    
    call checkpoint(msg = ' testing IJK')
    call matmul_loop_ijk(m1,m2,m3)

    call checkpoint(level=6, &
                    msg = "You see that there are many messages... What about errors?")


    ! LOW VERBOSITY tests -----------------------
    DEBUG_LEVEL = 2
    call checkpoint(msg = ' testing JIK (doing INTENTIONALLY an error)')
    call matmul_loop_jik(m2,m1,m3)  ! Making an intentional error here!

    ! exit
    deallocate(m1)
    deallocate(m2)
    deallocate(m3)

    call checkpoint(msg = 'bye bye')

end program main
