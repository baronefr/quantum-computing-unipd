! ========================================
!  QIC  >>  assignment 2 / exr 2
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 28 October 2022
! ========================================






module matmul_loops
!  DOC =============================================================
! 
!    This module implements matrix products with nested loops.
!  -----------------------------------------------------------------
!  TYPES:  none
! 
!  FUNCTIONS: none
! 
!  SUBROUTINES:
!    -  [family] matmul_loop_*    [all of them with same I/O conditions]
! 
!                 Input   |   m1, m2      input matrices to be multiplied
!                         |   m3          output matrix with the result of m1*m2
! 
!           Requirements  |   m1, m2, m3 must be allocated before call
! 
!                 Output  |   none        (implicit result returns to arg m3)
! 
!        Post-conditions  |   m3 will be overwritten
!  
!           Err handlers  |   I will check the correct matrix shapes and use the
!                         |   debugger to notice this error. The product is not
!                         |   executed in this case.
! 
!       +  list of available routines in this family:
!          - matmul_loop_ikj  ->  loop order is IKJ (worse loop ever)
!          - matmul_loop_ijk  ->  loop order is IJK (Rows by Columns)
!          - matmul_loop_jik  ->  loop order is JIK (Columns by Rows)
!          - matmul_loop_jki  ->  loop order is JKI (optimizable loop)
! 
! ==================================================================
    use mydebugger  ! use my checkpoints
    implicit none
contains

! ------------------------------------------------------------------
!   SUBROUTINE of matmul_loop family:  use IKJ loop order
! ------------------------------------------------------------------
subroutine matmul_loop_ikj(m1,m2,m3) ! the slowest one
    implicit none
    real, dimension(:,:) :: m1, m2, m3
    integer ::  I,  K, J
    integer :: ii, kk, jj

    I = size(m1, 1);  K = size(m1, 2);  J = size(m2, 2);
    
    ! check everything
    if ( K .ne. size(m2,1)) then ! 1) check m1 and m2 sizes 
        call checkpoint(level=DL_FATAL, msg = 'wrong input shapes for matrix product')
        stop '(fatal error)'
    end if
    if ( I .ne. size(m3,1)) then ! 2) check dim 1 of target matrix
        call checkpoint(level=DL_FATAL, msg = 'wrong target shape (1) for matrix product')
        stop '(fatal error)'
    end if
    if ( J .ne. size(m3,2)) then ! 3) check dim 2 of target matrix
        call checkpoint(level=DL_FATAL, msg = 'wrong target shape (2) for matrix product')
        stop '(fatal error)'
    end if
    
    m3=0

    call checkpoint(level=DL_PEDANTIC, msg = 'start matrix product')
    do ii=1,I
        do kk=1,K 
            do jj=1,J 
                m3(ii,jj)=m3(ii,jj)+m1(ii,kk)*m2(kk,jj)
            end do
        end do
    end do
    call checkpoint(level=DL_INFOP, msg = 'matrix product completed')
end subroutine

! ------------------------------------------------------------------
!   SUBROUTINE of matmul_loop family:  use IJK loop order
! ------------------------------------------------------------------
subroutine matmul_loop_ijk(m1,m2,m3) ! this is by row
    implicit none
    real, dimension(:,:) :: m1, m2, m3
    integer ::  I,  K, J
    integer :: ii, kk, jj

    I = size(m1, 1);  K = size(m1, 2);  J = size(m2, 2);
    
    ! check everything
    if ( K .ne. size(m2,1)) then ! 1) check m1 and m2 sizes 
        call checkpoint(level=DL_FATAL, msg = 'wrong input shapes for matrix product')
        stop '(fatal error)'
    end if
    if ( I .ne. size(m3,1)) then ! 2) check dim 1 of target matrix
        call checkpoint(level=DL_FATAL, msg = 'wrong target shape (1) for matrix product')
        stop '(fatal error)'
    end if
    if ( J .ne. size(m3,2)) then ! 3) check dim 2 of target matrix
        call checkpoint(level=DL_FATAL, msg = 'wrong target shape (2) for matrix product')
        stop '(fatal error)'
    end if

    m3=0

    call checkpoint(level=DL_PEDANTIC, msg = 'start matrix product')
    do ii=1,I
        do jj=1,J 
            do kk=1,K 
                m3(ii,jj)=m3(ii,jj)+m1(ii,kk)*m2(kk,jj)
            end do
        end do
    end do
    call checkpoint(level=DL_INFOP, msg = 'matrix product completed')
end subroutine


! ------------------------------------------------------------------
!   SUBROUTINE of matmul_loop family:  use JIK loop order
! ------------------------------------------------------------------
subroutine matmul_loop_jik(m1,m2,m3) ! this is by columns
    implicit none
    real, dimension(:,:) :: m1, m2, m3
    integer ::  I,  K, J
    integer :: ii, kk, jj

    I = size(m1, 1);  K = size(m1, 2);  J = size(m2, 2);

    ! check everything
    if ( K .ne. size(m2,1)) then ! 1) check m1 and m2 sizes 
        call checkpoint(level=DL_FATAL, msg = 'wrong input shapes for matrix product')
        stop '(fatal error)'
    end if
    if ( I .ne. size(m3,1)) then ! 2) check dim 1 of target matrix
        call checkpoint(level=DL_FATAL, msg = 'wrong target shape (1) for matrix product')
        stop '(fatal error)'
    end if
    if ( J .ne. size(m3,2)) then ! 3) check dim 2 of target matrix
        call checkpoint(level=DL_FATAL, msg = 'wrong target shape (2) for matrix product')
        stop '(fatal error)'
    end if

    m3=0

    call checkpoint(level=DL_PEDANTIC, msg = 'start matrix product')
    do jj=1,J 
        do ii=1,I
            do kk=1,K 
                m3(ii,jj)=m3(ii,jj)+m1(ii,kk)*m2(kk,jj)
            end do
        end do
    end do
    call checkpoint(level=DL_INFOP, msg = 'matrix product completed')
end subroutine


! ------------------------------------------------------------------
!   SUBROUTINE of matmul_loop family:  use JKI loop order
! ------------------------------------------------------------------
subroutine matmul_loop_jki(m1,m2,m3) ! this is optimizable by compiler
    implicit none
    real, dimension(:,:) :: m1, m2, m3
    integer ::  I,  K, J
    integer :: ii, kk, jj

    I = size(m1, 1);  K = size(m1, 2);  J = size(m2, 2);
    
    ! check everything
    if ( K .ne. size(m2,1)) then ! 1) check m1 and m2 sizes 
        call checkpoint(level=DL_FATAL, msg = 'wrong input shapes for matrix product')
        stop '(fatal error)'
    end if
    if ( I .ne. size(m3,1)) then ! 2) check dim 1 of target matrix
        call checkpoint(level=DL_FATAL, msg = 'wrong target shape (1) for matrix product')
        stop '(fatal error)'
    end if
    if ( J .ne. size(m3,2)) then ! 3) check dim 2 of target matrix
        call checkpoint(level=DL_FATAL, msg = 'wrong target shape (2) for matrix product')
        stop '(fatal error)'
    end if
    
    m3=0
    
    call checkpoint(level=DL_PEDANTIC, msg = 'start matrix product')
    do jj=1,J 
        do kk=1,K
            do ii=1,I
                m3(ii,jj)=m3(ii,jj)+m1(ii,kk)*m2(kk,jj)
            end do
        end do
    end do
    call checkpoint(level=DL_INFOP, msg = 'matrix product completed')
end subroutine

end module matmul_loops ! ==========================================
