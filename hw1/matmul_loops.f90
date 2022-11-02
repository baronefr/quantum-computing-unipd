! ========================================
!  QIC  >>  assignment 1 / exr 3
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 21 October 2022
! ========================================

module matmul_loops
    implicit none
contains

subroutine matmul_loop_ikj(m1,m2,m3) ! the one that should be the worse...
    implicit none
    real, dimension(:,:) :: m1, m2, m3
    integer ::  I,  K, J
    integer :: ii, kk, jj

    I = size(m1, 1);  K = size(m1, 2);
    J = size(m2, 2);

    do ii=1,I
        do kk=1,K 
            do jj=1,J 
                m3(ii,jj)=m3(ii,jj)+m1(ii,kk)*m2(kk,jj)
            end do
        end do
    end do

end subroutine

subroutine matmul_loop_ijk(m1,m2,m3) ! this is by row
    implicit none
    real, dimension(:,:) :: m1, m2, m3
    integer ::  I,  K, J
    integer :: ii, kk, jj

    I = size(m1, 1);  K = size(m1, 2);
    J = size(m2, 2);

    do ii=1,I
        do jj=1,J 
            do kk=1,K 
                m3(ii,jj)=m3(ii,jj)+m1(ii,kk)*m2(kk,jj)
            end do
        end do
    end do

end subroutine

subroutine matmul_loop_jik(m1,m2,m3) ! this is by columns
    implicit none
    real, dimension(:,:) :: m1, m2, m3
    integer ::  I,  K, J
    integer :: ii, kk, jj

    I = size(m1, 1);  K = size(m1, 2);
    J = size(m2, 2);

    do jj=1,J 
        do ii=1,I
            do kk=1,K 
                m3(ii,jj)=m3(ii,jj)+m1(ii,kk)*m2(kk,jj)
            end do
        end do
    end do

end subroutine


subroutine matmul_loop_jki(m1,m2,m3) ! this is optimizable by the compiler
    implicit none
    real, dimension(:,:) :: m1, m2, m3
    integer ::  I,  K, J
    integer :: ii, kk, jj

    I = size(m1, 1);  K = size(m1, 2);
    J = size(m2, 2);

    
    do jj=1,J 
        do kk=1,K
            do ii=1,I
                m3(ii,jj)=m3(ii,jj)+m1(ii,kk)*m2(kk,jj)
            end do
        end do
    end do

end subroutine

end module matmul_loops
