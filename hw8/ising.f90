! ========================================
!  QIC  >>  assignment 8
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 17 December 2022
! ========================================


module ising

    use matrix_r8
    implicit none

contains


    function ising_hamiltonian(N, lambda) result(H)
        type(real8_matrix) :: H
        real*8 :: lambda

        integer*8 :: N, ii
        type(real8_matrix) :: tmp1
        type(real8_matrix) :: sigmax, sigmaz, ssig_z1, ssig_xx


        ! create some useful building blocks
        sigmax = InitZero( (/2_8, 2_8/) )
        sigmax%data(2,1) = 1d0;       sigmax%data(1,2) = 1d0;
        sigmaz = InitZero( (/2_8, 2_8/) )
        sigmaz%data(1,1) = lambda;    sigmaz%data(2,2) = -lambda;

        ssig_xx = kron( sigmax, sigmax )
        ssig_z1 = kron_r1( sigmaz, 2_8 )
        !ssig_1z = kron_l1( 2_8, sigmaz )

        H = InitZero( (/ 2**N, 2**N /) )
        !print *, 'allocated GB :', sizeof(H%data)/(1024**3)

        ! BAK ----
        !do ii = 1, N
            !tmp1 = kron_l1( 2_8**(ii-1), sigmaz )
            !call kron_r1_inplace( tmp1, 2_8**(N-ii), H )
            !call Delete(tmp1)
        !end do

        do ii = 1, N - 1
            tmp1 = kron_l1( 2_8**(ii-1), ssig_z1 )
            call kron_l1_inplace( 2_8**(ii-1), ssig_xx, tmp1)
            call kron_r1_inplace( tmp1, 2_8**(N-ii-1), H )
            call Delete(tmp1)
        end do
        call kron_l1_inplace( 2_8**(N-1), sigmaz, H)

        call Delete(sigmax)
        call Delete(sigmaz)
        call Delete(ssig_xx)
        call Delete(ssig_z1)
    end function



end module