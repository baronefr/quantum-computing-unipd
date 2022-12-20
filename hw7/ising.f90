! ========================================
!  QIC  >>  assignment 7
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 15 December 2022
! ========================================

module ising

    use matrix_r8_pack
    implicit none

    contains

    function generate_packed_hamiltonian(N, lambda) result(H)
        type(real8_pack_matrix) :: H
        real*8 :: lambda

        integer*8 :: N, ii
        type(real8_pack_matrix) :: tmp1, tmp2
        type(real8_pack_matrix) :: sigmax, sigmaz, ssig_z1, ssig_xx


        ! create some useful building blocks
        sigmax = InitZerop( 2_8 )
        sigmax%data(2) = 1d0
        sigmaz = InitZerop( 2_8 )
        sigmaz%data(1) = lambda;    sigmaz%data(3) = -lambda;

        ssig_xx = kronp( sigmax, sigmax )
        ssig_z1 = kronp_r1( sigmaz, 2_8 )

        H = InitZerop( 2**N )
        print *, ' -> allocated GB :', real(sizeof(H%data))/(1024d0**3)

        ! inplace algo
        do ii = 1, N - 1
            tmp1 = kronp_l1( 2_8**(ii-1), ssig_z1 )
            call kronp_l1_inplace( 2_8**(ii-1), ssig_xx, tmp1)
            call kronp_r1_inplace( tmp1, 2_8**(N-ii-1), H )
            call Deletep(tmp1)
        end do
        call kronp_l1_inplace( 2_8**(N-1), sigmaz, H)

        call Deletep(sigmax)
        call Deletep(sigmaz)
        call Deletep(ssig_xx)
        call Deletep(ssig_z1)
    end function


end module