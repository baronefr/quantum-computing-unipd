! ========================================
!  QIC  >>  assignment 6
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 10 December 2022
! ========================================




module quantum_states_module
    implicit none

    type :: qstate
        integer :: dim, nn
        logical :: is_general
        complex(8), dimension(:), allocatable :: data
    end type

    ! global seed for LAPACK generator
    integer :: iseed(4) 

contains

    subroutine allocate_separable(dim, nn, psi)
        integer :: dim, nn
        type(qstate) :: psi

        allocate( psi%data(dim*nn) )

        ! generate random values with LAPACK & normalize
        call zlarnv(2, iseed, dim*nn, psi%data)
        psi%data = psi%data/sqrt( sum( real( (psi%data)*conjg(psi%data) ) ) )

        ! housekeeping
        psi%dim = dim
        psi%nn = nn
        psi%is_general = .false.
    end subroutine


    subroutine allocate_general(dim, nn, psi)
        integer :: dim, nn
        type(qstate) :: psi

        print *, '    i.e.', dim**nn, 'blochs'
        allocate( psi%data(dim**nn) )

        
        ! generate random values with LAPACK & normalize
        call zlarnv(2, iseed, dim**nn, psi%data)
        psi%data = psi%data/sqrt( sum( real( psi%data*conjg(psi%data) ) ) )

        ! housekeeping
        psi%dim = dim
        psi%nn = nn
        psi%is_general = .true.
    end subroutine



    subroutine qdeall(psi)
        type(qstate) :: psi

        deallocate( psi%data )

    end subroutine


end module