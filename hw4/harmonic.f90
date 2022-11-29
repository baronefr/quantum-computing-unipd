! ========================================
!  QIC  >>  assignment 4
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 22 November 2022
! ========================================


module harmonic
    use matrix_r8
    use mydebugger
    implicit none
    contains

    ! ------------------------------------------------------------------
    !  Init H for tridiagonal finite difference method (packed format)
    ! ------------------------------------------------------------------
    function harmonic_tridiagonal_packed(N, omega, xmax) result(H)
        integer :: N, ii
        double precision :: x, dx, xmax, omega
        type(real8_matrix) :: H

        ! check N greater than 2
        if(N.le.2) then
            call checkpoint(level=DL_FATAL, msg = 'N must be greater than 2')
            stop '(fatal error)'
        end if

        ! init packed hamiltonian to zero
        H = InitZero( (/2,N/) )

        ! set spacing & check if finite representation null
        dx = 2d0*xmax/(N-1)
        if(dx.eq.0d0) then
            call checkpoint(level=DL_FATAL, msg = 'spacing dx is null', var=dx)
            stop '(fatal error)'
        end if
        call checkpoint(level=DL_INFO, msg = 'spacing dx =', var=dx)


        ! set the matrix values ----------------------------------
        ! 
        !  [what's happening?]   the values are placed directly in 
        !                        upper-band-packed storage mode
        !
        !  ref: https://www.ibm.com/docs/en/essl/6.2?topic=representation-upper-band-packed-storage-mode
        !
        ! ----------------------------------------------------------

        do ii = 1, N
            x = xmax*(-1d0 + (ii-1)*2d0/(N-1))
            H%val(2,ii) = 2/(dx**2) + (omega**2)*(x**2)
        end do
        H%val(1,2:N) = -1/(dx**2)
        H%val = H%val * DBLE(0.5)
    end function



    ! ------------------------------------------------------------------
    !  Init H for pentadiagonal finite difference method (packed format)
    ! ------------------------------------------------------------------
    function harmonic_pentadiagonal_packed(N, omega, xmax) result(H)
        integer :: N, ii
        double precision :: x, dx, xmax, omega
        type(real8_matrix) :: H

        ! check N greater than 2
        if(N.le.2) then
            call checkpoint(level=DL_FATAL, msg = 'N must be greater than 2')
            stop '(fatal error)'
        end if

        ! init packed hamiltonian to zero
        H = InitZero( (/3,N/) )

        ! set spacing & check if finite representation null
        dx = 2d0*xmax/(N-1)
        if(dx.eq.0) then
            call checkpoint(level=DL_FATAL, msg = 'spacing dx is null', var=dx)
            stop '(fatal error)'
        end if
        call checkpoint(level=DL_INFO, msg = 'spacing dx =', var=dx)


        ! set the matrix values ----------------------------------
        ! 
        !  [what's happening?]   the values are placed directly in 
        !                        upper-band-packed storage mode
        !
        !  ref: https://www.ibm.com/docs/en/essl/6.2?topic=representation-upper-band-packed-storage-mode
        !
        ! ----------------------------------------------------------

        do ii = 1, N
            x = xmax*(-1d0 + (ii-1)*2d0/(N-1))
            H%val(3,ii) =  5/(2*(dx**2)) + (omega**2)*(x**2)
        end do
        H%val(2,2:N) = -4/(3*(dx**2))
        H%val(1,3:N) = 1/(12*(dx**2))
        H%val = H%val * DBLE(0.5)

    end function



    ! ------------------------------------------------------------------
    !  compute first k eigenvalues from packed n-banded simmetric matrix
    ! ------------------------------------------------------------------
    subroutine eigenvv_first_banded_packed(k, rmat, w, z)
        implicit none

        type(real8_matrix), intent(in) :: rmat
        integer, intent(in) :: k
        double precision, intent(out), allocatable :: w(:)
        type(real8_matrix), intent(out) :: z
        !  w  ->  array of eigenvalues
        !  z  ->  matrix  n x k  of first eigenvectors

        ! locals
        double precision :: abstol, vl, vu
        integer :: il, info, iu, j, kd, ldab, ldq, ldz, m, n
        double precision, allocatable :: q(:, :), work(:)
        integer, allocatable :: iwork(:), jfail(:)

        ! parameters
        character (1), parameter :: uplo = 'U'
        !  ->   define range of eigenvalues to find
        il = 1         ! recall: 1 <= IL <= IU <= N
        iu = k

        ldab = rmat%size(1)
        n = rmat%size(2)


        kd = ldab - 1
        ldq = n
        ldz = n
        m = iu-il+1 ! number of eigenvalues to be found (det by range)

        allocate(w(n), z%val(ldz,m), q(ldq,n), work(7*n), iwork(5*n), jfail(n))

        z%size = (/N,k/) ! manually override init of obj  
        ! Set the absolute error tolerance for eigenvalues.  
        ! With ABSTOL set to zero, the default value is used instead.
        !abstol = 0D0
        abstol = 1.0D-10

  
        ! solve the band symmetric eigenvalue problem
        call checkpoint(level = DL_PEDANTIC, msg = 'dsbevx (solver)')
        call dsbevx('V', 'I', uplo, n, kd, rmat%val, ldab, q, ldq, &
                    vl, vu, il, iu, abstol, &
                    m, w, z%val, ldz, work, iwork, jfail, info)
        
        ! [remark] from LAPACK documentation -------------
        ! about the 2nd arg of dsbevx:
        !   use 'V' to find eigenvalues in range (vl, vu]
        !   use 'I' to find eigenvalues of index  il -> iu (included)
        ! ----------------------------------------------------------

        ! [remark] from LAPACK documentation -------------
        !
        ! W is DOUBLE PRECISION array, dimension (N)
        !     The first M elements contain the selected eigenvalues in
        !     ascending order.
        !
        ! Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M))
        !     If JOBZ = 'V', then if INFO = 0, the first M columns of Z
        !     contain the orthonormal eigenvectors of the matrix A
        !     corresponding to the selected eigenvalues, with the i-th
        !     column of Z holding the eigenvector associated with W(i).
        !     If an eigenvector fails to converge, then that column of Z
        !     contains the latest approximation to the eigenvector, and the
        !     index of the eigenvector is returned in IFAIL.
        ! ----------------------------------------------------------
                    
        if(info == 0) then
            call checkpoint(level=DL_INFOP, msg = 'eigenv* ok')
        else
            call checkpoint(level=DL_FATAL, msg = 'eigenv* not computed')
            print *, " [info =", info, "]"
            stop '(fatal error)'
        end if

        deallocate(q, work, iwork, jfail)
    end subroutine

end module