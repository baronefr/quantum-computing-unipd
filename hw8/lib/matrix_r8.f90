! ========================================
!  Quantum Information and Computing
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   last update : 15 December 2022
! ========================================

module matrix_r8
!  DOC =============================================================
! 
!    This module implements a matrix with real*8 
!    (double precision) numbers.
! 
!  -----------------------------------------------------------------
! 
!  TYPES            |    real8_matrix  
!                   |       - size :  shape of the matrix
!                   |       - data :  allocated elements of matrix
! 
!  INTERFACE        |    .zeroInit.     init with zeros
! 
!  Requirements     |    use mydebugger class
! ==================================================================

    ! requirements:
    use mydebugger  ! use my checkpoints

    implicit none

    ! global seed for LAPACK generator
    integer :: iseed(4) 
    integer, parameter, private :: kkind = 8 ! sets real*8 type

    ! -------------------------------------------
    type real8_matrix
        ! store the dimension of the matrix
        integer*8, dimension(2) :: size
        ! to store the values of matrix
        real(kind=kkind), dimension(:,:), allocatable :: data
    end type

    interface operator(.zeroInit.)
        module procedure InitZero
    end interface

contains
    
    ! ------------------------------------------------------------------
    !   INITIALIZER : init a Real Matrix to zeros
    ! ------------------------------------------------------------------
    !        Input   |   shape      a 2-element array of positive 
    !                |              integers: the number of
    !                |              columns and rows
    !         
    !  Requirements  |   (be careful to not overwrite allocated objects
    !                |    with expression assignment)
    ! 
    !        Output  |   rmx        the initialized matrix
    ! 
    !    Post-cond.  |   rmx is allocated
    !  
    !  Err handlers  |   none
    ! ------------------------------------------------------------------
    function InitZero( shape ) result(rmx)
        integer*8, dimension(2), intent(in) :: shape
        type(real8_matrix) :: rmx

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        ! init properties of matrix
        allocate( rmx%data(shape(1), shape(2)) )
        rmx%size = shape

        ! writing values
        rmx%data = 0d0
    end function


    function InitIdentity(size) result(rmx)
        type(real8_matrix) :: rmx
        integer*8 :: ii, size

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        allocate( rmx%data(size, size) )
        rmx%size = (/ size, size /)

        rmx%data = 0d0 ! init to zero
        do ii = 1, rmx%size(1)
            rmx%data(ii,ii) = 1d0
        end do
    end function


    function copy_rmat(origin) result(dest)
        type(real8_matrix) :: origin, dest

        dest%size = origin%size
        allocate( dest%data( dest%size(1), dest%size(2)) )
        dest%data = origin%data
    end function


    ! ------------------------------------------------------------------
    !   DESTRUCTOR : deallocate Real Matrix
    ! ------------------------------------------------------------------
    !        Input   |   rmx      the matrix to delete from memory
    !         
    !  Requirements  |   rmx must be initialized
    ! 
    !        Output  |   none
    ! 
    !    Post-cond.  |   rmx is deallocated
    !  
    !  Err handlers  |   check rmx allocation (fatal if already dealloc)
    ! ------------------------------------------------------------------
    subroutine Delete(rmx)
        type(real8_matrix) :: rmx

        if( allocated(rmx%data) ) then
            call checkpoint(level=DL_PEDANTIC, msg = 'deallocating memory')
            deallocate( rmx%data )
            rmx%size(1) = 0;   rmx%size(2) = 0;
        else
            call checkpoint(level=DL_FATAL, msg = 'cannot deallocate rmat, already cleared')
            stop '(fatal error)'
        end if
    end subroutine







    ! ------------------------------------------------------------------
    !   MATH |  TENSOR PRODUCT
    ! ------------------------------------------------------------------
    function kron(rmxx, rmyy) result(rmo)
        type(real8_matrix) :: rmxx, rmyy
        type(real8_matrix) :: rmo
        integer*8 :: ii, jj

        rmo%size(1) = rmxx%size(1)*rmyy%size(1)
        rmo%size(2) = rmxx%size(2)*rmyy%size(2)

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        allocate( rmo%data( rmo%size(1), rmo%size(2) ) )

        rmo%data = (0d0, 0d0)

        do ii = 1, rmxx%size(1)
            do jj = 1, rmxx%size(2)
                rmo%data( (ii-1)*rmyy%size(1) + 1 : ii*rmyy%size(1) , &
                          (jj-1)*rmyy%size(2) + 1 : jj*rmyy%size(2) ) &
                        = rmxx%data(ii,jj) * rmyy%data
            end do
        end do
    end function



    ! ------------------------------------------------------------------
    !   MATH |  TENSOR PRODUCT
    ! ------------------------------------------------------------------
    subroutine kron_inplace(rmxx, rmyy, rmo)
        type(real8_matrix) :: rmxx, rmyy
        type(real8_matrix) :: rmo
        integer*8 :: ii, jj

        
        ! check rmo size
        call checkpoint(level=DL_PEDANTIC, msg = 'check inplace dimensions')
        if( rmo%size(1) .ne. rmxx%size(1)*rmyy%size(1) ) then
            call checkpoint(level=DL_FATAL, msg = 'inplace dimension 1 does not match')
            stop 'err'
        end if
        if( rmo%size(2) .ne. rmxx%size(2)*rmyy%size(2) ) then
            call checkpoint(level=DL_FATAL, msg = 'inplace dimension 1 does not match')
            stop 'err'
        end if


        do ii = 1, rmxx%size(1)
            do jj = 1, rmxx%size(2)
                rmo%data( (ii-1)*rmyy%size(1) + 1 : ii*rmyy%size(1) , &
                          (jj-1)*rmyy%size(2) + 1 : jj*rmyy%size(2) ) &
                        =  rmo%data( (ii-1)*rmyy%size(1) + 1 : ii*rmyy%size(1) , &
                        (jj-1)*rmyy%size(2) + 1 : jj*rmyy%size(2) ) &
                        + rmxx%data(ii,jj) * rmyy%data
            end do
        end do
    end subroutine


    ! ------------------------------------------------------------------
    !   MATH |  TENSOR PRODUCT with left IDENTITY
    ! ------------------------------------------------------------------
    function kron_l1(idsize, rmyy) result(rmo)
        integer*8, intent(in) :: idsize
        type(real8_matrix), intent(in) :: rmyy
        type(real8_matrix) :: rmo
        integer*8 :: ii, jj

        rmo%size(1) = idsize*rmyy%size(1)
        rmo%size(2) = idsize*rmyy%size(2)

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        allocate( rmo%data( rmo%size(1), rmo%size(2) ) )

        rmo%data = (0d0, 0d0)

        do ii = 1, idsize
            rmo%data( (ii-1)*rmyy%size(1) + 1 : ii*rmyy%size(1) , &
                      (ii-1)*rmyy%size(2) + 1 : ii*rmyy%size(2) ) &
                    = rmyy%data
        end do
    end function

    ! ------------------------------------------------------------------
    !   MATH |  TENSOR PRODUCT with left IDENTITY
    ! ------------------------------------------------------------------
    subroutine kron_l1_inplace(idsize, rmyy, rmo)
        integer*8, intent(in) :: idsize
        type(real8_matrix), intent(in) :: rmyy
        type(real8_matrix) :: rmo
        integer*8 :: ii

        ! check rmo size
        call checkpoint(level=DL_PEDANTIC, msg = 'check inplace dimensions')
        if( rmo%size(1) .ne. rmyy%size(1)*idsize ) then
            call checkpoint(level=DL_FATAL, msg = 'inplace dimension 1 does not match')
            stop 'err'
        end if

        if( rmo%size(2) .ne. rmyy%size(2)*idsize ) then
            call checkpoint(level=DL_FATAL, msg = 'inplace dimension 2 does not match')
            stop 'err'
        end if

        do ii = 1, idsize
            rmo%data( (ii-1)*rmyy%size(1) + 1 : ii*rmyy%size(1) , &
                      (ii-1)*rmyy%size(2) + 1 : ii*rmyy%size(2) ) &
                    =  rmo%data( (ii-1)*rmyy%size(1) + 1 : ii*rmyy%size(1) , &
                    (ii-1)*rmyy%size(2) + 1 : ii*rmyy%size(2) ) + rmyy%data
        end do
    end subroutine

    ! ------------------------------------------------------------------
    !   MATH |  TENSOR PRODUCT with right IDENTITY
    ! ------------------------------------------------------------------
    function kron_r1(rmxx, idsize) result(rmo)
        type(real8_matrix), intent(in) :: rmxx
        integer*8, intent(in) :: idsize
        type(real8_matrix) :: rmo

        integer*8 :: ii, jj, kk

        rmo%size(1) = rmxx%size(1)*idsize
        rmo%size(2) = rmxx%size(2)*idsize

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        allocate( rmo%data( rmo%size(1), rmo%size(2) ) )

        rmo%data = (0d0, 0d0)

        do ii = 1, rmxx%size(1)
            do jj = 1, rmxx%size(2)
                do kk = 1, idsize
                    rmo%data( (ii-1)*idsize + kk, (jj-1)*idsize + kk ) = rmxx%data(ii,jj)
                end do
            end do
        end do
    end function




    ! ------------------------------------------------------------------
    !   MATH |  TENSOR PRODUCT with right IDENTITY
    ! ------------------------------------------------------------------
    subroutine kron_r1_inplace(rmxx, idsize, rmo)
        type(real8_matrix), intent(in) :: rmxx
        integer*8, intent(in) :: idsize
        type(real8_matrix) :: rmo

        integer*8 :: ii, jj, kk

        ! check rmo size
        call checkpoint(level=DL_PEDANTIC, msg = 'check inplace dimensions')
        if( rmo%size(1) .ne. rmxx%size(1)*idsize ) then
            call checkpoint(level=DL_FATAL, msg = 'inplace dimension 1 does not match')
            stop 'err'
        end if

        if( rmo%size(2) .ne. rmxx%size(2)*idsize ) then
            call checkpoint(level=DL_FATAL, msg = 'inplace dimension 2 does not match')
            stop 'err'
        end if

        do ii = 1, rmxx%size(1)
            do jj = 1, rmxx%size(2)
                do kk = 1, idsize
                    rmo%data( (ii-1)*idsize + kk, (jj-1)*idsize + kk ) = &
                        rmo%data( (ii-1)*idsize + kk, (jj-1)*idsize + kk ) + rmxx%data(ii,jj)
                end do
            end do
        end do
    end subroutine





    ! ------------------------------------------------------------------
    !   MATH |  TRACE    Compute the trace of a real Matrix
    ! ------------------------------------------------------------------
    !        Input   |   rmx    real matrix
    !         
    !  Requirements  |   rmx must be square
    ! 
    !        Output  |   tr     (real)   trace of rmx
    ! 
    !    Post-cond.  |   none
    !  
    !  Err handlers  |   check if rmx is square (error)
    ! ------------------------------------------------------------------
    function Trace(rmx) result(tr)
        type(real8_matrix), intent(in) :: rmx
        real(kind=kkind) :: tr
        integer*8 :: ii

        tr = 0d0 ! init to zero before loop

        if( rmx%size(1) .eq. rmx%size(2) ) then
            do ii = 1, rmx%size(1)
                tr = tr + rmx%data(ii,ii)
            end do
        else
            call checkpoint(level=DL_ERROR, msg = 'input matrices bust be square to compute trace')
        end if
    end function

    



    ! ------------------------------------------------------------------
    !   MATH |  compute EIGENVALUES and EIGENVECTORS of tridiag matrix
    ! ------------------------------------------------------------------
    !        Input   |   rmx         real tridiagonal matrix (square)
    ! 
    !  Requirements  |   rmx must be initialized
    ! 
    !        Output  |   eigval      the real array of eigenvalues
    !                |               in ascending order
    ! 
    !    Post-cond.  |   an array of eigenvalues is allocated
    !  
    !  Err handlers  |   function throws fatal exception if eigenvalues
    !                |   are not computed (any reason in dstev() doc)
    ! ------------------------------------------------------------------
    subroutine eigenvv_symmetric(k, rmat, w, z, do_eigvec)
        implicit none
        type(real8_matrix), intent(in) :: rmat
        integer, intent(in) :: k
        double precision, intent(out), allocatable :: w(:)
        type(real8_matrix), intent(inout) :: z
        logical, optional :: do_eigvec ! if true, compute eigenvectors
        !  w  ->  array of eigenvalues
        !  z  ->  matrix  n x k  of first eigenvectors

        ! parameters
        integer, parameter :: nb = 64
        double precision :: abstol, vl, vu
        integer :: il, info, iu, lda, ldz, lwork, m, n

        ! tmp variables
        double precision, allocatable :: work(:)
        double precision :: dummy(1)
        integer, allocatable :: iwork(:), jfail(:)
        intrinsic :: max, nint
        character(1) :: jobz

        if (present( do_eigvec )) then
            if( do_eigvec ) then
                jobz = 'V'
            else
                call checkpoint(level = DL_INFO, msg = 'no eigenvectors will be computed')
                jobz = 'N'
            end if
        else
            jobz = 'V' ! default: do eigenvectors
        end if

        !  ->   define range of eigenvalues to find
        il = 1         ! recall: 1 <= IL <= IU <= N
        iu = k

        n   = rmat%size(2)
        lda = rmat%size(1)
        ldz = n
        m = iu-il+1 ! number of eigenvalues to be found (det by range)

        call checkpoint(level = DL_PEDANTIC, msg = 'allocate tmp for dsyevx')
        if ( jobz == 'V') then
            allocate( z%data(rmat%size(2), m) )
        end if
        allocate( w(n), iwork(5*n), jfail(n) )
        z%size = (/rmat%size(2), int8(m)/)     ! manually override init of obj  
        

        ! Use routine workspace query to get optimal workspace.
        lwork = -1                   ! to make a query
        call dsyevx(jobz, 'I', 'U', n, rmat%data, lda, &
                    vl, vu, il, iu, abstol, &
                    m, w, z%data, ldz, dummy, lwork, iwork, jfail, info)
        lwork = max((nb+3)*n, nint(dummy(1)))
        allocate(  work(lwork)  )
    
        ! Set the absolute error tolerance for eigenvalues.  With ABSTOL
        ! set to zero, the default value is used instead.
        abstol = 0d0

        ! eigenvalue problem solver
        call checkpoint(level = DL_PEDANTIC, msg = 'dsyevx (solver)')
        call dsyevx(jobz, 'I', 'U', n, rmat%data, lda, &
                    vl, vu, il, iu, abstol, &
                    m, w, z%data, ldz, work, lwork, iwork, jfail, info)

        ! [ ref ] ----------------- 
        ! https://netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga68612cdf4ed1051c08f0b0735b8dfdea.html
        
        ! check solver status
        if (info>=0) then
            if (info>0) then
                call checkpoint(level=DL_WARN, msg = 'eigenv* not converged')
                print *, 'INFO eigenvectors failed to converge, INFO =', info
                print *, 'Indices of eigenvectors that did not converge'
                print *, jfail(1:m)

            else
                call checkpoint(level=DL_INFOP, msg = 'eigenv* ok')
                !print *,  'Number of eigenvalues found =', m
            end if
        else
            call checkpoint(level=DL_FATAL, msg = 'eigenv* not computed')
            print *, " [info =", info, "]"
            stop '(fatal error)'
        end if

        deallocate(work, iwork, jfail)
    end subroutine

    
    ! ------------------------------------------------------------------
    !   I/O |  DUMP    Write Real Matrix to binary file (with records)
    ! ------------------------------------------------------------------
    !        Input   |   rmx        real matrix
    !                |
    !                |   filename   (char*) name of the output file
    !         
    !  Requirements  |   none, but check permission to write file
    ! 
    !        Output  |   none
    ! 
    !    Post-cond.  |   input matrix is written to file
    !  
    !  Err handlers  |   check if file is opened (error)
    ! ------------------------------------------------------------------
    subroutine RMatDump(rmx, filename)
        type(real8_matrix), intent(in) :: rmx
        character(*), intent(in) :: filename
        integer, parameter :: fileunit = 22
        logical :: file_isopened

        open(fileunit, file=filename, action='write', form="unformatted")
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then
            call checkpoint(level=DL_INFO, msg = 'writing RMat to (binary)'//filename)

            write(fileunit) rmx%size
            write(fileunit) rmx%data
        else
            call checkpoint(level=DL_ERROR, msg = 'file '//filename//' cannot be opened')
            return
        end if
 
        close(fileunit)
        call checkpoint(level=DL_INFOP, msg = 'file written '//filename)
    end subroutine



    ! ------------------------------------------------------------------
    !   I/O |  DUMP    Read Real Matrix from binary file (with records)
    ! ------------------------------------------------------------------
    function RMatRead(filename) result(rmx)
        type(real8_matrix) :: rmx
        character(*), intent(in) :: filename
        integer, parameter :: fileunit = 22
        logical :: file_isopened
        integer*8, dimension(2) :: shape

        open(fileunit, file=filename, action='read', form="unformatted")
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then
            call checkpoint(level=DL_INFO, msg = 'read RMat from (binary)'//filename)

            read(fileunit) shape

            allocate( rmx%data(shape(1), shape(2)) )
            rmx%size = shape
    
            ! reading matrix values
            read(fileunit) rmx%data
        else
            call checkpoint(level=DL_ERROR, msg = 'file '//filename//' cannot be opened')
            return
        end if
 
        close(fileunit)
        call checkpoint(level=DL_INFOP, msg = 'file closed '//filename)
    end function


    ! ------------------------------------------------------------------
    !   I/O |  DUMPTXT    Write Real Matrix to file (csv format)
    ! ------------------------------------------------------------------
    !        Input   |   rmx        real matrix
    !                |
    !                |   filename   (char*) name of the output file
    !         
    !  Requirements  |   none, but check permission to write file
    ! 
    !        Output  |   none
    ! 
    !    Post-cond.  |   input matrix is written to file
    !  
    !  Err handlers  |   check if file is opened (error)
    ! ------------------------------------------------------------------
    subroutine RMatDumpTXT(rmx, filename)
        type(real8_matrix), intent(in) :: rmx
        character(*), intent(in) :: filename
        integer, parameter :: fileunit = 22
        logical :: file_isopened
        integer*8 :: ii

        open(fileunit, file=filename, action='write')
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then
            call checkpoint(level=DL_INFO, msg = 'writing RMat to '//filename)

            ! write file header:  size of matrix
            write(fileunit, *) '# [real8_matrix] %size:', rmx%size

            ! write to file
            do ii = 1, rmx%size(1)
                write (fileunit,"(*(G0,:,','))") rmx%data(:,ii)
                ! write line by line with comma separation
            end do
        else
            call checkpoint(level=DL_ERROR, msg = 'file '//filename//' cannot be opened')
            return
        end if
 
        close(fileunit)
        call checkpoint(level=DL_INFOP, msg = 'file written '//filename)
    end subroutine


    ! ------------------------------------------------------------------
    !   I/O |  Print Head    Print the first nrow and ncol of a matrix
    ! ------------------------------------------------------------------
    !        Input   |   rmx        real matrix
    !                |
    !                |   nrow   (integer) index of last row to print
    !                |
    !                !   ncol   (integer) index of lareal8_matrixst column to print
    !         
    !  Requirements  |   rmx must be initialized
    ! 
    !        Output  |   none, just print in stdout
    ! 
    !    Post-cond.  |   none
    !  
    !  Err handlers  |   none
    ! ------------------------------------------------------------------
    subroutine RMatPrintHead(rmx, nrow, ncol)
        type(real8_matrix), intent(in) :: rmx
        integer*8, optional :: nrow, ncol
        integer*8 :: nrow_, ncol_
        integer*8 :: ii

        if( present(nrow) ) then
            nrow_ = nrow
        else
            nrow_ = rmx%size(1)
        end if

        if( present(ncol) ) then
            ncol_ = ncol
        else
            ncol_ = rmx%size(2)
        end if

        do ii = 1, nrow_
            print "(*(F0.3,:,'  '))", rmx%data(ii,:ncol_)
        end do

        ! compact format for complexes...
        !do ii = 1,nrow_
        !    write(*, "(*('('sf6.2xspf6.2x'i)':x))") rmx%data(ii,:ncol_)
        !end do
    end subroutine


    
end module matrix_r8
