! ========================================
!  Quantum Information and Computing
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   last update : 21 November 2022
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
!                   |       - val  :  allocated elements of matrix
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
        integer, dimension(2) :: size
        ! to store the values of matrix
        real(kind=kkind), dimension(:,:), allocatable :: val
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
        integer, dimension(2), intent(in) :: shape
        type(real8_matrix) :: rmx

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        ! init properties of matrix
        allocate( rmx%val(shape(1), shape(2)) )
        rmx%size = shape

        ! writing values
        rmx%val = 0d0
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

        if( allocated(rmx%val) ) then
            call checkpoint(level=DL_PEDANTIC, msg = 'deallocating memory')
            deallocate( rmx%val )
            rmx%size(1) = 0;   rmx%size(2) = 0;
        else
            call checkpoint(level=DL_FATAL, msg = 'cannot deallocate rmat, already cleared')
            stop '(fatal error)'
        end if
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
        integer :: ii

        tr = 0d0 ! init to zero before loop

        if( rmx%size(1) .eq. rmx%size(2) ) then
            do ii = 1, rmx%size(1)
                tr = tr + rmx%val(ii,ii)
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
    function eigenvv_tridiag(rmx, eigvec) result(eigval)
        implicit none
        ! args
        type(real8_matrix) :: rmx
        logical, optional :: eigvec
        real(kind=kkind), dimension(:), allocatable :: eigval

        ! zheev parameters & other stuff
        integer :: info, ldz, n, ii
        real(kind=kkind), dimension(:), allocatable :: work, e

        ! ARG HANDLERS ---------------------------------------------
        character*1 :: eigvec_ = 'N'  ! default: do not compute eigenvectors

        ! arg handler: if eigvec is True, eigenvectors are computed (default: False)
        if(present(eigvec)) then
            if(eigvec) then
                eigvec_ = 'V' ! compute eigenvalues and eigenvectors
            else
                eigvec_ = 'N' ! compute eigenvalues only
            end if
        end if

        call checkpoint(level=DL_PEDANTIC, msg = 'eigvec_ =', var=eigvec_)
        
        ! [remark]  this procedure is shamelessly inspired from LAPACK documentation
        ! https://numericalalgorithmsgroup.github.io/LAPACK_Examples/examples/doc/dstev_example.html
        ! ----------------------------------------------------------

        if( rmx%size(1) .eq. rmx%size(2) ) then
            n = rmx%size(1);   ldz = rmx%size(2);
            allocate( eigval(n), e(n-1), work( max(1,2*n-2) ) )

            ! init eigval (D in LAPACK doc) | the n diagonal elements of the tridiagonal matrix
            do ii = 1,n
                eigval(ii) = rmx%val(ii,ii)
            end do

            ! init e | the (n-1) subdiagonal elements of the tridiagonal matrix
            do ii = 1,n-1
                e(ii) = rmx%val(ii,ii+1)
            end do

            ! solve the Hermitian eigenvalue problem
            call checkpoint(level=DL_PEDANTIC, msg = 'dstev (solver)')
            call dstev(eigvec_, n, eigval, e, rmx%val, ldz, work, info)
            ! ref: https://netlib.org/lapack/explore-html/dc/dd2/group__double_o_t_h_e_reigen_gaaa6df51cfd92c4ab08d41a54bf05c3ab.html

            ! [remark] from LAPACK documentation (adapted) -------------
            !
            !  If eigvec_ = 'V', then if INFO = 0, RMX contains the orthonormal
            !  eigenvectors of the matrix RMX, with the i-th column of RMX
            !  holding the eigenvector associated with eigval(i).
            !
            ! ----------------------------------------------------------

            if(info == 0) then
                call checkpoint(level=DL_INFOP, msg = 'eigenv* ok')
            else
                call checkpoint(level=DL_FATAL, msg = 'eigenv* not computed')
                print *, " [info =", info, "]"
                stop '(fatal error)'
            end if

            deallocate(e, work)
        else
            call checkpoint(level=DL_ERROR, msg = 'input matrices bust be square to compute eigenvalues')
        end if
    end function

    
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
        integer :: ii

        open(fileunit, file=filename, action='write', form="unformatted")
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then
            call checkpoint(level=DL_INFO, msg = 'writing RMat to (binary)'//filename)

            write(fileunit) rmx%size
            write(fileunit) rmx%val
        else
            call checkpoint(level=DL_ERROR, msg = 'file '//filename//' cannot be opened')
            return
        end if
 
        close(fileunit)
        call checkpoint(level=DL_INFOP, msg = 'file written '//filename)
    end subroutine


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
        integer :: ii

        open(fileunit, file=filename, action='write')
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then
            call checkpoint(level=DL_INFO, msg = 'writing RMat to '//filename)

            ! write file header:  size of matrix
            write(fileunit, *) '# [real8_matrix] %size:', rmx%size

            ! write to file
            do ii = 1, rmx%size(1)
                write (fileunit,"(*(G0,:,','))") rmx%val(:,ii)
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
        integer :: ii, nrow, ncol

        do ii = 1,nrow
            print *, rmx%val(ii,:ncol)
        end do
    end subroutine


    
end module matrix_r8
