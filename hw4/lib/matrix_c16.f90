! ========================================
!  Quantum Information and Computing
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   last update : 10 November 2022
! ========================================

module matrix_c16
!  DOC =============================================================
! 
!    This module implements a matrix with complex*16 numbers.
! 
!  -----------------------------------------------------------------
! 
!  TYPES            |    complex16_matrix  
!                   |       - size :  shape of the matrix
!                   |       - val  :  allocated elements of matrix
! 
!  INTERFACE        |    .randInit.     init with rands
!                   |    .zeroInit.     init with zeros
!                   |    .Adj.          compute the adjoint
!                   |    .Tr.           compute the trace
! 
!  Requirements     |    use mydebugger class
! ==================================================================

    ! requirements:
    use mydebugger  ! use my checkpoints

    implicit none

    ! global seed for LAPACK generator
    integer :: iseed(4) 

    ! -------------------------------------------
    type complex16_matrix
        ! store the dimension of the matrix
        integer, dimension(2) :: size
        ! to store the values of matrix
        complex*16, dimension(:,:), allocatable :: val
    end type

    interface operator(.randInit.)
        module procedure InitRand
    end interface

    interface operator(.zeroInit.)
        module procedure InitZero
    end interface

    interface operator(.Adj.)
        module procedure Adjoint
    end interface

    interface operator(.Tr.)
        module procedure Trace
    end interface

contains
    
    ! ------------------------------------------------------------------
    !   INITIALIZER : init a Complex Matrix to Rand values
    ! ------------------------------------------------------------------
    !        Input   |   shape      a 2-element array of positive 
    !                |              integers: the number of
    !                |              columns and rows
    !         
    !  Requirements  |   (be careful to not overwrite allocated objects
    !                |    with expression assignment)
    ! 
    !        Output  |   cmx        the initialized matrix
    ! 
    !    Post-cond.  |   cmx is allocated
    !  
    !  Err handlers  |   shape has to be two positive integers
    ! ------------------------------------------------------------------
    function InitRand( shape ) result(cmx)
        integer, dimension(2), intent(in) :: shape
        type(complex16_matrix) :: cmx
        integer :: ii, jj

        if( .not.((shape(1).gt.0) .and. (shape(2).gt.0)) ) then
            call checkpoint(level=DL_FATAL, msg = 'shapes must be positive')
            stop '(fatal error)'
        end if

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        allocate( cmx%val(shape(1), shape(2)) )
        cmx%size = shape

        do jj=1,cmx%size(2)
            do ii=1,cmx%size(1)
                cmx%val(jj,ii) = complex( rand(), rand() )
            end do
        end do
    end function

    ! ------------------------------------------------------------------
    !   INITIALIZER : init a diagonal Complex Matrix to rand real values
    ! ------------------------------------------------------------------
    function InitRandDiagonal( size ) result(cmx)
        integer, intent(in) :: size
        type(complex16_matrix) :: cmx
        integer :: jj

        if( size.le.1 ) then
            call checkpoint(level=DL_FATAL, msg = 'size must be > 1')
            stop '(fatal error)'
        end if

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        allocate( cmx%val(size, size) )
        cmx%size = (/size, size/)

        cmx%val = complex(0d0, 0d0)
        do jj=1,size
            cmx%val(jj,jj) = complex( rand(), 0d0)
        end do
    end function

    ! ------------------------------------------------------------------
    !   INITIALIZER : init a Complex Matrix to 0+0i
    ! ------------------------------------------------------------------
    !        Input   |   shape      a 2-element array of positive 
    !                |              integers: the number of
    !                |              columns and rows
    !         
    !  Requirements  |   (be careful to not overwrite allocated objects
    !                |    with expression assignment)
    ! 
    !        Output  |   cmx        the initialized matrix
    ! 
    !    Post-cond.  |   cmx is allocated
    !  
    !  Err handlers  |   shape has to be two positive integers
    ! ------------------------------------------------------------------
    function InitZero( shape ) result(cmx)
        integer, dimension(2), intent(in) :: shape
        type(complex16_matrix) :: cmx

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        ! init properties of cmatrix
        allocate( cmx%val(shape(1), shape(2)) )
        cmx%size = shape

        ! writing values
        cmx%val = complex( 0d0, 0d0)
    end function

    ! ------------------------------------------------------------------
    !   INITIALIZER : init a diagonal Matrix with real values 
    ! ------------------------------------------------------------------
    !        Input   |   size  [int]   the number of columns and rows
    !                |                 of the matrix
    !                |   dist [char]   type of distribution to sample
    !                |                  if 'u' -> uniform in (0,1)
    !                |                  if 's' -> uniform in (-1,1)
    !                |                  if 'n' -> normal in 0 and std 1
    ! 
    !  Requirements  |   (be careful to not overwrite allocated objects
    !                |    with expression assignment)
    ! 
    !        Output  |   cmx        the initialized matrix
    ! 
    !    Post-cond.  |   cmx is allocated
    !  
    !  Err handlers  |   shape has to be two positive integers
    ! ------------------------------------------------------------------
    function InitDiagonalReal(size, dist) result(cmx)
        ! I/O arguments
        integer, intent(in) :: size
        character(1) :: dist
        type(complex16_matrix) :: cmx
        
        integer :: idist, jj
        double precision, dimension(:), allocatable :: tmp

        select case (dist)
            case ('u') ! uniform in (0,1)
                idist = 1
            case ('s') ! uniform in (-1,1)
                idist = 2
            case ('n') ! normal in 0 with std 1
                idist = 3
            case DEFAULT
                idist = 1
        end select
        call checkpoint(level=DL_PEDANTIC, msg = 'handle dist selection', var=idist)

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating cmat')
        ! init properties of cmatrix
        allocate( cmx%val(size, size), tmp(size))
        cmx%val = complex( 0d0, 0d0)
        cmx%size = (/size, size/)

        ! generating rand
        call checkpoint(level=DL_PEDANTIC, msg = 'generate rand values with LAPACK')
        call dlarnv(idist, iseed, size, tmp)

        ! placing values in diagonal
        do jj=1,size
            cmx%val(jj,jj) = dcmplx(tmp(jj))
        end do

        deallocate(tmp)
        call checkpoint(level=DL_INFOP, msg = 'rand done')
    end function

    ! ------------------------------------------------------------------
    !   INITIALIZER : init an Hermitian Complex Matrix
    ! ------------------------------------------------------------------
    !        Input   |   size    [integer]   dimension of square matrix
    !                |   dist    [char]      type of distribution to sample
    !                |                      if 'u' -> uniform in (0,1)
    !                |                      if 's' -> uniform in (-1,1)
    !                |                      if 'n' -> normal in 0 and std 1
    !                |                      if 'd' -> disc abs(z) < 1
    !                |                      if 'c' -> circle abs(z) = 1
    ! 
    !  Requirements  |   (be careful to not overwrite allocated objects
    !                |    with expression assignment)
    ! 
    !        Output  |   cmx        the initialized matrix
    ! 
    !    Post-cond.  |   cmx is allocated
    !  
    !  Err handlers  |   size has to be a positive integer
    ! ------------------------------------------------------------------
    function InitHermitian(size, dist) result(cmx)
        ! I/O arguments
        integer, intent(in) :: size
        type(complex16_matrix) :: cmx
        character(1) :: dist

        complex*16, dimension(:), allocatable :: tmp
        integer :: idist, jj

        select case (dist)
        case ('u') ! uniform in (0,1)
            idist = 1
        case ('s') ! uniform in (-1,1)
            idist = 2
        case ('n') ! normal in (0,1)
            idist = 3
        case ('d') ! disc abs(z) < 1
            idist = 4
        case ('c') ! circle abs(z) = 1
            idist = 5
        case DEFAULT
            idist = 1
        end select

        call checkpoint(level=DL_PEDANTIC, msg = 'handle dist selection', var=idist)

        if(size.le.1) then
            call checkpoint(level=DL_FATAL, msg = 'size must be 2 or greater')
            stop '(fatal error)'
        end if

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        allocate( cmx%val(size, size), tmp(size) )
        cmx%size = (/size, size/)

        call checkpoint(level=DL_INFO, msg = 'generate rand values with LAPACK')
        do jj = 1,size
            ! generate random numbers to tmp array
            call zlarnv(idist, iseed, jj, tmp)

            ! copy from tmp to columns & rows (conjugate)
            cmx%val(1:jj,jj)   =        tmp( 1:jj )
            cmx%val(jj,1:jj-1) = conjg( tmp( 1:jj ) )

            ! fix diagonal to real
            cmx%val(jj,jj) = realpart(cmx%val(jj,jj))
        end do

        deallocate(tmp)
        call checkpoint(level=DL_INFOP, msg = 'rand done')
    end function



    ! ------------------------------------------------------------------
    !   DESTRUCTOR : deallocate Complex Matrix
    ! ------------------------------------------------------------------
    !        Input   |   cmx      the matrix to delete from memory
    !         
    !  Requirements  |   cmx must be initialized
    ! 
    !        Output  |   none
    ! 
    !    Post-cond.  |   cmx is deallocated
    !  
    !  Err handlers  |   check cmx allocation (fatal)
    ! ------------------------------------------------------------------
    subroutine Delete(cmx)
        type(complex16_matrix) :: cmx

        if( allocated(cmx%val) ) then
            call checkpoint(level=DL_PEDANTIC, msg = 'deallocating memory')
            deallocate( cmx%val )
            cmx%size(1) = 0;   cmx%size(2) = 0;
        else
            call checkpoint(level=DL_FATAL, msg = 'cannot deallocate cmat, already cleared')
            stop '(fatal error)'
        end if
    end subroutine


    ! ------------------------------------------------------------------
    !   MATH |  ADJOINT    Compute the adjoint of a Complex Matrix
    ! ------------------------------------------------------------------
    !        Input   |   cmx        complex matrix
    !         
    !  Requirements  |   cmx must be initialized
    ! 
    !        Output  |   cmxadj     the adjoint matrix
    ! 
    !    Post-cond.  |   none
    !  
    !  Err handlers  |   none
    ! ------------------------------------------------------------------
    function Adjoint(cmx) result(cmxadj)
        type(complex16_matrix), intent(in) :: cmx
        type(complex16_matrix) :: cmxadj

        cmxadj%size(1) = cmx%size(2);    cmxadj%size(2) = cmx%size(1);

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        allocate( cmxadj%val(cmxadj%size(1),cmxadj%size(2)) )

        cmxadj%val = conjg(transpose(cmx%val))
    end function

    ! ------------------------------------------------------------------
    !   MATH |  TRACE    Compute the trace of a Complex Matrix
    ! ------------------------------------------------------------------
    !        Input   |   cmx    complex matrix
    !         
    !  Requirements  |   cmx must be square
    ! 
    !        Output  |   tr     (complex*16)   trace of cmx
    ! 
    !    Post-cond.  |   none
    !  
    !  Err handlers  |   check if cmx is square (error)
    ! ------------------------------------------------------------------
    function Trace(cmx) result(tr)
        type(complex16_matrix), intent(in) :: cmx
        complex*16 :: tr
        integer :: ii

        tr = complex(0d0,0d0) ! init to zero before loop

        if( cmx%size(1) .eq. cmx%size(2) ) then
            do ii = 1, cmx%size(1)
                tr = tr + cmx%val(ii,ii)
            end do
        else
            call checkpoint(level=DL_ERROR, msg = 'input matrices bust be square to compute trace')
        end if
    end function

    ! ------------------------------------------------------------------
    !   MATH |  compute EIGENVALUES and EIGENVECTORS of hermitian matrix
    ! ------------------------------------------------------------------
    !        Input   |   cmx         complex hermitian matrix (square)
    !                |   uselower    [logical]    if True, zheev() will
    !                |                            use the upper 
    !                |                            triangular matrix
    ! 
    !  Requirements  |   cmx must be initialized
    ! 
    !        Output  |   eigval        the double array of eigenvalues
    ! 
    !    Post-cond.  |   an array of eigenvalues is allocated
    !  
    !  Err handlers  |   function throws fatal exception if eigenvalues
    !                |   are not computed (any reason in zheev() doc)
    ! ------------------------------------------------------------------
    function eigenvv_hermitian(cmx, eigvec, uselower) result(eigval)
        implicit none
        ! args
        type(complex16_matrix) :: cmx
        logical, optional :: eigvec
        logical, optional :: uselower
        double precision, dimension(:), allocatable :: eigval

        ! zheev parameters & other stuff
        integer, parameter :: nb = 64
        integer :: info, lda, lwork, n
        complex*16, allocatable :: work(:)
        complex*16 :: dummy(1)
        double precision, allocatable :: rwork(:)
        

        ! ARG HANDLERS ---------------------------------------------
        character*1 :: uplo_ = 'U'    ! default: use upper triangular
        character*1 :: eigvec_ = 'N'  ! default: do not compute eigenvectors

        ! arg handler: if eigvec is True, eigenvectors are computed (default: False)
        if(present(eigvec)) then
            if(eigvec) then
                eigvec_ = 'V'
            else
                eigvec_ = 'N'
            end if
        end if

        ! arg handler: if eigvec is False, either upper/lower triangular is destroyed
        if(present(uselower)) then
            if(uselower) then
                uplo_ = 'L'
            else
                uplo_ = 'U'
            end if
        end if

        ! [remark] from LAPACK documentation (adapted) -------------
        !
        !   With A being cmx%val,
        !
        !  On exit, if eigvec_ = 'V', then if INFO = 0, A contains the
        !  orthonormal eigenvectors of the matrix A.
        !  If eigvec_ = 'N', then on exit the lower triangle (if uplo_='L')
        !  or the upper triangle (if uplo_='U') of A, including the
        !  diagonal, is destroyed.
        ! 
        ! ----------------------------------------------------------
        call checkpoint(level=DL_PEDANTIC, msg = 'eigvec_ =', var=eigvec_)
        call checkpoint(level=DL_PEDANTIC, msg = 'uplo_ =', var=uplo_)

        
        ! [remark]  this procedure is shamelessly inspired from LAPACK documentation
        ! https://numericalalgorithmsgroup.github.io/LAPACK_Examples/examples/doc/zheev_example.html
        ! ----------------------------------------------------------

        if( cmx%size(1) .eq. cmx%size(2) ) then
            n = cmx%size(1);   lda = n;
            allocate( rwork(3*n-2), eigval(n))
            
            ! use routine workspace query to get optimal workspace
            call checkpoint(level=DL_PEDANTIC, msg = 'zheev (optimal workspace)')
            lwork = -1
            call zheev(eigvec_, uplo_, n, cmx%val, lda, eigval, dummy, lwork, rwork, info)
            
            ! make sure that there is enough workspace for block size nb
            lwork = max((nb+1)*n, nint(real(dummy(1))))
            allocate(work(lwork))
            
            ! solve the Hermitian eigenvalue problem
            call checkpoint(level=DL_PEDANTIC, msg = 'zheev (solver)')
            call zheev(eigvec_, uplo_, n, cmx%val, lda, eigval, work, lwork, rwork, info)

            if(info==0) then
                call checkpoint(level=DL_INFOP, msg = 'eigenv* ok')
            else
                call checkpoint(level=DL_FATAL, msg = 'eigenv* not computed')
                print *, " [info =", info, "]"
                stop '(fatal error)'
            end if

            deallocate(rwork, work)
        else
            call checkpoint(level=DL_ERROR, msg = 'input matrices bust be square to compute eigenvalues')
        end if
    end function


    ! ------------------------------------------------------------------
    !   I/O |  DUMPTXT    Write Complex Matrix to file
    ! ------------------------------------------------------------------
    !        Input   |   cmx        complex matrix
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
    subroutine CMatDumpTXT(cmx, filename)
        type(complex16_matrix), intent(in) :: cmx
        character(*), intent(in) :: filename
        integer, parameter :: fileunit = 22
        logical :: file_isopened
        integer :: ii

        open(fileunit, file=filename, action='write')
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then
            call checkpoint(level=DL_INFO, msg = 'writing CMat to '//filename)

            ! write file header:  size of matrix
            write(fileunit, *) '# [complex16_matrix] %size:', cmx%size

            ! write to file
            do ii = 1, cmx%size(1)

                write (fileunit,"(*(G0,:,','))") cmx%val(:,ii) ! full size
                ! [ remark] ----------------------------------------
                ! each complex value is written as:   realpart,imgpart
                ! so the file looks like a matrix of 2*size(2) columns and size(1) rows
                ! --------------------------------------------------
            end do
        else
            call checkpoint(level=DL_ERROR, msg = 'file '//filename//' cannot be opened')
            return
        end if
 
        close(fileunit)
        call checkpoint(level=DL_INFOP, msg = 'file written '//filename)
    end subroutine


    ! ------------------------------------------------------------------
    !   I/O |  DUMPTXT    Write Real part of Complex Matrix to file
    ! ------------------------------------------------------------------
    !        Input   |   cmx        complex matrix
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
    subroutine CMatDumpTXT_real(cmx, filename)
        type(complex16_matrix), intent(in) :: cmx
        character(*), intent(in) :: filename
        integer, parameter :: fileunit = 22
        logical :: file_isopened
        integer :: ii

        open(fileunit, file=filename, action='write')
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then
            call checkpoint(level=DL_INFO, msg = 'writing CMat to '//filename)

            ! write file header:  size of matrix
            write(fileunit, *) '# [real8_matrix] %size:', cmx%size

            ! write to file
            do ii = 1, cmx%size(1)
                write (fileunit,"(*(G0,:,','))") real(cmx%val(:,ii))
            end do
        else
            call checkpoint(level=DL_ERROR, msg = 'file '//filename//' cannot be opened')
            return
        end if

        close(fileunit)
        call checkpoint(level=DL_INFOP, msg = 'file written '//filename)
    end subroutine

    ! ------------------------------------------------------------------
    !   I/O |  LOADTXT    Read Complex Matrix from file
    ! ------------------------------------------------------------------
    !        Input   |   filename   (char*) name of the input file
    !         
    !  Requirements  |   none, but check permission to read file
    ! 
    !        Output  |   cmx        the matrix parsed from file
    ! 
    !    Post-cond.  |   memory is allocated
    !  
    !  Err handlers  |   - file is opened (error)
    !                |   - file has correct header (warning)
    !                |   - unexpected end of file (fatal)
    ! ------------------------------------------------------------------
    function CMatLoadTXT(filename) result(cmx)
        character(*), intent(in) :: filename
        type(complex16_matrix) :: cmx
        integer, parameter :: fileunit = 22
        logical :: file_isopened
        integer :: ii, ios
        character(len=300) :: line  ! buffering line

        open(fileunit, file=filename, action='read')
        inquire(unit=fileunit, opened=file_isopened) 

        if(file_isopened) then
            call checkpoint(level=DL_INFO, msg = 'reading CMat from '//filename)

            ! parse the matrix size (must be in first line, after the delimiter "%size:")
            read(fileunit, '(A)', iostat=ios) line
            ii = index(line, "%size:")
            read (line( (ii+len("%size:")):), *) cmx%size
            write(line, *) cmx%size
            call checkpoint(level=DL_INFO, msg = 'parsed matrix of size ->'//line)
            
            call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
            allocate( cmx%val(cmx%size(1),cmx%size(2)) )

            call checkpoint(level=DL_INFO, msg = 'reading values')
            
            stop 'incomplete function' ! TODO fix this function

            do ii = 1, cmx%size(1)
                read(fileunit, *, iostat=ios) cmx%val(:,ii)
                if (ios /= 0) then
                    call checkpoint(level=DL_FATAL, msg = 'unexpected end of file')
                    print *, ii
                    call CMatPrintHead(cmx, 3, 3)
                    stop '(fatal error)'
                end if
            end do

            call checkpoint(level=DL_INFOP, msg = 'CMat reading completed')
        else
            call checkpoint(level=DL_ERROR, msg = 'file '//filename//' cannot be opened')
            return
        end if

        close(fileunit)
    end function


    ! ------------------------------------------------------------------
    !   I/O |  Print Head    Print the first nrow and ncol of a matrix
    ! ------------------------------------------------------------------
    !        Input   |   cmx        complex matrix
    !                |
    !                |   nrow   (integer) index of last row to print
    !                |
    !                !   ncol   (integer) index of lacomplex16_matrixst column to print
    !         
    !  Requirements  |   cmx must be initialized
    ! 
    !        Output  |   none, just print in stdout
    ! 
    !    Post-cond.  |   none
    !  
    !  Err handlers  |   none
    ! ------------------------------------------------------------------
    subroutine CMatPrintHead(cmx, nrow, ncol)
        type(complex16_matrix), intent(in) :: cmx
        integer :: ii, nrow, ncol

        do ii = 1,nrow
            print *, cmx%val(ii,:ncol)
        end do
    end subroutine


    
end module matrix_c16
