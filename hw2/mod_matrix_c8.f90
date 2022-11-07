! ========================================
!  QIC  >>  assignment 2 / exr 3
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 28 October 2022
! ========================================


module mod_matrix_c8
!  DOC =============================================================
! 
!    This module implements a matrix with complex*8 numbers.
!
!  -----------------------------------------------------------------
!
!  TYPES            |    complex8_matrix  
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
    use mydebugger  ! use my checkpoints

    type complex8_matrix
        ! store the dimension of the matrix
        integer, dimension(2) :: size
        ! to store the values of matrix
        complex*8, dimension(:,:), allocatable :: val
    end type

    interface operator(.randInit.)
        module procedure CMatInitRand
    end interface

    interface operator(.zeroInit.)
        module procedure CMatInitZero
    end interface

    interface operator(.Adj.)
        module procedure CMatAdjoint
    end interface

    interface operator(.Tr.)
        module procedure CMatTrace
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
    function CMatInitRand( shape ) result(cmx)
        integer, dimension(2), intent(in) :: shape
        type(complex8_matrix) :: cmx
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
    !   INITIALIZER : init a Complex Matrix to  (0 + i0)
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
    function CMatInitZero( shape ) result(cmx)
        integer, dimension(2), intent(in) :: shape
        type(complex8_matrix) :: cmx
        integer :: ii, jj

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        allocate( cmx%val(shape(1), shape(2)) )
        cmx%size = shape

        do jj=1,cmx%size(2)
            do ii=1,cmx%size(1)
                cmx%val(jj,ii) = complex( 0d0, 0d0)
            end do
        end do
    end function

    ! ------------------------------------------------------------------
    !   DESTRUCTOR : unallocate Complex Matrix
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
    subroutine CMatDelete(cmx)
        type(complex8_matrix) :: cmx

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
    function CMatAdjoint(cmx) result(cmxadj)
        type(complex8_matrix), intent(in) :: cmx
        type(complex8_matrix) :: cmxadj

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
    !        Output  |   tr     (complex*8)   trace of cmx
    ! 
    !    Post-cond.  |   none
    !  
    !  Err handlers  |   check if cmx is square (error)
    ! ------------------------------------------------------------------
    function CMatTrace(cmx) result(tr)
        type(complex8_matrix), intent(in) :: cmx
        complex*8 :: tr
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
        type(complex8_matrix), intent(in) :: cmx
        character(*), intent(in) :: filename
        integer, parameter :: fileunit = 22
        logical :: file_isopened
        integer :: ii

        open(fileunit, file=filename, action='write')
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then
            call checkpoint(level=DL_INFO, msg = 'writing CMat to '//filename)

            ! write file header:  size of matrix
            write(fileunit, *) '# [complex8_matrix] %size:', cmx%size
            write(fileunit, *) '# %val:'

            ! write to file
            do ii = 1, cmx%size(1)
                write(fileunit, *) cmx%val(ii, :)
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
        type(complex8_matrix) :: cmx
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

            ! read & skip second line, where I expect to find "%val:"
            read(fileunit, '(A)', iostat=ios) line
            ii = index(line, "%val:")
            if ( ii .eq. 0) then
                call checkpoint(level=DL_WARN, msg = 'delimiter %val not found')
            end if

            call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
            allocate( cmx%val(cmx%size(1),cmx%size(2)) )

            call checkpoint(level=DL_INFO, msg = 'reading values')
            do ii = 1, cmx%size(1)
                read(fileunit, *, iostat=ios) cmx%val(ii, :)
                if (ios /= 0) then
                    call checkpoint(level=DL_FATAL, msg = 'unexpected end of file')
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
    !                !   ncol   (integer) index of last column to print
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
        type(complex8_matrix), intent(in) :: cmx
        integer :: ii

        do ii = 1,nrow
            print *, cmx%val(ii,:ncol)
        end do
    end subroutine
    
end module mod_matrix_c8