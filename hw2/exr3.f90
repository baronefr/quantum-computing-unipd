! ========================================
!  QIC  >>  assignment 2 / exr 3
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 29 October 2022
! ========================================

!
! This program tests the module of complex matrices.
! 

program main
    use mod_matrix_c8 ! matrix of complex*8 numbers
    type(complex8_matrix) :: A, B, C, D
    complex*8 :: x

    integer, dimension(2) :: nice_shape = (/100, 150/)

    ! set an higher verbosity
    DEBUG_LEVEL = DL_PEDANTIC

    ! ---------------------------------
    !   INIT
    ! ---------------------------------
    A = CMatInitZero(shape=nice_shape)
    B = .randInit.(/100, 100/)  ! square matrix

    ! ---------------------------------
    !   MATH OPERATIONS
    ! ---------------------------------
    ! testing Trace (square matrix only)
    x = .Tr.B
    print *, "Trace of B is", x
    call CMatPrintHead(B, 3, 3)

    ! testing Adjoint
    C = .Adj.B
    
    ! ---------------------------------
    !   I/O
    ! ---------------------------------

    ! dumping a matrix to file
    call CMatDumpTXT(C, 'data/matrix.txt')
    print *, "The matrix written to file"
    call CMatPrintHead(C, 4, 4)

    ! loading the same matrix from file
    D = CMatLoadTXT('data/matrix.txt')
    print *, "The matrix taken from file"
    call CMatPrintHead(D, 4, 4)



    !! TEST1 : Check if matrix loaded from file is
    !          the same matrix already stored in memory.
    if(  all( C%val .eq. D%val ) ) then
        print *, "The matrices are the same ..."
    else
        call checkpoint(level=DL_ERROR, msg='matrices are NOT equal')
    end if



    !! TEST2: Check if an altered matrix is not equal to
    !         the one stored in memory.
    print *, " ... but now I change an element"
    D%val(1,1) = complex(1d0,0d0)
    if(  all( C%val .eq. D%val ) ) then
        print *, "The matrices are the same!"
    else
        call checkpoint(level=DL_ERROR, msg='matrices are NOT equal')
    end if

    
    ! ---------------------------------
    !   DEALLOCATING
    ! ---------------------------------
    call CMatDelete(A)
    call CMatDelete(B)
    !call CMatDelete(A)    ! to test the fatal error
end program