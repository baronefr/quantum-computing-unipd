! ========================================
!  Quantum Information and Computing
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   last update : 15 December 2022
! ========================================

module matrix_r8_pack
!  DOC =============================================================
! 
!    This module implements a symmetric matrix with real*8 
!    (double precision) numbers in packed mode.
! 
!  -----------------------------------------------------------------
! 
!  TYPES            |    real8_pack_matrix  
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
    type real8_pack_matrix
        ! store the dimension of the matrix
        integer*8 :: size
        ! to store the values of matrix
        real(kind=kkind), dimension(:), allocatable :: data

        character(1) :: uplo
    end type

    interface operator(.zeroInit.)
        module procedure InitZerop
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
    function InitZerop( shape ) result(rmx)
        integer*8, intent(in) :: shape
        type(real8_pack_matrix) :: rmx

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        ! init properties of matrix
        allocate( rmx%data( (shape*(shape+1)/2) ) )
        rmx%size = shape

        ! writing values
        rmx%data = 0d0
        rmx%uplo = 'U'     ! A(i,j) = AP(i + (j-1)*j/2)   for j =< i
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
    subroutine deletep(rmx)
        type(real8_pack_matrix) :: rmx

        if( allocated(rmx%data) ) then
            call checkpoint(level=DL_PEDANTIC, msg = 'deallocating memory')
            deallocate( rmx%data )
            rmx%size = 0;
        else
            call checkpoint(level=DL_FATAL, msg = 'cannot deallocate prmat, already cleared')
            stop '(fatal error)'
        end if
    end subroutine





    ! ------------------------------------------------------------------
    !   MATH |  TENSOR PRODUCT
    ! ------------------------------------------------------------------
    function kronp_bak(rmxx, rmyy) result(rmo)
        type(real8_pack_matrix) :: rmxx, rmyy
        type(real8_pack_matrix) :: rmo
        integer*8 :: ii, jj, aa, bb, oii, ojj, broad_i, stream_i, out_i

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        rmo = InitZerop( rmxx%size*rmyy%size )

        do jj = 1_8, rmxx%size
            do ii = 1_8, rmxx%size

                if( ii.gt.jj) then
                    stream_i = jj + (ii-1_8)*ii/2_8
                else
                    stream_i = ii + (jj-1_8)*jj/2_8 
                end if

                do bb = 1_8, rmyy%size
                    do aa = 1_8, rmyy%size

                        ! fetch output indexes
                        oii = (ii-1_8)*rmyy%size + aa
                        ojj = (jj-1_8)*rmyy%size + bb

                        ! fetch broadcast & stream indexes
                        if( aa.gt.bb) then
                            broad_i = bb + ((aa-1_8)*aa/2_8)
                        else
                            broad_i = aa + ((bb-1_8)*bb/2_8)
                        end if

                        if( oii.le.ojj) then
                            out_i = oii + (ojj-1_8)*ojj/2_8
                            rmo%data( out_i ) = rmxx%data( stream_i ) * rmyy%data( broad_i )
                        end if
                        
                    end do
                end do
            end do
        end do
    end function


    ! ------------------------------------------------------------------
    !   MATH |  TENSOR PRODUCT
    ! ------------------------------------------------------------------
    function kronp(rmxx, rmyy) result(rmo)
        type(real8_pack_matrix) :: rmxx, rmyy
        type(real8_pack_matrix) :: rmo
        integer*8 :: ii, jj, aa, bb, oii, ojj, broad_i, stream_i, out_i

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        rmo = InitZerop( rmxx%size*rmyy%size )

        do ii = 1_8, rmxx%size
        do jj = 1_8, rmxx%size

            if( jj.le.ii ) then
                stream_i = jj + ((ii-1_8)*ii/2_8)
            else
                stream_i = ii + ((jj-1_8)*jj/2_8)
            end if
                
            do aa = 1_8, rmyy%size
            do bb = 1_8, rmyy%size

                oii = (ii-1_8)*rmyy%size + aa
                ojj = (jj-1_8)*rmyy%size + bb

                if( bb.le.aa ) then
                    broad_i = bb + ((aa-1_8)*aa/2_8)
                else
                    broad_i = aa + ((bb-1_8)*bb/2_8)
                end if

                if( ojj.le.oii ) then
                    rmo%data( ojj + (oii*(oii-1_8)/2_8) ) = rmxx%data( stream_i ) * rmyy%data( broad_i )
                else
                    rmo%data( oii + (ojj*(ojj-1_8)/2_8)) = rmxx%data( stream_i ) * rmyy%data( broad_i )
                end if

            end do
            end do
        end do
        end do
    end function


    ! ------------------------------------------------------------------
    !   MATH |  TENSOR PRODUCT with left identity
    ! ------------------------------------------------------------------
    function kronp_l1(idsize, rmyy) result(rmo)
        implicit none
        integer*8, intent(in) :: idsize
        type(real8_pack_matrix) :: rmyy
        type(real8_pack_matrix) :: rmo
        integer*8 :: ii, jj, aa, bb, oii, ojj, broad_i, stream_i, out_i

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        rmo = InitZerop( idsize*rmyy%size )

        do ii = 1_8, idsize
        jj = ii
    
            stream_i = jj + ((ii-1_8)*ii/2_8)
                    
            do aa = 1_8, rmyy%size
            do bb = 1_8, rmyy%size
    
                oii = (ii-1_8)*rmyy%size + aa
                ojj = (jj-1_8)*rmyy%size + bb
    
                if( bb.le.aa ) then
                    broad_i = bb + ((aa-1_8)*aa/2_8)
                else
                    broad_i = aa + ((bb-1_8)*bb/2_8)
                end if
    
                if( ojj.le.oii ) then
                    rmo%data( ojj + (oii*(oii-1_8)/2_8) ) = rmyy%data( broad_i )
                else
                    !rmo%data( oii + (ojj*(ojj-1_8)/2_8) ) = rmyy%data( broad_i )
                end if
    
            end do
            end do
            
        end do
    end function

    ! ------------------------------------------------------------------
    !   MATH |  TENSOR PRODUCT with left identity
    ! ------------------------------------------------------------------
    subroutine kronp_l1_inplace(idsize, rmyy, rmo)
        implicit none
        integer*8, intent(in) :: idsize
        type(real8_pack_matrix) :: rmyy
        type(real8_pack_matrix), intent(inout) :: rmo
        integer*8 :: ii, jj, aa, bb, oii, ojj, broad_i, stream_i, out_i

        do ii = 1_8, idsize
        jj = ii
        
                stream_i = jj + ((ii-1_8)*ii/2_8)
                        
                do aa = 1_8, rmyy%size
                do bb = 1_8, rmyy%size
        
                    oii = (ii-1_8)*rmyy%size + aa
                    ojj = (jj-1_8)*rmyy%size + bb
        
                    if( bb.le.aa ) then
                        broad_i = bb + ((aa-1_8)*aa/2_8)
                    else
                        broad_i = aa + ((bb-1_8)*bb/2_8)
                    end if
        
                    if( ojj.le.oii ) then
                        rmo%data( ojj + (oii*(oii-1_8)/2_8) ) = rmo%data( ojj + (oii*(oii-1_8)/2_8) ) + rmyy%data( broad_i )
                    else
                        !rmo%data( oii + (ojj*(ojj-1_8)/2_8) ) = rmo%data( oii + (ojj*(ojj-1_8)/2_8) ) + rmyy%data( broad_i )
                    end if
        
                end do
                end do
                
        end do
    end subroutine

    ! ------------------------------------------------------------------
    !   MATH |  TENSOR PRODUCT with right identity
    ! ------------------------------------------------------------------
    function kronp_r1(rmxx, idsize) result(rmo)
        type(real8_pack_matrix) :: rmxx
        integer*8, intent(in) :: idsize
        type(real8_pack_matrix) :: rmo
        integer*8 :: ii, jj, aa, bb, oii, ojj, broad_i, stream_i, out_i

        call checkpoint(level=DL_PEDANTIC, msg = 'allocating memory')
        rmo = InitZerop( rmxx%size*idsize )

        do ii = 1_8, rmxx%size
        do jj = 1_8, rmxx%size
    
            if( jj.le.ii ) then
                stream_i = jj + ((ii-1_8)*ii/2_8)
            else
                stream_i = ii + ((jj-1_8)*jj/2_8)
            end if
                    
            do aa = 1_8, idsize
            bb = aa
    
                oii = (ii-1_8)*idsize + aa
                ojj = (jj-1_8)*idsize + bb
    
                broad_i = bb + ((aa-1_8)*aa/2_8)

    
                if( ojj.le.oii ) then
                    rmo%data( ojj + (oii*(oii-1_8)/2_8) ) = rmxx%data( stream_i )
                else
                    !rmo%data( oii + (ojj*(ojj-1_8)/2_8)) = rmxx%data( stream_i )
                end if
    
            end do
        end do
        end do
    end function

    ! ------------------------------------------------------------------
    !   MATH |  TENSOR PRODUCT with right identity
    ! ------------------------------------------------------------------
    subroutine kronp_r1_inplace(rmxx, idsize, rmo)
        type(real8_pack_matrix) :: rmxx
        integer*8, intent(in) :: idsize
        type(real8_pack_matrix), intent(inout) :: rmo
        integer*8 :: ii, jj, aa, bb, oii, ojj, broad_i, stream_i, out_i


        do ii = 1_8, rmxx%size
        do jj = 1_8, rmxx%size
        
            if( jj.le.ii ) then
                stream_i = jj + ((ii-1_8)*ii/2_8)
            else
                stream_i = ii + ((jj-1_8)*jj/2_8)
            end if
                        
            do aa = 1_8, idsize
            bb = aa
        
                oii = (ii-1_8)*idsize + aa
                ojj = (jj-1_8)*idsize + bb
        
                broad_i = bb + ((aa-1_8)*aa/2_8)
    
        
                if( ojj.le.oii ) then
                    rmo%data( ojj + (oii*(oii-1_8)/2_8) ) = rmxx%data( stream_i ) + rmo%data( ojj + (oii*(oii-1_8)/2_8) )
                else
                    !rmo%data( oii + (ojj*(ojj-1_8)/2_8) ) = rmxx%data( stream_i ) + rmo%data( ojj + (oii*(oii-1_8)/2_8) )
                end if
        
            end do
        end do
        end do
    end subroutine





    ! ------------------------------------------------------------------
    !   MATH |  compute EIGENVALUES and EIGENVECTORS of packed matrix
    ! ------------------------------------------------------------------
    subroutine eigenvv_pack(k, rmat, w, z, do_eigvec)
        implicit none
        type(real8_pack_matrix), intent(in) :: rmat
        integer, intent(in) :: k
        double precision, intent(out), allocatable :: w(:)
        !type(real8_matrix), intent(inout) :: z
        logical, optional :: do_eigvec ! if true, compute eigenvectors


        ! locals
        double precision, parameter :: zero = 0d0
        double precision :: abstol, r, vl, vu
        integer :: i, ifail, il, info, iu, j, ldz, m, n
        !double precision, allocatable :: ap(:), w(:), work(:), z(:, :)
        double precision, allocatable :: work(:), z(:, :)
        integer, allocatable :: iwork(:), jfail(:)
        character(1) :: jobz

        if( present( do_eigvec ) ) then
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
        n = rmat%size
        ldz = n
        m = iu-il+1

        call checkpoint(level = DL_PEDANTIC, msg = 'allocate tmp for dspevx')
        if ( jobz == 'V') then
            allocate( z(ldz,m) )
        end if
        allocate( w(n), work(8*n), iwork(5*n), jfail(n) )

        !     Set the absolute error tolerance for eigenvalues.  With ABSTOL
        !     set to zero, the default value is used instead
        abstol = zero

        ! eigenvalue problem solver
        call checkpoint(level = DL_PEDANTIC, msg = 'dspevx (solver)')
        call dspevx(jobz, 'I', rmat%uplo, n, rmat%data, vl, vu, il, iu, &
                    abstol, m, w, z, ldz, work, iwork, jfail, info)

        ! [ ref ] ----------------- 
        ! https://netlib.org/lapack/explore-html/dc/dd2/group__double_o_t_h_e_reigen_ga07ee2c397b1b0f73e296f20f8d36990a.html


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
    subroutine PMatDump(rmx, filename)
        type(real8_pack_matrix), intent(in) :: rmx
        character(*), intent(in) :: filename
        integer, parameter :: fileunit = 22
        logical :: file_isopened
        integer*8 :: ii

        open(fileunit, file=filename, action='write', form="unformatted")
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then
            call checkpoint(level=DL_INFO, msg = 'writing PMat to (binary)'//filename)

            !write(fileunit) rmx%size
            !write(fileunit) rmx%uplo
            write(fileunit) rmx%data
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
    !                !   ncol   (integer) index of lareal8_pack_matrixst column to print
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
        type(real8_pack_matrix), intent(in) :: rmx
        integer*8, optional :: nrow, ncol
        integer*8 :: nrow_, ncol_
        integer*8 :: ii, jj

        if( present(nrow) ) then
            nrow_ = nrow
        else
            nrow_ = rmx%size
        end if

        if( present(ncol) ) then
            ncol_ = ncol
        else
            ncol_ = rmx%size
        end if

        do ii = 1, nrow_
            print "(*(F0.3,:,'  '))",  ( rmx%data(ii + (jj-1)*jj/2), jj=1,ncol_)
        end do

        ! compact format for complexes...
        !do ii = 1,nrow_
        !    write(*, "(*('('sf6.2xspf6.2x'i)':x))") rmx%data(ii,:ncol_)
        !end do
    end subroutine


    
end module matrix_r8_pack
