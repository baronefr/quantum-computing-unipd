! ========================================
!  QIC  >>  assignment 6
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 9 December 2022
! ========================================




module density_matrix_module
    implicit none

    type :: density_matrix
        integer :: dim, nn, size
        complex(8), dimension(:,:), allocatable :: data
    end type

contains

    ! ------------------------------------------------------------------
    !   INITIALIZER : allocate a density matrix
    ! ------------------------------------------------------------------
    !        Input   |   dim      dimensions of a single quantum system
    !                |   nn       number of elements in the system
    !         
    !  Requirements  |   (be careful to not overwrite allocated objects
    !                |    with expression assignment)
    ! 
    !        Output  |   dm        the initialized densitymatrix
    ! 
    !    Post-cond.  |   dm is a square complex matrix with side dim**nn 
    !  
    !  Err handlers  |   none
    ! ------------------------------------------------------------------
    function InitDM(dim, nn) result(dm)
        integer, intent(in) :: dim, nn
        type(density_matrix) :: dm

        dm%dim = dim
        dm%nn = nn
        dm%size = dm%dim**dm%nn ! compute number of elements
        allocate( dm%data( dm%size, dm%size ) )
    end function


    ! ------------------------------------------------------------------
    !   INITIALIZER : generate a density matrix from a state vector psi
    ! ------------------------------------------------------------------
    function dm_from_state(psi, dim, nn) result(dm)
        complex(8), dimension(:) :: psi
        integer, intent(in) :: dim, nn
        type(density_matrix) :: dm

        complex(8), dimension(:,:), allocatable :: bra, ket

        dm = InitDM(dim, nn)

        allocate( bra(1,size(psi,1)), ket(size(psi,1),1) )
        ! to use matmul, I need these 2D objects:
        ket(:,1) = psi        ! as column
        bra(1,:) = CONJG(psi) ! as row

        dm%data = matmul( ket, bra )

        deallocate(bra, ket)
    end function

    

    ! ------------------------------------------------------------------
    !   PARTIAL TRACE |  given a system [ A, B ]
    !                 |  returns the density matrix over B
    ! ------------------------------------------------------------------
    subroutine ptr_remove_left(dm, nsys, trdm)
        type(density_matrix) :: dm    ! input density matrix
        type(density_matrix) :: trdm  ! traced density matrix
        integer :: nsys               ! system number to remove

        integer :: ii, jj, kk, lss, rss

        ! TODO check nsys > 0

        trdm = InitDM(dm%dim, dm%nn-nsys)
        lss = (dm%dim)**(dm%nn-nsys)
        rss = (dm%dim)**nsys

        trdm%data = (0.d0,0.d0)
        do ii = 1, lss
            do jj = 1, lss
                do kk = 1, rss
                    trdm%data(ii,jj) = trdm%data(ii,jj) + dm%data( rss*(ii-1) + kk, rss*(jj-1) + kk)
                enddo
            enddo
        enddo

        ! returns  [ -  B ]
    end subroutine



    ! ------------------------------------------------------------------
    !   PARTIAL TRACE |  given a system [ A, B ]
    !                 |  returns the density matrix over A
    ! ------------------------------------------------------------------
    subroutine ptr_remove_right(dm, nsys, trdm)
        type(density_matrix) :: dm    ! input density matrix
        type(density_matrix) :: trdm  ! traced density matrix
        integer :: nsys               ! system number to remove

        integer :: ii, jj, kk, lss, rss

        ! TODO check nsys > 0

        trdm = InitDM(dm%dim, dm%nn-nsys)
        lss = (dm%dim)**nsys
        rss = (dm%dim)**(dm%nn-nsys)

        trdm%data = (0.d0,0.d0)
        do ii = 1, rss
            do jj = 1, rss
                do kk = 1, lss
                    trdm%data(ii,jj) = trdm%data(ii,jj) + dm%data( rss*(kk-1) + ii, rss*(kk-1) + jj )
                enddo
            enddo
        enddo

        ! returns  [ A  - ]
    end subroutine





    ! ------------------------------------------------------------------
    !   PARTIAL TRACE |  given the density matrix dm, over n systems of
    !                 |  dimension D each,
    !                 |  returns the density matrix over system nsys
    ! ------------------------------------------------------------------
    function partial_trace(dm, system) result(trdm)
        type(density_matrix) :: dm, trdm
        integer :: system
        type(density_matrix) :: buffer
        logical :: two_steps

        ! TODO check systems has at least one element

        if( (1<system) .and. (system < dm%nn) ) then
            two_steps = .true.
        else
            two_steps = .false.
        end if

        if( system > 1 ) then
            if(two_steps) then
                ! place in buffer density matrix if RxTrace is planned ...
                call ptr_remove_left(dm, system-1, buffer)
            else
                call ptr_remove_left(dm, system-1, trdm)
            end if
        end if

        if( system < dm%nn ) then
            if(two_steps) then
                ! ... and ref buffer here
                call ptr_remove_right(buffer, dm%nn-system, trdm)
                deallocate(buffer%data)
            else
                call ptr_remove_right(dm, dm%nn-system, trdm)
            end if
        end if

    end function





    ! ------------------------------------------------------------------
    !   I/O |  DUMP    Read density matrix from binary file
    ! ------------------------------------------------------------------
    function dmread(filename) result(dm)
        type(density_matrix) :: dm
        character(*), intent(in) :: filename
        integer, parameter :: fileunit = 22
        logical :: file_isopened

        open(fileunit, file=filename, action='read', form="unformatted") !, access='stream')
        
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then

            ! manual init procedure .....
            read(fileunit) dm%dim
            read(fileunit) dm%nn
            dm%size = dm%dim**dm%nn ! compute number of elements
            print "(AI0AI0A)", 'reading [ dim = ', dm%dim, ', nn = ', dm%nn, ' ] density matrix'

            allocate( dm%data( dm%size, dm%size ) )
            read(fileunit) dm%data

            close(fileunit)
        else
            print *, 'ERR cannot open file'
            return
        end if

    end function


    ! ------------------------------------------------------------------
    !   I/O |  DUMP    Write density matrix to binary file
    ! ------------------------------------------------------------------
    subroutine dmwrite(dm, filename)
        type(density_matrix), intent(in) :: dm
        character(*), intent(in) :: filename
        integer, parameter :: fileunit = 22
        logical :: file_isopened

        open(fileunit, file=filename, action='write', form="unformatted")
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then
            write(fileunit) dm%dim
            write(fileunit) dm%nn
            write(fileunit) dm%data
            close(fileunit)
        else
            print *, 'ERR cannot write file'
            return
        end if        
    end subroutine

    ! ------------------------------------------------------------------
    !   I/O |  PRINT   print a density matrix in compact form
    ! ------------------------------------------------------------------
    subroutine dmprint(dm, nrow, ncol)
        type(density_matrix), intent(in) :: dm
        integer :: ii, nrow, ncol
        do ii = 1,nrow
            write(*, "(*('('sf6.2xspf6.2x'i)':x))") dm%data(ii,:ncol)
        end do
    end subroutine

end module


