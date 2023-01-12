! ========================================
!  QIC  >>  assignment 8
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 17 December 2022
! ========================================



module rg_algo

    ! requirements:
    use mydebugger  ! use my checkpoints
    use matrix_r8

    ! sys
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN

    implicit none


    ! -------------------------------------------
    type rga
        integer*8 :: N
        double precision :: thres = 1e-10
        integer :: max_iter = 1000
        integer :: counter
        logical :: converged
        type(real8_matrix) :: H, A, B
        integer :: energy_level = 1  ! select eigenvector to ref

        double precision, dimension(:), allocatable :: hist
        double precision :: cputime = 0
    end type




contains




    ! ------------------------------------------------------------------
    !   INITIALIZER |  init a RG algorithm type
    ! ------------------------------------------------------------------
    function rg_init(N, thres, max_iter) result(this)
        type(rga) :: this
        integer*8 :: N
        double precision, optional :: thres
        integer, optional :: max_iter

        this%N = N 

        if( present(thres) ) then
            this%thres = thres
        end if
        if( present(max_iter) ) then
            this%max_iter = max_iter
        end if

        this%counter = 1
        this%converged = .false.
    end function




    ! ------------------------------------------------------------------
    !   RG ITERATION    |  step of RG algorithm
    ! ------------------------------------------------------------------
    subroutine rg_exe(this)
        type(rga) :: this
        type(real8_matrix) :: H2N, tmp, H2N_copy
        type(real8_matrix) :: eigvec
        double precision, dimension(:), allocatable :: eigval
        double precision :: energy, energy_old

        double precision :: time_tic, time_toc

        if( this%converged ) then
            print *, ' [war] this RG already converged!'
        end if

        ! init the history keeper
        allocate( this%hist(this%max_iter) )
        this%hist = IEEE_VALUE(1d0, IEEE_QUIET_NAN)
        
        ! algo param init
        energy_old = huge(energy_old)

        print *, '+++ EXE SIMULATION +++'
        call CPU_TIME(time_tic)
        do while( .not.this%converged ) ! not referenced, but requested
            print *, 'iteration', this%counter
    
            ! build 2N system
            H2N = kron(this%A, this%B)
            call kron_r1_inplace(this%H, 2_8**(this%N), H2N)
            call kron_l1_inplace(2_8**(this%N), this%H, H2N)
    
    
            ! compute eigenvv
            H2N_copy = copy_rmat(H2N)   ! copy H2N, because it gets destroyed
            call eigenvv_symmetric( int(2_8**(this%N)), H2N_copy, eigval, eigvec)
            call Delete(H2N_copy)

            ! evaluate ground state ...
            energy = eigval( this%energy_level )/(2_8*(this%N))
            this%hist(this%counter) = energy  ! keep track of energy at each step

            if( abs( energy - energy_old ) < this%thres ) then        ! exit criteria: below abs thres
                this%converged = .true.   ! exit condition 
                print *, ' [exit] converged at iteration', this%counter

            else if( this%counter.gt.this%max_iter ) then    ! exit criteria: above allowed iterations
                this%converged = .true.   ! exit condition
                this%counter = -1         ! mark exit by iteration count
                print *, ' [exit] reached max iter'
            
            ! ... and getting ready for next iteration
            else
                ! project H
                this%H%data = matmul(  matmul( transpose(eigvec%data), H2N%data), eigvec%data ) 
                call Delete(H2N)
        
                ! project A
                tmp = kron_l1(2_8**(this%N), this%A)
                this%A%data = matmul(  matmul( transpose(eigvec%data), tmp%data), eigvec%data ) 
                call Delete(tmp)
    
                ! project B
                tmp = kron_r1(this%B, 2_8**(this%N))
                this%B%data = matmul(  matmul( transpose(eigvec%data), tmp%data), eigvec%data ) 
                call Delete(tmp)
    
                ! getting ready for next iteration
                this%counter = this%counter + 1
                energy_old = energy
                this%H%data = this%H%data / 2.0d0
                this%A%data = this%A%data / sqrt(2.0d0)
                this%B%data = this%B%data / sqrt(2.0d0)
            end if
    
            call Delete(eigvec)
            deallocate(eigval)
    
        end do
        call CPU_TIME(time_toc)
        this%cputime = time_toc-time_tic
    end subroutine


    ! ------------------------------------------------------------------
    !   I/O   |  write to file the parameters of RG object
    ! ------------------------------------------------------------------
    subroutine write_header(this, filename)
        type(rga) :: this
        character(*), intent(in) :: filename

        integer, parameter :: fileunit = 22
        logical :: file_isopened

        open(fileunit, file=filename, action='write', form="unformatted")
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then
            call checkpoint(level=DL_INFO, msg = 'writing RG simulation to (binary)'//filename)

            write(fileunit) this%N
            write(fileunit) this%thres
            write(fileunit) this%max_iter

            close(fileunit)
            call checkpoint(level=DL_INFOP, msg = 'file written '//filename)
        else
            call checkpoint(level=DL_ERROR, msg = 'file '//filename//' cannot be opened')
            return
        end if
    end subroutine



    ! ------------------------------------------------------------------
    !   I/O   |  write to file the results of RG run
    ! ------------------------------------------------------------------
    subroutine write_entry(this, par, filename)
        type(rga) :: this
        real*8 :: par  ! a value to be stored along this entry

        character(*), intent(in) :: filename

        integer, parameter :: fileunit = 22
        logical :: file_isopened

        open(fileunit, file=filename, status="old", action="write", &
                       form='unformatted', position='append')
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then
            call checkpoint(level=DL_INFO, msg = 'writing RG simulation entry to (binary)'//filename)

            write(fileunit) par
            write(fileunit) this%counter
            write(fileunit) this%hist
            write(fileunit) this%cputime
            
            close(fileunit)
            call checkpoint(level=DL_INFOP, msg = 'file written '//filename)
        else
            call checkpoint(level=DL_ERROR, msg = 'file '//filename//' cannot be opened')
            return
        end if
    end subroutine


end module rg_algo