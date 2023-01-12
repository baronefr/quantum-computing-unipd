! ========================================
!  QIC  >>  assignment 8
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 23 December 2022
! ========================================





module idmrg_algo

    ! requirements:
    use mydebugger  ! use my checkpoints
    use matrix_r8

    ! sys
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN

    implicit none

    ! -------------------------------------------

    type block
        type(real8_matrix) :: HH, Hx, Hz
        integer :: mm  ! size of block
    end type

    type idmrga
        integer :: max_m = 4        ! max size of block
        integer :: max_iter = 1000  ! max iterations allowed

        integer :: counter = 1
        logical :: converged = .false.
        logical :: save_hist = .true.
        double precision, dimension(:), allocatable :: hist, hist_err

        real*8 :: lambda
        double precision :: energy
        integer :: energy_level = 1  ! select eigenvector to ref
        
        type(block) :: lblock
        ! no right block, I am using lblock itself

        double precision :: cputime = 0
    end type



contains




    ! ------------------------------------------------------------------
    !   INITIALIZER |  init a DMRG algorithm type
    ! ------------------------------------------------------------------
    function idmrg_init(lambda, max_diag, max_m, max_iter, save_hist) result(this)
        type(idmrga) :: this
        real*8 :: lambda
        integer, optional :: max_diag, max_m, max_iter
        logical, optional :: save_hist
        
        ! dummy        
        type(real8_matrix) :: sigmax

        ! optional args:  max_diag ------------------------
        !    max_diag   |  The maximum size of matrix to be diagonalized
        !               |  at each step. Must be even.
        !               |  If present, then max_m = (max_diag - 2)/2 .
        !
        !    max_x      |  The maximum size of enlarged block of iDMRG.
        !               |  If max_diag and max_m are not present, then
        !               |  max_m = 4 by default.
        if( present(max_diag) ) then
            ! checking max_diag to be even
            if( mod(max_diag,2) .ne. 0 ) then
                print*, 'max_diag is odd, rounding up '
                max_diag = max_diag + 1
            end if
            this%max_m = (max_diag - 2)/2
            call checkpoint(level = DL_INFO, msg = 'set max_m by max_diag', var= this%max_m)

        else if( present(max_m) ) then
            this%max_m = max_m
            call checkpoint(level = DL_INFO, msg = 'set max_m by arg', var= this%max_m)

        else
            !this%max_m = 4  ! use default type init values
            call checkpoint(level = DL_INFO, msg = 'set max_m to default', var= this%max_m)
        end if

        ! optional args:  max_iter ------------------------
        !    max_iter     |   Max number of iterations for DMRG object.
        !                 |   If max_iter < 0, then no limits are enforced
        !                 |   on DMRG iterations.
        if( present(max_iter) ) then
            this%max_iter = max_iter
            call checkpoint(level = DL_INFO, msg = 'set max_iter by arg', var= this%max_iter)
        else
            !this%max_iter = 1000 ! use default type init values
            call checkpoint(level = DL_INFO, msg = 'set max_iter to default', var= this%max_iter)
        end if

        ! optional args:  save_hist -----------------------
        !    save_hist    |   If true, keeps track of energy through iterations.
        !                 !   History is not saved anyway if  max_iter <= 0 .
        if( present(save_hist) ) then
            this%save_hist = save_hist
        end if

        if( this%max_iter .le. 0 ) then
            this%save_hist = .false.
        end if

        ! parameters of DMRG object -----------------------
        this%counter = 1
        this%converged = .false.

        ! history object ----------------------------------
        if( this%save_hist ) then
            call checkpoint(level = DL_INFO, msg = 'hist -> on')
            allocate( this%hist(this%max_iter) )
            allocate( this%hist_err(this%max_iter) )
        else
            call checkpoint(level = DL_INFO, msg = 'hist -> off')
            allocate( this%hist( 1 ) )
            allocate( this%hist_err( 1 ) )
        end if
        this%hist = IEEE_VALUE(1d0, IEEE_QUIET_NAN) ! init array to NaN


        ! system properties & init ------------------------
        this%lambda = lambda
        this%lblock%mm = 1

        ! init an Hamiltonian for one site
        this%lblock%HH = InitZero( (/2_8, 2_8/) )
        this%lblock%HH%data(1,1) = lambda;    this%lblock%HH%data(2,2) = -lambda;

        ! interaction terms
        this%lblock%Hx = InitZero( (/2_8, 2_8/) )
        this%lblock%Hx%data(2,1) = 1d0;       this%lblock%Hx%data(1,2) = 1d0;

        ! spin term of site m + 1 (referenced when creating a superblock)
        this%lblock%Hz = InitZero( (/2_8, 2_8/) )
        this%lblock%Hz%data(1,1) = 1d0;       this%lblock%Hz%data(2,2) = 1d0;

        call checkpoint(level=DL_PEDANTIC, msg = 'DMRG initialized')
    end function



    ! ------------------------------------------------------------------
    !   DESTRUCTOR |  clear the allocated memory of a DMRG object
    ! ------------------------------------------------------------------
    subroutine idmrg_delete(this)
        type(idmrga) :: this

        deallocate( this%hist )
        call clear_block( this%lblock )

    end subroutine



    ! ------------------------------------------------------------------
    !   ENLARGE   ------------------------------------------------------
    !
    !   Add a site to the left of a block of current size m
    ! ------------------------------------------------------------------
    function enlarge(this) result(enl)
        type(idmrga) :: this
        type(block)  :: enl
        type(real8_matrix) :: sigmax, sigmaz

        ! some requirements
        sigmax = InitZero( (/2_8, 2_8/) )
        sigmax%data(2,1) = 1d0;       sigmax%data(1,2) = 1d0;

        sigmaz = InitZero( (/2_8, 2_8/) )
        sigmaz%data(1,1) = this%lambda;    sigmaz%data(2,2) = -this%lambda;

        ! enlarge H  ->  H x 1(d)
        enl%HH = kron_r1(this%lblock%HH, 2_8)
        call kron_inplace(this%lblock%Hx, sigmax, enl%HH)
        call kron_inplace(this%lblock%Hz, sigmaz, enl%HH)

        ! enlarge interaction Hamiltonian
        enl%Hx = kron_l1( 2_8**(this%lblock%mm), sigmax )
        enl%Hz = InitIdentity( 2_8**(this%lblock%mm + 1) )

        ! enlarge current block size
        enl%mm = this%lblock%mm + 1

        call Delete(sigmax)
        call Delete(sigmaz)
    end function


    ! ------------------------------------------------------------------
    !   CLEAR BLOCK  |  Deallocate data in block object.
    ! ------------------------------------------------------------------
    subroutine clear_block(blk)
        type(block) :: blk

        call Delete( blk%HH )
        call Delete( blk%Hx )
        call Delete( blk%Hz )
    end subroutine



    ! ------------------------------------------------------------------
    !   BUILD SUPERBLOCK HAMILTONIAN  ----------------------------------
    !
    !   Creates the full Hamiltonian given an enlarged block <enl>.
    ! ------------------------------------------------------------------
    function superblock_hamiltonian(enl) result(H)
        type(block) :: enl
        type(real8_matrix) :: H !, middle, tmp, sigmax

        ! allocate L/R block spin
        H = kron_r1( enl%HH, 2_8**(enl%mm) )
        call kron_l1_inplace(2_8**(enl%mm), enl%HH, H)

        call kron_inplace(enl%Hx, enl%Hx, H)
        !call kron_inplace(enl%Hz, enl%Hz, H)
    end function




    ! ------------------------------------------------------------------
    !   DMRG ITERATION  |  step of DMRG algorithm
    ! ------------------------------------------------------------------
    subroutine idmrg_exe(this)
        type(idmrga) :: this

        type(real8_matrix) :: H, rho, eigvec

        double precision, dimension(:), allocatable :: eigval
        double precision :: energy
        type(block) :: enlarged_block

        double precision, dimension(:, :), allocatable :: proj, proj_dag

        double precision :: eigsum, eigtmp
        double precision :: time_tic, time_toc
        type(real8_matrix) :: sigmax
        integer :: ii, crop_eigenvalues

        if( this%converged ) then
            call checkpoint(level=DL_ERROR, msg = 'DMRG marked as converged, skip iter')
            return
        end if


        ! tic the chronometer (performance monitor) -------
        call CPU_TIME(time_tic)



        ! enlarge the system and build the full Hamiltonian
        enlarged_block = enlarge( this )
        H = superblock_hamiltonian( enlarged_block )
        call clear_block( this%lblock )

        ! compute eigenvalues
        call eigenvv_symmetric( this%energy_level, H, eigval, eigvec, do_eigvec=.true.)
        !   remark: eigvec is of size (H%size(2), 1)

        ! take the eigenvalue (1 for ground state)
        this%energy = eigval( this%energy_level )
        call Delete(H)

        ! compute reduced density matrix
        proj = reshape( eigvec%data(:, this%energy_level), &
                        (/ 2**(enlarged_block%mm), 2**(enlarged_block%mm) /)  )
        !   remark: this reshape allows an easy contraction to compute
        !           the reduced density matrix directly
        rho%data = matmul( proj, transpose(proj) )
        rho%size = (/ size(rho%data, 1), size(rho%data, 2) /)

        ! free eigval, eigvec
        deallocate(proj)
        call Delete(eigvec)
        deallocate(eigval)


        ! compute eigenvalues of reduced density matrix
        call eigenvv_symmetric( int(rho%size(1)), rho, eigval, eigvec, do_eigvec=.true.)
        !   remark:  Eigenvalues are listed in increasing algebraic order
        !           but we need the largest (decreasing algebraic),
        !           so we need to invert the order!

        ! set the crop factor wrt max m allowed
        crop_eigenvalues = min( rho%size(2), 2_8**this%max_m )

        ! I will use the cropped eigval to quantify the error
        eigsum = sum(abs(  eigval(1:int(rho%size(1)))  ))        
        eigtmp = 0d0

        ! reverse the order of eigenvectors to build the projector
        allocate( proj(rho%size(1), crop_eigenvalues )  )
        do ii = 1, crop_eigenvalues
            proj(:, ii) = eigvec%data(:, rho%size(2) + 1 - ii)
            eigtmp = eigtmp + abs(eigval(rho%size(2) + 1 - ii))
        end do
        call Delete( eigvec )
        deallocate( eigval )
        proj_dag = transpose( proj )


        ! housekeeping: I need this to allocate the new blocks
        sigmax = InitZero( (/2_8, 2_8/) )
        sigmax%data(2,1) = 1d0;       sigmax%data(1,2) = 1d0;


        ! project the enlarged building blocks with cropped eigenvalues
        this%lblock%HH%data = matmul( proj_dag, matmul(enlarged_block%HH%data, proj) )
        this%lblock%HH%size = (/ size(this%lblock%HH%data, 1), size(this%lblock%HH%data, 2) /)

        this%lblock%Hx%data = matmul( proj_dag, matmul(enlarged_block%Hx%data, proj) )
        this%lblock%Hx%size = (/ size(this%lblock%Hx%data, 1), size(this%lblock%Hx%data, 2) /)

        this%lblock%Hz%data = matmul( proj_dag, matmul(enlarged_block%Hz%data, proj))
        this%lblock%Hz%size = (/ size(this%lblock%Hz%data, 1), size(this%lblock%Hz%data, 2) /)


        deallocate( proj, proj_dag )
        call Delete( sigmax )


        ! toc the chronometer -----------------------------
        call CPU_TIME(time_toc)
        this%cputime = this%cputime + time_toc - time_tic

        ! increment counter & check max_iter --------------
        this%counter = this%counter + 1
        if ( this%max_iter .gt. 0 ) then
            if ( this%counter .gt. this%max_iter) then
                this%converged = .true.
                call checkpoint(level=DL_INFO, msg = 'DMRG converged by max iterations')
            end if
        end if

        ! set block size according to crop rule -----------
        this%lblock%mm = min( enlarged_block%mm, this%max_m)

        ! return energy density ---------------------------
        this%energy = this%energy/( 2*(this%counter) )
        if( this%save_hist ) then
            this%hist( this%counter-1 ) = this%energy
            this%hist_err( this%counter-1 ) = 1 - (eigtmp/eigsum)
        end if

    end subroutine






    ! ------------------------------------------------------------------
    !   I/O   |  write to file the parameters of DMRG object
    ! ------------------------------------------------------------------
    subroutine write_header_dmrg(this, filename)
        type(idmrga) :: this
        character(*), intent(in) :: filename

        integer, parameter :: fileunit = 22
        logical :: file_isopened

        open(fileunit, file=filename, action='write', form="unformatted")
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then
            call checkpoint(level=DL_INFO, msg = 'writing iDMRG config to (binary)'//filename)

            write(fileunit) this%max_m
            write(fileunit) this%max_iter

            close(fileunit)
            call checkpoint(level=DL_INFOP, msg = 'file written '//filename)
        else
            call checkpoint(level=DL_ERROR, msg = 'file '//filename//' cannot be opened')
            return
        end if
    end subroutine



    ! ------------------------------------------------------------------
    !   I/O   |  write to file the results of DMRG run
    ! ------------------------------------------------------------------
    subroutine write_entry_dmrg(this, filename)
        type(idmrga) :: this

        character(*), intent(in) :: filename

        integer, parameter :: fileunit = 22
        logical :: file_isopened

        open(fileunit, file=filename, status="old", action="write", &
                       form='unformatted', position='append')
        inquire(unit=fileunit, opened=file_isopened)

        if(file_isopened) then
            call checkpoint(level=DL_INFO, msg = 'writing IDMRG simulation entry to (binary)'//filename)

            write(fileunit) this%lambda
            write(fileunit) this%cputime
            write(fileunit) this%counter
            write(fileunit) this%energy
            write(fileunit) this%hist
            write(fileunit) this%hist_err
            
            close(fileunit)
            call checkpoint(level=DL_INFOP, msg = 'file written '//filename)
        else
            call checkpoint(level=DL_ERROR, msg = 'file '//filename//' cannot be opened')
            return
        end if
    end subroutine





end module idmrg_algo