! ========================================
!  QIC  >>  assignment 5
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 28 November 2022
! ========================================



module trotter_suzuki
    use, intrinsic :: iso_c_binding ! to make fftw3 work correctly
    use mydebugger

    implicit none


    ! ------------------------------------------------------------------
    !  TYPE |   Data structure for Trotter-Suzuki simulation
    ! ------------------------------------------------------------------
    type :: trotter_suzuki_sim
        ! [parameters] space-time mesh
        double precision :: xmax, tmax
        integer :: Nx, Nt

        ! [parameters] harmonic oscillator & q(t)
        real*8 :: omega, tau
        ! we take m = 1

        ! [ globals ]
        double complex, dimension(:), allocatable :: psi
        real*8 :: time

        ! [ internals ]
        integer :: time_counter = 0
        real*8, dimension(:), allocatable, private :: HV ! potential
        real*8, dimension(:), allocatable, private :: HT ! kinetic

    end type trotter_suzuki_sim


    ! in case you wish to use a custom constructor
    interface trotter_suzuki_sim
        module procedure init_trotter_suzuki
    end interface trotter_suzuki_sim
    
contains
    
    ! ------------------------------------------------------------------
    !  EASY FFT |   Computes the Fourier transform of a complex vector
    ! ------------------------------------------------------------------
    function FFT(in, N) result(out)
        double complex, dimension(:) :: in
        double complex, dimension( N ) :: out
        integer :: N
        integer*8 :: plan

        integer, parameter :: FFTW_ESTIMATE = 64

        call dfftw_plan_dft_1d(plan, N, in, out, -1, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, in, out)

        ! [ ref https://www.fftw.org/fftw3.pdf ] -------------------
        ! The fourth argument, sign, can be either FFTW_FORWARD (-1) or FFTW_BACKWARD (+1), and
        ! indicates the direction of the transform you are interested in; technically, it is the sign of
        ! the exponent in the transform.
        ! The flags argument is usually either FFTW_MEASURE or FFTW_ESTIMATE. FFTW_MEASURE
        ! instructs FFTW to run and measure the execution time of several FFTs in order to find the
        ! best way to compute the transform of size n. This process takes some time (usually a few
        ! seconds), depending on your machine and on the size of the transform. FFTW_ESTIMATE,
        ! on the contrary, does not run any computation and just builds a reasonable plan that is
        ! probably sub-optimal.
        ! ----------------------------------------------------------

        ! about normalization --------------------------------------
        ! FFTW computes an unnormalized DFT. Thus, computing a forward followed by a backward
        ! transform (or vice versa) results in the original array scaled by n.
        out = out/sqrt( cmplx(N) )
    
        call dfftw_destroy_plan(plan)
    end function


    ! ------------------------------------------------------------------
    !  EASY FFT |   Computes the inverse Fourier transform
    ! ------------------------------------------------------------------
    !     notes : the ordering of input vector must be such that
    ! 
    !       -k/2          0            k/2
    !        --------------------------->  discretized p_i = 2 PI i / (2 *xmax)
    !                     [ ==== A ==== ]
    !        [ ==== B === )
    !  ------------------------------------------------------------------
    function invFFT(in, N) result(out)
        double complex, dimension(:) :: in
        double complex, dimension( N ) :: out
        integer :: N
        integer*8 :: plan

        integer, parameter :: FFTW_ESTIMATE = 64

        call dfftw_plan_dft_1d(plan, N, in, out, +1, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, in, out)

        out = out/sqrt( cmplx(N) )

        call dfftw_destroy_plan(plan)
    end function




    ! ------------------------------------------------------------------
    !  LINEAR UTILS | normalize to 1 the module square of a wavefunction
    ! ------------------------------------------------------------------
    subroutine normalize(psi)
        double complex, dimension(:), intent(inout) :: psi
        psi = psi/sqrt( sum( real( psi*conjg(psi) ) ) )
    end subroutine

    ! ------------------------------------------------------------------
    !  LINEAR UTILS | print the module square of a wavefunction
    ! ------------------------------------------------------------------
    subroutine check_norm(psi)
        double complex, dimension(:) :: psi
        real*8 :: norm

        norm = sum( real( psi*conjg(psi) ) )
        print *, 'norm =', norm
    end subroutine



    ! ------------------------------------------------------------------
    !  HAMILTONIAN | compute the potential (V) at given time
    ! ------------------------------------------------------------------
    subroutine simV(this, time)
        type(trotter_suzuki_sim), intent(inout) :: this
        real*8 :: time, x
        integer ii

        do ii = 1, this%Nx
            x = this%xmax*(-1d0 + (ii-1)*2d0/(this%Nx-1))
            this%HV(ii) = 0.5d0*(this%omega**2) * ((x-(time/this%tau))**2)
        end do
    end subroutine

    ! ------------------------------------------------------------------
    !  HAMILTONIAN | compute the kinetic term (T) at given time
    ! ------------------------------------------------------------------
    subroutine simT(this, time)
        type(trotter_suzuki_sim), intent(inout) :: this
        real*8 :: time
        real*8, parameter :: pi = 4.0d0*datan(1.0d0)
        integer ii

        ! be careful for ordering after FFT ------------------------
        !  from https://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html
        ! For those who like to think in terms of positive and negative frequencies, this means
        ! that the positive frequencies are stored in the first half of the output and the 
        ! negative frequencies are stored in backwards order in the second half of the output.
        ! (The frequency -k/n is the same as the frequency (n-k)/n.)

        do ii = 1, this%Nx
            ! p_ii = 2 Pi ii / (2 *xmax) = (Pi/xmax) * ii
            !  then  T = (p_ii)^2 / 2m
            if(ii .le. this%Nx/2 ) then
                this%HT(ii) = 0.5d0*( (ii*pi/this%xmax)**2 )
            else
                this%HT(ii) = 0.5d0*( ((this%Nx-ii)*pi/this%xmax)**2 )
            end if

        end do

    end subroutine



    ! ------------------------------------------------------------------
    !  TSS | initialize a Trotter-Suzuki simulation
    ! ------------------------------------------------------------------
    function init_trotter_suzuki(xmax, tmax, tau, Nx, Nt, psi0) result(this)
        type(trotter_suzuki_sim) :: this
        double precision :: xmax, tmax, tau
        integer :: Nx, Nt
        real*8, dimension(:) :: psi0

        call checkpoint(level=DL_PEDANTIC, msg = 'called trotter suzuki init')

        this%xmax = xmax;   this%tmax = tmax;
        this%Nx = Nx;       this%Nt = Nt;

        this%time_counter = 0

        allocate( this%HV(Nx), this%HT(Nx), this%psi(Nx) )

        this%omega = 1d0
        this%tau = tau ! or this%tmax, as written in the assignment...

        call checkpoint(level=DL_INFO, msg = 'parsing psi zero')
        ! parse input state from real to complex vector
        this%psi = cmplx( psi0, kind=4)
    end function


    ! ------------------------------------------------------------------
    !  TSS | print the most important parameters of a TSS object
    ! ------------------------------------------------------------------
    subroutine tellme_ts(this)
        type(trotter_suzuki_sim) :: this

        print *, "Trotter-Suzuki implementation ---------"
        print *, " xmax =", this%xmax
        print *, " tmax =", this%tmax
        print *, " Nx   =", this%Nx
        print *, " Nt   =", this%Nt
        print *, " - - - - - - - - - - -"
        print *, " omega =", this%omega
        print *, " tau   =", this%tau
    end subroutine


    ! ------------------------------------------------------------------
    !   TSS   |   deallocate a Trotter-Suzuki simulation
    ! ------------------------------------------------------------------
    subroutine destroy_ts(this)
        type(trotter_suzuki_sim) :: this
        deallocate( this%HV, this%HT, this%psi )
    end subroutine


    ! ------------------------------------------------------------------
    !   TSS   |   compute a time step of a TS simulation
    ! ------------------------------------------------------------------
    subroutine time_evolution(this)
        type(trotter_suzuki_sim), intent(inout) :: this
        real*8 :: time

        ! computing current time
        time = (this%time_counter)*(this%tmax)/(this%Nt)
        
        ! update V and T for current time
        call simV(this, time)
        call simT(this, time)

        ! Trotter-Suzuki approximation -----------------------------
        this%psi = this%psi * zexp( cmplx(0d0, -0.5d0*(this%HV)*(this%tmax)/(this%Nt), kind=8) )

        this%psi = FFT(this%psi, this%Nx)

        this%psi = this%psi * zexp( cmplx(0d0, -1d0*(this%HT)*(this%tmax)/(this%Nt), kind=8) )

        this%psi = invFFT(this%psi, this%Nx)

        this%psi = this%psi * zexp( cmplx(0d0, -0.5d0*(this%HV)*(this%tmax)/(this%Nt), kind=8) )
        
        ! flush
        this%time_counter = this%time_counter + 1
        this%time = time
    end subroutine


end module trotter_suzuki
