! ========================================
!  Quantum Information and Computing
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   last update : 09 November 2022
! ========================================


! ========================================
! QUICK REFERENCE: verbosity levels ------
!  0   CHECKPOINT (minimal)
!  1   FATAL
!  2   ERROR
!  3   WARNING
!  4   INFO (blue)
!  5   INFO (green)
!  6   CHECKPOINT (pedantic)
! ========================================


module mydebugger
!  DOC =============================================================
! 
!    This module implements a checkpoint function with
!    many levels of verbosity.
!  -----------------------------------------------------------------
!
!  GLOBAL VARIABLES    |  DEBUG       (logical)  If true, the checkpoints 
!                      |       are printed in stdout
!                      |  DEBUG_LEVEL (integer) Set the verbosity level
!                      |       of the checkpoints. If DEBUG_LEVEL is 
!                      |       smaller than the checkpoint verbosity level,
!                      |       the message is not printed in stdout.
!
!  TYPES                 none
! 
!  FUNCTIONS             no public methods
! 
!  SUBROUTINES:
!   -  debug_status()
!         Prints the values of the global variables DEBUG and DEBUG_LEVEL.
!
!   -  list_debug_verbosity()
!         Prints the verbosity levels & the corresponding integer value.
!    
!   -  checkpoint(level, flag, msg, var, wait)
!
!         Creates a checkpoint, which is printed only if DEBUG_LEVEL is
!         greater than the checkpoint level.
!
!           Inputs  |  level   (integer) (optional)
!                   |     The verbosity level of this checkpoint.
!                   |
!                   |  flag    (string) (optional)
!                   |     The text to display at beginning of the message
!                   |     between square brackets. If not provided, the
!                   |     debugger will print a default value associated
!                   |     to the checkpoint level.
!                   |
!                   |  msg     (string) (optional)
!                   |     Optional message to be printed.
!                   |
!                   |  var     (---) (optional)
!                   |     Optional variable to be printed.
!                   |     Supported types:  character, logical,
!                   |                       integer*1, *2, *4,
!                   |                       real, double,
!                   |                       complex*4, complex*8
!                   |  wait    (logical) (optional) (default: FALSE)
!                   |     If True, the execution of the code is paused
!                   |     until the user hits Enter.
!
!            Requirements  |   DEBUG and DEBUG_LEVEL globals
! 
!                  Output  |   none   (just print in stdout)
! 
!         Post-conditions  |   none
! 
! ==================================================================
    implicit none
    
    ! globally accessible variables
    logical :: DEBUG = .TRUE.
    integer :: DEBUG_LEVEL = 1

    ! enumerate the verbosity levels, for convenience
    enum, bind(c)
        enumerator :: DL_CHECK = 0
        enumerator :: DL_FATAL, DL_ERROR
        enumerator :: DL_WARN
        enumerator :: DL_INFO, DL_INFOP
        enumerator :: DL_PEDANTIC
    endenum

    ! private methods
    private :: manage_default_flag
contains

    ! ------------------------------------------------------------------
    !   Prints the values of the global variables DEBUG and DEBUG_LEVEL
    ! ------------------------------------------------------------------
    subroutine debug_status()
        print *, "DEBUG :", DEBUG
        print *, "LEVEL :", DEBUG_LEVEL
    end subroutine

    ! ------------------------------------------------------------------
    !   Prints the verbosity levels & the corresponding integer value.
    ! ------------------------------------------------------------------
    subroutine list_debug_verbosity()
        print *, "DL_CHECK    :", DL_CHECK
        print *, "DL_FATAL    :", DL_FATAL
        print *, "DL_ERROR    :", DL_ERROR
        print *, "DL_WARN     :", DL_WARN
        print *, "DL_INFO     :", DL_INFO
        print *, "DL_INFOP    :", DL_INFOP
        print *, "DL_PEDANTIC :", DL_PEDANTIC
    end subroutine

    ! ------------------------------------------------------------------
    !   Creates a checkpoint, which is printed only if DEBUG_LEVEL is
    !  greater than the checkpoint level.
    ! ------------------------------------------------------------------
    subroutine checkpoint(level, flag, msg, var, wait)

        ! args --------------------------
        integer, optional, intent(in) :: level
        character(*), optional, intent(in) :: msg, flag
        class(*), optional, intent(in) :: var
        logical, optional, intent(in) :: wait
        
        ! default values ----------------
        integer :: level_
        logical :: wait_
        if(present(level)) then
            level_ = level;
        else ! use the else to prevent memory bound out-reference
            level_ = 0;
        end if
        if(present(wait)) then
            wait_ = wait
        else
            wait_ = .false.
        end if

        ! check if debug is active ------
        if(DEBUG .eqv. .FALSE.) then
            return
        end if
        
        ! SKIP if checkpoint level is higher than DEBUG_LEVEL
        if(level_ .gt. DEBUG_LEVEL) then
            return
        end if

        !  print flag
        if( present(flag) ) then
            ! use the flag in arg ...
            write (*, fmt='(A)', advance='NO') '['//flag//'] '
        else
            ! ... or pick a default value if not present
            write (*, fmt='(A)', advance='NO') manage_default_flag(level_)
        end if

        !  print message (if provided)
        if( present(msg) ) then
            print *, msg
        end if

        ! print variable
        if(present(var)) then
            ! compiler requires an interface to print the variable!
            select type(var)
                type is (character(*))
                    print *, "  (char) ", var
                type is (  logical )
                    print *, "  (logical) ", var
                type is (integer(1))
                    print *, "  (int1) ", var
                type is (integer(2))
                    print *, "  (int2) ", var
                type is (integer(4))
                    print *, "  (int4) ", var
                type is (  real(4) )
                    print *, "  (real) ", var
                type is (  real(8) )
                    print *, "  (double) ", var
                type is (complex(4))
                    print *, "  (complex4) ", var
                type is (complex(8))
                    print *, "  (complex8) ", var
            end select
        end if

        if(wait_) then
            print *, "-> ["//achar(27)//"[0;36m"//"pause"//achar(27)//"[0m] press enter to continue"
            read(*,*)
            print *, '-> resume'
        end if

    end subroutine

    function manage_default_flag(level) result(flag)
        integer :: level
        character(:) ,allocatable :: flag

        select case (level)
            case (DL_FATAL)
                flag = "["//achar(27)//"[1;31m"//" fatal "//achar(27)//"[0m]"
            case (DL_ERROR)
                flag = "["//achar(27)//"[1;31m"//" error "//achar(27)//"[0m]"
            case (DL_WARN)
                flag = "["//achar(27)//"[1;33m"//"warning"//achar(27)//"[0m]"
            case (DL_INFO)
                flag = "["//achar(27)//"[1;34m"//" info  "//achar(27)//"[0m]"
            case (DL_INFOP)
                flag = "["//achar(27)//"[1;32m"//" info  "//achar(27)//"[0m]"
            case (DL_PEDANTIC)
                flag = "["//achar(27)//"[0;35m"//"  dbg  "//achar(27)//"[0m]"
            case default
                flag = '[checkpoint]'
        end select
    end function

end module mydebugger

