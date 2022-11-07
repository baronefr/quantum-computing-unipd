! ========================================
!  QIC  >>  assignment 2 / exr 1
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 28 October 2022
! ========================================

!  This programs tests the 'mydebugger' module, which
! implements a checkpoint() subroutine. 
program main
    use mydebugger
    implicit none

    ! set debug & debug level
    DEBUG = .true.
    DEBUG_LEVEL = DL_CHECK  ! btw, this is the default value

    ! check the debug & debug level used
    call debug_status()
    print *, "-- list of verbosity levels --"
    call list_debug_verbosity()
    print *, "-----------"


    ! test messages
    call checkpoint(msg = 'your message')
    call checkpoint(flag = ' -!- ', var = 1.0)

    ! this is not executed because DEBUG_LEVEL is lower than 3
    call checkpoint(level=2, msg = 'I will not be printed')

    ! wait for user imput
    call checkpoint(msg = 'now we will test higher debug levels', wait=.true.) 

    print *, ""
    print *, ""

    ! change the debug level
    DEBUG_LEVEL = DL_PEDANTIC

    ! now you see this ...
    call checkpoint(level=2, msg = 'print me (2) error message')
    call checkpoint(level=DL_WARN, msg = 'print me (3) warning')
    call checkpoint(level=DL_INFO, msg = 'print me (4) info')
    call checkpoint(level=DL_INFOP, msg = 'print me (5) info but green')
    call checkpoint(level=DL_PEDANTIC, msg = 'print me (6) pedantic', var=DEBUG_LEVEL)

    ! ... but not this one
    call checkpoint(level=7, msg = '!! do NOT print me!! (7)')

    ! say goodbye
    call checkpoint(msg = 'last test, have a good day')

end program main
