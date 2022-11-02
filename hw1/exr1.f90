! ========================================
!  QIC  >>  assignment 1 / exr 1
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 21 October 2022
! ========================================

! This program is pointless... but it implements the 4th order
! central difference.

program test
  implicit none

  ! number of values to generate & dx step
  integer, parameter :: nn = 20000
  real, parameter :: dx = 0.001
  ! loop variable
  integer :: ii

  ! arrays
  real, dimension(:), allocatable :: y, dy

  ! allocate dynamic memory
  allocate(  y(nn)   )
  allocate( dy(nn-4) )

  ! initialize the vector
  do ii = 0,nn
      y(ii) = sin(ii*dx)
  end do

  ! compute 4th order central difference
  dy = 0 + y(1:nn-4)
  dy = dy - 8*y(2:nn-3)
  dy = dy + 8*y(4:nn-1)
  dy = dy - y(5:nn)
  dy = dy/(12*dx)

  !   print to file, if you wish ...
  !open(1, file = 'data/dy.txt')
  !do ii = 1, nn-4
  !  write(1,*) dy(ii)
  !end do
  !close(1)

  print *, dy(:100)

  ! say bye bye
  deallocate(y, dy)

end program test
