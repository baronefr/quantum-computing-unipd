! ========================================
!  Quantum Information and Computing
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   last update : 10 November 2022
! ========================================

module mytools
    use mydebugger

    implicit none

contains

    ! set the random iseed for LAPACK library from given values
    subroutine iseed_manual(iseed, s1, s2, s3, s4)
        implicit none
        integer :: iseed(4)
        integer, intent(in) :: s1, s2, s3, s4

        if( (s1.ge.0).and.(s1.le.4095) ) then
            iseed(1) = s1
        else
            call checkpoint(level=DL_ERROR, msg = 'iseed(1) not in [0, 4095] ... setting =0')
            iseed(1) = 0
        end if

        if( (s2.ge.0).and.(s2.le.4095) ) then;
            iseed(2) = s2
        else
            call checkpoint(level=DL_ERROR, msg = 'iseed(2) not in [0, 4095] ... setting =0')
            iseed(2) = 0
        end if

        if( (s3.ge.0).and.(s3.le.4095) ) then
            iseed(3) = s3
        else
            call checkpoint(level=DL_ERROR, msg = 'iseed(3) not in [0, 4095] ... setting =0')
            iseed(3) = 0
        end if

        if( (s4.ge.0).and.(s4.le.4095) ) then
            if( Mod(s4,2).eq.1 ) then
                iseed(4) = s4
            else
                call checkpoint(level=DL_ERROR, msg = 'seed(4) not odd ... setting =1')
                iseed(4) = 1
            end if
        else
            call checkpoint(level=DL_ERROR, msg = 'iseed(4) not in [0, 4095] ... setting =1')
            iseed(4) = 0
        end if
    end subroutine iseed_manual


    ! set the random iseed for LAPACK library from random values
    subroutine iseed_now(iseed)
        real :: val(4)
        integer :: iseed(4)
        integer :: seed(4)

        ! take some rands
        call random_number(val)
        seed = FLOOR(val*4095)

        ! if last value is even, make it odd
        if( Mod(seed(4),2).eq.0 ) seed(4) = seed(4) + 1

        call checkpoint(level=DL_INFO, msg = 'setting iseed:', var = seed(1))
        call iseed_manual(iseed, seed(1),seed(2),seed(3),seed(4))
    end subroutine

    
    ! write an array of real*8 to txt file (line by line)
    subroutine dump_real8_array_txt(filename, array, first_index, last_index, append)
        implicit none
        character(len=256) :: filename
        double precision, dimension(:) :: array
        logical :: append, file_isopened
        integer :: first_index, last_index

        integer :: ii
        integer, parameter :: fileunit = 22

        if (trim(filename) == '') then
            call checkpoint(level=DL_ERROR, msg = "skipping file operation, arg is empty")
        else

            if(append) then
                call checkpoint(level=DL_INFO, msg = "opening file "//trim(filename)//" -> append")
                open(fileunit, file=filename, action='write', position='append')
            else
                call checkpoint(level=DL_INFO, msg = "opening file "//trim(filename)//" -> write")
                open(fileunit, file=filename, action='write')
            end if

            inquire(unit=fileunit, opened=file_isopened)
            if(file_isopened) then
                call checkpoint(level=DL_INFO, msg = 'writing file '//trim(filename))
                do ii = first_index,last_index
                    write(fileunit, '(1G0)') array(ii)
                end do
            else
                call checkpoint(level=DL_FATAL, msg = 'file '//trim(filename)//' cannot be opened')
            end if

        end if
    end subroutine

end module mytools
