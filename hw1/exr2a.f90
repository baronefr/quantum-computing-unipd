! ========================================
!  QIC  >>  assignment 1 / exr 2A
!  UniPD, AY 2022/23, Physics of Data
! ----------------------------------------
!   coder : Barone Francesco
!   dated : 21 October 2022
! ========================================

program int_overflow

    implicit none

    ! variables for int*2 test
    integer*2 :: big_int2 = 2000000
    integer*2 :: res_int2
    ! variables for int*4 test
    integer*4 :: big_int4 = 2000000
    integer*4 :: res_int4

    
    print *, "twomillion increment test -------"
    ! doing the actual test
    res_int2 = big_int2 + 1
    print *, " sum with int*2 =", res_int2

    res_int4 = big_int4 + 1
    print *, " sum with int*4 =", res_int4


    ! useful things to know ---------------------------
    print *, "[i] direct overflow test ------"
    big_int2 = huge(big_int2) + 1
    big_int4 = huge(big_int4) + 1
    print *, " int*2 :", huge(big_int2), "+ 1  != ", big_int2
    print *, " int*4 :", huge(big_int4), "+ 1  != ", big_int4

end program int_overflow
