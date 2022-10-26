module initialize_random_mod
    implicit none

    !This module initializes random number generater 

    integer, parameter :: double = kind(0d0)
    integer, parameter :: iseed = 1912382169
    integer,allocatable:: seed(:)
    integer seedsize
    private iseed, seed, seedsize

contains
    !!!!!!!!!! 乱数発生機構の設定 !!!!!!!!!

    subroutine initialize_random
    implicit none
    
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    seed = iseed
    call random_seed(put=seed)

    return
    end subroutine initialize_random

    real(double) function rnd()

    call random_number(rnd)

    return
    end function rnd

end module initialize_random_mod