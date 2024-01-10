module GlobalTSConstants 
    use, intrinsic :: iso_fortran_env
    IMPLICIT NONE

    ! integer :: 
    ! integer,allocatable :: 
    ! !便于后续切换为双精度
     integer,parameter :: TS_DOUBLE = 8
     real(TS_DOUBLE) :: convergence_limit=0.00001
    ! real(TS_KOUBLE),allocatable :: 
     real(TS_DOUBLE),parameter :: pi=3.1415926


end module