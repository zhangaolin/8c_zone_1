module GlobalTSConstants 
!定义全局常量
    IMPLICIT NONE

    ! integer :: 
    ! integer,allocatable :: 
    ! !便于后续切换为双精度
     integer,parameter :: TS_DOUBLE = 4 
     real(TS_DOUBLE) :: convergence_limit=0.00001
    ! real(TS_KOUBLE),allocatable :: 
    ! real(TS_KOUBLE),parameter :: 


end module