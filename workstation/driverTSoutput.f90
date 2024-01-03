module DriverTSoutput
    use GlobalTSVariables
    use GlobalTSConstants
    IMPLICIT NONE
    !生成输出文件
    contains
    subroutine print_ts_output()    
    write(iunit,*)'output for heat conduction module'
    write(iunit,*)!
    end subroutine
end module