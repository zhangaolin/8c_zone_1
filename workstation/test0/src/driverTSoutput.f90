module DriverTSoutput
    use GlobalTSVariables
    use GlobalTSConstants
    !use driverthinput
    use driver_input_test
    use DataProcess
    use DriverSteadyTSSolver
    IMPLICIT NONE
    !生成输出文件

contains

    subroutine print_ts_output()    
    write(iunit,*)'output for heat conduction module'
    write(iunit,*)'================================='
    write(iunit,*)'t_solid_average:',ts_ave(:,:,1)
    end subroutine

end module