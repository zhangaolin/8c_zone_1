module DriverTSoutput
    use GlobalTSVariables
    use GlobalTSConstants
    use driverthinput
    IMPLICIT NONE
    !生成输出文件
contains

    subroutine print_ts_output()    
    write(iunit,*)'output for heat conduction module'
    write(iunit,*)'num_z=',num_z
    write(iunit,*)'num_r=',num_r
    do k=1,num_z
        do j=1,num_theta
            write(iunit,*)'config:',cfg(:,j,k)
        end do
    end do
    write(iunit,*)'r=',r
    end subroutine

end module