module DriverTSoutput
    use GlobalTSVariables
    use GlobalTSConstants
    !use driverthinput
    use driver_input_test
    use DataProcess
    IMPLICIT NONE
    !生成输出文件
contains

    subroutine print_ts_output()    
    write(iunit,*)'output for heat conduction module'

    !num_z_output_test
    write(iunit,*)'num_z=',num_z
    write(iunit,*)'num_r=',num_r

    !config_output_test
    do k=1,num_z
        do j=1,num_theta
            write(iunit,*)'config:',cfg(:,j,k)
        end do
    end do

    !qn_output_test

        do k=1,num_z
            do j=1,num_theta
                write(iunit,*)'qn:',power_den(:,j,k)
            end do
        end do
    !v_output_test

    do k=1,num_z
        do j=1,num_theta
            write(iunit,*)'v:',volume_node(:,j,k)
        end do
    end do
    
    end subroutine

end module