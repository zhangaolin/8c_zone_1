module DriverTSoutput
    use GlobalTSVariables
    use GlobalTSConstants
    use driver_input_test
    use DataProcess
    use DriverSteadyTSSolver
    IMPLICIT NONE
    contains

    subroutine print_ts_output() 
        write(iunit,*)'================================='  
        write(iunit,*)'tests for heat conduction module'
        write(iunit,*)'================================='   
        write(iunit,*)'timestep=',timestep_solid
        write(iunit,*)'error=',error_Tsolid
        write(iunit,*)'coeff_b_r1=',coeff_b_r1
        write(iunit,*)'=================================' 
        write(iunit,*)'input and intermediate variables checking:' 

        do k=1,num_z
            do j=1,num_theta
                write(iunit,*)'q_vir:',q_vir(:,j,k)
            end do
        end do

        do k=1,num_z
            do j=1,num_theta
                write(iunit,*)'alpha:',alpha(:,j,k)
            end do
        end do

        write(iunit,*)'================================='  
        write(iunit,*)'output for heat conduction module'
        write(iunit,*)'================================='
        do k=1,num_z
            do j=1,num_theta
                write(iunit,*)'t_solid_average:',ts_ave(:,j,k)
            end do
        end do
        !call lu_solve_test()
    end subroutine

end module