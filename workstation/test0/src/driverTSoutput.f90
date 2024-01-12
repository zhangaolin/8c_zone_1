module DriverTSoutput
    use GlobalTSVariables
    use GlobalTSConstants
    use driver_input_test
    use DataProcess
    use DriverSteadyTSSolver
    IMPLICIT NONE
    contains

    subroutine print_ts_output() 

        ! write(iunit,*)'================================='  
        ! write(iunit,*)'output for heat conduction module'
        ! write(iunit,*)'================================='
        do k=1,num_z
            do j=1,num_theta
                write(iunit,*)'t_solid_average:',ts_ave(:,j,k)
            end do
        end do
        ! write(iunit,*)'=================================' 
        
        
        write(iunit,*)'tests for heat conduction module'
        write(iunit,*)'================================='   
        write(iunit,*)'timestep=',timestep_solid
        write(iunit,*)'=================================' 
        write(iunit,*)'error=',error_Tsolid
        write(iunit,*)'=================================' 
        write(iunit,*)'input and intermediate variables checking:' 
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_b_r1=',coeff_b_r1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_b_r2=',coeff_b_r2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_b_r3=',coeff_b_r3
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_b_z1=',coeff_b_z1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_b_z2=',coeff_b_z2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_b_z3=',coeff_b_z3
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_b_th1=',coeff_b_th1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_b_th2=',coeff_b_th2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_b_th3=',coeff_b_th3
        write(iunit,*)'=================================' 

        write(iunit,*)'coeff_c_r1=',coeff_c_r1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_c_r2=',coeff_c_r2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_c_r3=',coeff_c_r3
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_c_z1=',coeff_c_z1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_c_z2=',coeff_c_z2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_c_z3=',coeff_c_z3
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_c_th1=',coeff_c_th1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_c_th2=',coeff_c_th2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_c_th3=',coeff_c_th3
        write(iunit,*)'=================================' 

        write(iunit,*)'coeff_d_r1=',coeff_d_r1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_d_r2=',coeff_d_r2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_d_r3=',coeff_d_r3
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_d_z1=',coeff_d_z1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_d_z2=',coeff_d_z2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_d_z3=',coeff_d_z3
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_d_th1=',coeff_d_th1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_d_th2=',coeff_d_th2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_d_th3=',coeff_d_th3
        write(iunit,*)'=================================' 

        write(iunit,*)'coeff_e_r1=',coeff_e_r1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_e_r2=',coeff_e_r2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_e_r3=',coeff_e_r3
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_e_z1=',coeff_e_z1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_e_z2=',coeff_e_z2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_e_z3=',coeff_e_z3
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_e_th1=',coeff_e_th1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_e_th2=',coeff_e_th2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_e_th3=',coeff_e_th3
        write(iunit,*)'=================================' 

        write(iunit,*)'coeff_f_r1=',coeff_f_r1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_f_r2=',coeff_f_r2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_f_r3=',coeff_f_r3
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_f_z1=',coeff_f_z1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_f_z2=',coeff_f_z2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_f_z3=',coeff_f_z3
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_f_th1=',coeff_f_th1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_f_th2=',coeff_f_th2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_f_th3=',coeff_f_th3
        write(iunit,*)'=================================' 

        write(iunit,*)'coeff_g_r1=',coeff_g_r1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_g_r2=',coeff_g_r2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_g_r3=',coeff_g_r3
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_g_z1=',coeff_g_z1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_g_z2=',coeff_g_z2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_g_z3=',coeff_g_z3
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_g_th1=',coeff_g_th1
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_g_th2=',coeff_g_th2
        write(iunit,*)'=================================' 
        write(iunit,*)'coeff_g_th3=',coeff_g_th3
        write(iunit,*)'=================================' 
        ! do k=1,num_z
        !     do j=1,num_theta
        !         write(iunit,*)'q_vir:',q_vir(:,j,k)
        !     end do
        ! end do
        ! write(iunit,*)'=================================' 
        ! do k=1,num_z
        !     do j=1,num_theta
        !         write(iunit,*)'coeff_q_r:',coeff_q_r(:,j,k)
        !     end do
        ! end do
        ! write(iunit,*)'=================================' 
        ! do k=1,num_z
        !     do j=1,num_theta
        !         write(iunit,*)'coeff_q_theta:',coeff_q_th(:,j,k)
        !     end do
        ! end do
        ! write(iunit,*)'=================================' 
        ! do k=1,num_z
        !     do j=1,num_theta
        !         write(iunit,*)'coeff_q_z:',coeff_q_z(:,j,k)
        !     end do
        ! end do
        ! write(iunit,*)'=======================================' 
        ! do k=1,num_z
        !     do j=1,num_theta
        !         write(iunit,*)'alpha:',alpha(:,j,k)
        !     end do
        ! end do
        ! write(iunit,*)'=======================================' 
        ! do k=1,num_z
        !     write(iunit,*)'sub_z:',c(k)
        ! end do
        ! write(iunit,*)'=======================================' 
        ! do j=1,num_theta
        !     write(iunit,*)'sub_theta:',b(j)
        ! end do
        ! write(iunit,*)'=======================================' 
        ! do i=1,num_r
        !     write(iunit,*)'sub_r:',a(i)
        ! end do
        ! write(iunit,*)'=======================================' 

        ! write(iunit,*)'=======================================' 
        ! do k=1,num_z
        !     do j=1,num_theta
        !         write(iunit,*)'lambda:',lambda(:,j,k)
        !     end do
        ! end do
        ! write(iunit,*)'=======================================' 
        ! do k=1,num_z
        !     do j=1,num_theta
        !         write(iunit,*)'leakage:',z_leakage(:,j,k)
        !     end do
        ! end do
        ! write(iunit,*)'=======================================' 

    end subroutine

end module