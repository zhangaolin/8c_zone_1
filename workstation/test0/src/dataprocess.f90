module DataProcess !data processing

use GlobalTSVariables
use GlobalTSConstants
use driver_input_test

IMPLICIT NONE

contains

! subroutine channel_approximation() !流道近似，调整节块尺寸

! end subroutine

! subroutine alpha_calculation() !计算每个节块的固气换热系数
!     IMPLICIT NONE
!     real(TS_KDUBLE),allocatable :: alpha(:,:,:)
!     allocate(alpha(z_num,r_num,theta_num)) 
!     do i=1,z_num
!         do j=1,r_num
!             do k=1,theta_num
!                 if(composition_solidmodule(node_config(i,j,k))%phase_type == 3.)then
!                     alpha(i,j,k)=h_solidmodule(i,j,k)*6*(1-composition_solidmodule(node_config(i,j,k))%epsilon)/pebble_diameter 
!                 else 
!                     alpha(i,j,k)=0
!                 end if
!             end do
!         end do
!     end do
! end subroutine


subroutine geo_process() !计算节块的r值和体积
    !r_calculate
    allocate(r(num_r,num_z,num_theta))  
        do j=1,num_theta
            do k=1,num_z
                any_real=0.0d0
                do i=1,num_r             
                    any_real=any_real+sub_r(i)                              
                    r(i,j,k)=any_real-0.5*sub_r(i)                 
            end do
        end do
    end do
    !v_calculate
    allocate(volume_node(num_r,num_z,num_theta))  
    do j=1,num_theta
        do k=1,num_z
            do i=1,num_r 
            volume_node(i,j,k)=pi*sub_r(i)*sub_theta(j)*sub_z(k)*r(i,j,k)/180.0d0
        end do
    end do
end do
end subroutine


subroutine alpha_test0() !用于test0,直接给每种材料的alpha赋0值
    allocate(alpha(num_r,num_z,num_theta))
    do i=1,num_r
        do j=1,num_theta
            do k=1,num_z
                if(cfg(i,j,k)==1.)then
                    alpha(i,j,k)=0
                else if(cfg(i,j,k)==2.)then
                    alpha(i,j,k)=0
                end if
            end do
        end do
    end do
end subroutine

subroutine q_calculate() !用于test0,直接给每种材料的alpha赋0值
    allocate(power_den(num_r,num_z,num_theta))
    allocate(q_vir(num_r,num_z,num_theta))
    do i=1,num_r
        do j=1,num_theta
            do k=1,num_z
                if(cfg(i,j,k)==1.)then
                    power_den(i,j,k)=0
                else if(cfg(i,j,k)==2.)then
                    power_den(i,j,k)=0   
                end if 
                q_vir(i,j,k)=power_den(i,j,k)+alpha(i,j,k)*tg(i,j,k)
            end do
        end do
    end do
end subroutine

subroutine data_processing() 
    IMPLICIT NONE
    call geo_process()
    ! call alpha_test0()
    ! call q_calculate()
    !call channel_approximation()
   ! call alpha_calculation()
end subroutine

end  module