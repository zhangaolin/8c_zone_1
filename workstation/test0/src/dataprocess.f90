module DataProcess !data processing

use GlobalTSVariables
use GlobalTSConstants
use driver_input_test

IMPLICIT NONE

contains

! subroutine channel_approximation() !流道近似，调整节块尺寸
!implicit none
!!================undone==================================
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


subroutine geo_process() !calculate radial length and volume for every node
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


subroutine alpha_calculate() !only for test0,give alpha as 0 for every node
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

subroutine q_calculate() 
    allocate(power_den(num_r,num_z,num_theta))
    allocate(q_vir(num_r,num_z,num_theta))
    do i=1,num_r
        do j=1,num_theta
            do k=1,num_z
                if(cfg(i,j,k)==1.)then ! 1D problem, material A
                    power_den(i,j,k)=1500000
                else if(cfg(i,j,k)==2.)then ! 1D problem, material B
                    power_den(i,j,k)=0   
                end if 
                q_vir(i,j,k)=power_den(i,j,k)+alpha(i,j,k)*tg(i,j,k)
            end do
        end do
    end do
end subroutine

subroutine chase()
    implicit none
    u_ch(1)=b_ch(1)
    y_ch(1)=d_ch(1)
    do i_ch=2,num_ch
        l_ch(i_ch)=a_ch(i_ch)/u_ch(i_ch-1)
        u_ch(i_ch)=b_ch(i_ch)-l_ch(i_ch)*c_ch(i_ch-1)
        y_ch(i_ch)=d_ch(i_ch)-l_ch(i_ch)*y_ch(i_ch-1)
    end do
    x_ch(num_ch)=y_ch(num_ch)/u_ch(num_ch)
    do i_ch=num_ch-1,1,-1
        x_ch(i_ch)=(y_ch(i_ch)-c_ch(i_ch)*x_ch(i_ch+1))/u_ch(i_ch)
    end do    
end subroutine


subroutine lu_solve()
    implicit none
    l_lu = 0.0d0
    u_lu = 0.0d0
    y_lu = 0.0d0
    x_lu = 0.0d0
    do j_lu = 1, n_lu
        u_lu(j_lu, 1) = A_lu(j_lu, 1)
    end do
    do j_lu = 2, n_lu
        l_lu(1, j_lu) = A_lu(1, j_lu) / u_lu(1, 1)
    end do
    do i_lu = 2, n_lu-1
        do j_lu = 1, i_lu-1
            u_lu(i_lu, i_lu) = u_lu(i_lu, i_lu) - l_lu(j_lu, i_lu) * u_lu(i_lu, j_lu)
        end do
        u_lu(i_lu, i_lu) = u_lu(i_lu, i_lu) + a_lu(i_lu, i_lu)
        do j_lu = i_lu+1, n_lu
            do k_lu = 1, i_lu-1
                u_lu(j_lu, i_lu) = u_lu(j_lu, i_lu) - l_lu(k_lu, i_lu) * u_lu(j_lu, k_lu)
                l_lu(i_lu, j_lu) = l_lu(i_lu, j_lu) - l_lu(k_lu, j_lu) * u_lu(i_lu, k_lu)
            end do
            u_lu(j_lu, i_lu) = u_lu(j_lu, i_lu) + a_lu(j_lu, i_lu)
            l_lu(i_lu, j_lu) = (l_lu(i_lu, j_lu) + a_lu(i_lu, j_lu)) / u_lu(i_lu, i_lu)
        end do
    end do
    do i_lu = 1, n_lu-1
        u_lu(n_lu, n_lu) = u_lu(n_lu, n_lu) - l_lu(i_lu, n_lu) * u_lu(n_lu, i_lu)
    end do
    u_lu(n_lu, n_lu) = u_lu(n_lu, n_lu) + a_lu(n_lu, n_lu)
    do i_lu = 1, n_lu
        l_lu(i_lu, i_lu) = 1.0d0
    end do
    y_lu(1) = b_lu(1)
    do i_lu = 2, n_lu
        do j_lu = 1, i_lu-1
            y_lu(i_lu) = y_lu(i_lu) - y_lu(j_lu) * l_lu(j_lu, i_lu)
        end do
        y_lu(i_lu) = y_lu(i_lu) + b_lu(i_lu)
    end do
    x_lu(n_lu) = y_lu(n_lu) / u_lu(n_lu, n_lu)
    do i_lu = n_lu-1, 1, -1
        do j_lu = i_lu+1, n_lu
            x_lu(i_lu) = x_lu(i_lu) - x_lu(j_lu) * u_lu(j_lu, i_lu)
        end do
        x_lu(i_lu) = (x_lu(i_lu) + y_lu(i_lu)) / u_lu(i_lu, i_lu)
    end do
end subroutine



subroutine data_processing() 
    IMPLICIT NONE
    call geo_process()
    call alpha_calculate()
    call q_calculate()
    !call channel_approximation()
   ! call alpha_calculation()
end subroutine

end  module