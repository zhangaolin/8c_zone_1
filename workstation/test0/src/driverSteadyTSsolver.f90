module DriverSteadyTSSolver
use GlobalTSVariables
use GlobalTSConstants
use driver_input_test
use DataProcess
IMPLICIT NONE
!integer,allocatable :: 
!integer,parameter :: 
!real(TS_DOUBLE) :: 
real(TS_DOUBLE),allocatable :: j_radial_l(:,:,:),j_radial_r(:,:,:),j_circum_l(:,:,:),&
                                &j_circum_r(:,:,:),j_axial_l(:,:,:), j_axial_r(:,:,:)
real(TS_DOUBLE),allocatable :: t_radial_l(:,:,:),t_radial_r(:,:,:),t_circum_l(:,:,:),&
                                &t_circum_r(:,:,:),t_axial_l(:,:,:), t_axial_r(:,:,:)
real(TS_DOUBLE),allocatable :: z_leakage(:,:,:),r_leakage(:,:,:),theta_leakage(:,:,:)
real(TS_DOUBLE),allocatable :: a_ch(:),b_ch(:),c_ch(:),d_ch(:)&
                                &,x_ch(:),y_ch(:),u_ch(:),l_ch(:)
integer :: num_ch,i_ch
!real(TS_DOUBLE),parameter :: 
contains 

! subroutine initialize() !initialize
!     allocate(j_axial_l(num_r,num_theta,num_z))
!     allocate(j_axial_r(num_r,num_theta,num_z))
!     allocate(j_radial_l(num_r,num_theta,num_z))
!     allocate(j_radial_r(num_r,num_theta,num_z))
!     allocate(j_circum_l(num_r,num_theta,num_z))
!     allocate(j_circum_r(num_r,num_theta,num_z))
!     allocate(t_axial_l(num_r,num_theta,num_z))
!     allocate(t_axial_r(num_r,num_theta,num_z))
!     allocate(t_radial_l(num_r,num_theta,num_z))
!     allocate(t_radial_r(num_r,num_theta,num_z))
!     allocate(t_circum_l(num_r,num_theta,num_z))
!     allocate(t_circum_r(num_r,num_theta,num_z))
!     allocate(z_leakage(num_r,num_theta,num_z))
!     allocate(r_leakage(num_r,num_theta,num_z))
!     allocate(theta_leakage(num_r,num_theta,num_z))
!     allocate(lambda(num_r,num_theta,num_z))

! allocate(a_ch(num_ch))
! allocate(b_ch(num_ch))
! allocate(c_ch(num_ch))
! allocate(d_ch(num_ch))
! allocate(x_ch(num_ch))
! allocate(y_ch(num_ch))
! allocate(u_ch(num_ch))
! allocate(l_ch(num_ch))
!     do k=1,num_z
!         do j=1,num_theta
!             do i=1,num_r
!                 j_axial_l(i,j,k)=1
!                 j_axial_r(i,j,k)=1
!                 j_radial_l(i,j,k)=1
!                 j_radial_r(i,j,k)=1
!                 j_circum_l(i,j,k)=1
!                 j_circum_r(i,j,k)=1
!             end do
!         end do
!     end do
!     ! allocate(r(num_r))
!     ! r(1)=0.5*sub_r(1)
    
!     ! do j=2,num_r 
!     !     r(j)=r(1)
!     !     do i=2,j
!     !         r(j)=r(j)+sub_r(i)
!     !     end do
!     ! end do

! end subroutine

! subroutine leakage_calculation() !计算横向泄漏
!     IMPLICIT NONE
!     do k=1,num_z
!         do j=1,num_theta
!             do i=1,num_r
!                 z_leakage(i,j,k)=((r(i)+sub_r(i))*j_radial_r(i,j,k)-((r(i)-sub_r(i))*j_radial_l(i,j,k)))+(j_circum_r(i,j,k)-j_circum_l(i,j,k))/(2*sub_theta(j))
!                 r_leakage(i,j,k)=((r(i)+sub_r(i))*j_radial_r(i,j,k)-((r(i)-sub_r(i))*j_radial_l(i,j,k)))+(j_axial_r(i,j,k)-j_axial_l(i,j,k))/(2*sub_z(k))
!                 theta_leakage(i,j,k)=(j_circum_r(i,j,k)-j_circum_l(i,j,k))/(2*sub_theta(j))+(j_axial_r(i,j,k)-j_axial_l(i,j,k))/(2*sub_z(k))
!             end do
!         end do
!     end do
! end subroutine

! subroutine lambda_calculation()!计算每种材料的导热系数
!     IMPLICIT NONE
!     do k=1,num_z
!         do j=1,num_theta
!             do i=1,num_r
!                 if(cfg(i,j,k)==1.)then
!                     lambda(i,j,k)=75
!                 else if(cfg(i,j,k)==2.)then
!                     lambda(i,j,k)=150
!                 end if
!             end do
!         end do
!     end do
! end subroutine

subroutine chase_test()
    implicit none
    num_ch=5
    allocate(a_ch(num_ch), b_ch(num_ch), c_ch(num_ch), d_ch(num_ch))
    allocate(u_ch(num_ch), l_ch(num_ch), x_ch(num_ch), y_ch(num_ch))
    a_ch = (/0,2,-3,4,-5/)
    b_ch = (/1,3,4,7,6/)
    c_ch = (/2,1,2,1,0/)
    d_ch = (/5,9,2,19,-4/)
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
    if (all(x_ch == (/1,2,1,2,1/))) then
        write(*,*) 'test_chase passed!!!!'
    else
        write(*,*) 'test_chase failed, x_ch=', x_ch
    end if
end subroutine

! subroutine heatflux_update_z() !更新边界热流
! implicit none
! allocate(a_ch(num_ch))
! allocate(b_ch(num_ch))
! allocate(c_ch(num_ch))
! allocate(d_ch(num_ch))
! allocate(l_ch(num_ch))
! allocate(u_ch(num_ch))
! allocate(x_ch(num_ch))
! allocate(y_ch(num_ch))
! num_ch=num_z+1
! do k=2,num_ch-1
! a_ch(k)=coeff_f_z2(i,j,k)
! end do
! deallocate(a_ch)
! deallocate(b_ch)
! deallocate(c_ch)
! deallocate(d_ch)
! deallocate(l_ch)
! deallocate(u_ch)
! deallocate(x_ch)
! deallocate(y_ch)
! call chase()
! end subroutine

! subroutine heatflux_update_r() 
!     implicit none
! num_ch=num_r+1
! call chase()
! end subroutine

! subroutine heatflux_update_theta()
!     implicit none
!     num_ch=num_theta+1
!     call chase()
! end subroutine

! subroutine T_soild_average() !计算平均温度
!     implicit none
!     real(TS_DOUBLE),allocatable :: temp_1(:,:,:),temp_2(:,:,:)
!     allocate(temp_1(num_r,num_theta,num_z))
!     allocate(temp_2(num_r,num_theta,num_z))
!     do k=1,num_z
!         do j=1,num_theta
!             do i=1,num_r
!                 temp_1(i,j,k)=alpha(i,j,k)+3*lambda(i,j,k)*(1/(sub_r(i)*sub_r(i))+1/((sub_theta(k)*r(i))**2)+1/(sub_z(k)*sub_z(k)))
!                 temp_2(i,j,k)=(sub_r(i)*(t_radial_r(i,j,k)-t_radial_l(i,j,k))/6+r(i)*((0.5+sub_r(i)/(6*r(i)))*t_radial_r(i,j,k)+&
!                 &(0.5-sub_r(i)/(6*r(i)))*t_radial_l(i,j,k)))/(r(i)*sub_r(i)*sub_r(i))+&
!                 &(t_axial_l(i,j,k)+t_axial_r(i,j,k))/(2*sub_z(k)*sub_z(k))+(t_circum_l(i,j,k)+t_circum_r(i,j,k))/(2*((sub_theta(k)*r(i))**2))
!                 ts_ave(i,j,k)=(q_vir(i,j,k)+3*lambda(i,j,k)*temp_2(i,j,k))/temp_1(i,j,k)
!             end do
!         end do
!     end do
! end subroutine
! subroutine driver_steady_ts_solver() !调用各个子例程求解温度
! IMPLICIT NONE
! integer :: timestep_solid
! real(TS_DOUBLE),allocatable :: error_Tsolid(:,:,:),old_T_solid(:,:,:)
! allocate(error_Tsolid(num_r,num_theta,num_z))
! allocate(old_T_solid(num_r,num_theta,num_z))
! call initialize()
! timestep_solid=0
! do while(any(error_Tsolid > convergence_limit) .and. timestep_solid<1000.)
!     do k=1,num_z
!         do j=1,num_theta
!             do i=1,num_r
!                 old_T_solid(i,j,k)=ts_ave(i,j,k)
!             end do
!         end do
!     end do  
!     call leakage_calculation()
!     call lambda_calculation()
!     call heatflux_update_z()
!     call heatflux_update_r()
!     call heatflux_update_theta()
!     call T_soild_average()
!     do k=1,num_z
!         do j=1,num_theta
!             do i=1,num_r
!                 error_Tsolid(i,j,k)=(ts_ave(i,j,k)-old_T_solid(i,j,k))/ts_ave(i,j,k)
!             end do
!         end do
!     end do
!     timestep_solid=timestep_solid+1  
! end do

! end subroutine

end  module