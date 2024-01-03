module DriverSteadyTSSolver
use GlobalTSVariables
use GlobalTSConstants
IMPLICIT NONE
integer :: 
integer,allocatable :: 
!integer,parameter :: 
real(TS_DOUBLE) :: 
real(TS_DOUBLE),allocatable :: j_radial_l(:,:,:),j_radial_r(:,:,:),j_circum_l(:,:,:),j_circum_r(:,:,:)
real(TS_DOUBLE),allocatable :: j_axial_l(:,:,:), j_axial_r(:,:,:)
real(TS_DOUBLE),allocatable :: z_leakage(:,:,:),r_leakage(:,:,:),theta_leakage(:,:,:),lambda(:,:,:)
!real(TS_DOUBLE),parameter :: 
contains 

subroutine initialize() !初始化，给定初值
    allocate(j_axial_l(num_r,num_theta,num_z))
    allocate(j_axial_r(num_r,num_theta,num_z))
    allocate(j_radial_l(num_r,num_theta,num_z))
    allocate(j_radial_r(num_r,num_theta,num_z))
    allocate(j_circum_l(num_r,num_theta,num_z))
    allocate(j_circum_r(num_r,num_theta,num_z))
    do k=1,num_z
        do j=1,num_theta
            do i=1,num_r
                j_axial_l(i,j,k)=1
                j_axial_r(i,j,k)=1
                j_radial_l(i,j,k)=1
                j_radial_r(i,j,k)=1
                j_circum_l(i,j,k)=1
                j_circum_r(i,j,k)=1
            end do
        end do
    end do
    allocate(r(num_r))
    r(1)=0.5*sub_r(1)
    do j=2,num_r 
        r(j)=r(1)
        do i=2,j
            r(j)=r(j)+sub_r(i)
        end do
    end do

end subroutine

subroutine leakage_calculation() !计算横向泄漏
IMPLICIT NONE

allocate(z_leakage(num_r,num_theta,num_z))
allocate(r_leakage(num_r,num_theta,num_z))
allocate(theta_leakage(num_r,num_theta,num_z))
do k=1,num_z
    do j=1,num_theta
        do i=1,num_r
            z_leakage(i,j,k)=(())
            r_leakage(i,j,k)=
            theta_leakage(i,j,k)=
        end do
    end do
end do
end subroutine

subroutine lambda_calculation()!计算每种材料的导热系数
    IMPLICIT NONE
    allocate(lambda(num_r,num_theta,num_z))
    do k=1,num_z
        do j=1,num_theta
            do i=1,num_r


            end do
        end do
    end do
    end subroutine

subroutine heatflux_update_z() !更新边界热流
end subroutine
subroutine heatflux_update_r() 
end subroutine
subroutine heatflux_update_theta()
end subroutine

subroutine T_soild_average() !计算平均温度

end subroutine
subroutine driver_steady_ts_solver() !调用各个子例程求解温度
IMPLICIT NONE
integer :: timestep_solid
real(TS_KDUBLE),allocatable :: error_Tsolid(:,:,:),old_T_solid(:,:,:)
allocate(error_Tsolid(z_num,r_num,theta_num))
allocate(old_T_solid(z_num,r_num,theta_num))
call initialize()
do while(any(error_Tsolid) > convergence_limit .and. timestep_solid<1000.)
    do i=1,z_num
        do j=1,r_num
            do k=1,theta_num
                old_T_solid=T_solid
            end do
        end do
    end do  
    call leakage_calculation()
    call lambda_calculation()
    call heatflux_update_z()
    call heatflux_update_r()
    call heatflux_update_theta()
    call T_soild_average()
    do i=1,z_num
        do j=1,r_num
            do k=1,theta_num
                error_Tsolid(i,j,k)=(T_solid(i,j,k)-old_T_solid(i,j,k))/T_solid(i,j,k)
            end do
        end do
    end do  
end do

end subroutine

end  module