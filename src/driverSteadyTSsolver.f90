module DriverSteadyTSSolver
use GlobalTSVariables
use GlobalTSConstants
IMPLICIT NONE
contains 

subroutine initialize() !初始化，给定初值

end subroutine

subroutine leakage_calculation() !计算横向泄漏
    IMPLICIT NONE
    real(TS_KDUBLE),allocatable :: z_leakage(:,:,:),r_leakage(:,:,:),theta_leakage(:,:,:)
    allocate(z_leakage(z_num,r_num,theta_num))
    allocate(r_leakage(z_num,r_num,theta_num))
    allocate(theta_leakage(z_num,r_num,theta_num))
    do i=1,z_num
        do j=1,r_num
            do k=1,theta_num
                z_leakage(i,j,k)=!
                r_leakage(i,j,k)=!
                theta_leakage(i,j,k)=!
            end do
        end do
    end do
end subroutine

subroutine lambda_calculation()!计算每种材料的导热系数
    IMPLICIT NONE
    real(TS_KDUBLE),allocatable ::met(:,:,:)
    allocate(met(z_num,r_num,theta_num))
    do i=1,z_num
        do j=1,r_num
            do k=1,theta_num
                met(i,j,k)=composition_solidmodule(node_config(i,j,k))%method
                !具体的计算方法需要在确定流道近似的方案后填写
                if(met(i,j,k) == 1.)then
                    heat_conductivity=!
                else if(met(i,j,k) == 2.)then
                    heat_conductivity=!
                else if(met(i,j,k) == 3.)then
                    heat_conductivity=!
                end if
                lambda_solidmodule(i,j,k)=heat_conductivity
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