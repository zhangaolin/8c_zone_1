module DriverSteadyTSSolver
use GlobalTSVariables
use GlobalTSConstants
use driver_input_test
use DataProcess
IMPLICIT NONE
!integer,allocatable :: 
!integer,parameter :: 
real(TS_DOUBLE) :: xi
real(TS_DOUBLE),allocatable :: j_radial_l(:,:,:),j_radial_r(:,:,:),j_circum_l(:,:,:),&
&j_circum_r(:,:,:),j_axial_l(:,:,:), j_axial_r(:,:,:)

real(TS_DOUBLE),allocatable :: t_radial_l(:,:,:),t_radial_r(:,:,:),t_circum_l(:,:,:),&
&t_circum_r(:,:,:),t_axial_l(:,:,:), t_axial_r(:,:,:)

real(TS_DOUBLE),allocatable :: z_leakage(:,:,:),r_leakage(:,:,:),theta_leakage(:,:,:)

real(TS_DOUBLE),allocatable :: r_leakage_temp(:,:,:),theta_leakage_temp(:,:,:),&
&z_leakage_temp(:,:,:)

real(TS_DOUBLE),allocatable :: a(:),b(:),c(:) 

!real(TS_DOUBLE),allocatable :: coeff_a_z1(:,:,:),coeff_a_z2(:,:,:),coeff_a_z3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_b_z1(:,:,:),coeff_b_z2(:,:,:),coeff_b_z3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_c_z1(:,:,:),coeff_c_z2(:,:,:),coeff_c_z3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_d_z1(:,:,:),coeff_d_z2(:,:,:),coeff_d_z3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_e_z1(:,:,:),coeff_e_z2(:,:,:),coeff_e_z3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_f_z1(:,:,:),coeff_f_z2(:,:,:),coeff_f_z3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_g_z1(:,:,:),coeff_g_z2(:,:,:),coeff_g_z3(:,:,:)

!real(TS_DOUBLE),allocatable :: coeff_a_r1(:,:,:),coeff_a_r2(:,:,:),coeff_a_r3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_b_r1(:,:,:),coeff_b_r2(:,:,:),coeff_b_r3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_c_r1(:,:,:),coeff_c_r2(:,:,:),coeff_c_r3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_d_r1(:,:,:),coeff_d_r2(:,:,:),coeff_d_r3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_e_r1(:,:,:),coeff_e_r2(:,:,:),coeff_e_r3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_f_r1(:,:,:),coeff_f_r2(:,:,:),coeff_f_r3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_g_r1(:,:,:),coeff_g_r2(:,:,:),coeff_g_r3(:,:,:)

!real(TS_DOUBLE),allocatable :: coeff_a_th1(:,:,:),coeff_a_th2(:,:,:),coeff_a_th3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_b_th1(:,:,:),coeff_b_th2(:,:,:),coeff_b_th3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_c_th1(:,:,:),coeff_c_th2(:,:,:),coeff_c_th3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_d_th1(:,:,:),coeff_d_th2(:,:,:),coeff_d_th3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_e_th1(:,:,:),coeff_e_th2(:,:,:),coeff_e_th3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_f_th1(:,:,:),coeff_f_th2(:,:,:),coeff_f_th3(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_g_th1(:,:,:),coeff_g_th2(:,:,:),coeff_g_th3(:,:,:)

real(TS_DOUBLE),allocatable :: heatflux_temp1(:,:,:)
real(TS_DOUBLE),allocatable :: coeff_q_z(:,:,:),coeff_q_r(:,:,:),coeff_q_th(:,:,:)

integer :: timestep_solid
real(TS_DOUBLE),allocatable :: error_Tsolid(:,:,:),old_T_solid(:,:,:)
!real(TS_DOUBLE),parameter :: 
contains 

subroutine initialize() !initialize
    allocate(j_axial_l(num_r,num_theta,num_z))
    allocate(j_axial_r(num_r,num_theta,num_z))
    allocate(j_radial_l(num_r,num_theta,num_z))
    allocate(j_radial_r(num_r,num_theta,num_z))
    allocate(j_circum_l(num_r,num_theta,num_z))
    allocate(j_circum_r(num_r,num_theta,num_z))

    allocate(t_axial_l(num_r,num_theta,num_z))
    allocate(t_axial_r(num_r,num_theta,num_z))
    allocate(t_radial_l(num_r,num_theta,num_z))
    allocate(t_radial_r(num_r,num_theta,num_z))
    allocate(t_circum_l(num_r,num_theta,num_z))
    allocate(t_circum_r(num_r,num_theta,num_z))

    allocate(z_leakage(num_r,num_theta,num_z))
    allocate(r_leakage(num_r,num_theta,num_z))
    allocate(theta_leakage(num_r,num_theta,num_z))

    allocate(z_leakage_temp(num_r,num_theta,num_z))
    allocate(r_leakage_temp(num_r,num_theta,num_z))
    allocate(theta_leakage_temp(num_r,num_theta,num_z))

    ! allocate(coeff_a_r1(num_r,num_theta,num_z))
    ! allocate(coeff_a_r2(num_r,num_theta,num_z))
    ! allocate(coeff_a_r3(num_r,num_theta,num_z))
    allocate(coeff_b_r1(num_r,num_theta,num_z))
    allocate(coeff_b_r2(num_r,num_theta,num_z))
    allocate(coeff_b_r3(num_r,num_theta,num_z))
    allocate(coeff_c_r1(num_r,num_theta,num_z))
    allocate(coeff_c_r2(num_r,num_theta,num_z))
    allocate(coeff_c_r3(num_r,num_theta,num_z))
    allocate(coeff_d_r1(num_r,num_theta,num_z))
    allocate(coeff_d_r2(num_r,num_theta,num_z))
    allocate(coeff_d_r3(num_r,num_theta,num_z))
    allocate(coeff_e_r1(num_r,num_theta,num_z))
    allocate(coeff_e_r2(num_r,num_theta,num_z))
    allocate(coeff_e_r3(num_r,num_theta,num_z))
    allocate(coeff_f_r1(num_r,num_theta,num_z))
    allocate(coeff_f_r2(num_r,num_theta,num_z))
    allocate(coeff_f_r3(num_r,num_theta,num_z))
    allocate(coeff_g_r1(num_r,num_theta,num_z))
    allocate(coeff_g_r2(num_r,num_theta,num_z))
    allocate(coeff_g_r3(num_r,num_theta,num_z))

    ! allocate(coeff_a_z1(num_r,num_theta,num_z))
    ! allocate(coeff_a_z2(num_r,num_theta,num_z))
    ! allocate(coeff_a_z3(num_r,num_theta,num_z))
    allocate(coeff_b_z1(num_r,num_theta,num_z))
    allocate(coeff_b_z2(num_r,num_theta,num_z))
    allocate(coeff_b_z3(num_r,num_theta,num_z))
    allocate(coeff_c_z1(num_r,num_theta,num_z))
    allocate(coeff_c_z2(num_r,num_theta,num_z))
    allocate(coeff_c_z3(num_r,num_theta,num_z))
    allocate(coeff_d_z1(num_r,num_theta,num_z))
    allocate(coeff_d_z2(num_r,num_theta,num_z))
    allocate(coeff_d_z3(num_r,num_theta,num_z))
    allocate(coeff_e_z1(num_r,num_theta,num_z))
    allocate(coeff_e_z2(num_r,num_theta,num_z))
    allocate(coeff_e_z3(num_r,num_theta,num_z))
    allocate(coeff_f_z1(num_r,num_theta,num_z))
    allocate(coeff_f_z2(num_r,num_theta,num_z))
    allocate(coeff_f_z3(num_r,num_theta,num_z))
    allocate(coeff_g_z1(num_r,num_theta,num_z))
    allocate(coeff_g_z2(num_r,num_theta,num_z))
    allocate(coeff_g_z3(num_r,num_theta,num_z))

    ! allocate(coeff_a_th1(num_r,num_theta,num_z))
    ! allocate(coeff_a_th2(num_r,num_theta,num_z))
    ! allocate(coeff_a_th3(num_r,num_theta,num_z))
    allocate(coeff_b_th1(num_r,num_theta,num_z))
    allocate(coeff_b_th2(num_r,num_theta,num_z))
    allocate(coeff_b_th3(num_r,num_theta,num_z))
    allocate(coeff_c_th1(num_r,num_theta,num_z))
    allocate(coeff_c_th2(num_r,num_theta,num_z))
    allocate(coeff_c_th3(num_r,num_theta,num_z))
    allocate(coeff_d_th1(num_r,num_theta,num_z))
    allocate(coeff_d_th2(num_r,num_theta,num_z))
    allocate(coeff_d_th3(num_r,num_theta,num_z))
    allocate(coeff_e_th1(num_r,num_theta,num_z))
    allocate(coeff_e_th2(num_r,num_theta,num_z))
    allocate(coeff_e_th3(num_r,num_theta,num_z))
    allocate(coeff_f_th1(num_r,num_theta,num_z))
    allocate(coeff_f_th2(num_r,num_theta,num_z))
    allocate(coeff_f_th3(num_r,num_theta,num_z))
    allocate(coeff_g_th1(num_r,num_theta,num_z))
    allocate(coeff_g_th2(num_r,num_theta,num_z))
    allocate(coeff_g_th3(num_r,num_theta,num_z))

    allocate(coeff_q_r(num_r,num_theta,num_z))
    allocate(coeff_q_z(num_r,num_theta,num_z))
    allocate(coeff_q_th(num_r,num_theta,num_z))

    allocate(lambda(num_r,num_theta,num_z))
    allocate(error_Tsolid(num_r,num_theta,num_z))
    allocate(old_T_solid(num_r,num_theta,num_z))
    allocate(heatflux_temp1(num_r,num_theta,num_z))
    allocate(ts_ave(num_r,num_theta,num_z))
    allocate(a(num_r),b(num_theta),c(num_z))
    do k=1,num_z
        do i=1,num_r
            do j=1,num_theta
                a(i)=0.5*sub_r(i)
                b(j)=0.5*sub_theta(j)
                c(k)=0.5*sub_z(k)
            end do
        end do
    end do
    j_axial_l =1
    j_axial_r =1
    j_radial_l =1
    j_radial_r =1
    j_circum_l =1
    j_circum_r =1

    coeff_b_r1 =1
    coeff_b_r2 =1
    coeff_b_r3 =1
    coeff_c_r1 =1
    coeff_c_r2 =1
    coeff_c_r3 =1
    coeff_d_r1 =1
    coeff_d_r2 =1
    coeff_d_r3 =1
    coeff_e_r1 =1
    coeff_e_r2 =1
    coeff_e_r3 =1
    coeff_f_r1 =1
    coeff_f_r2 =1
    coeff_f_r3 =1
    coeff_g_r1 =1
    coeff_g_r2 =1
    coeff_g_r3 =1

    coeff_b_z1 =1
    coeff_b_z2 =1
    coeff_b_z3 =1
    coeff_c_z1 =1
    coeff_c_z2 =1
    coeff_c_z3 =1
    coeff_d_z1 =1
    coeff_d_z2 =1
    coeff_d_z3 =1
    coeff_e_z1 =1
    coeff_e_z2 =1
    coeff_e_z3 =1
    coeff_f_z1 =1
    coeff_f_z2 =1
    coeff_f_z3 =1
    coeff_g_z1 =1
    coeff_g_z2 =1
    coeff_g_z3 =1

    coeff_b_th1 =1
    coeff_b_th2 =1
    coeff_b_th3 =1
    coeff_c_th1 =1
    coeff_c_th2 =1
    coeff_c_th3 =1
    coeff_d_th1 =1
    coeff_d_th2 =1
    coeff_d_th3 =1
    coeff_e_th1 =1
    coeff_e_th2 =1
    coeff_e_th3 =1
    coeff_f_th1 =1
    coeff_f_th2 =1
    coeff_f_th3 =1
    coeff_g_th1 =1
    coeff_g_th2 =1
    coeff_g_th3 =1

    coeff_q_r =0
    coeff_q_th =0
    coeff_q_z =0
end subroutine

subroutine leakage_calculation() !计算横向泄漏
    IMPLICIT NONE

    do k=1,num_z
        do i=1,num_r
            do j=1,num_theta
                z_leakage_temp(i,j,k)=(j_axial_r(i,j,k)-j_axial_l(i,j,k))/(2*c(k))
                r_leakage_temp(i,j,k)=((r(i,j,k)+a(i))*j_radial_l(i,j,k)-(r(i,j,k)-a(i))*j_radial_r(i,j,k))/(2*a(k)*r(i,j,k))
                theta_leakage_temp(i,j,k)=(j_circum_r(i,j,k)-j_circum_l(i,j,k))/(2*b(k))

                z_leakage(i,j,k)=r_leakage_temp(i,j,k)+theta_leakage_temp(i,j,k)
                r_leakage(i,j,k)=z_leakage_temp(i,j,k)+theta_leakage_temp(i,j,k)
                theta_leakage(i,j,k)=r_leakage_temp(i,j,k)+z_leakage_temp(i,j,k)
            end do
        end do
    end do
end subroutine

subroutine lambda_calculation()!计算每种材料的导热系数
    IMPLICIT NONE
    do k=1,num_z
        do i=1,num_r
            do j=1,num_theta
                if(cfg(i,j,k)==1.)then
                    lambda(i,j,k)=75
                else if(cfg(i,j,k)==2.)then
                    lambda(i,j,k)=150
                end if
            end do
        end do
    end do
end subroutine

subroutine t_bound_calculation()!calculate boundary temperature for every node
    IMPLICIT NONE
    do k=1,num_z
        do i=1,num_r
            do j=1,num_theta
                t_axial_r(i,j,k)=coeff_f_z1(i,j,k)*j_axial_r(i,j,k)+coeff_f_z2(i,j,k)*&
                &j_axial_l(i,j,k)+coeff_f_z3(i,j,k)*coeff_q_z(i,j,k)
                t_axial_l(i,j,k)=coeff_g_z1(i,j,k)*j_axial_r(i,j,k)+&
                &coeff_g_z2(i,j,k)*j_axial_l(i,j,k)+coeff_g_z3(i,j,k)*coeff_q_z(i,j,k)

                t_radial_r(i,j,k)=coeff_f_r1(i,j,k)*j_axial_r(i,j,k)+&
                &coeff_f_r2(i,j,k)*j_axial_l(i,j,k)+coeff_f_r3(i,j,k)*coeff_q_r(i,j,k)
                t_radial_l(i,j,k)=coeff_g_r1(i,j,k)*j_axial_r(i,j,k)+&
                &coeff_g_r2(i,j,k)*j_axial_l(i,j,k)+coeff_g_r3(i,j,k)*coeff_q_r(i,j,k)

                t_circum_r(i,j,k)=coeff_f_th1(i,j,k)*j_axial_r(i,j,k)+&
                &coeff_f_th2(i,j,k)*j_axial_l(i,j,k)+coeff_f_th3(i,j,k)*coeff_q_th(i,j,k)
                t_circum_l(i,j,k)=coeff_g_th1(i,j,k)*j_axial_r(i,j,k)+&
                &coeff_g_th2(i,j,k)*j_axial_l(i,j,k)+coeff_g_th3(i,j,k)*coeff_q_th(i,j,k)
            end do
        end do
    end do
end subroutine




subroutine heatflux_update_z() !更新边界热流
    implicit none

    !calculation of coeffs
        do k=1,num_z
            do i=1,num_r
                do j=1,num_theta
                    coeff_b_z1(i,j,k)=0.5d0
                    coeff_b_z2(i,j,k)=-0.5d0
                    coeff_b_z3(i,j,k)=0.0d0

                    coeff_c_z1(i,j,k)=0.5*alpha(i,j,k)/(alpha(i,j,k)+3*lambda(i,j,k)/(c(k)**2))
                    coeff_c_z2(i,j,k)=0.5*alpha(i,j,k)/(alpha(i,j,k)+3*lambda(i,j,k)/(c(k)**2))
                    coeff_c_z3(i,j,k)=-1.0/(alpha(i,j,k)+3*lambda(i,j,k)/(c(k)**2))


                    coeff_d_z1(i,j,k)=-1*lambda(i,j,k)/c(k)*(coeff_b_z1(i,j,k)+3*coeff_c_z1(i,j,k))
                    coeff_d_z2(i,j,k)=-1*lambda(i,j,k)/c(k)*(coeff_b_z2(i,j,k)+3*coeff_c_z2(i,j,k))
                    coeff_d_z3(i,j,k)=-1*lambda(i,j,k)/c(k)*(coeff_b_z3(i,j,k)+3*coeff_c_z3(i,j,k))

                    coeff_e_z1(i,j,k)=-1*lambda(i,j,k)/c(k)*(coeff_b_z1(i,j,k)-3*coeff_c_z1(i,j,k))
                    coeff_e_z2(i,j,k)=-1*lambda(i,j,k)/c(k)*(coeff_b_z2(i,j,k)-3*coeff_c_z2(i,j,k))
                    coeff_e_z3(i,j,k)=-1*lambda(i,j,k)/c(k)*(coeff_b_z3(i,j,k)-3*coeff_c_z3(i,j,k))

                    heatflux_temp1=coeff_d_z1(i,j,k)*coeff_e_z2(i,j,k)-coeff_d_z2(i,j,k)*coeff_e_z1(i,j,k)

                    coeff_f_z1(i,j,k)=coeff_e_z2(i,j,k)/heatflux_temp1(i,j,k)
                    coeff_f_z2(i,j,k)=-1*coeff_d_z2(i,j,k)/heatflux_temp1(i,j,k)
                    coeff_f_z3(i,j,k)=(coeff_d_z2(i,j,k)*coeff_e_z3(i,j,k)-coeff_d_z3(i,j,k)*&
                    &coeff_e_z2(i,j,k))/heatflux_temp1(i,j,k)

                    coeff_g_z1(i,j,k)=-1*coeff_e_z1(i,j,k)/heatflux_temp1(i,j,k)
                    coeff_g_z2(i,j,k)=coeff_d_z1(i,j,k)/heatflux_temp1(i,j,k)
                    coeff_g_z3(i,j,k)=(coeff_d_z3(i,j,k)*coeff_e_z1(i,j,k)-coeff_d_z1(i,j,k)*&
                    &coeff_e_z3(i,j,k))/heatflux_temp1(i,j,k)

                    coeff_q_z(i,j,k)=q_vir(i,j,k)-z_leakage(i,j,k)

                end do
            end do
        end do
    !matrix filling
        num_ch=num_z+1
            !allocate chase parameters
        allocate(a_ch(num_ch), b_ch(num_ch), c_ch(num_ch), d_ch(num_ch))
        allocate(u_ch(num_ch), l_ch(num_ch), x_ch(num_ch), y_ch(num_ch))
        do i=1,num_r
            do j=1,num_theta
                do i_ch=2,num_ch-1
                    a_ch(i_ch)=coeff_f_z2(i,j,i_ch)
                    b_ch(i_ch)=coeff_f_z1(i,j,i_ch)-coeff_g_z2(i,j,i_ch+1)
                    c_ch(i_ch)=-1*coeff_g_z1(i,j,i_ch+1)
                    d_ch(i_ch)=coeff_g_z3(i,j,i_ch+1)*coeff_q_z(i,j,i_ch+1)-coeff_f_z3(i,j,i_ch)*coeff_q_z(i,j,i_ch)
                end do
                !boundary condition
                    a_ch(1)=0
                    b_ch(1)=1
                    c_ch(1)=0
                    d_ch(1)=0
                    a_ch(num_ch)=coeff_f_r2(i,j,1)/303.15
                    b_ch(num_ch)=coeff_f_r1(i,j,1)/303.15-1/(1000*303.15)
                    c_ch(num_ch)=0
                    d_ch(num_ch)=1-coeff_f_r3(i,j,num_z)/303.15*coeff_q_z(i,j,num_z)
                call chase()
                do k=1,num_z
                    j_axial_r(i,j,k)=x_ch(k+1)
                    j_axial_l(i,j,k)=x_ch(k)
                end do
            end do
        end do

    !deallocate chase parameters
        deallocate(a_ch)
        deallocate(b_ch)
        deallocate(c_ch)
        deallocate(d_ch)
        deallocate(l_ch)
        deallocate(u_ch)
        deallocate(x_ch)
        deallocate(y_ch)
end subroutine

subroutine heatflux_update_r() !更新边界热流
    implicit none

    !calculation of coeffs
        do k=1,num_z
            do i=1,num_r
                do j=1,num_theta
                    coeff_b_r1(i,j,k)=1.0d0/6.0d0
                    coeff_b_r2(i,j,k)=1.0d0/6.0d0
                    coeff_b_r3(i,j,k)=0.0d0

                    coeff_c_r1(i,j,k)=(3*alpha(i,j,k)*r(i,j,k)*(a(i)**2)+(alpha(i,j,k)*a(i)**3-3*a(i)*lambda(i,j,k)))&
                    &/(6*alpha(i,j,k)*r(i,j,k)*(a(i)**2)+18*lambda(i,j,k)*r(i,j,k))
                    coeff_c_r2(i,j,k)=(3*alpha(i,j,k)*r(i,j,k)*(a(i)**2)-(alpha(i,j,k)*a(i)**3-3*a(i)*lambda(i,j,k)))&
                    &/(6*alpha(i,j,k)*r(i,j,k)*(a(i)**2)+18*lambda(i,j,k)*r(i,j,k))
                    coeff_c_r3(i,j,k)=-1.0d0*(a(i)**2)/(alpha(i,j,k)*a(i)**2+3*lambda(i,j,k))


                    coeff_d_r1(i,j,k)=-3.0d0*lambda(i,j,k)/a(i)*(coeff_b_r1(i,j,k)+coeff_c_r1(i,j,k))
                    coeff_d_r2(i,j,k)=-3.0d0*lambda(i,j,k)/a(i)*(coeff_b_r2(i,j,k)+coeff_c_r2(i,j,k))
                    coeff_d_r3(i,j,k)=-3.0d0*lambda(i,j,k)/a(i)*(coeff_b_r3(i,j,k)+coeff_c_r3(i,j,k))

                    coeff_e_r1(i,j,k)=-3.0d0*lambda(i,j,k)/a(i)*(coeff_b_r1(i,j,k)-coeff_c_r1(i,j,k))
                    coeff_e_r2(i,j,k)=-3.0d0*lambda(i,j,k)/a(i)*(coeff_b_r2(i,j,k)-coeff_c_r2(i,j,k))
                    coeff_e_r3(i,j,k)=-3.0d0*lambda(i,j,k)/a(i)*(coeff_b_r3(i,j,k)-coeff_c_r3(i,j,k))

                    heatflux_temp1=coeff_d_r1(i,j,k)*coeff_e_r2(i,j,k)-coeff_d_r2(i,j,k)*coeff_e_r1(i,j,k)

                    coeff_f_r1(i,j,k)=coeff_e_r2(i,j,k)/heatflux_temp1(i,j,k)
                    coeff_f_r2(i,j,k)=-1*coeff_d_r2(i,j,k)/heatflux_temp1(i,j,k)
                    coeff_f_r3(i,j,k)=(coeff_d_r2(i,j,k)*coeff_e_r3(i,j,k)-coeff_d_r3(i,j,k)*&
                    &coeff_e_r2(i,j,k))/heatflux_temp1(i,j,k)

                    coeff_g_r1(i,j,k)=-1*coeff_e_r1(i,j,k)/heatflux_temp1(i,j,k)
                    coeff_g_r2(i,j,k)=coeff_d_r1(i,j,k)/heatflux_temp1(i,j,k)
                    coeff_g_r3(i,j,k)=(coeff_d_r3(i,j,k)*coeff_e_r1(i,j,k)-coeff_d_r1(i,j,k)*&
                    &coeff_e_r3(i,j,k))/heatflux_temp1(i,j,k)

                    coeff_q_r(i,j,k)=q_vir(i,j,k)-r_leakage(i,j,k)
                end do
            end do
        end do
    !matrix filling
        num_ch=num_r+1
            !allocate chase parameters
        allocate(a_ch(num_ch), b_ch(num_ch), c_ch(num_ch), d_ch(num_ch))
        allocate(u_ch(num_ch), l_ch(num_ch), x_ch(num_ch), y_ch(num_ch))
        do k=1,num_z
            do j=1,num_theta
                do i_ch=2,num_ch-1
                    a_ch(i_ch)=coeff_f_r2(i_ch,j,k)
                    b_ch(i_ch)=coeff_f_r1(i_ch,j,k)-coeff_g_r2(i_ch+1,j,k)
                    c_ch(i_ch)=-1*coeff_g_r1(i_ch+1,j,k)
                    d_ch(i_ch)=coeff_g_r3(i_ch+1,j,k)*coeff_q_r(i_ch+1,j,k)-coeff_f_r3(i_ch,j,k)*coeff_q_r(i_ch,j,k)
                end do
                !boundary condition
                    a_ch(1)=0
                    b_ch(1)=1
                    c_ch(1)=0
                    d_ch(1)=0
                    a_ch(num_ch)=0
                    b_ch(num_ch)=1
                    c_ch(num_ch)=0
                    d_ch(num_ch)=0
                call chase()
                do i=1,num_r
                    j_radial_r(i,j,k)=x_ch(i+1)
                    j_radial_l(i,j,k)=x_ch(i)
                end do
            end do
        end do

    !deallocate chase parameters
        deallocate(a_ch)
        deallocate(b_ch)
        deallocate(c_ch)
        deallocate(d_ch)
        deallocate(l_ch)
        deallocate(u_ch)
        deallocate(x_ch)
        deallocate(y_ch)
end subroutine


subroutine heatflux_update_theta() !更新边界热流
    implicit none

    !calculation of coeffs
        do k=1,num_z
            do i=1,num_r
                do j=1,num_theta
                    coeff_b_th1(i,j,k)=0.5d0
                    coeff_b_th2(i,j,k)=-0.5d0
                    coeff_b_th3(i,j,k)=0.0d0

                    coeff_c_th1(i,j,k)=0.5d0*alpha(i,j,k)/(alpha(i,j,k)+3/(r(i,j,k)**2)*lambda(i,j,k)/(c(k)**2))
                    coeff_c_th2(i,j,k)=0.5d0*alpha(i,j,k)/(alpha(i,j,k)+3/(r(i,j,k)**2)*lambda(i,j,k)/(c(k)**2))
                    coeff_c_th3(i,j,k)=-1.0d0/(alpha(i,j,k)+3*lambda(i,j,k)/((r(i,j,k)*c(k))**2))


                    coeff_d_th1(i,j,k)=-1*lambda(i,j,k)/c(k)*(coeff_b_th1(i,j,k)+3*coeff_c_th1(i,j,k))
                    coeff_d_th2(i,j,k)=-1*lambda(i,j,k)/c(k)*(coeff_b_th2(i,j,k)+3*coeff_c_th2(i,j,k))
                    coeff_d_th3(i,j,k)=-1*lambda(i,j,k)/c(k)*(coeff_b_th3(i,j,k)+3*coeff_c_th3(i,j,k))

                    coeff_e_th1(i,j,k)=-1*lambda(i,j,k)/c(k)*(coeff_b_th1(i,j,k)-3*coeff_c_th1(i,j,k))
                    coeff_e_th2(i,j,k)=-1*lambda(i,j,k)/c(k)*(coeff_b_th2(i,j,k)-3*coeff_c_th2(i,j,k))
                    coeff_e_th3(i,j,k)=-1*lambda(i,j,k)/c(k)*(coeff_b_th3(i,j,k)-3*coeff_c_th3(i,j,k))

                    heatflux_temp1=coeff_d_th1(i,j,k)*coeff_e_th2(i,j,k)-coeff_d_th2(i,j,k)*coeff_e_th1(i,j,k)

                    coeff_f_th1(i,j,k)=coeff_e_th2(i,j,k)/heatflux_temp1(i,j,k)
                    coeff_f_th2(i,j,k)=-1*coeff_d_th2(i,j,k)/heatflux_temp1(i,j,k)
                    coeff_f_th3(i,j,k)=(coeff_d_th2(i,j,k)*coeff_e_th3(i,j,k)-coeff_d_th3(i,j,k)*&
                    &coeff_e_th2(i,j,k))/heatflux_temp1(i,j,k)

                    coeff_g_th1(i,j,k)=-1*coeff_e_th1(i,j,k)/heatflux_temp1(i,j,k)
                    coeff_g_th2(i,j,k)=coeff_d_th1(i,j,k)/heatflux_temp1(i,j,k)
                    coeff_g_th3(i,j,k)=(coeff_d_th3(i,j,k)*coeff_e_th1(i,j,k)-coeff_d_th1(i,j,k)*&
                    &coeff_e_th3(i,j,k))/heatflux_temp1(i,j,k)

                    coeff_q_th(i,j,k)=q_vir(i,j,k)-theta_leakage(i,j,k)
                end do
            end do
        end do
    !matrix filling

    !!matrix filling
        ! num_ch=num_theta+1
        !     !allocate chase parameters
        ! allocate(a_ch(num_ch), b_ch(num_ch), c_ch(num_ch), d_ch(num_ch))
        ! allocate(u_ch(num_ch), l_ch(num_ch), x_ch(num_ch), y_ch(num_ch))
        ! do i=1,num_r
        !     do k=1,num_z
        !         do i_ch=2,num_ch-1
        !             a_ch(i_ch)=coeff_f_th2(i,i_ch,k)
        !             b_ch(i_ch)=coeff_f_th1(i,i_ch,k)-coeff_g_th2(i,i_ch+1,k)
        !             c_ch(i_ch)=-1*coeff_g_th1(i,i_ch+1,k)
        !             d_ch(i_ch)=coeff_g_th3(i,i_ch+1,k)*coeff_q_th(i,i_ch+1,k)-coeff_f_th3(i,i_ch,k)*coeff_q_th(i,j,i_ch)
        !         end do
        !         call chase()
        !         !boundary condition
        !             !=============================!
        !             !===========undone============!
        !             !=============================!
        !         do j=1,num_theta
        !             j_circum_r(i,j,k)=x_ch(j+1)
        !             j_circum_l(i,j,k)=x_ch(j)
        !         end do
        !     end do
        ! end do
    !!deallocate chase parameters
        ! deallocate(a_ch)
        ! deallocate(b_ch)
        ! deallocate(c_ch)
        ! deallocate(d_ch)
        ! deallocate(l_ch)
        ! deallocate(u_ch)
        ! deallocate(x_ch)
        ! deallocate(y_ch)
end subroutine




subroutine T_soild_average() !计算平均温度
    implicit none
    real(TS_DOUBLE),allocatable :: temp_1(:,:,:),temp_2(:,:,:)
    allocate(temp_1(num_r,num_theta,num_z))
    allocate(temp_2(num_r,num_theta,num_z))

    do k=1,num_z
        do j=1,num_theta
            do i=1,num_r
                temp_1(i,j,k)=alpha(i,j,k)+3*lambda(i,j,k)*(1/(a(i)*a(i))+1/((b(j)*r(i,j,k))**2)+1/(c(k)*c(k)))
                temp_2(i,j,k)=(a(i)*(t_radial_r(i,j,k)-t_radial_l(i,j,k))/6+r(i,j,k)*((0.5+a(i)/(6*r(i,j,k)))*t_radial_r(i,j,k)+&
                &(0.5-a(i)/(6*r(i,j,k)))*t_radial_l(i,j,k)))/(r(i,j,k)*a(i)*a(i))+&
                &(t_axial_l(i,j,k)+t_axial_r(i,j,k))/(2*c(k)*c(k))+(t_circum_l(i,j,k)+t_circum_r(i,j,k))/(2*((b(j)*r(i,j,k))**2))
                ts_ave(i,j,k)=(q_vir(i,j,k)+3*lambda(i,j,k)*temp_2(i,j,k))/temp_1(i,j,k)
            end do
        end do
    end do

end subroutine

subroutine driver_steady_ts_solver() !调用各个子例程求解温度
    IMPLICIT NONE

    call initialize()
    timestep_solid=0
    call data_processing()
    do while(any(error_Tsolid > convergence_limit) .and. timestep_solid<100.)
        call t_bound_calculation()
        call T_soild_average()
        !record old_t
            do k=1,num_z
                do j=1,num_theta
                    do i=1,num_r
                        old_T_solid(i,j,k)=ts_ave(i,j,k)
                    end do
                end do
            end do  
        call leakage_calculation()
        call lambda_calculation()
        call heatflux_update_z()
        call heatflux_update_r()
        !call heatflux_update_theta()
        !calculate error_t
            do k=1,num_z
                do j=1,num_theta
                    do i=1,num_r
                        error_Tsolid(i,j,k)=(ts_ave(i,j,k)-old_T_solid(i,j,k))/ts_ave(i,j,k)
                    end do
                end do
            end do
        timestep_solid=timestep_solid+1  
    end do

end subroutine

end  module