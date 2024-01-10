module GlobalTSVariables
use GlobalTSConstants
use, intrinsic :: iso_fortran_env
IMPLICIT NONE
integer :: num_r,num_z,num_theta
integer :: bc_u,bc_d,bc_s
integer :: i,j,k,any_int,i_lu,j_lu, k_lu,n_lu,num_ch,i_ch
integer,allocatable :: cfg(:,:,:)
!integer,parameter :: 
real(TS_DOUBLE) :: dp,any_real
real(TS_DOUBLE),allocatable :: r(:,:,:),sub_r(:),sub_z(:),sub_theta(:),volume_node(:,:,:)
real(TS_DOUBLE),allocatable :: tg(:,:,:),ts_ave(:,:,:),h(:,:,:),alpha(:,:,:)
real(TS_DOUBLE),allocatable :: lambda(:,:,:),power_den(:,:,:),q_vir(:,:,:)
real(TS_DOUBLE),allocatable :: a_ch(:),b_ch(:),c_ch(:),d_ch(:)&
&,x_ch(:),y_ch(:),u_ch(:),l_ch(:)
real(TS_DOUBLE),allocatable :: a_lu(:,:),b_lu(:),l_lu(:,:),u_lu(:,:),y_lu(:),x_lu(:)

!real(TS_DOUBLE),parameter :: 

end module