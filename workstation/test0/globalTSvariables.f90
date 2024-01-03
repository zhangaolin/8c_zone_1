module GlobalTSVariables
use GlobalTSConstants
use, intrinsic :: iso_fortran_env
IMPLICIT NONE
integer :: num_r,num_z,num_theta,bc_l,bc_r,bc_s,i,j,k
integer,allocatable :: cfg(:,:,:)
!integer,parameter :: 
real(TS_DOUBLE) :: dp
real(TS_DOUBLE),allocatable :: sub_r(:),sub_z(:),sub_theta(:),tg(:,:,:),h(:,:,:),alpha(:,:,:)
real(TS_DOUBLE),allocatable :: r(:),theta(:),z(:),lambda(:,:,:),power_den(:,:,:),ts_ave(:,:,:),q_vir(:,:,:)
!real(TS_DOUBLE),parameter :: 

end module