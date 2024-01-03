module driverthinput
use GlobalTSConstants
use GlobalTSVariables
use, intrinsic :: iso_fortran_env
Implicit None
contains

subroutine read_TH_Input()
character(100) :: aline,aword,keyword,filename
integer :: iunit,io_error
open(unit=8,file="input.dat")
do
    read(unit=8,fmt='(A)',iostat=io_error) aline
    if(aline(1:100).eq. ' ') then
        cycle
    end if  
    if (io_error == IOSTAT_END)EXIT 
    read(unit=aline, fmt=*, iostat=io_error) aword
    select case(trim(adjustl(aword)))
        case('PebbleDiameter')
            read(unit=aline,fmt=*,iostat=io_error) keyword,dp
        case('Axial_node_num')
            read(unit=aline,fmt=*,iostat=io_error) keyword,num_r
        case('Radial_node_num')
            read(unit=aline,fmt=*,iostat=io_error) keyword,num_z
        case('Circumferential_node_num')
            read(unit=aline,fmt=*,iostat=io_error) keyword,num_theta
        case('Axial_node_size')
            allocate(sub_r(num_r))
            read(unit=aline,fmt=*,iostat=io_error) keyword,sub_r(1:num_r)
        case('Radial_node_size')
            allocate(sub_z(num_z))
            read(unit=aline,fmt=*,iostat=io_error) keyword,sub_z(1:num_z)
        case('Circumferential_node_size')
            allocate(sub_theta(num_theta))
            read(unit=aline,fmt=*,iostat=io_error) keyword,sub_theta(1:num_theta)
        case('bc_r_left')
            read(unit=aline,fmt=*,iostat=io_error) keyword,bc_l
        case('bc_r_right')
            read(unit=aline,fmt=*,iostat=io_error) keyword,bc_r
        case('bc_r_side')
            read(unit=aline,fmt=*,iostat=io_error) keyword,bc_s

        case('config')
            allocate(cfg(num_r,num_z,num_theta))
            do j=1,num_z
                do i = 1,num_theta
                    if(i .eq. 1) then
                        read(unit=aline, fmt=*, iostat=io_error) keyword,cfg(1:num_r,1,j)
                    else
                        read(unit=8, fmt=*, iostat=io_error) cfg(1:num_r, i,j)
                    end if
                end do 
            end do
    end select
end do
filename = "output.dat"
open(unit=iunit, file=filename)
end subroutine

end module
!write(iunit,*)"output:"