program test_solver    
implicit none


call chase_test()
end program

subroutine chase_test()
    implicit none
    integer :: num_ch,i_ch
    real(4),allocatable :: a_ch(:),b_ch(:),c_ch(:),d_ch(:)&
                                &,x_ch(:),y_ch(:),u_ch(:),l_ch(:)
    
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
        write(*,*) 'from unittest'
        write(*,*) 'test_chase passed!!!'
    else
        write(*,*) 'test_chase failed, x_ch=', x_ch
    end if
end subroutine



