program test_solver    
    use DataProcess
    use GlobalTSVariables
    use GlobalTSConstants
    use driver_input_test
    use DriverTSoutput
    use DriverSteadyTSSolver
    implicit none

    call chase_test()
   ! call lu_test()
    !call lu_test_from_src()
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
    write(*,*) 'from unittest' 
    if (all(x_ch == (/1,2,1,2,1/))) then
        write(*,*) 'test_chase passed!!!x_ch=',x_ch
    else
        write(*,*) 'test_chase failed, x_ch=', x_ch
    end if
end subroutine

subroutine lu_test()
    implicit none
    integer :: n, i, j, k
    real(8), dimension(:,:), allocatable ::  L, U
    REAL(8) :: a(4,4)
    real(8),allocatable :: y_lu(:),x_lu(:),b_lu(:)
    n=4
    allocate(x_lu(n),y_lu(n),b_lu(n))
    data A /9,18,9,-27,&
    18,45,0,-45, &
    9,0,126,9,&
    -27,-45,9,135/
    b_lu = (/1,2,16,8/)
    allocate(L(n,n), U(n,n))
    l=0.
    u=0.
    y_lu=0.
    x_lu=0.
    do j = 1, n
        U(j,1) = A(j,1)
    end do
    do j=2,n
        L(1,j) = A(1,j) / U(1,1)
    end do
    do i = 2, n-1
        do j=1,i-1
            u(i,i)=u(i,i)-l(j,i)*u(i,j)
        end do
        u(i,i)=u(i,i)+a(i,i)
        do j=i+1,n
            do k=1,i-1
                u(j,i)=u(j,i)-l(k,i)*u(j,k)
                l(i,j)=l(i,j)-l(k,j)*u(i,k)
            end do
            u(j,i)=u(j,i)+a(j,i)
            l(i,j)=(l(i,j)+a(i,j))/u(i,i)
        end do
    end do
    do i=1,n-1
    u(n,n)=u(n,n)-l(i,n)*u(n,i)
    end do
    u(n,n)=u(n,n)+a(n,n)
    do i=1,n
        l(i,i)=1.0d0
    end do

    y_lu(1)=b_lu(1)
    do i=2,n
        do j=1,i-1
            y_lu(i)=y_lu(i)-y_lu(j)*l(j,i)
        end do
        y_lu(i)=y_lu(i)+b_lu(i)
    end do

    x_lu(n)=y_lu(n)/u(n,n)
    do i=n-1,1,-1
        do j=i+1,n
            x_lu(i)=x_lu(i)-x_lu(j)*u(j,i)
        end do
        x_lu(i)=(x_lu(i)+y_lu(i))/u(i,i)
    end do
    write(*,*)'test_lu:x_lu=',x_lu
end subroutine

! subroutine lu_test_from_src()
!     implicit none
!     real(8),allocatable :: a_lu(:,:),b_lu(:),l_lu(:,:),u_lu(:,:),y_lu(:),x_lu(:)
!     integer :: i_lu,j_lu, k_lu,n_lu
!     n_lu=4
!     allocate(x_lu(n_lu), y_lu(n_lu), b_lu(n_lu))
!     allocate(l_lu(n_lu, n_lu), u_lu(n_lu, n_lu),a_lu(n_lu, n_lu))
!     a_lu(1,1) = 9.0d0
!     a_lu(1,2) = 18.0d0
!     a_lu(1,3) = 9.0d0
!     a_lu(1,4) = -27.0d0

!     a_lu(2,1) = 18.0d0
!     a_lu(2,2) = 45.0d0
!     a_lu(2,3) = 0.0d0
!     a_lu(2,4) = -45.0d0

!     a_lu(3,1) = 9.0d0
!     a_lu(3,2) = 0.0d0
!     a_lu(3,3) = 126.0d0
!     a_lu(3,4) = 9.0d0

!     a_lu(4,1) = -27.0d0
!     a_lu(4,2) = -45.0d0
!     a_lu(4,3) = 9.0d0
!     a_lu(4,4) = 135.0d0

!     b_lu = (/1,2,16,8/)
!     call lu_solve()
!     write(*,*)'test_lu_from_src :x_lu=',x_lu
!    end subroutine



