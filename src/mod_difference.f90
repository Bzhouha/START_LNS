module mod_difference
! ----------------------------------------------------
! 
!  这个模块是计算导数的函数
! 
!       1.call fd1(out_array,mo,no,po,qo,ro,so,in_array,mi,ni,pi,qi,ri,si,flg,dof) 计算一整块数据的一阶导数
! 
!       2.call fd2(out_array,mo,no,po,qo,ro,so,in_array,mi,ni,pi,qi,ri,si,flg,dof) 计算一整块数据的二阶导数
! 
! ----------------------------------------------------
    use penf, only: R_P
    use mod_parameters
    implicit none
    public :: fd1,fd2
    contains

    subroutine fd1(out_array,mo,no,po,qo,ro,so,&
                    in_array,mi,ni,pi,qi,ri,si,flg,dof)
        implicit none
        integer :: i,j,k
        integer,intent(in) :: flg,mi,ni,pi,qi,ri,si,mo,no,po,qo,ro,so,dof
        real(R_P), dimension(dof,mi:ni,pi:qi,ri:si), intent(in) :: in_array
        real(R_P), dimension(dof,mo:no,po:qo,ro:so), intent(out) :: out_array
        select case (flg)
        case(1)
            do k=ro,so 
                do j=po,qo 
                    do i=mo,no 
                        if(i > 1 .and. i < (in-2))then
                            out_array(:,i,j,k) = (in_array(:,i-2,j,k)-8.0d0*in_array(:,i-1,j,k) &
                                                &+8.0d0*in_array(:,i+1,j,k)-in_array(:,i+2,j,k))/12.0d0
                        else if (i == 0) then
                            out_array(:,i,j,k) = ((-3.0d0)*in_array(:,i,j,k)+4.0d0*in_array(:,i+1,j,k)-in_array(:,i+2,j,k))/2.0d0
                        else if (i == 1) then
                            out_array(:,i,j,k) = ((-2.0d0)*in_array(:,i-1,j,k)-3.0d0*in_array(:,i,j,k) &
                                                &+6.0d0*in_array(:,i+1,j,k)-in_array(:,i+2,j,k))/6.0d0
                        else if (i == (in-2)) then
                            out_array(:,i,j,k) = (in_array(:,i-2,j,k)-6.0d0*in_array(:,i-1,j,k) &
                                                &+3.0d0*in_array(:,i,j,k)+2.0d0*in_array(:,i+1,j,k))/6.0d0
                        else if (i == (in-1)) then
                            out_array(:,i,j,k) = (in_array(:,i-2,j,k)-4.0d0*in_array(:,i-1,j,k)+3.0d0*in_array(:,i,j,k))/2.0d0
                        end if
                    enddo
                enddo
            enddo
        case(2)
            do k=ro,so 
                do j=po,qo 
                    do i=mo,no 
                        if(j > 1 .and. j < (jn-2))then
                            out_array(:,i,j,k) = (in_array(:,i,j-2,k)-8.0d0*in_array(:,i,j-1,k) &
                                                &+8.0d0*in_array(:,i,j+1,k)-in_array(:,i,j+2,k))/12.0d0
                        else if (j == 0) then
                            out_array(:,i,j,k) = ((-3.0d0)*in_array(:,i,j,k)+4.0d0*in_array(:,i,j+1,k)-in_array(:,i,j+2,k))/2.0d0
                        else if (j == 1) then
                            out_array(:,i,j,k) = ((-2.0d0)*in_array(:,i,j-1,k)-3.0d0*in_array(:,i,j,k) &
                                                &+6.0d0*in_array(:,i,j+1,k)-in_array(:,i,j+2,k))/6.0d0
                        else if (j == (jn-2)) then
                            out_array(:,i,j,k) = (in_array(:,i,j-2,k)-6.0d0*in_array(:,i,j-1,k) &
                                                &+3.0d0*in_array(:,i,j,k)+2.0d0*in_array(:,i,j+1,k))/6.0d0
                        else if (j == (jn-1)) then
                            out_array(:,i,j,k) = (in_array(:,i,j-2,k)-4.0d0*in_array(:,i,j-1,k)+3.0d0*in_array(:,i,j,k))/2.0d0
                        end if
                    enddo
                enddo
            enddo
        case(3)
            do k=ro,so 
                do j=po,qo 
                    do i=mo,no 
                        if(k > 1 .and. k < (kn-2))then
                            out_array(:,i,j,k) = (in_array(:,i,j,k-2)-8.0d0*in_array(:,i,j,k-1) &
                                                &+8.0d0*in_array(:,i,j,k+1)-in_array(:,i,j,k+2))/12.0d0
                        else if (k == 0) then
                            out_array(:,i,j,k) = ((-3.0d0)*in_array(:,i,j,k)+4.0d0*in_array(:,i,j,k+1)-in_array(:,i,j,k+2))/2.0d0
                        else if (k == 1) then
                            out_array(:,i,j,k) = ((-2.0d0)*in_array(:,i,j,k-1)-3.0d0*in_array(:,i,j,k) &
                                                &+6.0d0*in_array(:,i,j,k+1)-in_array(:,i,j,k+2))/6.0d0
                        else if (k == (kn-2)) then
                            out_array(:,i,j,k) = (in_array(:,i,j,k-2)-6.0d0*in_array(:,i,j,k-1) &
                                                &+3.0d0*in_array(:,i,j,k)+2.0d0*in_array(:,i,j,k+1))/6.0d0
                        else if (k == (kn-1)) then
                            out_array(:,i,j,k) = (in_array(:,i,j,k-2)-4.0d0*in_array(:,i,j,k-1)+3.0d0*in_array(:,i,j,k))/2.0d0
                        end if
                    enddo
                enddo
            enddo
        end select
    end subroutine fd1

    subroutine fd2(out_array,mo,no,po,qo,ro,so,&
                    in_array,mi,ni,pi,qi,ri,si,flg,dof)
        implicit none
        integer :: i,j,k
        integer,intent(in) :: flg,mi,ni,pi,qi,ri,si,mo,no,po,qo,ro,so,dof
        real(R_P), dimension(dof,mi:ni,pi:qi,ri:si), intent(in) :: in_array
        real(R_P), dimension(dof,mo:no,po:qo,ro:so), intent(out) :: out_array
        select case (flg)
        case(1)
            do k=ro,so 
                do j=po,qo 
                    do i=mo,no 
                        if(i > 1 .and. i < (in-2))then
                            out_array(:,i,j,k) = ((-1.0d0)*in_array(:,i-2,j,k)+16.0d0*in_array(:,i-1,j,k) &
                                                &-30.0d0*in_array(:,i,j,k)+16.0d0*in_array(:,i+1,j,k)-in_array(:,i+2,j,k))/12.0d0
                        else if (i == 0) then
                            out_array(:,i,j,k) = (in_array(:,i,j,k)-2.0d0*in_array(:,i+1,j,k)+in_array(:,i+2,j,k))
                        else if (i == 1 .or. i == (in-2)) then
                            out_array(:,i,j,k) = (in_array(:,i-1,j,k)-2.0d0*in_array(:,i,j,k)+in_array(:,i+1,j,k))
                        else if (i == (in-1)) then
                            out_array(:,i,j,k) = (in_array(:,i-2,j,k)-2.0d0*in_array(:,i-1,j,k)+in_array(:,i,j,k))
                        end if
                    enddo
                enddo
            enddo
        case(2)
            do k=ro,so 
                do j=po,qo 
                    do i=mo,no 
                        if(j > 1 .and. j < (jn-2))then
                            out_array(:,i,j,k) = ((-1.0d0)*in_array(:,i,j-2,k)+16.0d0*in_array(:,i,j-1,k) &
                                                &-30.0d0*in_array(:,i,j,k)+16.0d0*in_array(:,i,j+1,k)-in_array(:,i,j+2,k))/12.0d0
                        else if (j == 0) then
                            out_array(:,i,j,k) = (in_array(:,i,j,k)-2.0d0*in_array(:,i,j+1,k)+in_array(:,i,j+2,k))
                        else if (j == 1 .or. j == (jn-2)) then
                            out_array(:,i,j,k) = (in_array(:,i,j-1,k)-2.0d0*in_array(:,i,j,k)+in_array(:,i,j+1,k))
                        else if (j == (jn-1)) then
                            out_array(:,i,j,k) = (in_array(:,i,j-2,k)-2.0d0*in_array(:,i,j-1,k)+in_array(:,i,j,k))
                        end if
                    enddo
                enddo
            enddo
        case(3)
            do k=ro,so 
                do j=po,qo 
                    do i=mo,no 
                        if(k > 1 .and. k < (kn-2))then
                            out_array(:,i,j,k) = ((-1.0d0)*in_array(:,i,j,k-2)+16.0d0*in_array(:,i,j,k-1) &
                                                &-30.0d0*in_array(:,i,j,k)+16.0d0*in_array(:,i,j,k+1)-in_array(:,i,j,k+2))/12.0d0
                        else if (k == 0) then
                            out_array(:,i,j,k) = (in_array(:,i,j,k)-2.0d0*in_array(:,i,j,k+1)+in_array(:,i,j,k+2))
                        else if (k == 1 .or. k == (kn-2)) then
                            out_array(:,i,j,k) = (in_array(:,i,j,k-1)-2.0d0*in_array(:,i,j,k)+in_array(:,i,j,k+1))
                        else if (k == (kn-1)) then
                            out_array(:,i,j,k) = (in_array(:,i,j,k-2)-2.0d0*in_array(:,i,j,k-1)+in_array(:,i,j,k))
                        end if
                    enddo
                enddo
            enddo
        end select
    end subroutine fd2
end module mod_difference

