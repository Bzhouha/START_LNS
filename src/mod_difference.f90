module mod_difference
! ----------------------------------------------------
! 
!  这个模块是计算导数的函数
! 
!       1.call Cal_dif(result, array, Delta, istart, iend, ighost_start, ighost_end, ln) 计算某一行数据的一阶导数，不太使用了
! 
!       2.call fd1(out_array,mo,no,po,qo,ro,so,in_array,mi,ni,pi,qi,ri,si,flg,tag) 计算一整块数据的一阶导数
! 
!       3.call fd2(out_array,mo,no,po,qo,ro,so,in_array,mi,ni,pi,qi,ri,si,flg,tag) 计算一整块数据的二阶导数
! 
! ----------------------------------------------------
    use penf, only: R_P
    use global_parameters
    implicit none
    public :: fd1,fd2
    contains
    subroutine Cal_dif(result, array, Delta, istart, iend, ighost_start, ighost_end, ln)
        implicit None
        real(R_P), intent(out),dimension(istart:iend) :: result
        real(R_P), intent(in), dimension(ighost_start:ighost_end) :: array 
        real(R_P), intent(in) :: Delta 
        integer, intent(in) :: istart, iend, ighost_start, ighost_end, ln
        integer :: l
        
        do l=istart,iend
            if ( l > 1 .and. l < (ln-2) ) then
                result(l) = (array(l-2) - 8.0d0*array(l-1) + 8.0d0*array(l+1) - array(l+2))/(12.0d0*Delta)
            else if (l == 0) then
                result(l) = ((-3.0d0)*array(l) + 4.0d0*array(l+1) - array(l+2))/(2.0d0*Delta)
            else if (l == 1) then
                result(l) = ((-2.0d0)*array(l-1) - 3.0d0*array(l) + 6.0d0*array(l+1) - array(l+2))/(6.0d0*Delta)
            else if (l == (ln-2)) then
                result(l) = (array(l-2) - 6.0d0*array(l-1) + 3.0d0*array(l) + 2.0d0*array(l+1))/(6.0d0*Delta)
            else if (l == (ln-1)) then
                result(l) = (array(l-2) - 4.0d0*array(l-1) + 3.0d0*array(l))/(2.0d0*Delta)
            end if
        enddo
    end subroutine Cal_dif

    subroutine fd1(out_array,mo,no,po,qo,ro,so,&
                    in_array,mi,ni,pi,qi,ri,si,&
                    flg,tag)
        implicit none
        integer :: i,j,k
        integer,intent(in) :: flg,mi,ni,pi,qi,ri,si,mo,no,po,qo,ro,so,tag
        real(R_P), dimension(tag,mi:ni,pi:qi,ri:si), intent(in) :: in_array
        real(R_P), dimension(tag,mo:no,po:qo,ro:so), intent(out) :: out_array
        select case (flg)
        case(1)
            do k=ro,so 
                do j=po,qo 
                    do i=mo,no 
                        if(i > 1 .and. i < (in-2))then
                            out_array(:,i,j,k)=(in_array(:,i-2,j,k)-8.0d0*in_array(:,i-1,j,k) &
                            +8.0d0*in_array(:,i+1,j,k)-in_array(:,i+2,j,k))/12.0d0
                        else if (i == 0) then
                            out_array(:,i,j,k) = ((-3.0d0)*in_array(:,i,j,k) + 4.0d0*in_array(:,i+1,j,k) - in_array(:,i+2,j,k))/2.0d0
                        else if (i == 1) then
                            out_array(:,i,j,k) = ((-2.0d0)*in_array(:,i-1,j,k) - 3.0d0*in_array(:,i,j,k) &
                             + 6.0d0*in_array(:,i+1,j,k) - in_array(:,i+2,j,k))/6.0d0
                        else if (i == (in-2)) then
                            out_array(:,i,j,k) = (in_array(:,i-2,j,k) - 6.0d0*in_array(:,i-1,j,k) &
                            + 3.0d0*in_array(:,i,j,k) + 2.0d0*in_array(:,i+1,j,k))/6.0d0
                        else if (i == (in-1)) then
                            out_array(:,i,j,k) = (in_array(:,i-2,j,k) - 4.0d0*in_array(:,i-1,j,k) + 3.0d0*in_array(:,i,j,k))/2.0d0
                        end if
                    enddo
                enddo
            enddo
        case(2)
            do k=ro,so 
                do j=po,qo 
                    do i=mo,no 
                        if(j > 1 .and. j < (jn-2))then
                            out_array(:,i,j,k)=( in_array(:,i,j-2,k)-8.0d0*in_array(:,i,j-1,k)&
                            +8.0d0*in_array(:,i,j+1,k)-in_array(:,i,j+2,k) )/12.0d0
                        else if (j == 0) then
                            out_array(:,i,j,k) = ((-3.0d0)*in_array(:,i,j,k) + 4.0d0*in_array(:,i,j+1,k) - in_array(:,i,j+2,k))/2.0d0
                        else if (j == 1) then
                            out_array(:,i,j,k) = ((-2.0d0)*in_array(:,i,j-1,k) - 3.0d0*in_array(:,i,j,k) &
                            + 6.0d0*in_array(:,i,j+1,k) - in_array(:,i,j+2,k))/6.0d0
                        else if (j == (jn-2)) then
                            out_array(:,i,j,k) = (in_array(:,i,j-2,k) - 6.0d0*in_array(:,i,j-1,k) &
                             + 3.0d0*in_array(:,i,j,k) + 2.0d0*in_array(:,i,j+1,k))/6.0d0
                        else if (j == (jn-1)) then
                            out_array(:,i,j,k) = (in_array(:,i,j-2,k) - 4.0d0*in_array(:,i,j-1,k) + 3.0d0*in_array(:,i,j,k))/2.0d0
                        end if
                    enddo
                enddo
            enddo
        case(3)
            do k=ro,so 
                do j=po,qo 
                    do i=mo,no 
                        if(k > 1 .and. k < (kn-2))then
                            out_array(:,i,j,k)=( in_array(:,i,j,k-2)-8.0d0*in_array(:,i,j,k-1) &
                            +8.0d0*in_array(:,i,j,k+1)-in_array(:,i,j,k+2) )/12.0d0
                        else if (k == 0) then
                            out_array(:,i,j,k) = ((-3.0d0)*in_array(:,i,j,k) + 4.0d0*in_array(:,i,j,k+1) - in_array(:,i,j,k+2))/2.0d0
                        else if (k == 1) then
                            out_array(:,i,j,k) = ((-2.0d0)*in_array(:,i,j,k-1) - 3.0d0*in_array(:,i,j,k) &
                             + 6.0d0*in_array(:,i,j,k+1) - in_array(:,i,j,k+2))/6.0d0
                        else if (k == (kn-2)) then
                            out_array(:,i,j,k) = (in_array(:,i,j,k-2) - 6.0d0*in_array(:,i,j,k-1) &
                             + 3.0d0*in_array(:,i,j,k) + 2.0d0*in_array(:,i,j,k+1))/6.0d0
                        else if (k == (kn-1)) then
                            out_array(:,i,j,k) = (in_array(:,i,j,k-2) - 4.0d0*in_array(:,i,j,k-1) + 3.0d0*in_array(:,i,j,k))/2.0d0
                        end if
                    enddo
                enddo
            enddo
        end select
    end subroutine fd1

    subroutine fd2(out_array,mo,no,po,qo,ro,so,&
                    in_array,mi,ni,pi,qi,ri,si,&
                    flg,tag)
        implicit none
        integer :: i,j,k
        integer,intent(in) :: flg,mi,ni,pi,qi,ri,si,mo,no,po,qo,ro,so,tag
        real(R_P), dimension(tag,mi:ni,pi:qi,ri:si), intent(in) :: in_array
        real(R_P), dimension(tag,mo:no,po:qo,ro:so), intent(out) :: out_array
        select case (flg)
        case(1)
            do k=ro,so 
                do j=po,qo 
                    do i=mo,no 
                        if(i > 1 .and. i < (in-2))then
                            out_array(:,i,j,k)=( (-1.0d0)*in_array(:,i-2,j,k)+16.0d0*in_array(:,i-1,j,k) &
                            -30.0d0*in_array(:,i,j,k)+16.0d0*in_array(:,i+1,j,k)-in_array(:,i+2,j,k) )/12.0d0
                        else if (i == 0) then
                            out_array(:,i,j,k) = (in_array(:,i,j,k) - 2.0d0*in_array(:,i+1,j,k) + in_array(:,i+2,j,k))
                        else if (i == 1 .or. i == (in-2)) then
                            out_array(:,i,j,k) = (in_array(:,i-1,j,k) - 2.0d0*in_array(:,i,j,k) + in_array(:,i+1,j,k))
                        else if (i == (in-1)) then
                            out_array(:,i,j,k) = (in_array(:,i-2,j,k) - 2.0d0*in_array(:,i-1,j,k) + in_array(:,i,j,k))
                        end if
                    enddo
                enddo
            enddo
        case(2)
            do k=ro,so 
                do j=po,qo 
                    do i=mo,no 
                       if(j > 1 .and. j < (jn-2))then
                           out_array(:,i,j,k)=( (-1.0d0)*in_array(:,i,j-2,k)+16.0d0*in_array(:,i,j-1,k) &
                           -30.0d0*in_array(:,i,j,k)+16.0d0*in_array(:,i,j+1,k)-in_array(:,i,j+2,k) )/12.0d0
                       else if (j == 0) then
                           out_array(:,i,j,k) = (in_array(:,i,j,k) - 2.0d0*in_array(:,i,j+1,k) + in_array(:,i,j+2,k))
                       else if (j == 1 .or. j == (jn-2)) then
                           out_array(:,i,j,k) = (in_array(:,i,j-1,k) - 2.0d0*in_array(:,i,j,k) + in_array(:,i,j+1,k))
                       else if (j == (jn-1)) then
                           out_array(:,i,j,k) = (in_array(:,i,j-2,k) - 2.0d0*in_array(:,i,j-1,k) + in_array(:,i,j,k))
                       end if
                    enddo
                enddo
            enddo
        case(3)
            do k=ro,so 
                do j=po,qo 
                    do i=mo,no 
                       if(k > 1 .and. k < (kn-2))then
                           out_array(:,i,j,k)=( (-1.0d0)*in_array(:,i,j,k-2)+16.0d0*in_array(:,i,j,k-1) &
                           -30.0d0*in_array(:,i,j,k)+16.0d0*in_array(:,i,j,k+1)-in_array(:,i,j,k+2) )/12.0d0
                       else if (k == 0) then
                           out_array(:,i,j,k) = (in_array(:,i,j,k) - 2.0d0*in_array(:,i,j,k+1) + in_array(:,i,j,k+2))
                       else if (k == 1 .or. k == (kn-2)) then
                           out_array(:,i,j,k) = (in_array(:,i,j,k-1) - 2.0d0*in_array(:,i,j,k) + in_array(:,i,j,k+1))
                       else if (k == (kn-1)) then
                           out_array(:,i,j,k) = (in_array(:,i,j,k-2) - 2.0d0*in_array(:,i,j,k-1) + in_array(:,i,j,k))
                       end if
                    enddo
                enddo
            enddo
        end select
    end subroutine fd2

end module mod_difference

