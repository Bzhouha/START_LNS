module mod_mftools
    use mod_parameters, only: in,jn,kn
    use penf, only: R_P
    implicit none
    public
    contains
    subroutine f5d1(out_array,f)
        implicit none
        complex(R_P),dimension(5),intent(out) :: out_array
        complex(R_P),dimension(5,5),intent(in) :: f
            out_array(:)=(1.0d0*f(:,1)-8.0d0*f(:,2)+8.0d0*f(:,4)-1.0d0*f(:,5))/12.0d0
    end subroutine f5d1

    subroutine f5d11(out_array,f)
        implicit none
        complex(R_P),dimension(5),intent(out) :: out_array
        complex(R_P),dimension(5,5),intent(in) :: f
            out_array(:)=(-1.0d0*f(:,1)+16.0d0*f(:,2)-30.0d0*f(:,3)+16.0d0*f(:,4)-1.0d0*f(:,5))/12.0d0 
    end subroutine f5d11

    subroutine f5d12(out_array,f)
        implicit none
        complex(R_P),dimension(5),intent(out) :: out_array
        complex(R_P),dimension(5,5,5),intent(in) :: f
        complex(R_P),dimension(5,5) :: tmp
            call f5d1(tmp(:,1),f(:,1,:))
            call f5d1(tmp(:,2),f(:,2,:))
            call f5d1(tmp(:,4),f(:,4,:))
            call f5d1(tmp(:,5),f(:,5,:))
            call f5d1(out_array,tmp)
    end subroutine f5d12

    subroutine f4d1(out_array,f,i,j,k,flg)
        implicit none
        complex(R_P),dimension(5),intent(out) :: out_array
        complex(R_P),dimension(5,4),intent(in) :: f
        integer,intent(in) :: i,j,k
        integer,intent(in) :: flg 
        integer :: l,ln 
        select case (flg) 
        case(1)
            l=i;ln=in 
        case(2)
            l=j;ln=jn
        case(3)
            l=k;ln=kn
        end select
        if(l>=1 .and. l<=(ln-3))then
            out_array(:)=((-2.0d0)*f(:,1)-3.0d0*f(:,2)+6.0d0*f(:,3)-1.0d0*f(:,4))/6.0d0
        else if(l==(ln-2))then
            out_array(:)=(1.0d0*f(:,1)-6.0d0*f(:,2)+3.0d0*f(:,3)+2.0d0*f(:,4))/6.0d0
        endif
    end subroutine f4d1

    subroutine f4d11(out_array,f)
        implicit none
        complex(R_P),dimension(5),intent(out) :: out_array
        complex(R_P),dimension(5,3),intent(in) :: f
        out_array(:)=1.0d0*f(:,1)-2.0d0*f(:,2)+1.0d0*f(:,3)
    end subroutine f4d11

    subroutine f4d12(out_array,f,i,j,k,flg)
        implicit none
        complex(R_P),dimension(5),intent(out) :: out_array
        complex(R_P),dimension(5,4,4),intent(in) :: f
        complex(R_P),dimension(5,4) :: tmp 
        integer,intent(in) :: i,j,k 
        integer,intent(in) :: flg
        integer :: flg1,flg2
        select case (flg)
        case(12)
            flg1=1;flg2=2
        case(13)
            flg1=1;flg2=3
        case(23)
            flg1=2;flg2=3
        end select
        call f4d1(tmp(:,1),f(:,1,:),i,j,k,flg2)
        call f4d1(tmp(:,2),f(:,2,:),i,j,k,flg2)
        call f4d1(tmp(:,3),f(:,3,:),i,j,k,flg2)
        call f4d1(tmp(:,4),f(:,4,:),i,j,k,flg2)
        call f4d1(out_array(:),tmp(:,:),i,j,k,flg1)
    end subroutine f4d12 

    subroutine f45d12(out_array,f,i,j,k,flg)
        implicit none
        complex(R_P),dimension(5),intent(out) :: out_array
        complex(R_P),dimension(5,4,5),intent(in) :: f
        complex(R_P),dimension(5,4) :: tmp 
        integer,intent(in) :: i,j,k 
        integer,intent(in) :: flg
        call f5d1(tmp(:,1),f(:,1,:))
        call f5d1(tmp(:,2),f(:,2,:))
        call f5d1(tmp(:,3),f(:,3,:))
        call f5d1(tmp(:,4),f(:,4,:))
        call f4d1(out_array(:),tmp(:,:),i,j,k,flg)
    end subroutine f45d12 

    subroutine f4_index(i1,i2,j1,j2,i,j)
        implicit none
        integer,intent(out) :: i1,i2,j1,j2
        integer,intent(in) :: i,j
        if(i>=1 .and. i<=(in-3))then
            i1=i-1;i2=i+2
        else if(i==(in-2))then
            i1=i-2;i2=i+1
        endif
        if(j>=1 .and. j<=(jn-3))then
            j1=j-1;j2=j+2
        else if(j==(jn-2))then
            j1=j-2;j2=j+1
        endif
    end subroutine f4_index

    subroutine f3d1f(out_array,f)
        implicit none
        complex(R_P), dimension(5), intent(out) :: out_array
        complex(R_P),dimension(5,3),intent(in) :: f
        out_array(:)=(-3.0d0*f(:,1)+4.0d0*f(:,2)-1.0d0*f(:,3))/2.0d0 
    end subroutine f3d1f

    subroutine f3d1r(out_array,f)
        implicit none
        complex(R_P), dimension(5), intent(out) :: out_array
        complex(R_P),dimension(5,3),intent(in) :: f
        out_array(:)=(1.0d0*f(:,1)-4.0d0*f(:,2)+3.0d0*f(:,3))/2.0d0 
    end subroutine f3d1r

    subroutine f2d1f(out_array,f)
        implicit none
        complex(R_P), dimension(5), intent(out) :: out_array
        complex(R_P),dimension(5,2),intent(in) :: f
        out_array(:)=(-1.0d0*f(:,1)+1.0d0*f(:,2))
    end subroutine f2d1f

    subroutine f2d1r(out_array,f)
        implicit none
        complex(R_P), dimension(5), intent(out) :: out_array
        complex(R_P),dimension(5,2),intent(in) :: f
        out_array(:)=(-1.0d0*f(:,1)+1.0d0*f(:,2))
    end subroutine f2d1r

end module mod_mftools