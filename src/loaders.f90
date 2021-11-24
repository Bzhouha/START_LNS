module mod_loaders
! ----------------------------------------------------------------
!
!  这个模块在0号进程读取文件。
!
!       call plot3d_load() 读取Plot3D格式的基本流和网格和来流扰动
!
! ----------------------------------------------------------------
   use global_parameters
   implicit none
   public :: plot3d_load
   private
contains
   subroutine plot3d_load()
      select case(mode)
      case(0)
         call load2d()
      case(1)
         call load3d()
      end select
      call printinfo()
   end subroutine plot3d_load
   subroutine load2d()
      implicit none
      real(R_P),dimension(:,:,:),allocatable :: qq_0
      integer :: l,i,j,k 
      write(*,*) "-----------------------------------"
      write(*,*) "              读取数据              "
      write(*,*) "-----------------------------------"
      write(*,*) " "
      write(*,*) "开始读取网格数据..."
      kn=1
      open(11, file=trim(gridfile),action='read',form='unformatted')
      read(11)
      read(11) in, jn
      write(*,109) in,jn,kn
      109 format ('   流向的网格数in=',I5,/,'   法向的网格数jn=',I5,/,'   展向的网格数kn=',I5)
      allocate(xx(in, jn, kn), yy(in, jn, kn),zz(in, jn, kn))
      read(11) xx,yy
      close(11)
      write(*,*) '  网格数据读取结束。'
      write(*,*) ""
      ! 读取基本流数据
      write(*,*) "开始读取流场数据..."
      open(12, file=trim(flowfile),action='read',form='unformatted')
      read(12)
      read(12) in,jn,ln
      write(*,110) in,jn,kn,ln
      110 format ('   流向的网格数in=',I5,/,'   法向的网格数jn=',I5,/,&
      &'   展向的网格数kn=',I5,/,"   自由度ln=",I3)
      allocate(qq_0(in, jn, 5))
      read(12)((((qq_0(i, j, l), i=1, in), j=1, jn), k=1, kn), l=1, 5)
      close(12)
      ln=5;zz=0.0d0
      allocate(qq(5,in,jn,kn))
      do j=1,jn 
         do i=1,in 
            qq(:,i,j,1)=qq_0(i,j,:)
         enddo
      enddo
      deallocate(qq_0)
      write(*,*) '  流场信息读取结束。'
      write(*,*) ""
      ! 读取来流扰动
      write(*,*) "开始读取来流扰动..."
      open(13, file="..//files//in//Disturbance.dat",action='read')
      read(13,*)
      allocate(inflow(ln,0:jn-1,kn))
      do j=0,jn-1 
            read(13,*) inflow(:,j,1)
      enddo
      close(13)
      write(*,*) "  来流扰动读取结束。"
   end subroutine load2d
   subroutine load3d()
      implicit none
      integer :: i,j,k 
      write(*,*) "-----------------------------------"
      write(*,*) "              读取数据              "
      write(*,*) "-----------------------------------"
      write(*,*) " "
      write(*,*) "开始读取网格数据..."
      open(11, file=trim(gridfile),action='read',form='unformatted')
      read(11)
      read(11) in, jn, kn
      write(*,109) in,jn,kn
      109 format ('   流向的网格数in=',I5,/,'   法向的网格数jn=',I5,/,'   展向的网格数kn=',I5)
      allocate(xx(in, jn, kn), yy(in, jn, kn), zz(in, jn, kn))
      do k=1,kn
         do j=1,jn 
             do i=1,in
                 read(11) xx(i,j,k),yy(i,j,k),zz(i,j,k)
             enddo
         enddo
      enddo
      close(11)
      write(*,*) '  网格数据读取结束。'
      write(*,*) ""
      ! 读取基本流数据
      write(*,*) "开始读取流场数据..."
      open(12, file=trim(flowfile),action='read',form='unformatted')
      read(12)
      read(12) in, jn, kn, ln
      write(*,110) in,jn,kn,ln
      110 format ('   流向的网格数in=',I5,/,'   法向的网格数jn=',I5,/,&
      &'   展向的网格数kn=',I5,/,"   自由度ln=",I3)
      allocate(qq(ln, in, jn, kn))
      do k=1,kn
         do j=1,jn 
             do i=1,in 
                 read(12) qq(:,i,j,k)
             enddo
         enddo
      enddo
      close(12)
      write(*,*) '  流场信息读取结束。'
      write(*,*) ""
      ! 读取来流扰动
      write(*,*) "开始读取来流扰动..."
      open(13, file="..//files//in//Disturbance.dat",action='read')
      read(13,*)
      allocate(inflow(ln,0:jn-1,0:kn-1))
      do k=0,kn-1
         do j=0,jn-1 
               read(13,*) inflow(:,j,k)
         enddo
      enddo
      close(13)
      write(*,*) "  来流扰动读取结束。"
   end subroutine load3d

   subroutine printinfo()
      implicit none
      write(*,*) ""
      write(*,*) "输出部分信息："
      write(*,"(A,F10.5,F10.5)") "  流向起止位置: ",xx(1, 1, 1), xx(in, jn, kn)
      write(*,103) "  第一个数据是：",qq(1,1,1,1),qq(2,1,1,1),qq(3,1,1,1),qq(4,1,1,1),qq(5,1,1,1)
      103 format (A,5(F10.5))
      write(*,103) "  第二个数据是：",qq(1,2,1,1),qq(2,2,1,1),qq(3,2,1,1),qq(4,2,1,1),qq(5,2,1,1)
      write(*,103) "  第三个数据是：",qq(1,3,1,1),qq(2,3,1,1),qq(3,3,1,1),qq(4,3,1,1),qq(5,3,1,1)
      write(*,104) "  第一个坐标是：",xx(1,1,1),yy(1,1,1),zz(1,1,1)
      104 format (A,3(F10.5))
      write(*,104) "  第二个坐标是：",xx(2,1,1),yy(2,1,1),zz(2,1,1)
      write(*,104) "  第三个坐标是：",xx(3,1,1),yy(3,1,1),zz(3,1,1)
      write(*,*)
   end subroutine printinfo
end module mod_loaders
