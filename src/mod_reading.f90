#include <slepc/finclude/slepc.h>

module mod_reading
! ----------------------------------------------------------------
!
!  这个模块在0号进程读取文件。
!
!       call plot3d_load() 读取Plot3D格式的基本流和网格和来流扰动
!
! ----------------------------------------------------------------
    use petsc
    use mod_parameters
    implicit none
    public :: load
    private
    PetscScalar, pointer :: grid(:,:,:,:)
    PetscScalar, pointer :: flow(:,:,:,:)
    PetscErrorCode :: ierr
    PetscViewer :: Viewer
    Vec :: Multi_disturb
    Vec :: Flowfield
    DM :: singleDA 
    Vec :: Coord

contains

    subroutine load(comm)
        implicit none
        PetscInt, INTENT(in) :: comm

        call signal_loading(comm)
        call load_plot3d()

        call signal_convert()
        call set_seq_da(comm)
        call load_petsc_file(comm)
        call to_petsc_files(comm)

        call format_basic_file()
        call print_info()
        call deallocate_memory()
    end subroutine load

    subroutine signal_loading(comm)
        implicit none
        PetscInt,intent(in) :: comm
        call PetscPrintf(comm, "\n", ierr)
        call PetscPrintf(comm, " ===========================================================================\n", ierr)
        call PetscPrintf(comm, " =                                 读    取                                = \n", ierr)
        call PetscPrintf(comm, " ===========================================================================\n", ierr)
        call PetscPrintf(comm, " 「 主 进 程 」\n",ierr)
    end subroutine signal_loading

    ! 读取普通的数据文件
    subroutine load_plot3d()
        implicit none
        real(R_P),dimension(:,:,:,:),allocatable :: qq_0
        integer :: l,i,j,k 
        write(*,*) "-----------------------------------"
        write(*,*) "              读取数据              "
        write(*,*) "-----------------------------------"
        write(*,*) " "

        ! 读取网格信息
        write(*,*) "开始读取网格数据..."
        if(lns_mode==0) kn=1
        open(11, file=trim(gridfile),action='read',form='unformatted')
        read(11)
        select case (lns_mode)
        case(0)
            read(11) in,jn
            write(*,"(A,I5)") '   流向的网格数in=',in 
            write(*,"(A,I5)") '   法向的网格数jn=',jn 
            allocate(xx(in,jn,kn), yy(in,jn,kn), zz(in,jn,kn))
            read(11) xx,yy
            zz=0.0d0
        case(1)
            read(11) in,jn,kn
            write(*,"(A,I5)") '   流向的网格数in=',in
            write(*,"(A,I5)") '   法向的网格数jn=',jn
            write(*,"(A,I5)") '   展向的网格数kn=',kn
            allocate(xx(in,jn,kn), yy(in,jn,kn), zz(in,jn,kn))
            read(11) xx,yy,zz
        end select
        close(11)
        write(*,*) '  网格数据读取结束。'
        write(*,*) ""

        ! 读取基本流数据
        write(*,*) "开始读取流场数据..."
        open(12, file=trim(flowfile),action='read',form='unformatted')
        read(12)
        select case (lns_mode)
        case(0)
            read(12) in,jn,ln
            write(*,"(A,I5)") '   流向的网格数in=',in 
            write(*,"(A,I5)") '   法向的网格数jn=',jn 
            write(*,"(A,I5)") '   自由度ln=',ln 
        case(1)
            read(12) in,jn,kn,ln
            write(*,"(A,I5)") '   流向的网格数in=',in 
            write(*,"(A,I5)") '   法向的网格数jn=',jn 
            write(*,"(A,I5)") '   展向的网格数kn=',kn 
            write(*,"(A,I5)") '   自由度ln=',ln 
        end select
        allocate(qq_0(in,jn,kn,5))
        read(12)((((qq_0(i,j,k,l), i=1,in), j=1,jn), k=1,kn), l=1,5)
        close(12)
        ln=5
        allocate(qq(5,in,jn,kn))
        do k=1,kn
            do j=1,jn 
                do i=1,in 
                    qq(:,i,j,k)=qq_0(i,j,k,:)
                enddo
            enddo
        enddo
        deallocate(qq_0)
        write(*,*) '  流场信息读取结束。'
        write(*,*) ""

    end subroutine load_plot3d

    subroutine signal_convert()
        implicit none 
        write(*,*) "-----------------------------------"
        write(*,*) "         转换数据并生成文件          "
        write(*,*) "-----------------------------------"
        write(*,*)
    end subroutine signal_convert

    subroutine set_seq_da(comm)
        implicit none
        PetscInt, INTENT(in) :: comm

        call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
        &                 DMDA_STENCIL_BOX, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,&
        &                 3, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, coordDA, ierr)
        call DMSetUp(coordDA, ierr)

        call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
        &                 DMDA_STENCIL_STAR, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,&
        &                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, meshDA, ierr)
        call DMSetUp(meshDA, ierr)

        call DMDACreate2d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
        &                 DMDA_STENCIL_BOX, jn, kn, PETSC_DECIDE, PETSC_DECIDE, &
        &                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, singleDA, ierr)
        call DMSetUp(singleDA, ierr)

        call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
        &                 DMDA_STENCIL_BOX, size, jn, kn, PETSC_DECIDE, 1, 1,&
        &                 5, 0, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, disturbDA, ierr)
        call DMSetUp(disturbDA, ierr)

    end subroutine set_seq_da

    subroutine load_petsc_file(comm)
        implicit none
        PetscScalar, pointer :: multi(:,:,:,:)
        PetscScalar, pointer :: grid(:,:,:,:)
        PetscScalar, pointer :: flow(:,:,:,:)
        PetscScalar, pointer :: single(:,:,:)
        integer :: xs, ys, zs, xl, yl, zl
        PetscInt, INTENT(in) :: comm
        integer :: i, j, k, l
        Vec :: Single_disturb
        
        write(*,*) "开始转换为PetsC数据类型..."

        call DMGetGlobalVector(coordDA, Coord, ierr)
        call DMDAGetCorners(coordDA,xs,ys,zs,xl,yl,zl,ierr)
        call DMDAVecGetArrayF90(coordDA, coord, grid, ierr)
        do i=xs,xs+xl-1
            do j=ys,ys+yl-1
                do k=zs,zs+zl-1
                    grid(0, i, j, k) = xx(i+1, j+1, k+1)
                    grid(1, i, j, k) = yy(i+1, j+1, k+1)
                    grid(2, i, j, k) = zz(i+1, j+1, k+1)
                enddo
            enddo
        enddo
        call DMDAVecRestoreArrayF90(coordDA, coord, grid,ierr)
        write(*,*) '  网格信息转换结束。'

        call DMGetGlobalVector(MeshDA, Flowfield, ierr)
        call DMDAGetCorners(MeshDA,xs,ys,zs,xl,yl,zl,ierr)
        call DMDAVecGetArrayF90(MeshDA, Flowfield, flow, ierr)
        do i=xs,xs+xl-1
            do j=ys,ys+yl-1
                do k=zs,zs+zl-1
                    do l=0, 4
                        flow(l,i,j,k) = qq(l+1, i+1, j+1, k+1)
                    enddo
                enddo
            enddo
        enddo
        call DMDAVecRestoreArrayF90(MeshDA, Flowfield, flow,ierr)
        write(*,*) '  流场信息转换结束。'

        call DMGetGlobalVector(singleDA, Single_disturb, ierr)
        call PetscViewerBinaryOpen(comm, "in/LPSE_disturbance+1.petsc", FILE_MODE_READ, Viewer, ierr)
        call VecLoad(Single_disturb, Viewer, ierr)
        call PetscViewerDestroy(Viewer, ierr)

        call DMDAGetCorners(disturbDA,xs,ys,zs,xl,yl,zl,ierr)
        call DMGetGlobalVector(disturbDA, Multi_disturb, ierr)
        call DMDAVecGetArrayF90(disturbDA, Multi_disturb, multi, ierr)
        call DMDAVecGetArrayReadF90(singleDA, Single_disturb, single, ierr)
        do i=xs,xs+xl-1
            do j=ys,ys+yl-1
                do k=zs,zs+zl-1
                    multi(:,i,j,k)=single(:,j,k)
                enddo
            enddo
        enddo
        call DMDAVecRestoreArrayReadF90(singleDA, Single_disturb, single, ierr)
        call DMDAVecRestoreArrayF90(disturbDA, Multi_disturb, multi, ierr)
        write(*,*) '  来流信息转换结束。'

        write(*,*)
    end subroutine load_petsc_file

    subroutine to_petsc_files(comm)
        implicit none
        PetscInt, intent(in) :: comm

        write(*,*) "开始生成文件..."

        call PetscViewerBinaryOpen(comm, "in/grid.petsc",FILE_MODE_WRITE, Viewer, ierr)
        call VecView(Coord, Viewer, ierr)
        call PetscViewerDestroy(Viewer, ierr)
        write(*,*) '  网格文件已生成。'

        call PetscViewerBinaryOpen(comm, "in/flow.petsc",FILE_MODE_WRITE, Viewer, ierr)
        call VecView(Flowfield, Viewer, ierr)
        call PetscViewerDestroy(Viewer, ierr)
        write(*,*) '  流场文件已生成。'

        call PetscViewerBinaryOpen(comm, "in/disturb.petsc",FILE_MODE_WRITE, Viewer, ierr)
        call VecView(Multi_disturb, Viewer, ierr)
        call PetscViewerDestroy(Viewer, ierr)
        write(*,*) '  来流文件已生成。'

        write(*,*) ""
    end subroutine to_petsc_files

    subroutine print_info()
        implicit none
        write(*,*) "输出部分信息："
        write(*,*)
        select case (lns_mode)
        case(0)
            write(*,*) "  Ma =",Ma 
            write(*,*) "  Re =",Re 
            write(*,*) "  Te =",Te
            write(*,"(3X,A,2(F20.15))") "Alpha =",Alpha  
            write(*,"(3X,A,2(F20.15))") "Beta  =",Beta
            write(*,"(3X,A,2(F20.15))") "Omega =",Omega
        case(1)
            write(*,*) "  Ma =",Ma 
            write(*,*) "  Re =",Re 
            write(*,*) "  Te =",Te
            write(*,"(3X,A,2(F20.15))") "Omega =",Omega
        end select
        write(*,"(A,F10.5,' ->',F10.5)") "   流向起止位置: ",xx(1, 1, 1), xx(in, 1, 1)
        write(*,"(A,F10.5,' ->',F10.5)") "   法向起止位置: ",yy(1, 1, 1), yy(1, jn, 1)
        write(*,"(A,F10.5,' ->',F10.5)") "   展向起止位置: ",zz(1, 1, 1), zz(1, 1, kn)
        write(*,103) "  第一个数据是：",qq(1,1,1,1),qq(2,1,1,1),qq(3,1,1,1),qq(4,1,1,1),qq(5,1,1,1)
        103 format (1X,A,5(F10.5))
        write(*,103) "  第二个数据是：",qq(1,2,1,1),qq(2,2,1,1),qq(3,2,1,1),qq(4,2,1,1),qq(5,2,1,1)
        write(*,103) "  第三个数据是：",qq(1,3,1,1),qq(2,3,1,1),qq(3,3,1,1),qq(4,3,1,1),qq(5,3,1,1)
        write(*,104) "  第一个坐标是：",xx(1,1,1),yy(1,1,1),zz(1,1,1)
        104 format (1X,A,3(F10.5))
        write(*,104) "  第二个坐标是：",xx(2,1,1),yy(2,1,1),zz(2,1,1)
        write(*,104) "  第三个坐标是：",xx(3,1,1),yy(3,1,1),zz(3,1,1)
        write(*,*)
    end subroutine print_info

    subroutine format_basic_file()
        implicit none
        integer :: i,j,k
        open(30,file='out/hlns_info.csv',action='write',status='replace')
        write(30,*) "In,Jn,Kn,Alpha_r,Alpha_i,Beta_r,Beta_i,Omega_r,Omega_i"
        write(30,*) in,',',jn,',',kn,',',real(Alpha),',',aimag(Alpha),',',&
                    & real(Beta),',',aimag(Beta),',',real(Omega),',',aimag(Omega) 
        close(30)
        open(31, file='out/grid.csv',action='write',status='replace')
        write(31,*) "xx,yy,zz"
        do i=1,in 
            do j=1,jn 
                do k=1,kn
                    write(31,*) xx(i,j,k),',',yy(i,j,k),',',zz(i,j,k)
                enddo
            enddo
        enddo
        close(31)
        open(32,file='out/flow.csv',action='write',status='replace')
        write(32,*) 'rho,u,v,w,T'
        do i=1,in 
            do j=1,jn 
                do k=1,kn
                    write(32,*) qq(1,i,j,k),',',qq(2,i,j,k),',',qq(3,i,j,k),',',qq(4,i,j,k),',',qq(5,i,j,k)
                enddo
            enddo
        enddo
        close(32)
    end subroutine format_basic_file

    subroutine deallocate_memory()
        implicit none 
        deallocate(xx)
        deallocate(yy)
        deallocate(zz)
        deallocate(qq)
        call DMDestroy(disturbDA,ierr)
        call DMDestroy(coordDA,ierr)
        call DMDestroy(meshDA,ierr)
    end subroutine deallocate_memory
end module mod_reading
