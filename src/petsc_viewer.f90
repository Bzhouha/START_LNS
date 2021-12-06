!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_petsc_viewer
! ----------------------------------------------------------------
!
!  这个模块在0号进程生成PetsC格式文件。
!
!       call petsc_viewer(comm) 
!
! ----------------------------------------------------------------
  use petsc
  use global_parameters
  implicit none
  private
  PetscScalar, pointer :: grid(:,:,:,:)
  PetscScalar, pointer :: flow(:,:,:,:)
  public :: petsc_viewer
  PetscErrorCode :: ierr
  Vec :: Coord, Flowfield
  PetscViewer :: Viewer
contains
  subroutine petsc_viewer(comm)
    implicit none
    PetscInt, INTENT(in) :: comm
    call signal_viewer()
    call create_mesh_3d(comm)
    call mesh_output(comm)
    call deallocate_memory()
  end subroutine petsc_viewer

  subroutine signal_viewer()
    implicit none 
    write(*,*) "-----------------------------------"
    write(*,*) "         转换数据并生成文件          "
    write(*,*) "-----------------------------------"
    write(*,*)
  end subroutine signal_viewer

  subroutine create_mesh_3d(comm)
    implicit none
    PetscInt, INTENT(in) :: comm
    integer :: xs, ys, zs, xl, yl, zl
    integer :: i, j, k, l
    write(*,*) "开始转换为PetsC数据类型..."
    ! 转换网格
    call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
    &                 DMDA_STENCIL_BOX, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,&
    &                 3, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, coordDA, ierr)
    call DMSetUp(coordDA, ierr)
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
    ! 转换流场
    call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
    &                 DMDA_STENCIL_STAR, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,&
    &                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, meshDA, ierr)
    call DMSetUp(meshDA, ierr)
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
    write(*,*)
  end subroutine create_mesh_3d

  subroutine mesh_output(comm)
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
  end subroutine mesh_output

  subroutine deallocate_memory()
    implicit none 
    deallocate(xx)
    deallocate(yy)
    deallocate(zz)
    deallocate(qq)
    call DMDestroy(coordDA,ierr)
    call DMDestroy(meshDA,ierr)
  end subroutine deallocate_memory
end module mod_petsc_viewer
