!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_petsc_loader ! 读入并分发数据
  use petsc
  use global_parameters
  implicit none
  private
  PetscScalar, pointer :: grid(:,:,:,:)
  PetscScalar, pointer :: flow(:,:,:,:)
  Vec :: Coord_local, Flowfield_local
  PetscErrorCode  :: ierr
  Vec :: Coord, Flowfield
  PetscViewer :: Viewer
  public :: load_info
  contains
	subroutine load_info(comm)
		implicit none
		PetscInt, intent(in) :: comm
		call read_mesh_3d(comm)
		call get_layout()
		call load_mesh_info()
		call load_flow_info()
		call OhWeHaveThese(comm)
	end subroutine load_info

	subroutine OhWeHaveThese(comm)
		implicit none 
		PetscInt,intent(in) :: comm 
		PetscErrorCode :: ierr  
		call PetscPrintf(comm," -----------------------------------\n",ierr)
		call PetscPrintf(comm,"          读入数据信息...      \n",ierr)
		call PetscPrintf(comm," -----------------------------------\n",ierr)
		if(rank==0)then
			write(*,*)
			write(*,*) "输出部分信息："
			write(*,*) 
			write(*,9) in,jn,kn
			9 format ('   流向的网格数in=',I5,/,'   法向的网格数jn=',I5,/,'   展向的网格数kn=',I5)
			write(*,113) "   第一个数据是：",qq(1,0,0,0),qq(2,0,0,0),qq(3,0,0,0),qq(4,0,0,0),qq(5,0,0,0)
			113 format (A,5(F10.5))
			write(*,113) "   第二个数据是：",qq(1,1,0,0),qq(2,1,0,0),qq(3,1,0,0),qq(4,1,0,0),qq(5,1,0,0)
			write(*,113) "   第三个数据是：",qq(1,2,0,0),qq(2,2,0,0),qq(3,2,0,0),qq(4,2,0,0),qq(5,2,0,0)
			write(*,114) "   第一个坐标是：",xx(0,0,0),yy(0,0,0),zz(0,0,0)
			114 format (A,3(F10.5))
			write(*,114) "   第二个坐标是：",xx(1,0,0),yy(1,0,0),zz(1,0,0)
			write(*,114) "   第三个坐标是：",xx(2,0,0),yy(2,0,0),zz(2,0,0)
			write(*,*) 
		endif
		call PetscPrintf(comm," -----------------------------------\n",ierr)
		call PetscPrintf(comm,"        流场和网格数据已录入。      \n",ierr)
		call PetscPrintf(comm," -----------------------------------\n",ierr)
	end subroutine OhWeHaveThese

	subroutine read_mesh_3d(comm)
		implicit none
		PetscInt, INTENT(in) :: comm
		call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
		&                 DMDA_STENCIL_BOX, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, &
		&                 1, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, DA, ierr)
		call DMSetFromOptions(DA,ierr)
		call DMSetUp(DA, ierr)
		call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
		&                 DMDA_STENCIL_BOX, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,&
		&                 3, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, coordDA, ierr)
		call DMSetFromOptions(coordDA,ierr)
		call DMSetUp(coordDA, ierr)
		call DMGetGlobalVector(coordDA, Coord, ierr)
		call DMGetLocalVector(coordDA, Coord_local, ierr)
		call PetscViewerBinaryOpen(comm, trim(FileLocation)//"in//"//"grid.petsc",FILE_MODE_READ, Viewer, ierr)
		call VecLoad(Coord, Viewer, ierr)
		call DMGlobalToLocalBegin(coordDA, Coord, INSERT_VALUES, Coord_local, ierr)
		call DMGlobalToLocalEnd(coordDA, Coord, INSERT_VALUES, Coord_local, ierr)
		call PetscViewerDestroy(Viewer, ierr)
		call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC,& !DM_BOUNDARY_PERIODIC
		&                 DMDA_STENCIL_BOX, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,&
		&                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, meshDA, ierr)
		call DMSetFromOptions(meshDA,ierr)
		call DMSetUp(meshDA, ierr)
		call DMGetGlobalVector(meshDA, Flowfield, ierr)
		call DMGetLocalVector(meshDA, Flowfield_local, ierr)
		call PetscViewerBinaryOpen(comm, trim(FileLocation)//"in//"//"flow.petsc",FILE_MODE_READ, Viewer, ierr)
		call VecLoad(Flowfield, Viewer, ierr)
		call DMGlobalToLocalBegin(meshDA, Flowfield, INSERT_VALUES, Flowfield_local, ierr)
		call DMGlobalToLocalEnd(meshDA, Flowfield, INSERT_VALUES, Flowfield_local, ierr)
		call PetscViewerDestroy(Viewer, ierr)
	end subroutine read_mesh_3d

	subroutine get_layout()
		implicit none
		call DMDAGetGhostCorners(DA,igs,jgs,kgs,igl,jgl,kgl,ierr)
		ige=igs+igl-1; jge=jgs+jgl-1; kge=kgs+kgl-1
		call DMDAGetCorners(DA,is,js,ks,il,jl,kl,ierr)
		ie=is+il-1; je=js+jl-1; ke=ks+kl-1
	end subroutine get_layout

	subroutine load_mesh_info()
		implicit none
		integer :: i, j, k
		allocate(xx(igs:ige, jgs:jge, kgs:kge))
		allocate(yy(igs:ige, jgs:jge, kgs:kge))
		allocate(zz(igs:ige, jgs:jge, kgs:kge))
		call DMDAVecGetArrayReadF90(coordDA, coord_local, grid, ierr)
		do k=kgs, kge
			do j=jgs, jge
				do i=igs, ige
					xx(i, j, k)=real(grid(0, i, j, k))
					yy(i, j, k)=real(grid(1, i, j, k))
					zz(i, j, k)=real(grid(2, i, j, k))
				enddo
		  enddo
		enddo
		call DMDAVecRestoreArrayReadF90(coordDA, coord_local, grid, ierr)
	end subroutine load_mesh_info

	subroutine load_flow_info()
		  implicit none
		  integer :: i, j, k, l
		  allocate(qq(5, igs:ige, jgs:jge, kgs:kge))
		  call DMDAVecGetArrayReadF90(meshDA, Flowfield_local, flow, ierr)
		  do k=kgs, kge
			  do j=jgs, jge
				do i=igs, ige
				  do l=0,4
					qq(l+1, i, j, k)=real(flow(l, i, j, k))
				  enddo
				enddo
			  enddo
		  enddo
		  call DMDAVecRestoreArrayReadF90(meshDA, Flowfield_local, flow, ierr)
	end subroutine load_flow_info
end module mod_petsc_loader