!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_loading ! 读入并分发数据
	use petsc
	use mod_loaders
	use mod_petsc_viewer
	use mod_parameters
	use mod_cfgio_adapter
	implicit none
	private
	PetscScalar, pointer :: grid(:,:,:,:)
	PetscScalar, pointer :: flow(:,:,:,:)
	Vec :: Coord_local, Flowfield_local
	Vec :: Coord, Flowfield
	PetscErrorCode  :: ierr
	public :: loading_data
	PetscViewer :: Viewer
	contains
	subroutine loading_data(comm)
		implicit none
		PetscInt, intent(in) :: comm
		call read_from_file(comm)
		call bcast_parameters(comm)
		!call pack_parameters(comm)
		call read_mesh_3d(comm)
		call load_disturbance(comm)
		call get_layout()
		call load_mesh_info()
		call load_flow_info()
		call print_info(comm)
		call MPI_Barrier(comm,ierr)
	end subroutine loading_data

	subroutine read_from_file(comm)
		use mod_loaders
		use mod_petsc_viewer
		use mod_cfgio_adapter
		implicit none
		PetscInt,intent(in) :: comm
		character(len=256) :: cfg_file
		PetscBool :: set
		call signal_loading(comm)
		call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-f',cfg_file,set,ierr)
		if(.not. set) then
			write(*,*) 'should use -f option to determin the config file.'
			stop
		endif
		call cfg_loader(trim(cfg_file))
		if(rank==0)then
			call plot3d_load()
			call petsc_viewer(PETSC_COMM_SELF)
			call cfg_writer(trim(cfg_file))
		endif
		call MPI_Barrier(comm,ierr)
		call signal_starting(comm)
	end subroutine read_from_file

	subroutine bcast_parameters(comm)
		use mod_parameters
		implicit none 
		integer(KIND=MPI_ADDRESS_KIND) :: address_in,address_jn,address_kn,address_ln
		integer(KIND=MPI_ADDRESS_KIND) :: address_mode,address_Ma,address_Re,address_Te
		integer(KIND=MPI_ADDRESS_KIND) :: address_Alpha,address_Omega,address_Beta
		integer(KIND=MPI_ADDRESS_KIND) :: displacement(11)
		PetscInt,intent(in) :: comm
		integer :: block_lengths(11)
		integer :: pack_type
		integer :: types(11)
		call MPI_Get_address(in,address_in,ierr)
		call MPI_Get_address(jn,address_jn,ierr)
		call MPI_Get_address(kn,address_kn,ierr)
		call MPI_Get_address(ln,address_ln,ierr)
		call MPI_Get_address(lns_mode,address_mode,ierr)
		call MPI_Get_address(Ma,address_Ma,ierr)
		call MPI_Get_address(Re,address_Re,ierr)
		call MPI_Get_address(Te,address_Te,ierr)
		call MPI_Get_address(Alpha,address_Alpha,ierr)
		call MPI_Get_address(Beta,address_Beta,ierr)
		call MPI_Get_address(Omega,address_Omega,ierr)
		block_lengths=1
		displacement(1)=0
		displacement(2)=address_jn-address_in
		displacement(3)=address_kn-address_in
		displacement(4)=address_ln-address_in
		displacement(5)=address_mode-address_in
		displacement(6)=address_Ma-address_in
		displacement(7)=address_Re-address_in
		displacement(8)=address_Te-address_in
		displacement(9)=address_Alpha-address_in
		displacement(10)=address_Beta-address_in
		displacement(11)=address_Omega-address_in
		types=(/MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,&
			MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,&
			MPI_DOUBLE_COMPLEX,MPI_DOUBLE_COMPLEX,MPI_DOUBLE_COMPLEX/)
		call MPI_Type_create_struct(11,block_lengths,displacement,types,pack_type,ierr)
		call MPI_Type_commit(pack_type,ierr)
		call MPI_Bcast(in,1,pack_type,0,comm,ierr)
		call MPI_Barrier(comm,ierr)
		call MPI_Type_free(pack_type,ierr)
	end subroutine bcast_parameters

	subroutine pack_parameters(comm)
		use mod_parameters
		implicit none 
		integer,intent(in) :: comm 
		integer :: packsize,position
		character(len=120) :: packbuf
		if(rank==0)then 
			position = 0
			call MPI_Pack(in,1,MPI_INT,packbuf,120,position,comm,ierr)
			call MPI_Pack(jn,1,MPI_INT,packbuf,120,position,comm,ierr)
			call MPI_Pack(kn,1,MPI_INT,packbuf,120,position,comm,ierr)
			call MPI_Pack(lns_mode,1,MPI_INT,packbuf,120,position,comm,ierr)
			call MPI_Pack(Ma,1,MPI_DOUBLE,packbuf,120,position,comm,ierr)
			call MPI_Pack(Re,1,MPI_DOUBLE,packbuf,120,position,comm,ierr)
			call MPI_Pack(Te,1,MPI_DOUBLE,packbuf,120,position,comm,ierr)
			call MPI_Pack(Alpha,1,MPI_DOUBLE_COMPLEX,packbuf,120,position,comm,ierr)
			call MPI_Pack(Beta,1,MPI_DOUBLE_COMPLEX,packbuf,120,position,comm,ierr)
			call MPI_Pack(Omega,1,MPI_DOUBLE_COMPLEX,packbuf,120,position,comm,ierr)
		endif 
		call MPI_Bcast(packbuf,120,MPI_PACKED,0,comm,ierr)
		if(rank/=0)then 
			position = 0
			call MPI_Unpack(packbuf,120,position,in,1,MPI_INT,comm,ierr)
			call MPI_Unpack(packbuf,120,position,jn,1,MPI_INT,comm,ierr)
			call MPI_Unpack(packbuf,120,position,kn,1,MPI_INT,comm,ierr)
			call MPI_Unpack(packbuf,120,position,lns_mode,1,MPI_INT,comm,ierr)
			call MPI_Unpack(packbuf,120,position,Ma,1,MPI_DOUBLE,comm,ierr)
			call MPI_Unpack(packbuf,120,position,Re,1,MPI_DOUBLE,comm,ierr)
			call MPI_Unpack(packbuf,120,position,Te,1,MPI_DOUBLE,comm,ierr)
			call MPI_Unpack(packbuf,120,position,Alpha,1,MPI_DOUBLE_COMPLEX,comm,ierr)
			call MPI_Unpack(packbuf,120,position,Beta,1,MPI_DOUBLE_COMPLEX,comm,ierr)
			call MPI_Unpack(packbuf,120,position,Omega,1,MPI_DOUBLE_COMPLEX,comm,ierr)
		endif 
		call MPI_Barrier(comm,ierr)
	end subroutine pack_parameters

	subroutine read_mesh_3d(comm)
		implicit none
		PetscInt, intent(in) :: comm
		call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
		&                 DMDA_STENCIL_BOX, in, jn, kn, 1, PETSC_DECIDE, PETSC_DECIDE, &
		&                 1, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, DA, ierr)
		call DMSetFromOptions(DA,ierr)
		call DMSetUp(DA, ierr)
		call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
		&                 DMDA_STENCIL_BOX, in, jn, kn, 1, PETSC_DECIDE, PETSC_DECIDE,&
		&                 3, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, coordDA, ierr)
		call DMSetFromOptions(coordDA,ierr)
		call DMSetUp(coordDA, ierr)
		call DMGetGlobalVector(coordDA, Coord, ierr)
		call DMGetLocalVector(coordDA, Coord_local, ierr)
		call PetscViewerBinaryOpen(comm, "in/grid.petsc",FILE_MODE_READ, Viewer, ierr)
		call VecLoad(Coord, Viewer, ierr)
		call DMGlobalToLocalBegin(coordDA, Coord, INSERT_VALUES, Coord_local, ierr)
		call DMGlobalToLocalEnd(coordDA, Coord, INSERT_VALUES, Coord_local, ierr)
		call PetscViewerDestroy(Viewer, ierr)
		call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC,&
		&                 DMDA_STENCIL_BOX, in, jn, kn, 1, PETSC_DECIDE, PETSC_DECIDE,&
		&                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, meshDA, ierr)
		call DMSetFromOptions(meshDA,ierr)
		call DMSetUp(meshDA, ierr)
		call DMGetGlobalVector(meshDA, Flowfield, ierr)
		call DMGetLocalVector(meshDA, Flowfield_local, ierr)
		call PetscViewerBinaryOpen(comm, "in/flow.petsc",FILE_MODE_READ, Viewer, ierr)
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

	subroutine load_disturbance(comm)
		implicit none 
		PetscInt,intent(in) :: comm 
		select case (lns_mode)
		case(0)
			call load_disturbance_2d(comm)
		case(1)
			call load_disturbance_3d(comm)
		end select
	end subroutine load_disturbance

	subroutine load_disturbance_2d(comm)
		implicit none 
		PetscInt,intent(in) :: comm
		call DMDACreate1d(comm, DM_BOUNDARY_NONE, jn, 5, 2, PETSC_NULL_INTEGER, disturbDA, ierr)
		call DMSetFromOptions(disturbDA,ierr)
		call DMSetUp(disturbDA, ierr)
		call DMGetGlobalVector(disturbDA, disturb, ierr)
        call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "in//LPSE_disturbance+1.petsc", FILE_MODE_READ, Viewer, ierr)
        call VecLoad(disturb, Viewer, ierr)
        call PetscViewerDestroy(Viewer, ierr)
	end subroutine load_disturbance_2d

	subroutine load_disturbance_3d(comm)
		implicit none 
		PetscInt,intent(in) :: comm 
		call DMDACreate2d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
		&                 DMDA_STENCIL_BOX, jn, kn, PETSC_DECIDE, PETSC_DECIDE, &
		&                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, disturbDA, ierr)
		call DMSetFromOptions(disturbDA,ierr)
		call DMSetUp(disturbDA, ierr)
		call DMGetGlobalVector(disturbDA, disturb, ierr)
        call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "in//LPSE_disturbance+1.petsc", FILE_MODE_READ, Viewer, ierr)
        call VecLoad(disturb, Viewer, ierr)
        call PetscViewerDestroy(Viewer, ierr)
	end subroutine load_disturbance_3d

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

	subroutine signal_loading(comm)
		implicit none
		PetscInt,intent(in) :: comm
		call PetscPrintf(comm, "\n", ierr)
		call PetscPrintf(comm, " ===========================================================================\n", ierr)
		call PetscPrintf(comm, " =                                 读    取                                = \n", ierr)
		call PetscPrintf(comm, " ===========================================================================\n", ierr)
		call PetscPrintf(comm, " 「 单 进 程 」\n",ierr)
	end subroutine signal_loading

	subroutine signal_starting(comm)
		implicit none
		PetscInt,intent(in) :: comm 
		call PetscPrintf(comm, "\n", ierr)
		call PetscPrintf(comm, " ===========================================================================\n", ierr)
		call PetscPrintf(comm, " =                                 计    算                                = \n", ierr)
		call PetscPrintf(comm, " ===========================================================================\n", ierr)
		call PetscPrintf(comm, " 「 多 进 程 」\n",ierr)
	end subroutine signal_starting

	subroutine print_info(comm)
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
			write(*,9) in,jn,kn,ln
			9 format ('   流向的网格数in=',I5,/,'   法向的网格数jn=',I5,/,&
				'   展向的网格数kn=',I5/,'   自由度ln=',I5)
			write(*,113) "   第一个数据是：",qq(1,0,0,0),qq(2,0,0,0),qq(3,0,0,0),qq(4,0,0,0),qq(5,0,0,0)
			113 format (A,5(F10.5))
			write(*,113) "   第二个数据是：",qq(1,1,0,0),qq(2,1,0,0),qq(3,1,0,0),qq(4,1,0,0),qq(5,1,0,0)
			write(*,113) "   第三个数据是：",qq(1,2,0,0),qq(2,2,0,0),qq(3,2,0,0),qq(4,2,0,0),qq(5,2,0,0)
			write(*,114) "   第一个坐标是：",xx(0,0,0),yy(0,0,0),zz(0,0,0)
			114 format (A,3(F10.5))
			write(*,114) "   第二个坐标是：",xx(1,0,0),yy(1,0,0),zz(1,0,0)
			write(*,114) "   第三个坐标是：",xx(2,0,0),yy(2,0,0),zz(2,0,0)
		endif
		call PetscPrintf(comm,"\n",ierr)
		call PetscPrintf(comm," -----------------------------------\n",ierr)
		call PetscPrintf(comm,"        流场和网格数据已录入。      \n",ierr)
		call PetscPrintf(comm," -----------------------------------\n",ierr)
	end subroutine print_info
end module mod_loading
