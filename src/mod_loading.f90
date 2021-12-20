!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_loading ! 读入并分发数据
	use mod_parameters
	use mod_reading
	use petsc
	public :: loading_data
	private
	PetscErrorCode  :: ierr
	Vec :: Flowfield_local
	Vec :: Coord,Flowfield	
	PetscViewer :: Viewer
	Vec :: Multi_disturb
	Vec :: Coord_local
	contains
	subroutine loading_data(comm)
		implicit none
		PetscInt, intent(in) :: comm

		call read_argv_and_file(comm)

		call signal_mpiing(comm)
		call bcast_parameters(comm)
		! call pack_parameters(comm)
		call set_mpi_da(comm)

		call signal_loading(comm)
		call load_petsc_file(comm)
		call get_layout()
		call load_disturb_mesh_flow()
		call deallocate_memory()
		call MPI_Barrier(comm,ierr)

		call signal_printing(comm)
		call print_info()
	end subroutine loading_data

	subroutine read_argv_and_file(comm)
		use mod_cfgio_adapter
		implicit none
		PetscInt,intent(in) :: comm
		character(len=256) :: cfg_file
		PetscBool :: set

		call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-f',cfg_file,set,ierr)
		if(.not. set) then
			write(*,*) 'should use -f option to determin the config file.'
			stop
		endif

		call cfg_loader(trim(cfg_file))

		if(rank==0)then
			call load(PETSC_COMM_SELF)
			call cfg_writer(trim(cfg_file))
		endif

		call MPI_Barrier(comm,ierr)
	end subroutine read_argv_and_file

	subroutine signal_mpiing(comm)
		implicit none
		PetscInt,intent(in) :: comm 
		call PetscPrintf(comm, "\n", ierr)
		call PetscPrintf(comm, " ===========================================================================\n", ierr)
		call PetscPrintf(comm, " =                                 计    算                                = \n", ierr)
		call PetscPrintf(comm, " ===========================================================================\n", ierr)
		call PetscPrintf(comm, " 「 M P I 」\n",ierr)
	end subroutine signal_mpiing

	subroutine bcast_parameters(comm)
		implicit none 
		integer(KIND=MPI_ADDRESS_KIND) :: address_in,address_jn,address_kn,address_ln
		integer(KIND=MPI_ADDRESS_KIND) :: address_mode,address_Ma,address_Re,address_Te
		integer(KIND=MPI_ADDRESS_KIND) :: address_Alpha,address_Omega,address_Beta
		integer(KIND=MPI_ADDRESS_KIND) :: address_initguess
		integer(KIND=MPI_ADDRESS_KIND) :: displacement(12)
		integer :: block_lengths(12)
		PetscInt,intent(in) :: comm
		integer :: pack_type
		integer :: types(12)

		call MPI_Get_address(in,address_in,ierr)
		call MPI_Get_address(jn,address_jn,ierr)
		call MPI_Get_address(kn,address_kn,ierr)
		call MPI_Get_address(ln,address_ln,ierr)
		call MPI_Get_address(lns_mode,address_mode,ierr)
		call MPI_Get_address(initial_guess,address_initguess,ierr)
		call MPI_Get_address(Ma,address_Ma,ierr)
		call MPI_Get_address(Re,address_Re,ierr)
		call MPI_Get_address(Te,address_Te,ierr)
		call MPI_Get_address(Alpha,address_Alpha,ierr)
		call MPI_Get_address(Beta,address_Beta,ierr)
		call MPI_Get_address(Omega,address_Omega,ierr)

		displacement(1)=0
		displacement(2)=address_jn-address_in
		displacement(3)=address_kn-address_in
		displacement(4)=address_ln-address_in
		displacement(5)=address_mode-address_in
		displacement(6)=address_initguess-address_in
		displacement(7)=address_Ma-address_in
		displacement(8)=address_Re-address_in
		displacement(9)=address_Te-address_in
		displacement(10)=address_Alpha-address_in
		displacement(11)=address_Beta-address_in
		displacement(12)=address_Omega-address_in

		block_lengths=1

		types=(/MPI_INTEGER4,MPI_INTEGER4,MPI_INTEGER4,MPI_INTEGER4,MPI_INTEGER4,MPI_LOGICAL,&
				MPI_REAL8,MPI_REAL8,MPI_REAL8,&
				MPI_COMPLEX16,MPI_COMPLEX16,MPI_COMPLEX16/)

		call MPI_Type_create_struct(12,block_lengths,displacement,types,pack_type,ierr)
		call MPI_Type_commit(pack_type,ierr)
		call MPI_Bcast(in,1,pack_type,0,comm,ierr)
		call MPI_Barrier(comm,ierr)
		call MPI_Type_free(pack_type,ierr)

	end subroutine bcast_parameters

	subroutine pack_parameters(comm)
		implicit none 
		character(len=120) :: packbuf
		integer :: packsize,position
		integer,intent(in) :: comm 

		if(rank==0)then 
			position = 0
			call MPI_Pack(in,1,MPI_INTEGER4,packbuf,120,position,comm,ierr)
			call MPI_Pack(jn,1,MPI_INTEGER4,packbuf,120,position,comm,ierr)
			call MPI_Pack(kn,1,MPI_INTEGER4,packbuf,120,position,comm,ierr)
			call MPI_Pack(ln,1,MPI_INTEGER4,packbuf,120,position,comm,ierr)
			call MPI_Pack(lns_mode,1,MPI_INTEGER4,packbuf,120,position,comm,ierr)
			call MPI_Pack(initial_guess,1,MPI_LOGICAL,packbuf,120,position,comm,ierr)
			call MPI_Pack(Ma,1,MPI_REAL8,packbuf,120,position,comm,ierr)
			call MPI_Pack(Re,1,MPI_REAL8,packbuf,120,position,comm,ierr)
			call MPI_Pack(Te,1,MPI_REAL8,packbuf,120,position,comm,ierr)
			call MPI_Pack(Alpha,1,MPI_COMPLEX16,packbuf,120,position,comm,ierr)
			call MPI_Pack(Beta,1,MPI_COMPLEX16,packbuf,120,position,comm,ierr)
			call MPI_Pack(Omega,1,MPI_COMPLEX16,packbuf,120,position,comm,ierr)
		endif 
		call MPI_Bcast(packbuf,120,MPI_PACKED,0,comm,ierr)

		if(rank/=0)then 
			position = 0
			call MPI_Unpack(packbuf,120,position,in,1,MPI_INTEGER4,comm,ierr)
			call MPI_Unpack(packbuf,120,position,jn,1,MPI_INTEGER4,comm,ierr)
			call MPI_Unpack(packbuf,120,position,kn,1,MPI_INTEGER4,comm,ierr)
			call MPI_Unpack(packbuf,120,position,ln,1,MPI_INTEGER4,comm,ierr)
			call MPI_Unpack(packbuf,120,position,lns_mode,1,MPI_INTEGER4,comm,ierr)
			call MPI_Unpack(packbuf,120,position,initial_guess,1,MPI_LOGICAL,comm,ierr)
			call MPI_Unpack(packbuf,120,position,Ma,1,MPI_REAL8,comm,ierr)
			call MPI_Unpack(packbuf,120,position,Re,1,MPI_REAL8,comm,ierr)
			call MPI_Unpack(packbuf,120,position,Te,1,MPI_REAL8,comm,ierr)
			call MPI_Unpack(packbuf,120,position,Alpha,1,MPI_COMPLEX16,comm,ierr)
			call MPI_Unpack(packbuf,120,position,Beta,1,MPI_COMPLEX16,comm,ierr)
			call MPI_Unpack(packbuf,120,position,Omega,1,MPI_COMPLEX16,comm,ierr)
		endif 
		call MPI_Barrier(comm,ierr)

	end subroutine pack_parameters

	subroutine set_mpi_da(comm)
		implicit none
		PetscInt, intent(in) :: comm

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

		call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC,&
		&                 DMDA_STENCIL_BOX, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,&
		&                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, meshDA, ierr)
		call DMSetFromOptions(meshDA,ierr)
		call DMSetUp(meshDA, ierr)

        call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
        &                 DMDA_STENCIL_BOX, sink, jn, kn, PETSC_DECIDE, 1, 1,&
        &                 5, 0, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, disturbDA, ierr)
        call DMSetUp(disturbDA, ierr)

	end subroutine set_mpi_da

	subroutine signal_loading(comm)
		implicit none 
		PetscInt,intent(in) :: comm 
		call PetscPrintf(comm," -----------------------------------\n",ierr)
		call PetscPrintf(comm,"               开始分发              \n",ierr)
		call PetscPrintf(comm," -----------------------------------\n",ierr)
	end subroutine signal_loading

	subroutine load_petsc_file(comm)
		implicit none 
		PetscInt, intent(in) :: comm

		call DMGetGlobalVector(coordDA, Coord, ierr)
		call DMGetLocalVector(coordDA, Coord_local, ierr)
		call PetscViewerBinaryOpen(comm, "in/grid.petsc",FILE_MODE_READ, Viewer, ierr)
		call VecLoad(Coord, Viewer, ierr)
		call DMGlobalToLocalBegin(coordDA, Coord, INSERT_VALUES, Coord_local, ierr)
		call DMGlobalToLocalEnd(coordDA, Coord, INSERT_VALUES, Coord_local, ierr)
		call PetscViewerDestroy(Viewer, ierr)

		call DMGetGlobalVector(meshDA, Flowfield, ierr)
		call DMGetLocalVector(meshDA, Flowfield_local, ierr)
		call PetscViewerBinaryOpen(comm, "in/flow.petsc",FILE_MODE_READ, Viewer, ierr)
		call VecLoad(Flowfield, Viewer, ierr)
		call DMGlobalToLocalBegin(meshDA, Flowfield, INSERT_VALUES, Flowfield_local, ierr)
		call DMGlobalToLocalEnd(meshDA, Flowfield, INSERT_VALUES, Flowfield_local, ierr)
		call PetscViewerDestroy(Viewer, ierr)

        call DMGetGlobalVector(disturbDA, Multi_disturb, ierr)
		call PetscViewerBinaryOpen(comm, "in/disturb.petsc",FILE_MODE_READ, Viewer, ierr)
		call VecLoad(Multi_disturb, Viewer, ierr)
		call PetscViewerDestroy(Viewer, ierr)

		call DMGetGlobalVector(meshDA, Turtle, ierr)
		select case (initial_guess)
		case(.True.)
			call VecZeroEntries(Turtle,ierr)
			call PetscViewerBinaryOpen(comm, trim(initfile),FILE_MODE_READ, Viewer, ierr)
			call VecLoad(Turtle, Viewer, ierr)
			call PetscViewerDestroy(Viewer, ierr)
		case(.False.)
			call VecZeroEntries(Turtle,ierr)
		end select

	end subroutine load_petsc_file

	subroutine get_layout()
		implicit none

		call DMDAGetGhostCorners(DA,igs,jgs,kgs,igl,jgl,kgl,ierr)
		call DMDAGetCorners(DA,is,js,ks,il,jl,kl,ierr)

		ige=igs+igl-1; jge=jgs+jgl-1; kge=kgs+kgl-1
		ie=is+il-1; je=js+jl-1; ke=ks+kl-1

	end subroutine get_layout

	subroutine load_disturb_mesh_flow()
		implicit none 
		integer :: xs,ys,zs,xl,yl,zl,xe,ye,ze
		PetscScalar, pointer :: multi(:,:,:,:)
		PetscScalar, pointer :: grid(:,:,:,:)
		PetscScalar, pointer :: flow(:,:,:,:)
		integer :: l,i,j,k 

		call DMDAGetCorners(disturbDA,xs,ys,zs,xl,yl,zl,ierr)
		xe=xs+xl-1;ye=ys+yl-1;ze=zs+zl-1
		allocate(disturb(0:4, ys:ye, zs:ze))
		call DMDAVecGetArrayReadF90(disturbDA, Multi_disturb, multi, ierr)
		do i=xs,xe
			do j=ys,ye
				do k=zs,ze
					disturb(:,j,k) = multi(:,i,j,k)
				enddo
			enddo
		enddo
		call DMDAVecRestoreArrayReadF90(disturbDA, Multi_disturb, multi, ierr)

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
	end subroutine load_disturb_mesh_flow

	subroutine deallocate_memory()
		call VecDestroy(Flowfield_local,ierr)
		call VecDestroy(Multi_disturb,ierr)
		call VecDestroy(Coord_local,ierr)
		call VecDestroy(Flowfield,ierr)
		call VecDestroy(Coord,ierr)
	end subroutine deallocate_memory

	subroutine signal_printing(comm)
		implicit none 
		PetscInt,intent(in) :: comm 
		call PetscPrintf(comm," -----------------------------------\n",ierr)
		call PetscPrintf(comm,"               分发结束              \n",ierr)
		call PetscPrintf(comm," -----------------------------------\n",ierr)
		call PetscPrintf(comm," -----------------------------------\n",ierr)
		call PetscPrintf(comm,"               查对信息              \n",ierr)
		call PetscPrintf(comm," -----------------------------------\n",ierr)
		call PetscPrintf(comm,"\n",ierr)
	end subroutine signal_printing

	subroutine print_info()
		implicit none 
		if(rank==0)then
			write(*,"(A,I5)") "   进程数 =",sink
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
	end subroutine print_info
end module mod_loading
