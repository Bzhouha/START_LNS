!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_files
!
! 这个模块处理文件流，包括:
!
!                 1).从config配置文件读入参数；
!
!                 2).支持 raw Plot3d格式、binary格式、hdf5格式文件的输入输出
!
!   call istream(comm)  读入配置文件，初始化PETSc以及从文件载入流场、网格、边界、初值等。
!
!   call ostream(comm) 将计算结果以Binary、HDF5格式输出。
!
    use mod_parameters
    use petsc
    public :: istream,ostream
    private
    PetscErrorCode  :: ierr
    PetscViewer :: viewer
    contains
    subroutine istream(comm)
        implicit none
        PetscInt, intent(in) :: comm

        call config(comm)

        call set_DM(comm)

        call load(comm)

    end subroutine istream

    subroutine config(comm)
        use mod_parameters,only : fk
        use mod_cfgio_adapter
        implicit none
        PetscInt,intent(in) :: comm
        character(len=256) :: cfg_file
        logical :: ksp_flg,snes_flg
        PetscBool :: set

        call mpi_comm_rank(comm,rank,ierr)
        call mpi_comm_size(comm,sink,ierr)
        call PetscPrintf(comm, "\n         S T A R T - L N S\n", ierr)
        call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-f',cfg_file,set,ierr)
        if(.not. set) then
            write(*,*) 'should use -f option to determin the config file.'
            stop
        endif

        call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-ksp',ksp_flg,ierr)
        call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-snes',snes_flg,ierr)
        if(ksp_flg)then
            solver_mode='ksp'; split_mode=0
        endif
        if(snes_flg)then
            solver_mode='snes';split_mode=1
        endif

        call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-raw',set,ierr)
        if(set) io_type="raw"
        call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-binary',set,ierr)
        if(set) io_type="binary"
        call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-hdf5',set,ierr)
        if(set) io_type="hdf5"
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-lk',fk,set,ierr)

        if(ksp_flg)  call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_monitor",PETSC_NULL_CHARACTER,ierr)
        if(snes_flg) call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-snes_monitor",PETSC_NULL_CHARACTER,ierr)
        ! if(snes_flg) call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_monitor",PETSC_NULL_CHARACTER,ierr)
        ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-sub_pc_factor_in_place",PETSC_NULL_CHARACTER,ierr)
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-pc_asm_sub_mat_type","baij",ierr)

        call cfg_loader(trim(cfg_file))
        call MPI_Barrier(comm,ierr)

        call bcast(comm)

    end subroutine config

    subroutine bcast(comm)
        implicit none
        integer(KIND=MPI_ADDRESS_KIND) :: address_in,address_jn,address_kn,address_ln
        integer(KIND=MPI_ADDRESS_KIND) :: address_mode,address_Ma,address_Re,address_Te
        integer(KIND=MPI_ADDRESS_KIND) :: address_Alpha,address_Omega,address_Beta
        integer(KIND=MPI_ADDRESS_KIND) :: address_initguess,address_h5fn
        integer(KIND=MPI_ADDRESS_KIND) :: displacement(13)
        integer :: block_lengths(13)
        PetscInt,intent(in) :: comm
        integer :: pack_type
        integer :: types(13)

        call MPI_Get_address(in,address_in,ierr)
        call MPI_Get_address(jn,address_jn,ierr)
        call MPI_Get_address(kn,address_kn,ierr)
        call MPI_Get_address(ln,address_ln,ierr)
        call MPI_Get_address(lns_mode,address_mode,ierr)
        call MPI_Get_address(init_guess_flg,address_initguess,ierr)
        call MPI_Get_address(Ma,address_Ma,ierr)
        call MPI_Get_address(Re,address_Re,ierr)
        call MPI_Get_address(Te,address_Te,ierr)
        call MPI_Get_address(Alpha,address_Alpha,ierr)
        call MPI_Get_address(Beta,address_Beta,ierr)
        call MPI_Get_address(Omega,address_Omega,ierr)
        call MPI_Get_address(hdf5file,address_h5fn,ierr)
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
        displacement(13)=address_h5fn-address_in
        block_lengths=1
        block_lengths(13)=256
        types=(/MPI_INTEGER4,MPI_INTEGER4,MPI_INTEGER4,MPI_INTEGER4,MPI_INTEGER4,&
                MPI_LOGICAL,MPI_REAL8,MPI_REAL8,MPI_REAL8,&
                MPI_COMPLEX16,MPI_COMPLEX16,MPI_COMPLEX16,MPI_CHAR/)
        call MPI_Type_create_struct(13,block_lengths,displacement,types,pack_type,ierr)
        call MPI_Type_commit(pack_type,ierr)
        call MPI_Bcast(in,1,pack_type,0,comm,ierr)
        call MPI_Barrier(comm,ierr)
        call MPI_Type_free(pack_type,ierr)
    end subroutine bcast

    subroutine set_DM(comm)
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
        call DMSetMatType(meshDA,MATBAIJ,ierr)
        call DMSetFromOptions(meshDA,ierr)
        call DMSetUp(meshDA, ierr)

        call DMGetGlobalVector(meshDA, turtle, ierr)
        call DMDAGetGhostCorners(DA,igs,jgs,kgs,igl,jgl,kgl,ierr)
        call DMDAGetCorners(DA,is,js,ks,il,jl,kl,ierr)
        ige=igs+igl-1; jge=jgs+jgl-1; kge=kgs+kgl-1
        ie=is+il-1; je=js+jl-1; ke=ks+kl-1
    end subroutine set_DM

    subroutine load(comm)
        implicit none
        integer,intent(in) :: comm

        call PetscPrintf(comm, "\n ----------------------------------\n", ierr)
        call PetscPrintf(comm, "              IStream            \n", ierr)

        select case(io_type)
            case("raw")
                call PetscPrintf(comm, "\n   I/O type : raw\n\n", ierr)
                call load_raw_files(comm)
            case("binary")
                call PetscPrintf(comm, "\n   I/O type : binary\n\n", ierr)
                call load_binary_files(comm)
            case("hdf5")
                call PetscPrintf(comm, "\n   I/O type : hdf5\n\n", ierr)
                call load_hdf5_files(comm)
        end select

        call set_disturb(comm)
        call set_init_guess(comm)
        call check(comm)
    end subroutine load

    subroutine load_raw_files(comm)
        implicit none
        integer,intent(in) :: comm

        bigridfile = "./data/grid.pet"
        biflowfile = "./data/flow.pet"
        if(rank==0)then
            call raw_to_binary()
        endif
        call MPI_Barrier(comm,ierr)
        call load_binary_files(comm)
    end subroutine load_raw_files

    subroutine raw_to_binary()
        implicit none
        real(R_P),dimension(:,:,:,:),allocatable :: qq_0
        PetscScalar, pointer :: disturbs(:,:,:,:)
        PetscScalar, pointer :: grid(:,:,:,:)
        PetscScalar, pointer :: flow(:,:,:,:)
        PetscScalar, pointer :: slice(:,:,:)
        integer :: xs,ys,zs,xl,yl,zl
        DM :: uni_coordDA,uni_meshDA
        Vec :: coord,flowfield
        integer :: l,i,j,k

        ! 读取网格信息
        open(11, file=trim(gridfile),action='read',form='unformatted')
        read(11)
        select case (lns_mode)
        case(2)
            read(11) in,jn
            allocate(xx(in,jn,kn), yy(in,jn,kn), zz(in,jn,kn))
            read(11) xx,yy
            zz=0.0d0
        case(3)
            read(11) in,jn,kn
            allocate(xx(in,jn,kn), yy(in,jn,kn), zz(in,jn,kn))
            read(11) xx,yy,zz
        end select
        close(11)

        ! 读取基本流数据
        open(12, file=trim(flowfile),action='read',form='unformatted')
        read(12)
        select case (lns_mode)
        case(2)
            read(12) in,jn,ln
        case(3)
            read(12) in,jn,kn,ln
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

        call DMDACreate3d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
        &                 DMDA_STENCIL_BOX, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,&
        &                 3, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, uni_coordDA, ierr)
        call DMSetUp(uni_coordDA, ierr)

        call DMDACreate3d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
        &                 DMDA_STENCIL_BOX, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,&
        &                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, uni_meshDA, ierr)
        call DMSetUp(uni_meshDA, ierr)

        call DMGetGlobalVector(uni_coordDA,coord,ierr)
        call DMDAGetCorners(uni_coordDA,xs,ys,zs,xl,yl,zl,ierr)
        call DMDAVecGetArrayF90(uni_coordDA,coord,grid,ierr)
        do i=xs,xs+xl-1
            do j=ys,ys+yl-1
                do k=zs,zs+zl-1
                    grid(0, i, j, k) = xx(i+1, j+1, k+1)
                    grid(1, i, j, k) = yy(i+1, j+1, k+1)
                    grid(2, i, j, k) = zz(i+1, j+1, k+1)
                enddo
            enddo
        enddo
        call DMDAVecRestoreArrayF90(uni_coordDA,coord,grid,ierr)

        call DMGetGlobalVector(uni_meshDA,flowfield,ierr)
        call DMDAGetCorners(uni_meshDA,xs,ys,zs,xl,yl,zl,ierr)
        call DMDAVecGetArrayF90(uni_meshDA,flowfield,flow,ierr)
        do i=xs,xs+xl-1
            do j=ys,ys+yl-1
                do k=zs,zs+zl-1
                    do l=0, 4
                        flow(l,i,j,k) = qq(l+1, i+1, j+1, k+1)
                    enddo
                enddo
            enddo
        enddo
        call DMDAVecRestoreArrayF90(uni_meshDA,flowfield,flow,ierr)

        call PetscViewerBinaryOpen(PETSC_COMM_SELF, trim(bigridfile),FILE_MODE_WRITE, viewer, ierr)
        call VecView(coord, viewer, ierr)
        call PetscViewerDestroy(viewer, ierr)

        call PetscViewerBinaryOpen(PETSC_COMM_SELF, trim(biflowfile),FILE_MODE_WRITE, viewer, ierr)
        call VecView(flowfield, viewer, ierr)
        call PetscViewerDestroy(viewer, ierr)

        deallocate(xx)
        deallocate(yy)
        deallocate(zz)
        deallocate(qq)
        call DMRestoreGlobalVector(uni_meshDA,flowfield,ierr)
        call DMRestoreGlobalVector(uni_coordDA,coord,ierr)
        call DMDestroy(uni_coordDA,ierr)
        call DMDestroy(uni_meshDA,ierr)

    end subroutine raw_to_binary

    subroutine load_binary_files(comm)
        implicit none
        PetscScalar, pointer :: grid(:,:,:,:)
        PetscScalar, pointer :: flow(:,:,:,:)
        Vec :: flowfield_local,coord_local
        PetscInt, intent(in) :: comm
        Vec :: coord,flowfield
        integer :: l,i,j,k

        call DMGetGlobalVector(coordDA, coord, ierr)
        call DMGetLocalVector(coordDA, coord_local, ierr)
        call PetscViewerBinaryOpen(comm, trim(bigridfile),FILE_MODE_READ, viewer, ierr)
        call VecLoad(coord, viewer, ierr)
        call DMGlobalToLocalBegin(coordDA, coord, INSERT_VALUES, coord_local, ierr)
        call DMGlobalToLocalEnd(coordDA, coord, INSERT_VALUES, coord_local, ierr)
        call PetscViewerDestroy(viewer, ierr)

        call DMGetGlobalVector(meshDA, flowfield, ierr)
        call DMGetLocalVector(meshDA, flowfield_local, ierr)
        call PetscViewerBinaryOpen(comm, trim(biflowfile),FILE_MODE_READ, viewer, ierr)
        call VecLoad(flowfield, viewer, ierr)
        call DMGlobalToLocalBegin(meshDA, flowfield, INSERT_VALUES, flowfield_local, ierr)
        call DMGlobalToLocalEnd(meshDA, flowfield, INSERT_VALUES, flowfield_local, ierr)
        call PetscViewerDestroy(viewer, ierr)

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
        call DMDAVecGetArrayReadF90(meshDA, flowfield_local, flow, ierr)
        do k=kgs, kge
            do j=jgs, jge
                do i=igs, ige
                    do l=0,4
                        qq(l+1, i, j, k)=real(flow(l, i, j, k))
                    enddo
                enddo
            enddo
        enddo
        call DMDAVecRestoreArrayReadF90(meshDA, flowfield_local, flow, ierr)

        call DMRestoreLocalVector(meshDA,flowfield_local,ierr)
        call DMRestoreLocalVector(coordDA,coord_local,ierr)
        call DMRestoreGlobalVector(meshDA,flowfield,ierr)
        call DMRestoreGlobalVector(coordDA,coord,ierr)

    end subroutine load_binary_files

    subroutine load_hdf5_files(comm)
        implicit none
        PetscScalar, pointer :: grid(:,:,:,:)
        PetscScalar, pointer :: flow(:,:,:,:)
        Vec :: flowfield_local,coord_local
        PetscInt, intent(in) :: comm
        Vec :: coord,flowfield
        integer :: l,i,j,k

        call DMGetGlobalVector(meshDA, flowfield, ierr)
        call PetscObjectSetName(flowfield,"baseflow",ierr)
        call DMGetGlobalVector(coordDA, coord, ierr)
        call PetscObjectSetName(coord,"grid",ierr)

        call PetscViewerHDF5Open(comm,trim(hdf5file),FILE_MODE_READ,viewer,ierr)
        call VecLoad(flowfield,viewer,ierr)
        call VecLoad(coord,viewer,ierr)
        call PetscViewerDestroy(viewer, ierr)

        call DMGetLocalVector(coordDA, coord_local, ierr)
        call DMGetLocalVector(meshDA, flowfield_local, ierr)
        call DMGlobalToLocalBegin(coordDA, coord, INSERT_VALUES, coord_local, ierr)
        call DMGlobalToLocalBegin(meshDA, flowfield, INSERT_VALUES, flowfield_local, ierr)
        call DMGlobalToLocalEnd(coordDA, coord, INSERT_VALUES, coord_local, ierr)
        call DMGlobalToLocalEnd(meshDA, flowfield, INSERT_VALUES, flowfield_local, ierr)

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
        call DMDAVecGetArrayReadF90(meshDA, flowfield_local, flow, ierr)
        do k=kgs, kge
            do j=jgs, jge
                do i=igs, ige
                    do l=0,4
                        qq(l+1, i, j, k)=real(flow(l, i, j, k))
                    enddo
                enddo
            enddo
        enddo
        call DMDAVecRestoreArrayReadF90(meshDA, flowfield_local, flow, ierr)

        call DMRestoreLocalVector(meshDA,flowfield_local,ierr)
        call DMRestoreLocalVector(coordDA,coord_local,ierr)
        call DMRestoreGlobalVector(meshDA,flowfield,ierr)
        call DMRestoreGlobalVector(coordDA,coord,ierr)

    end subroutine load_hdf5_files

    subroutine set_disturb(comm)
        implicit none
        PetscScalar, pointer :: disturbs(:,:,:,:)
        integer :: xs,ys,zs,xl,yl,zl,xe,ye,ze
        PetscScalar, pointer :: slice(:,:,:)
        integer,intent(in) :: comm
        Vec :: disturb_gather
        Vec :: disturb_slice
        DM :: disturbDA
        integer :: i,j,k
        DM :: sliceDA

        turbfiles = "./data/disturbs.pet"
        if(rank==0)then
            call DMDACreate2d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
            &                 DMDA_STENCIL_BOX, jn, kn, PETSC_DECIDE, PETSC_DECIDE, &
            &                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, sliceDA, ierr)
            call DMSetUp(sliceDA, ierr)

            call DMDACreate3d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
            &                 DMDA_STENCIL_BOX, sink, jn, kn, PETSC_DECIDE, 1, 1,&
            &                 5, 0, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, disturbDA, ierr)
            call DMSetUp(disturbDA, ierr)

            call DMGetGlobalVector(sliceDA, disturb_slice, ierr)
            call PetscViewerBinaryOpen(PETSC_COMM_SELF, trim(turbfile), FILE_MODE_READ, viewer, ierr)
            call VecLoad(disturb_slice, viewer, ierr)
            call PetscViewerDestroy(viewer, ierr)

            call DMDAGetCorners(disturbDA,xs,ys,zs,xl,yl,zl,ierr)
            call DMGetGlobalVector(disturbDA, disturb_gather, ierr)
            call DMDAVecGetArrayF90(disturbDA, disturb_gather, disturbs, ierr)
            call DMDAVecGetArrayReadF90(sliceDA, disturb_slice, slice, ierr)
            do i=xs,xs+xl-1
                do j=ys,ys+yl-1
                    do k=zs,zs+zl-1
                        disturbs(:,i,j,k)=slice(:,j,k)
                    enddo
                enddo
            enddo
            call DMDAVecRestoreArrayReadF90(sliceDA, disturb_slice, slice, ierr)
            call DMDAVecRestoreArrayF90(disturbDA, disturb_gather, disturbs, ierr)

            call PetscViewerBinaryOpen(PETSC_COMM_SELF, trim(turbfiles),FILE_MODE_WRITE, viewer, ierr)
            call VecView(disturb_gather, viewer, ierr)
            call PetscViewerDestroy(viewer, ierr)

            call DMRestoreGlobalVector(disturbDA,disturb_gather,ierr)
            call DMRestoreGlobalVector(sliceDA,disturb_slice,ierr)
            call DMDestroy(sliceDA,ierr)
            call DMDestroy(disturbDA,ierr)
        endif
        call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
        &                 DMDA_STENCIL_BOX, sink, jn, kn, PETSC_DECIDE, 1, 1,&
        &                 5, 0, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, disturbDA, ierr)
        call DMSetUp(disturbDA, ierr)
        call DMGetGlobalVector(disturbDA, disturb_gather, ierr)
        call PetscViewerBinaryOpen(comm, trim(turbfiles),FILE_MODE_READ, viewer, ierr)
        call VecLoad(disturb_gather, viewer, ierr)
        call PetscViewerDestroy(viewer, ierr)

        call DMDAGetCorners(disturbDA,xs,ys,zs,xl,yl,zl,ierr)
        xe=xs+xl-1;ye=ys+yl-1;ze=zs+zl-1
        allocate(disturb(0:4, ys:ye, zs:ze))
        call DMDAVecGetArrayReadF90(disturbDA, disturb_gather, disturbs, ierr)
        do i=xs,xe
            do j=ys,ye
                do k=zs,ze
                    disturb(:,j,k) = disturbs(:,i,j,k)
                enddo
            enddo
        enddo
        call DMDAVecRestoreArrayReadF90(disturbDA, disturb_gather, disturbs, ierr)

        if(rank==0) call system("rm -f "//trim(turbfiles)//"*")
        call DMRestoreGlobalVector(disturbDA,disturb_gather,ierr)
        call DMDestroy(disturbDA,ierr)

    end subroutine set_disturb

    subroutine set_init_guess(comm) ! binary
        implicit none
        integer,intent(in) :: comm
        select case (init_guess_flg)
            case(.True.)
                call VecZeroEntries(turtle,ierr)
                call PetscViewerBinaryOpen(comm, trim(initfile),FILE_MODE_READ,viewer,ierr)
                call VecLoad(turtle, viewer, ierr)
                call PetscViewerDestroy(viewer, ierr)
            case(.False.)
                call VecZeroEntries(turtle,ierr)
        end select
    end subroutine set_init_guess

    subroutine check(comm)
        implicit none
        PetscInt,intent(in) :: comm
        if(rank==0)then
            write(*,"(3X,A,I5)") "Process Count :",sink
            write(*,*)
        endif
        call MPI_Barrier(comm,ierr)
        if(rank==(sink-1))then
            write(*,"(3X,A,I5,1X,A,I5)") "Rank:",rank,"jn    =",jn
            write(*,"(3X,A,I5,1X,A,F20.10)") "Rank:",rank,"Re    =",Re
            write(*,"(3X,A,I5,1X,A,2(F20.15))") "Rank:",rank,"Alpha =",Alpha
            write(*,"(3X,A,I5,1X,A,2(F20.15))") "Rank:",rank,"Omega =",Omega
        endif
        call MPI_Barrier(comm,ierr)
        if(rank==0)then
            write(*,"(3X,A,I5,1X,A,I5)") "Rank:",rank,"jn    =",jn
            write(*,"(3X,A,I5,1X,A,F20.10)") "Rank:",rank,"Re    =",Re
            write(*,"(3X,A,I5,1X,A,2(F20.15))") "Rank:",rank,"Alpha =",Alpha
            write(*,"(3X,A,I5,1X,A,2(F20.15))") "Rank:",rank,"Omega =",Omega
        endif
        call MPI_Barrier(comm,ierr)
        if(rank==0)then
            write(*,*)
            write(*,113) "data[0,0,0]：",qq(1,0,0,0),qq(2,0,0,0),qq(3,0,0,0),qq(4,0,0,0),qq(5,0,0,0)
            113 format (3X,A,5(F10.5))
            write(*,113) "data[1,0,0]：",qq(1,1,0,0),qq(2,1,0,0),qq(3,1,0,0),qq(4,1,0,0),qq(5,1,0,0)
            write(*,113) "data[2,0,0]：",qq(1,2,0,0),qq(2,2,0,0),qq(3,2,0,0),qq(4,2,0,0),qq(5,2,0,0)
            write(*,114) "mesh[0,0,0]：",xx(0,0,0),yy(0,0,0),zz(0,0,0)
            114 format (3X,A,3(F10.5))
            write(*,114) "mesh[1,0,0]：",xx(1,0,0),yy(1,0,0),zz(1,0,0)
            write(*,114) "mesh[2,0,0]：",xx(2,0,0),yy(2,0,0),zz(2,0,0)
        endif
        call MPI_Barrier(comm,ierr)
    end subroutine check

    subroutine ostream(comm)
        implicit none
        character(len=256) :: resultfile
        PetscInt, intent(in) :: comm

        call PetscPrintf(comm, "\n ----------------------------------\n", ierr)
        call PetscPrintf(comm, "               Ostream              \n\n", ierr)

        resultfile = "./data/turtle.pet"
        call PetscViewerBinaryOpen(comm,trim(resultfile),FILE_MODE_WRITE,viewer,ierr)
        call VecView(turtle, viewer, ierr)
        call PetscViewerDestroy(viewer, ierr)
        call PetscPrintf(comm, "   Binary Result: "//resultfile//"\n", ierr)
        resultfile = hdf5file
        if(io_type/='hdf5') call preload_hdf5(comm,resultfile)
        call PetscViewerHDF5Open(comm,trim(resultfile),FILE_MODE_UPDATE,viewer,ierr)
        call PetscObjectSetName(turtle,"hlns",ierr)
        call VecView(turtle,viewer,ierr)
        call PetscViewerHDF5WriteAttribute(viewer,"hlns","disturb.Alpha.r",PETSC_DOUBLE,real(Alpha),ierr)
        call PetscViewerHDF5WriteAttribute(viewer,"hlns","disturb.Alpha.i",PETSC_DOUBLE,aimag(Alpha),ierr)
        call PetscViewerHDF5WriteAttribute(viewer,"hlns","disturb.Beta.r",PETSC_DOUBLE,real(Beta),ierr)
        call PetscViewerHDF5WriteAttribute(viewer,"hlns","disturb.Beta.i",PETSC_DOUBLE,aimag(Beta),ierr)
        call PetscViewerHDF5WriteAttribute(viewer,"hlns","disturb.Omega.r",PETSC_DOUBLE,real(Omega),ierr)
        call PetscViewerHDF5WriteAttribute(viewer,"hlns","disturb.Omega.i",PETSC_DOUBLE,aimag(Omega),ierr)
        call PetscViewerDestroy(viewer, ierr)
        call PetscPrintf(comm, "   HDF5 Result: "//resultfile//"\n", ierr)

    end subroutine ostream

    subroutine preload_hdf5(comm,resultfile)
        implicit none
        character(len=256),intent(in) :: resultfile
        DM :: uni_coordDA,uni_meshDA
        integer,intent(in) :: comm
        Vec :: coord,flowfield
        PetscViewer :: Sviewer

        if(rank==0)then
            call DMDACreate3d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
            &                 DMDA_STENCIL_BOX, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,&
            &                 3, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, uni_coordDA, ierr)
            call DMSetUp(uni_coordDA, ierr)

            call DMDACreate3d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
            &                 DMDA_STENCIL_BOX, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,&
            &                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, uni_meshDA, ierr)
            call DMSetUp(uni_meshDA, ierr)
            call DMGetGlobalVector(uni_coordDA,coord,ierr)
            call DMGetGlobalVector(uni_meshDA,flowfield,ierr)

            call PetscViewerHDF5Open(PETSC_COMM_SELF,trim(resultfile),FILE_MODE_WRITE,viewer,ierr)

            call PetscViewerBinaryOpen(PETSC_COMM_SELF, trim(biflowfile),FILE_MODE_READ, Sviewer, ierr)
            call VecLoad(flowfield, Sviewer, ierr)
            call PetscViewerDestroy(Sviewer, ierr)
            call PetscObjectSetName(flowfield,"baseflow",ierr)
            call VecView(flowfield,viewer,ierr)

            call PetscViewerBinaryOpen(PETSC_COMM_SELF, trim(bigridfile),FILE_MODE_READ, Sviewer, ierr)
            call VecLoad(coord, Sviewer, ierr)
            call PetscViewerDestroy(Sviewer, ierr)
            call PetscObjectSetName(coord,"grid",ierr)
            call VecView(coord,viewer,ierr)
            call PetscViewerHDF5WriteAttribute(viewer,"grid","grid.In",PETSC_INT,in,ierr)
            call PetscViewerHDF5WriteAttribute(viewer,"grid","grid.Jn",PETSC_INT,jn,ierr)
            call PetscViewerHDF5WriteAttribute(viewer,"grid","grid.Kn",PETSC_INT,kn,ierr)
            call PetscViewerDestroy(viewer, ierr)

            call DMRestoreGlobalVector(uni_meshDA,flowfield,ierr)
            call DMRestoreGlobalVector(uni_coordDA,coord,ierr)
            call DMDestroy(uni_coordDA,ierr)
            call DMDestroy(uni_meshDA,ierr)
        endif

        call MPI_Barrier(comm,ierr)

    end subroutine preload_hdf5

end module mod_files
