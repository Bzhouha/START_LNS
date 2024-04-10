!------------------------------------------------------------------------------
!
! Copyright (c) 2019-2024 Bzhouha 
! All rights reserved.
! 
! This file is part of START_LNS.
! START_LNS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! START_LNS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>. 
!
!------------------------------------------------------------------------------
!
! This work used the PETSc library, which is developed by the PETSc Development Team. See https://petsc.org/ for more information.
! Copyright (c) 2023, PETSc Development Team
! All rights reserved.
! This file is subject to the terms and conditions of the BSD 2-Clause License. See the file LICENSE in the top-level directory for more information.
!
!------------------------------------------------------------------------------
!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_files 
!
! 这个模块处理文件流，包括:
!
!                 1).从config配置文件读入参数；
!
!                 2).支持Plot3d格式、binary格式、hdf5格式文件的输入输出
!
!   call istream(comm)  读入配置文件，初始化PETSc以及从文件载入流场、网格、边界、初值等。
!
!   call ostream(comm) 将计算结果以Binary、HDF5格式输出。
!
    use mod_parameters
    use petsc
    public :: istream
    public :: ostream
    private
    PetscErrorCode  :: ierr
    PetscViewer :: viewer
    DM :: coordDA
    contains
    subroutine istream(comm)
        implicit none
        PetscInt, intent(in) :: comm

        call config(comm)

        call set_DM(comm)

        call load(comm)

    end subroutine istream

    subroutine config(comm)
        use mod_cfgio_adapter
        implicit none
        PetscBool :: ksp_flg,snes_flg,newt_flg,newtsub_flg,mg_flg
        character(len=256) :: cfg_file
        PetscInt,intent(in) :: comm
        PetscBool :: set

        time0 = MPI_Wtime()

        call mpi_comm_rank(comm,rank,ierr)
        call mpi_comm_size(comm,sink,ierr)
        ! call PetscPrintf(comm, "\n"//char(27)//"[0m"//"           "//char(27)//"[0;1;4;37;44m"// &
        ! &    "S T A R T - L N S"// char(27)//"[0m"//"\n\n\n", ierr)
        call PetscPrintf(comm, "\n\n"//char(27)//"[0;1;31m"//" > "//char(27)//"[0m"//" "//char(27)//"[0;1;4;37;44m"// &
        &    "S T A R T - L N S"// char(27)//"[0m"//"\n\n\n", ierr)

        call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-f',cfg_file,set,ierr)
        if(.not. set) then
            write(*,*) 'should use -f option to determin the config file.'
            stop
        endif

        call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-asf',ksp_flg,ierr)
        call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-snes',snes_flg,ierr)
        call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nasf',newt_flg,ierr)
        call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-lnasf',newtsub_flg,ierr)
        if(ksp_flg)then
            solver_mode='asf'; split_mode=0
        endif
        if(snes_flg)then
            solver_mode='snes';split_mode=1
        endif
        if(newt_flg)then
            solver_mode='nasf';split_mode=0
        endif
        if(newtsub_flg)then
            solver_mode='lnasf';split_mode=0
        endif

        call PetscOptionsGetInt (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-lv',levels,mg_flg,ierr)
        if(.not.mg_flg) levels=0

        call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-plt',set,ierr)
        if(set) io_type="plt"
        call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-binary',set,ierr)
        if(set) io_type="binary"
        call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-hdf5',set,ierr)
        if(set) io_type="hdf5"

        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-cfl',cfl,set,ierr)
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-dt',dt,set,ierr)
        call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-dt',usedt,ierr)
        if(usedt)then
            if(abs(dt-TAG)<1e-5)then
                calculate_dt=.True.
            else
                calculate_dt=.False.
            endif
        else
            calculate_dt=.False.
        endif

        call PetscOptionsGetInt (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nx',nx,set,ierr)
        call PetscOptionsGetInt (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-ny',ny,set,ierr)
        call PetscOptionsGetInt (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nz',nz,set,ierr)

        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-lm',lm,set,ierr)

        if(ksp_flg)  call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_monitor",PETSC_NULL_CHARACTER,ierr)

        call cfg_loader(trim(cfg_file))
        call MPI_Barrier(comm,ierr)

        call bcast(comm)

    end subroutine config

    subroutine bcast(comm)
        implicit none
        integer(KIND=MPI_ADDRESS_KIND) :: adres_initguess,adres_inlet,adres_h5fn
        integer(KIND=MPI_ADDRESS_KIND) :: adres_mode,adres_Ma,adres_Re,adres_Te
        integer(KIND=MPI_ADDRESS_KIND) :: adres_in,adres_jn,adres_kn,adres_ln
        integer(KIND=MPI_ADDRESS_KIND) :: adres_Alpha,adres_Omega,adres_Beta
        integer(KIND=MPI_ADDRESS_KIND) :: displacement(14)
        integer :: block_lengths(14)
        PetscInt,intent(in) :: comm
        integer :: pack_type
        integer :: types(14)

        call MPI_Get_address(in,adres_in,ierr)
        call MPI_Get_address(jn,adres_jn,ierr)
        call MPI_Get_address(kn,adres_kn,ierr)
        call MPI_Get_address(ln,adres_ln,ierr)
        call MPI_Get_address(lns_mode,adres_mode,ierr)
        call MPI_Get_address(ex_ini_gus_flg,adres_initguess,ierr)
        call MPI_Get_address(inlet_file_flg,adres_inlet,ierr)
        call MPI_Get_address(Ma,adres_Ma,ierr)
        call MPI_Get_address(Re,adres_Re,ierr)
        call MPI_Get_address(Te,adres_Te,ierr)
        call MPI_Get_address(Alpha,adres_Alpha,ierr)
        call MPI_Get_address(Beta,adres_Beta,ierr)
        call MPI_Get_address(Omega,adres_Omega,ierr)
        call MPI_Get_address(hdf5file,adres_h5fn,ierr)
        displacement(1)  = 0
        displacement(2)  = adres_jn-adres_in
        displacement(3)  = adres_kn-adres_in
        displacement(4)  = adres_ln-adres_in
        displacement(5)  = adres_mode-adres_in
        displacement(6)  = adres_initguess-adres_in
        displacement(7)  = adres_inlet-adres_in
        displacement(8)  = adres_Ma-adres_in
        displacement(9)  = adres_Re-adres_in
        displacement(10) = adres_Te-adres_in
        displacement(11) = adres_Alpha-adres_in
        displacement(12) = adres_Beta-adres_in
        displacement(13) = adres_Omega-adres_in
        displacement(14) = adres_h5fn-adres_in
        block_lengths=1
        block_lengths(14)=256
        types=(/MPI_INTEGER4,MPI_INTEGER4,MPI_INTEGER4,MPI_INTEGER4,MPI_INTEGER4,&
                MPI_LOGICAL,MPI_LOGICAL,MPI_REAL8,MPI_REAL8,MPI_REAL8,&
                MPI_COMPLEX16,MPI_COMPLEX16,MPI_COMPLEX16,MPI_CHAR/)
        call MPI_Type_create_struct(14,block_lengths,displacement,types,pack_type,ierr)
        call MPI_Type_commit(pack_type,ierr)
        call MPI_Bcast(in,1,pack_type,0,comm,ierr)
        call MPI_Barrier(comm,ierr)
        call MPI_Type_free(pack_type,ierr)
    end subroutine bcast

    subroutine set_DM(comm)
        implicit none
        PetscInt, intent(in) :: comm

        call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
        &                 DMDA_STENCIL_BOX, in, jn, kn, nx, ny, nz, &
        &                 1, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, DA, ierr)
        call DMSetFromOptions(DA,ierr)
        call DMSetUp(DA, ierr)

        call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
        &                 DMDA_STENCIL_BOX, in, jn, kn, nx, ny, nz,&
        &                 3, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, coordDA, ierr)
        call DMSetFromOptions(coordDA,ierr)
        call DMSetUp(coordDA, ierr)

        call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC,&
        &                 DMDA_STENCIL_BOX, in, jn, kn, nx, ny, nz,&
        &                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, meshDA, ierr)
        call DMSetMatType(meshDA,MATBAIJ,ierr)
        call DMSetFromOptions(meshDA,ierr)
        call DMSetUp(meshDA, ierr)

        call DMDAGetInfo(DA, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
        &   nx, ny, nz, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
        &   PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

        call DMGetGlobalVector(meshDA, turtle, ierr)
        call VecZeroEntries(turtle,ierr)
        call DMDAGetGhostCorners(DA,igs,jgs,kgs,igl,jgl,kgl,ierr)
        call DMDAGetCorners(DA,is,js,ks,il,jl,kl,ierr)
        ige=igs+igl-1; jge=jgs+jgl-1; kge=kgs+kgl-1
        ie=is+il-1; je=js+jl-1; ke=ks+kl-1
    end subroutine set_DM

    subroutine load(comm)
        implicit none
        integer,intent(in) :: comm

        call PetscPrintf(comm, char(27)//"[0;36m"//" IStream            \n"//char(27)//"[0m", ierr)
        call PetscPrintf(comm, char(27)//"[0;1;36m"//" -----------------------------------\n"//char(27)//"[0m", ierr)
        select case(io_type)
            case("plt")
                call PetscPrintf(comm, "\n   I/O type      -> "//char(27)//"[0;33m"//"Plot3d"//char(27)//"[0m"//"\n\n", ierr)
                call load_raw_files(comm)
            case("binary")
                call PetscPrintf(comm, "\n   I/O type      -> "//char(27)//"[0;33m"//"Binary"//char(27)//"[0m"//"\n\n", ierr)
                call load_binary_files(comm)
            case("hdf5")
                call PetscPrintf(comm, "\n   I/O type      -> "//char(27)//"[0;33m"//"HDF5"//char(27)//"[0m"//"\n\n", ierr)
                call load_hdf5_files(comm)
        end select

        call set_init_guess(comm)
        call load_inlet(comm)
        if(wall_bc==797)then 
            call load_wallbc_from_ini_gus(comm)
        endif
        call check(comm)
        call DMDestroy(coordDA,ierr)
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
        PetscScalar, pointer :: inlets(:,:,:,:)
        PetscScalar, pointer :: grid(:,:,:,:)
        PetscScalar, pointer :: flow(:,:,:,:)
        PetscScalar, pointer :: slice(:,:,:)
        integer :: xs,ys,zs,xl,yl,zl
        DM :: unicordDA,unimeshDA
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
        ! read(12) ((((qq(l,i,j,k), i=1,in), j=1,jn), k=1,kn), l=1,5)
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
        &                 3, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, unicordDA, ierr)
        call DMSetUp(unicordDA, ierr)

        call DMDACreate3d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
        &                 DMDA_STENCIL_BOX, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,&
        &                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, unimeshDA, ierr)
        call DMSetUp(unimeshDA, ierr)

        call DMGetGlobalVector(unicordDA,coord,ierr)
        call DMDAGetCorners(unicordDA,xs,ys,zs,xl,yl,zl,ierr)
        call DMDAVecGetArrayF90(unicordDA,coord,grid,ierr)
        do i=xs,xs+xl-1
            do j=ys,ys+yl-1
                do k=zs,zs+zl-1
                    grid(0, i, j, k) = xx(i+1, j+1, k+1)
                    grid(1, i, j, k) = yy(i+1, j+1, k+1)
                    grid(2, i, j, k) = zz(i+1, j+1, k+1)
                enddo
            enddo
        enddo
        call DMDAVecRestoreArrayF90(unicordDA,coord,grid,ierr)

        call DMGetGlobalVector(unimeshDA,flowfield,ierr)
        call DMDAGetCorners(unimeshDA,xs,ys,zs,xl,yl,zl,ierr)
        call DMDAVecGetArrayF90(unimeshDA,flowfield,flow,ierr)
        do i=xs,xs+xl-1
            do j=ys,ys+yl-1
                do k=zs,zs+zl-1
                    do l=0, 4
                        flow(l,i,j,k) = qq(l+1, i+1, j+1, k+1)
                    enddo
                enddo
            enddo
        enddo
        call DMDAVecRestoreArrayF90(unimeshDA,flowfield,flow,ierr)

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
        call DMRestoreGlobalVector(unimeshDA,flowfield,ierr)
        call DMRestoreGlobalVector(unicordDA,coord,ierr)
        call DMDestroy(unicordDA,ierr)
        call DMDestroy(unimeshDA,ierr)

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
        PetscScalar, pointer :: tmp(:,:,:,:)
        Vec :: flowfield_local,coord_local
        PetscInt, intent(in) :: comm
        Vec :: coord,flowfield
        integer :: l,i,j,k
        PetscBool :: has

        call PetscPrintf(comm,"   Input File    -> "//char(27)//"[0;33m"//trim(hdf5file)//char(27)//"[0m"//"\n\n",ierr)

        call DMGetGlobalVector(meshDA, flowfield, ierr)
        call PetscObjectSetName(flowfield,"flow",ierr)
        call DMGetGlobalVector(coordDA, coord, ierr)
        call PetscObjectSetName(coord,"grid",ierr)
        call PetscObjectSetName(turtle,"shapefunc",ierr)

        call PetscViewerHDF5Open(comm,trim(hdf5file),FILE_MODE_READ,viewer,ierr)
        call VecLoad(flowfield,viewer,ierr)
        call VecLoad(coord,viewer,ierr)
        call PetscViewerHDF5HasDataset(viewer,'shapefunc',has,ierr)
        if(has) call VecLoad(turtle,viewer,ierr)
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

    subroutine set_init_guess(comm) ! binary
        implicit none
        integer,intent(in) :: comm
        if(ex_ini_gus_flg)then
            ini_gus_flg=.True.
            call VecZeroEntries(turtle,ierr)
            call PetscViewerBinaryOpen(comm, trim(initfile),FILE_MODE_READ,viewer,ierr)
            call VecLoad(turtle, viewer, ierr)
            call PetscViewerDestroy(viewer, ierr)
        endif
    end subroutine set_init_guess

    subroutine load_inlet(comm)
        implicit none
        integer,intent(in) :: comm
        select case(inlet_file_flg)
        case(.True.)
            call load_inlet_file(comm)
        case(.False.)
            call load_inlet_from_ini_gus(comm)
        end select
    end subroutine load_inlet

    subroutine load_inlet_from_ini_gus(comm)
        implicit none
        PetscScalar,pointer :: xr(:,:,:,:)
        integer,intent(in) :: comm

        allocate(inlet(0:4, js:je, ks:ke))
        inlet = 0.0d0
        call DMDAVecGetArrayReadF90(meshDA, turtle, xr, ierr)
        if(is==0)then
            inlet(:,:,:)=xr(:,is,:,:)
        endif
        call DMDAVecRestoreArrayReadF90(meshDA, turtle, xr, ierr)
        call MPI_Barrier(comm,ierr)
    end subroutine load_inlet_from_ini_gus

    subroutine load_wallbc_from_ini_gus(comm)
        implicit none
        PetscScalar,pointer :: xr(:,:,:,:)
        integer,intent(in) :: comm

        allocate(wall(0:4,is:ie,ks:ke))
        wall = 0.0d0 
        call DMDAVecGetArrayReadF90(meshDA, turtle, xr, ierr)
        if(js==0)then 
            wall(:,:,:) = xr(:,:,js,:)
        endif
        call DMDAVecRestoreArrayReadF90(meshDA, turtle, xr, ierr)
        call MPI_Barrier(comm,ierr)
    end subroutine load_wallbc_from_ini_gus

    subroutine load_inlet_file(comm)
        implicit none
        PetscScalar, pointer :: inlets(:,:,:,:)
        integer :: xs,ys,zs,xl,yl,zl,xe,ye,ze
        PetscScalar, pointer :: slice(:,:,:)
        character(len=256) :: inletfiles ! 文件名：边界文件
        integer,intent(in) :: comm
        Vec :: inlet_gather
        Vec :: inlet_slice
        integer :: i,j,k
        DM :: inletDA
        DM :: sliceDA

        inletfiles = "./data/inlets.pet"
        if(rank==0)then
            call DMDACreate2d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
            &                 DMDA_STENCIL_BOX, jn, kn, PETSC_DECIDE, PETSC_DECIDE, &
            &                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, sliceDA, ierr)
            call DMSetUp(sliceDA, ierr)

            call DMDACreate3d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
            &                 DMDA_STENCIL_BOX, sink, jn, kn, PETSC_DECIDE, 1, 1,&
            &                 5, 0, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, inletDA, ierr)
            call DMSetUp(inletDA, ierr)

            call DMGetGlobalVector(sliceDA, inlet_slice, ierr)
            call PetscViewerBinaryOpen(PETSC_COMM_SELF, trim(inletfile), FILE_MODE_READ, viewer, ierr)
            call VecLoad(inlet_slice, viewer, ierr)
            call PetscViewerDestroy(viewer, ierr)

            call DMDAGetCorners(inletDA,xs,ys,zs,xl,yl,zl,ierr)
            call DMGetGlobalVector(inletDA, inlet_gather, ierr)
            call DMDAVecGetArrayF90(inletDA, inlet_gather, inlets, ierr)
            call DMDAVecGetArrayReadF90(sliceDA, inlet_slice, slice, ierr)
            do i=xs,xs+xl-1
                do j=ys,ys+yl-1
                    do k=zs,zs+zl-1
                        inlets(:,i,j,k)=slice(:,j,k)
                    enddo
                enddo
            enddo
            call DMDAVecRestoreArrayReadF90(sliceDA, inlet_slice, slice, ierr)
            call DMDAVecRestoreArrayF90(inletDA, inlet_gather, inlets, ierr)

            call PetscViewerBinaryOpen(PETSC_COMM_SELF, trim(inletfiles),FILE_MODE_WRITE, viewer, ierr)
            call VecView(inlet_gather, viewer, ierr)
            call PetscViewerDestroy(viewer, ierr)

            call DMRestoreGlobalVector(inletDA,inlet_gather,ierr)
            call DMRestoreGlobalVector(sliceDA,inlet_slice,ierr)
            call DMDestroy(sliceDA,ierr)
            call DMDestroy(inletDA,ierr)
        endif
        call MPI_Barrier(comm,ierr)
        call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
        &                 DMDA_STENCIL_BOX, sink, jn, kn, PETSC_DECIDE, 1, 1,&
        &                 5, 0, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, inletDA, ierr)
        call DMSetUp(inletDA, ierr)
        call DMGetGlobalVector(inletDA, inlet_gather, ierr)
        call PetscViewerBinaryOpen(comm, trim(inletfiles),FILE_MODE_READ, viewer, ierr)
        call VecLoad(inlet_gather, viewer, ierr)
        call PetscViewerDestroy(viewer, ierr)

        call DMDAGetCorners(inletDA,xs,ys,zs,xl,yl,zl,ierr)
        xe=xs+xl-1;ye=ys+yl-1;ze=zs+zl-1
        allocate(inlet(0:4, ys:ye, zs:ze))
        inlet = 0.0d0
        call DMDAVecGetArrayReadF90(inletDA, inlet_gather, inlets, ierr)
        do i=xs,xe
            do j=ys,ye
                do k=zs,ze
                    inlet(:,j,k) = inlets(:,i,j,k)
                enddo
            enddo
        enddo
        call DMDAVecRestoreArrayReadF90(inletDA, inlet_gather, inlets, ierr)

        if(rank==0) call system("rm -f "//trim(inletfiles)//"*")
        call DMRestoreGlobalVector(inletDA,inlet_gather,ierr)
        call DMDestroy(inletDA,ierr)
        call MPI_Barrier(comm,ierr)

    end subroutine load_inlet_file

    subroutine check(comm)
        implicit none
        PetscInt,intent(in) :: comm
        character(len=45) :: str_45
        character(len=34) :: str_34
        character(len=30) :: str_30
        character(len=27) :: str_27
        character(len=18) :: str_18
        character(len=5) :: str_5
        character(len=1) :: str_1

        write(str_5,"(I5)") sink
        call PetscPrintf(comm,"   Process Count -> "//char(27)//"[0;33m"//str_5//char(27)//"[0m"//"\n\n",ierr)

        call PetscPrintf(comm,"   Solver Type   -> "//char(27)//"[0;33m"//solver_mode//char(27)//"[0m"//"\n\n",ierr)

        write(str_1,"(I1)") lns_mode
        call PetscPrintf(comm,"   LNS Dimension -> "//char(27)//"[0;33m"//str_1//"D-LNS"//char(27)//"[0m"//"\n\n",ierr)

        write(str_18,"(3I6)") in,jn,kn
        call PetscPrintf(comm,"   Grid Size     -> "//char(27)//"[0;33m"//str_18//char(27)//"[0m"//"\n\n",ierr)

        write(str_18,"(3I6)") nx,ny,nz
        call PetscPrintf(comm,"   Partition     -> "//char(27)//"[0;33m"//str_18//char(27)//"[0m"//"\n\n",ierr)

        7888 format (F17.9)
        7777 format (2F17.13) 
        7666 format (3X,'Grid[',2(I5,','),I5,' ] ->')
        7555 format (3(F9.3))
        7444 format (3X,'Flow[',2(I5,','),I5,' ] ->')
        7333 format (5(F9.3))

        write(str_18,7888) Re 
        call PetscPrintf(comm,"   Re    -> "//char(27)//"[0;33m"//str_18//char(27)//"[0m"//"\n",ierr)
        write(str_18,7888) Ma 
        call PetscPrintf(comm,"   Ma    -> "//char(27)//"[0;33m"//str_18//char(27)//"[0m"//"\n",ierr)
        write(str_18,7888) Te
        call PetscPrintf(comm,"   Te    -> "//char(27)//"[0;33m"//str_18//char(27)//"[0m"//"\n",ierr)

        write(str_34,7777) Alpha
        call PetscPrintf(comm,"   Alpha -> "//char(27)//"[0;33m"//str_34//char(27)//"[0m"//"\n",ierr)
        write(str_34,7777) Beta
        call PetscPrintf(comm,"   Beta  -> "//char(27)//"[0;33m"//str_34//char(27)//"[0m"//"\n",ierr)
        write(str_34,7777) Omega
        call PetscPrintf(comm,"   Omega -> "//char(27)//"[0;33m"//str_34//char(27)//"[0m"//"\n\n",ierr)

        write(str_30,7666) igs,jgs,kgs; write(str_27,7555) xx(igs,jgs,kgs),yy(igs,jgs,kgs),zz(igs,jgs,kgs)
        call PetscPrintf(comm,str_30//char(27)//"[0;33m"//str_27//char(27)//"[0m"//"\n",ierr)
        write(str_30,7666) ige,jge,kge; write(str_27,7555) xx(ige,jge,kge),yy(ige,jge,kge),zz(ige,jge,kge)
        call PetscPrintf(comm,str_30//char(27)//"[0;33m"//str_27//char(27)//"[0m"//"\n",ierr)

        write(str_30,7444) igs,jgs,kgs
        write(str_45,7333) qq(1,igs,jgs,kgs),qq(2,igs,jgs,kgs),    &
        &                  qq(3,igs,jgs,kgs),qq(4,igs,jgs,kgs),qq(5,igs,jgs,kgs)
        call PetscPrintf(comm,str_30//char(27)//"[0;33m"//str_45//char(27)//"[0m"//"\n",ierr)

        write(str_30,7444) ige,jge,kge
        write(str_45,7333) qq(1,ige,jge,kge),qq(2,ige,jge,kge),    &
        &                  qq(3,ige,jge,kge),qq(4,ige,jge,kge),qq(5,ige,jge,kge)
        call PetscPrintf(comm,str_30//char(27)//"[0;33m"//str_45//char(27)//"[0m"//"\n",ierr)

        call PetscPrintf(comm,"\n\n",ierr)

        call MPI_Barrier(comm,ierr)

    end subroutine check

    subroutine ostream(comm)
        implicit none
        complex(R_P),parameter :: Ci = cmplx(0.0d0,1.0d0,R_P)
        PetscScalar,pointer :: tmp(:,:,:,:),grid(:,:,:,:)
        real(R_P),allocatable,dimension(:,:) :: x0
        character(len=96) :: hdf5out,pltout
        character(len=20) :: str_time
        PetscInt, intent(in) :: comm
        Vec :: VecT,coord,flowfield
        DM :: unimeshDA,unicordDA
        PetscViewer :: sviewer
        integer :: i,j,k,l
        Vec :: VecX

        call PetscPrintf(comm, "\n\n", ierr)
        call PetscPrintf(comm, char(27)//"[0;36m"//" OStream            \n"//char(27)//"[0m", ierr)
        call PetscPrintf(comm, char(27)//"[0;1;36m"//" -----------------------------------\n"//char(27)//"[0m", ierr)
        call PetscPrintf(comm, "\n", ierr)

        biresfile = "./data/shapefunc.pet"
        call PetscViewerBinaryOpen(comm,trim(biresfile),FILE_MODE_WRITE,viewer,ierr)
        call VecView(turtle, viewer, ierr)
        call PetscViewerDestroy(viewer, ierr)
        ! call PetscPrintf(comm, "   Binary Result -> "//trim(biresfile)//"\n", ierr)
        call MPI_Barrier(comm,ierr)
        call DMDestroy(meshDA,ierr)
        call MPI_Barrier(comm,ierr)

        if(rank==0)then
            hdf5out = trim(output_prefix)//".h5"
            ! 生成单进程DM
            call DMDACreate3d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
            &                 DMDA_STENCIL_BOX, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,&
            &                 5, 0, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, unimeshDA, ierr)
            call DMSetUp(unimeshDA, ierr)
            call DMDACreate3d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
            &                 DMDA_STENCIL_BOX, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,&
            &                 3, 0, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, unicordDA, ierr)
            call DMSetUp(unicordDA, ierr)
            ! 单进程读入大乌龟
            call DMGetGlobalVector(unimeshDA,VecT,ierr)
            call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(biresfile),FILE_MODE_READ,viewer,ierr)
            call VecLoad(VecT, viewer, ierr)
            call PetscViewerDestroy(viewer, ierr)

            call DMGetGlobalVector(unicordDA, VecX, ierr)
            call PetscObjectSetName(VecX,"grid",ierr)
            call PetscViewerHDF5Open(PETSC_COMM_SELF,trim(hdf5file),FILE_MODE_READ,viewer,ierr)
            call VecLoad(VecX,viewer,ierr)
            call PetscViewerDestroy(viewer, ierr)

            allocate(xx(0:in-1,0:jn-1,0:kn-1))
            call DMDAVecGetArrayReadF90(unicordDA, VecX, grid, ierr)
            xx = real(grid(0, :, :, :))
            call DMDAVecRestoreArrayReadF90(unicordDA, VecX, grid, ierr)

            ! 后处理
            allocate(x0(0:jn-1,0:kn-1))
            x0(:,:) = xx(0,:,:)
            call DMDAVecGetArrayF90(unimeshDA, VecT, tmp, ierr)
            do k=0,kn-1
                do j=0,jn-1
                    do i=0,in-1
                        tmp(:,i,j,k) = tmp(:,i,j,k)*exp(Ci*Alpha*(xx(i,j,k)-x0(j,k)))
                    enddo 
                enddo
            enddo
            call DMDAVecRestoreArrayF90(unimeshDA, VecT, tmp, ierr)

            ! 把大乌龟写入Plot3D格式文件

            ! pltout = trim(output_prefix)//".p3d"
            ! call DMDAVecGetArrayReadF90(unimeshDA, VecT, tmp, ierr)
            ! open(18, file=trim(pltout),action='write',form='unformatted')
            ! write(18)
            ! select case (lns_mode)
            ! case(2)
            !     write(18) in,jn,ln
            ! case(3)
            !     write(18) in,jn,kn,ln
            ! end select
            ! write(18) ((((tmp(l,i,j,k), i=0,in-1), j=0,jn-1), k=0,kn-1), l=0,4)
            ! close(18)
            ! call DMDAVecRestoreArrayReadF90(unimeshDA, VecT, tmp, ierr)
            ! call PetscPrintf(comm, "   Plot3D Result -> "//trim(pltout)//"\n", ierr)

            ! 把大乌龟写入HDF5格式文件

            ! call DMGetGlobalVector(unicordDA,coord,ierr)
            ! call DMGetGlobalVector(unimeshDA,flowfield,ierr)
            call PetscViewerHDF5Open(PETSC_COMM_SELF,trim(hdf5out),FILE_MODE_WRITE,viewer,ierr)
            ! call PetscViewerBinaryOpen(PETSC_COMM_SELF, trim(biflowfile),FILE_MODE_READ, sviewer, ierr)
            ! call VecLoad(flowfield, sviewer, ierr)
            ! call PetscViewerDestroy(sviewer, ierr)
            ! call PetscObjectSetName(flowfield,"flow",ierr)
            ! call VecView(flowfield,viewer,ierr)
            ! call DMRestoreGlobalVector(unimeshDA,flowfield,ierr)

            ! call PetscViewerBinaryOpen(PETSC_COMM_SELF, trim(bigridfile),FILE_MODE_READ, sviewer, ierr)
            ! call VecLoad(coord, sviewer, ierr)
            ! call PetscViewerDestroy(sviewer, ierr)
            ! call PetscObjectSetName(coord,"grid",ierr)
            ! call VecView(coord,viewer,ierr)
            ! call PetscViewerHDF5WriteAttribute(viewer,"grid","In",PETSC_INT,in,ierr)
            ! call PetscViewerHDF5WriteAttribute(viewer,"grid","Jn",PETSC_INT,jn,ierr)
            ! call PetscViewerHDF5WriteAttribute(viewer,"grid","Kn",PETSC_INT,kn,ierr)
            ! call DMRestoreGlobalVector(unicordDA,coord,ierr)
            ! call DMDestroy(unicordDA,ierr)

            call PetscObjectSetName(VecT,"shapefunc",ierr)
            call VecView(VecT,viewer,ierr)
            call PetscViewerHDF5WriteAttribute(viewer,"shapefunc","disturb.Alpha.r",PETSC_DOUBLE,real(Alpha),ierr)
            call PetscViewerHDF5WriteAttribute(viewer,"shapefunc","disturb.Alpha.i",PETSC_DOUBLE,aimag(Alpha),ierr)
            call PetscViewerHDF5WriteAttribute(viewer,"shapefunc","disturb.Beta.r",PETSC_DOUBLE,real(Beta),ierr)
            call PetscViewerHDF5WriteAttribute(viewer,"shapefunc","disturb.Beta.i",PETSC_DOUBLE,aimag(Beta),ierr)
            call PetscViewerHDF5WriteAttribute(viewer,"shapefunc","disturb.Omega.r",PETSC_DOUBLE,real(Omega),ierr)
            call PetscViewerHDF5WriteAttribute(viewer,"shapefunc","disturb.Omega.i",PETSC_DOUBLE,aimag(Omega),ierr)
            call DMRestoreGlobalVector(unimeshDA,VecT,ierr)
            call DMDestroy(unimeshDA,ierr)
            call PetscViewerDestroy(viewer, ierr)
            call PetscPrintf(comm, "   HDF5 Result -> "//char(27)//"[0;33m"//trim(hdf5out)//char(27)//"[0m"//"\n", ierr)
            call PetscPrintf(comm, "\n", ierr)

            call PetscViewerDestroy(viewer, ierr)
            call system("rm "//trim(biresfile)//"*")
        endif
        call MPI_Barrier(comm,ierr)
        time = MPI_Wtime()-time0
        write(str_time,"(F20.1)") time
        call PetscPrintf(comm, char(27)//"[0;1m"//"   Total time"//char(27)//"[0m"//"  -> "//char(27)//"[0;32m"//trim(str_time)//"s"//char(27)//"[0m"//"\n", ierr)
        call PetscPrintf(comm, "\n\n", ierr)
        call PetscBarrier(PETSC_NULL_VEC, ierr)

    end subroutine ostream

end module mod_files
