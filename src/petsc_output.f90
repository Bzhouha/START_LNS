!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_petsc_output
    use petsc
    use global_parameters
    implicit none
    private
    public :: result_to_file
    PetscErrorCode :: ierr
    PetscViewer :: Viewer
    contains
    subroutine result_to_file(comm)
        implicit none
        PetscInt, INTENT(in) :: comm
        logical :: test 
        test = .True.
        call signal_ending(comm)
        call petsc_file(comm)
        if (rank==0) then
        call format_file(PETSC_COMM_SELF)
        if(test) call test_format_file(PETSC_COMM_SELF)
        endif
        call print_info(comm)
    end subroutine result_to_file

    subroutine petsc_file(comm)
        implicit none
        PetscInt,intent(in) :: comm
        call PetscViewerBinaryOpen(comm, "out/Turtle.petsc", FILE_MODE_WRITE, Viewer, ierr)
        call VecView(Turtle, Viewer, ierr)
        call PetscViewerDestroy(Viewer, ierr)
    end subroutine petsc_file

    subroutine format_file(comm)
        use penf,only:R_P
        implicit none
        complex(R_P),dimension(:,:,:,:),allocatable :: resrray
        real(R_P),dimension(:,:,:,:),allocatable :: modrray
        PetscScalar,dimension(:,:,:,:),pointer :: tmp 
        integer :: xs,ys,zs,xl,yl,zl,xe,ye,ze
        PetscInt,intent(in) :: comm
        integer :: i,j,k,l
        DM :: ResDA
        Vec :: Res
        call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
        &                 DMDA_STENCIL_STAR, in, jn, kn, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,&
        &                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ResDA, ierr)
        call DMSetUp(ResDA, ierr)
        call DMGetGlobalVector(ResDA, Res, ierr)
        call DMDAGetCorners(ResDA,xs,ys,zs,xl,yl,zl,ierr)
        xe=xs+xl-1;ye=ys+yl-1;ze=zs+zl-1
        call PetscViewerBinaryOpen(comm, "out/Turtle.petsc",FILE_MODE_READ, Viewer, ierr)
        call VecLoad(Res, Viewer, ierr)
        call PetscViewerDestroy(Viewer, ierr)
        allocate(resrray(0:4,xs:xe,ys:ye,zs:ze))
        allocate(modrray(0:4,xs:xe,ys:ye,zs:ze))
        call DMDAVecGetArrayReadF90(ResDA, Res, tmp, ierr)
        resrray=tmp
        call DMDAVecRestoreArrayReadF90(ResDA, Res, tmp, ierr)
        modrray=abs(resrray)
        open(37, file='out/turtle.csv',action='write',status='replace')
        write(37,*) "rho_r,rho_i,rho_m,u_r,u_i,u_m,v_r,v_i,v_m,w_r,w_i,w_m,T_r,T_i,T_m"
        do i=xs,xe 
            do j=ys,ye 
                do k=zs,ze
                    write(37,*)  real(resrray(0,i,j,k)),',',aimag(resrray(0,i,j,k)),',',modrray(0,i,j,k),',',&
                                &real(resrray(1,i,j,k)),',',aimag(resrray(1,i,j,k)),',',modrray(1,i,j,k),',',&
                                &real(resrray(2,i,j,k)),',',aimag(resrray(2,i,j,k)),',',modrray(2,i,j,k),',',&
                                &real(resrray(3,i,j,k)),',',aimag(resrray(3,i,j,k)),',',modrray(3,i,j,k),',',&
                                &real(resrray(4,i,j,k)),',',aimag(resrray(4,i,j,k)),',',modrray(4,i,j,k)
                enddo
            enddo
        enddo
        close(37)
        deallocate(resrray)
        deallocate(modrray)
    end subroutine format_file

    subroutine test_format_file(comm)
        use penf,only:R_P
        implicit none
        complex(R_P),dimension(:,:,:),allocatable :: resrray
        real(R_P),dimension(:,:,:),allocatable :: modrray
        PetscScalar,dimension(:,:,:),pointer :: tmp 
        integer :: xs,ys,zs,xl,yl,zl,xe,ye,ze
        PetscInt,intent(in) :: comm
        integer :: i,j,k,l
        DM :: ResDA
        Vec :: Res
        call DMDACreate2d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
        &                 DMDA_STENCIL_BOX, in, jn,  PETSC_DECIDE, PETSC_DECIDE, &
        &                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ResDA, ierr)
        call DMSetUp(ResDA, ierr)
        call DMGetGlobalVector(ResDA, Res, ierr)
        call DMDAGetCorners(ResDA,xs,ys,zs,xl,yl,zl,ierr)
        xe=xs+xl-1;ye=ys+yl-1;ze=zs+zl-1
        zs=0;ze=0
        call PetscViewerBinaryOpen(comm, "out/LPSE_istart=+1_iend=+301_phi.petsc",FILE_MODE_READ, Viewer, ierr)
        call VecLoad(Res, Viewer, ierr)
        call PetscViewerDestroy(Viewer, ierr)
        allocate(resrray(0:4,xs:xe,ys:ye))
        allocate(modrray(0:4,xs:xe,ys:ye))
        call DMDAVecGetArrayReadF90(ResDA, Res, tmp, ierr)
        resrray=tmp
        call DMDAVecRestoreArrayReadF90(ResDA, Res, tmp, ierr)
        modrray=abs(resrray)
        open(38, file='out/lpse.csv',action='write',status='replace')
        write(38,*) "rho_r,rho_i,rho_m,u_r,u_i,u_m,v_r,v_i,v_m,w_r,w_i,w_m,T_r,T_i,T_m"
        do i=xs,xe 
            do j=ys,ye 
                do k=zs,ze
                    write(38,*)  real(resrray(0,i,j)),',',aimag(resrray(0,i,j)),',',modrray(0,i,j),',',&
                                &real(resrray(1,i,j)),',',aimag(resrray(1,i,j)),',',modrray(1,i,j),',',&
                                &real(resrray(2,i,j)),',',aimag(resrray(2,i,j)),',',modrray(2,i,j),',',&
                                &real(resrray(3,i,j)),',',aimag(resrray(3,i,j)),',',modrray(3,i,j),',',&
                                &real(resrray(4,i,j)),',',aimag(resrray(4,i,j)),',',modrray(4,i,j)
                enddo
            enddo
        enddo
        close(38)
        deallocate(resrray)
        deallocate(modrray)
    end subroutine test_format_file
    
    subroutine signal_ending(comm)
        implicit none 
        PetscInt,intent(in) :: comm 
        call PetscPrintf(comm, "\n", ierr)
        call PetscPrintf(comm, " ===========================================================================\n", ierr)
        call PetscPrintf(comm, " =                                 输    出                                = \n", ierr)
        call PetscPrintf(comm, " ===========================================================================\n", ierr)
        call PetscPrintf(comm, " ----------------------------------\n", ierr)
        call PetscPrintf(comm, "              输出结果               \n", ierr)
        call PetscPrintf(comm, " ----------------------------------\n", ierr)
    end subroutine signal_ending

    subroutine print_info(comm)
        implicit none
        PetscInt,INTENT(in) :: comm
        call PetscPrintf(comm," \n", ierr)
        call PetscPrintf(comm," 输出解向量...                                     ooo    ooo\n", ierr)
        call PetscPrintf(comm,"   输出结束。                                     o   o  o   o\n", ierr)
        call PetscPrintf(comm,"                                            ooo   o   o  o   o   ooo\n", ierr)
        call PetscPrintf(comm," ========================================  o   o   ooo    ooo   o   o\n", ierr)
        call PetscPrintf(comm,"                                           o   o                o   o\n",ierr)
        call PetscPrintf(comm,"                                            ooo    oooooooooo    ooo\n",ierr)
        call PetscPrintf(comm,"               ooo    ooo                        o            o\n",ierr)
        call PetscPrintf(comm,"              o   o  o   o                      o              o\n",ierr)
        call PetscPrintf(comm,"        ooo   o   o  o   o   ooo                 o            o\n",ierr)
        call PetscPrintf(comm,"       o   o   ooo    ooo   o   o                  oooooooooo\n",ierr)
        call PetscPrintf(comm,"       o   o                o   o\n",ierr)
        call PetscPrintf(comm,"        ooo    oooooooooo    ooo\n",ierr)
        call PetscPrintf(comm,"             o            o\n",ierr)
        call PetscPrintf(comm,"            o              o\n",ierr)
        call PetscPrintf(comm,"             o            o\n",ierr)
        call PetscPrintf(comm,"               oooooooooo\n",ierr)
    end subroutine print_info
end module mod_petsc_output  
