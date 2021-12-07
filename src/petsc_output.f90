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
        integer :: xs, ys, zs, xl, yl, zl
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
        call PetscViewerBinaryOpen(comm, "out/Turtle.petsc",FILE_MODE_READ, Viewer, ierr)
        call VecLoad(Res, Viewer, ierr)
        call PetscViewerDestroy(Viewer, ierr)
        allocate(resrray(5,in,jn,kn))
        allocate(modrray(5,in,jn,kn))
        call DMDAVecGetArrayReadF90(ResDA, Res, tmp, ierr)
        resrray=tmp
        call DMDAVecRestoreArrayReadF90(ResDA, Res, tmp, ierr)
        modrray=abs(resrray)
        open(37, file='out/turtle.csv',action='write',status='replace')
        write(37,*) "rho_r,rho_i,rho_m,u_r,u_i,u_m,v_r,v_i,v_m,w_r,w_i,w_m,T_r,T_i,T_m"
        do i=1,in 
            do j=1,jn 
                do k=1,kn
                    write(37,*)  real(resrray(1,i,j,k)),',',aimag(resrray(1,i,j,k)),',',modrray(1,i,j,k),',',&
                                &real(resrray(2,i,j,k)),',',aimag(resrray(2,i,j,k)),',',modrray(2,i,j,k),',',&
                                &real(resrray(3,i,j,k)),',',aimag(resrray(3,i,j,k)),',',modrray(3,i,j,k),',',&
                                &real(resrray(4,i,j,k)),',',aimag(resrray(4,i,j,k)),',',modrray(4,i,j,k),',',&
                                &real(resrray(5,i,j,k)),',',aimag(resrray(5,i,j,k)),',',modrray(5,i,j,k)
                enddo
            enddo
        enddo
        close(37)
    end subroutine format_file

    subroutine test_format_file(comm)
        use penf,only:R_P
        implicit none
        complex(R_P),dimension(:,:,:),allocatable :: resrray
        real(R_P),dimension(:,:,:),allocatable :: modrray
        PetscScalar,dimension(:,:,:),pointer :: tmp 
        integer :: xs, ys, zs, xl, yl, zl
        PetscInt,intent(in) :: comm
        integer :: i,j,k,l
        DM :: ResDA
        Vec :: Res
        call DMDACreate2d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
        &                 DMDA_STENCIL_BOX, in, jn,  1, PETSC_DECIDE, &
        &                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ResDA, ierr)
        call DMSetUp(ResDA, ierr)
        call DMGetGlobalVector(ResDA, Res, ierr)
        call DMDAGetCorners(ResDA,xs,ys,zs,xl,yl,zl,ierr)
        call PetscViewerBinaryOpen(comm, "lpse/lpse.petsc",FILE_MODE_READ, Viewer, ierr)
        call VecLoad(Res, Viewer, ierr)
        call PetscViewerDestroy(Viewer, ierr)
        allocate(resrray(5,in,jn))
        allocate(modrray(5,in,jn))
        call DMDAVecGetArrayReadF90(ResDA, Res, tmp, ierr)
        resrray=tmp
        call DMDAVecRestoreArrayReadF90(ResDA, Res, tmp, ierr)
        modrray=abs(resrray)
        open(38, file='out/lpse.csv',action='write',status='replace')
        write(38,*) "rho_r,rho_i,rho_m,u_r,u_i,u_m,v_r,v_i,v_m,w_r,w_i,w_m,T_r,T_i,T_m"
        do i=1,in 
            do j=1,jn 
                write(38,*)  real(resrray(1,i,j)),',',aimag(resrray(1,i,j)),',',modrray(1,i,j),',',&
                            &real(resrray(2,i,j)),',',aimag(resrray(2,i,j)),',',modrray(2,i,j),',',&
                            &real(resrray(3,i,j)),',',aimag(resrray(3,i,j)),',',modrray(3,i,j),',',&
                            &real(resrray(4,i,j)),',',aimag(resrray(4,i,j)),',',modrray(4,i,j),',',&
                            &real(resrray(5,i,j)),',',aimag(resrray(5,i,j)),',',modrray(5,i,j)
            enddo
        enddo
        close(38)
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
        call PetscPrintf(comm," ========================================  o   o   ooo    ooo   o   o  =====\n", ierr)
        call PetscPrintf(comm,"                                           o   o                o   o \n",ierr)
        call PetscPrintf(comm,"                                            ooo     oooooooo     ooo\n",ierr)
        call PetscPrintf(comm,"               ooo    ooo                        o            o\n",ierr)
        call PetscPrintf(comm,"              o   o  o   o                      o              o\n",ierr)
        call PetscPrintf(comm,"        ooo   o   o  o   o   ooo                 o            o\n",ierr)
        call PetscPrintf(comm,"       o   o   ooo    ooo   o   o                   oooooooo\n",ierr)
        call PetscPrintf(comm,"       o   o                o   o\n",ierr)
        call PetscPrintf(comm,"        ooo     oooooooo     ooo\n",ierr)
        call PetscPrintf(comm,"             o            o\n",ierr)
        call PetscPrintf(comm,"            o              o\n",ierr)
        call PetscPrintf(comm,"             o            o\n",ierr)
        call PetscPrintf(comm,"                oooooooo\n",ierr)
    end subroutine print_info
end module mod_petsc_output  
