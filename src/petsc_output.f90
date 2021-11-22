!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_petsc_output
    use petsc
    use global_parameters
    implicit none
    private
    public :: petsc_output
    contains
    subroutine petsc_output(comm)
        implicit none
        PetscInt, INTENT(in) :: comm
        PetscViewer :: Viewer
        PetscErrorCode  :: ierr
        call LightTurningOff(comm)
        call PetscViewerBinaryOpen(comm, trim(FileLocation)//"out//"//"Turtle.petsc", FILE_MODE_WRITE, Viewer, ierr)
        call VecView(turtle, Viewer, ierr)
        call PetscViewerDestroy(Viewer, ierr)
        call CurtainCall(comm)
    end subroutine petsc_output
    subroutine LightTurningOff(comm)
        implicit none
        PetscInt, INTENT(in) :: comm
        PetscErrorCode  :: ierr
        call PetscPrintf(comm, " ----------------------------------\n", ierr)
        call PetscPrintf(comm, "              输出结果               \n", ierr)
        call PetscPrintf(comm, " ----------------------------------\n", ierr)
        call PetscPrintf(comm, "\n", ierr)
        call PetscPrintf(comm, " 输出解向量...\n", ierr)
        end subroutine LightTurningOff
    subroutine CurtainCall(comm)
        implicit none
        PetscInt,INTENT(in) :: comm
        PetscErrorCode :: ierr
        call PetscPrintf(comm, "   输出结束。                                      ooo    ooo\n", ierr)
        call PetscPrintf(comm, "                                                  o   o  o   o\n", ierr)
        call PetscPrintf(comm, " ========================================   ooo   o   o  o   o   ooo   =====\n", ierr)
        call PetscPrintf(comm,"                                           o   o   ooo    ooo   o   o\n",ierr)
        call PetscPrintf(comm,"                                           o   o                o   o\n",ierr)
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
    end subroutine CurtainCall
end module mod_petsc_output  
