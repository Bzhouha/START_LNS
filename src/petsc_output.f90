!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_petsc_output
    use petsc
    use global_parameters
    implicit none
    private
    public :: result_to_file
    contains
    subroutine result_to_file(comm)
        implicit none
        PetscInt, INTENT(in) :: comm
        PetscViewer :: Viewer
        PetscErrorCode  :: ierr
        call signal_ending(comm)
        call PetscViewerBinaryOpen(comm, "out//Turtle.petsc", FILE_MODE_WRITE, Viewer, ierr)
        call VecView(turtle, Viewer, ierr)
        call PetscViewerDestroy(Viewer, ierr)
        call print_info(comm)
    end subroutine result_to_file
    
    subroutine signal_ending(comm)
        implicit none 
        PetscInt,intent(in) :: comm 
        PetscErrorCode :: ierr 
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
        PetscErrorCode :: ierr
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
