!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

program main
    use mod_files
    use mod_solving
    use petsc
    implicit none
    PetscErrorCode :: ierr

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    if (ierr /= 0) then
        write(*,*) 'PetscInitialize failed'
        stop
    endif

    call istream(PETSC_COMM_WORLD)
    call dstream(PETSC_COMM_WORLD)
    call ostream(PETSC_COMM_WORLD)

    call PetscFinalize(ierr)
end program main
