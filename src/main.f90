!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

program main
	use petsc
	use mod_parameters,only:rank,size
	use mod_loading
	use mod_solving
	use mod_output
	implicit none
	PetscErrorCode :: ierr
    
	call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
	if (ierr /= 0) then
		write(*,*) 'PetscInitialize failed'
		stop
	endif
	call mpi_comm_rank(PETSC_COMM_WORLD,rank,ierr) 
	call mpi_comm_size(PETSC_COMM_WORLD,size,ierr) 
    
	call loading_data(PETSC_COMM_WORLD)
	call working(PETSC_COMM_WORLD)              
	call result_to_file(PETSC_COMM_WORLD)
     
	call PetscFinalize(ierr)
end program main
