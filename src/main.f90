!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

program main
	use petsc
	use mod_loaders
	use mod_petsc_viewer
	use mod_petsc_loader
	use mod_petsc_output
	use mod_solving
	use global_parameters,only:rank,size
	use global_parameters,only:gridfile,flowfile
	use mod_cfgio_adapter
	implicit none
	character(len=256) :: cfg_file
	PetscErrorCode :: ierr
	PetscBool :: set

	call PetscInitialize(PETSC_NULL_CHARACTER,ierr) ! 初始化PetsC
	if (ierr /= 0) then
		write(*,*) 'PetscInitialize failed'
		stop
	endif

	call mpi_comm_rank(PETSC_COMM_WORLD,rank,ierr) 
	call mpi_comm_size(PETSC_COMM_WORLD,size,ierr) 

	call LoadingAndSayHi(PETSC_COMM_WORLD)
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

	call WeBcastSoHard(PETSC_COMM_WORLD)
	call MPI_Barrier(PETSC_COMM_WORLD,ierr)

	call StartWorkingNow(PETSC_COMM_WORLD)
	call load_info(PETSC_COMM_WORLD)
	call MPI_Barrier(PETSC_COMM_WORLD,ierr)
	call solve(PETSC_COMM_WORLD)
	call MPI_Barrier(PETSC_COMM_WORLD,ierr)
	
	call AndWeAreDone(PETSC_COMM_WORLD)
	call petsc_output(PETSC_COMM_WORLD)
	call MPI_Barrier(PETSC_COMM_WORLD,ierr)

	call PetscFinalize(ierr)
end program main

subroutine LoadingAndSayHi(comm)
	implicit none
	PetscInt,intent(in) :: comm
	PetscErrorCode :: ierr 
	call PetscPrintf(comm, "\n", ierr)
	call PetscPrintf(comm, " ===========================================================================\n", ierr)
	call PetscPrintf(comm, " =                                 读    取                                = \n", ierr)
	call PetscPrintf(comm, " ===========================================================================\n", ierr)
	call PetscPrintf(comm, " 「 单 进 程 」\n",ierr)
end subroutine LoadingAndSayHi

subroutine WeBcastSoHard(comm)
	use global_parameters
	implicit none 
	PetscInt,intent(in) :: comm
	PetscErrorCode :: ierr
	call MPI_Bcast(in,1,MPI_INT,0,comm,ierr)
	call MPI_Bcast(jn,1,MPI_INT,0,comm,ierr)
	call MPI_Bcast(kn,1,MPI_INT,0,comm,ierr)
	call MPI_Bcast(BC_type,1,MPI_INT,0,comm,ierr)
	call MPI_Bcast(Ma,1,MPI_DOUBLE,0,comm,ierr)
	call MPI_Bcast(Re,1,MPI_DOUBLE,0,comm,ierr)
	call MPI_Bcast(Te,1,MPI_DOUBLE,0,comm,ierr)
	call MPI_Bcast(Omega,1,MPI_DOUBLE_COMPLEX,0,comm,ierr)
end subroutine WeBcastSoHard

subroutine StartWorkingNow(comm)
	implicit none
	PetscInt,intent(in) :: comm 
	PetscErrorCode :: ierr 
	call PetscPrintf(comm, "\n", ierr)
	call PetscPrintf(comm, " ===========================================================================\n", ierr)
	call PetscPrintf(comm, " =                                 计    算                                = \n", ierr)
	call PetscPrintf(comm, " ===========================================================================\n", ierr)
	call PetscPrintf(comm, " 「 多 进 程 」\n",ierr)
end subroutine StartWorkingNow

subroutine AndWeAreDone(comm)
	implicit none 
	PetscInt,intent(in) :: comm 
	PetscErrorCode :: ierr 
	call PetscPrintf(comm, "\n", ierr)
	call PetscPrintf(comm, " ===========================================================================\n", ierr)
	call PetscPrintf(comm, " =                                 输    出                                = \n", ierr)
	call PetscPrintf(comm, " ===========================================================================\n", ierr)
end subroutine AndWeAreDone
