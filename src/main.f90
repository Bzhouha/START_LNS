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
	use mod_cfgio_adapter
	implicit none
	character(len=256) :: cfg_file
	PetscErrorCode :: ierr
	PetscBool :: set

	!                   初始化PetsC                  
	call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
	if (ierr /= 0) then
		write(*,*) 'PetscInitialize failed'
		stop
	endif
	call mpi_comm_rank(PETSC_COMM_WORLD,rank,ierr) 
	call mpi_comm_size(PETSC_COMM_WORLD,size,ierr) 
	call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-f',cfg_file,set,ierr)
	if(.not. set) then
		write(*,*) 'should use -f option to determin the config file.'
		stop
	endif
	!                读取信息               
	call Signal_Loading(PETSC_COMM_WORLD)
	call cfg_loader(trim(cfg_file))
	if(rank==0)then
		call plot3d_load()
		call petsc_viewer(PETSC_COMM_SELF)
		call cfg_writer(trim(cfg_file))
	endif

	!               分发并计算              
	call Signal_Starting(PETSC_COMM_WORLD)
	call LoadingData(PETSC_COMM_WORLD)
	!call Working(PETSC_COMM_WORLD)

	!               输出文件              
	call Signal_Ending(PETSC_COMM_WORLD)
	!call ResultToFile(PETSC_COMM_WORLD)

	!        终止PetsC       
	call PetscFinalize(ierr)
end program main

subroutine Signal_Loading(comm)
	implicit none
	PetscInt,intent(in) :: comm
	PetscErrorCode :: ierr 
	call PetscPrintf(comm, "\n", ierr)
	call PetscPrintf(comm, " ===========================================================================\n", ierr)
	call PetscPrintf(comm, " =                                 读    取                                = \n", ierr)
	call PetscPrintf(comm, " ===========================================================================\n", ierr)
	call PetscPrintf(comm, " 「 单 进 程 」\n",ierr)
end subroutine Signal_Loading

subroutine Signal_Starting(comm)
	implicit none
	PetscInt,intent(in) :: comm 
	PetscErrorCode :: ierr 
	call PetscPrintf(comm, "\n", ierr)
	call PetscPrintf(comm, " ===========================================================================\n", ierr)
	call PetscPrintf(comm, " =                                 计    算                                = \n", ierr)
	call PetscPrintf(comm, " ===========================================================================\n", ierr)
	call PetscPrintf(comm, " 「 多 进 程 」\n",ierr)
end subroutine Signal_Starting

subroutine Signal_Ending(comm)
	implicit none 
	PetscInt,intent(in) :: comm 
	PetscErrorCode :: ierr 
	call PetscPrintf(comm, "\n", ierr)
	call PetscPrintf(comm, " ===========================================================================\n", ierr)
	call PetscPrintf(comm, " =                                 输    出                                = \n", ierr)
	call PetscPrintf(comm, " ===========================================================================\n", ierr)
end subroutine Signal_Ending
