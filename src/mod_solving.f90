#include <slepc/finclude/slepc.h>
!#include <petsc/finclude/petscksp.h>

module mod_solving
	use petsc
	use mod_metrics
	use bf_points
	use mod_forming
	use global_parameters
	private
	PetscErrorCode :: ierr
	public :: Working
	Vec :: rhs
	contains
	subroutine Working(comm)
		PetscInt,intent(in) :: comm
		PetscErrorCode :: ierr 
		call DMCreateGlobalVector(meshDA,turtle,ierr)
		call Flatfish(comm)
		call MPI_Barrier(comm,ierr)
		call Tunas(comm)
		call MPI_Barrier(comm,ierr)
		!call DolphinComing(comm)
		call WhaleComing(comm)
		call MPI_Barrier(comm,ierr)
		call solving(comm)
		call DropTheWaste()
	end subroutine Working

	subroutine solving(comm)
		implicit none
		PetscScalar,pointer :: tmp(:,:,:,:)
		PetscInt,intent(in) :: comm
		integer :: l,i,j,k

		call WhaleReady(comm,0)
		
		! call DMDAVecGetArrayReadF90(meshDA,turtle,tmp,ierr)
		! do l=0,4
		! 	do k=ks,ke
		! 		do j=js,je
		! 			write(*,*) (real(tmp(l,i,j,k)),i=is,ie)
		! 		enddo
		! 		write(*,*) 
		! 	enddo
		! 	write(*,*) "---"
		! enddo
		! call DMDAVecRestoreArrayReadF90(meshDA,turtle,tmp,ierr)

	end subroutine solving 

	subroutine WhaleReady(comm,level)
		implicit none
		integer,intent(in) :: level
		PetscInt,intent(in) :: comm
		PetscScalar :: one
		Vec :: b,u
		KSP :: ksp 
		PC :: pc
		one=1.0
		call VecDuplicate(turtle,b,ierr)
		call VecDuplicate(turtle,u,ierr)
		call VecSet(b,one,ierr)
		call MatMult(Whale,b,u,ierr)
		if(level==0)then ! 如果level是0，那么不使用多重网格
			call KSPCreate(comm,ksp,ierr)
			call KSPSetOperators(ksp,Whale,Whale,ierr)
			call KSPSetType(ksp,KSPFGMRES,ierr)
			call KSPSetFromOptions(ksp,ierr)
			call KSPSetUp(ksp,ierr)
			call PetscPrintf(comm, " -----------------------------------\n", ierr)
			call PetscPrintf(comm, "               计算中...      \n", ierr)
			call PetscPrintf(comm, " -----------------------------------\n", ierr)
			call KSPSolve(ksp,u,turtle,ierr)
			call PetscPrintf(comm, " -----------------------------------\n", ierr)
			call PetscPrintf(comm, "              计算完毕。      \n", ierr)
			call PetscPrintf(comm, " -----------------------------------\n", ierr)
			call PetscPrintf(comm, " -----------------------------------\n", ierr)
			call PetscPrintf(comm, "             K S P 信 息           \n", ierr)
			call PetscPrintf(comm, " -----------------------------------\n", ierr)
			call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
		elseif(level>0)then
			call KSPCreate(comm,ksp,ierr)
			call KSPSetOperators(ksp,whale,whale,ierr)
			call KSPSetType(ksp,KSPFGMRES,ierr)
			call KSPGetPC(ksp,pc,ierr)
			call PCSetType(pc,PCMG,ierr)
			call KSPSetFromOptions(ksp,ierr)
			call KSPSetUp(ksp,ierr)
			! #####设置网格层数 
			call PCMGSetLevels(pc,level,comm,ierr)
			call PCMGSetType(pc,PC_MG_MULTIPLICATIVE,ierr)
			call PCMGSetCycleType(pc,PC_MG_CYCLE_W,ierr)
			! #####设置网格粗化 
			! #####设置光滑子 
			! #####设置每层迭代矩阵
		endif
	end subroutine WhaleReady

	subroutine SetRightValues(comm)
		implicit none 
		PetscScalar,pointer :: tmp(:,:,:,:)
		PetscInt,intent(in) :: comm
		integer :: j,k 
		call VecDuplicate(turtle,rhs,ierr)
		if(is==0)then
		call DMDAVecGetArrayF90(meshDA,rhs,tmp,ierr)
		do k=ks,ke 
			do j=js,je 
				tmp(:,0,j,k)=wave(:,j,k)
			enddo
		enddo
		call DMDAVecRestoreArrayF90(meshDA,rhs,tmp,ierr)
		endif 
		call MPI_Barrier(comm,ierr)
		deallocate(wave)
	end subroutine SetRightValues

	subroutine DropTheWaste()
		implicit none
		deallocate(bf)
	end subroutine DropTheWaste
end module mod_solving
