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
	public :: working
	Vec :: RHS
	contains
	subroutine working(comm)
		PetscInt,intent(in) :: comm
		PetscErrorCode :: ierr 
		call allocate_memory()
		call metric_coefficient(comm)
		call partial_derivatives(comm)
		call linear_equations(comm)
		call deallocate_memory()
	end subroutine working

	subroutine allocate_memory()
		implicit none 
		call DMCreateGlobalVector(meshDA,Turtle,ierr)
		call VecDuplicate(Turtle,RHS,ierr)
		call VecZeroEntries(Turtle,ierr)
		call VecZeroEntries(RHS,ierr)
	end subroutine allocate_memory

	subroutine linear_equations(comm)
		implicit none
		integer,intent(in) :: comm
		logical :: Matrix_Free
		Matrix_Free=.False.
		call set_right_hand_side(comm)
		select case (Matrix_Free)
			case (.True.)
				call dolphin_coming(comm)
				!call DolphinReady(comm,0)
			case (.False.)
				call whale_coming(comm)
				call whale_ready(comm,0)
				! call print_result()
		end select
	end subroutine linear_equations

	subroutine whale_ready(comm,level)
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
			call PetscPrintf(comm, "               计算中...             \n", ierr)
			call PetscPrintf(comm, " -----------------------------------\n", ierr)
			call KSPSolve(ksp,RHS,turtle,ierr)
			call PetscPrintf(comm, " -----------------------------------\n", ierr)
			call PetscPrintf(comm, "              计算完毕。              \n", ierr)
			call PetscPrintf(comm, " -----------------------------------\n", ierr)
			call PetscPrintf(comm, " -----------------------------------\n", ierr)
			call PetscPrintf(comm, "             K S P 信 息             \n", ierr)
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
		call KSPDestroy(ksp,ierr)
		call VecDestroy(b,ierr)
		call VecDestroy(u,ierr)
	end subroutine whale_ready

	subroutine set_right_hand_side(comm)
		implicit none 
		PetscScalar,pointer :: disturb_array(:,:,:)
		PetscScalar,pointer :: RHS_array(:,:,:,:)
		PetscInt,intent(in) :: comm
		integer :: j,k 
		call DMDAVecGetArrayF90(meshDA,RHS,RHS_array,ierr)
		call DMDAVecGetArrayReadF90(disturbDA,disturb,disturb_array,ierr)
		do k=ks,ke 
			do j=js,je 
				RHS_array(:,0,j,k)=disturb_array(:,j,k)
			enddo
		enddo
		call DMDAVecRestoreArrayF90(meshDA,RHS,RHS_array,ierr)
		call DMDAVecRestoreArrayReadF90(disturbDA,disturb,disturb_array,ierr)
		call MPI_Barrier(comm,ierr)
	end subroutine set_right_hand_side

	subroutine deallocate_memory()
		implicit none
		deallocate(bf)
		call VecDestroy(RHS,ierr)
	end subroutine deallocate_memory

	subroutine print_result()
		implicit none
		PetscScalar,pointer :: tmp(:,:,:,:)
		PetscErrorCode :: ierr
		integer :: l,i,j,k
		call DMDAVecGetArrayReadF90(meshDA,turtle,tmp,ierr)
		do l=0,4
			do k=ks,ke
				do j=js,je
					write(*,*) (real(tmp(l,i,j,k)),i=is,ie)
				enddo
				write(*,*) 
			enddo
			write(*,*) "---"
		enddo
		call DMDAVecRestoreArrayReadF90(meshDA,turtle,tmp,ierr)
	end subroutine print_result
end module mod_solving
