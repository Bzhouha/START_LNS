#include <slepc/finclude/slepc.h>

module mod_solving
! ----------------------------------------------------
! 
!  这个模块是主工作流。
! 
!       1.call working(comm) 主工作流。
! 
! 			1).call allocate_memory() 分配内存。
!
!			2).call linear_equations(comm) 求解线性系统。
!
!				a).call dolphin_ready(comm,level) 设置免矩阵求解方法。
!
!				b).call whale_ready(comm,level) 设置显式矩阵求解方法。
!
!			3).call set_right_hand_side(comm) 设置右边量，即来流。
!
! 			4).call deallocate_memory() 释放内存。
! 
! ----------------------------------------------------
	use petsc
	use mod_metrics
	use mod_points
	use mod_forming
	use mod_parameters
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
				!call dolphin_ready(comm,0)
			case (.False.)
				call whale_coming(comm)
				call whale_ready(comm,0)
		end select
	end subroutine linear_equations

	subroutine whale_ready(comm,level)
		implicit none
		integer,intent(in) :: level
		PetscInt,intent(in) :: comm
		KSP :: ksp 
		PC :: pc
		if(level==0)then ! 如果level是0，那么不使用多重网格
			call KSPCreate(comm,ksp,ierr)
			call KSPSetOperators(ksp,Whale,Whale,ierr)
			call KSPSetType(ksp,KSPFGMRES,ierr)
			call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)
			call KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization,ierr)
			call KspGetPC(ksp,pc,ierr)
			call PCSetType(pc,PCASM,ierr)
			call PCASMSetOverlap(pc,10,ierr)
			! call PCFactorSetUseInPlace(pc,PETSC_TRUE,ierr)
			call PCSetFromOptions(pc,ierr)
			call PCSetUp(pc,ierr)
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
			call PetscPrintf(comm, "\n", ierr)
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
	end subroutine whale_ready

	subroutine set_right_hand_side(comm)
		implicit none 
		PetscScalar,pointer :: RHS_array(:,:,:,:)
		PetscInt,intent(in) :: comm
		integer :: j,k 
		if (is==0) then
			call DMDAVecGetArrayF90(meshDA,RHS,RHS_array,ierr)
				do k=ks,ke 
					do j=js,je 
						RHS_array(:,0,j,k)=disturb(:,j,k)
					enddo
				enddo
			call DMDAVecRestoreArrayF90(meshDA,RHS,RHS_array,ierr)
		endif
		call MPI_Barrier(comm,ierr)
		deallocate(disturb)
	end subroutine set_right_hand_side

	subroutine deallocate_memory()
		implicit none
		deallocate(bf)
		call VecDestroy(RHS,ierr)
	end subroutine deallocate_memory

end module mod_solving
