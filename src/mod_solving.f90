#include <slepc/finclude/slepc.h>

module mod_solving
! ----------------------------------------------------
! 
!  这个模块是主工作流。
! 
!       1.call working(comm) 主工作流。
! 
!           1).call allocate_memory() 分配内存。
!
!           2).call linear_equations(comm) 求解线性系统。
!
!               a).call dolphin_ready(comm,level) 设置免矩阵求解方法。
!
!               b).call whale_ready(comm,level) 设置显式矩阵求解方法。
!
!           3).call deallocate_memory() 释放内存。
! 
! ----------------------------------------------------
    use mod_parameters
    use mod_metrics
    use mod_forming
    use mod_points
    use petsc
    public :: working
    private
    PetscErrorCode :: ierr
    contains
    subroutine working(comm)
        PetscInt,intent(in) :: comm
        call metric_coefficient(comm)
        call partial_derivatives(comm)
        select case (solver_mode)
            case(0)
                call linear_equations(comm)
            case(1)
                call nonlinear_equations(comm)
        end select
        deallocate(bf)
    end subroutine working

    subroutine linear_equations(comm)
        implicit none
        integer,intent(in) :: comm
        logical :: Matrix_Free
        Matrix_Free=.False.
        call VecDuplicate(Turtle,RHS,ierr)
        call VecZeroEntries(RHS,ierr)
        call set_right_hand_side(comm)
        select case (Matrix_Free)
            case (.True.)
                call dolphin_coming(comm)
                !call dolphin_ready(comm,0)
            case (.False.)
                call whale_coming(comm)
                call whale_ready(comm,0)
        end select
        call VecDestroy(RHS,ierr)
    end subroutine linear_equations

    subroutine whale_ready(comm,level)
        implicit none
        integer,intent(in) :: level
        PetscInt,intent(in) :: comm
        real(8) :: rtol
        KSP :: ksp 
        PC :: pc
        rtol = 1e-8
        if(level==0)then ! 如果level是0，那么不使用多重网格
            call KSPCreate(comm,ksp,ierr)
            call KSPSetOperators(ksp,Whale,Whale,ierr)
            call KSPSetType(ksp,KSPFGMRES,ierr)
            call KSPSetInitialGuessNonzero(ksp,initial_guess,ierr)
            call KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization,ierr)
            call KSPGMRESSetRestart(ksp,40,ierr)
            call KSPSetTolerances(ksp,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
            ! call KSPSetDiagonalScale(ksp,PETSC_TRUE,ierr)
            ! call KSPSetDiagonalScaleFix(ksp,PETSC_TRUE,ierr)
            call KspGetPC(ksp,pc,ierr)
            call PCSetType(pc,PCASM,ierr)
            call PCASMSetOverlap(pc,10,ierr)
            ! call PCASMSetSubMatType(pc,MATBAIJ,ierr)
            call PCFactorSetUseInPlace(pc,PETSC_TRUE,ierr)
            call PCSetFromOptions(pc,ierr)
            call PCSetUp(pc,ierr)
            call KSPSetFromOptions(ksp,ierr)
            call KSPSetUp(ksp,ierr)
            call PetscPrintf(comm, " -----------------------------------\n", ierr)
            call PetscPrintf(comm, "                计算中              \n", ierr)
            call PetscPrintf(comm, " -----------------------------------\n", ierr)
            call KSPSolve(ksp,RHS,turtle,ierr)
            call PetscPrintf(comm, " -----------------------------------\n", ierr)
            call PetscPrintf(comm, "               计算完毕              \n", ierr)
            call PetscPrintf(comm, " -----------------------------------\n", ierr)
            call PetscPrintf(comm, " -----------------------------------\n", ierr)
            call PetscPrintf(comm, "               KSP 信息            \n", ierr)
            call PetscPrintf(comm, " -----------------------------------\n", ierr)
            call PetscPrintf(comm, "\n", ierr)
            call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
        elseif(level>0)then
            call KSPCreate(comm,ksp,ierr)
            call KSPSetOperators(ksp,Whale,Whale,ierr)
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

    subroutine nonlinear_equations(comm)
        implicit none 
        integer, intent(in) :: comm
        call shark_coming(comm)
        call shark_ready(comm)
    end subroutine nonlinear_equations

    subroutine shark_ready(comm)
        implicit none
        integer, intent(in) :: comm
        real(8) :: rtol
        SNES :: snes 
        KSP :: ksp 
        PC :: pc
        rtol = 1e-8
        call SNESCreate(comm,snes,ierr)
        call SNESSetFunction(snes,PETSC_NULL_VEC,RHS_with_BC,0,ierr)
        call SNESSetJacobian(snes,Shark,Shark,shark_growing_up,0,ierr)
        call SNESGetKSP(snes,ksp,ierr)
        call KSPSetTolerances(ksp,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)    
        call KSPGetPC(ksp,pc,ierr)
        ! call PCSetType(pc,PCASM,ierr)
        call PCSetFromOptions(pc,ierr)
        call PCSetUp(pc,ierr)
        call KSPSetFromOptions(ksp,ierr)
        call KSPSetUp(ksp,ierr)
        call SNESSetType(snes,SNESNEWTONLS,ierr)
        call SNESSetFromOptions(snes,ierr)
        call SNESSetUp(snes,ierr)
        call SNESSolve(snes,PETSC_NULL_VEC,Turtle,ierr)
        call SNESDestroy(snes,ierr)
    end subroutine shark_ready

end module mod_solving
