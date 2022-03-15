#include <slepc/finclude/slepc.h>

module mod_solving
! ----------------------------------------------------
! 
!  这个模块是主工作流。
! 
!       call working(comm) 主工作流。
!
!           1.call linear_equations(comm) 基于KSP求解方程。
!
!               1).call dolphin_ready(comm,level) 设置KSP免矩阵求解方法。
!
!               2).call whale_ready(comm,level) 设置KSP显式矩阵求解方法。
!
!           2.call nonlinear_equations(comm) 基于SNES求解方程。
!
!               call shark_ready(comm) 设置SNES求解参数。
! 
! ----------------------------------------------------
    use mod_metrics
    use mod_forming
    use mod_points
    use petsc
    public :: working
    private
    PetscErrorCode :: ierr
    contains
    subroutine working(comm)
        use mod_parameters,only : solver_mode
        implicit none
        PetscInt,intent(in) :: comm
        call metric_coefficient(comm)
        call partial_derivatives(comm)
        select case (solver_mode)
            case(0)
                call linear_equations(comm)
            case(1)
                call nonlinear_equations(comm)
        end select
    end subroutine working

    subroutine linear_equations(comm)
        use mod_parameters,only : ksp_mat_free_flg
        implicit none
        integer,intent(in) :: comm
        select case (ksp_mat_free_flg)
            case (.True.)
                call dolphin_coming(comm)
                call dolphin_ready(comm,0)
            case (.False.)
                call whale_coming(comm)
                call whale_ready(comm,0)
        end select
    end subroutine linear_equations

    subroutine nonlinear_equations(comm)
        implicit none 
        integer, intent(in) :: comm
        call shark_coming(comm)
        call shark_ready(comm)
    end subroutine nonlinear_equations

    subroutine whale_ready(comm,level)
        use mod_parameters,only : whale,turtle,RHS,init_guess_flg
        implicit none
        integer,intent(in) :: level
        PetscInt,intent(in) :: comm
        real(8) :: rtol
        KSP :: ksp 
        PC :: pc
        rtol = 1e-8
        if(level==0)then ! 如果level是0，那么不使用多重网格
            call KSPCreate(comm,ksp,ierr)
            call KSPSetOperators(ksp,whale,whale,ierr)
            call KSPSetType(ksp,KSPFGMRES,ierr)
            call KSPSetInitialGuessNonzero(ksp,init_guess_flg,ierr)
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
        call VecDestroy(RHS,ierr)
        call KSPDestroy(ksp,ierr)
    end subroutine whale_ready

    subroutine dolphin_ready(comm,level)
        use mod_parameters,only : turtle,RHS,meshDA,tinkle_bell
        implicit none
        integer,intent(in) :: level
        PetscInt,intent(in) :: comm
        call DMRestoreLocalVector(meshDA,tinkle_bell,ierr)
        call deallocate_bfinfo_and_metrics()
        call VecDestroy(RHS,ierr)
    end subroutine dolphin_ready

    subroutine shark_ready(comm)
        use mod_parameters,only : shark,turtle,meshDA,tinkle_bell
        implicit none
        integer, intent(in) :: comm
        real(8) :: rtol
        SNES :: snes 
        KSP :: ksp 
        PC :: pc
        rtol = 1e-8
        call SNESCreate(comm,snes,ierr)
        call SNESSetFunction(snes,PETSC_NULL_VEC,RHS_with_BC,0,ierr)
        call SNESSetJacobian(snes,shark,shark,shark_growing_up,0,ierr)
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
        call SNESSolve(snes,PETSC_NULL_VEC,turtle,ierr)
        call SNESDestroy(snes,ierr)
        call deallocate_bfinfo_and_metrics()
        call DMRestoreLocalVector(meshDA,tinkle_bell,ierr)
    end subroutine shark_ready
end module mod_solving
