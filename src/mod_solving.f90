#include <slepc/finclude/slepc.h>

module mod_solving
! ----------------------------------------------------
!
!  这个模块是主工作流。
!
!       call dstream(comm) 数据流。
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
    public :: dstream
    private
    PetscErrorCode :: ierr
    contains
    subroutine dstream(comm)
        use mod_parameters,only : solver_mode
        implicit none
        PetscInt,intent(in) :: comm
        call PetscPrintf(comm, "\n ----------------------------------\n", ierr)
        call PetscPrintf(comm, "              DStream            \n", ierr)
        call PetscPrintf(comm,"\n   Data :: Preparation\n",ierr)
        call metric_coefficient(comm)
        call partial_derivatives(comm)
        select case (solver_mode)
            case('ksp')
                call linear_equations(comm)
            case('snes')
                call nonlinear_equations(comm)
        end select
    end subroutine dstream

    subroutine linear_equations(comm)
        use mod_parameters,only : ksp_mat_free_flg
        implicit none
        integer,intent(in) :: comm
        ksp_mat_free_flg = .False.
        select case (ksp_mat_free_flg)
            case (.True.)
                call dolphin_coming(comm,ierr)
                call dolphin_ready(comm,0)
            case (.False.)
                call whale_coming(comm,ierr)
                call whale_ready(comm,0)
        end select
    end subroutine linear_equations

    subroutine nonlinear_equations(comm)
        implicit none
        integer, intent(in) :: comm
        call shark_coming(comm,ierr)
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
            call PetscPrintf(comm, "\n   KSP :: Solve\n\n", ierr)
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
            call KSPSolve(ksp,RHS,turtle,ierr)
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
        use mod_parameters,only : turtle,RHS,meshDA,bell
        implicit none
        integer,intent(in) :: level
        PetscInt,intent(in) :: comm
        call DMRestoreLocalVector(meshDA,bell,ierr)
        call cleanup()
        call VecDestroy(RHS,ierr)
    end subroutine dolphin_ready

    subroutine shark_ready(comm)
        use mod_parameters,only : shark,turtle,meshDA,bell,whale
        implicit none
        integer,intent(in) :: comm
        real(8) :: rtol
        SNES :: snes
        KSP :: ksp
        Vec :: r
        PC :: pc
        rtol = 1e-8
        call PetscPrintf(comm, "\n   SNES :: Solve\n\n", ierr)
        call VecDuplicate(turtle,r,ierr)
        call SNESCreate(comm,snes,ierr)
        call SNESSetFunction(snes,r,snes_fx,0,ierr)
        ! call SNESSetJacobian(snes,shark,shark,snes_jac,0,ierr)
        call SNESSetJacobian(snes,shark,shark,PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER,ierr)
        call SNESGetKSP(snes,ksp,ierr)
        call KSPSetType(ksp,KSPFGMRES,ierr)
        call KSPSetTolerances(ksp,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
        call KSPGetPC(ksp,pc,ierr)
        call PCSetType(pc,PCASM,ierr)
        call PCSetFromOptions(pc,ierr)
        call KSPSetFromOptions(ksp,ierr)
        ! call SNESSetType(snes,SNESNEWTONLS,ierr)
        call SNESSetTolerances(snes,PETSC_DEFAULT_REAL,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_INTEGER,ierr)
        call SNESSetFromOptions(snes,ierr)
        call SNESSetUp(snes,ierr)
        call SNESSolve(snes,PETSC_NULL_VEC,turtle,ierr)
        call PetscPrintf(comm, "\n", ierr)
        call SNESView(snes,PETSC_VIEWER_STDOUT_WORLD,ierr)
        call SNESDestroy(snes,ierr)
        call cleanup()
        call DMRestoreLocalVector(meshDA,bell,ierr)
    end subroutine shark_ready
end module mod_solving
