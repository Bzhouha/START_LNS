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
!               1).call solve_ksp_mf(comm,level) 设置KSP免矩阵求解方法。
!
!               2).call solve_ksp(comm,level) 设置KSP显式矩阵求解方法。
!
!           2.call nonlinear_equations(comm) 基于SNES求解方程。
!
!               call solve_snes(comm) 设置SNES求解参数。
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
        use mod_parameters,only : solver_mode,whale,shark,dolphin,turtle,RHS
        implicit none
        PetscInt,intent(in) :: comm
        call PetscPrintf(comm, "\n ----------------------------------\n", ierr)
        call PetscPrintf(comm, "              DStream            \n", ierr)
        call PetscPrintf(comm,"\n   Data :: Preparation\n",ierr)
        call metric_coefficient(comm)
        call partial_derivatives(comm)
        select case (solver_mode)
            case('ksp')
                call linear_equations(comm,whale,turtle,RHS)
            case('snes')
                call nonlinear_equations(comm,shark,turtle,snes_fx4o,RHS)
        end select
    end subroutine dstream

    subroutine linear_equations(comm,mat,x,r)
        use mod_parameters,only : dolphin,whale
        implicit none
        integer,intent(in) :: comm
        Vec :: x,r
        Mat :: mat
        call PetscPrintf(comm,"\n   KSP :: Matrix\n",ierr)
        call form_mat(comm,mat)
        call ksp_rhs(comm,r,ierr)
        if(mat==dolphin)then
            call solve_ksp_mf(comm,mat,x,r,0)
        else
            call solve_ksp(comm,mat,x,r,0)
        endif
    end subroutine linear_equations

    subroutine nonlinear_equations(comm,mat,x,fx,r)
        implicit none
        integer, intent(in) :: comm
        external :: fx
        Mat :: mat
        Vec :: x,r
        call PetscPrintf(comm,"\n   SNES :: Jacobi&fx\n",ierr)
        call form_mat(comm,mat)
        call ksp_rhs(comm,r,ierr)
        call solve_snes(comm,mat,x,fx,r)
    end subroutine nonlinear_equations

    subroutine solve_ksp_mf(comm,mat,x,r,level)
        implicit none
        integer,intent(in) :: level
        PetscInt,intent(in) :: comm
        Mat,intent(in) :: mat
        Vec,intent(in) :: r
        Vec,intent(inout) :: x
        call cleanup()
        call VecDestroy(r,ierr)
    end subroutine solve_ksp_mf

    subroutine solve_ksp(comm,mat,x,r,level)
        use mod_parameters,only : init_guess_flg,meshDA
        implicit none
        integer,intent(in) :: level
        PetscInt,intent(in) :: comm
        Mat,intent(in) :: mat
        Vec,intent(in) :: r
        Vec,intent(inout) :: x
        real(8) :: rtol
        DM,pointer :: da_list(:)
        KSP :: ksp
        PC :: pc
        rtol = 1e-8
        if(level==0)then ! 如果level是0，那么不使用多重网格
            call PetscPrintf(comm, "\n   KSP :: Solve\n\n", ierr)
            call KSPCreate(comm,ksp,ierr)
            call KSPSetOperators(ksp,mat,mat,ierr)
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
            call KSPSolve(ksp,r,x,ierr)
            call PetscPrintf(comm, "\n", ierr)
            call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
        elseif(level>0)then
            call KSPCreate(comm,ksp,ierr)
            call KSPSetOperators(ksp,mat,mat,ierr)
            ! call KSPSetType(ksp,KSPFGMRES,ierr)
            call KSPGetPC(ksp,pc,ierr)
            call PCSetType(pc,PCMG,ierr)
            call KSPSetFromOptions(ksp,ierr)
            call KSPSetUp(ksp,ierr)
            ! #####设置网格层数
            ! call PCMGSetLevels(pc,level,comm,ierr)
            ! call PCMGSetType(pc,PC_MG_MULTIPLICATIVE,ierr)
            ! call PCMGSetCycleType(pc,PC_MG_CYCLE_W,ierr)
            ! #####设置网格粗化
            allocate(da_list(level))
            ! da_list(1)=meshDA
            ! call DMCoarsen(meshDA,comm,da_list(2),ierr)
            ! call DMCoarsen(da_list(2),comm,da_list(3),ierr)
            !
            ! call PCMGSetLevels(pc,level,comm,ierr)
            ! call PCMGSetGalerkin(pc,PC_MG_GALERKIN_PMAT,ierr)
            ! #####设置光滑子
            ! #####设置每层迭代矩阵
        endif
        call VecDestroy(r,ierr)
        call KSPDestroy(ksp,ierr)
    end subroutine solve_ksp

    subroutine solve_snes(comm,jac,x,fx,r)
        implicit none
        integer,intent(in) :: comm
        Mat,intent(in) :: jac
        Vec,intent(inout) :: x
        Vec,intent(in) :: r
        character(len=256) :: reason
        external :: fx
        real(8) :: rtol
        SNES :: snes
        KSP :: ksp
        Vec :: t
        PC :: pc
        rtol = 1e-8
        call PetscPrintf(comm, "\n   SNES :: Solve\n\n", ierr)
        call VecDuplicate(x,t,ierr)
        call SNESCreate(comm,snes,ierr)
        call SNESSetFunction(snes,t,fx,0,ierr)
        call SNESSetJacobian(snes,jac,jac,PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER,ierr)
        ! call SNESSetType(snes,SNESNEWTONTR,ierr)
        call SNESGetKSP(snes,ksp,ierr)
        call KSPSetType(ksp,KSPFGMRES,ierr)
        call KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization,ierr)
        ! call KSPGMRESSetRestart(ksp,40,ierr)
        call KSPSetTolerances(ksp,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
        call KSPGetPC(ksp,pc,ierr)
        ! call PCSetType(pc,PCASM,ierr)
        ! call PCASMSetOverlap(pc,10,ierr)
        ! call PCFactorSetUseInPlace(pc,PETSC_TRUE,ierr)
        call PCSetFromOptions(pc,ierr)
        call KSPSetFromOptions(ksp,ierr)
        ! call SNESSetTolerances(snes,PETSC_DEFAULT_REAL,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_INTEGER,ierr)
        call SNESSetFromOptions(snes,ierr)
        call SNESSetUp(snes,ierr)
        call SNESSolve(snes,r,x,ierr)
        call SNESGetConvergedReasonString(snes,reason,ierr)
        call PetscPrintf(comm, "\n", ierr)
        call PetscPrintf(comm, "The Reason of Converged is : "//reason//"\n", ierr)
        call SNESView(snes,PETSC_VIEWER_STDOUT_WORLD,ierr)
        call VecDestroy(t,ierr)
        call SNESDestroy(snes,ierr)
        call cleanup()
    end subroutine solve_snes

end module mod_solving
