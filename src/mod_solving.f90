#include <slepc/finclude/slepc.h>

module mod_solving
! ----------------------------------------------------
!
!  这个模块是主工作流。
!
!       call dstream(comm) 数据流。
!
!           1.call ksp_equations(comm,mat,x,r) 迭代格式一：标准线性求解器 Ax=b
!
!               1).call solve_ksp_mf(comm,mat,x,r,level) 设置KSP免矩阵求解过程。
!
!               2).call solve_ksp(comm,mat,x,r,level) 设置KSP显式矩阵求解过程。
!
!           2.call snes_equations(comm,mat,x,fx,r) 迭代格式二：借用SNES模块 Jac x = - F(x)
!
!               1).call solve_snes(comm,jac,x,fx,r) 设置SNES求解过程。
!
!               2).call MySNESConverged(snes,it,xnorm,snorm,fnorm,reason,dummy,ierr) SNES收敛判断函数。
!
!           3.call ksps_equations(comm,mat,x,ksps_fx) 迭代格式三：借用KSP模块 A x = - (Ax-b)
!
!               call solve_ksps(comm,mat,x,fx_rhs) 设置KSPs求解过程。
!
! ----------------------------------------------------
    use mod_metrics
    use mod_forming
    use mod_points
    use petsc

    public :: dstream
    private
    PetscErrorCode :: ierr
    integer, parameter :: NEWTON_LIKE=1001, TIME_DISCRETE=1002

    contains

    subroutine dstream(comm)
        use mod_parameters,only : meshDA,whale,turtle,RHS,solver_mode,shark
        implicit none
        PetscInt,intent(in) :: comm

        call PetscPrintf(comm, "\n ----------------------------------\n", ierr)
        call PetscPrintf(comm, "              DStream            \n", ierr)
        call PetscPrintf(comm, "\n   Data :: Preparation\n", ierr)
        call metric_coefficient(comm)
        call partial_derivatives(comm)
        select case (solver_mode)
            case('ksp')
                call ksp_equations(comm,meshDA,whale,turtle,RHS)
            case('snes')
                call snes_equations(comm,meshDA,whale,turtle,snes_fx4o,RHS)
            case('ksps')
                call ksps_equations(comm,meshDA,whale,turtle,ksps_rhs_fx_b_Ax,NEWTON_LIKE)
        end select
    end subroutine dstream

    subroutine ksp_equations(comm,da,mat,x,r)
        implicit none
        integer,intent(in) :: comm
        Vec :: x,r
        Mat :: mat
        DM :: da
        call PetscPrintf(comm,"\n   KSP :: Matrix\n",ierr)
        call initialize_mat_from_da(comm,da,mat)
        call form_mat_4_precision(mat)
        call set_rhs(comm,da,x,r)
        call solve_ksp(comm,mat,x,r,0)
    end subroutine ksp_equations

    subroutine snes_equations(comm,da,mat,x,fx,r)
        implicit none
        integer, intent(in) :: comm
        external :: fx
        Mat :: mat
        Vec :: x,r
        DM :: da
        call PetscPrintf(comm,"\n   SNES :: Jacobi\n",ierr)
        call initialize_mat_from_da(comm,da,mat)
        call form_mat_2_precision(mat)
        call set_rhs(comm,da,x,r)
        call solve_snes(comm,mat,x,fx,r)
    end subroutine snes_equations

    subroutine ksps_equations(comm,da,mat,x,ksps_fx,type)
        implicit none
        integer,intent(in) :: type
        integer,intent(in) :: comm
        external :: ksps_fx
        Mat :: mat
        Vec :: x
        DM :: da

        call PetscPrintf(comm,"\n   KSPs :: Matrix\n",ierr)
        call initialize_mat_from_da(comm,da,mat)
        call form_mat_2_precision(mat)
        call solve_ksps(comm,mat,x,ksps_fx,type)
    end subroutine ksps_equations

    ! -----------------------------------------------------------------------------------------------------
    !   迭代格式一：标准线性求解器 Ax=b

    subroutine solve_ksp_mf(comm,mat,x,r,level)
        implicit none
        integer,intent(in) :: level
        PetscInt,intent(in) :: comm
        Mat,intent(in) :: mat
        Vec,intent(in) :: r
        Vec,intent(inout) :: x
    end subroutine solve_ksp_mf

    subroutine solve_ksp(comm,mat,x,r,level)
        use mod_parameters,only : init_guess_flg
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

        call PetscPrintf(comm, "\n   KSP :: Solve\n\n", ierr)
        ! Set parameter
        rtol = 1e-8
        if(level==0)then ! 如果level是0，那么不使用多重网格
            ! Create & set KSP
            call KSPCreate(comm,ksp,ierr)
            call KSPSetOperators(ksp,mat,mat,ierr)
            call KSPSetType(ksp,KSPFGMRES,ierr)
            call KSPSetInitialGuessNonzero(ksp,init_guess_flg,ierr)
            call KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization,ierr)
            call KSPGMRESSetRestart(ksp,40,ierr)
            call KSPSetTolerances(ksp,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
            ! Get & set PC
            call KspGetPC(ksp,pc,ierr)
            call PCSetType(pc,PCASM,ierr)
            call PCASMSetOverlap(pc,10,ierr)
            call PCFactorSetUseInPlace(pc,PETSC_TRUE,ierr)
            call PCSetFromOptions(pc,ierr)
            call PCSetUp(pc,ierr)
            call KSPSetFromOptions(ksp,ierr)
            call KSPSetUp(ksp,ierr)
            ! Solve KSP
            call KSPSolve(ksp,r,x,ierr)
            ! View KSP
            call PetscPrintf(comm, "\n", ierr)
            call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
        elseif(level>0)then
            ! Create & set KSP
            call KSPCreate(comm,ksp,ierr)
            call KSPSetOperators(ksp,mat,mat,ierr)
            call KSPSetType(ksp,KSPFGMRES,ierr)
            call KSPGetPC(ksp,pc,ierr)
            call PCSetType(pc,PCMG,ierr)
            call KSPSetFromOptions(ksp,ierr)
            call KSPSetUp(ksp,ierr)
            ! ### 设置网格层数
            ! call PCMGSetLevels(pc,level,comm,ierr)
            ! call PCMGSetType(pc,PC_MG_MULTIPLICATIVE,ierr)
            ! call PCMGSetCycleType(pc,PC_MG_CYCLE_W,ierr)
            ! ### 设置网格粗化
            ! allocate(da_list(level))
            ! da_list(1)=meshDA
            ! call DMCoarsen(meshDA,comm,da_list(2),ierr)
            ! call DMCoarsen(da_list(2),comm,da_list(3),ierr)
            ! call PCMGSetLevels(pc,level,comm,ierr)
            ! call PCMGSetGalerkin(pc,PC_MG_GALERKIN_PMAT,ierr)
            ! ### 设置光滑子
            ! ### 设置每层迭代矩阵
            ! Solve KSP
            call KSPSolve(ksp,r,x,ierr)
            ! View KSP
            call PetscPrintf(comm, "\n", ierr)
            call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
        endif
        ! Destory variables
        call VecDestroy(r,ierr)
        call KSPDestroy(ksp,ierr)
        call cleanup()

    end subroutine solve_ksp

    ! -----------------------------------------------------------------------------------------------------
    !   迭代格式二：借用SNES模块 Jac x = - F(x)

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

        call PetscPrintf(comm, "\n   SNES :: Solve\n\n", ierr)
        ! Set parameters
        rtol = 1e-8
        ! Initialize Vecs
        call VecDuplicate(x,t,ierr)
        call VecZeroEntries(t,ierr)
        ! Create & set SNES
        call SNESCreate(comm,snes,ierr)
        call SNESSetConvergenceTest(snes,MySNESConverged,0,PETSC_NULL_FUNCTION,ierr)
        call SNESSetFunction(snes,t,fx,0,ierr)
        call SNESSetJacobian(snes,jac,jac,PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER,ierr)
        ! Get & set KSP
        call SNESGetKSP(snes,ksp,ierr)
        call KSPSetType(ksp,KSPFGMRES,ierr)
        call KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization,ierr)
        call KSPGMRESSetRestart(ksp,40,ierr)
        call KSPSetTolerances(ksp,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
        ! Get & set PC
        call KSPGetPC(ksp,pc,ierr)
        call PCSetType(pc,PCASM,ierr)
        call PCASMSetOverlap(pc,10,ierr)
        call PCFactorSetUseInPlace(pc,PETSC_TRUE,ierr)
        call PCSetFromOptions(pc,ierr)
        call KSPSetFromOptions(ksp,ierr)
        call SNESSetTolerances(snes,PETSC_DEFAULT_REAL,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_INTEGER,ierr)
        call SNESSetFromOptions(snes,ierr)
        call SNESSetUp(snes,ierr)
        ! Solve
        call SNESSolve(snes,r,x,ierr)
        ! Get convergedr reason
        call SNESGetConvergedReasonString(snes,reason,ierr)
        call PetscPrintf(comm, "\n", ierr)
        call PetscPrintf(comm, "The Reason of Converged is : "//reason//"\n", ierr)
        ! View SNES
        call SNESView(snes,PETSC_VIEWER_STDOUT_WORLD,ierr)
        ! Destory variables
        call VecDestroy(t,ierr)
        call SNESDestroy(snes,ierr)
        call cleanup()

    end subroutine solve_snes

    subroutine MySNESConverged(snes,it,xnorm,snorm,fnorm,reason,dummy,ierr)
        implicit none
        PetscReal :: xnorm,snorm,fnorm,nrm
        SNESConvergedReason :: reason
        PetscErrorCode :: ierr
        PetscInt :: it,dummy
        SNES :: snes
        Vec :: f

        call SNESGetFunction(snes,f,PETSC_NULL_FUNCTION,dummy,ierr)
        call VecNorm(f,NORM_INFINITY,nrm,ierr)
        if (nrm .le. 1.e-5) reason = SNES_CONVERGED_FNORM_ABS
    end subroutine MySNESConverged

    ! -----------------------------------------------------------------------------------------------------
    !   迭代格式三：借用KSP模块 Newton-Like

    subroutine solve_ksps(comm,mat,x,fx_rhs,type)
        implicit none
        character(len=20) :: str_norm
        character(len=6) :: str_count
        integer,intent(in) :: comm
        integer, intent(in) :: type
        Vec,intent(inout) :: x
        PetscScalar :: one,ine
        PetscErrorCode :: ierr
        Mat,intent(in) :: mat
        external :: fx_rhs
        PetscInt :: count
        PetscReal :: rtol
        PetscReal :: nrm
        Vec :: f,res
        KSP :: ksp
        PC :: pc

        call PetscPrintf(comm, "\n   KSPs :: Solve\n\n", ierr)
        ! Set parameters
        one = 1.0d0
        ine = -1.0d0
        rtol = 1e-8
        count = 0
        ! Initialize Vecs
        call VecDuplicate(x,f,ierr)
        call VecZeroEntries(f,ierr)
        call VecDuplicate(x,res,ierr)
        call VecZeroEntries(res,ierr)
        ! Create & set KSP
        call KSPCreate(comm,ksp,ierr)
        call KSPSetOperators(ksp,mat,mat,ierr)
        call KSPSetType(ksp,KSPFGMRES,ierr)
        call KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization,ierr)
        call KSPGMRESSetRestart(ksp,40,ierr)
        call KSPSetTolerances(ksp,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
        ! Get & Set PC
        call KspGetPC(ksp,pc,ierr)
        call PCSetType(pc,PCASM,ierr)
        call PCASMSetOverlap(pc,10,ierr)
        call PCFactorSetUseInPlace(pc,PETSC_TRUE,ierr)
        call PCSetFromOptions(pc,ierr)
        call PCSetUp(pc,ierr)
        call KSPSetFromOptions(ksp,ierr)
        call KSPSetUp(ksp,ierr)
        ! Iteration loops

        select case (type)

            case(NEWTON_LIKE)

                call PetscPrintf(comm,"    using *Newtom-Like* method",ierr)

                do while(.True.)
                    ! Get rhs
                    call fx_rhs(x,f)
                    ! Solve
                    call KSPSolve(ksp,f,res,ierr)
                    ! Get residual
                    call VecNorm(res,NORM_INFINITY,nrm,ierr)
                    ! Print residual
                    write(str_count,"(I5)") count
                    write(str_norm,"(ES20.12)") nrm
                    call PetscPrintf(comm," "//str_count//" < residual i-Norm > "//str_norm//"\n",ierr)
                    ! Get solution
                    call VecAXPY(x,one,res,ierr)
                    ! If iterated too much times
                    if(count>50)then
                        call PetscPrintf(comm,"   Maximum number of iterations reached.\n",ierr)
                        exit
                    endif
                    ! If converged
                    if(nrm<1e-5)then
                        call PetscPrintf(comm,"   Converged.\n",ierr)
                        exit
                    endif
                    ! Count
                    count = count + 1
                enddo

            case(TIME_DISCRETE)

                call PetscPrintf(comm,"    using *Time-Discrete* method",ierr)

                do while(.True.)
                    ! Get rhs
                    call fx_rhs(x,f)
                    ! Solve
                    call KSPSolve(ksp,f,res,ierr)
                    ! Get residual
                    call VecAYPX(x,ine,res,ierr)
                    call VecNorm(x,NORM_INFINITY,nrm,ierr)
                    ! Print residual
                    write(str_count,"(I5)") count
                    write(str_norm,"(ES20.12)") nrm
                    call PetscPrintf(comm," "//str_count//" < residual i-Norm > "//str_norm//"\n",ierr)
                    ! Get solution
                    call VecCopy(res,x,ierr)
                    ! If iterated too much times
                    if(count>50)then
                        call PetscPrintf(comm,"   Maximum number of iterations reached.\n",ierr)
                        exit
                    endif
                    ! If converged
                    if(nrm<1e-5)then
                        call PetscPrintf(comm,"   Converged.\n",ierr)
                        exit
                    endif
                    ! Count
                    count = count + 1
                enddo

        end select

        ! Destory variables
        call KSPDestroy(ksp,ierr)
        call VecDestroy(f,ierr)
        call VecDestroy(res,ierr)
        call cleanup()

    end subroutine solve_ksps

end module mod_solving
