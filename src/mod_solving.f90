#include <slepc/finclude/slepc.h>

module mod_solving
! ----------------------------------------------------
!
!  这个模块是主工作流。
!
!       call dstream(comm) 数据流。
!
!           1.call ksp_equation(comm,mat,x,r) 迭代格式一：标准线性求解器 Ax=b
!
!               1).call solve_ksp_mf(comm,mat,x,r,level) 设置KSP免矩阵求解过程。
!
!               2).call solve_ksp(comm,mat,x,r,level) 设置KSP显式矩阵求解过程。
!
!           2.call snes_equation(comm,mat,x,fx,r) 迭代格式二：借用SNES模块 Jac x = - F(x)
!
!               1).call solve_snes(comm,jac,x,fx,r) 设置SNES求解过程。
!
!               2).call snes_converged_test(snes,it,xnorm,snorm,fnorm,reason,dummy,ierr) SNES收敛判断函数。
!
!           3.call newt_equation(comm,mat,x,ksps_fx) 迭代格式三：借用KSP模块 A x = - (Ax-b)
!
!               call solve_newt(comm,mat,x,fx_rhs) 设置KSPs求解过程。
!
! ----------------------------------------------------
    use mod_parameters
    use mod_metrics
    use mod_forming
    use mod_iterate
    use mod_points
    use petsc
    public :: dstream
    private
    PetscErrorCode :: ierr

    contains

    subroutine dstream(comm)
        implicit none
        PetscInt,intent(in) :: comm

        call PetscPrintf(comm, "\n -----------------------------------\n", ierr)
        call PetscPrintf(comm, "               DStream            \n", ierr)
        call PetscPrintf(comm, "\n   「 Data Preparation 」\n", ierr)
        call metric_coefficient(comm)
        call partial_derivatives(comm)
        select case (solver_mode)
            case('ksp')
                call ksp_equation(comm,levels)
            case('snes')
                call snes_equation(comm,snes_rhs_fx_4ord)
            case('newt')
                call newt_equation(comm,fx_rhs_Ax_4ord)
            case('newt_sub')
                call newtsub_equation(comm,fx_rhs_Ax_4ord)
        end select
    end subroutine dstream

    subroutine ksp_equation(comm,level)
        implicit none
        integer,intent(in) :: level
        integer,intent(in) :: comm

        call PetscPrintf(comm,"\n   「 Forming Matrix 」\n",ierr)
        call init_mat_from_da(comm,meshDA,whale)
        call form_global_mat_4ord(whale)
        call duplicate_vec(turtle,RHS)
        call linear_rhs(comm,RHS)
        call solve_ksp(comm,whale,turtle,RHS,level)
    end subroutine ksp_equation

    subroutine snes_equation(comm,fx)
        implicit none
        integer, intent(in) :: comm
        external :: fx

        call PetscPrintf(comm,"\n   「 Forming Jacobi Matrix 」\n",ierr)
        call init_mat_from_da(comm,meshDA,whale)
        call form_global_mat_2ord(whale)
        call duplicate_vec(turtle,RHS)
        call linear_rhs(comm,RHS)
        call solve_snes(comm,whale,turtle,fx,RHS)
    end subroutine snes_equation

    subroutine newt_equation(comm,fx_rhs)
        implicit none
        integer,intent(in) :: comm
        external :: fx_rhs

        call PetscPrintf(comm,"\n   「 Forming Jacobi Matrix 」\n",ierr)
        call set_medDA(comm,med1DA)
        call init_mat_from_da(comm,med1DA,whale)
        call form_global_mat_2ord(whale)
        call solve_newt(comm,whale,turtle,fx_rhs)
    end subroutine newt_equation

    subroutine newtsub_equation(comm,fx_rhs)
        implicit none
        integer,intent(in) :: comm
        external :: fx_rhs

        call PetscPrintf(comm,"\n   「 Forming sub Jacobi Matrix 」\n",ierr)
        call set_subDA(subDA)
        call init_sub_vecs()
        call init_mat_from_da(comm,subDA,whale)
        call form_sub_mat_2_ord(whale)
        call solve_newtsub(comm,whale,turtle,fx_rhs)
    end subroutine newtsub_equation

    ! -----------------------------------------------------------------------------------------------------
    !   迭代格式一：标准线性求解器 Ax=b

    subroutine solve_ksp(comm,mat,x,r,level)
        implicit none
        integer,intent(in) :: level
        PetscInt,intent(in) :: comm
        Mat,intent(in) :: mat
        Vec,intent(in) :: r
        DM,pointer :: da_list(:)
        Vec,intent(inout) :: x
        real(8) :: rtol
        integer :: i
        KSP :: ksp
        Mat :: Ref
        DM :: fDA
        PC :: pc

        ! Set parameter
        rtol = 1e-8
        if(level==0)then
            call PetscPrintf(comm, "\n   「 Krylov Subspace Method 」\n\n", ierr)
            ! Create & set KSP
            call KSPCreate(comm,ksp,ierr)
            call KSPSetOperators(ksp,mat,mat,ierr)
            call KSPSetType(ksp,KSPFGMRES,ierr)
            call KSPSetInitialGuessNonzero(ksp,ex_ini_gus_flg,ierr)
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
            call KSPSetFromOptions(ksp,ierr)
            call KSPSetType(ksp,KSPFGMRES,ierr)
            call KSPGetPC(ksp,pc,ierr)
            call PCSetType(pc,PCMG,ierr)
            call PCSetDM(pc,meshDA,ierr)

            ! ### 设置网格层数
            call PCMGSetLevels(pc,level,PETSC_NULL_MPI_COMM,ierr)
            call PCMGSetGalerkin(pc,PC_MG_GALERKIN_PMAT,ierr)
            call PCMGSetType(pc,PC_MG_MULTIPLICATIVE,ierr)
            ! call PCMGSetCycleType(pc,PC_MG_CYCLE_W,ierr)
            ! ### 设置网格粗化
            allocate(da_list(level))
            if(lns_mode==2) call DMDASetRefinementFactor(meshDA,PETSC_DECIDE,PETSC_DECIDE,1,ierr)
            call DMRefine(meshDA,comm,fDA,ierr)
            call DMCoarsen(fDA,comm,da_list(1),ierr)
            do i=1,level-1
                call DMCoarsen(da_list(i),comm,da_list(i+1),ierr)
            enddo

            do i=1,level-1
                call DMCreateInterpolation(da_list(i+1),da_list(i),Ref,PETSC_NULL_VEC,ierr)
                ! call PCMGSetInterpolation(pc,i,Ref,ierr)
                ! call MatPtAP(mat,Ref,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,INTP,ierr)
                ! call MatRARt(mat,Ref,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,INTP,ierr)
                ! call PCMGSetOperators(pc,i,INTP,INTP,ierr)
                call MatDestroy(Ref,ierr)
            enddo

            deallocate(da_list)
            ! ### 设置光滑子
            ! ### 设置每层迭代矩阵
            ! call KSPSetUp(ksp,ierr)
            ! Solve KSP
            call KSPSolve(ksp,r,x,ierr)
            ! View KSP
            call PetscPrintf(comm, "\n", ierr)
            ! call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
        endif
        ! Destory variables
        call VecDestroy(r,ierr)
        call KSPDestroy(ksp,ierr)
        call MatDestroy(mat,ierr)
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
        SNES :: snes
        KSP :: ksp
        Vec :: t
        PC :: pc

        call PetscPrintf(comm, "\n    「 Nonlinear Solver 」\n\n", ierr)

        ! 初始化向量
        call VecDuplicate(x,t,ierr)
        call VecZeroEntries(t,ierr)
        ! 设置SNES
        call SNESCreate(comm,snes,ierr)
        call SNESSetConvergenceTest(snes,snes_converged_test,0,PETSC_NULL_FUNCTION,ierr)
        call SNESSetFunction(snes,t,fx,0,ierr)
        call SNESSetJacobian(snes,jac,jac,PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER,ierr)
        ! 设置KSP
        call SNESGetKSP(snes,ksp,ierr)
        call KSPSetType(ksp,KSPFGMRES,ierr)
        call KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization,ierr)
        call KSPGMRESSetRestart(ksp,40,ierr)
        ! 设置PC
        call KSPGetPC(ksp,pc,ierr)
        call PCSetType(pc,PCASM,ierr)
        call PCASMSetOverlap(pc,10,ierr)
        call PCFactorSetUseInPlace(pc,PETSC_TRUE,ierr)
        call PCSetFromOptions(pc,ierr)
        call KSPSetFromOptions(ksp,ierr)
        call SNESSetFromOptions(snes,ierr)
        call SNESSetUp(snes,ierr)
        ! 求解
        call SNESSolve(snes,r,x,ierr)
        ! 输出收敛状况
        call SNESGetConvergedReasonString(snes,reason,ierr)
        call PetscPrintf(comm, "\n", ierr)
        call PetscPrintf(comm, "The Reason of Converged is : "//reason//"\n", ierr)
        ! 查看SNES
        call SNESView(snes,PETSC_VIEWER_STDOUT_WORLD,ierr)
        ! 清理
        call VecDestroy(t,ierr)
        call SNESDestroy(snes,ierr)
        call VecDestroy(r,ierr)
        call MatDestroy(jac,ierr)
        call cleanup()

    end subroutine solve_snes

    ! -----------------------------------------------------------------------------------------------------
    !   迭代格式三：借用KSP模块 Newton-Like

    subroutine solve_newt(comm,mat,x,fx_rhs)
        implicit none
        character(len=20) :: str_norm
        character(len=6) :: str_count
        integer,intent(in) :: comm
        Vec,intent(inout) :: x
        PetscScalar :: one,ine
        PetscErrorCode :: ierr
        Mat,intent(in) :: mat
        external :: fx_rhs
        PetscInt :: count
        PetscReal :: nrm
        Vec :: f,res
        KSP :: ksp
        PC :: pc

        call PetscPrintf(comm, "\n   「 Newton-Like Method 」\n\n", ierr)
        ! Set parameters
        one = 1.0d0
        ine = -1.0d0
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
        ! call PetscPrintf(comm,"    using *Newtom-Like* method\n",ierr)

        do while(.True.)
            ! Count
            count = count + 1
            ! Update Boundary Conditions
            call push_bc(comm,x)
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
            if(count>5000)then
                call PetscPrintf(comm,"\n   < Maximum number of iterations reached. >\n",ierr)
                exit
            endif
            ! If converged
            if(nrm<1e-5)then
                call PetscPrintf(comm,"\n   < Converged. >\n",ierr)
                exit
            endif
        enddo

        ! Destory variables
        call KSPDestroy(ksp,ierr)
        call VecDestroy(f,ierr)
        call VecDestroy(res,ierr)
        call MatDestroy(mat,ierr)
        call cleanup()

    end subroutine solve_newt

    ! -----------------------------------------------------------------------------------------------------
    !   迭代格式四：分块求解 Newton-Like

    subroutine solve_newtsub(comm,mat,x,fx_rhs)
        implicit none
        character(len=20) :: str_norm
        character(len=6) :: str_count
        integer,intent(in) :: comm
        Vec,intent(inout) :: x
        PetscErrorCode :: ierr
        Mat,intent(in) :: mat
        Vec :: subres,subf
        PetscScalar :: one
        external :: fx_rhs
        PetscInt :: count
        PetscReal :: nrm
        Vec :: f,res
        KSP :: ksp
        PC :: pc

        call PetscPrintf(comm, "\n   「 Newton-Like Subs Method 」\n\n", ierr)

        ! 初始化变量
        one = 1.0d0
        count = 0
        call VecDuplicate(x,f,ierr)
        call VecZeroEntries(f,ierr)
        call VecDuplicate(x,res,ierr)
        call VecZeroEntries(res,ierr)
        call VecDuplicate(subx,subf,ierr)
        call VecZeroEntries(subf,ierr)
        call VecDuplicate(subx,subres,ierr)
        call VecZeroEntries(subres,ierr)

        ! 设置KSP
        call KSPCreate(PETSC_COMM_SELF,ksp,ierr)
        call KSPSetOperators(ksp,mat,mat,ierr)
        call KSPSetType(ksp,KSPFGMRES,ierr)
        call KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization,ierr)
        call KSPGMRESSetRestart(ksp,40,ierr)
        call KspGetPC(ksp,pc,ierr)
        call PCSetType(pc,PCILU,ierr)
        call PCSetFromOptions(pc,ierr)
        call PCSetUp(pc,ierr)
        call KSPSetFromOptions(ksp,ierr)
        call KSPSetUp(ksp,ierr)

        ! 设置迭代过程
        do while(.True.)
            ! 迭代计数
            count = count + 1
            ! 载入边界条件
            call push_bc(comm,x)
            ! 计算右端项
            call fx_rhs(x,f)
            ! 分发右端项
            call get_subf(f,subf)
            ! 求解
            call KSPSolve(ksp,subf,subres,ierr)
            ! 拼接残差
            call merge_res(subres,res)
            ! 计算残差范数
            call VecNorm(res,NORM_INFINITY,nrm,ierr)
            ! 输出残差
            write(str_count,"(I5)") count
            write(str_norm,"(ES20.12)") nrm
            call PetscPrintf(comm," "//str_count//" < residual i-Norm > "//str_norm//"\n",ierr)
            ! 刷新近似解
            call VecAXPY(x,one,res,ierr)
            ! 如果达到最大迭代数
            if(count>5000)then
                call PetscPrintf(comm,"\n   < Maximum number of iterations reached. >\n",ierr)
                exit
            endif
            ! 如果收敛
            if(nrm<1e-5)then
                call PetscPrintf(comm,"\n   < Converged. >\n",ierr)
                exit
            endif
        enddo
        call DMRestoreGlobalVector(subDA,subx,ierr)
        call DMDestroy(subDA,ierr)
        call KSPDestroy(ksp,ierr)
        call MatDestroy(mat,ierr)
        call VecDestroy(subf,ierr)
        call VecDestroy(subres,ierr)
        call VecDestroy(f,ierr)
        call VecDestroy(res,ierr)
        call cleanup()
    end subroutine solve_newtsub

    subroutine cleanup()

        use mod_parameters
        deallocate(bf)
        deallocate(inlet)
        deallocate(xi_x,xi_y,xi_z)
        deallocate(eta_x,eta_y,eta_z)
        deallocate(phi_x,phi_y,phi_z)
        deallocate(xi_xx,xi_yy,xi_zz)
        deallocate(eta_xx,eta_yy,eta_zz)
        deallocate(phi_xx,phi_yy,phi_zz)
        deallocate(xi_xy,xi_xz,xi_yz)
        deallocate(eta_xy,eta_yz,eta_xz)
        deallocate(phi_xy,phi_yz,phi_xz)
        call DMDestroy(DA,ierr)

    end subroutine cleanup

end module mod_solving
