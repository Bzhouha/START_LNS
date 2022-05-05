!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_forming
! ------------------------------------------------------------------
!
!  这个模块包含调用xxxSolve前的所有准备步骤，包括矩阵、右端量的生成函数，清理函数等。
!
!   这个模块是混乱的。
!
!   在线性求解系统 KSP 框架下，生成左端矩阵；
!   在非线性求解系统 SNES 框架下，生成雅各比矩阵和右端项函数。
!
!       1.call set_matfree_mat(comm,ierr) 免矩阵形式的矩阵生成函数
!
!           call mat_mult_4ord(A, X, F, ierr) 免矩阵需要的矩阵向量乘法函数
!
!       2.call get_mat_from_dm(comm,ierr) 显式矩阵形式的矩阵生成函数
!
!           call form_mat_4ord(ierr) 填充数据函数
!
!               (1).call mat_set_boundary_conditions(i,j,k) 边界部分的填充数据函数
!
!               (2).call mat_insert_values_4_ord(i,j,k) 内部部分的填充数据函数
!
!       3.call set_rhs(comm,ierr) 设置右边量，即边界。
!
!       4.call cleanup() 清理
!
!       5.call shark_coming(comm,ierr) 一阶精度分裂的矩阵
!
!           call form_mat_2ord(ierr) 生成一阶精度分裂的矩阵
!
!       6.call snes_fx1o(snes,x,f,null_int,ierr) 一阶精度右端项函数
!
!       7.call snes_rhs_fx_4ord(snes,x,f,null_int,ierr) 四阶精度右端项函数
!
! ------------------------------------------------------------------
    use penf, only: R_P
    use mod_cubes
    use petsc
    implicit none
    public :: set_matfree_mat
    public :: mat_mult_4ord
    public :: init_mat_from_da
    public :: form_mat_2ord
    public :: form_mat_4ord
    public :: duplicate_vec
    public :: set_rhs
    public :: snes_rhs_fx
    public :: snes_rhs_fx_4ord
    public :: snes_converged_test
    public :: push_bc
    public :: rhs_fx_Ax
    public :: ksps_rhs_fx_b_Gtx
    public :: form_sub_mat_2_ord
    public :: set_subDA
    public :: init_sub_vecs
    public :: get_subf
    public :: merge_res
    public :: cleanup
    private
    type(lns_OP_point_type) :: Jor
    real(R_P), parameter, dimension(-2:2,-2:2) :: FDM_1nd_4ORD_CENTER=reshape( [&
        0.0d0       ,        0.0d0, -3.0d0/2.0d0, 2.0d0      ,  -1.0d0/2.0d0, &
        0.0d0       , -1.0d0/3.0d0, -1.0d0/2.0d0, 1.0d0      ,  -1.0d0/6.0d0, &
        1.0d0/12.0d0, -2.0d0/3.0d0,        0.0d0, 2.0d0/3.0d0, -1.0d0/12.0d0, &
        1.0d0/6.0d0 ,       -1.0d0,  1.0d0/2.0d0, 1.0d0/3.0d0,         0.0d0, &
        1.0d0/2.0d0 ,       -2.0d0,  3.0d0/2.0d0, 0.0d0      ,         0.0d0  &
        ], [5,5])
    real(R_P), parameter, dimension(-2:2,-2:2) :: FDM_2nd_4ORD_CENTER=reshape( [&
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0, &
        0.0d0        ,       1.0d0,        -2.0d0,       1.0d0,         0.0d0, &
        -1.0d0/12.0d0, 4.0d0/3.0d0, -15.0d0/6.0d0, 4.0d0/3.0d0, -1.0d0/12.0d0, &
        0.0d0        ,       1.0d0,        -2.0d0,       1.0d0,         0.0d0, &
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0  &
        ], [5,5])
    real(R_P), parameter, dimension(-2:2,-2:2) :: FDM_1nd_4ORD_Backward=reshape( [&
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0, &
        0.0d0        ,      -1.0d0,         1.0d0,       0.0d0,         0.0d0, &
        1.0d0/2.0d0  ,      -2.0d0,   3.0d0/2.0d0,       0.0d0,         0.0d0, &
        0.0d0        ,      -1.0d0,         1.0d0,       0.0d0,         0.0d0, &
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0  &
        ],[5,5])
    real(R_P), parameter, dimension(-2:2,-2:2) :: FDM_1nd_4ORD_Forward=reshape( [&
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0, &
        0.0d0        ,       0.0d0,        -1.0d0,       1.0d0,         0.0d0, &
        0.0d0        ,       0.0d0,  -3.0d0/2.0d0,       2.0d0,  -1.0d0/2.0d0, &
        0.0d0        ,       0.0d0,        -1.0d0,       1.0d0,         0.0d0, &
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0  &
        ],[5,5])
    real(R_P), parameter, dimension(-1:1,-1:1) :: FDM_1nd_1ORD_Backward=reshape( [&
           0.0d0,           -1.0d0,                1.0d0, &
          -1.0d0,            1.0d0,                0.0d0, &
          -1.0d0,            1.0d0,                0.0d0  &
        ],[3,3])
    real(R_P), parameter, dimension(-1:1,-1:1) :: FDM_1nd_1ORD_Forward=reshape( [&
           0.0d0,           -1.0d0,                1.0d0, &
           0.0d0,           -1.0d0,                1.0d0, &
          -1.0d0,            1.0d0,                0.0d0  &
        ],[3,3])
    real(R_P), parameter, dimension(-1:1,-1:1) :: FDM_2nd_2ORD_CENTER=reshape( [&
           0.0d0,            0.0d0,                0.0d0, &
           1.0d0,           -2.0d0,                1.0d0, &
           0.0d0,            0.0d0,                0.0d0  &
        ], [3,3])
    real(R_P), parameter, dimension(-1:1,-1:1) :: FDM_1nd_2ORD_CENTER=reshape( [&
           0.0d0,          -1.0d0,                 1.0d0, &
    -1.0d0/2.0d0,           0.0d0,           1.0d0/2.0d0, &
          -1.0d0,           1.0d0,                 0.0d0  &
        ], [3,3])
    real(R_P), parameter :: delta_i(-2:2)=[0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]
    real(R_P), parameter :: delta_j(-2:2)=[0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]
    real(R_P), parameter :: delta_k(-2:2)=[0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]
    contains

    ! -----------------------------------------------------------------------------------------------------
    !   迭代格式一：标准线性求解器 Ax=b

    subroutine set_matfree_mat(comm,mat)
        use mod_parameters,only:turtle
        implicit none
        PetscInt,intent(in) :: comm
        Mat,intent(inout) :: mat
        PetscErrorCode :: ierr
        PetscInt :: ls

        call VecGetLocalSize(turtle,ls,ierr)
        call MatCreateShell(comm,ls,ls,PETSC_DETERMINE,PETSC_DETERMINE,PETSC_NULL_INTEGER,mat,ierr)
        call MatShellSetOperation(mat,MATOP_MULT,mat_mult_4ord,ierr)
        call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)

    end subroutine set_matfree_mat

    subroutine mat_mult_4ord(A,X,F,ierr)

        use mod_parameters,only : meshDA,lns_mode,in,jn,kn,is,ie,js,je,ks,ke
        implicit none
        PetscScalar,pointer :: fr(:,:,:,:),xr(:,:,:,:)
        integer :: ic_index,jc_index,kc_index
        integer :: lib,lie,ljb,lje,lkb,lke
        PetscScalar,dimension(5,16) :: yi
        PetscErrorCode :: ierr
        integer :: li,lj,lk
        integer :: i,j,k
        Vec :: bell
        Vec :: X,F
        Mat :: A

        call DMGetLocalVector(meshDA,bell,ierr)
        call VecZeroEntries(bell,ierr)
        call DMGlobalToLocalBegin(meshDA,X,INSERT_VALUES,bell,ierr)
        call DMGlobalToLocalEnd(meshDA,X,INSERT_VALUES,bell,ierr)
        call DMDAVecGetArrayReadF90(meshDA,bell,xr,ierr)
        call DMDAVecGetArrayF90(meshDA,F,fr,ierr)
        associate ( &
        &   coef_c4  => FDM_1nd_4ORD_CENTER,   &
        &   coef_d4  => FDM_2nd_4ORD_CENTER,   &
        &   coef_c4f => FDM_1nd_4ORD_Forward,  &
        &   coef_c4b => FDM_1nd_4ORD_Backward, &
        &   f   => fr,        x => xr,         &
        &   G   => Jor%G,     D => Jor%D,      &
        &   A   => Jor%A,     B => Jor%B,     C => jor%C,   &
        &   A_p => Jor%A_p, A_m => Jor%A_m, A_v => Jor%A_v, &
        &   B_p => Jor%B_p, B_m => Jor%B_m, B_v => Jor%B_v, &
        &   C_p => Jor%C_p, C_m => Jor%C_m, C_v => Jor%C_v, &
        &   Vxx => Jor%Vxx, Vyy => Jor%Vyy, Vzz => Jor%Vzz, &
        &   Vxy => Jor%Vxy, Vxz => Jor%Vxz, Vyz => Jor%Vyz)
        do k=ks,ke
            do j=js,je
                do i=is,ie
                    if(i==0 .or. j==(jn-1))then
                        f(:,i,j,k)=x(:,i,j,k)
                    elseif(i==(in-1) .and. j/=0 .and. j/=(jn-1))then
                        f(:,i,j,k)=(x(:,i-2,j,k)-4*x(:,i-1,j,k)+3*x(:,i,j,k))/2.0d0
                    elseif(j==0 .and. i/=0)then
                        call Jor%get_transed_cubes(i,j,k)
                        do lj=1,5
                            do li=2,5
                                D(li,lj)=0.0d0
                                B(li,lj)=0.0d0
                            enddo
                        enddo
                        D(2,2)=1.0d0;D(3,3)=1.0d0
                        D(4,4)=1.0d0;D(5,5)=1.0d0
                        f(:,i,j,k)=matmul(D,x(:,i,j,k))+matmul(B,x(:,i,j+1,k))-matmul(B,x(:,i,j,k))
                    else
                        yi=0.0d0
                        call Jor%get_adorned_cubes(i,j,k)
                        if(i==0)then
                            lib=0; lie=2
                            ic_index=-2
                        elseif(i==1)then
                            lib=-1; lie=2
                            ic_index=-1
                        elseif(i==(in-2))then
                            lib=-2; lie=1
                            ic_index=1
                        elseif(i==(in-1))then
                            lib=-2; lie=0
                            ic_index=2
                        else
                            lib=-2; lie=2
                            ic_index=0
                        endif
                        if(j==0)then
                            ljb=0; lje=2
                            jc_index=-2
                        elseif(j==1)then
                            ljb=-1; lje=2
                            jc_index=-1
                        elseif(j==(jn-2))then
                            ljb=-2; lje=1
                            jc_index=1
                        elseif(j==(jn-1))then
                            ljb=-2; lje=0
                            jc_index=2
                        else
                            ljb=-2; lje=2
                            jc_index=0
                        endif
                        select case (lns_mode)
                            case(2)
                                lkb=0; lke=0; kc_index=0
                            case(3)
                                lkb=-2; lke=2; kc_index=0
                        end select

                        yi(:,1) = x(:,i,j,k) ! D

                        do li=lib,lie
                            yi(:,2)  = yi(:,2)  + x(:,i+li,j,k)*coef_c4b(li,ic_index) ! A_p
                            yi(:,3)  = yi(:,3)  + x(:,i+li,j,k)*coef_c4f(li,ic_index) ! A_m
                            yi(:,4)  = yi(:,4)  + x(:,i+li,j,k)*coef_c4(li,ic_index)  ! A_v
                            yi(:,11) = yi(:,11) + x(:,i+li,j,k)*coef_d4(li,ic_index)  ! Vxx
                        enddo

                        do lj=ljb,lje
                            yi(:,5)  = yi(:,5)  + x(:,i,j+lj,k)*coef_c4b(lj,jc_index) ! B_p
                            yi(:,6)  = yi(:,6)  + x(:,i,j+lj,k)*coef_c4f(lj,jc_index) ! B_m
                            yi(:,7)  = yi(:,7)  + x(:,i,j+lj,k)*coef_c4(lj,jc_index)  ! B_v
                            yi(:,12) = yi(:,12) + x(:,i,j+lj,k)*coef_d4(lj,jc_index)  ! Vxx
                        enddo

                        do lk=lkb,lke
                            yi(:,8)  = yi(:,8)  + x(:,i,j,k+lk)*coef_c4b(lk,kc_index) ! C_p
                            yi(:,9)  = yi(:,9)  + x(:,i,j,k+lk)*coef_c4f(lk,kc_index) ! C_m
                            yi(:,10) = yi(:,10) + x(:,i,j,k+lk)*coef_c4(lk,kc_index)  ! C_v
                            yi(:,13) = yi(:,13) + x(:,i,j,k+lk)*coef_d4(lk,kc_index)  ! Vxx
                        enddo

                        do lj=ljb,lje
                            do li=lib,lie
                                yi(:,14)=yi(:,14)+x(:,i+li,j+lj,k)*coef_c4(li,ic_index)*coef_c4(lj,jc_index) ! Vxy
                            enddo
                        enddo

                        do lk=lkb,lke
                            do li=lib,lie
                                yi(:,15)=yi(:,15)+x(:,i+li,j,k+lk)*coef_c4(li,ic_index)*coef_c4(lk,kc_index) ! Vxz
                            enddo
                        enddo

                        do lk=lkb,lke
                            do lj=ljb,lje
                                yi(:,16)=yi(:,16)+x(:,i,j+lj,k+lk)*coef_c4(lj,jc_index)*coef_c4(lk,kc_index) ! Vyz
                            enddo
                        enddo

                        f(:,i,j,k)=matmul(D,yi(:,1))+                                              &
                        &          matmul(A_p,yi(:,2))+matmul(A_m,yi(:,3))+matmul(A_v,yi(:,4))+    &
                        &          matmul(B_p,yi(:,5))+matmul(B_m,yi(:,6))+matmul(B_v,yi(:,7))+    &
                        &          matmul(C_p,yi(:,8))+matmul(C_m,yi(:,9))+matmul(C_v,yi(:,10))-   &
                        &          matmul(Vxx,yi(:,11))-matmul(Vyy,yi(:,12))-matmul(Vzz,yi(:,13))- &
                        &          matmul(Vxy,yi(:,14))-matmul(Vxz,yi(:,15))-matmul(Vyz,yi(:,16))
                    endif
                enddo
            enddo
        enddo
        end associate
        call DMDAVecRestoreArrayReadF90(meshDA,bell,xr,ierr)
        call DMDAVecRestoreArrayF90(meshDA,F,fr,ierr)
        call DMRestoreLocalVector(meshDA,bell,ierr)

    end subroutine mat_mult_4ord

    subroutine init_mat_from_da(comm,da,mat)

        implicit none
        PetscInt,intent(in) :: comm
        PetscErrorCode :: ierr
        Mat :: mat
        DM :: da

        call DMCreateMatrix(da,mat,ierr)
        call MatZeroEntries(mat,ierr)

    end subroutine init_mat_from_da

    subroutine form_mat_4ord(mat)

        use mod_parameters,only : in,jn,kn,is,ie,js,je,ks,ke
        implicit none
        Mat,intent(inout) :: mat
        PetscErrorCode :: ierr
        integer :: i,j,k

        do k=ks,ke
            do j=js,je
                do i=is,ie
                    if(i==0 .or. i==(in-1) .or. j==0 .or. j==(jn-1))then
                        call mat_set_boundary_conditions(mat,i,j,k)
                    else
                        call mat_insert_values_4_ord(mat,i,j,k)
                    endif
                enddo
            enddo
        enddo

        call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)

    end subroutine form_mat_4ord

    subroutine mat_set_boundary_conditions(mat,i,j,k)

        use mod_parameters,only : in,jn,kn
        implicit none
        PetscScalar :: box(5,5),trans(5,5)
        MatStencil :: idxm(4,1),idxn(4,1)
        integer :: ic_index, jc_index
        integer :: lib, lie, ljb, lje
        integer,intent(in) :: i,j,k
        Mat,intent(inout) :: mat
        PetscErrorCode :: ierr
        integer :: li,lj

        associate( &
            coef_c4f=>FDM_1nd_4ORD_Forward, &
            coef_c4b=>FDM_1nd_4ORD_Backward)
        if(i==0 .or. j==(jn-1))then
            box=0.0d0
            box(1,1)=1.0d0;box(2,2)=1.0d0;box(3,3)=1.0d0;box(4,4)=1.0d0;box(5,5)=1.0d0
            idxm(MatStencil_i, 1)=i; idxm(MatStencil_j, 1)=j; idxm(MatStencil_k, 1)=k
            idxn(MatStencil_i, 1)=i; idxn(MatStencil_j, 1)=j; idxn(MatStencil_k, 1)=k
            call MatSetValuesBlockedStencil(mat, 1, idxm, 1, idxn, box, INSERT_VALUES, ierr)
        elseif(j==0 .and. i/=0)then
            ljb=0;lje=1
            jc_index=1
            idxm(MatStencil_i, 1)=i; idxm(MatStencil_j, 1)=j; idxm(MatStencil_k, 1)=k
            call Jor%get_transed_cubes(i,j,k)
            do lj=1,5
                do li=2,5
                    Jor%D(li,lj)=0.0d0
                    Jor%B(li,lj)=0.0d0
                enddo
            enddo
            Jor%D(2,2)=1.0d0;Jor%D(3,3)=1.0d0
            Jor%D(4,4)=1.0d0;Jor%D(5,5)=1.0d0
            do lj=ljb,lje
                idxn(MatStencil_i, 1)=i
                idxn(MatStencil_j, 1)=j+lj
                idxn(MatStencil_k, 1)=k
                box=0.0d0;trans=0.0d0
                box=delta_j(lj)*Jor%D+coef_c4f(lj,jc_index)*Jor%B
                trans=transpose(box)
                call MatSetValuesBlockedStencil(mat, 1, idxm, 1, idxn, trans, INSERT_VALUES, ierr)
            enddo
        elseif(i==(in-1) .and. j/=0 .and. j/=(jn-1))then
            lib=-2;lie=0
            ic_index=0
            idxm(MatStencil_i, 1)=i; idxm(MatStencil_j, 1)=j; idxm(MatStencil_k, 1)=k
            do li=lib,lie
                idxn(MatStencil_i, 1)=i+li
                idxn(MatStencil_j, 1)=j
                idxn(MatStencil_k, 1)=k
                box=0.0d0
                box(1,1)=1.0d0;box(2,2)=1.0d0;box(3,3)=1.0d0;box(4,4)=1.0d0;box(5,5)=1.0d0
                box=coef_c4b(li,ic_index)*box
                call MatSetValuesBlockedStencil(mat, 1, idxm, 1, idxn, box, INSERT_VALUES, ierr)
            enddo
        endif
        end associate

    end subroutine mat_set_boundary_conditions

    subroutine mat_bc_set_I(mat,i,j,k)

        implicit none
        MatStencil :: idxm(4,1),idxn(4,1)
        integer,intent(in) :: i,j,k
        Mat,intent(inout) :: mat
        PetscScalar :: box(5,5)
        PetscErrorCode :: ierr

        box=0.0d0
        box(1,1)=1.0d0;box(2,2)=1.0d0;box(3,3)=1.0d0;box(4,4)=1.0d0;box(5,5)=1.0d0
        idxm(MatStencil_i, 1)=i; idxm(MatStencil_j, 1)=j; idxm(MatStencil_k, 1)=k
        idxn(MatStencil_i, 1)=i; idxn(MatStencil_j, 1)=j; idxn(MatStencil_k, 1)=k
        call MatSetValuesBlockedStencil(mat, 1, idxm, 1, idxn, box, INSERT_VALUES, ierr)

    end subroutine mat_bc_set_I

    subroutine mat_insert_values_4_ord(mat,i,j,k)

        use mod_parameters,only : lns_mode,in,jn,kn
        implicit none
        integer :: ic_index, jc_index, kc_index
        integer :: lib, lie, ljb, lje, lkb, lke
        PetscScalar :: box(5,5),trans(5,5)
        MatStencil :: idxm(4,1),idxn(4,1)
        integer,intent(in) :: i,j,k
        Mat,intent(inout) :: mat
        PetscErrorCode :: ierr
        integer :: li, lj, lk

        call Jor%get_adorned_cubes(i,j,k)
        if(i==0)then
            lib=0; lie=2
            ic_index=-2
        elseif(i==1)then
            lib=-1; lie=2
            ic_index=-1
        elseif(i==(in-2))then
            lib=-2; lie=1
            ic_index=1
        elseif(i==(in-1))then
            lib=-2; lie=0
            ic_index=2
        else
            lib=-2; lie=2
            ic_index=0
        endif
        if(j==0)then
            ljb=0; lje=2
            jc_index=-2
        elseif(j==1)then
            ljb=-1; lje=2
            jc_index=-1
        elseif(j==(jn-2))then
            ljb=-2; lje=1
            jc_index=1
        elseif(j==(jn-1))then
            ljb=-2; lje=0
            jc_index=2
        else
            ljb=-2; lje=2
            jc_index=0
        endif
        select case (lns_mode)
            case(2)
                lkb=0; lke=0; kc_index=0
            case(3)
                lkb=-2; lke=2; kc_index=0
        end select
        idxm=0;
        idxm(MatStencil_i, 1)=i; idxm(MatStencil_j, 1)=j; idxm(MatStencil_k, 1)=k
        associate( &
        &   coef_c4  => FDM_1nd_4ORD_CENTER,   &
        &   coef_d4  => FDM_2nd_4ORD_CENTER,   &
        &   coef_c4f => FDM_1nd_4ORD_Forward,  &
        &   coef_c4b => FDM_1nd_4ORD_Backward, &
        &   G   => Jor%G,   D   => Jor%D,      &
        &   A_p => Jor%A_p, A_m => Jor%A_m, A_v => Jor%A_v, &
        &   B_p => Jor%B_p, B_m => Jor%B_m, B_v => Jor%B_v, &
        &   C_p => Jor%C_p, C_m => Jor%C_m, C_v => Jor%C_v, &
        &   Vxx => Jor%Vxx, Vyy => Jor%Vyy, Vzz => Jor%Vzz, &
        &   Vxy => Jor%Vxy, Vxz => Jor%Vxz, Vyz => Jor%Vyz)
        do lk = lkb, lke
            do lj = ljb, lje
                do li = lib, lie
                    idxn=0;box=0.0d0;trans=0.0d0
                    idxn(MatStencil_i, 1)=i+li
                    idxn(MatStencil_j, 1)=j+lj
                    idxn(MatStencil_k, 1)=k+lk
                    !if(idxn(MatStencil_k, 1)>(kn-1)) idxn(MatStencil_k, 1)=idxn(MatStencil_k, 1)-kn
                    !if(idxn(MatStencil_k, 1)<0)  idxn(MatStencil_k, 1)=idxn(MatStencil_k, 1)+kn  !! kn为展向一个周期的点数
                    box=delta_i(li)*delta_j(lj)*delta_k(lk)*D+                     &
                    &   delta_j(lj)*delta_k(lk)*A_v*coef_c4(li,ic_index)+          &
                    &   delta_j(lj)*delta_k(lk)*A_m*coef_c4f(li,ic_index)+         &
                    &   delta_j(lj)*delta_k(lk)*A_p*coef_c4b(li,ic_index)+         &
                    &   delta_i(li)*delta_k(lk)*B_v*coef_c4(lj,jc_index)+          &
                    &   delta_i(li)*delta_k(lk)*B_m*coef_c4f(lj,jc_index)+         &
                    &   delta_i(li)*delta_k(lk)*B_p*coef_c4b(lj,jc_index)+         &
                    &   delta_i(li)*delta_j(lj)*C_v*coef_c4(lk,kc_index)+          &
                    &   delta_i(li)*delta_j(lj)*C_m*coef_c4f(lk,kc_index)+         &
                    &   delta_i(li)*delta_j(lj)*C_p*coef_c4b(lk,kc_index)-         &
                    &   delta_j(lj)*delta_k(lk)*Vxx*coef_d4(li,ic_index)-          &
                    &   delta_i(li)*delta_k(lk)*Vyy*coef_d4(lj,jc_index)-          &
                    &   delta_i(li)*delta_j(lj)*Vzz*coef_d4(lk,kc_index)-          &
                    &   delta_k(lk)*Vxy*coef_c4(li,ic_index)*coef_c4(lj,jc_index)- &
                    &   delta_j(lj)*Vxz*coef_c4(li,ic_index)*coef_c4(lk,kc_index)- &
                    &   delta_i(li)*Vyz*coef_c4(lj,jc_index)*coef_c4(lk,kc_index)
                    trans=transpose(box)
                    call MatSetValuesBlockedStencil(mat, 1, idxm, 1, idxn, trans, INSERT_VALUES, ierr)
                end do
            end do
        end do
        end associate

    end subroutine mat_insert_values_4_ord

    subroutine duplicate_vec(x,r)
        implicit none
        Vec,intent(inout) :: r
        PetscErrorCode :: ierr
        Vec,intent(in) :: x

        call VecDuplicate(x,r,ierr)
        call VecZeroEntries(r,ierr)
    end subroutine duplicate_vec

    subroutine set_rhs(comm,r)

        use mod_parameters,only : disturb,meshDA,is,ie,js,je,ks,ke
        implicit none
        PetscScalar,pointer :: r_array(:,:,:,:)
        PetscInt,intent(in) :: comm
        Vec,intent(inout) :: r
        PetscErrorCode :: ierr
        integer :: j,k

        if (is==0) then
            call DMDAVecGetArrayF90(meshDA,r,r_array,ierr)
                do k=ks,ke
                    do j=js,je
                        r_array(:,0,j,k)=disturb(:,j,k)
                    enddo
                enddo
            call DMDAVecRestoreArrayF90(meshDA,r,r_array,ierr)
        endif
        call MPI_Barrier(comm,ierr)

    end subroutine set_rhs

    subroutine form_mat_2ord(mat)

        use mod_parameters,only : in,jn,kn,is,ie,js,je,ks,ke
        implicit none
        Mat,intent(inout) :: mat
        PetscErrorCode :: ierr
        integer :: i,j,k

        do k=ks,ke
            do j=js,je
                do i=is,ie
                    if(i==0 .or. i==(in-1) .or. j==0 .or. j==(jn-1))then
                        call mat_bc_set_I(mat,i,j,k)
                    else
                        call mat_insert_values_2_ord(mat,i,j,k)
                    endif
                enddo
            enddo
        enddo

        call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)

    end subroutine form_mat_2ord

    subroutine mat_insert_values_2_ord(mat,i,j,k)

        use mod_parameters,only : lns_mode,in,jn,kn
        implicit none
        integer :: ic_index, jc_index, kc_index
        integer :: lib, lie, ljb, lje, lkb, lke
        PetscScalar :: box(5,5),trans(5,5)
        MatStencil :: idxm(4,1),idxn(4,1)
        integer,intent(in) :: i,j,k
        Mat,intent(inout) :: mat
        PetscErrorCode :: ierr
        integer :: li, lj, lk

        call Jor%get_adorned_cubes(i,j,k)
        if(i==0)then
            lib=0; lie=1
            ic_index=-1
        elseif(i==(in-1))then
            lib=-1; lie=0
            ic_index=1
        else
            lib=-1; lie=1
            ic_index=0
        endif
        if(j==0)then
            ljb=0; lje=1
            jc_index=-1
        elseif(j==(jn-1))then
            ljb=-1; lje=0
            jc_index=1
        else
            ljb=-1; lje=1
            jc_index=0
        endif
        select case (lns_mode)
            case(2)
                lkb=0; lke=0; kc_index=0
            case(3)
                lkb=-1; lke=1; kc_index=0
        end select
        idxm=0;
        idxm(MatStencil_i, 1)=i; idxm(MatStencil_j, 1)=j; idxm(MatStencil_k, 1)=k
        associate( &
            &   coef_c1f=>FDM_1nd_1ORD_Forward,  &
            &   coef_c1b=>FDM_1nd_1ORD_Backward, &
            &   coef_d2 =>FDM_2nd_2ORD_CENTER,   &
            &   coef_c2 =>FDM_1nd_2ORD_CENTER,   &
            &   G   => Jor%G,   D   => Jor%D,    &
            &   A_p => Jor%A_p, A_m => Jor%A_m, A_v => Jor%A_v, &
            &   B_p => Jor%B_p, B_m => Jor%B_m, B_v => Jor%B_v, &
            &   C_p => Jor%C_p, C_m => Jor%C_m, C_v => Jor%C_v, &
            &   Vxx => Jor%Vxx, Vyy => Jor%Vyy, Vzz => Jor%Vzz, &
            &   Vxy => Jor%Vxy, Vxz => Jor%Vxz, Vyz => Jor%Vyz)
            do lk = lkb, lke
                do lj = ljb, lje
                    do li = lib, lie
                        idxn=0;box=0.0d0;trans=0.0d0
                        idxn(MatStencil_i, 1)=i+li
                        idxn(MatStencil_j, 1)=j+lj
                        idxn(MatStencil_k, 1)=k+lk
                        box=delta_i(li)*delta_j(lj)*delta_k(lk)*D+                     &
                        &   delta_j(lj)*delta_k(lk)*A_v*coef_c2(li,ic_index)+          &
                        &   delta_j(lj)*delta_k(lk)*A_m*coef_c1f(li,ic_index)+         &
                        &   delta_j(lj)*delta_k(lk)*A_p*coef_c1b(li,ic_index)+         &
                        &   delta_i(li)*delta_k(lk)*B_v*coef_c2(lj,jc_index)+          &
                        &   delta_i(li)*delta_k(lk)*B_m*coef_c1f(lj,jc_index)+         &
                        &   delta_i(li)*delta_k(lk)*B_p*coef_c1b(lj,jc_index)+         &
                        &   delta_i(li)*delta_j(lj)*C_v*coef_c2(lk,kc_index)+          &
                        &   delta_i(li)*delta_j(lj)*C_m*coef_c1f(lk,kc_index)+         &
                        &   delta_i(li)*delta_j(lj)*C_p*coef_c1b(lk,kc_index)-         &
                        &   delta_j(lj)*delta_k(lk)*Vxx*coef_d2(li,ic_index)-          &
                        &   delta_i(li)*delta_k(lk)*Vyy*coef_d2(lj,jc_index)-          &
                        &   delta_i(li)*delta_j(lj)*Vzz*coef_d2(lk,kc_index)-          &
                        &   delta_k(lk)*Vxy*coef_c2(li,ic_index)*coef_c2(lj,jc_index)- &
                        &   delta_j(lj)*Vxz*coef_c2(li,ic_index)*coef_c2(lk,kc_index)- &
                        &   delta_i(li)*Vyz*coef_c2(lj,jc_index)*coef_c2(lk,kc_index)
                        trans=transpose(box)
                        call MatSetValuesBlockedStencil(mat, 1, idxm, 1, idxn, trans, INSERT_VALUES, ierr)
                    enddo
                enddo
            enddo
        end associate

    end subroutine mat_insert_values_2_ord

    subroutine snes_rhs_fx(snes,x,f,null_int,ierr)

        use mod_parameters,only : meshDA,lns_mode,in,jn,kn,is,ie,js,je,ks,ke
        implicit none
        PetscScalar,pointer :: fr(:,:,:,:),xr(:,:,:,:)
        integer :: ic_index,jc_index,kc_index
        integer :: lib,lie,ljb,lje,lkb,lke
        PetscScalar,dimension(5,10) :: yi
        PetscErrorCode :: ierr
        integer :: null_int(*)
        integer :: li,lj,lk
        integer :: i,j,k
        SNES :: snes
        Vec :: bell
        Vec :: x,f

        call DMGetLocalVector(meshDA,bell,ierr)
        call VecZeroEntries(bell,ierr)
        call DMGlobalToLocalBegin(meshDA,X,INSERT_VALUES,bell,ierr)
        call DMGlobalToLocalEnd(meshDA,X,INSERT_VALUES,bell,ierr)
        call DMDAVecGetArrayReadF90(meshDA,bell,xr,ierr)
        call DMDAVecGetArrayF90(meshDA,F,fr,ierr)
        associate ( &
        &   coef_c4 =>FDM_1nd_4ORD_CENTER, &
        &   coef_d4 =>FDM_2nd_4ORD_CENTER, &
        &   f   => fr,    x => xr,         &
        &   G   => Jor%G, D => Jor%D,      &
        &   A   => Jor%A, B => Jor%B, C => jor%C,           &
        &   Vxx => Jor%Vxx, Vyy => Jor%Vyy, Vzz => Jor%Vzz, &
        &   Vxy => Jor%Vxy, Vxz => Jor%Vxz, Vyz => Jor%Vyz)
        do k=ks,ke
            do j=js,je
                do i=is,ie
                    if(i==0)then
                        f(:,i,j,k)=x(:,i,j,k)
                    elseif(j==(jn-1))then
                        f(:,i,j,k)=x(:,i,j,k)
                    elseif(i==(in-1) .and. j/=0 .and. j/=(jn-1))then
                        f(:,i,j,k)=(x(:,i-2,j,k)-4*x(:,i-1,j,k)+3*x(:,i,j,k))/2.0d0
                    elseif(j==0 .and. i/=0)then
                        call Jor%get_transed_cubes(i,j,k)
                        do lj=1,5
                            do li=2,5
                                D(li,lj)=0.0d0
                                B(li,lj)=0.0d0
                            enddo
                        enddo
                        D(2,2)=1.0d0;D(3,3)=1.0d0
                        D(4,4)=1.0d0;D(5,5)=1.0d0
                        f(:,i,j,k)=matmul(D,x(:,i,j,k))+matmul(B,x(:,i,j+1,k))-matmul(B,x(:,i,j,k))
                    else
                        yi=0.0d0
                        call Jor%get_transed_cubes(i,j,k)
                        if(i==0)then
                            lib=0; lie=2
                            ic_index=-2
                        elseif(i==1)then
                            lib=-1; lie=2
                            ic_index=-1
                        elseif(i==(in-2))then
                            lib=-2; lie=1
                            ic_index=1
                        elseif(i==(in-1))then
                            lib=-2; lie=0
                            ic_index=2
                        else
                            lib=-2; lie=2
                            ic_index=0
                        endif
                        if(j==0)then
                            ljb=0; lje=2
                            jc_index=-2
                        elseif(j==1)then
                            ljb=-1; lje=2
                            jc_index=-1
                        elseif(j==(jn-2))then
                            ljb=-2; lje=1
                            jc_index=1
                        elseif(j==(jn-1))then
                            ljb=-2; lje=0
                            jc_index=2
                        else
                            ljb=-2; lje=2
                            jc_index=0
                        endif
                        select case (lns_mode)
                            case(2)
                                lkb=0; lke=0; kc_index=0
                            case(3)
                                lkb=-2; lke=2; kc_index=0
                        end select

                        yi(:,1) = x(:,i,j,k)  ! D

                        do li=lib,lie
                            yi(:,2) = yi(:,2) + x(:,i+li,j,k)*coef_c4(li,ic_index)  ! A
                            yi(:,5) = yi(:,5) + x(:,i+li,j,k)*coef_d4(li,ic_index)  ! Vxx
                        enddo

                        do lj=ljb,lje
                            yi(:,3) = yi(:,3) + x(:,i,j+lj,k)*coef_c4(lj,jc_index)  ! B
                            yi(:,6) = yi(:,6) + x(:,i,j+lj,k)*coef_d4(lj,jc_index)  ! Vxx
                        enddo

                        do lk=lkb,lke
                            yi(:,4) = yi(:,4) + x(:,i,j,k+lk)*coef_c4(lk,kc_index)  ! C
                            yi(:,7) = yi(:,7) + x(:,i,j,k+lk)*coef_d4(lk,kc_index)  ! Vxx
                        enddo

                        do lj=ljb,lje
                            do li=lib,lie
                                yi(:,8)=yi(:,8)+x(:,i+li,j+lj,k)*coef_c4(li,ic_index)*coef_c4(lj,jc_index) ! Vxy
                            enddo
                        enddo

                        do lk=lkb,lke
                            do li=lib,lie
                                yi(:,9)=yi(:,9)+x(:,i+li,j,k+lk)*coef_c4(li,ic_index)*coef_c4(lk,kc_index) ! Vxz
                            enddo
                        enddo

                        do lk=lkb,lke
                            do lj=ljb,lje
                                yi(:,10)=yi(:,10)+x(:,i,j+lj,k+lk)*coef_c4(lj,jc_index)*coef_c4(lk,kc_index) ! Vyz
                            enddo
                        enddo

                        f(:,i,j,k)=matmul(D,yi(:,1))+                                           &
                        &          matmul(A,yi(:,2))+matmul(B,yi(:,3))+matmul(C,yi(:,4))-       &
                        &          matmul(Vxx,yi(:,5))-matmul(Vyy,yi(:,6))-matmul(Vzz,yi(:,7))- &
                        &          matmul(Vxy,yi(:,8))-matmul(Vxz,yi(:,9))-matmul(Vyz,yi(:,10))

                    endif
                enddo
            enddo
        enddo
        end associate
        call DMDAVecRestoreArrayReadF90(meshDA,bell,xr,ierr)
        call DMDAVecRestoreArrayF90(meshDA,F,fr,ierr)
        call DMRestoreLocalVector(meshDA,bell,ierr)

    end subroutine snes_rhs_fx

    subroutine snes_rhs_fx_4ord(snes,x,f,null_int,ierr)

        use mod_parameters,only : meshDA,lns_mode,in,jn,kn,is,ie,js,je,ks,ke
        implicit none
        PetscScalar,pointer :: fr(:,:,:,:),xr(:,:,:,:)
        integer :: ic_index,jc_index,kc_index
        integer :: lib,lie,ljb,lje,lkb,lke
        PetscScalar,dimension(5,16) :: yi
        PetscErrorCode :: ierr
        integer :: null_int(*)
        integer :: li,lj,lk
        integer :: i,j,k
        SNES :: snes
        Vec :: bell
        Vec :: x,f
        call DMGetLocalVector(meshDA,bell,ierr)
        call VecZeroEntries(bell,ierr)
        call DMGlobalToLocalBegin(meshDA,X,INSERT_VALUES,bell,ierr)
        call DMGlobalToLocalEnd(meshDA,X,INSERT_VALUES,bell,ierr)
        call DMDAVecGetArrayReadF90(meshDA,bell,xr,ierr)
        call DMDAVecGetArrayF90(meshDA,F,fr,ierr)
        associate ( &
            &   coef_c4  => FDM_1nd_4ORD_CENTER,   &
            &   coef_d4  => FDM_2nd_4ORD_CENTER,   &
            &   coef_c4f => FDM_1nd_4ORD_Forward,  &
            &   coef_c4b => FDM_1nd_4ORD_Backward, &
            &   f   => fr,        x => xr,         &
            &   G   => Jor%G,     D => Jor%D,      &
            &   A   => Jor%A,     B => Jor%B,     C => jor%C,   &
            &   A_p => Jor%A_p, A_m => Jor%A_m, A_v => Jor%A_v, &
            &   B_p => Jor%B_p, B_m => Jor%B_m, B_v => Jor%B_v, &
            &   C_p => Jor%C_p, C_m => Jor%C_m, C_v => Jor%C_v, &
            &   Vxx => Jor%Vxx, Vyy => Jor%Vyy, Vzz => Jor%Vzz, &
            &   Vxy => Jor%Vxy, Vxz => Jor%Vxz, Vyz => Jor%Vyz)
        do k=ks,ke
            do j=js,je
                do i=is,ie
                    if(i==0)then
                        f(:,i,j,k)=x(:,i,j,k)
                    elseif(j==(jn-1))then
                        f(:,i,j,k)=x(:,i,j,k)
                    elseif(i==(in-1) .and. j/=0 .and. j/=(jn-1))then
                        f(:,i,j,k)=(x(:,i-2,j,k)-4*x(:,i-1,j,k)+3*x(:,i,j,k))/2.0d0
                    elseif(j==0 .and. i/=0)then
                        call Jor%get_transed_cubes(i,j,k)
                        do lj=1,5
                            do li=2,5
                                D(li,lj)=0.0d0
                                B(li,lj)=0.0d0
                            enddo
                        enddo
                        D(2,2)=1.0d0;D(3,3)=1.0d0
                        D(4,4)=1.0d0;D(5,5)=1.0d0
                        f(:,i,j,k)=matmul(D,x(:,i,j,k))+matmul(B,x(:,i,j+1,k))-matmul(B,x(:,i,j,k))
                    else
                        yi=0.0d0
                        call Jor%get_transed_cubes(i,j,k)
                        if(i==0)then
                            lib=0; lie=2
                            ic_index=-2
                        elseif(i==1)then
                            lib=-1; lie=2
                            ic_index=-1
                        elseif(i==(in-2))then
                            lib=-2; lie=1
                            ic_index=1
                        elseif(i==(in-1))then
                            lib=-2; lie=0
                            ic_index=2
                        else
                            lib=-2; lie=2
                            ic_index=0
                        endif
                        if(j==0)then
                            ljb=0; lje=2
                            jc_index=-2
                        elseif(j==1)then
                            ljb=-1; lje=2
                            jc_index=-1
                        elseif(j==(jn-2))then
                            ljb=-2; lje=1
                            jc_index=1
                        elseif(j==(jn-1))then
                            ljb=-2; lje=0
                            jc_index=2
                        else
                            ljb=-2; lje=2
                            jc_index=0
                        endif
                        select case (lns_mode)
                            case(2)
                                lkb=0; lke=0; kc_index=0
                            case(3)
                                lkb=-2; lke=2; kc_index=0
                        end select

                        yi(:,1) = x(:,i,j,k) ! D

                        do li=lib,lie
                            yi(:,2)  = yi(:,2)  + x(:,i+li,j,k)*coef_c4b(li,ic_index) ! A_p
                            yi(:,3)  = yi(:,3)  + x(:,i+li,j,k)*coef_c4f(li,ic_index) ! A_m
                            yi(:,4)  = yi(:,4)  + x(:,i+li,j,k)*coef_c4(li,ic_index)  ! A_v
                            yi(:,11) = yi(:,11) + x(:,i+li,j,k)*coef_d4(li,ic_index)  ! Vxx
                        enddo

                        do lj=ljb,lje
                            yi(:,5)  = yi(:,5)  + x(:,i,j+lj,k)*coef_c4b(lj,jc_index) ! B_p
                            yi(:,6)  = yi(:,6)  + x(:,i,j+lj,k)*coef_c4f(lj,jc_index) ! B_m
                            yi(:,7)  = yi(:,7)  + x(:,i,j+lj,k)*coef_c4(lj,jc_index)  ! B_v
                            yi(:,12) = yi(:,12) + x(:,i,j+lj,k)*coef_d4(lj,jc_index)  ! Vxx
                        enddo

                        do lk=lkb,lke
                            yi(:,8)  = yi(:,8)  + x(:,i,j,k+lk)*coef_c4b(lk,kc_index) ! C_p
                            yi(:,9)  = yi(:,9)  + x(:,i,j,k+lk)*coef_c4f(lk,kc_index) ! C_m
                            yi(:,10) = yi(:,10) + x(:,i,j,k+lk)*coef_c4(lk,kc_index)  ! C_v
                            yi(:,13) = yi(:,13) + x(:,i,j,k+lk)*coef_d4(lk,kc_index)  ! Vxx
                        enddo

                        do lj=ljb,lje
                            do li=lib,lie
                                yi(:,14)=yi(:,14)+x(:,i+li,j+lj,k)*coef_c4(li,ic_index)*coef_c4(lj,jc_index) ! Vxy
                            enddo
                        enddo

                        do lk=lkb,lke
                            do li=lib,lie
                                yi(:,15)=yi(:,15)+x(:,i+li,j,k+lk)*coef_c4(li,ic_index)*coef_c4(lk,kc_index) ! Vxz
                            enddo
                        enddo

                        do lk=lkb,lke
                            do lj=ljb,lje
                                yi(:,16)=yi(:,16)+x(:,i,j+lj,k+lk)*coef_c4(lj,jc_index)*coef_c4(lk,kc_index) ! Vyz
                            enddo
                        enddo

                        f(:,i,j,k)=matmul(D,yi(:,1))+                                              &
                        &          matmul(A_p,yi(:,2))+matmul(A_m,yi(:,3))+matmul(A_v,yi(:,4))+    &
                        &          matmul(B_p,yi(:,5))+matmul(B_m,yi(:,6))+matmul(B_v,yi(:,7))+    &
                        &          matmul(C_p,yi(:,8))+matmul(C_m,yi(:,9))+matmul(C_v,yi(:,10))-   &
                        &          matmul(Vxx,yi(:,11))-matmul(Vyy,yi(:,12))-matmul(Vzz,yi(:,13))- &
                        &          matmul(Vxy,yi(:,14))-matmul(Vxz,yi(:,15))-matmul(Vyz,yi(:,16))

                    endif
                enddo
            enddo
        enddo
        end associate
        call DMDAVecRestoreArrayReadF90(meshDA,bell,xr,ierr)
        call DMDAVecRestoreArrayF90(meshDA,F,fr,ierr)
        call DMRestoreLocalVector(meshDA,bell,ierr)

    end subroutine snes_rhs_fx_4ord

    subroutine snes_converged_test(snes,it,xnorm,snorm,fnorm,reason,dummy,ierr)
        implicit none
        PetscReal :: xnorm,snorm,fnorm,nrm
        SNESConvergedReason :: reason
        PetscErrorCode :: ierr
        PetscInt :: it,dummy
        SNES :: snes
        Vec :: f
        character(len=20) :: str_norm
        character(len=6) :: str_it

        call SNESGetFunction(snes,f,PETSC_NULL_FUNCTION,dummy,ierr)
        call VecNorm(f,NORM_INFINITY,nrm,ierr)
        write(str_it,"(I5)") it
        write(str_norm,"(ES20.12)") nrm
        call PetscPrintf(PETSC_COMM_WORLD," "//str_it//" < residual i-Norm > "//str_norm//"\n",ierr)
        if (nrm .le. 1.e-5) reason = SNES_CONVERGED_FNORM_ABS
    end subroutine snes_converged_test

    ! -----------------------------------------------------------------------------------------------------
    !   迭代格式三：借用KSP模块

    subroutine push_bc(comm,x)
        use mod_parameters,only:meshDA,disturb,bf,in,jn,kn,is,ie,js,je,ks,ke
        implicit none
        PetscReal,parameter :: d4d3=4.0d0/3.0d0,d1d3=1.0d0/3.0d0
        PetscScalar,pointer :: xr(:,:,:,:)
        integer,intent(in) :: comm
        Vec,intent(inout) :: x
        PetscErrorCode :: ierr
        integer :: i,j,k

        call DMDAVecGetArrayF90(meshDA,x,xr,ierr)

        associate( x => xr )

        if(is==0)then
            do k=ks,ke
                do j=js,je
                    do i=is,is
                        x(:,i,j,k)=disturb(:,j,k)
                    enddo
                enddo
            enddo
        endif

        if(ie==(in-1))then
            do k=ks,ke
                do j=js,je
                    do i=ie,ie
                        x(:,i,j,k)=d4d3*x(:,i-1,j,k)-d1d3*x(:,i-2,j,k)
                    enddo
                enddo
            enddo
        endif

        if(js==0)then
            do k=ks,ke
                do j=js,js
                    do i=is,ie
                        x(:,i,j,k)=0.0d0
                        x(0,i,j,k)=(bf(i,j+1,k)%BF%rho*x(4,i,j+1,k)+bf(i,j+1,k)%BF%T*x(0,i,j+1,k))/bf(i,j,k)%BF%T
                    enddo
                enddo
            enddo
        endif

        if(je==(jn-1))then
            do k=ks,ke
                do j=je,je
                    do i=is,ie
                        x(:,i,j,k)=0.0d0
                    enddo
                enddo
            enddo
        endif

        end associate

        call DMDAVecRestoreArrayF90(meshDA,x,xr,ierr)

        call MPI_Barrier(comm,ierr)

    end subroutine push_bc

    subroutine rhs_fx_Ax(x,f)

        use mod_parameters,only : meshDA,lns_mode,in,jn,kn,is,ie,js,je,ks,ke
        implicit none
        PetscScalar,pointer :: fr(:,:,:,:),xr(:,:,:,:)
        integer :: ic_index,jc_index,kc_index
        integer :: lib,lie,ljb,lje,lkb,lke
        PetscScalar,dimension(5,16) :: yi
        PetscErrorCode :: ierr
        Vec,intent(inout) :: f
        Vec,intent(in) :: x
        integer :: li,lj,lk
        integer :: i,j,k
        Vec :: bell

        call DMGetLocalVector(meshDA,bell,ierr)
        call VecZeroEntries(bell,ierr)
        call DMGlobalToLocalBegin(meshDA,x,INSERT_VALUES,bell,ierr)
        call DMGlobalToLocalEnd(meshDA,x,INSERT_VALUES,bell,ierr)
        call DMDAVecGetArrayReadF90(meshDA,bell,xr,ierr)
        call DMDAVecGetArrayF90(meshDA,f,fr,ierr)
        associate ( &
        &   coef_c4  => FDM_1nd_4ORD_CENTER,   &
        &   coef_d4  => FDM_2nd_4ORD_CENTER,   &
        &   coef_c4f => FDM_1nd_4ORD_Forward,  &
        &   coef_c4b => FDM_1nd_4ORD_Backward, &
        &   f   => fr,        x => xr,         &
        &   G   => Jor%G,     D => Jor%D,      &
        &   A   => Jor%A,     B => Jor%B,     C => jor%C,   &
        &   A_p => Jor%A_p, A_m => Jor%A_m, A_v => Jor%A_v, &
        &   B_p => Jor%B_p, B_m => Jor%B_m, B_v => Jor%B_v, &
        &   C_p => Jor%C_p, C_m => Jor%C_m, C_v => Jor%C_v, &
        &   Vxx => Jor%Vxx, Vyy => Jor%Vyy, Vzz => Jor%Vzz, &
        &   Vxy => Jor%Vxy, Vxz => Jor%Vxz, Vyz => Jor%Vyz)
        do k=ks,ke
            do j=js,je
                do i=is,ie
                    if(i==0 .or. i==(in-1) .or. j==0 .or. j==(jn-1))then
                        f(:,i,j,k)=0.0d0
                    else
                        yi=0.0d0
                        call Jor%get_transed_cubes(i,j,k)
                        if(i==0)then
                            lib=0; lie=2
                            ic_index=-2
                        elseif(i==1)then
                            lib=-1; lie=2
                            ic_index=-1
                        elseif(i==(in-2))then
                            lib=-2; lie=1
                            ic_index=1
                        elseif(i==(in-1))then
                            lib=-2; lie=0
                            ic_index=2
                        else
                            lib=-2; lie=2
                            ic_index=0
                        endif
                        if(j==0)then
                            ljb=0; lje=2
                            jc_index=-2
                        elseif(j==1)then
                            ljb=-1; lje=2
                            jc_index=-1
                        elseif(j==(jn-2))then
                            ljb=-2; lje=1
                            jc_index=1
                        elseif(j==(jn-1))then
                            ljb=-2; lje=0
                            jc_index=2
                        else
                            ljb=-2; lje=2
                            jc_index=0
                        endif
                        select case (lns_mode)
                            case(2)
                                lkb=0; lke=0; kc_index=0
                            case(3)
                                lkb=-2; lke=2; kc_index=0
                        end select

                        yi(:,1) = x(:,i,j,k) ! D

                        do li=lib,lie
                            yi(:,2)  = yi(:,2)  + x(:,i+li,j,k)*coef_c4b(li,ic_index) ! A_p
                            yi(:,3)  = yi(:,3)  + x(:,i+li,j,k)*coef_c4f(li,ic_index) ! A_m
                            yi(:,4)  = yi(:,4)  + x(:,i+li,j,k)*coef_c4(li,ic_index)  ! A_v
                            yi(:,11) = yi(:,11) + x(:,i+li,j,k)*coef_d4(li,ic_index)  ! Vxx
                        enddo

                        do lj=ljb,lje
                            yi(:,5)  = yi(:,5)  + x(:,i,j+lj,k)*coef_c4b(lj,jc_index) ! B_p
                            yi(:,6)  = yi(:,6)  + x(:,i,j+lj,k)*coef_c4f(lj,jc_index) ! B_m
                            yi(:,7)  = yi(:,7)  + x(:,i,j+lj,k)*coef_c4(lj,jc_index)  ! B_v
                            yi(:,12) = yi(:,12) + x(:,i,j+lj,k)*coef_d4(lj,jc_index)  ! Vxx
                        enddo

                        do lk=lkb,lke
                            yi(:,8)  = yi(:,8)  + x(:,i,j,k+lk)*coef_c4b(lk,kc_index) ! C_p
                            yi(:,9)  = yi(:,9)  + x(:,i,j,k+lk)*coef_c4f(lk,kc_index) ! C_m
                            yi(:,10) = yi(:,10) + x(:,i,j,k+lk)*coef_c4(lk,kc_index)  ! C_v
                            yi(:,13) = yi(:,13) + x(:,i,j,k+lk)*coef_d4(lk,kc_index)  ! Vxx
                        enddo

                        do lj=ljb,lje
                            do li=lib,lie
                                yi(:,14)=yi(:,14)+x(:,i+li,j+lj,k)*coef_c4(li,ic_index)*coef_c4(lj,jc_index) ! Vxy
                            enddo
                        enddo

                        do lk=lkb,lke
                            do li=lib,lie
                                yi(:,15)=yi(:,15)+x(:,i+li,j,k+lk)*coef_c4(li,ic_index)*coef_c4(lk,kc_index) ! Vxz
                            enddo
                        enddo

                        do lk=lkb,lke
                            do lj=ljb,lje
                                yi(:,16)=yi(:,16)+x(:,i,j+lj,k+lk)*coef_c4(lj,jc_index)*coef_c4(lk,kc_index) ! Vyz
                            enddo
                        enddo

                        f(:,i,j,k)=-1.0d0*(matmul(D,yi(:,1))+                                      &
                        &          matmul(A_p,yi(:,2))+matmul(A_m,yi(:,3))+matmul(A_v,yi(:,4))+    &
                        &          matmul(B_p,yi(:,5))+matmul(B_m,yi(:,6))+matmul(B_v,yi(:,7))+    &
                        &          matmul(C_p,yi(:,8))+matmul(C_m,yi(:,9))+matmul(C_v,yi(:,10))-   &
                        &          matmul(Vxx,yi(:,11))-matmul(Vyy,yi(:,12))-matmul(Vzz,yi(:,13))- &
                        &          matmul(Vxy,yi(:,14))-matmul(Vxz,yi(:,15))-matmul(Vyz,yi(:,16)))

                    endif
                enddo
            enddo
        enddo
        end associate

        call DMDAVecRestoreArrayReadF90(meshDA,bell,xr,ierr)
        call DMDAVecRestoreArrayF90(meshDA,f,fr,ierr)
        call DMRestoreLocalVector(meshDA,bell,ierr)

    end subroutine rhs_fx_Ax

    subroutine ksps_rhs_fx_b_Gtx(x,f)
       use mod_parameters,only : meshDA,is,ie,js,je,ks,ke,disturb,dt
       implicit none
       PetscScalar,pointer :: fr(:,:,:,:),xr(:,:,:,:)
       PetscErrorCode :: ierr
       Vec,intent(inout) :: f
       Vec,intent(in) :: x
       integer :: i,j,k
       PetscReal :: idt

       idt = 1.0d0/dt

       call DMDAVecGetArrayReadF90(meshDA,x,xr,ierr)
       call DMDAVecGetArrayF90(meshDA,f,fr,ierr)

       associate(    &
       &    x => xr, &
       &    f => fr)
       do k=ks,ke
           do j=js,je
               do i=is,ie
                   call Jor%get_adorned_cubes(i,j,k)
                   if(i==0)then
                       f(:,i,j,k)=disturb(:,j,k)+matmul(idt*Jor%G,x(:,i,j,k))
                   else
                       f(:,i,j,k)=matmul(idt*Jor%G,x(:,i,j,k))
                   endif
               enddo
           enddo
       enddo
       end associate

       call DMDAVecRestoreArrayReadF90(meshDA,x,xr,ierr)
       call DMDAVecRestoreArrayF90(meshDA,f,fr,ierr)
    end subroutine ksps_rhs_fx_b_Gtx

    ! --------------------------------------------------------------------------
    ! KSPs

    subroutine form_sub_mat_2_ord(mat)  !!设置子块当地矩阵

        use mod_parameters,only : igs,ige,jgs,jge,kgs,kge,rank
        implicit none
        Mat,intent(inout) :: mat
        PetscErrorCode :: ierr
        integer :: i,j,k

        do k=kgs,kge
            do j=jgs,jge
                do i=igs,ige
                    if(i==igs .or. i==ige .or. j==jgs .or. j==jge)then !!暂时不考虑展向
                        call sub_mat_set_boundary_conditions(mat,i,j,k)  !!子矩阵边界点的\delta\hat\phi=0
                    else
                        call sub_mat_insert_values_2_ord(mat,i,j,k)  !!设置子矩阵内点的值
                    endif
                enddo
            enddo
        enddo

        call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)

    end subroutine form_sub_mat_2_ord

    subroutine sub_mat_set_boundary_conditions(mat,i,j,k)
        use mod_parameters,only:igs,jgs,kgs
        implicit none
        MatStencil :: idxm(4,1),idxn(4,1)
        integer,intent(in) :: i,j,k
        Mat,intent(inout) :: mat
        PetscScalar :: box(5,5)
        PetscErrorCode :: ierr

        box=0.0d0
        box(1,1)=1.0d0;box(2,2)=1.0d0;box(3,3)=1.0d0;box(4,4)=1.0d0;box(5,5)=1.0d0
        idxm(MatStencil_i, 1)=i-igs; idxm(MatStencil_j, 1)=j-jgs; idxm(MatStencil_k, 1)=k-kgs
        idxn(MatStencil_i, 1)=i-igs; idxn(MatStencil_j, 1)=j-jgs; idxn(MatStencil_k, 1)=k-kgs
        call MatSetValuesBlockedStencil(mat, 1, idxm, 1, idxn, box, INSERT_VALUES, ierr)

    end subroutine sub_mat_set_boundary_conditions

    subroutine sub_mat_insert_values_2_ord(mat,i,j,k)

        use mod_parameters,only : lns_mode,igs,jgs,kgs,bf
        implicit none
        integer :: ic_index, jc_index, kc_index
        integer :: lib, lie, ljb, lje, lkb, lke
        PetscScalar :: box(5,5),trans(5,5)
        MatStencil :: idxm(4,1),idxn(4,1)
        integer,intent(in) :: i,j,k
        Mat,intent(inout) :: mat
        PetscErrorCode :: ierr
        integer :: li,lj,lk
        integer :: is,js,ks

        call Jor%get_adorned_cubes(i,j,k)
        lib=-1; lie=1
        ic_index=0
        ljb=-1; lje=1
        jc_index=0
        select case (lns_mode)
            case(2)
                lkb=0; lke=0; kc_index=0
            case(3)
                lkb=-1; lke=1; kc_index=0
        end select
        idxm=0;
        idxm(MatStencil_i, 1)=i-igs; idxm(MatStencil_j, 1)=j-jgs; idxm(MatStencil_k, 1)=k-kgs
        associate( &
            &   coef_c1f=>FDM_1nd_1ORD_Forward,  &
            &   coef_c1b=>FDM_1nd_1ORD_Backward, &
            &   coef_d2 =>FDM_2nd_2ORD_CENTER,   &
            &   coef_c2 =>FDM_1nd_2ORD_CENTER,   &
            &   G   => Jor%G,   D   => Jor%D,    &
            &   A_p => Jor%A_p, A_m => Jor%A_m, A_v => Jor%A_v, &
            &   B_p => Jor%B_p, B_m => Jor%B_m, B_v => Jor%B_v, &
            &   C_p => Jor%C_p, C_m => Jor%C_m, C_v => Jor%C_v, &
            &   Vxx => Jor%Vxx, Vyy => Jor%Vyy, Vzz => Jor%Vzz, &
            &   Vxy => Jor%Vxy, Vxz => Jor%Vxz, Vyz => Jor%Vyz)
            do lk = lkb, lke
                do lj = ljb, lje
                    do li = lib, lie
                        idxn=0;box=0.0d0;trans=0.0d0
                        idxn(MatStencil_i, 1)=i+li-igs
                        idxn(MatStencil_j, 1)=j+lj-jgs
                        idxn(MatStencil_k, 1)=k+lk-kgs
                        box=delta_i(li)*delta_j(lj)*delta_k(lk)*D+                     &
                        &   delta_j(lj)*delta_k(lk)*A_v*coef_c2(li,ic_index)+          &
                        &   delta_j(lj)*delta_k(lk)*A_m*coef_c1f(li,ic_index)+         &
                        &   delta_j(lj)*delta_k(lk)*A_p*coef_c1b(li,ic_index)+         &
                        &   delta_i(li)*delta_k(lk)*B_v*coef_c2(lj,jc_index)+          &
                        &   delta_i(li)*delta_k(lk)*B_m*coef_c1f(lj,jc_index)+         &
                        &   delta_i(li)*delta_k(lk)*B_p*coef_c1b(lj,jc_index)+         &
                        &   delta_i(li)*delta_j(lj)*C_v*coef_c2(lk,kc_index)+          &
                        &   delta_i(li)*delta_j(lj)*C_m*coef_c1f(lk,kc_index)+         &
                        &   delta_i(li)*delta_j(lj)*C_p*coef_c1b(lk,kc_index)-         &
                        &   delta_j(lj)*delta_k(lk)*Vxx*coef_d2(li,ic_index)-          &
                        &   delta_i(li)*delta_k(lk)*Vyy*coef_d2(lj,jc_index)-          &
                        &   delta_i(li)*delta_j(lj)*Vzz*coef_d2(lk,kc_index)-          &
                        &   delta_k(lk)*Vxy*coef_c2(li,ic_index)*coef_c2(lj,jc_index)- &
                        &   delta_j(lj)*Vxz*coef_c2(li,ic_index)*coef_c2(lk,kc_index)- &
                        &   delta_i(li)*Vyz*coef_c2(lj,jc_index)*coef_c2(lk,kc_index)
                        trans=transpose(box)
                        call MatSetValuesBlockedStencil(mat, 1, idxm, 1, idxn, trans, INSERT_VALUES, ierr)
                    enddo
                enddo
            enddo
        end associate

    end subroutine sub_mat_insert_values_2_ord

    subroutine set_subDA(da)
        use mod_parameters,only:igl,jgl,kgl
        implicit none
        DM,intent(inout) :: da
        PetscErrorCode :: ierr

        call DMDACreate3d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC,&
        &                 DMDA_STENCIL_BOX, igl, jgl, kgl, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, &
        &                 5, 1, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, da, ierr)
        call DMSetMatType(da,MATBAIJ,ierr)
        call DMSetFromOptions(da,ierr)
        call DMSetUp(da, ierr)
    end subroutine set_subDA

    subroutine init_sub_vecs()
        use mod_parameters,only:meshDA,subDA,localx,subx
        implicit none
        PetscErrorCode :: ierr

        call DMGetLocalVector(meshDA,localx,ierr)
        call VecZeroEntries(localx,ierr)
        call DMGetGlobalVector(subDA,subx,ierr)
        call VecZeroEntries(subx,ierr)
    end subroutine init_sub_vecs

    subroutine get_subf(f,subf)
        use mod_parameters,only:lns_mode,igs,ige,jgs,jge,kgs,kge,localx,meshDA,subDA
        implicit none
        PetscScalar,pointer :: fr(:,:,:,:),sfr(:,:,:,:)
        Vec,intent(inout) :: subf
        PetscErrorCode :: ierr
        Vec,intent(in) :: f
        integer :: i,j,k
        Vec :: localf

        call VecDuplicate(localx,localf,ierr)
        call VecZeroEntries(localf,ierr)

        call DMGlobalToLocalBegin(meshDA, f, INSERT_VALUES, localf, ierr)
        call DMGlobalToLocalEnd(meshDA, f, INSERT_VALUES, localf, ierr)
        call DMDAVecGetArrayF90(meshDA, localf, fr, ierr)
        call DMDAVecGetArrayF90(subDA, subf, sfr, ierr)

        select case(lns_mode)
        case(2)
            fr(:, igs:igs, :, :)=0.0d0
            fr(:, ige:ige, :, :)=0.0d0
            fr(:, :, jgs:jgs, :)=0.0d0
            fr(:, :, jge:jge, :)=0.0d0
        case(3)
            fr(:, igs:igs, :, :)=0.0d0
            fr(:, ige:ige, :, :)=0.0d0
            fr(:, :, jgs:jgs, :)=0.0d0
            fr(:, :, jge:jge, :)=0.0d0
            fr(:, :, :, kgs:kgs)=0.0d0
            fr(:, :, :, kge:kge)=0.0d0
        end select

        do k=kgs, kge
            do j=jgs, jge
                do i=igs, ige
                    sfr(:, i-igs, j-jgs, k-kgs)=fr(:, i, j, k)
                enddo
            enddo
        enddo

        call DMDAVecRestoreArrayF90(meshDA, localf, fr, ierr)
        call DMDAVecRestoreArrayF90(subDA, subf, sfr, ierr)
        call VecDestroy(localf,ierr)
    end subroutine get_subf

    subroutine merge_res(subres,res)
        use mod_parameters,only:subDA,meshDA,is,ie,js,je,ks,ke,igs,jgs,kgs
        implicit none
        PetscScalar,pointer :: sr(:,:,:,:),r(:,:,:,:)
        Vec,intent(inout) :: res
        Vec,intent(in) :: subres
        PetscErrorCode :: ierr
        integer :: i,j,k

        call DMDAVecGetArrayReadF90(subDA,subres,sr,ierr)
        call DMDAVecGetArrayF90(meshDA,res,r,ierr)

        do k=ks, ke
            do j=js, je
                do i=is, ie
                    r(:,i,j,k)=sr(:,i-igs,j-jgs,k-kgs)
                enddo
            enddo
        enddo

        call DMDAVecRestoreArrayF90(meshDA,res,r,ierr)
        call DMDAVecRestoreArrayReadF90(subDA,subres,sr,ierr)

    end subroutine merge_res

    ! --------------------------------------------------------------------------
    ! 清理函数

    subroutine cleanup()

        use mod_parameters
        deallocate(bf)
        deallocate(disturb)
        deallocate(xi_x,xi_y,xi_z)
        deallocate(eta_x,eta_y,eta_z)
        deallocate(phi_x,phi_y,phi_z)
        deallocate(xi_xx,xi_yy,xi_zz)
        deallocate(eta_xx,eta_yy,eta_zz)
        deallocate(phi_xx,phi_yy,phi_zz)
        deallocate(xi_xy,xi_xz,xi_yz)
        deallocate(eta_xy,eta_yz,eta_xz)
        deallocate(phi_xy,phi_yz,phi_xz)

    end subroutine cleanup

end module mod_forming
