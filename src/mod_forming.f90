!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_forming
! ------------------------------------------------------------------
!
!  这个模块包含调用xxxSolve前的所有准备步骤，包括矩阵、右端量的生成函数，清理函数等。
!
!   在线性求解系统 KSP 框架下，生成左端矩阵；
!   在非线性求解系统 SNES 框架下，生成雅各比矩阵和右端项函数。
!
!       1.call set_matfree_mat(comm,ierr) 免矩阵形式的矩阵生成函数
!
!           call mat_mult_4_precision(A, X, F, ierr) 免矩阵需要的矩阵向量乘法函数
!
!       2.call get_mat_from_dm(comm,ierr) 显式矩阵形式的矩阵生成函数
!
!           call form_mat_4_precision(ierr) 填充数据函数
!
!               (1).call mat_set_boundary_conditions(i,j,k) 边界部分的填充数据函数
!
!               (2).call mat_insert_values_4_precision(i,j,k) 内部部分的填充数据函数
!
!       3.call set_rhs(comm,ierr) 设置右边量，即边界。
!
!       4.call cleanup() 清理
!
!       5.call shark_coming(comm,ierr) 一阶精度分裂的矩阵
!
!           call form_mat_2_precision(ierr) 生成一阶精度分裂的矩阵
!
!       6.call snes_fx1o(snes,x,f,null_int,ierr) 一阶精度右端项函数
!
!       7.call snes_fx4o(snes,x,f,null_int,ierr) 四阶精度右端项函数
!
! ------------------------------------------------------------------
    use penf, only: R_P
    use mod_cubes
    use petsc
    implicit none
    public :: set_matfree_mat
    public :: mat_mult_4_precision
    public :: initialize_mat_from_da
    public :: form_mat_4_precision
    public :: form_mat_2_precision
    public :: form_mat_4_precision_Gtx
    public :: form_mat_2_precision_Gtx
    public :: snes_fx,snes_fx4o
    public :: ksps_rhs_fx_b_Ax
    public :: ksps_rhs_fx_b_Gtx
    public :: set_rhs
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

    subroutine set_matfree_mat(comm,x,mat)

        implicit none
        PetscInt,intent(in) :: comm
        Mat,intent(inout) :: mat
        PetscErrorCode :: ierr
        Vec,intent(in) :: x
        PetscInt :: ls

        call VecGetLocalSize(x,ls,ierr)
        call MatCreateShell(comm,ls,ls,PETSC_DETERMINE,PETSC_DETERMINE,PETSC_NULL_INTEGER,mat,ierr)
        call MatShellSetOperation(mat,MATOP_MULT,mat_mult_4_precision,ierr)
        call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)

    end subroutine set_matfree_mat

    subroutine mat_mult_4_precision(A,X,F,ierr)

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
                        call Jor%get_adorned_cubes(i,j,k)
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

    end subroutine mat_mult_4_precision

    subroutine initialize_mat_from_da(comm,da,mat)

        implicit none
        PetscInt,intent(in) :: comm
        PetscErrorCode :: ierr
        Mat :: mat
        DM :: da

        call DMCreateMatrix(da,mat,ierr)
        call MatZeroEntries(mat,ierr)

    end subroutine initialize_mat_from_da

    subroutine form_mat_4_precision(mat)

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
                        call mat_insert_values_4_precision(mat,i,j,k)
                    endif
                enddo
            enddo
        enddo

        call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)

    end subroutine form_mat_4_precision

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
            call Jor%get_adorned_cubes(i,j,k)
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

    subroutine mat_insert_values_4_precision(mat,i,j,k)

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

    end subroutine mat_insert_values_4_precision

    subroutine set_rhs(comm,da,x,r)

        use mod_parameters,only : disturb,is,ie,js,je,ks,ke
        implicit none
        PetscScalar,pointer :: r_array(:,:,:,:)
        PetscInt,intent(in) :: comm
        Vec,intent(inout) :: r
        PetscErrorCode :: ierr
        Vec,intent(in) :: x
        DM,intent(in) :: da
        integer :: j,k

        call VecDuplicate(x,r,ierr)
        call VecZeroEntries(r,ierr)
        if (is==0) then
            call DMDAVecGetArrayF90(da,r,r_array,ierr)
                do k=ks,ke
                    do j=js,je
                        r_array(:,0,j,k)=disturb(:,j,k)
                    enddo
                enddo
            call DMDAVecRestoreArrayF90(da,r,r_array,ierr)
        endif
        call MPI_Barrier(comm,ierr)

    end subroutine set_rhs

    ! -----------------------------------------------------------------------------------------------------
    !   迭代格式二：借用SNES模块 Jac x = - F(x)

    subroutine form_mat_2_precision(mat)

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
                        call mat_insert_values_2_precision(mat,i,j,k)
                    endif
                enddo
            enddo
        enddo

        call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)

    end subroutine form_mat_2_precision

    subroutine mat_insert_values_2_precision(mat,i,j,k)

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

    end subroutine mat_insert_values_2_precision

    subroutine snes_fx(snes,x,f,null_int,ierr)

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
                        call Jor%get_adorned_cubes(i,j,k)
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

    end subroutine snes_fx

    subroutine snes_fx4o(snes,x,f,null_int,ierr)

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
                        call Jor%get_adorned_cubes(i,j,k)
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

    end subroutine snes_fx4o

    ! -----------------------------------------------------------------------------------------------------
    !   迭代格式三：借用KSP模块

    subroutine form_mat_4_precision_Gtx(mat)

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
                        call mat_insert_values_4_precision_handin(mat,i,j,k)
                    endif
                enddo
            enddo
        enddo

        call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)

    end subroutine form_mat_4_precision_Gtx

    subroutine mat_insert_values_4_precision_handin(mat,i,j,k)

        use mod_parameters,only : lns_mode,in,jn,kn,ck
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
                    box=delta_i(li)*delta_j(lj)*delta_k(lk)/ck*G+                  &
                    &   delta_i(li)*delta_j(lj)*delta_k(lk)*D+                     &
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

    end subroutine mat_insert_values_4_precision_handin

    subroutine ksps_rhs_fx_b_Ax(x,f)

        use mod_parameters,only : meshDA,lns_mode,in,jn,kn,is,ie,js,je,ks,ke,disturb
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
                    if(i==0)then
                        f(:,i,j,k)=disturb(:,j,k)-x(:,i,j,k)
                    elseif(j==(jn-1))then
                        f(:,i,j,k)=-1.0d0*x(:,i,j,k)
                    elseif(i==(in-1) .and. j/=0 .and. j/=(jn-1))then
                        f(:,i,j,k)=-1.0d0*(x(:,i-2,j,k)-4*x(:,i-1,j,k)+3*x(:,i,j,k))/2.0d0
                    elseif(j==0 .and. i/=0)then
                        call Jor%get_adorned_cubes(i,j,k)
                        do lj=1,5
                            do li=2,5
                                D(li,lj)=0.0d0
                                B(li,lj)=0.0d0
                            enddo
                        enddo
                        D(2,2)=1.0d0;D(3,3)=1.0d0
                        D(4,4)=1.0d0;D(5,5)=1.0d0
                        f(:,i,j,k)=-1.0d0*(matmul(D,x(:,i,j,k))+matmul(B,x(:,i,j+1,k))-matmul(B,x(:,i,j,k)))
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

    end subroutine ksps_rhs_fx_b_Ax

    subroutine ksps_rhs_fx_b_Gtx(x,f)
       use mod_parameters,only : meshDA,is,ie,js,je,ks,ke,disturb,ck
       implicit none
       PetscScalar,pointer :: fr(:,:,:,:),xr(:,:,:,:)
       PetscErrorCode :: ierr
       Vec,intent(inout) :: f
       Vec,intent(in) :: x
       integer :: i,j,k
       PetscReal :: ick

       ick = 1.0d0/ck

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
                       f(:,i,j,k)=disturb(:,j,k)+matmul(ick*Jor%G,x(:,i,j,k))
                   else
                       f(:,i,j,k)=matmul(ick*Jor%G,x(:,i,j,k))
                   endif
               enddo
           enddo
       enddo
       end associate

       call DMDAVecRestoreArrayReadF90(meshDA,x,xr,ierr)
       call DMDAVecRestoreArrayF90(meshDA,f,fr,ierr)
   end subroutine ksps_rhs_fx_b_Gtx

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

    subroutine form_mat_2_precision_Gtx(mat)

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
                        call mat_insert_values_2_precision_handin(mat,i,j,k)
                    endif
                enddo
            enddo
        enddo

        call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)

    end subroutine form_mat_2_precision_Gtx

    subroutine mat_insert_values_2_precision_handin(mat,i,j,k)

        use mod_parameters,only : lns_mode,in,jn,kn,ck
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
                        box=delta_i(li)*delta_j(lj)*delta_k(lk)/ck*G+                  &
                        &   delta_i(li)*delta_j(lj)*delta_k(lk)*D+                     &
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

    end subroutine mat_insert_values_2_precision_handin

end module mod_forming
