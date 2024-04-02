!------------------------------------------------------------------------------
!
! Copyright (C) 2019-2024 Bzhouha
! 
! This file is part of START_LNS.
!
! START_LNS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! START_LNS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>. 
!
!------------------------------------------------------------------------------
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
!           call form_global_mat_4ord(ierr) 填充数据函数
!
!               (1).call mat_set_boundary_conditions(i,j,k) 边界部分的填充数据函数
!
!               (2).call mat_insert_values_4_ord(i,j,k) 内部部分的填充数据函数
!
!       3.call linear_rhs(comm,ierr) 设置右边量，即边界。
!
!       4.call cleanup() 清理
!
!       5.call shark_coming(comm,ierr) 一阶精度分裂的矩阵
!
!           call form_global_mat_2ord(ierr) 生成一阶精度分裂的矩阵
!
!       6.call snes_fx1o(snes,x,f,null_int,ierr) 一阶精度右端项函数
!
!       7.call snes_rhs_fx_4ord(snes,x,f,null_int,ierr) 四阶精度右端项函数
!
! ------------------------------------------------------------------
    use penf, only: R_P
    use mod_parameters
    use mod_cubes
    use petsc
    implicit none
    public :: set_matfree_mat
    public :: mat_mult_4ord
    public :: init_mat_from_da
    public :: form_global_mat_2ord
    public :: form_global_mat_4ord
    public :: form_sub_mat_2_ord
    private
    type(lns_OP_point_type) :: Jor
    contains

    ! -----------------------------------------------------------------------------------------------------
    !   迭代格式一：标准线性求解器 Ax=b

    subroutine set_matfree_mat(comm,mat)
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

    subroutine form_global_mat_4ord(mat)

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

    end subroutine form_global_mat_4ord

    subroutine mat_set_boundary_conditions(mat,i,j,k)

        implicit none
        PetscScalar :: box(5,5),trans(5,5)
        MatStencil :: idxm(4,1),idxn(4,1)
        integer :: ic_index, jc_index
        integer :: lib, lie, ljb, lje
        integer,intent(in) :: i,j,k
        Mat,intent(inout) :: mat
        PetscErrorCode :: ierr
        integer :: li,lj

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
            idxn(MatStencil_i, 1)=i;idxn(MatStencil_j, 1)=j;idxn(MatStencil_k, 1)=k
            box=0.0d0;trans=0.0d0
            box=Jor%D-1.0d0*Jor%B
            trans=transpose(box)
            call MatSetValuesBlockedStencil(mat, 1, idxm, 1, idxn, trans, INSERT_VALUES, ierr)
            idxn(MatStencil_i, 1)=i;idxn(MatStencil_j, 1)=j+1;idxn(MatStencil_k, 1)=k
            box=0.0d0;trans=0.0d0
            box=Jor%B
            trans=transpose(box)
            call MatSetValuesBlockedStencil(mat, 1, idxm, 1, idxn, trans, INSERT_VALUES, ierr)
        elseif(i==(in-1) .and. j/=0 .and. j/=(jn-1))then
            lib=-2;lie=0
            ic_index=0
            idxm(MatStencil_i, 1)=i; idxm(MatStencil_j, 1)=j; idxm(MatStencil_k, 1)=k
            idxn(MatStencil_i, 1)=i-2;idxn(MatStencil_j, 1)=j;idxn(MatStencil_k, 1)=k
            box=0.0d0;box(1,1)=1.0d0;box(2,2)=1.0d0;box(3,3)=1.0d0;box(4,4)=1.0d0;box(5,5)=1.0d0
            box=1.0d0/2.0d0*box
            call MatSetValuesBlockedStencil(mat, 1, idxm, 1, idxn, box, INSERT_VALUES, ierr)
            idxn(MatStencil_i, 1)=i-1;idxn(MatStencil_j, 1)=j;idxn(MatStencil_k, 1)=k
            box=0.0d0;box(1,1)=1.0d0;box(2,2)=1.0d0;box(3,3)=1.0d0;box(4,4)=1.0d0;box(5,5)=1.0d0
            box=-2.0d0*box
            call MatSetValuesBlockedStencil(mat, 1, idxm, 1, idxn, box, INSERT_VALUES, ierr)
            idxn(MatStencil_i, 1)=i;idxn(MatStencil_j, 1)=j;idxn(MatStencil_k, 1)=k
            box=0.0d0;box(1,1)=1.0d0;box(2,2)=1.0d0;box(3,3)=1.0d0;box(4,4)=1.0d0;box(5,5)=1.0d0
            box=3.0d0/2.0d0*box
            call MatSetValuesBlockedStencil(mat, 1, idxm, 1, idxn, box, INSERT_VALUES, ierr)
        endif

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

        implicit none
        integer :: ic_index, jc_index, kc_index
        integer :: lib, lie, ljb, lje, lkb, lke
        PetscScalar :: box(5,5),trans(5,5)
        MatStencil :: idxm(4,1),idxn(4,1)
        integer,intent(in) :: i,j,k
        Mat,intent(inout) :: mat
        PetscErrorCode :: ierr
        integer :: li, lj, lk

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

    subroutine form_global_mat_2ord(mat)

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

    end subroutine form_global_mat_2ord

    subroutine mat_insert_values_2_ord(mat,i,j,k)

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

    subroutine form_sub_mat_2_ord(mat)  !!设置子块当地矩阵

        implicit none
        Mat,intent(inout) :: mat
        PetscErrorCode :: ierr
        integer :: i,j,k

        do k=kgs,kge
            do j=jgs,jge
                do i=igs,ige
                    if(i==igs .or. i==ige .or. j==jgs .or. j==jge .or. &
                    & ((k==kgs .or. k==kge) .and. k/=0 .and. k/=(kn-1)))then
                        call sub_mat_set_boundary_conditions(mat,i,j,k)  !!子矩阵边界点的\delta\hat\phi=0
                    else
                        call sub_mat_insert_values_2_ord(mat,i,j,k)  !!设置子矩阵内点的值
                    endif
                enddo
            enddo
        enddo

        ! select case(lns_mode)
        ! case(2)
        !     do k=kgs,kge
        !         do j=jgs,jge
        !             do i=igs,ige
        !                 if(i==igs .or. i==ige .or. j==jgs .or. j==jge)then
        !                     call sub_mat_set_boundary_conditions(mat,i,j,k)  !!子矩阵边界点的\delta\hat\phi=0
        !                 else
        !                     call sub_mat_insert_values_2_ord(mat,i,j,k)  !!设置子矩阵内点的值
        !                 endif
        !             enddo
        !         enddo
        !     enddo
        ! case(3)
        !     do k=kgs,kge
        !         do j=jgs,jge
        !             do i=igs,ige
        !                 if(i==igs .or. i==ige .or. j==jgs .or. j==jge .or. &
        !                 &  k==kgs .or. k==kge)then
        !                     call sub_mat_set_boundary_conditions(mat,i,j,k)  !!子矩阵边界点的\delta\hat\phi=0
        !                 else
        !                     call sub_mat_insert_values_2_ord(mat,i,j,k)  !!设置子矩阵内点的值
        !                 endif
        !             enddo
        !         enddo
        !     enddo
        ! end select

        call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)

    end subroutine form_sub_mat_2_ord

    subroutine sub_mat_set_boundary_conditions(mat,i,j,k)

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

end module mod_forming
