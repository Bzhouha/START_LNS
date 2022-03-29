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
!   for KSP :: Linear System Solvers
!
!       1.call dolphin_coming(comm) 免矩阵形式的矩阵生成函数
!
!           1).call dolphin_growing_up(A, X, F, ierr) 免矩阵需要的矩阵向量乘法函数
!
!           2).call dolphin_say_hi(comm) 进度标记：免矩阵生效中
!
!       2.call whale_coming(comm) 显式矩阵形式的矩阵生成函数
!
!           1).call whale_is_born(comm) 进度标记：前绪开始生成矩阵
!
!           2).call whale_growing_up() 填充数据函数
!
!               (1).call whale_eat_shrimps(i,j,k) 边界部分的填充数据函数
!
!               (2).call whale_eat_sardine(i,j,k) 内部部分的填充数据函数
!
!           3).call whale_say_hi(comm)
!
!       3.call set_right_hand_side(comm) 设置右边量，即边界。
!
!   for SNES :: Nonlinear Solvers
!
!       1.call shark_coming(comm) SNES所需的矩阵、数组分配
!
!       2.call Jacobi(snes,x,jac,B,null_int,ierr) 雅各比矩阵函数
!
!       3.call RHS_with_BC(snes,x,f,null_int,ierr) 右端项函数
!
!   call deallocate_bfinfo_and_metrics() 释放基本流类和度量系数数组内存
!
! ------------------------------------------------------------------
    use penf, only: R_P
    use mod_parameters
    use mod_flowtype
    use mod_cubes
    use petsc
    implicit none
    public :: dolphin_coming, whale_coming, shark_coming
    public :: shark_growing_up, RHS_with_BC, deallocate_bfinfo_and_metrics
    private
    type(lns_OP_point_type) :: Jor
    ! 注：这里是列优先，存储在内存中的样子是下面形式的转置，所以实际使用时需要将行列调换，如C1(li,c_index)。
    real(R_P), parameter, dimension(-2:2,-2:2) :: FDM_1nd_4ORD_CENTER=reshape( [&
        0.0d0       ,        0.0d0, -3.0d0/2.0d0, 2.0d0      ,  -1.0d0/2.0d0, &
        0.0d0       , -1.0d0/3.0d0, -1.0d0/2.0d0, 1.0d0      ,  -1.0d0/6.0d0, &
        1.0d0/12.0d0, -2.0d0/3.0d0,        0.0d0, 2.0d0/3.0d0, -1.0d0/12.0d0, &
        1.0d0/6.0d0 ,       -1.0d0,  1.0d0/2.0d0, 1.0d0/3.0d0,         0.0d0, &
        1.0d0/2.0d0 ,       -2.0d0,  3.0d0/2.0d0, 0.0d0      ,         0.0d0  &
        ], [5, 5])
    real(R_P), parameter, dimension(-2:2,-2:2) :: FDM_2nd_4ORD_CENTER=reshape( [&
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0, &
        0.0d0        ,       1.0d0,        -2.0d0,       1.0d0,         0.0d0, &
        -1.0d0/12.0d0, 4.0d0/3.0d0, -15.0d0/6.0d0, 4.0d0/3.0d0, -1.0d0/12.0d0, &
        0.0d0        ,       1.0d0,        -2.0d0,       1.0d0,         0.0d0, &
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0  &
        ], [5, 5])
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
    real(R_P), parameter, dimension(-2:2,-2:2) :: FDM_1nd_1ORD_Backward=reshape( [&
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0, &
        0.0d0        ,      -1.0d0,         1.0d0,       0.0d0,         0.0d0, &
        0.0d0        ,      -1.0d0,         1.0d0,       0.0d0,         0.0d0, &
        0.0d0        ,      -1.0d0,         1.0d0,       0.0d0,         0.0d0, &
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0  &
        ],[5,5])
    real(R_P), parameter, dimension(-2:2,-2:2) :: FDM_1nd_1ORD_Forward=reshape( [&
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0, &
        0.0d0        ,       0.0d0,        -1.0d0,       1.0d0,         0.0d0, &
        0.0d0        ,       0.0d0,        -1.0d0,       1.0d0,         0.0d0, &
        0.0d0        ,       0.0d0,        -1.0d0,       1.0d0,         0.0d0, &
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0  &
        ],[5,5])
    real(R_P), parameter :: delta_i(-2:2)=[0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]
    real(R_P), parameter :: delta_j(-2:2)=[0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]
    real(R_P), parameter :: delta_k(-2:2)=[0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]
    contains

    !   KSP :: Linear System Solvers

    subroutine dolphin_coming(comm)
        implicit none
        PetscInt,intent(in) :: comm
        PetscErrorCode :: ierr
        PetscInt :: ls
        call DMGetLocalVector(meshDA,tinkle_bell,ierr)
        call VecZeroEntries(tinkle_bell,ierr)
        call VecGetLocalSize(turtle,ls,ierr)
        call MatCreateShell(comm,ls,ls,PETSC_DETERMINE,PETSC_DETERMINE,PETSC_NULL_INTEGER,dolphin,ierr)
        call MatShellSetOperation(dolphin,MATOP_MULT,dolphin_growing_up,ierr)
        call MatAssemblyBegin(dolphin,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(dolphin,MAT_FINAL_ASSEMBLY,ierr)
        call VecDuplicate(turtle,RHS,ierr)
        call VecZeroEntries(RHS,ierr)
        call set_right_hand_side(comm)
        call PetscPrintf(comm,"\n   KSP :: Matrix-Free\n",ierr)
    end subroutine dolphin_coming

    subroutine dolphin_growing_up(A, X, F, ierr)
        implicit none
        PetscScalar,pointer :: fr(:,:,:,:),xr(:,:,:,:)
        PetscScalar,dimension(5,16) :: liy
        integer :: ic_index,jc_index,kc_index
        integer :: lib,lie,ljb,lje,lkb,lke
        PetscErrorCode :: ierr
        integer :: li,lj,lk
        integer :: i,j,k
        Vec :: X, F
        Mat :: A
        call DMGlobalToLocalBegin(meshDA,X,INSERT_VALUES,tinkle_bell,ierr)
        call DMGlobalToLocalEnd(meshDA,X,INSERT_VALUES,tinkle_bell,ierr)
        call DMDAVecGetArrayReadF90(meshDA,tinkle_bell,xr,ierr)
        call DMDAVecGetArrayF90(meshDA,F,fr,ierr)
        associate ( &
            coef_c4 =>FDM_1nd_4ORD_CENTER,   &
            coef_d4 =>FDM_2nd_4ORD_CENTER,   &
            coef_c4f=>FDM_1nd_4ORD_Forward,  &
            coef_c4b=>FDM_1nd_4ORD_Backward, &
            f   => fr,    x => xr, &
            G   => Jor%G, D => Jor%D, &
            A   => Jor%A, B => Jor%B, C => jor%C, &
            A_p => Jor%A_p, A_m => Jor%A_m, A_v => Jor%A_v, &
            B_p => Jor%B_p, B_m => Jor%B_m, B_v => Jor%B_v, &
            C_p => Jor%C_p, C_m => Jor%C_m, C_v => Jor%C_v, &
            Vxx => Jor%Vxx, Vyy => Jor%Vyy, Vzz => Jor%Vzz, &
            Vxy => Jor%Vxy, Vxz => Jor%Vxz, Vyz => Jor%Vyz)
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
                        liy=0.0d0
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

                        liy(:,1)=x(:,i,j,k) !D

                        do li=lib,lie
                            liy(:,2)  = liy(:,2)  + x(:,i+li,j,k)*coef_c4b(li,ic_index) ! A_p
                            liy(:,3)  = liy(:,3)  + x(:,i+li,j,k)*coef_c4f(li,ic_index) ! A_m
                            liy(:,4)  = liy(:,4)  + x(:,i+li,j,k)*coef_c4(li,ic_index)  ! A_v
                            liy(:,11) = liy(:,11) + x(:,i+li,j,k)*coef_d4(li,ic_index)  ! Vxx
                        enddo

                        do lj=ljb,lje
                            liy(:,5)  = liy(:,5)  + x(:,i,j+lj,k)*coef_c4b(lj,jc_index) ! B_p
                            liy(:,6)  = liy(:,6)  + x(:,i,j+lj,k)*coef_c4f(lj,jc_index) ! B_m
                            liy(:,7)  = liy(:,7)  + x(:,i,j+lj,k)*coef_c4(lj,jc_index)  ! B_v
                            liy(:,12) = liy(:,12) + x(:,i,j+lj,k)*coef_d4(lj,jc_index)  ! Vxx
                        enddo

                        do lk=lkb,lke
                            liy(:,8)  = liy(:,8)  + x(:,i,j,k+lk)*coef_c4b(lk,kc_index) ! C_p
                            liy(:,9)  = liy(:,9)  + x(:,i,j,k+lk)*coef_c4f(lk,kc_index) ! C_m
                            liy(:,10) = liy(:,10) + x(:,i,j,k+lk)*coef_c4(lk,kc_index)  ! C_v
                            liy(:,13) = liy(:,13) + x(:,i,j,k+lk)*coef_d4(lk,kc_index)  ! Vxx
                        enddo

                        do lj=ljb,lje
                            do li=lib,lie
                                liy(:,14)=liy(:,14)+x(:,i+li,j+lj,k)*coef_c4(li,ic_index)*coef_c4(lj,jc_index) ! Vxy
                            enddo
                        enddo

                        do lk=lkb,lke
                            do li=lib,lie
                                liy(:,15)=liy(:,15)+x(:,i+li,j,k+lk)*coef_c4(li,ic_index)*coef_c4(lk,kc_index) ! Vxz
                            enddo
                        enddo

                        do lk=lkb,lke
                            do lj=ljb,lje
                                liy(:,16)=liy(:,16)+x(:,i,j+lj,k+lk)*coef_c4(lj,jc_index)*coef_c4(lk,kc_index) ! Vyz
                            enddo
                        enddo

                        f(:,i,j,k)=matmul(D,liy(:,1)) &
                                & +matmul(A_p,liy(:,2))+matmul(A_m,liy(:,3))+matmul(A_v,liy(:,4)) &
                                & +matmul(B_p,liy(:,5))+matmul(B_m,liy(:,6))+matmul(B_v,liy(:,7)) &
                                & +matmul(C_p,liy(:,8))+matmul(C_m,liy(:,9))+matmul(C_v,liy(:,10)) &
                                & -matmul(Vxx,liy(:,11))-matmul(Vyy,liy(:,12))-matmul(Vzz,liy(:,13)) &
                                & -matmul(Vxy,liy(:,14))-matmul(Vxz,liy(:,15))-matmul(Vyz,liy(:,16))
                    endif
                enddo
            enddo
        enddo
        end associate
        call DMDAVecRestoreArrayReadF90(meshDA,tinkle_bell,xr,ierr)
        call DMDAVecRestoreArrayF90(meshDA,F,fr,ierr)
    end subroutine dolphin_growing_up

    subroutine whale_coming(comm)
        implicit none
        PetscInt,intent(in) :: comm
        PetscErrorCode :: ierr
        call DMCreateMatrix(meshDA, whale, ierr)
        call MatZeroEntries(whale,ierr)
        call whale_growing_up()
        call VecDuplicate(turtle,RHS,ierr)
        call VecZeroEntries(RHS,ierr)
        call set_right_hand_side(comm)
        call deallocate_bfinfo_and_metrics()
        call PetscPrintf(comm,"\n   KSP :: Matrix-Assembled\n",ierr)
    end subroutine whale_coming

    subroutine whale_growing_up()
        implicit none
        PetscErrorCode :: ierr
        integer :: i,j,k
        do k=ks,ke
            do j=js,je
                do i=is,ie
                    if(i==0 .or. i==(in-1) .or. j==0 .or. j==(jn-1))then
                        call whale_eat_shrimps(i,j,k)
                    else
                        call whale_eat_sardine(i,j,k)
                    endif
                enddo
            enddo
        enddo
        call MatAssemblyBegin(whale,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(whale,MAT_FINAL_ASSEMBLY,ierr)
    end subroutine whale_growing_up

    subroutine whale_eat_shrimps(i,j,k)
        implicit none
        PetscScalar :: box(5,5),trans(5,5)
        MatStencil :: idxm(4,1),idxn(4,1)
        integer :: ic_index, jc_index
        integer :: lib, lie, ljb, lje
        integer,intent(in) :: i,j,k
        integer :: li,lj,ii,jj
        PetscErrorCode :: ierr
        associate( &
            coef_c4f=>FDM_1nd_4ORD_Forward, &
            coef_c4b=>FDM_1nd_4ORD_Backward)
        if(i==0 .or. j==(jn-1))then
            box=0.0d0
            box(1,1)=1.0d0;box(2,2)=1.0d0;box(3,3)=1.0d0;box(4,4)=1.0d0;box(5,5)=1.0d0
            idxm(MatStencil_i, 1)=i; idxm(MatStencil_j, 1)=j; idxm(MatStencil_k, 1)=k
            idxn(MatStencil_i, 1)=i; idxn(MatStencil_j, 1)=j; idxn(MatStencil_k, 1)=k
            call MatSetValuesBlockedStencil(whale, 1, idxm, 1, idxn, box, INSERT_VALUES, ierr)
        elseif(j==0 .and. i/=0)then
            ljb=0;lje=1
            jc_index=1
            idxm(MatStencil_i, 1)=i; idxm(MatStencil_j, 1)=j; idxm(MatStencil_k, 1)=k
            call Jor%get_adorned_cubes(i,j,k)
            do jj=1,5
                do ii=2,5
                    Jor%D(ii,jj)=0.0d0
                    Jor%B(ii,jj)=0.0d0
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
                call MatSetValuesBlockedStencil(whale, 1, idxm, 1, idxn, trans, INSERT_VALUES, ierr)
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
                call MatSetValuesBlockedStencil(whale, 1, idxm, 1, idxn, box, INSERT_VALUES, ierr)
            enddo
        endif
        end associate
    end subroutine whale_eat_shrimps

    subroutine whale_eat_sardine(i,j,k)
        implicit none
        integer :: ic_index, jc_index, kc_index
        integer :: lib, lie, ljb, lje, lkb, lke
        PetscScalar :: box(5,5),trans(5,5)
        MatStencil :: idxm(4,1),idxn(4,1)
        integer,intent(in) :: i,j,k
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
        coef_c4 =>FDM_1nd_4ORD_CENTER,   &
        coef_d4 =>FDM_2nd_4ORD_CENTER,   &
        coef_c4f=>FDM_1nd_4ORD_Forward,  &
        coef_c4b=>FDM_1nd_4ORD_Backward, &
        G   => Jor%G,   D   => Jor%D,    &
        A_p => Jor%A_p, A_m => Jor%A_m, A_v => Jor%A_v, &
        B_p => Jor%B_p, B_m => Jor%B_m, B_v => Jor%B_v, &
        C_p => Jor%C_p, C_m => Jor%C_m, C_v => Jor%C_v, &
        Vxx => Jor%Vxx, Vyy => Jor%Vyy, Vzz => Jor%Vzz, &
        Vxy => Jor%Vxy, Vxz => Jor%Vxz, Vyz => Jor%Vyz)
        do lk = lkb, lke
            do lj = ljb, lje
                do li = lib, lie
                    idxn=0;box=0.0d0;trans=0.0d0
                    idxn(MatStencil_i, 1)=i+li
                    idxn(MatStencil_j, 1)=j+lj
                    idxn(MatStencil_k, 1)=k+lk
                    !if(idxn(MatStencil_k, 1)>(kn-1)) idxn(MatStencil_k, 1)=idxn(MatStencil_k, 1)-kn
                    !if(idxn(MatStencil_k, 1)<0)  idxn(MatStencil_k, 1)=idxn(MatStencil_k, 1)+kn  !! kn为展向一个周期的点数
                    box=delta_i(li)*delta_j(lj)*delta_k(lk)*D + &
                        delta_j(lj)*delta_k(lk)*A_v*coef_c4(li,ic_index)+ &
                        delta_j(lj)*delta_k(lk)*A_m*coef_c4f(li,ic_index)+ &
                        delta_j(lj)*delta_k(lk)*A_p*coef_c4b(li,ic_index)+ &
                        delta_i(li)*delta_k(lk)*B_v*coef_c4(lj,jc_index)+ &
                        delta_i(li)*delta_k(lk)*B_m*coef_c4f(lj,jc_index)+ &
                        delta_i(li)*delta_k(lk)*B_p*coef_c4b(lj,jc_index)+ &
                        delta_i(li)*delta_j(lj)*C_v*coef_c4(lk,kc_index)+ &
                        delta_i(li)*delta_j(lj)*C_m*coef_c4f(lk,kc_index)+ &
                        delta_i(li)*delta_j(lj)*C_p*coef_c4b(lk,kc_index)- &
                        delta_j(lj)*delta_k(lk)*Vxx*coef_d4(li,ic_index)- &
                        delta_i(li)*delta_k(lk)*Vyy*coef_d4(lj,jc_index)- &
                        delta_i(li)*delta_j(lj)*Vzz*coef_d4(lk,kc_index)- &
                        delta_k(lk)*Vxy*coef_c4(li,ic_index)*coef_c4(lj,jc_index)- &
                        delta_j(lj)*Vxz*coef_c4(li,ic_index)*coef_c4(lk,kc_index)- &
                        delta_i(li)*Vyz*coef_c4(lj,jc_index)*coef_c4(lk,kc_index)
                    trans=transpose(box)
                    call MatSetValuesBlockedStencil(whale, 1, idxm, 1, idxn, trans, INSERT_VALUES, ierr)
                end do
            end do
        end do
        end associate
    end subroutine whale_eat_sardine

    subroutine set_right_hand_side(comm)
        implicit none
        PetscScalar,pointer :: RHS_array(:,:,:,:)
        PetscInt,intent(in) :: comm
        PetscErrorCode :: ierr
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

    !   SNES :: Nonlinear Solvers

    subroutine shark_coming(comm)
        implicit none
        PetscInt,intent(in) :: comm
        PetscErrorCode :: ierr
        call shark_is_born(comm)
        call DMGetLocalVector(meshDA,tinkle_bell,ierr)
        call VecZeroEntries(tinkle_bell,ierr)
        call DMCreateMatrix(meshDA,shark,ierr)
        call MatZeroEntries(shark,ierr)
        call shark_say_hi(comm)
    end subroutine shark_coming

    subroutine shark_is_born(comm)
        implicit none
        PetscInt,intent(in) :: comm
        PetscErrorCode :: ierr
        call PetscPrintf(comm," -----------------------------------\n",ierr)
        call PetscPrintf(comm,"          Jacobi :: be born.\n",ierr)
        call PetscPrintf(comm," -----------------------------------\n",ierr)
    end subroutine shark_is_born

    subroutine shark_growing_up(snes,x,jac,B,null_int,ierr)
        implicit none
        PetscErrorCode :: ierr
        integer :: null_int(*)
        integer :: i,j,k
        Mat :: jac, B
        SNES :: snes
        Vec :: x
        do k=ks,ke
            do j=js,je
                do i=is,ie
                    if(i==0 .or. i==(in-1) .or. j==0 .or. j==(jn-1))then
                        call shark_eat_shrimps(jac,i,j,k)
                    else
                        call shark_eat_sardine(jac,i,j,k)
                    endif
                enddo
            enddo
        enddo
        call MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr)
    end subroutine shark_growing_up

    subroutine shark_eat_shrimps(jac,i,j,k)
        implicit none
        PetscScalar :: box(5,5),trans(5,5)
        MatStencil :: idxm(4,1),idxn(4,1)
        integer :: ic_index, jc_index
        integer :: lib, lie, ljb, lje
        integer,intent(in) :: i,j,k
        integer :: li,lj,ii,jj
        PetscErrorCode :: ierr
        Mat :: jac
        associate( &
            coef_c4f=>FDM_1nd_4ORD_Forward, &
            coef_c4b=>FDM_1nd_4ORD_Backward)
        if(i==0 .or. j==(jn-1))then
            box=0.0d0
            box(1,1)=1.0d0;box(2,2)=1.0d0;box(3,3)=1.0d0;box(4,4)=1.0d0;box(5,5)=1.0d0
            idxm(MatStencil_i, 1)=i; idxm(MatStencil_j, 1)=j; idxm(MatStencil_k, 1)=k
            idxn(MatStencil_i, 1)=i; idxn(MatStencil_j, 1)=j; idxn(MatStencil_k, 1)=k
            call MatSetValuesBlockedStencil(jac, 1, idxm, 1, idxn, box, INSERT_VALUES, ierr)
        elseif(j==0 .and. i/=0)then
            ljb=0;lje=1
            jc_index=1
            idxm(MatStencil_i, 1)=i; idxm(MatStencil_j, 1)=j; idxm(MatStencil_k, 1)=k
            call Jor%get_adorned_cubes(i,j,k)
            do jj=1,5
                do ii=2,5
                    Jor%D(ii,jj)=0.0d0
                    Jor%B(ii,jj)=0.0d0
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
                call MatSetValuesBlockedStencil(jac, 1, idxm, 1, idxn, trans, INSERT_VALUES, ierr)
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
                call MatSetValuesBlockedStencil(jac, 1, idxm, 1, idxn, box, INSERT_VALUES, ierr)
            enddo
        endif
        end associate
    end subroutine shark_eat_shrimps

    subroutine shark_eat_sardine(jac,i,j,k)
        implicit none
        integer :: ic_index, jc_index, kc_index
        integer :: lib, lie, ljb, lje, lkb, lke
        PetscScalar :: box(5,5),trans(5,5)
        MatStencil :: idxm(4,1),idxn(4,1)
        integer,intent(in) :: i,j,k
        PetscErrorCode :: ierr
        integer :: li,lj,lk
        Mat :: jac
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
        coef_c4 =>FDM_1nd_4ORD_CENTER,   &
        coef_d4 =>FDM_2nd_4ORD_CENTER,   &
        coef_c1f=>FDM_1nd_1ORD_Forward,  &
        coef_c1b=>FDM_1nd_1ORD_Backward, &
        G   => Jor%G,   D   => Jor%D,    &
        A_p => Jor%A_p, A_m => Jor%A_m, A_v => Jor%A_v, &
        B_p => Jor%B_p, B_m => Jor%B_m, B_v => Jor%B_v, &
        C_p => Jor%C_p, C_m => Jor%C_m, C_v => Jor%C_v, &
        Vxx => Jor%Vxx, Vyy => Jor%Vyy, Vzz => Jor%Vzz, &
        Vxy => Jor%Vxy, Vxz => Jor%Vxz, Vyz => Jor%Vyz)
        do lk = lkb, lke
            do lj = ljb, lje
                do li = lib, lie
                    idxn=0;box=0.0d0;trans=0.0d0
                    idxn(MatStencil_i, 1)=i+li
                    idxn(MatStencil_j, 1)=j+lj
                    idxn(MatStencil_k, 1)=k+lk
                    box=delta_i(li)*delta_j(lj)*delta_k(lk)*D + &
                        delta_j(lj)*delta_k(lk)*A_v*coef_c4(li,ic_index)+ &
                        delta_j(lj)*delta_k(lk)*A_m*coef_c1f(li,ic_index)+ &
                        delta_j(lj)*delta_k(lk)*A_p*coef_c1b(li,ic_index)+ &
                        delta_i(li)*delta_k(lk)*B_v*coef_c4(lj,jc_index)+ &
                        delta_i(li)*delta_k(lk)*B_m*coef_c1f(lj,jc_index)+ &
                        delta_i(li)*delta_k(lk)*B_p*coef_c1b(lj,jc_index)+ &
                        delta_i(li)*delta_j(lj)*C_v*coef_c4(lk,kc_index)+ &
                        delta_i(li)*delta_j(lj)*C_m*coef_c1f(lk,kc_index)+ &
                        delta_i(li)*delta_j(lj)*C_p*coef_c1b(lk,kc_index)- &
                        delta_j(lj)*delta_k(lk)*Vxx*coef_d4(li,ic_index)- &
                        delta_i(li)*delta_k(lk)*Vyy*coef_d4(lj,jc_index)- &
                        delta_i(li)*delta_j(lj)*Vzz*coef_d4(lk,kc_index)- &
                        delta_k(lk)*Vxy*coef_c4(li,ic_index)*coef_c4(lj,jc_index)- &
                        delta_j(lj)*Vxz*coef_c4(li,ic_index)*coef_c4(lk,kc_index)- &
                        delta_i(li)*Vyz*coef_c4(lj,jc_index)*coef_c4(lk,kc_index)
                    trans=transpose(box)
                    call MatSetValuesBlockedStencil(jac, 1, idxm, 1, idxn, trans, INSERT_VALUES, ierr)
                end do
            end do
        end do
        end associate
    end subroutine shark_eat_sardine

    subroutine shark_say_hi(comm)
        implicit none
        PetscInt,intent(in) :: comm
        PetscErrorCode :: ierr
        call PetscPrintf(comm," -----------------------------------\n",ierr)
        call PetscPrintf(comm,"          Jacobi :: be grown.\n",ierr)
        call PetscPrintf(comm," -----------------------------------\n",ierr)
    end subroutine shark_say_hi

    subroutine RHS_with_BC(snes,x,f,null_int,ierr)
        implicit none
        PetscScalar,pointer :: f_local(:,:,:,:),x_local(:,:,:,:)
        complex(R_P), dimension(5,9) :: cab
        PetscScalar :: M1(5,5),M2(5,5)
        PetscErrorCode :: ierr
        integer :: i1,i2,j1,j2
        integer :: null_int(*)
        integer :: i,j,k
        SNES :: snes
        Vec :: x, f
    end subroutine RHS_with_BC

    subroutine deallocate_bfinfo_and_metrics()
        ! 免矩阵版本需在全部计算结束之后才可释放内存
        implicit none
        deallocate(bf)
        deallocate(xi_x,xi_y,xi_z)
        deallocate(eta_x,eta_y,eta_z)
        deallocate(phi_x,phi_y,phi_z)
        deallocate(xi_xx,xi_yy,xi_zz)
        deallocate(eta_xx,eta_yy,eta_zz)
        deallocate(phi_xx,phi_yy,phi_zz)
        deallocate(xi_xy,xi_xz,xi_yz)
        deallocate(eta_xy,eta_yz,eta_xz)
        deallocate(phi_xy,phi_yz,phi_xz)
    end subroutine deallocate_bfinfo_and_metrics
end module mod_forming
