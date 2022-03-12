!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_forming
! ------------------------------------------------------------------
!
!  这个模块:
!
!   在线性求解系统 KSP 框架下，生成最终的大矩阵；
!   在非线性求解系统 SNES 框架下，生成雅各比矩阵和右端项函数。
!
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
!           1).call we_hear_a_sound(comm) 进度标记：前绪开始生成矩阵
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
!
!   for SNES :: Nonlinear Solvers  
!
!       1.call Jacobi(snes,x,jac,B,null_int,ierr) 雅各比矩阵函数
!
!       2.call RHS_with_BC(snes,x,f,null_int,ierr) 右端项函数
!
! ------------------------------------------------------------------
    use penf, only: R_P
    use mod_parameters
    use mod_flowtype
    use mod_cubes
    use petsc
    implicit none
    public :: dolphin_coming, whale_coming, set_right_hand_side
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
        call VecGetLocalSize(turtle,ls,ierr)
        call MatCreateShell(comm,ls,ls,PETSC_DETERMINE,PETSC_DETERMINE,PETSC_NULL_INTEGER,Dolphin,ierr)
        call MatShellSetOperation(Dolphin,MATOP_MULT,dolphin_growing_up,ierr)
        call MatAssemblyBegin(Dolphin,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(Dolphin,MAT_FINAL_ASSEMBLY,ierr)
        call dolphin_say_hi(comm)
    end subroutine dolphin_coming

    subroutine dolphin_growing_up(A, X, F, ierr)
        use mod_mftools
        implicit none
        PetscScalar,pointer :: F_r(:,:,:,:),X_r(:,:,:,:)
        complex(R_P), dimension(5,15) :: cab
        PetscScalar :: M1(5,5),M2(5,5)
        integer :: i1,i2,j1,j2
        PetscErrorCode :: ierr
        integer :: i,j,k
        Vec :: Clams
        Vec :: X, F
        Mat :: A 
        call DMGetLocalVector(meshDA,Clams,ierr)
        call VecZeroEntries(Clams,ierr)
        call DMGlobalToLocalBegin(meshDA,X,INSERT_VALUES,Clams,ierr)
        call DMGlobalToLocalEnd(meshDA,X,INSERT_VALUES,Clams,ierr)
        call DMDAVecGetArrayReadF90(meshDA,Clams,X_r,ierr)
        call DMDAVecGetArrayF90(meshDA,F,F_r,ierr)
        associate ( F => F_r, x => X_r, &
            G   => Jor%G, D => Jor%D, &
            A_p => Jor%A_p, A_m => Jor%A_m, A_v => Jor%A_v, &
            B_p => Jor%B_p, B_m => Jor%B_m, B_v => Jor%B_v, &
            C_p => Jor%C_p, C_m => Jor%C_m, C_v => Jor%C_v, &
            Vxx => Jor%Vxx, Vyy => Jor%Vyy, Vzz => Jor%Vzz, &
            Vxy => Jor%Vxy, Vxz => Jor%Vxz, Vyz => Jor%Vyz) 
        do k=ks,ke
            do j=js,je
                do i=is,ie
                    if(i==0 .or. j==(jn-1))then
                        F(:,i,j,k)=x(:,i,j,k)
                    elseif(i==(in-1) .and. j/=0 .and. j/=(jn-1))then
                        F(:,i,j,k)=(x(:,i-2,j,k)-4*x(:,i-1,j,k)+3*x(:,i,j,k))/2.0d0
                    elseif(j==0 .and. i/=0)then
                        M1=0.0d0;M2=0.0d0 
                        M1(1,3)=bf(i,j,k)%BF%rho
                        M2(1,1)=bf(i,j,k)%BFDy%y-cmplx(0.0d0,1.0d0,R_P)*Omega
                        M2(2,2)=1.0d0;M2(3,3)=1.0d0;M2(4,4)=1.0d0;M2(5,5)=1.0d0
                        F(:,i,j,k)=matmul(M2,x(:,i,j,k))+matmul(M1,x(:,i,j+1,k))-matmul(M1,x(:,i,j,k))
                    elseif(i>1 .and. j>1 .and. i<(in-2) .and. j<(jn-2))then
                        cab=0.0d0
                        call Jor%get_adorned_cubes(i,j,k)
                        call f3d1r(cab(:,1),x(:,i-2:i,j,k)) ! A1p
                        call f3d1f(cab(:,2),x(:,i:i+2,j,k)) ! A1m
                        call f5d1(cab(:,3),x(:,i-2:i+2,j,k)) ! A2
                        call f3d1r(cab(:,4),x(:,i,j-2:j,k)) ! B1p
                        call f3d1f(cab(:,5),x(:,i,j:j+2,k)) ! B1m
                        call f5d1(cab(:,6),x(:,i,j-2:j+2,k)) ! B2
                        call f3d1r(cab(:,7),x(:,i,j,k-2:k)) ! C1p
                        call f3d1f(cab(:,8),x(:,i,j,k:k+2)) ! C1m
                        call f5d1(cab(:,9),x(:,i,j,k-2:k+2)) ! C2
                        call f5d11(cab(:,10),x(:,i-2:i+2,j,k)) ! Vxx
                        call f5d11(cab(:,11),x(:,i,j-2:j+2,k)) ! Vyy
                        call f5d11(cab(:,12),x(:,i,j,k-2:k+2)) ! Vzz
                        call f5d12(cab(:,13),x(:,i-2:i+2,j-2:j+2,k)) ! Vxy
                        call f5d12(cab(:,14),x(:,i-2:i+2,j,k-2:k+2)) ! Vxz
                        call f5d12(cab(:,15),x(:,i,j-2:j+2,k-2:k+2)) ! Vyz
                        F(:,i,j,k) = matmul(A_p,cab(:,1))+matmul(A_m,cab(:,2))+matmul(A_v,cab(:,3)) &
                                    +matmul(B_p,cab(:,4))+matmul(B_m,cab(:,5))+matmul(B_v,cab(:,6)) &
                                    +matmul(C_p,cab(:,7))+matmul(C_m,cab(:,8))+matmul(C_v,cab(:,9)) &
                                    +matmul(D,x(:,i,j,k)) &
                                    -matmul(Vxx,cab(:,10))-matmul(Vyy,cab(:,11))-matmul(Vzz,cab(:,12)) &
                                    -matmul(Vxy,cab(:,13))-matmul(Vxz,cab(:,14))-matmul(Vyz,cab(:,15))
                    else 
                        cab=0.0d0
                        call Jor%get_adorned_cubes(i,j,k)
                        call f4_index(i1,i2,j1,j2,i,j)
                        call f2d1r(cab(:,1),x(:,i-1:i,j,k)) ! A1p
                        call f2d1f(cab(:,2),x(:,i:i+1,j,k)) ! A1m
                        call f4d1(cab(:,3),x(:,i1:i2,j,k),i,j,k,1) ! A2
                        call f2d1r(cab(:,4),x(:,i,j-1:j,k)) ! B1p
                        call f2d1f(cab(:,5),x(:,i,j:j+1,k)) ! B1m
                        call f4d1(cab(:,6),x(:,i,j1:j2,k),i,j,k,2) ! B2
                        call f3d1r(cab(:,7),x(:,i,j,k-2:k)) ! C1p
                        call f3d1f(cab(:,8),x(:,i,j,k:k+2)) ! C1m
                        call f5d1(cab(:,9),x(:,i,j,k-2:k+2)) ! C2
                        call f4d11(cab(:,10),x(:,i-1:i+1,j,k)) ! Vxx
                        call f4d11(cab(:,11),x(:,i,j-1:j+1,k)) ! Vyy
                        call f5d11(cab(:,12),x(:,i,j,k-2:k+2)) ! Vzz
                        call f4d12(cab(:,13),x(:,i1:i2,j1:j2,k),i,j,k,12) ! Vxy
                        call f45d12(cab(:,14),x(:,i1:i2,j,k-2:k+2),i,j,k,1) ! Vxz
                        call f45d12(cab(:,15),x(:,i,j1:j2,k-2:k+2),i,j,k,2) ! Vyz
                        F(:,i,j,k) = matmul(A_p,cab(:,1))+matmul(A_m,cab(:,2))+matmul(A_v,cab(:,3)) &
                                    +matmul(B_p,cab(:,4))+matmul(B_m,cab(:,5))+matmul(B_v,cab(:,6)) &
                                    +matmul(C_p,cab(:,7))+matmul(C_m,cab(:,8))+matmul(C_v,cab(:,9)) &
                                    +matmul(D,x(:,i,j,k)) &
                                    -matmul(Vxx,cab(:,10))-matmul(Vyy,cab(:,11))-matmul(Vzz,cab(:,12)) &
                                    -matmul(Vxy,cab(:,13))-matmul(Vxz,cab(:,14))-matmul(Vyz,cab(:,15))
                    endif 
                enddo
            enddo
        enddo
        end associate
        call DMDAVecRestoreArrayReadF90(meshDA,Clams,X_r,ierr)
        call DMDAVecRestoreArrayF90(meshDA,F,F_r,ierr)
    end subroutine dolphin_growing_up

    subroutine dolphin_say_hi(comm)
        implicit none
        PetscInt,intent(in) :: comm
        PetscErrorCode :: ierr
        call PetscPrintf(comm," -----------------------------------\n",ierr)
        call PetscPrintf(comm,"             免矩阵生效中             \n",ierr)
        call PetscPrintf(comm," -----------------------------------\n",ierr)
    end subroutine dolphin_say_hi

    subroutine whale_coming(comm)
        implicit none
        PetscInt,intent(in) :: comm
        PetscErrorCode :: ierr 
        call we_hear_a_sound(comm)
        ! call DMSetMatType(meshDA, MATBAIJ, ierr)
        call DMCreateMatrix(meshDA, Whale, ierr)
        call MatZeroEntries(Whale,ierr)
        call whale_growing_up()
        call whale_say_hi(comm)
    end subroutine whale_coming

    subroutine we_hear_a_sound(comm)
        implicit none
        PetscInt,intent(in) :: comm
        PetscErrorCode :: ierr
        call PetscPrintf(comm," -----------------------------------\n",ierr)
        call PetscPrintf(comm,"             开始生成矩阵             \n",ierr)
        call PetscPrintf(comm," -----------------------------------\n",ierr)
    end subroutine we_hear_a_sound
 
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
        call MatAssemblyBegin(Whale,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(Whale,MAT_FINAL_ASSEMBLY,ierr)
    end subroutine whale_growing_up

    subroutine whale_eat_shrimps(i,j,k)
        implicit none
        PetscScalar :: box(5,5),trans(5,5)
        MatStencil :: idxm(4,1),idxn(4,1)
        integer :: ic_index, jc_index
        integer :: lib, lie, ljb, lje
        integer,intent(in) :: i,j,k
        PetscErrorCode :: ierr
        integer :: ii,jj
        integer :: li,lj
        associate( &
            coef_c1f=>FDM_1nd_4ORD_Forward, &
            coef_c1b=>FDM_1nd_4ORD_Backward)
        if(i==0 .or. j==(jn-1))then
            box=0.0d0
            box(1,1)=1.0d0;box(2,2)=1.0d0;box(3,3)=1.0d0;box(4,4)=1.0d0;box(5,5)=1.0d0
            idxm(MatStencil_i, 1)=i; idxm(MatStencil_j, 1)=j; idxm(MatStencil_k, 1)=k
            idxn(MatStencil_i, 1)=i; idxn(MatStencil_j, 1)=j; idxn(MatStencil_k, 1)=k
            call MatSetValuesBlockedStencil(Whale, 1, idxm, 1, idxn, box, INSERT_VALUES, ierr)
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
                box=delta_j(lj)*Jor%D+coef_c1f(lj,jc_index)*Jor%B
                trans=transpose(box)
                call MatSetValuesBlockedStencil(Whale, 1, idxm, 1, idxn, trans, INSERT_VALUES, ierr)
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
                box=coef_c1b(li,ic_index)*box
                call MatSetValuesBlockedStencil(Whale, 1, idxm, 1, idxn, box, INSERT_VALUES, ierr)
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
            case(0)
                lkb=0; lke=0; kc_index=0
            case(1)
                lkb=-2; lke=2; kc_index=0
        end select
        idxm=0; 
        idxm(MatStencil_i, 1)=i; idxm(MatStencil_j, 1)=j; idxm(MatStencil_k, 1)=k
        associate( &
        coef_c1 =>FDM_1nd_4ORD_CENTER,   &
        coef_d1 =>FDM_2nd_4ORD_CENTER,   &
        coef_c1b=>FDM_1nd_4ORD_Backward, &
        coef_c1f=>FDM_1nd_4ORD_Forward,  &
          G => Jor%G,     D => Jor%D,    &
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
                        delta_j(lj)*delta_k(lk)*A_v*coef_c1(li,ic_index)+ &
                        delta_j(lj)*delta_k(lk)*A_m*coef_c1f(li,ic_index)+ &
                        delta_j(lj)*delta_k(lk)*A_p*coef_c1b(li,ic_index)+ &
                        delta_i(li)*delta_k(lk)*B_v*coef_c1(lj,jc_index)+ &
                        delta_i(li)*delta_k(lk)*B_m*coef_c1f(lj,jc_index)+ &
                        delta_i(li)*delta_k(lk)*B_p*coef_c1b(lj,jc_index)+ &
                        delta_i(li)*delta_j(lj)*C_v*coef_c1(lk,kc_index)+ &
                        delta_i(li)*delta_j(lj)*C_m*coef_c1f(lk,kc_index)+ &
                        delta_i(li)*delta_j(lj)*C_p*coef_c1b(lk,kc_index)- &
                        delta_j(lj)*delta_k(lk)*Vxx*coef_d1(li,ic_index)- &
                        delta_i(li)*delta_k(lk)*Vyy*coef_d1(lj,jc_index)- &
                        delta_i(li)*delta_j(lj)*Vzz*coef_d1(lk,kc_index)- &
                        delta_k(lk)*Vxy*coef_c1(li,ic_index)*coef_c1(lj,jc_index)- &
                        delta_j(lj)*Vxz*coef_c1(li,ic_index)*coef_c1(lk,kc_index)- &
                        delta_i(li)*Vyz*coef_c1(lj,jc_index)*coef_c1(lk,kc_index)
                    trans=transpose(box)
                    call MatSetValuesBlockedStencil(Whale, 1, idxm, 1, idxn, trans, INSERT_VALUES, ierr)
                end do
            end do
        end do
        end associate
    end subroutine whale_eat_sardine

    subroutine whale_say_hi(comm)
        implicit none
        PetscInt,intent(in) :: comm
        PetscErrorCode :: ierr
        call PetscPrintf(comm," -----------------------------------\n",ierr)
        call PetscPrintf(comm,"             矩阵生成结束             \n",ierr)
        call PetscPrintf(comm," -----------------------------------\n",ierr)
    end subroutine whale_say_hi

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

    subroutine Jacobi(snes,x,jac,B,null_int,ierr)
        implicit none
        PetscErrorCode :: ierr
        integer :: null_int(*)
        Mat :: jac, B
        SNES :: snes  
        Vec :: x 
    end subroutine Jacobi

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

        ! tinkle_bell should be already allocated before calling this function

        call DMGlobalToLocalBegin(meshDA,x,INSERT_VALUES,tinkle_bell,ierr)
        call DMGlobalToLocalEnd(meshDA,x,INSERT_VALUES,tinkle_bell,ierr)

        call DMDAVecGetArrayReadF90(meshDA,tinkle_bell,x_local,ierr)
        call DMDAVecGetArrayF90(meshDA,f,f_local,ierr)

        associate ( f => f_local, x => x_local, &
            G   => Jor%G, D => Jor%D, &
            A   => Jor%A, B => Jor%B, C => Jor%C, &
            Vxx => Jor%Vxx, Vyy => Jor%Vyy, Vzz => Jor%Vzz, &
            Vxy => Jor%Vxy, Vxz => Jor%Vxz, Vyz => Jor%Vyz) 
            do k=ks,ke
                do j=js,je
                    do i=is,ie
                        if(i==0)then
                            f(:,i,j,k)=disturb(:,j,k)-x(:,i,j,k)
                        elseif(j==(jn-1))then
                            f(:,i,j,k)=(-1.0d0)*x(:,i,j,k)
                        elseif(i==(in-1) .and. j/=0 .and. j/=(jn-1))then
                            f(:,i,j,k)=(-1.0d0)*(x(:,i-2,j,k)-4*x(:,i-1,j,k)+3*x(:,i,j,k))/2.0d0
                        elseif(j==0 .and. i/=0)then
                            M1=0.0d0;M2=0.0d0 
                            M1(1,3)=bf(i,j,k)%BF%rho
                            M2(1,1)=bf(i,j,k)%BFDy%y-cmplx(0.0d0,1.0d0,R_P)*Omega
                            M2(2,2)=1.0d0;M2(3,3)=1.0d0;M2(4,4)=1.0d0;M2(5,5)=1.0d0
                            f(:,i,j,k)=(-1.0d0)*(matmul(M2,x(:,i,j,k))+matmul(M1,x(:,i,j+1,k))-matmul(M1,x(:,i,j,k)))
                        elseif(i>1 .and. j>1 .and. i<(in-2) .and. j<(jn-2))then
                            cab=0.0d0
                            call Jor%get_adorned_cubes(i,j,k)
                            call f5d1(cab(:,1),x(:,i-2:i+2,j,k)) ! A
                            call f5d1(cab(:,2),x(:,i,j-2:j+2,k)) ! B
                            call f5d1(cab(:,3),x(:,i,j,k-2:k+2)) ! C
                            call f5d11(cab(:,4),x(:,i-2:i+2,j,k)) ! Vxx
                            call f5d11(cab(:,5),x(:,i,j-2:j+2,k)) ! Vyy
                            call f5d11(cab(:,6),x(:,i,j,k-2:k+2)) ! Vzz
                            call f5d12(cab(:,7),x(:,i-2:i+2,j-2:j+2,k)) ! Vxy
                            call f5d12(cab(:,8),x(:,i-2:i+2,j,k-2:k+2)) ! Vxz
                            call f5d12(cab(:,9),x(:,i,j-2:j+2,k-2:k+2)) ! Vyz
                            f(:,i,j,k) = (-1.0d0)*(matmul(A,cab(:,1))+matmul(B,cab(:,2))+matmul(C,cab(:,3)) &
                                        +matmul(D,x(:,i,j,k)) &
                                        -matmul(Vxx,cab(:,4))-matmul(Vyy,cab(:,5))-matmul(Vzz,cab(:,6)) &
                                        -matmul(Vxy,cab(:,7))-matmul(Vxz,cab(:,8))-matmul(Vyz,cab(:,9)))
                        else 
                            cab=0.0d0
                            call Jor%get_adorned_cubes(i,j,k)
                            call f4_index(i1,i2,j1,j2,i,j)
                            call f4d1(cab(:,1),x(:,i1:i2,j,k),i,j,k,1) ! A
                            call f4d1(cab(:,2),x(:,i,j1:j2,k),i,j,k,2) ! B
                            call f5d1(cab(:,3),x(:,i,j,k-2:k+2)) ! C
                            call f4d11(cab(:,4),x(:,i-1:i+1,j,k)) ! Vxx
                            call f4d11(cab(:,5),x(:,i,j-1:j+1,k)) ! Vyy
                            call f5d11(cab(:,6),x(:,i,j,k-2:k+2)) ! Vzz
                            call f4d12(cab(:,7),x(:,i1:i2,j1:j2,k),i,j,k,12) ! Vxy
                            call f45d12(cab(:,8),x(:,i1:i2,j,k-2:k+2),i,j,k,1) ! Vxz
                            call f45d12(cab(:,9),x(:,i,j1:j2,k-2:k+2),i,j,k,2) ! Vyz
                            f(:,i,j,k) = (-1.0d0)*(matmul(A,cab(:,1))+matmul(B,cab(:,2))+matmul(C,cab(:,3)) &
                                        +matmul(D,x(:,i,j,k)) &
                                        -matmul(Vxx,cab(:,4))-matmul(Vyy,cab(:,5))-matmul(Vzz,cab(:,6)) &
                                        -matmul(Vxy,cab(:,7))-matmul(Vxz,cab(:,8))-matmul(Vyz,cab(:,9)))
                        endif
                    enddo
                enddo
            enddo
        end associate

        call DMDAVecRestoreArrayReadF90(meshDA,tinkle_bell,x_local,ierr)
        call DMDAVecRestoreArrayF90(meshDA,f,f_local,ierr)

    end subroutine RHS_with_BC

end module mod_forming
