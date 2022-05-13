!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_iterate
    use penf, only: R_P
    use mod_parameters
    use mod_cubes
    use petsc
    implicit none
    public :: duplicate_vec
    public :: linear_rhs
    public :: snes_rhs_fx_center
    public :: snes_rhs_fx_4ord
    public :: snes_converged_test
    public :: push_bc
    public :: fx_rhs_Ax_4ord
    public :: fx_rhs_Ax_center
    public :: set_medDA
    public :: set_subDA
    public :: init_sub_vecs
    public :: get_subf
    public :: merge_res

    private
    type(lns_OP_point_type) :: Jor

    contains

    subroutine duplicate_vec(x,r)
        implicit none
        Vec,intent(inout) :: r
        PetscErrorCode :: ierr
        Vec,intent(in) :: x

        call VecDuplicate(x,r,ierr)
        call VecZeroEntries(r,ierr)
    end subroutine duplicate_vec

    subroutine linear_rhs(comm,r)

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
                        r_array(:,0,j,k)=inlet(:,j,k)
                    enddo
                enddo
            call DMDAVecRestoreArrayF90(meshDA,r,r_array,ierr)
        endif
        call MPI_Barrier(comm,ierr)

    end subroutine linear_rhs

    subroutine snes_rhs_fx_center(snes,x,f,null_int,ierr)

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

    end subroutine snes_rhs_fx_center

    subroutine snes_rhs_fx_4ord(snes,x,f,null_int,ierr)

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
        character(len=20) :: str_norm
        character(len=6) :: str_it
        PetscErrorCode :: ierr
        PetscInt :: it,dummy
        SNES :: snes
        Vec :: f

        call SNESGetFunction(snes,f,PETSC_NULL_FUNCTION,dummy,ierr)
        call VecNorm(f,NORM_INFINITY,nrm,ierr)
        write(str_it,"(I5)") it
        write(str_norm,"(ES20.12)") nrm
        call PetscPrintf(PETSC_COMM_WORLD," "//str_it//" < residual i-Norm > "//str_norm//"\n",ierr)
        if (nrm .le. 1.e-5) reason = SNES_CONVERGED_FNORM_ABS

    end subroutine snes_converged_test

    subroutine push_bc(comm,x)

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
                        x(:,i,j,k)=inlet(:,j,k)
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

    subroutine fx_rhs_Ax_4ord(x,f)

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

    end subroutine fx_rhs_Ax_4ord

    subroutine fx_rhs_Ax_center(x,f)

        implicit none
        PetscScalar,pointer :: fr(:,:,:,:),xr(:,:,:,:)
        integer :: ic_index,jc_index,kc_index
        integer :: lib,lie,ljb,lje,lkb,lke
        PetscScalar,dimension(5,10) :: yi
        PetscErrorCode :: ierr
        Vec,intent(inout) :: f
        Vec,intent(in) :: x
        integer :: li,lj,lk
        integer :: i,j,k
        Vec :: bell

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

                        f(:,i,j,k)=-1.0d0*(matmul(D,yi(:,1))+                                   &
                        &          matmul(A,yi(:,2))+matmul(B,yi(:,3))+matmul(C,yi(:,4))-       &
                        &          matmul(Vxx,yi(:,5))-matmul(Vyy,yi(:,6))-matmul(Vzz,yi(:,7))- &
                        &          matmul(Vxy,yi(:,8))-matmul(Vxz,yi(:,9))-matmul(Vyz,yi(:,10)))

                    endif
                enddo
            enddo
        enddo
        end associate
        call DMDAVecRestoreArrayReadF90(meshDA,bell,xr,ierr)
        call DMDAVecRestoreArrayF90(meshDA,F,fr,ierr)
        call DMRestoreLocalVector(meshDA,bell,ierr)

    end subroutine fx_rhs_Ax_center

    subroutine set_medDA(comm,da)
        implicit none
        integer,intent(in) :: comm
        DM,intent(inout) :: da
        PetscErrorCode :: ierr

        call DMDACreate3d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC,&
        &                 DMDA_STENCIL_BOX, in, jn, kn, nx, ny, nz,&
        &                 5, 1, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, da, ierr)
        call DMSetMatType(da,MATBAIJ,ierr)
        call DMSetUp(da, ierr)
    end subroutine set_medDA

    subroutine set_subDA(da)

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

        implicit none
        PetscErrorCode :: ierr

        call DMGetLocalVector(meshDA,localx,ierr)
        call VecZeroEntries(localx,ierr)
        call DMGetGlobalVector(subDA,subx,ierr)
        call VecZeroEntries(subx,ierr)
    end subroutine init_sub_vecs

    subroutine get_subf(f,subf)

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

end module mod_iterate
