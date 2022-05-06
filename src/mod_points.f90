#include <slepc/finclude/slepc.h>

module mod_points
! -----------------------------------------------------------
!
!   这个模块计算流场中基本流的物理量对坐标的导数
!
!       call partial_derivatives(comm) 计算基本流的物理量对坐标的导数的函数
!
!           1).call allocate_memoru() 分配内存
!
!           2).call partial_derivatives_on_IJK() 计算在计算坐标下的导数的函数
!
!           3).call from_IJK_to_XYZ() 将计算坐标转换到物理坐标下的函数
!
!               call turnItoX(fx,fi,fj,fk,ix,jx,kx) 一阶导数坐标转换函数
!
!           4).call insert_to_BF() 将数据整合到BF结构中
!
!               call insert(obj,arr) 将数组形式的数据插入bf结构的函数
!
!           5).call deallocate_memoru() 释放内存
!
! -----------------------------------------------------------
    use penf, only: R_P
    use mod_parameters
    use petsc
    implicit none
    public :: partial_derivatives
    private
    real(R_P), allocatable, dimension(:, :, :, :) :: qq_x_local_array
    real(R_P), allocatable, dimension(:, :, :, :) :: qq_y_local_array
    real(R_P), allocatable, dimension(:, :, :, :) :: qq_z_local_array
    real(R_P), allocatable, dimension(:, :, :, :) :: qqxx,qqyy,qqzz
    real(R_P), allocatable, dimension(:, :, :, :) :: qqxy,qqxz,qqyz
    real(R_P), allocatable, dimension(:, :, :, :) :: qqxi,qqyi,qqzi
    real(R_P), allocatable, dimension(:, :, :, :) :: qqxj,qqyj,qqzj
    real(R_P), allocatable, dimension(:, :, :, :) :: qqxk,qqyk,qqzk
    real(R_P), allocatable, dimension(:, :, :, :) :: qqi,qqj,qqk
    real(R_P), allocatable, dimension(:, :, :, :) :: qqx,qqy,qqz
    Vec :: QQ_X_local,QQ_Y_local,QQ_Z_local
    PetscErrorCode :: ierr
    Vec :: QQ_X,QQ_Y,QQ_Z
    contains
    subroutine partial_derivatives(comm)
        implicit none
        PetscInt,intent(in) :: comm
        call allocate_memory()
        call partial_derivatives_on_IJK()
        call from_IJK_to_XYZ()
        call insert_to_BF()
        call deallocate_memory()
        call MPI_Barrier(comm,ierr)
    end subroutine partial_derivatives

    subroutine allocate_memory()
        implicit none
        allocate(bf(is:ie,js:je,ks:ke))
        call DMGetGlobalVector(meshDA, QQ_X, ierr)
        call VecDuplicate(QQ_X, QQ_Y, ierr)
        call VecDuplicate(QQ_X, QQ_Z, ierr)
        call DMGetLocalVector(meshDA, QQ_X_local, ierr)
        call VecDuplicate(QQ_X_local, QQ_Y_local, ierr)
        call VecDuplicate(QQ_X_local, QQ_Z_local, ierr)
        call VecZeroEntries(QQ_X,ierr)
        call VecZeroEntries(QQ_Y,ierr)
        call VecZeroEntries(QQ_Z,ierr)
        call VecZeroEntries(QQ_X_local,ierr)
        call VecZeroEntries(QQ_Y_local,ierr)
        call VecZeroEntries(QQ_Z_local,ierr)
        allocate(qqi(5,is:ie,js:je,ks:ke))
        allocate(qqj(5,is:ie,js:je,ks:ke))
        allocate(qqk(5,is:ie,js:je,ks:ke))
        allocate(qqx(5,is:ie,js:je,ks:ke))
        allocate(qqy(5,is:ie,js:je,ks:ke))
        allocate(qqz(5,is:ie,js:je,ks:ke))
        allocate(qqxi(5,is:ie,js:je,ks:ke))
        allocate(qqxj(5,is:ie,js:je,ks:ke))
        allocate(qqxk(5,is:ie,js:je,ks:ke))
        allocate(qqyi(5,is:ie,js:je,ks:ke))
        allocate(qqyj(5,is:ie,js:je,ks:ke))
        allocate(qqyk(5,is:ie,js:je,ks:ke))
        allocate(qqzi(5,is:ie,js:je,ks:ke))
        allocate(qqzj(5,is:ie,js:je,ks:ke))
        allocate(qqzk(5,is:ie,js:je,ks:ke))
        allocate(qqxx(5,is:ie,js:je,ks:ke))
        allocate(qqyy(5,is:ie,js:je,ks:ke))
        allocate(qqzz(5,is:ie,js:je,ks:ke))
        allocate(qqxy(5,is:ie,js:je,ks:ke))
        allocate(qqxz(5,is:ie,js:je,ks:ke))
        allocate(qqyz(5,is:ie,js:je,ks:ke))
        allocate(qq_x_local_array(5,igs:ige,jgs:jge,kgs:kge))
        allocate(qq_y_local_array(5,igs:ige,jgs:jge,kgs:kge))
        allocate(qq_z_local_array(5,igs:ige,jgs:jge,kgs:kge))
    end subroutine allocate_memory

    subroutine partial_derivatives_on_IJK()
        use mod_difference,only:fd1
        implicit none
        select case (lns_mode)
        case(2)
            call fd1(qqi,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,1,5)
            call fd1(qqj,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,2,5)
            qqk=0.0d0
        case(3)
            call fd1(qqi,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,1,5)
            call fd1(qqj,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,2,5)
            call fd1(qqk,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,3,5)
        end select
    end subroutine partial_derivatives_on_IJK

    subroutine from_IJK_to_XYZ()
        use mod_difference,only:fd1
        implicit none
        PetscScalar,pointer :: tmp(:,:,:,:)
        integer :: l,i,j,k
        select case (lns_mode)
        case(2)
            do l=1,5
                do i=is,ie
                    do j=js,je
                        do k=ks,ke
                            call turnItoX(qqx(l,i,j,k),qqi(l,i,j,k),qqj(l,i,j,k),qqk(l,i,j,k), &
                            &             xi_x(i,j,k),eta_x(i,j,k),phi_x(i,j,k))
                            call turnItoX(qqy(l,i,j,k),qqi(l,i,j,k),qqj(l,i,j,k),qqk(l,i,j,k), &
                            &             xi_y(i,j,k),eta_y(i,j,k),phi_y(i,j,k))
                        enddo
                    enddo
                enddo
            enddo
            qqz=0.0d0
        case(3)
            do l=1,5
                do i=is,ie
                    do j=js,je
                        do k=ks,ke
                            call turnItoX(qqx(l,i,j,k),qqi(l,i,j,k),qqj(l,i,j,k),qqk(l,i,j,k), &
                            &             xi_x(i,j,k),eta_x(i,j,k),phi_x(i,j,k))
                            call turnItoX(qqy(l,i,j,k),qqi(l,i,j,k),qqj(l,i,j,k),qqk(l,i,j,k), &
                            &             xi_y(i,j,k),eta_y(i,j,k),phi_y(i,j,k))
                            call turnItoX(qqz(l,i,j,k),qqi(l,i,j,k),qqj(l,i,j,k),qqk(l,i,j,k), &
                            &             xi_z(i,j,k),eta_z(i,j,k),phi_z(i,j,k))
                        enddo
                    enddo
                enddo
            enddo
        end select
        call DMDAVecGetArrayF90(meshDA, QQ_X, tmp, ierr)
        tmp(:,is:ie,js:je,ks:ke) = qqx(:,is:ie,js:je,ks:ke)
        call DMDAVecRestoreArrayF90(meshDA, QQ_X, tmp, ierr)
        call DMGlobalToLocalBegin(meshDA, QQ_X, INSERT_VALUES, QQ_X_local, ierr)
        call DMGlobalToLocalEnd(meshDA, QQ_X, INSERT_VALUES, QQ_X_local, ierr)
        call DMDAVecGetArrayF90(meshDA, QQ_Y, tmp, ierr)
        tmp(:,is:ie,js:je,ks:ke) = qqy(:,is:ie,js:je,ks:ke)
        call DMDAVecRestoreArrayF90(meshDA, QQ_Y, tmp, ierr)
        call DMGlobalToLocalBegin(meshDA, QQ_Y, INSERT_VALUES, QQ_Y_local, ierr)
        call DMGlobalToLocalEnd(meshDA, QQ_Y, INSERT_VALUES, QQ_Y_local, ierr)
        call DMDAVecGetArrayF90(meshDA, QQ_Z, tmp, ierr)
        tmp(:,is:ie,js:je,ks:ke) = qqz(:,is:ie,js:je,ks:ke)
        call DMDAVecRestoreArrayF90(meshDA, QQ_Z, tmp, ierr)
        call DMGlobalToLocalBegin(meshDA, QQ_Z, INSERT_VALUES, QQ_Z_local, ierr)
        call DMGlobalToLocalEnd(meshDA, QQ_Z, INSERT_VALUES, QQ_Z_local, ierr)
        call DMDAVecGetArrayReadF90(meshDA, QQ_X_local, tmp, ierr)
        qq_x_local_array(:,igs:ige,jgs:jge,kgs:kge)=tmp(:,igs:ige,jgs:jge,kgs:kge)
        call DMDAVecRestoreArrayReadF90(meshDA, QQ_X_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(meshDA, QQ_Y_local, tmp, ierr)
        qq_y_local_array(:,igs:ige,jgs:jge,kgs:kge)=tmp(:,igs:ige,jgs:jge,kgs:kge)
        call DMDAVecRestoreArrayReadF90(meshDA, QQ_Y_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(meshDA, QQ_Z_local, tmp, ierr)
        qq_z_local_array(:,igs:ige,jgs:jge,kgs:kge)=tmp(:,igs:ige,jgs:jge,kgs:kge)
        call DMDAVecRestoreArrayReadF90(meshDA, QQ_Z_local, tmp, ierr)
        select case (lns_mode)
        case(2)
            call fd1(qqxi,is,ie,js,je,ks,ke,qq_x_local_array,igs,ige,jgs,jge,kgs,kge,1,5)
            call fd1(qqxj,is,ie,js,je,ks,ke,qq_x_local_array,igs,ige,jgs,jge,kgs,kge,2,5)
            call fd1(qqyi,is,ie,js,je,ks,ke,qq_y_local_array,igs,ige,jgs,jge,kgs,kge,1,5)
            call fd1(qqyj,is,ie,js,je,ks,ke,qq_y_local_array,igs,ige,jgs,jge,kgs,kge,2,5)
            call fd1(qqzi,is,ie,js,je,ks,ke,qq_z_local_array,igs,ige,jgs,jge,kgs,kge,1,5)
            call fd1(qqzj,is,ie,js,je,ks,ke,qq_z_local_array,igs,ige,jgs,jge,kgs,kge,2,5)
            qqxk=0.0d0;qqyk=0.0d0;qqzk=0.0d0
            do l=1,5
                do i=is,ie
                    do j=js,je
                        do k=ks,ke
                            call turnItoX(qqxx(l,i,j,k),qqxi(l,i,j,k),qqxj(l,i,j,k),qqxk(l,i,j,k), &
                            &             xi_x(i,j,k),eta_x(i,j,k),phi_x(i,j,k))
                            call turnItoX(qqyy(l,i,j,k),qqyi(l,i,j,k),qqyj(l,i,j,k),qqyk(l,i,j,k), &
                            &             xi_y(i,j,k),eta_y(i,j,k),phi_y(i,j,k))
                            call turnItoX(qqxy(l,i,j,k),qqxi(l,i,j,k),qqxj(l,i,j,k),qqxk(l,i,j,k), &
                            &             xi_y(i,j,k),eta_y(i,j,k),phi_y(i,j,k))
                        enddo
                    enddo
                enddo
            enddo
            qqzz=0.0d0;qqxz=0.0d0;qqyz=0.0d0
        case(3)
            call fd1(qqxi,is,ie,js,je,ks,ke,qq_x_local_array,igs,ige,jgs,jge,kgs,kge,1,5)
            call fd1(qqxj,is,ie,js,je,ks,ke,qq_x_local_array,igs,ige,jgs,jge,kgs,kge,2,5)
            call fd1(qqxk,is,ie,js,je,ks,ke,qq_x_local_array,igs,ige,jgs,jge,kgs,kge,3,5)
            call fd1(qqyi,is,ie,js,je,ks,ke,qq_y_local_array,igs,ige,jgs,jge,kgs,kge,1,5)
            call fd1(qqyj,is,ie,js,je,ks,ke,qq_y_local_array,igs,ige,jgs,jge,kgs,kge,2,5)
            call fd1(qqyk,is,ie,js,je,ks,ke,qq_y_local_array,igs,ige,jgs,jge,kgs,kge,3,5)
            call fd1(qqzi,is,ie,js,je,ks,ke,qq_z_local_array,igs,ige,jgs,jge,kgs,kge,1,5)
            call fd1(qqzj,is,ie,js,je,ks,ke,qq_z_local_array,igs,ige,jgs,jge,kgs,kge,2,5)
            call fd1(qqzk,is,ie,js,je,ks,ke,qq_z_local_array,igs,ige,jgs,jge,kgs,kge,3,5)
            do l=1,5
                do i=is,ie
                    do j=js,je
                        do k=ks,ke
                            call turnItoX(qqxx(l,i,j,k),qqxi(l,i,j,k),qqxj(l,i,j,k),qqxk(l,i,j,k), &
                            &             xi_x(i,j,k),eta_x(i,j,k),phi_x(i,j,k))
                            call turnItoX(qqyy(l,i,j,k),qqyi(l,i,j,k),qqyj(l,i,j,k),qqyk(l,i,j,k), &
                            &             xi_y(i,j,k),eta_y(i,j,k),phi_y(i,j,k))
                            call turnItoX(qqzz(l,i,j,k),qqzi(l,i,j,k),qqzj(l,i,j,k),qqzk(l,i,j,k), &
                            &             xi_z(i,j,k),eta_z(i,j,k),phi_z(i,j,k))
                            call turnItoX(qqxy(l,i,j,k),qqxi(l,i,j,k),qqxj(l,i,j,k),qqxk(l,i,j,k), &
                            &             xi_y(i,j,k),eta_y(i,j,k),phi_y(i,j,k))
                            call turnItoX(qqxz(l,i,j,k),qqxi(l,i,j,k),qqxj(l,i,j,k),qqxk(l,i,j,k), &
                            &             xi_z(i,j,k),eta_z(i,j,k),phi_z(i,j,k))
                            call turnItoX(qqyz(l,i,j,k),qqyi(l,i,j,k),qqyj(l,i,j,k),qqyk(l,i,j,k), &
                            &             xi_z(i,j,k),eta_z(i,j,k),phi_z(i,j,k))
                        enddo
                    enddo
                enddo
            enddo
        end select
    end subroutine from_IJK_to_XYZ

    subroutine turnItoX(fx,fi,fj,fk,ix,jx,kx)
        implicit none
        real(R_P),intent(out) :: fx
        real(R_P),intent(in) :: fi,fj,fk
        real(R_P),intent(in) :: ix,jx,kx
        fx = ix*fi+jx*fj+kx*fk
    end subroutine turnItoX

    subroutine insert_to_BF()
        implicit none
        integer :: i,j,k
        do k=ks,ke
            do j=js,je
                do i=is,ie
                    call insert(bf(i,j,k)%BF,qq(:,i,j,k))
                    call insert(bf(i,j,k)%BFDx,qqx(:,i,j,k))
                    call insert(bf(i,j,k)%BFDy,qqy(:,i,j,k))
                    call insert(bf(i,j,k)%BFDz,qqz(:,i,j,k))
                    call insert(bf(i,j,k)%BFDxx,qqxx(:,i,j,k))
                    call insert(bf(i,j,k)%BFDyy,qqyy(:,i,j,k))
                    call insert(bf(i,j,k)%BFDzz,qqzz(:,i,j,k))
                    call insert(bf(i,j,k)%BFDxy,qqxy(:,i,j,k))
                    call insert(bf(i,j,k)%BFDxz,qqxz(:,i,j,k))
                    call insert(bf(i,j,k)%BFDyz,qqyz(:,i,j,k))
                enddo
            enddo
        enddo
    end subroutine insert_to_BF

    subroutine insert(obj,arr)
        use mod_flowtype
        implicit none
        real(R_P),dimension(5),intent(in) :: arr
        type(basetype),intent(out) :: obj
        obj%rho = arr(1)
        obj%x = arr(2)
        obj%y = arr(3)
        obj%z = arr(4)
        obj%T = arr(5)
    end subroutine insert

    subroutine deallocate_memory()
        implicit none
        deallocate(qq)
        deallocate(qqi);deallocate(qqj);deallocate(qqk)
        deallocate(qqx);deallocate(qqy);deallocate(qqz)
        deallocate(qqxx);deallocate(qqyy);deallocate(qqzz)
        deallocate(qqxy);deallocate(qqxz);deallocate(qqyz)
        deallocate(qqxi);deallocate(qqxj);deallocate(qqxk)
        deallocate(qqyi);deallocate(qqyj);deallocate(qqyk)
        deallocate(qqzi);deallocate(qqzj);deallocate(qqzk)
        deallocate(qq_x_local_array)
        deallocate(qq_y_local_array)
        deallocate(qq_z_local_array)
    end subroutine deallocate_memory

end module mod_points
