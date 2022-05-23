!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_metrics
! ----------------------------------------------------------------
!
!  这个模块计算度量系数。
!
!       call metric_coefficient(comm) 获得度量系数
!
!           1).call allocate_memory() 分配内存
!
!           2).call compute_contravariant_metrics() 计算坐标的一阶导数的函数
!
!           3).call compute_convariant_metrics() 计算度量系数的函数
!
!           4).call deallocate_memory() 释放内存
!
! ----------------------------------------------------------------
    use penf, only: R_P
    use mod_parameters
    use petsc
    implicit none
    public :: metric_coefficient
    private
    Vec :: XIX_local, XIY_local, XIZ_local, ETAX_local, ETAY_local, ETAZ_local, PHIX_local, PHIY_local, PHIZ_local
    real(R_P), dimension(:, :, :), allocatable :: x_xi,x_eta,x_phi,y_xi,y_eta,y_phi,z_xi,z_eta,z_phi
    Vec :: XIX, XIY, XIZ, ETAX, ETAY, ETAZ, PHIX, PHIY, PHIZ
    real(R_P), dimension(:, :, :), allocatable :: jacobi
    PetscErrorCode  :: ierr
    contains
    subroutine metric_coefficient(comm)
        implicit none
        PetscInt,intent(in) :: comm
        call allocate_memory()
        call compute_contravariant_metrics()
        call compute_convariant_metrics()
        call delivery_by_dmda()
        call deallocate_memory()
        call MPI_Barrier(comm,ierr)
    end subroutine metric_coefficient

    subroutine allocate_memory()
        implicit none
        call DMGetGlobalVector(DA, XIX, ierr)
        call VecDuplicate(XIX, XIY, ierr)
        call VecDuplicate(XIX, XIZ, ierr)
        call VecDuplicate(XIX, ETAX, ierr)
        call VecDuplicate(XIX, ETAY, ierr)
        call VecDuplicate(XIX, ETAZ, ierr)
        call VecDuplicate(XIX, PHIX, ierr)
        call VecDuplicate(XIX, PHIY, ierr)
        call VecDuplicate(XIX, PHIZ, ierr)

        call DMGetLocalVector(DA, XIX_local, ierr)
        call VecDuplicate(XIX_local, XIY_local, ierr)
        call VecDuplicate(XIX_local, XIZ_local, ierr)
        call VecDuplicate(XIX_local, ETAX_local, ierr)
        call VecDuplicate(XIX_local, ETAY_local, ierr)
        call VecDuplicate(XIX_local, ETAZ_local, ierr)
        call VecDuplicate(XIX_local, PHIX_local, ierr)
        call VecDuplicate(XIX_local, PHIY_local, ierr)
        call VecDuplicate(XIX_local, PHIZ_local, ierr)

        allocate(x_xi(is:ie, js:je, ks:ke))
        allocate(y_xi(is:ie, js:je, ks:ke))
        allocate(z_xi(is:ie, js:je, ks:ke))
        allocate(x_eta(is:ie, js:je, ks:ke))
        allocate(y_eta(is:ie, js:je, ks:ke))
        allocate(z_eta(is:ie, js:je, ks:ke))
        allocate(x_phi(is:ie, js:je, ks:ke))
        allocate(y_phi(is:ie, js:je, ks:ke))
        allocate(z_phi(is:ie, js:je, ks:ke))

        allocate(xi_x(is:ie, js:je, ks:ke))
        allocate(xi_y(is:ie, js:je, ks:ke))
        allocate(xi_z(is:ie, js:je, ks:ke))
        allocate(eta_x(is:ie, js:je, ks:ke))
        allocate(eta_y(is:ie, js:je, ks:ke))
        allocate(eta_z(is:ie, js:je, ks:ke))
        allocate(phi_x(is:ie, js:je, ks:ke))
        allocate(phi_y(is:ie, js:je, ks:ke))
        allocate(phi_z(is:ie, js:je, ks:ke))

        allocate(xi_xx(is:ie, js:je, ks:ke))
        allocate(xi_yy(is:ie, js:je, ks:ke))
        allocate(xi_zz(is:ie, js:je, ks:ke))
        allocate(xi_xy(is:ie, js:je, ks:ke))
        allocate(xi_xz(is:ie, js:je, ks:ke))
        allocate(xi_yz(is:ie, js:je, ks:ke))
        allocate(eta_xx(is:ie, js:je, ks:ke))
        allocate(eta_yy(is:ie, js:je, ks:ke))
        allocate(eta_zz(is:ie, js:je, ks:ke))
        allocate(eta_xy(is:ie, js:je, ks:ke))
        allocate(eta_xz(is:ie, js:je, ks:ke))
        allocate(eta_yz(is:ie, js:je, ks:ke))
        allocate(phi_xx(is:ie, js:je, ks:ke))
        allocate(phi_yy(is:ie, js:je, ks:ke))
        allocate(phi_zz(is:ie, js:je, ks:ke))
        allocate(phi_xy(is:ie, js:je, ks:ke))
        allocate(phi_xz(is:ie, js:je, ks:ke))
        allocate(phi_yz(is:ie, js:je, ks:ke))

        allocate(jacobi(is:ie, js:je, ks:ke))
    end subroutine allocate_memory

    subroutine compute_contravariant_metrics()
        use mod_difference
        implicit none
        select case (lns_mode)
        case(2)
            call fd1(x_xi,is,ie,js,je,ks,ke,xx,igs,ige,jgs,jge,kgs,kge,1,1)
            call fd1(y_xi,is,ie,js,je,ks,ke,yy,igs,ige,jgs,jge,kgs,kge,1,1)
            call fd1(x_eta,is,ie,js,je,ks,ke,xx,igs,ige,jgs,jge,kgs,kge,2,1)
            call fd1(y_eta,is,ie,js,je,ks,ke,yy,igs,ige,jgs,jge,kgs,kge,2,1)
            x_phi=0.0d0;y_phi=0.0d0;z_xi=0.0d0;z_eta=0.0d0;z_phi=0.0d0
        case(3)
            call fd1(x_xi,is,ie,js,je,ks,ke,xx,igs,ige,jgs,jge,kgs,kge,1,1)
            call fd1(y_xi,is,ie,js,je,ks,ke,yy,igs,ige,jgs,jge,kgs,kge,1,1)
            call fd1(z_xi,is,ie,js,je,ks,ke,zz,igs,ige,jgs,jge,kgs,kge,1,1)
            call fd1(x_eta,is,ie,js,je,ks,ke,xx,igs,ige,jgs,jge,kgs,kge,2,1)
            call fd1(y_eta,is,ie,js,je,ks,ke,yy,igs,ige,jgs,jge,kgs,kge,2,1)
            call fd1(z_eta,is,ie,js,je,ks,ke,zz,igs,ige,jgs,jge,kgs,kge,2,1)
            call fd1(x_phi,is,ie,js,je,ks,ke,xx,igs,ige,jgs,jge,kgs,kge,3,1)
            call fd1(y_phi,is,ie,js,je,ks,ke,yy,igs,ige,jgs,jge,kgs,kge,3,1)
            call fd1(z_phi,is,ie,js,je,ks,ke,zz,igs,ige,jgs,jge,kgs,kge,3,1)
        end select
    end subroutine compute_contravariant_metrics

    subroutine compute_convariant_metrics()
        select case (lns_mode)
        case(2)
            call compute_convariant_metrics_2d()
        case(3)
            call compute_convariant_metrics_3d()
        end select
    end subroutine compute_convariant_metrics

    subroutine compute_convariant_metrics_2d()
        use mod_difference
        implicit none
        PetscScalar, pointer :: tmp(:, :, :)
        integer :: i, j, k
        real(R_P), allocatable, dimension(:,:,:) :: xi_x_local
        real(R_P), allocatable, dimension(:,:,:) :: xi_y_local
        real(R_P), allocatable, dimension(:,:,:) :: eta_x_local
        real(R_P), allocatable, dimension(:,:,:) :: eta_y_local
        real(R_P), allocatable, dimension(:,:,:) :: xi_x_xi
        real(R_P), allocatable, dimension(:,:,:) :: xi_x_eta
        real(R_P), allocatable, dimension(:,:,:) :: xi_y_xi
        real(R_P), allocatable, dimension(:,:,:) :: xi_y_eta
        real(R_P), allocatable, dimension(:,:,:) :: eta_x_xi
        real(R_P), allocatable, dimension(:,:,:) :: eta_x_eta
        real(R_P), allocatable, dimension(:,:,:) :: eta_y_xi
        real(R_P), allocatable, dimension(:,:,:) :: eta_y_eta
        do k=ks,ke
            do j=js,je
                do i=is,ie
                    jacobi(i,j,k)=1.0d0/(x_xi(i,j,k)*y_eta(i,j,k)-y_xi(i,j,k)*x_eta(i,j,k))
                enddo
            enddo
        enddo

        do k=ks,ke
            do j=js,je
                do i=is,ie
                    xi_x(i,j,k) =y_eta(i,j,k)*jacobi(i,j,k)
                    xi_y(i,j,k) =-1.0d0*x_eta(i,j,k)*jacobi(i,j,k)
                    eta_x(i,j,k)=-1.0d0*y_xi(i,j,k)*jacobi(i,j,k)
                    eta_y(i,j,k)=x_xi(i,j,k)*jacobi(i,j,k)
                enddo
            enddo
        enddo
        call DMDAVecGetArrayF90(DA, XIX, tmp, ierr)
        tmp(:,:,:)=xi_x(:,:,:)
        call DMDAVecRestoreArrayF90(DA, XIX, tmp, ierr)
        call DMGlobalToLocalBegin(DA, XIX, INSERT_VALUES, XIX_local, ierr)

        call DMDAVecGetArrayF90(DA, XIY, tmp, ierr)
        tmp(:,:,:)=xi_y(:,:,:)
        call DMDAVecRestoreArrayF90(DA, XIY, tmp, ierr)
        call DMGlobalToLocalEnd(DA, XIX, INSERT_VALUES, XIX_local, ierr)
        call DMGlobalToLocalBegin(DA, XIY, INSERT_VALUES, XIY_local, ierr)

        call DMDAVecGetArrayF90(DA, ETAX, tmp, ierr)
        tmp(:,:,:)=eta_x(:,:,:)
        call DMDAVecRestoreArrayF90(DA, ETAX, tmp, ierr)
        call DMGlobalToLocalEnd(DA, XIY, INSERT_VALUES, XIY_local, ierr)
        call DMGlobalToLocalBegin(DA, ETAX, INSERT_VALUES, ETAX_local, ierr)

        call DMDAVecGetArrayF90(DA, ETAY, tmp, ierr)
        tmp(:,:,:)=eta_y(:,:,:)
        call DMDAVecRestoreArrayF90(DA, ETAY, tmp, ierr)
        call DMGlobalToLocalEnd(DA, ETAX, INSERT_VALUES, ETAX_local, ierr)
        call DMGlobalToLocalBegin(DA, ETAY, INSERT_VALUES, ETAY_local, ierr)
        call DMGlobalToLocalEND(DA, ETAY, INSERT_VALUES, ETAY_local, ierr)

        allocate(xi_x_local(igs:ige,jgs:jge,kgs:kge))
        allocate(xi_y_local(igs:ige,jgs:jge,kgs:kge))
        allocate(eta_x_local(igs:ige,jgs:jge,kgs:kge))
        allocate(eta_y_local(igs:ige,jgs:jge,kgs:kge))

        call DMDAVecGetArrayReadF90(DA, XIX_local, tmp, ierr)
        xi_x_local=tmp
        call DMDAVecRestoreArrayReadF90(DA, XIX_local, tmp, ierr)

        call DMDAVecGetArrayReadF90(DA, XIY_local, tmp, ierr)
        xi_y_local=tmp
        call DMDAVecRestoreArrayReadF90(DA, XIY_local, tmp, ierr)

        call DMDAVecGetArrayReadF90(DA, ETAX_local, tmp, ierr)
        eta_x_local=tmp
        call DMDAVecRestoreArrayReadF90(DA, ETAX_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, ETAY_local, tmp, ierr)
        eta_y_local=tmp
        call DMDAVecRestoreArrayReadF90(DA, ETAY_local, tmp, ierr)

        allocate(xi_x_xi(is:ie,js:je,ks:ke))
        allocate(xi_x_eta(is:ie,js:je,ks:ke))
        allocate(xi_y_xi(is:ie,js:je,ks:ke))
        allocate(xi_y_eta(is:ie,js:je,ks:ke))
        allocate(eta_x_xi(is:ie,js:je,ks:ke))
        allocate(eta_x_eta(is:ie,js:je,ks:ke))
        allocate(eta_y_xi(is:ie,js:je,ks:ke))
        allocate(eta_y_eta(is:ie,js:je,ks:ke))

        call fd1(xi_x_xi,is,ie,js,je,ks,ke,xi_x_local,igs,ige,jgs,jge,kgs,kge,1,1)
        call fd1(xi_y_xi,is,ie,js,je,ks,ke,xi_y_local,igs,ige,jgs,jge,kgs,kge,1,1)
        call fd1(eta_x_xi,is,ie,js,je,ks,ke,eta_x_local,igs,ige,jgs,jge,kgs,kge,1,1)
        call fd1(eta_y_xi,is,ie,js,je,ks,ke,eta_y_local,igs,ige,jgs,jge,kgs,kge,1,1)

        call fd1(xi_x_eta,is,ie,js,je,ks,ke,xi_x_local,igs,ige,jgs,jge,kgs,kge,2,1)
        call fd1(xi_y_eta,is,ie,js,je,ks,ke,xi_y_local,igs,ige,jgs,jge,kgs,kge,2,1)
        call fd1(eta_x_eta,is,ie,js,je,ks,ke,eta_x_local,igs,ige,jgs,jge,kgs,kge,2,1)
        call fd1(eta_y_eta,is,ie,js,je,ks,ke,eta_y_local,igs,ige,jgs,jge,kgs,kge,2,1)

        do k=ks,ke
            do j=js,je
                do i=is,ie
                    xi_xx(i,j,k)=xi_x(i,j,k)*xi_x_xi(i,j,k)+eta_x(i,j,k)*xi_x_eta(i,j,k)
                    xi_yy(i,j,k)=xi_y(i,j,k)*xi_y_xi(i,j,k)+eta_y(i,j,k)*xi_y_eta(i,j,k)
                    xi_xy(i,j,k)=xi_y(i,j,k)*xi_x_xi(i,j,k)+eta_y(i,j,k)*xi_x_eta(i,j,k)
                    eta_xx(i,j,k)=xi_x(i,j,k)*eta_x_xi(i,j,k)+eta_x(i,j,k)*eta_x_eta(i,j,k)
                    eta_yy(i,j,k)=xi_y(i,j,k)*eta_y_xi(i,j,k)+eta_y(i,j,k)*eta_y_eta(i,j,k)
                    eta_xy(i,j,k)=xi_y(i,j,k)*eta_x_xi(i,j,k)+eta_y(i,j,k)*eta_x_eta(i,j,k)
                enddo
            enddo
        enddo

        xi_zz=0.0d0;xi_xz=0.0d0;xi_yz=0.0d0;eta_zz=0.0d0;eta_xz=0.0d0;eta_yz=0.0d0
        phi_xx=0.0d0;phi_yy=0.0d0;phi_zz=0.0d0;phi_xy=0.0d0;phi_xz=0.0d0;phi_yz=0.0d0

        deallocate(xi_x_xi)
        deallocate(xi_x_eta)
        deallocate(xi_y_xi)
        deallocate(xi_y_eta)
        deallocate(eta_x_xi)
        deallocate(eta_x_eta)
        deallocate(eta_y_xi)
        deallocate(eta_y_eta)
        deallocate(xi_x_local)
        deallocate(xi_y_local)
        deallocate(eta_x_local)
        deallocate(eta_y_local)
    end subroutine compute_convariant_metrics_2d

    subroutine compute_convariant_metrics_3d()
        use mod_difference
        implicit none
        PetscScalar, pointer :: tmp(:, :, :)
        integer :: i, j, k
        real(R_P), allocatable, dimension(:,:,:) :: xi_x_local
        real(R_P), allocatable, dimension(:,:,:) :: xi_y_local
        real(R_P), allocatable, dimension(:,:,:) :: xi_z_local
        real(R_P), allocatable, dimension(:,:,:) :: eta_x_local
        real(R_P), allocatable, dimension(:,:,:) :: eta_y_local
        real(R_P), allocatable, dimension(:,:,:) :: eta_z_local
        real(R_P), allocatable, dimension(:,:,:) :: phi_x_local
        real(R_P), allocatable, dimension(:,:,:) :: phi_y_local
        real(R_P), allocatable, dimension(:,:,:) :: phi_z_local
        real(R_P), allocatable, dimension(:,:,:) :: xi_x_xi
        real(R_P), allocatable, dimension(:,:,:) :: xi_x_eta
        real(R_P), allocatable, dimension(:,:,:) :: xi_x_phi
        real(R_P), allocatable, dimension(:,:,:) :: xi_y_xi
        real(R_P), allocatable, dimension(:,:,:) :: xi_y_eta
        real(R_P), allocatable, dimension(:,:,:) :: xi_y_phi
        real(R_P), allocatable, dimension(:,:,:) :: xi_z_xi
        real(R_P), allocatable, dimension(:,:,:) :: xi_z_eta
        real(R_P), allocatable, dimension(:,:,:) :: xi_z_phi
        real(R_P), allocatable, dimension(:,:,:) :: eta_x_xi
        real(R_P), allocatable, dimension(:,:,:) :: eta_x_eta
        real(R_P), allocatable, dimension(:,:,:) :: eta_x_phi
        real(R_P), allocatable, dimension(:,:,:) :: eta_y_xi
        real(R_P), allocatable, dimension(:,:,:) :: eta_y_eta
        real(R_P), allocatable, dimension(:,:,:) :: eta_y_phi
        real(R_P), allocatable, dimension(:,:,:) :: eta_z_xi
        real(R_P), allocatable, dimension(:,:,:) :: eta_z_eta
        real(R_P), allocatable, dimension(:,:,:) :: eta_z_phi
        real(R_P), allocatable, dimension(:,:,:) :: phi_x_xi
        real(R_P), allocatable, dimension(:,:,:) :: phi_x_eta
        real(R_P), allocatable, dimension(:,:,:) :: phi_x_phi
        real(R_P), allocatable, dimension(:,:,:) :: phi_y_xi
        real(R_P), allocatable, dimension(:,:,:) :: phi_y_eta
        real(R_P), allocatable, dimension(:,:,:) :: phi_y_phi
        real(R_P), allocatable, dimension(:,:,:) :: phi_z_xi
        real(R_P), allocatable, dimension(:,:,:) :: phi_z_eta
        real(R_P), allocatable, dimension(:,:,:) :: phi_z_phi
        do k=ks,ke
            do j=js,je
                do i=is,ie
                    jacobi(i, j, k) = 1.0d0/(x_xi(i,j,k)*(y_eta(i,j,k)*z_phi(i,j,k)-z_eta(i,j,k)*y_phi(i,j,k)) &
                                            -x_eta(i,j,k)*(y_xi(i,j,k)*z_phi(i,j,k)-z_xi(i,j,k)*y_phi(i,j,k)) &
                                            +x_phi(i,j,k)*(y_xi(i,j,k)*z_eta(i,j,k)-z_xi(i,j,k)*y_eta(i,j,k)))
                enddo
            enddo
        enddo
        do k=ks,ke
            do j=js,je
                do i=is,ie
                    xi_x(i,j,k)=jacobi(i,j,k)*(y_eta(i,j,k)*z_phi(i,j,k)-z_eta(i,j,k)*y_phi(i,j,k))
                    eta_x(i,j,k)=jacobi(i,j,k)*(y_phi(i,j,k)*z_xi(i,j,k)-z_phi(i,j,k)*y_xi(i,j,k))
                    phi_x(i,j,k)=jacobi(i,j,k)*(y_xi(i,j,k)*z_eta(i,j,k)-z_xi(i,j,k)*y_eta(i,j,k))
                    xi_y(i,j,k)=jacobi(i,j,k)*(z_eta(i,j,k)*x_phi(i,j,k)-x_eta(i,j,k)*z_phi(i,j,k))
                    eta_y(i,j,k)=jacobi(i,j,k)*(z_phi(i,j,k)*x_xi(i,j,k)-x_phi(i,j,k)*z_xi(i,j,k))
                    phi_y(i,j,k)=jacobi(i,j,k)*(z_xi(i,j,k)*x_eta(i,j,k)-x_xi(i,j,k)*z_eta(i,j,k))
                    xi_z(i,j,k)=jacobi(i,j,k)*(x_eta(i,j,k)*y_phi(i,j,k)-y_eta(i,j,k)*x_phi(i,j,k))
                    eta_z(i,j,k)=jacobi(i,j,k)*(x_phi(i,j,k)*y_xi(i,j,k)-y_phi(i,j,k)*x_xi(i,j,k))
                    phi_z(i,j,k)=jacobi(i,j,k)*(x_xi(i,j,k)*y_eta(i,j,k)-y_xi(i,j,k)*x_eta(i,j,k))
                enddo
            enddo
        enddo
        call DMDAVecGetArrayF90(DA, XIX, tmp, ierr)
        tmp(:,:,:) = xi_x(:,:,:)
        call DMDAVecRestoreArrayF90(DA, XIX, tmp, ierr)
        call DMGlobalToLocalBegin(DA, XIX, INSERT_VALUES, XIX_local, ierr)

        call DMDAVecGetArrayF90(DA, XIY, tmp, ierr)
        tmp(:,:,:) = xi_y(:,:,:)
        call DMDAVecRestoreArrayF90(DA, XIY, tmp, ierr)
        call DMGlobalToLocalEnd(DA, XIX, INSERT_VALUES, XIX_local, ierr)
        call DMGlobalToLocalBegin(DA, XIY, INSERT_VALUES, XIY_local, ierr)

        call DMDAVecGetArrayF90(DA, XIZ, tmp, ierr)
        tmp(:,:,:) = xi_z(:,:,:)
        call DMDAVecRestoreArrayF90(DA, XIZ, tmp, ierr)
        call DMGlobalToLocalEnd(DA, XIY, INSERT_VALUES, XIY_local, ierr)
        call DMGlobalToLocalBegin(DA, XIZ, INSERT_VALUES, XIZ_local, ierr)

        call DMDAVecGetArrayF90(DA, ETAX, tmp, ierr)
        tmp(:,:,:) = eta_x(:,:,:)
        call DMDAVecRestoreArrayF90(DA, ETAX, tmp, ierr)
        call DMGlobalToLocalEnd(DA, XIZ, INSERT_VALUES, XIZ_local, ierr)
        call DMGlobalToLocalBegin(DA, ETAX, INSERT_VALUES, ETAX_local, ierr)

        call DMDAVecGetArrayF90(DA, ETAY, tmp, ierr)
        tmp(:,:,:) = eta_y(:,:,:)
        call DMDAVecRestoreArrayF90(DA, ETAY, tmp, ierr)
        call DMGlobalToLocalEnd(DA, ETAX, INSERT_VALUES, ETAX_local, ierr)
        call DMGlobalToLocalBegin(DA, ETAY, INSERT_VALUES, ETAY_local, ierr)

        call DMDAVecGetArrayF90(DA, ETAZ, tmp, ierr)
        tmp(:,:,:) = eta_z(:,:,:)
        call DMDAVecRestoreArrayF90(DA, ETAZ, tmp, ierr)
        call DMGlobalToLocalEnd(DA, ETAY, INSERT_VALUES, ETAY_local, ierr)
        call DMGlobalToLocalBegin(DA, ETAZ, INSERT_VALUES, ETAZ_local, ierr)

        call DMDAVecGetArrayF90(DA, PHIX, tmp, ierr)
        tmp(:,:,:) = phi_x(:,:,:)
        call DMDAVecRestoreArrayF90(DA, PHIX, tmp, ierr)
        call DMGlobalToLocalEnd(DA, ETAZ, INSERT_VALUES, ETAZ_local, ierr)
        call DMGlobalToLocalBegin(DA, PHIX, INSERT_VALUES, PHIX_local, ierr)

        call DMDAVecGetArrayF90(DA, PHIY, tmp, ierr)
        tmp(:,:,:) = phi_y(:,:,:)
        call DMDAVecRestoreArrayF90(DA, PHIY, tmp, ierr)
        call DMGlobalToLocalEnd(DA, PHIX, INSERT_VALUES, PHIX_local, ierr)
        call DMGlobalToLocalBegin(DA, PHIY, INSERT_VALUES, PHIY_local, ierr)

        call DMDAVecGetArrayF90(DA, PHIZ, tmp, ierr)
        tmp(:,:,:) = phi_z(:,:,:)
        call DMDAVecRestoreArrayF90(DA, PHIZ, tmp, ierr)
        call DMGlobalToLocalEnd(DA, PHIY, INSERT_VALUES, PHIY_local, ierr)
        call DMGlobalToLocalBegin(DA, PHIZ, INSERT_VALUES, PHIZ_local, ierr)
        call DMGlobalToLocalEnd(DA, PHIZ, INSERT_VALUES, PHIZ_local, ierr)

        allocate(xi_x_local (igs:ige, jgs:jge, kgs:kge))
        allocate(xi_y_local (igs:ige, jgs:jge, kgs:kge))
        allocate(xi_z_local (igs:ige, jgs:jge, kgs:kge))
        allocate(eta_x_local(igs:ige, jgs:jge, kgs:kge))
        allocate(eta_y_local(igs:ige, jgs:jge, kgs:kge))
        allocate(eta_z_local(igs:ige, jgs:jge, kgs:kge))
        allocate(phi_x_local(igs:ige, jgs:jge, kgs:kge))
        allocate(phi_y_local(igs:ige, jgs:jge, kgs:kge))
        allocate(phi_z_local(igs:ige, jgs:jge, kgs:kge))

        allocate(xi_x_xi(is:ie, js:je, ks:ke))
        allocate(xi_x_eta(is:ie, js:je, ks:ke))
        allocate(xi_x_phi(is:ie, js:je, ks:ke))
        allocate(xi_y_xi(is:ie, js:je, ks:ke))
        allocate(xi_y_eta(is:ie, js:je, ks:ke))
        allocate(xi_y_phi(is:ie, js:je, ks:ke))
        allocate(xi_z_xi(is:ie, js:je, ks:ke))
        allocate(xi_z_eta(is:ie, js:je, ks:ke))
        allocate(xi_z_phi(is:ie, js:je, ks:ke))
        allocate(eta_x_xi(is:ie, js:je, ks:ke))
        allocate(eta_x_eta(is:ie, js:je, ks:ke))
        allocate(eta_x_phi(is:ie, js:je, ks:ke))
        allocate(eta_y_xi(is:ie, js:je, ks:ke))
        allocate(eta_y_eta(is:ie, js:je, ks:ke))
        allocate(eta_y_phi(is:ie, js:je, ks:ke))
        allocate(eta_z_xi(is:ie, js:je, ks:ke))
        allocate(eta_z_eta(is:ie, js:je, ks:ke))
        allocate(eta_z_phi(is:ie, js:je, ks:ke))
        allocate(phi_x_xi(is:ie, js:je, ks:ke))
        allocate(phi_x_eta(is:ie, js:je, ks:ke))
        allocate(phi_x_phi(is:ie, js:je, ks:ke))
        allocate(phi_y_xi(is:ie, js:je, ks:ke))
        allocate(phi_y_eta(is:ie, js:je, ks:ke))
        allocate(phi_y_phi(is:ie, js:je, ks:ke))
        allocate(phi_z_xi(is:ie, js:je, ks:ke))
        allocate(phi_z_eta(is:ie, js:je, ks:ke))
        allocate(phi_z_phi(is:ie, js:je, ks:ke))

        call DMDAVecGetArrayReadF90(DA, XIX_local, tmp, ierr)
        xi_x_local=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, XIX_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, XIY_local, tmp, ierr)
        xi_y_local=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, XIY_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, XIZ_local, tmp, ierr)
        xi_z_local=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, XIZ_local, tmp, ierr)

        call DMDAVecGetArrayReadF90(DA, ETAX_local, tmp, ierr)
        eta_x_local=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, ETAX_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, ETAY_local, tmp, ierr)
        eta_y_local=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, ETAY_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, ETAZ_local, tmp, ierr)
        eta_z_local=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, ETAZ_local, tmp, ierr)

        call DMDAVecGetArrayReadF90(DA, PHIX_local, tmp, ierr)
        phi_x_local=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, PHIX_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, PHIY_local, tmp, ierr)
        phi_y_local=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, PHIY_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, PHIZ_local, tmp, ierr)
        phi_z_local=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, PHIZ_local, tmp, ierr)

        call fd1(xi_x_xi,is,ie,js,je,ks,ke,xi_x_local,igs,ige,jgs,jge,kgs,kge,1,1)
        call fd1(xi_y_xi,is,ie,js,je,ks,ke,xi_y_local,igs,ige,jgs,jge,kgs,kge,1,1)
        call fd1(xi_z_xi,is,ie,js,je,ks,ke,xi_z_local,igs,ige,jgs,jge,kgs,kge,1,1)
        call fd1(eta_x_xi,is,ie,js,je,ks,ke,eta_x_local,igs,ige,jgs,jge,kgs,kge,1,1)
        call fd1(eta_y_xi,is,ie,js,je,ks,ke,eta_y_local,igs,ige,jgs,jge,kgs,kge,1,1)
        call fd1(eta_z_xi,is,ie,js,je,ks,ke,eta_z_local,igs,ige,jgs,jge,kgs,kge,1,1)
        call fd1(phi_x_xi,is,ie,js,je,ks,ke,phi_x_local,igs,ige,jgs,jge,kgs,kge,1,1)
        call fd1(phi_y_xi,is,ie,js,je,ks,ke,phi_y_local,igs,ige,jgs,jge,kgs,kge,1,1)
        call fd1(phi_z_xi,is,ie,js,je,ks,ke,phi_z_local,igs,ige,jgs,jge,kgs,kge,1,1)

        call fd1(xi_x_eta,is,ie,js,je,ks,ke,xi_x_local,igs,ige,jgs,jge,kgs,kge,2,1)
        call fd1(xi_y_eta,is,ie,js,je,ks,ke,xi_y_local,igs,ige,jgs,jge,kgs,kge,2,1)
        call fd1(xi_z_eta,is,ie,js,je,ks,ke,xi_z_local,igs,ige,jgs,jge,kgs,kge,2,1)
        call fd1(eta_x_eta,is,ie,js,je,ks,ke,eta_x_local,igs,ige,jgs,jge,kgs,kge,2,1)
        call fd1(eta_y_eta,is,ie,js,je,ks,ke,eta_y_local,igs,ige,jgs,jge,kgs,kge,2,1)
        call fd1(eta_z_eta,is,ie,js,je,ks,ke,eta_z_local,igs,ige,jgs,jge,kgs,kge,2,1)
        call fd1(phi_x_eta,is,ie,js,je,ks,ke,phi_x_local,igs,ige,jgs,jge,kgs,kge,2,1)
        call fd1(phi_y_eta,is,ie,js,je,ks,ke,phi_y_local,igs,ige,jgs,jge,kgs,kge,2,1)
        call fd1(phi_z_eta,is,ie,js,je,ks,ke,phi_z_local,igs,ige,jgs,jge,kgs,kge,2,1)

        call fd1(xi_x_phi,is,ie,js,je,ks,ke,xi_x_local,igs,ige,jgs,jge,kgs,kge,3,1)
        call fd1(xi_y_phi,is,ie,js,je,ks,ke,xi_y_local,igs,ige,jgs,jge,kgs,kge,3,1)
        call fd1(xi_z_phi,is,ie,js,je,ks,ke,xi_z_local,igs,ige,jgs,jge,kgs,kge,3,1)
        call fd1(eta_x_phi,is,ie,js,je,ks,ke,eta_x_local,igs,ige,jgs,jge,kgs,kge,3,1)
        call fd1(eta_y_phi,is,ie,js,je,ks,ke,eta_y_local,igs,ige,jgs,jge,kgs,kge,3,1)
        call fd1(eta_z_phi,is,ie,js,je,ks,ke,eta_z_local,igs,ige,jgs,jge,kgs,kge,3,1)
        call fd1(phi_x_phi,is,ie,js,je,ks,ke,phi_x_local,igs,ige,jgs,jge,kgs,kge,3,1)
        call fd1(phi_y_phi,is,ie,js,je,ks,ke,phi_y_local,igs,ige,jgs,jge,kgs,kge,3,1)
        call fd1(phi_z_phi,is,ie,js,je,ks,ke,phi_z_local,igs,ige,jgs,jge,kgs,kge,3,1)

        do k=ks,ke
            do j=js,je
                do i=is,ie
                    xi_xx(i,j,k)   = xi_x(i,j,k)*xi_x_xi(i,j,k) + eta_x(i,j,k)*xi_x_eta(i,j,k) + phi_x(i,j,k)*xi_x_phi(i,j,k)
                    xi_yy(i,j,k)   = xi_y(i,j,k)*xi_y_xi(i,j,k) + eta_y(i,j,k)*xi_y_eta(i,j,k) + phi_y(i,j,k)*xi_y_phi(i,j,k)
                    xi_zz(i,j,k)   = xi_z(i,j,k)*xi_z_xi(i,j,k) + eta_z(i,j,k)*xi_z_eta(i,j,k) + phi_z(i,j,k)*xi_z_phi(i,j,k)
                    xi_xy(i,j,k)   = xi_y(i,j,k)*xi_x_xi(i,j,k) + eta_y(i,j,k)*xi_x_eta(i,j,k) + phi_y(i,j,k)*xi_x_phi(i,j,k)
                    xi_yz(i,j,k)   = xi_z(i,j,k)*xi_y_xi(i,j,k) + eta_z(i,j,k)*xi_y_eta(i,j,k) + phi_z(i,j,k)*xi_y_phi(i,j,k)
                    xi_xz(i,j,k)   = xi_z(i,j,k)*xi_x_xi(i,j,k) + eta_z(i,j,k)*xi_x_eta(i,j,k) + phi_z(i,j,k)*xi_x_phi(i,j,k)
                    eta_xx(i,j,k)  = xi_x(i,j,k)*eta_x_xi(i,j,k) + eta_x(i,j,k)*eta_x_eta(i,j,k) + phi_x(i,j,k)*eta_x_phi(i,j,k)
                    eta_yy(i,j,k)  = xi_y(i,j,k)*eta_y_xi(i,j,k) + eta_y(i,j,k)*eta_y_eta(i,j,k) + phi_y(i,j,k)*eta_y_phi(i,j,k)
                    eta_zz(i,j,k)  = xi_z(i,j,k)*eta_z_xi(i,j,k) + eta_z(i,j,k)*eta_z_eta(i,j,k) + phi_z(i,j,k)*eta_z_phi(i,j,k)
                    eta_xy(i,j,k)  = xi_y(i,j,k)*eta_x_xi(i,j,k) + eta_y(i,j,k)*eta_x_eta(i,j,k) + phi_y(i,j,k)*eta_x_phi(i,j,k)
                    eta_yz(i,j,k)  = xi_z(i,j,k)*eta_y_xi(i,j,k) + eta_z(i,j,k)*eta_y_eta(i,j,k) + phi_z(i,j,k)*eta_y_phi(i,j,k)
                    eta_xz(i,j,k)  = xi_z(i,j,k)*eta_x_xi(i,j,k) + eta_z(i,j,k)*eta_x_eta(i,j,k) + phi_z(i,j,k)*eta_x_phi(i,j,k)
                    phi_xx(i,j,k)  = xi_x(i,j,k)*phi_x_xi(i,j,k) + eta_x(i,j,k)*phi_x_eta(i,j,k) + phi_x(i,j,k)*phi_x_phi(i,j,k)
                    phi_yy(i,j,k)  = xi_y(i,j,k)*phi_y_xi(i,j,k) + eta_y(i,j,k)*phi_y_eta(i,j,k) + phi_y(i,j,k)*phi_y_phi(i,j,k)
                    phi_zz(i,j,k)  = xi_z(i,j,k)*phi_z_xi(i,j,k) + eta_z(i,j,k)*phi_z_eta(i,j,k) + phi_z(i,j,k)*phi_z_phi(i,j,k)
                    phi_xy(i,j,k)  = xi_y(i,j,k)*phi_x_xi(i,j,k) + eta_y(i,j,k)*phi_x_eta(i,j,k) + phi_y(i,j,k)*phi_x_phi(i,j,k)
                    phi_yz(i,j,k)  = xi_z(i,j,k)*phi_y_xi(i,j,k) + eta_z(i,j,k)*phi_y_eta(i,j,k) + phi_z(i,j,k)*phi_y_phi(i,j,k)
                    phi_xz(i,j,k)  = xi_z(i,j,k)*phi_x_xi(i,j,k) + eta_z(i,j,k)*phi_x_eta(i,j,k) + phi_z(i,j,k)*phi_x_phi(i,j,k)
                enddo
            enddo
        enddo

        deallocate(xi_x_xi)
        deallocate(xi_x_eta)
        deallocate(xi_x_phi)
        deallocate(xi_y_xi)
        deallocate(xi_y_eta)
        deallocate(xi_y_phi)
        deallocate(xi_z_xi)
        deallocate(xi_z_eta)
        deallocate(xi_z_phi)
        deallocate(eta_x_xi)
        deallocate(eta_x_eta)
        deallocate(eta_x_phi)
        deallocate(eta_y_xi)
        deallocate(eta_y_eta)
        deallocate(eta_y_phi)
        deallocate(eta_z_xi)
        deallocate(eta_z_eta)
        deallocate(eta_z_phi)
        deallocate(phi_x_xi)
        deallocate(phi_x_eta)
        deallocate(phi_x_phi)
        deallocate(phi_y_xi)
        deallocate(phi_y_eta)
        deallocate(phi_y_phi)
        deallocate(phi_z_xi)
        deallocate(phi_z_eta)
        deallocate(phi_z_phi)
        deallocate(xi_x_local)
        deallocate(xi_y_local)
        deallocate(xi_z_local)
        deallocate(eta_x_local)
        deallocate(eta_y_local)
        deallocate(eta_z_local)
        deallocate(phi_x_local)
        deallocate(phi_y_local)
        deallocate(phi_z_local)
    end subroutine compute_convariant_metrics_3d

    subroutine delivery_by_dmda()
        implicit none
        Vec :: XIXX, XIYY, XIZZ, XIXY, XIXZ, XIYZ
        Vec :: ETAXX,ETAYY,ETAZZ,ETAXY,ETAXZ,ETAYZ
        Vec :: PHIXX,PHIYY,PHIZZ,PHIXY,PHIXZ,PHIYZ
        Vec :: XIXX_local, XIYY_local, XIZZ_local, XIXY_local, XIXZ_local, XIYZ_local
        Vec :: ETAXX_local,ETAYY_local,ETAZZ_local,ETAXY_local,ETAXZ_local,ETAYZ_local
        Vec :: PHIXX_local,PHIYY_local,PHIZZ_local,PHIXY_local,PHIXZ_local,PHIYZ_local
        PetscScalar, pointer :: tmp(:, :, :)

        call VecDuplicate(XIX, XIXX, ierr)
        call VecDuplicate(XIX, XIYY, ierr)
        call VecDuplicate(XIX, XIZZ, ierr)
        call VecDuplicate(XIX, XIXY, ierr)
        call VecDuplicate(XIX, XIXZ, ierr)
        call VecDuplicate(XIX, XIYZ, ierr)
        call VecDuplicate(XIX, ETAXX, ierr)
        call VecDuplicate(XIX, ETAYY, ierr)
        call VecDuplicate(XIX, ETAZZ, ierr)
        call VecDuplicate(XIX, ETAXY, ierr)
        call VecDuplicate(XIX, ETAXZ, ierr)
        call VecDuplicate(XIX, ETAYZ, ierr)
        call VecDuplicate(XIX, PHIXX, ierr)
        call VecDuplicate(XIX, PHIYY, ierr)
        call VecDuplicate(XIX, PHIZZ, ierr)
        call VecDuplicate(XIX, PHIXY, ierr)
        call VecDuplicate(XIX, PHIXZ, ierr)
        call VecDuplicate(XIX, PHIYZ, ierr)

        call VecDuplicate(XIX_local, XIXX_local, ierr)
        call VecDuplicate(XIX_local, XIYY_local, ierr)
        call VecDuplicate(XIX_local, XIZZ_local, ierr)
        call VecDuplicate(XIX_local, XIXY_local, ierr)
        call VecDuplicate(XIX_local, XIXZ_local, ierr)
        call VecDuplicate(XIX_local, XIYZ_local, ierr)
        call VecDuplicate(XIX_local, ETAXX_local, ierr)
        call VecDuplicate(XIX_local, ETAYY_local, ierr)
        call VecDuplicate(XIX_local, ETAZZ_local, ierr)
        call VecDuplicate(XIX_local, ETAXY_local, ierr)
        call VecDuplicate(XIX_local, ETAXZ_local, ierr)
        call VecDuplicate(XIX_local, ETAYZ_local, ierr)
        call VecDuplicate(XIX_local, PHIXX_local, ierr)
        call VecDuplicate(XIX_local, PHIYY_local, ierr)
        call VecDuplicate(XIX_local, PHIZZ_local, ierr)
        call VecDuplicate(XIX_local, PHIXY_local, ierr)
        call VecDuplicate(XIX_local, PHIXZ_local, ierr)
        call VecDuplicate(XIX_local, PHIYZ_local, ierr)

        ! 整合入Global Vector,并分配给Local Vector
        call DMDAVecGetArrayF90(DA, XIXX, tmp, ierr)
        tmp(:,:,:)=xi_xx(:,:,:)
        call DMDAVecRestoreArrayF90(DA, XIXX, tmp, ierr)
        call DMGlobalToLocalBegin(DA, XIXX, INSERT_VALUES, XIXX_local, ierr)

        call DMDAVecGetArrayF90(DA, XIYY, tmp, ierr)
        tmp(:,:,:)=xi_yy(:,:,:)
        call DMDAVecRestoreArrayF90(DA, XIYY, tmp, ierr)
        call DMGlobalToLocalEnd(DA, XIXX, INSERT_VALUES, XIXX_local, ierr)
        call DMGlobalToLocalBegin(DA, XIYY, INSERT_VALUES, XIYY_local, ierr)

        call DMDAVecGetArrayF90(DA, XIZZ, tmp, ierr)
        tmp(:,:,:)=xi_zz(:,:,:)
        call DMDAVecRestoreArrayF90(DA, XIZZ, tmp, ierr)
        call DMGlobalToLocalEnd(DA, XIYY, INSERT_VALUES, XIYY_local, ierr)
        call DMGlobalToLocalBegin(DA, XIZZ, INSERT_VALUES, XIZZ_local, ierr)

        call DMDAVecGetArrayF90(DA, XIXY, tmp, ierr)
        tmp(:,:,:)=xi_xy(:,:,:)
        call DMDAVecRestoreArrayF90(DA, XIXY, tmp, ierr)
        call DMGlobalToLocalEnd(DA, XIZZ, INSERT_VALUES, XIZZ_local, ierr)
        call DMGlobalToLocalBegin(DA, XIXY, INSERT_VALUES, XIXY_local, ierr)

        call DMDAVecGetArrayF90(DA, XIXZ, tmp, ierr)
        tmp(:,:,:)=xi_xz(:,:,:)
        call DMDAVecRestoreArrayF90(DA, XIXZ, tmp, ierr)
        call DMGlobalToLocalEnd(DA, XIXY, INSERT_VALUES, XIXY_local, ierr)
        call DMGlobalToLocalBegin(DA, XIXZ, INSERT_VALUES, XIXZ_local, ierr)

        call DMDAVecGetArrayF90(DA, XIYZ, tmp, ierr)
        tmp(:,:,:)=xi_yz(:,:,:)
        call DMDAVecRestoreArrayF90(DA, XIYZ, tmp, ierr)
        call DMGlobalToLocalEnd(DA, XIXZ, INSERT_VALUES, XIXZ_local, ierr)
        call DMGlobalToLocalBegin(DA, XIYZ, INSERT_VALUES, XIYZ_local, ierr)

        call DMDAVecGetArrayF90(DA, ETAXX, tmp, ierr)
        tmp(:,:,:)=eta_xx(:,:,:)
        call DMDAVecRestoreArrayF90(DA, ETAXX, tmp, ierr)
        call DMGlobalToLocalEnd(DA, XIYZ, INSERT_VALUES, XIYZ_local, ierr)
        call DMGlobalToLocalBegin(DA, ETAXX, INSERT_VALUES, ETAXX_local, ierr)

        call DMDAVecGetArrayF90(DA, ETAYY, tmp, ierr)
        tmp(:,:,:)=eta_yy(:,:,:)
        call DMDAVecRestoreArrayF90(DA, ETAYY, tmp, ierr)
        call DMGlobalToLocalEnd(DA, ETAXX, INSERT_VALUES, ETAXX_local, ierr)
        call DMGlobalToLocalBegin(DA, ETAYY, INSERT_VALUES, ETAYY_local, ierr)

        call DMDAVecGetArrayF90(DA, ETAZZ, tmp, ierr)
        tmp(:,:,:)=eta_zz(:,:,:)
        call DMDAVecRestoreArrayF90(DA, ETAZZ, tmp, ierr)
        call DMGlobalToLocalEnd(DA, ETAYY, INSERT_VALUES, ETAYY_local, ierr)
        call DMGlobalToLocalBegin(DA, ETAZZ, INSERT_VALUES, ETAZZ_local, ierr)

        call DMDAVecGetArrayF90(DA, ETAXY, tmp, ierr)
        tmp(:,:,:)=eta_xy(:,:,:)
        call DMDAVecRestoreArrayF90(DA, ETAXY, tmp, ierr)
        call DMGlobalToLocalEnd(DA, ETAZZ, INSERT_VALUES, ETAZZ_local, ierr)
        call DMGlobalToLocalBegin(DA, ETAXY, INSERT_VALUES, ETAXY_local, ierr)

        call DMDAVecGetArrayF90(DA, ETAXZ, tmp, ierr)
        tmp(:,:,:)=eta_xz(:,:,:)
        call DMDAVecRestoreArrayF90(DA, ETAXZ, tmp, ierr)
        call DMGlobalToLocalEnd(DA, ETAXY, INSERT_VALUES, ETAXY_local, ierr)
        call DMGlobalToLocalBegin(DA, ETAXZ, INSERT_VALUES, ETAXZ_local, ierr)

        call DMDAVecGetArrayF90(DA, ETAYZ, tmp, ierr)
        tmp(:,:,:)=eta_yz(:,:,:)
        call DMDAVecRestoreArrayF90(DA, ETAYZ, tmp, ierr)
        call DMGlobalToLocalEnd(DA, ETAXZ, INSERT_VALUES, ETAXZ_local, ierr)
        call DMGlobalToLocalBegin(DA, ETAYZ, INSERT_VALUES, ETAYZ_local, ierr)

        call DMDAVecGetArrayF90(DA, PHIXX, tmp, ierr)
        tmp(:,:,:)=phi_xx(:,:,:)
        call DMDAVecRestoreArrayF90(DA, PHIXX, tmp, ierr)
        call DMGlobalToLocalEnd(DA, ETAYZ, INSERT_VALUES, ETAYZ_local, ierr)
        call DMGlobalToLocalBegin(DA, PHIXX, INSERT_VALUES, PHIXX_local, ierr)

        call DMDAVecGetArrayF90(DA, PHIYY, tmp, ierr)
        tmp(:,:,:)=phi_yy(:,:,:)
        call DMDAVecRestoreArrayF90(DA, PHIYY, tmp, ierr)
        call DMGlobalToLocalEnd(DA, PHIXX, INSERT_VALUES, PHIXX_local, ierr)
        call DMGlobalToLocalBegin(DA, PHIYY, INSERT_VALUES, PHIYY_local, ierr)

        call DMDAVecGetArrayF90(DA, PHIZZ, tmp, ierr)
        tmp(:,:,:)=phi_zz(:,:,:)
        call DMDAVecRestoreArrayF90(DA, PHIZZ, tmp, ierr)
        call DMGlobalToLocalEnd(DA, PHIYY, INSERT_VALUES, PHIYY_local, ierr)
        call DMGlobalToLocalBegin(DA, PHIZZ, INSERT_VALUES, PHIZZ_local, ierr)

        call DMDAVecGetArrayF90(DA, PHIXY, tmp, ierr)
        tmp(:,:,:)=phi_xy(:,:,:)
        call DMDAVecRestoreArrayF90(DA, PHIXY, tmp, ierr)
        call DMGlobalToLocalEnd(DA, PHIZZ, INSERT_VALUES, PHIZZ_local, ierr)
        call DMGlobalToLocalBegin(DA, PHIXY, INSERT_VALUES, PHIXY_local, ierr)

        call DMDAVecGetArrayF90(DA, PHIXZ, tmp, ierr)
        tmp(:,:,:)=phi_xz(:,:,:)
        call DMDAVecRestoreArrayF90(DA, PHIXZ, tmp, ierr)
        call DMGlobalToLocalEnd(DA, PHIXY, INSERT_VALUES, PHIXY_local, ierr)
        call DMGlobalToLocalBegin(DA, PHIXZ, INSERT_VALUES, PHIXZ_local, ierr)

        call DMDAVecGetArrayF90(DA, PHIYZ, tmp, ierr)
        tmp(:,:,:)=phi_yz(:,:,:)
        call DMDAVecRestoreArrayF90(DA, PHIYZ, tmp, ierr)
        call DMGlobalToLocalEnd(DA, PHIXZ, INSERT_VALUES, PHIXZ_local, ierr)
        call DMGlobalToLocalBegin(DA, PHIYZ, INSERT_VALUES, PHIYZ_local, ierr)
        call DMGlobalToLocalEnd(DA, PHIYZ, INSERT_VALUES, PHIYZ_local, ierr)

        ! 重新分配数组
        deallocate(xi_x)
        deallocate(xi_y)
        deallocate(xi_z)
        deallocate(xi_xx)
        deallocate(xi_yy)
        deallocate(xi_zz)
        deallocate(xi_xy)
        deallocate(xi_xz)
        deallocate(xi_yz)
        deallocate(eta_x)
        deallocate(eta_y)
        deallocate(eta_z)
        deallocate(eta_xx)
        deallocate(eta_yy)
        deallocate(eta_zz)
        deallocate(eta_xy)
        deallocate(eta_xz)
        deallocate(eta_yz)
        deallocate(phi_x)
        deallocate(phi_y)
        deallocate(phi_z)
        deallocate(phi_xx)
        deallocate(phi_yy)
        deallocate(phi_zz)
        deallocate(phi_xy)
        deallocate(phi_xz)
        deallocate(phi_yz)

        allocate(xi_x(igs:ige, jgs:jge, kgs:kge))
        allocate(xi_y(igs:ige, jgs:jge, kgs:kge))
        allocate(xi_z(igs:ige, jgs:jge, kgs:kge))
        allocate(xi_xx(igs:ige, jgs:jge, kgs:kge))
        allocate(xi_yy(igs:ige, jgs:jge, kgs:kge))
        allocate(xi_zz(igs:ige, jgs:jge, kgs:kge))
        allocate(xi_xy(igs:ige, jgs:jge, kgs:kge))
        allocate(xi_xz(igs:ige, jgs:jge, kgs:kge))
        allocate(xi_yz(igs:ige, jgs:jge, kgs:kge))
        allocate(eta_x(igs:ige, jgs:jge, kgs:kge))
        allocate(eta_y(igs:ige, jgs:jge, kgs:kge))
        allocate(eta_z(igs:ige, jgs:jge, kgs:kge))
        allocate(eta_xx(igs:ige, jgs:jge, kgs:kge))
        allocate(eta_yy(igs:ige, jgs:jge, kgs:kge))
        allocate(eta_zz(igs:ige, jgs:jge, kgs:kge))
        allocate(eta_xy(igs:ige, jgs:jge, kgs:kge))
        allocate(eta_xz(igs:ige, jgs:jge, kgs:kge))
        allocate(eta_yz(igs:ige, jgs:jge, kgs:kge))
        allocate(phi_x(igs:ige, jgs:jge, kgs:kge))
        allocate(phi_y(igs:ige, jgs:jge, kgs:kge))
        allocate(phi_z(igs:ige, jgs:jge, kgs:kge))
        allocate(phi_xx(igs:ige, jgs:jge, kgs:kge))
        allocate(phi_yy(igs:ige, jgs:jge, kgs:kge))
        allocate(phi_zz(igs:ige, jgs:jge, kgs:kge))
        allocate(phi_xy(igs:ige, jgs:jge, kgs:kge))
        allocate(phi_xz(igs:ige, jgs:jge, kgs:kge))
        allocate(phi_yz(igs:ige, jgs:jge, kgs:kge))

        ! 从Local Vector获得数组

        call DMDAVecGetArrayReadF90(DA, XIX_local, tmp, ierr)
        xi_x=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, XIX_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, XIY_local, tmp, ierr)
        xi_y=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, XIY_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, XIZ_local, tmp, ierr)
        xi_z=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, XIZ_local, tmp, ierr)

        call DMDAVecGetArrayReadF90(DA, XIXX_local, tmp, ierr)
        xi_xx=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, XIXX_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, XIYY_local, tmp, ierr)
        xi_yy=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, XIYY_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, XIZZ_local, tmp, ierr)
        xi_zz=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, XIZZ_local, tmp, ierr)

        call DMDAVecGetArrayReadF90(DA, XIXY_local, tmp, ierr)
        xi_xy=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, XIXY_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, XIXZ_local, tmp, ierr)
        xi_xz=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, XIXZ_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, XIYZ_local, tmp, ierr)
        xi_yz=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, XIYZ_local, tmp, ierr)

        call DMDAVecGetArrayReadF90(DA, ETAX_local, tmp, ierr)
        eta_x=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, ETAX_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, ETAY_local, tmp, ierr)
        eta_y=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, ETAY_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, ETAZ_local, tmp, ierr)
        eta_z=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, ETAZ_local, tmp, ierr)

        call DMDAVecGetArrayReadF90(DA, ETAXX_local, tmp, ierr)
        eta_xx=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, ETAXX_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, ETAYY_local, tmp, ierr)
        eta_yy=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, ETAYY_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, ETAZZ_local, tmp, ierr)
        eta_zz=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, ETAZZ_local, tmp, ierr)

        call DMDAVecGetArrayReadF90(DA, ETAXY_local, tmp, ierr)
        eta_xy=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, ETAXY_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, ETAXZ_local, tmp, ierr)
        eta_xz=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, ETAXZ_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, ETAYZ_local, tmp, ierr)
        eta_yz=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, ETAYZ_local, tmp, ierr)

        call DMDAVecGetArrayReadF90(DA, PHIX_local, tmp, ierr)
        phi_x=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, PHIX_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, PHIY_local, tmp, ierr)
        phi_y=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, PHIY_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, PHIZ_local, tmp, ierr)
        phi_z=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, PHIZ_local, tmp, ierr)

        call DMDAVecGetArrayReadF90(DA, PHIXX_local, tmp, ierr)
        phi_xx=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, PHIXX_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, PHIYY_local, tmp, ierr)
        phi_yy=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, PHIYY_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, PHIZZ_local, tmp, ierr)
        phi_zz=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, PHIZZ_local, tmp, ierr)

        call DMDAVecGetArrayReadF90(DA, PHIXY_local, tmp, ierr)
        phi_xy=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, PHIXY_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, PHIXZ_local, tmp, ierr)
        phi_xz=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, PHIXZ_local, tmp, ierr)
        call DMDAVecGetArrayReadF90(DA, PHIYZ_local, tmp, ierr)
        phi_yz=real(tmp)
        call DMDAVecRestoreArrayReadF90(DA, PHIYZ_local, tmp, ierr)

        call VecDestroy(XIXX,ierr) ;call VecDestroy(XIYY,ierr) ;call VecDestroy(XIZZ,ierr)
        call VecDestroy(XIXY,ierr) ;call VecDestroy(XIXZ,ierr) ;call VecDestroy(XIYZ,ierr)

        call VecDestroy(ETAXX,ierr);call VecDestroy(ETAYY,ierr);call VecDestroy(ETAZZ,ierr)
        call VecDestroy(ETAXY,ierr);call VecDestroy(ETAXZ,ierr);call VecDestroy(ETAYZ,ierr)

        call VecDestroy(PHIXX,ierr);call VecDestroy(PHIYY,ierr);call VecDestroy(PHIZZ,ierr)
        call VecDestroy(PHIXY,ierr);call VecDestroy(PHIXZ,ierr);call VecDestroy(PHIYZ,ierr)

        call VecDestroy(XIXX_local,ierr) ;call VecDestroy(XIYY_local,ierr)
        call VecDestroy(XIZZ_local,ierr) ;call VecDestroy(XIXY_local,ierr)
        call VecDestroy(XIXZ_local,ierr) ;call VecDestroy(XIYZ_local,ierr)

        call VecDestroy(ETAXX_local,ierr);call VecDestroy(ETAYY_local,ierr)
        call VecDestroy(ETAZZ_local,ierr);call VecDestroy(ETAXY_local,ierr)
        call VecDestroy(ETAXZ_local,ierr);call VecDestroy(ETAYZ_local,ierr)

        call VecDestroy(PHIXX_local,ierr);call VecDestroy(PHIYY_local,ierr)
        call VecDestroy(PHIZZ_local,ierr);call VecDestroy(PHIXY_local,ierr)
        call VecDestroy(PHIXZ_local,ierr);call VecDestroy(PHIYZ_local,ierr)

    end subroutine  delivery_by_dmda

    subroutine deallocate_memory()
        implicit none
        deallocate(xx)
        deallocate(yy)
        deallocate(zz)
        deallocate(jacobi)
        deallocate(x_xi)
        deallocate(y_xi)
        deallocate(z_xi)
        deallocate(x_eta)
        deallocate(y_eta)
        deallocate(z_eta)
        deallocate(x_phi)
        deallocate(y_phi)
        deallocate(z_phi)
        call VecDestroy(XIX, ierr)
        call VecDestroy(XIY, ierr)
        call VecDestroy(XIZ, ierr)
        call VecDestroy(ETAX, ierr)
        call VecDestroy(ETAY, ierr)
        call VecDestroy(ETAZ, ierr)
        call VecDestroy(PHIX, ierr)
        call VecDestroy(PHIY, ierr)
        call VecDestroy(PHIZ, ierr)
        call VecDestroy(XIX_local, ierr)
        call VecDestroy(XIY_local, ierr)
        call VecDestroy(XIZ_local, ierr)
        call VecDestroy(ETAX_local, ierr)
        call VecDestroy(ETAY_local, ierr)
        call VecDestroy(ETAZ_local, ierr)
        call VecDestroy(PHIX_local, ierr)
        call VecDestroy(PHIY_local, ierr)
        call VecDestroy(PHIZ_local, ierr)
    end subroutine deallocate_memory

end module mod_metrics
