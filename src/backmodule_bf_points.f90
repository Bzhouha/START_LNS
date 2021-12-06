#include <slepc/finclude/slepc.h>

module bf_points
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
!           3).call fromIJKtoXYZ() 将计算坐标转换到物理坐标下的函数
!
!               a).call turnItoX(px,pi,pj,pk,ix,jx,kx) 一阶导数坐标转换函数
!
!               b).call turnIItoXX(pxx,pi,pj,pk,...) 二阶导数坐标转换函数
!
!               c).call turnIJtoXY(pxy,pi,pj,pk,...) 混合导数坐标转换函数
!
!           4).call insert_to_BF() 将数据整合到BF结构中
!
!               call insert(obj,arr) 将数组形式的数据插入bf结构的函数
!
!           5).call deallocate_memoru() 释放内存
!
! 			6).call print_info(comm) 输出本模块运行结束信息
!
! -----------------------------------------------------------
	use petsc
	use global_parameters
	use penf, only: R_P
	implicit none
	private
	real(R_P), allocatable, dimension(:, :, :, :) :: qq_i_local_array
	real(R_P), allocatable, dimension(:, :, :, :) :: qq_j_local_array
	real(R_P), allocatable, dimension(:, :, :, :) :: qq_k_local_array
	real(R_P), allocatable, dimension(:, :, :, :) :: qqii,qqjj,qqkk
	real(R_P), allocatable, dimension(:, :, :, :) :: qqij,qqik,qqjk
	real(R_P), allocatable, dimension(:, :, :, :) :: qqxx,qqyy,qqzz
	real(R_P), allocatable, dimension(:, :, :, :) :: qqxy,qqxz,qqyz
	real(R_P), allocatable, dimension(:, :, :, :) :: qqi,qqj,qqk
	real(R_P), allocatable, dimension(:, :, :, :) :: qqx,qqy,qqz
	Vec :: QQ_I_local,QQ_J_local,QQ_K_local
	PetscErrorCode :: ierr
	Vec :: QQ_I,QQ_J,QQ_K
	public :: partial_derivatives
	contains

	subroutine partial_derivatives(comm)
		implicit none
		PetscInt,intent(in) :: comm
		call allocate_memory()
		call partial_derivatives_on_IJK()
		call fromIJKtoXYZ()
		call insert_to_BF()
		call deallocate_memory()
		call print_info(comm)
		call MPI_Barrier(comm,ierr)
	end subroutine partial_derivatives

	subroutine allocate_memory()
		implicit none
		allocate(bf(is:ie,js:je,ks:ke))
		call DMGetGlobalVector(meshDA, QQ_I, ierr)
		call VecDuplicate(QQ_I, QQ_J, ierr)
		call VecDuplicate(QQ_I, QQ_K, ierr)
		call DMGetLocalVector(meshDA, QQ_I_local, ierr)
		call VecDuplicate(QQ_I_local, QQ_J_local, ierr)
		call VecDuplicate(QQ_I_local, QQ_K_local, ierr)
		call VecZeroEntries(QQ_I,ierr)
		call VecZeroEntries(QQ_J,ierr)
		call VecZeroEntries(QQ_K,ierr)
		call VecZeroEntries(QQ_I_local,ierr)
		call VecZeroEntries(QQ_J_local,ierr)
		call VecZeroEntries(QQ_K_local,ierr)
		allocate(qqi(5,is:ie,js:je,ks:ke))
		allocate(qqj(5,is:ie,js:je,ks:ke))
		allocate(qqk(5,is:ie,js:je,ks:ke))
		allocate(qqii(5,is:ie,js:je,ks:ke))
		allocate(qqjj(5,is:ie,js:je,ks:ke))
		allocate(qqkk(5,is:ie,js:je,ks:ke))
		allocate(qqij(5,is:ie,js:je,ks:ke))
		allocate(qqjk(5,is:ie,js:je,ks:ke))
		allocate(qqik(5,is:ie,js:je,ks:ke))
		allocate(qqx(5,is:ie,js:je,ks:ke))
		allocate(qqy(5,is:ie,js:je,ks:ke))
		allocate(qqz(5,is:ie,js:je,ks:ke))
		allocate(qqxx(5,is:ie,js:je,ks:ke))
		allocate(qqyy(5,is:ie,js:je,ks:ke))
		allocate(qqzz(5,is:ie,js:je,ks:ke))
		allocate(qqxy(5,is:ie,js:je,ks:ke))
		allocate(qqxz(5,is:ie,js:je,ks:ke))
		allocate(qqyz(5,is:ie,js:je,ks:ke))
		allocate(qq_i_local_array(5,igs:ige,jgs:jge,kgs:kge))
		allocate(qq_j_local_array(5,igs:ige,jgs:jge,kgs:kge))
		allocate(qq_k_local_array(5,igs:ige,jgs:jge,kgs:kge))
	end subroutine allocate_memory

	subroutine partial_derivatives_on_IJK()
		use mod_difference
		implicit none
		PetscScalar, pointer :: tmp(:, :, :, :)
		select case (lns_mode)
		case(0)
			call fd1(qqi,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,1,5)
			call fd1(qqj,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,2,5)
			call fd2(qqii,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,1,5)
			call fd2(qqjj,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,2,5)
			qqk=0.0d0;qqkk=0.0d0
		case(1)
			call fd1(qqi,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,1,5)
			call fd1(qqj,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,2,5)
			call fd1(qqk,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,3,5)
			call fd2(qqii,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,1,5)
			call fd2(qqjj,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,2,5)
			call fd2(qqkk,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,3,5)
		end select
		call DMDAVecGetArrayF90(meshDA, QQ_I, tmp, ierr)
		tmp(:,is:ie,js:je,ks:ke) = qqi(:,is:ie,js:je,ks:ke)
		call DMDAVecRestoreArrayF90(meshDA, QQ_I, tmp, ierr)
		call DMGlobalToLocalBegin(meshDA, QQ_I, INSERT_VALUES, QQ_I_local, ierr)
		call DMGlobalToLocalEnd(meshDA, QQ_I, INSERT_VALUES, QQ_I_local, ierr)
		call DMDAVecGetArrayF90(meshDA, QQ_J, tmp, ierr)
		tmp(:,is:ie,js:je,ks:ke) = qqj(:,is:ie,js:je,ks:ke)
		call DMDAVecRestoreArrayF90(meshDA, QQ_J, tmp, ierr)
		call DMGlobalToLocalBegin(meshDA, QQ_J, INSERT_VALUES, QQ_J_local, ierr)
		call DMGlobalToLocalEnd(meshDA, QQ_J, INSERT_VALUES, QQ_J_local, ierr)
		call DMDAVecGetArrayF90(meshDA, QQ_K, tmp, ierr)
		tmp(:,is:ie,js:je,ks:ke) = qqk(:,is:ie,js:je,ks:ke)
		call DMDAVecRestoreArrayF90(meshDA, QQ_K, tmp, ierr)
		call DMGlobalToLocalBegin(meshDA, QQ_K, INSERT_VALUES, QQ_K_local, ierr)
		call DMGlobalToLocalEnd(meshDA, QQ_K, INSERT_VALUES, QQ_K_local, ierr)
		call DMDAVecGetArrayReadF90(meshDA, QQ_I_local, tmp, ierr)
		qq_i_local_array(:,igs:ige,jgs:jge,kgs:kge)=tmp(:,igs:ige,jgs:jge,kgs:kge)
		call DMDAVecRestoreArrayReadF90(meshDA, QQ_I_local, tmp, ierr)
		call DMDAVecGetArrayReadF90(meshDA, QQ_J_local, tmp, ierr)
		qq_j_local_array(:,igs:ige,jgs:jge,kgs:kge)=tmp(:,igs:ige,jgs:jge,kgs:kge)
		call DMDAVecRestoreArrayReadF90(meshDA, QQ_J_local, tmp, ierr)
		call DMDAVecGetArrayReadF90(meshDA, QQ_K_local, tmp, ierr)
		qq_k_local_array(:,igs:ige,jgs:jge,kgs:kge)=tmp(:,igs:ige,jgs:jge,kgs:kge)
		call DMDAVecRestoreArrayReadF90(meshDA, QQ_K_local, tmp, ierr)
		select case (lns_mode)
		case(0)
			call fd1(qqij,is,ie,js,je,ks,ke,qq_i_local_array,igs,ige,jgs,jge,kgs,kge,2,5)
			qqik=0.0d0;qqjk=0.0d0
		case(1)
			call fd1(qqij,is,ie,js,je,ks,ke,qq_i_local_array,igs,ige,jgs,jge,kgs,kge,2,5)
			call fd1(qqik,is,ie,js,je,ks,ke,qq_k_local_array,igs,ige,jgs,jge,kgs,kge,1,5)
			call fd1(qqjk,is,ie,js,je,ks,ke,qq_j_local_array,igs,ige,jgs,jge,kgs,kge,3,5)
		end select
		deallocate(qq_i_local_array)
		deallocate(qq_j_local_array)
		deallocate(qq_k_local_array)
		call VecDestroy(QQ_I,ierr)
		call VecDestroy(QQ_J,ierr)
		call VecDestroy(QQ_K,ierr)
		call VecDestroy(QQ_I_local,ierr)
		call VecDestroy(QQ_J_local,ierr)
		call VecDestroy(QQ_K_local,ierr)
	end subroutine partial_derivatives_on_IJK

	subroutine fromIJKtoXYZ()
		implicit none
		integer :: l
		qqx=0.0d0;qqz=0.0d0;qqy=0.0d0
		qqxx=0.0d0;qqzz=0.0d0;qqyy=0.0d0
		qqxy=0.0d0;qqyz=0.0d0;qqxz=0.0d0
		do l=1,5
			call turnItoX(qqx(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),xi_x,eta_x,phi_x)
			call turnItoX(qqy(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),xi_y,eta_y,phi_y)
			call turnItoX(qqz(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),xi_z,eta_z,phi_z)
		enddo 
		do l=2,5
			call turnIItoXX(qqxx(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),&
			qqii(l,:,:,:),qqjj(l,:,:,:),qqkk(l,:,:,:),&
			qqij(l,:,:,:),qqik(l,:,:,:),qqjk(l,:,:,:),&
			xi_xx,eta_xx,phi_xx,xi_x,eta_x,phi_x)
			call turnIItoXX(qqyy(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),&
			qqii(l,:,:,:),qqjj(l,:,:,:),qqkk(l,:,:,:),&
			qqij(l,:,:,:),qqik(l,:,:,:),qqjk(l,:,:,:),&
			xi_yy,eta_yy,phi_yy,xi_y,eta_y,phi_y)
			call turnIItoXX(qqzz(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),&
			qqii(l,:,:,:),qqjj(l,:,:,:),qqkk(l,:,:,:),&
			qqij(l,:,:,:),qqik(l,:,:,:),qqjk(l,:,:,:),&
			xi_zz,eta_zz,phi_zz,xi_z,eta_z,phi_z)
		enddo
		do l=2,4
			call turnIJtoXY(qqxy(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),&
			qqii(l,:,:,:),qqjj(l,:,:,:),qqkk(l,:,:,:),&
			qqij(l,:,:,:),qqik(l,:,:,:),qqjk(l,:,:,:),&
			xi_xy,eta_xy,phi_xy,xi_x,xi_y,eta_x,eta_y,phi_x,phi_y)
			call turnIJtoXY(qqxz(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),&
			qqii(l,:,:,:),qqjj(l,:,:,:),qqkk(l,:,:,:),&
			qqij(l,:,:,:),qqik(l,:,:,:),qqjk(l,:,:,:),&
			xi_xz,eta_xz,phi_xz,xi_x,xi_z,eta_x,eta_z,phi_x,phi_z)
			call turnIJtoXY(qqyz(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),&
			qqii(l,:,:,:),qqjj(l,:,:,:),qqkk(l,:,:,:),&
			qqij(l,:,:,:),qqik(l,:,:,:),qqjk(l,:,:,:),&
			xi_yz,eta_yz,phi_yz,xi_y,xi_z,eta_y,eta_z,phi_y,phi_z)
		enddo
		select case (lns_mode)
		case(0)
			qqz=0.0d0;qqzz=0.0d0;qqxz=0.0d0;qqyz=0.0d0
		end select
	end subroutine fromIJKtoXYZ

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
		use mod_baseflow_org 
		implicit none
		real(R_P),dimension(5),intent(in) :: arr
		type(bf_flux_org_type),intent(out) :: obj
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
		deallocate(qqii);deallocate(qqjj);deallocate(qqkk)
		deallocate(qqij);deallocate(qqik);deallocate(qqjk)
		deallocate(qqx);deallocate(qqy);deallocate(qqz)
		deallocate(qqxx);deallocate(qqyy);deallocate(qqzz)
		deallocate(qqxy);deallocate(qqxz);deallocate(qqyz)
	end subroutine deallocate_memory

	subroutine print_info(comm)
		implicit none 
		PetscInt,intent(in) :: comm 
		PetscErrorCode :: ierr  
		call PetscPrintf(comm," -----------------------------------\n",ierr)
		call PetscPrintf(comm,"         导数信息计算结束。      \n",ierr)
		call PetscPrintf(comm," -----------------------------------\n",ierr)
	end subroutine print_info

	elemental subroutine turnItoX(fx,fi,fj,fk,ix,jx,kx)
		implicit none
		real(R_P),intent(out) :: fx
		real(R_P),intent(in) :: fi,fj,fk
		real(R_P),intent(in) :: ix,jx,kx
		fx = ix*fi+jx*fj+kx*fk 
	end subroutine turnItoX

	elemental subroutine turnIItoXX(fxx,fi,fj,fk,fii,fjj,fkk,fij,fik,fjk,&
		ixx,jxx,kxx,ix,jx,kx)
		implicit none
		real(R_P),intent(out) :: fxx 
		real(R_P),intent(in) :: fi,fj,fk,fii,fjj,fkk,fij,fik,fjk
		real(R_P),intent(in) :: ixx,jxx,kxx,ix,jx,kx 
		fxx = ixx*fi+jxx*fj+kxx*fk+ix*ix*fii+jx*jx*fjj+kx*kx*fkk+&
		2.0d0*ix*jx*fij+2.0d0*ix*kx*fik+2.0d0*jx*kx*fjk
	end subroutine turnIItoXX

	elemental subroutine turnIJtoXY(fxy,fi,fj,fk,fii,fjj,fkk,fij,fik,fjk,&
		ixy,jxy,kxy,ix,iy,jx,jy,kx,ky)
		implicit none
		real(R_P),intent(out) :: fxy  
		real(R_P),intent(in) :: fi,fj,fk,fii,fjj,fkk,fij,fik,fjk
		real(R_P),intent(in) :: ixy,jxy,kxy,ix,iy,jx,jy,kx,ky 
		fxy = ixy*fi+jxy*fj+kxy*fk+ix*iy*fii+jx*jy*fjj+kx*ky*fkk+&
		(ix*jy+iy*jx)*fij+(ix*ky+iy*kx)*fik+(jx*ky+jy*kx)*fjk
	end subroutine turnIJtoXY
	
end module bf_points