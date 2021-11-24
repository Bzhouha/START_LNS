#include <slepc/finclude/slepc.h>

module bf_points
! -----------------------------------------------------------
!
!   这个模块计算流场中基本流的物理量对计算坐标的导数
!
!       call Tunas(comm) 计算基本流的物理量对计算坐标的导数的函数
!
!           1).call PrepareNet() 分配内存
!
!           2).call CatchAllTunas() 计算在计算坐标下的导数的函数
!
!           3).call FromIJKToXYZ() 将计算坐标转换到物理坐标下的函数
!
!               a).call TurnIToX(px,pi,pj,pk,ix,jx,kx) 一阶导数坐标转换函数
!
!               b).call TurnIIToXX(pxx,pi,pj,pk,...) 二阶导数坐标转换函数
!
!               c).call TurnIJToXY(pxy,pi,pj,pk,...) 混合导数坐标转换函数
!
!           4).call InsertToBF() 将数据整合到BF结构中
!
!               call insert(obj,arr) 将数组形式的数据插入bf结构的函数
!
!           5).call DropTheWaste() 释放内存
!
! 			6).call TunasSayHi(comm) 输出本模块运行结束信息
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
	public :: Tunas
	contains
	subroutine Tunas(comm)
		implicit none
		PetscInt,intent(in) :: comm
		call PrepareNet()
		call CatchAllTunas()
		call FromIJKToXYZ()
		call InsertToBF()
		call DropTheWaste()
		call TunasSayHi(comm)
	end subroutine Tunas
	subroutine PrepareNet()
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
		allocate(qqi(0:4,is:ie,js:je,ks:ke))
		allocate(qqj(0:4,is:ie,js:je,ks:ke))
		allocate(qqk(0:4,is:ie,js:je,ks:ke))
		allocate(qqii(0:4,is:ie,js:je,ks:ke))
		allocate(qqjj(0:4,is:ie,js:je,ks:ke))
		allocate(qqkk(0:4,is:ie,js:je,ks:ke))
		allocate(qqij(0:4,is:ie,js:je,ks:ke))
		allocate(qqjk(0:4,is:ie,js:je,ks:ke))
		allocate(qqik(0:4,is:ie,js:je,ks:ke))
		allocate(qq_i_local_array(0:4,igs:ige,jgs:jge,kgs:kge))
		allocate(qq_j_local_array(0:4,igs:ige,jgs:jge,kgs:kge))
		allocate(qq_k_local_array(0:4,igs:ige,jgs:jge,kgs:kge))
	end subroutine PrepareNet
	subroutine CatchAllTunas()
		use mod_difference
		implicit none
		PetscScalar, pointer :: tmp(:, :, :, :)
		call fd1(qqi,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,1,5)
		call fd1(qqj,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,2,5)
		call fd1(qqk,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,3,5)
		call fd2(qqii,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,1,5)
		call fd2(qqjj,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,2,5)
		call fd2(qqkk,is,ie,js,je,ks,ke,qq,igs,ige,jgs,jge,kgs,kge,3,5)
		call DMDAVecGetArrayF90(meshDA, QQ_I, tmp, ierr)
		tmp = qqi
		call DMDAVecRestoreArrayF90(meshDA, QQ_I, tmp, ierr)
		call DMGlobalToLocalBegin(meshDA, QQ_I, INSERT_VALUES, QQ_I_local, ierr)
		call DMGlobalToLocalEnd(meshDA, QQ_I, INSERT_VALUES, QQ_I_local, ierr)
		call DMDAVecGetArrayF90(meshDA, QQ_J, tmp, ierr)
		tmp = qqj
		call DMDAVecRestoreArrayF90(meshDA, QQ_J, tmp, ierr)
		call DMGlobalToLocalBegin(meshDA, QQ_J, INSERT_VALUES, QQ_J_local, ierr)
		call DMGlobalToLocalEnd(meshDA, QQ_J, INSERT_VALUES, QQ_J_local, ierr)
		call DMDAVecGetArrayF90(meshDA, QQ_K, tmp, ierr)
		tmp = qqk
		call DMDAVecRestoreArrayF90(meshDA, QQ_K, tmp, ierr)
		call DMGlobalToLocalBegin(meshDA, QQ_K, INSERT_VALUES, QQ_K_local, ierr)
		call DMGlobalToLocalEnd(meshDA, QQ_K, INSERT_VALUES, QQ_K_local, ierr)
		call DMDAVecGetArrayReadF90(meshDA, QQ_I_local, tmp, ierr)
		qq_i_local_array=tmp(:,igs:ige,jgs:jge,kgs:kge)
		call DMDAVecRestoreArrayReadF90(meshDA, QQ_I_local, tmp, ierr)
		call DMDAVecGetArrayReadF90(meshDA, QQ_J_local, tmp, ierr)
		qq_j_local_array=tmp(:,igs:ige,jgs:jge,kgs:kge)
		call DMDAVecRestoreArrayReadF90(meshDA, QQ_J_local, tmp, ierr)
		call DMDAVecGetArrayReadF90(meshDA, QQ_K_local, tmp, ierr)
		qq_k_local_array=tmp(:,igs:ige,jgs:jge,kgs:kge)
		call DMDAVecRestoreArrayReadF90(meshDA, QQ_K_local, tmp, ierr)
		call fd1(qqij,is,ie,js,je,ks,ke,qq_i_local_array,igs,ige,jgs,jge,kgs,kge,2,5)
		call fd1(qqik,is,ie,js,je,ks,ke,qq_k_local_array,igs,ige,jgs,jge,kgs,kge,1,5)
		call fd1(qqjk,is,ie,js,je,ks,ke,qq_j_local_array,igs,ige,jgs,jge,kgs,kge,3,5)
		deallocate(qq_i_local_array)
		deallocate(qq_j_local_array)
		deallocate(qq_k_local_array)
		call VecDestroy(QQ_I,ierr)
		call VecDestroy(QQ_J,ierr)
		call VecDestroy(QQ_K,ierr)
		call VecDestroy(QQ_I_local,ierr)
		call VecDestroy(QQ_J_local,ierr)
		call VecDestroy(QQ_K_local,ierr)
	end subroutine CatchAllTunas
	subroutine FromIJKToXYZ()
		implicit none
		integer :: l
		allocate(qqx(0:4,is:ie,js:je,ks:ke))
		allocate(qqy(0:4,is:ie,js:je,ks:ke))
		allocate(qqz(0:4,is:ie,js:je,ks:ke))
		allocate(qqxx(0:4,is:ie,js:je,ks:ke))
		allocate(qqyy(0:4,is:ie,js:je,ks:ke))
		allocate(qqzz(0:4,is:ie,js:je,ks:ke))
		allocate(qqxy(0:4,is:ie,js:je,ks:ke))
		allocate(qqxz(0:4,is:ie,js:je,ks:ke))
		allocate(qqyz(0:4,is:ie,js:je,ks:ke))
		qqx=0.0d0;qqz=0.0d0;qqy=0.0d0
		qqxx=0.0d0;qqzz=0.0d0;qqyy=0.0d0
		qqxy=0.0d0;qqyz=0.0d0;qqxz=0.0d0
		do l=0,4
			call TurnIToX(qqx(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),xi_x,eta_x,phi_x)
			call TurnIToX(qqy(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),xi_y,eta_y,phi_y)
			call TurnIToX(qqz(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),xi_z,eta_z,phi_z)
		enddo 
		do l=1,4
			call TurnIIToXX(qqxx(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),&
			qqii(l,:,:,:),qqjj(l,:,:,:),qqkk(l,:,:,:),&
			qqij(l,:,:,:),qqik(l,:,:,:),qqjk(l,:,:,:),&
			xi_xx,eta_xx,phi_xx,xi_x,eta_x,phi_x)
			call TurnIIToXX(qqyy(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),&
			qqii(l,:,:,:),qqjj(l,:,:,:),qqkk(l,:,:,:),&
			qqij(l,:,:,:),qqik(l,:,:,:),qqjk(l,:,:,:),&
			xi_yy,eta_yy,phi_yy,xi_y,eta_y,phi_y)
			call TurnIIToXX(qqzz(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),&
			qqii(l,:,:,:),qqjj(l,:,:,:),qqkk(l,:,:,:),&
			qqij(l,:,:,:),qqik(l,:,:,:),qqjk(l,:,:,:),&
			xi_zz,eta_zz,phi_zz,xi_z,eta_z,phi_z)
		enddo
		do l=1,3
			call TurnIJToXY(qqxy(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),&
			qqii(l,:,:,:),qqjj(l,:,:,:),qqkk(l,:,:,:),&
			qqij(l,:,:,:),qqik(l,:,:,:),qqjk(l,:,:,:),&
			xi_xy,eta_xy,phi_xy,xi_x,xi_y,eta_x,eta_y,phi_x,phi_y)
			call TurnIJToXY(qqyz(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),&
			qqii(l,:,:,:),qqjj(l,:,:,:),qqkk(l,:,:,:),&
			qqij(l,:,:,:),qqik(l,:,:,:),qqjk(l,:,:,:),&
			xi_yz,eta_yz,phi_yz,xi_y,xi_z,eta_y,eta_z,phi_y,phi_z)
			call TurnIJToXY(qqxz(l,:,:,:),qqi(l,:,:,:),qqj(l,:,:,:),qqk(l,:,:,:),&
			qqii(l,:,:,:),qqjj(l,:,:,:),qqkk(l,:,:,:),&
			qqij(l,:,:,:),qqik(l,:,:,:),qqjk(l,:,:,:),&
			xi_xz,eta_xz,phi_xz,xi_x,xi_z,eta_x,eta_z,phi_x,phi_z)
		enddo
	end subroutine FromIJKToXYZ
	subroutine InsertToBF()
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
	end subroutine InsertToBF
	subroutine DropTheWaste()
		implicit none
		deallocate(qq)
		deallocate(qqi);deallocate(qqj);deallocate(qqk)
		deallocate(qqii);deallocate(qqjj);deallocate(qqkk)
		deallocate(qqij);deallocate(qqik);deallocate(qqjk)
		deallocate(qqx);deallocate(qqy);deallocate(qqz)
		deallocate(qqxx);deallocate(qqyy);deallocate(qqzz)
		deallocate(qqxy);deallocate(qqxz);deallocate(qqyz)
	end subroutine DropTheWaste
	subroutine TunasSayHi(comm)
		implicit none 
		PetscInt,intent(in) :: comm 
		PetscErrorCode :: ierr  
		call PetscPrintf(comm," -----------------------------------\n",ierr)
		call PetscPrintf(comm,"         导数信息计算结束。      \n",ierr)
		call PetscPrintf(comm," -----------------------------------\n",ierr)
	end subroutine TunasSayHi
	elemental subroutine TurnIToX(px,pi,pj,pk,ix,jx,kx)
		implicit none
		real(R_P),intent(out) :: px
		real(R_P),intent(in) :: pi,pj,pk
		real(R_P),intent(in) :: ix,jx,kx
		px = ix*pi+jx*pj+kx*pk 
	end subroutine TurnIToX
	elemental subroutine TurnIIToXX(pxx,pi,pj,pk,pii,pjj,pkk,pij,pik,pjk,&
		iixx,jjxx,kkxx,ix,jx,kx)
		implicit none
		real(R_P),intent(out) :: pxx 
		real(R_P),intent(in) :: pi,pj,pk,pii,pjj,pkk,pij,pik,pjk
		real(R_P),intent(in) :: iixx,jjxx,kkxx,ix,jx,kx 
		pxx = iixx*pi+jjxx*pj+kkxx*pk+ix*ix*pii+jx*jx*pjj+kx*kx*pkk+&
		2.0d0*ix*jx*pij+2.0d0*ix*kx*pik+2.0d0*jx*kx*pjk
	end subroutine TurnIIToXX
	elemental subroutine TurnIJToXY(pxy,pi,pj,pk,pii,pjj,pkk,pij,pik,pjk,&
		iixy,jjxy,kkxy,ix,iy,jx,jy,kx,ky)
		implicit none
		real(R_P),intent(out) :: pxy  
		real(R_P),intent(in) :: pi,pj,pk,pii,pjj,pkk,pij,pik,pjk
		real(R_P),intent(in) :: iixy,jjxy,kkxy,ix,iy,jx,jy,kx,ky 
		pxy = iixy*pi+jjxy*pj+kkxy*pk+ix*iy*pii+jx*jy*pjj+kx*ky*pkk+&
		(ix*jy+iy*jx)*pij+(ix*ky+iy*kx)*pik+(jx*ky+jy*kx)*pjk
	end subroutine TurnIJToXY
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
end module bf_points