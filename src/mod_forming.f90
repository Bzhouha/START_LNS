!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_forming
! ------------------------------------------------------------------
!
!  这个模块生成最终的大矩阵。
!
!       1.call dolphin_coming(comm) 免矩阵形式的矩阵生成函数
!       
!           call Mymult(A, X, F, ierr) 免矩阵需要的矩阵向量乘法函数
!
!       2.call shark_coming(comm) 矩阵形式的矩阵生成函数
!
!           call shark_growing_up() 填充数据函数
!
!               1).call shark_eat_shrimps(i,j,k) 边界部分的填充数据函数
!
!               2).call shark_eat_sardine(i,j,k) 内部部分的填充数据函数
!
! ------------------------------------------------------------------
	use penf, only: R_P
	use mod_parameters
	use mod_flowtype
	use mod_cubes
	use petsc
	implicit none
	public :: dolphin_coming, shark_coming
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
	subroutine dolphin_coming(comm)
		implicit none
		PetscInt,intent(in) :: comm
		PetscErrorCode :: ierr
		PetscInt :: ls
		call VecGetLocalSize(turtle,ls,ierr)
		call MatCreateShell(comm,ls,ls,PETSC_DETERMINE,PETSC_DETERMINE,PETSC_NULL_INTEGER,Dolphin,ierr)
		call MatShellSetOperation(Dolphin,MATOP_MULT,Mymult,ierr)
		call MatAssemblyBegin(Dolphin,MAT_FINAL_ASSEMBLY,ierr)
		call MatAssemblyEnd(Dolphin,MAT_FINAL_ASSEMBLY,ierr)
		call dolphin_say_hi(comm)
	end subroutine dolphin_coming

	subroutine Mymult(A, X, F, ierr)
		use mod_mftools
		implicit none
		PetscScalar,pointer :: F_r(:,:,:,:),X_r(:,:,:,:)
		complex(R_P), dimension(5,15) :: crab
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
						crab=0.0d0
						call Jor%get_adorned_cubes(i,j,k)
						call f3d1r(crab(:,1),x(:,i-2:i,j,k)) ! A1p
						call f3d1f(crab(:,2),x(:,i:i+2,j,k)) ! A1m
						call f5d1(crab(:,3),x(:,i-2:i+2,j,k)) ! A2
						call f3d1r(crab(:,4),x(:,i,j-2:j,k)) ! B1p
						call f3d1f(crab(:,5),x(:,i,j:j+2,k)) ! B1m
						call f5d1(crab(:,6),x(:,i,j-2:j+2,k)) ! B2
						call f3d1r(crab(:,7),x(:,i,j,k-2:k)) ! C1p
						call f3d1f(crab(:,8),x(:,i,j,k:k+2)) ! C1m
						call f5d1(crab(:,9),x(:,i,j,k-2:k+2)) ! C2
						call f5d11(crab(:,10),x(:,i-2:i+2,j,k)) ! Vxx
						call f5d11(crab(:,11),x(:,i,j-2:j+2,k)) ! Vyy
						call f5d11(crab(:,12),x(:,i,j,k-2:k+2)) ! Vzz
						call f5d12(crab(:,13),x(:,i-2:i+2,j-2:j+2,k)) ! Vxy
						call f5d12(crab(:,14),x(:,i-2:i+2,j,k-2:k+2)) ! Vxz
						call f5d12(crab(:,15),x(:,i,j-2:j+2,k-2:k+2)) ! Vyz
						F(:,i,j,k) = matmul(A_p,crab(:,1))+matmul(A_m,crab(:,2))+matmul(A_v,crab(:,3)) &
									+matmul(B_p,crab(:,4))+matmul(B_m,crab(:,5))+matmul(B_v,crab(:,6)) &
									+matmul(C_p,crab(:,7))+matmul(C_m,crab(:,8))+matmul(C_v,crab(:,9)) &
									+matmul(D,x(:,i,j,k)) &
									-matmul(Vxx,crab(:,10))-matmul(Vyy,crab(:,11))-matmul(Vzz,crab(:,12)) &
									-matmul(Vxy,crab(:,13))-matmul(Vxz,crab(:,14))-matmul(Vyz,crab(:,15))
					else 
						crab=0.0d0
						call Jor%get_adorned_cubes(i,j,k)
						call f4_index(i1,i2,j1,j2,i,j)
						call f2d1r(crab(:,1),x(:,i-1:i,j,k)) ! A1p
						call f2d1f(crab(:,2),x(:,i:i+1,j,k)) ! A1m
						call f4d1(crab(:,3),x(:,i1:i2,j,k),i,j,k,1) ! A2
						call f2d1r(crab(:,4),x(:,i,j-1:j,k)) ! B1p
						call f2d1f(crab(:,5),x(:,i,j:j+1,k)) ! B1m
						call f4d1(crab(:,6),x(:,i,j1:j2,k),i,j,k,2) ! B2
						call f3d1r(crab(:,7),x(:,i,j,k-2:k)) ! C1p
						call f3d1f(crab(:,8),x(:,i,j,k:k+2)) ! C1m
						call f5d1(crab(:,9),x(:,i,j,k-2:k+2)) ! C2
						call f4d11(crab(:,10),x(:,i-1:i+1,j,k)) ! Vxx
						call f4d11(crab(:,11),x(:,i,j-1:j+1,k)) ! Vyy
						call f5d11(crab(:,12),x(:,i,j,k-2:k+2)) ! Vzz
						call f4d12(crab(:,13),x(:,i1:i2,j1:j2,k),i,j,k,12) ! Vxy
						call f45d12(crab(:,14),x(:,i1:i2,j,k-2:k+2),i,j,k,1) ! Vxz
						call f45d12(crab(:,15),x(:,i,j1:j2,k-2:k+2),i,j,k,2) ! Vyz
						F(:,i,j,k) = matmul(A_p,crab(:,1))+matmul(A_m,crab(:,2))+matmul(A_v,crab(:,3)) &
									+matmul(B_p,crab(:,4))+matmul(B_m,crab(:,5))+matmul(B_v,crab(:,6)) &
									+matmul(C_p,crab(:,7))+matmul(C_m,crab(:,8))+matmul(C_v,crab(:,9)) &
									+matmul(D,x(:,i,j,k)) &
									-matmul(Vxx,crab(:,10))-matmul(Vyy,crab(:,11))-matmul(Vzz,crab(:,12)) &
									-matmul(Vxy,crab(:,13))-matmul(Vxz,crab(:,14))-matmul(Vyz,crab(:,15))
					endif 
				enddo
			enddo
		enddo
		end associate
		call DMDAVecRestoreArrayReadF90(meshDA,Clams,X_r,ierr)
		call DMDAVecRestoreArrayF90(meshDA,F,F_r,ierr)
	end subroutine Mymult

	subroutine shark_coming(comm)
		implicit none
		PetscInt,intent(in) :: comm
		PetscErrorCode :: ierr 
		call we_hear_a_sound(comm)
		call DMCreateMatrix(meshDA, Shark, ierr)
		call MatZeroEntries(Shark,ierr)
		call shark_growing_up()
		call whale_say_hi(comm)
	end subroutine shark_coming

	subroutine we_hear_a_sound(comm)
		implicit none
		PetscInt,intent(in) :: comm
		PetscErrorCode :: ierr
		call PetscPrintf(comm," -----------------------------------\n",ierr)
		call PetscPrintf(comm,"             开始生成矩阵             \n",ierr)
		call PetscPrintf(comm," -----------------------------------\n",ierr)
	end subroutine we_hear_a_sound
 
	subroutine shark_growing_up()
		implicit none
		PetscErrorCode :: ierr
		integer :: i,j,k
		select case (differential_scheme)
		case(0)
			do k=ks,ke
				do j=js,je
					do i=is,ie
						if(i==0 .or. i==(in-1) .or. j==0 .or. j==(jn-1))then
							call shark_eat_shrimps(i,j,k)
						else
							call shark_eat_tuna(i,j,k)
						endif
					enddo
				enddo
			enddo 
		case(1)
			do k=ks,ke
				do j=js,je
					do i=is,ie
						if(i==0 .or. i==(in-1) .or. j==0 .or. j==(jn-1))then
							call shark_eat_shrimps(i,j,k)
						else
							call shark_eat_sardine(i,j,k)
						endif
					enddo
				enddo
			enddo
		end select
		call MatAssemblyBegin(Shark,MAT_FINAL_ASSEMBLY,ierr)
		call MatAssemblyEnd(Shark,MAT_FINAL_ASSEMBLY,ierr)
	end subroutine shark_growing_up

	subroutine shark_eat_shrimps(i,j,k)
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
			call MatSetValuesBlockedStencil(Shark, 1, idxm, 1, idxn, box, INSERT_VALUES, ierr)
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
				call MatSetValuesBlockedStencil(Shark, 1, idxm, 1, idxn, trans, INSERT_VALUES, ierr)
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
				call MatSetValuesBlockedStencil(Shark, 1, idxm, 1, idxn, box, INSERT_VALUES, ierr)
			enddo
		endif
		end associate
	end subroutine shark_eat_shrimps

	subroutine shark_eat_sardine(i,j,k)
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
					call MatSetValuesBlockedStencil(Shark, 1, idxm, 1, idxn, trans, INSERT_VALUES, ierr)
				end do
			end do
		end do
		end associate
	end subroutine shark_eat_sardine

	subroutine shark_eat_tuna(i,j,k)
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
		  G => Jor%G,     D => Jor%D,    &
		  A => Jor%A,     B => Jor%B,     C => Jor%C,   &
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
						delta_j(lj)*delta_k(lk)*A*coef_c1(li,ic_index)+ &
						delta_i(li)*delta_k(lk)*B*coef_c1(lj,jc_index)+ &
						delta_i(li)*delta_j(lj)*C*coef_c1(lk,kc_index)+ &
						delta_j(lj)*delta_k(lk)*Vxx*coef_d1(li,ic_index)- &
						delta_i(li)*delta_k(lk)*Vyy*coef_d1(lj,jc_index)- &
						delta_i(li)*delta_j(lj)*Vzz*coef_d1(lk,kc_index)- &
						delta_k(lk)*Vxy*coef_c1(li,ic_index)*coef_c1(lj,jc_index)- &
						delta_j(lj)*Vxz*coef_c1(li,ic_index)*coef_c1(lk,kc_index)- &
						delta_i(li)*Vyz*coef_c1(lj,jc_index)*coef_c1(lk,kc_index)
					trans=transpose(box)
					call MatSetValuesBlockedStencil(Shark, 1, idxm, 1, idxn, trans, INSERT_VALUES, ierr)
				end do
			end do
		end do
		end associate
	end subroutine shark_eat_tuna

	subroutine whale_say_hi(comm)
		implicit none
		PetscInt,intent(in) :: comm
		PetscErrorCode :: ierr
		call PetscPrintf(comm," -----------------------------------\n",ierr)
		call PetscPrintf(comm,"             矩阵生成结束             \n",ierr)
		call PetscPrintf(comm," -----------------------------------\n",ierr)
	end subroutine whale_say_hi

	subroutine dolphin_say_hi(comm)
		implicit none
		PetscInt,intent(in) :: comm
		PetscErrorCode :: ierr
		call PetscPrintf(comm," -----------------------------------\n",ierr)
		call PetscPrintf(comm,"             免矩阵生效中             \n",ierr)
		call PetscPrintf(comm," -----------------------------------\n",ierr)
	end subroutine dolphin_say_hi
end module mod_forming
