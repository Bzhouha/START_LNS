module matrix_used_as_cofficient 
! -----------------------------------------------------------
!
!   这个模块生成系数小矩阵G、A、A_p、A_m、A_v、B、...、D、Vxx、Vxy等。
!
! 			(A_c:对流部分; A_v:粘性部分; A_p:正通量; A_m:负通量)
!
!       Type :: lns_OP_point_type 该点的系数小矩阵类
!
!           1.call get_unadorned_cubes(i,j,k) 获得原本的系数小矩阵
!
!               call split(A,G,Aplus,Aminus) 矢通量分裂函数
!
!           2.call get_adorned_cubes(i,j,k) 获得坐标变换后的系数小矩阵
!
! 			3.call colored_cubes(i,j,k) 获得对应方程的系数小矩阵
!
! 				1).call teal_cubes(i,j,k) 水鸭色的小块儿，2D-HLNS对应的系数小矩阵disturbance=f(x,y)*e^i(bz-wt)
!
! 				2).call mint_cubes(i,j,k) 薄荷色的小块儿，2D-HLNS对应的系数小矩阵,disturbance=f(x,y)*e^i(ax+bz-wt)
!
! 				3).call skyblue_cubes(i,j,k) 天空蓝色的小块儿，3D-HLNS对应的系数小矩阵,disturbance=f(x,y,z)*e^i(-wt)
!
! 				4).call lilac_cubes(i,j,k) 浅紫色的小块儿，3D-HLNS对应的系数小矩阵,disturbance=f(x,y,z)*e^i(ax-wt)
!
! -----------------------------------------------------------
	use penf, only: R_P
	implicit none
	private
	type,public :: lns_OP_point_type
		complex(R_P), dimension(5, 5) :: G=0.0d0
		complex(R_P), dimension(5, 5) :: A=0.0d0 
		complex(R_P), dimension(5, 5) :: A_p=0.0d0 
		complex(R_P), dimension(5, 5) :: A_m=0.0d0 
		complex(R_P), dimension(5, 5) :: A_v=0.0d0  
		complex(R_P), dimension(5, 5) :: B=0.0d0  
		complex(R_P), dimension(5, 5) :: B_p=0.0d0 
		complex(R_P), dimension(5, 5) :: B_m=0.0d0 
		complex(R_P), dimension(5, 5) :: B_v=0.0d0 
		complex(R_P), dimension(5, 5) :: C=0.0d0  
		complex(R_P), dimension(5, 5) :: C_p=0.0d0 
		complex(R_P), dimension(5, 5) :: C_m=0.0d0 
		complex(R_P), dimension(5, 5) :: C_v=0.0d0  
		complex(R_P), dimension(5, 5) :: D=0.0d0 
		complex(R_P), dimension(5, 5) :: Vxx=0.0d0 
		complex(R_P), dimension(5, 5) :: Vyy=0.0d0 
		complex(R_P), dimension(5, 5) :: Vzz=0.0d0 
		complex(R_P), dimension(5, 5) :: Vxy=0.0d0 
		complex(R_P), dimension(5, 5) :: Vxz=0.0d0 
		complex(R_P), dimension(5, 5) :: Vyz=0.0d0 
		Contains
		  procedure::get_unadorned_cubes,get_adorned_cubes
		  procedure::colored_cubes,teal_cubes,mint_cubes,skyblue_cubes,lilac_cubes
	end type lns_OP_point_type
	complex(R_P),parameter :: Li = cmplx(0.0d0,1.0d0,R_P)
	contains
	subroutine get_unadorned_cubes(this,i,j,k)
		use bf_point_org
		use global_parameters
		implicit none
		class(lns_OP_point_type),intent(inout) :: this 
		integer,intent(in) :: i,j,k
		real(R_P) :: Pe, gf, g1, g2, cm
		real(R_P) :: n_Miu, n_MiuT, n_MiuTT, n_Miux, n_Miuy, n_Miuz
		real(R_P) :: d1d3=1.0d0/3.0d0, d2d3=2.0d0/3.0d0, d4d3=4.0d0/3.0d0
		real(R_P) :: G(5, 5), A(5, 5), B(5, 5), C(5, 5), D(5, 5)
		real(R_P) :: A_p(5, 5), A_m(5, 5), A_c(5, 5), A_v(5, 5)
		real(R_P) :: B_p(5, 5), B_m(5, 5), B_c(5, 5), B_v(5, 5)
		real(R_P) :: C_p(5, 5), C_m(5, 5), C_c(5, 5), C_v(5, 5)
		real(R_P) :: Vxx(5, 5), Vyy(5, 5), Vzz(5, 5), Vxy(5, 5), Vxz(5, 5), Vyz(5, 5)
		real(R_P),parameter :: C1=110.4D0
		! 初始化矩阵
		G=0.0d0;A=0.0d0;B=0.0d0;C=0.0d0;D=0.0d0
		A_p=0.0d0;A_m=0.0d0;A_c=0.0d0;A_v=0.0d0
		B_p=0.0d0;B_m=0.0d0;B_c=0.0d0;B_v=0.0d0
		C_p=0.0d0;C_m=0.0d0;C_c=0.0d0;C_v=0.0d0
		Vxx=0.0d0;Vyy=0.0d0;Vzz=0.0d0
		Vxy=0.0d0;Vxz=0.0d0;Vyz=0.0d0
		! 设定参数
		Pe = 1.0d0/(GAMMA*Ma*Ma)
		gf = 1.0d0/GAMMA 
		g1 = (1.0d0-GAMMA)/GAMMA
		g2 = 1.0d0/((GAMMA-1.0d0)*Ma*Ma)
		associate ( &
			rho  => bf(i,j,k)%BF%rho, &
			U    => bf(i,j,k)%BF%x, &
			V    => bf(i,j,k)%BF%y, &
			W    => bf(i,j,k)%BF%z, &
			T    => bf(i,j,k)%BF%T, &
		
			rhox => bf(i,j,k)%BFDx%rho, &
			Ux   => bf(i,j,k)%BFDx%x, &
			Vx   => bf(i,j,k)%BFDx%y, &
			Wx   => bf(i,j,k)%BFDx%z, &
			Tx   => bf(i,j,k)%BFDx%T, &
		
			rhoy => bf(i,j,k)%BFDy%rho, &
			Uy   => bf(i,j,k)%BFDy%x, &
			Vy   => bf(i,j,k)%BFDy%y, &
			Wy   => bf(i,j,k)%BFDy%z, &
			Ty   => bf(i,j,k)%BFDy%T, &
		
			rhoz => bf(i,j,k)%BFDz%rho, &
			Uz   => bf(i,j,k)%BFDz%x, &
			Vz   => bf(i,j,k)%BFDz%y, &
			Wz   => bf(i,j,k)%BFDz%z, &
			Tz   => bf(i,j,k)%BFDz%T, &
		
			Uxx  => bf(i,j,k)%BFDxx%x, &
			nVxx => bf(i,j,k)%BFDxx%y, & ! 加前缀n与系数矩阵区分
			Wxx  => bf(i,j,k)%BFDxx%z, &
			Txx  => bf(i,j,k)%BFDxx%T, &
		
			Uyy  => bf(i,j,k)%BFDyy%x, &
			nVyy => bf(i,j,k)%BFDyy%y, & 
			Wyy  => bf(i,j,k)%BFDyy%z, &
			Tyy  => bf(i,j,k)%BFDyy%T, &
		
			Uzz  => bf(i,j,k)%BFDzz%x, &
			nVzz => bf(i,j,k)%BFDzz%y, & 
			Wzz  => bf(i,j,k)%BFDzz%z, &
			Tzz  => bf(i,j,k)%BFDzz%T, &
		
			Uxy  => bf(i,j,k)%BFDxy%x, &
			nVxy => bf(i,j,k)%BFDxy%y, & 
			Wxy  => bf(i,j,k)%BFDxy%z, &
		
			Uxz  => bf(i,j,k)%BFDxz%x, &
			nVxz => bf(i,j,k)%BFDxz%y, & 
			Wxz  => bf(i,j,k)%BFDxz%z, &
		
			Uyz  => bf(i,j,k)%BFDyz%x, &
			nVyz => bf(i,j,k)%BFDyz%y, & 
			Wyz  => bf(i,j,k)%BFDyz%z  )
			! Surthland公式相关
			cm=C1/Te
			n_Miu = T*sqrt(T)*(1.0d0+cm)/(T+cm)
			n_MiuT = n_Miu*(1.5d0/T-1.0d0/(T+cm))
			n_MiuTT = n_MiuT*(1.5d0/T-1.0d0/(T+cm))-n_Miu*(1.5d0/T**2-1.0d0/(T+cm)**2)
			n_Miux = n_MiuT*Tx
			n_Miuy = n_MiuT*Ty
			n_Miuz = n_MiuT*Tz
			! 下面开始制造矩阵
			G(1, 1) = 1.0d0
			G(2, 2) = rho
			G(3, 3) = rho
			G(4, 4) = rho
			G(5, 1) = -1.0d0*Pe*T
			G(5, 5) = rho*g2-Pe*rho
			! G(5, 1) = g1*T 
			! G(5, 5) = rho/GAMMA

			A(1, 1) = U
			A(1, 2) = rho
			A(2, 1) = Pe*T 
			A(2, 2) = rho*U - d4d3*n_Miux/Re
			A(2, 3) = -1.0d0*n_Miuy/Re
			A(2, 4) = -1.0d0*n_Miuz/Re
			A(2, 5) = Pe*rho - n_MiuT/Re*(d4d3*Ux-d2d3*Vy-d2d3*Wz)
			A(3, 2) = d2d3*n_Miuy/Re
			A(3, 3) = rho*U - n_Miux/Re
			A(3, 5) = -1.0d0*n_MiuT*(Uy+Vx)/Re
			A(4, 2) = d2d3*n_Miuz/Re
			A(4, 4) = rho*U - n_Miux/Re
			A(4, 5) = -1.0d0*n_MiuT*(Wx+Uz)/Re
			A(5, 1) = -1.0d0*Pe*U*T
			A(5, 2) = -2.0d0*n_Miu*(d4d3*Ux-d2d3*Vy-d2d3*Wz)/Re
			A(5, 3) = -2.0d0*n_Miu*(Uy+Vx)/Re
			A(5, 4) = -2.0d0*n_Miu*(Wx+Uz)/Re
			A(5, 5) = rho*U*g2-rho*U*Pe-2.0d0*Tx*n_MiuT/Re/Pr*g2
			! A(5, 1) = g1*T*U
			! A(5, 2) = -1.0d0*d4d3*g2*n_Miu*(2.0d0*Ux-Vy-Wz)/Re 
			! A(5, 3) = -2.0d0*g2*n_Miu*(Vx+Uy)/Re
			! A(5, 4) = -2.0d0*g2*n_Miu*(Wx+Uz)/Re
			! A(5, 5) = gf*rho*U - 2.0d0*n_Miux/Re/Pr

			A_c(1, 1) = U
			A_c(1, 2) = rho
			A_c(2, 1) = Pe*T 
			A_c(2, 2) = rho*U
			A_c(2, 5) = Pe*rho
			A_c(3, 3) = rho*U
			A_c(4, 4) = rho*U
			A_c(5, 1) = -1.0d0*Pe*U*T
			A_c(5, 5) = rho*U*g2-rho*U*Pe

			A_v(2, 2) = -1.0d0*d4d3*n_Miux/Re
			A_v(2, 3) = -1.0d0*n_Miuy/Re
			A_v(2, 4) = -1.0d0*n_Miuz/Re
			A_v(2, 5) = -1.0d0*n_MiuT/Re*(d4d3*Ux-d2d3*Vy-d2d3*Wz)
			A_v(3, 2) = d2d3*n_Miuy/Re
			A_v(3, 3) = -1.0d0*n_Miux/Re
			A_v(3, 5) = -1.0d0*n_MiuT*(Uy+Vx)/Re
			A_v(4, 2) = d2d3*n_Miuz/Re
			A_v(4, 4) = -1.0d0*n_Miux/Re
			A_v(4, 5) = -1.0d0*n_MiuT*(Wx+Uz)/Re
			A_v(5, 2) = -2.0d0*n_Miu*(d4d3*Ux-d2d3*Vy-d2d3*Wz)/Re
			A_v(5, 3) = -2.0d0*n_Miu*(Uy+Vx)/Re
			A_v(5, 4) = -2.0d0*n_Miu*(Wx+Uz)/Re
			A_v(5, 5) = -2.0d0*Tx*n_MiuT/Re/Pr*g2
			
			B(1, 1) = V
			B(1, 3) = rho
			B(2, 2) = rho*V - n_Miuy/Re
			B(2, 3) = d2d3*n_Miux/Re
			B(2, 5) = -1.0d0*n_MiuT*(Uy+Vx)/Re
			B(3, 1) = Pe*T
			B(3, 2) = -1.0d0*n_Miux/Re
			B(3, 3) = rho*V - d4d3*n_Miuy/Re
			B(3, 4) = -1.0d0*n_Miuz/Re
			B(3, 5) = Pe*rho - n_MiuT*(-d2d3*Ux+d4d3*Vy-d2d3*Wz)/Re
			B(4, 3) = d2d3*n_Miuz/Re
			B(4, 4) = rho*V - n_Miuy/Re
			B(4, 5) = -1.0d0*n_MiuT*(Vz+Wy)/Re
			B(5, 1) = -1.0d0*Pe*V*T
			B(5, 2) = -2.0d0*n_Miu*(Uy+Vx)/Re
			B(5, 3) = -2.0d0*n_Miu*(-d2d3*Ux+d4d3*Vy-d2d3*Wz)/Re
			B(5, 4) = -2.0d0*n_Miu*(Vz+Wy)/Re
			B(5, 5) = rho*V*g2-rho*V*Pe-2.0d0*Ty*n_MiuT*g2/Re/Pr
			! B(5, 1) = g1*T*V
			! B(5, 2) = -2.0d0*g2*n_Miu*(Uy+Vx)/Re
			! B(5, 3) = -1.0d0*d4d3*g2*n_Miu*(2.0d0*Vy-Ux-Wz)/Re
			! B(5, 4) = -2.0d0*g2*n_Miu*(Wy+Vz)/Re
			! B(5, 5) = gf*rho*V - 2.0d0*n_Miuy/Re/Pr 

			B_c(1, 1) = V
			B_c(1, 3) = rho
			B_c(2, 2) = rho*V
			B_c(3, 1) = Pe*T
			B_c(3, 3) = rho*V
			B_c(3, 5) = Pe*rho
			B_c(4, 4) = rho*V
			B_c(5, 1) = -1.0d0*Pe*V*T
			B_c(5, 5) = rho*V*g2-rho*V*Pe

			B_v(2, 2) = -1.0d0*n_Miuy/Re
			B_v(2, 3) = d2d3*n_Miux/Re
			B_v(2, 5) = -1.0d0*n_MiuT*(Uy+Vx)/Re
			B_v(3, 2) = -1.0d0*n_Miux/Re
			B_v(3, 3) = -1.0d0*d4d3*n_Miuy/Re
			B_v(3, 4) = -1.0d0*n_Miuz/Re
			B_v(3, 5) = -1.0d0*n_MiuT*(-d2d3*Ux+d4d3*Vy-d2d3*Wz)/Re
			B_v(4, 3) = d2d3*n_Miuz/Re
			B_v(4, 4) = -1.0d0*n_Miuy/Re
			B_v(4, 5) = -1.0d0*n_MiuT*(Vz+Wy)/Re
			B_v(5, 2) = -2.0d0*n_Miu*(Uy+Vx)/Re
			B_v(5, 3) = -2.0d0*n_Miu*(-d2d3*Ux+d4d3*Vy-d2d3*Wz)/Re
			B_v(5, 4) = -2.0d0*n_Miu*(Vz+Wy)/Re
			B_v(5, 5) = -2.0d0*Ty*n_MiuT*g2/Re/Pr 
			
			C(1, 1) = W
			C(1, 4) = rho
			C(2, 2) = rho*W - n_Miuz/Re
			C(2, 4) = d2d3*n_Miux/Re
			C(2, 5) = -1.0d0*n_MiuT*(Wx+Uz)/Re
			C(3, 3) = rho*W - n_Miuz/Re
			C(3, 4) = d2d3*n_Miuy/Re
			C(3, 5) = -1.0d0*n_MiuT*(Vz+Wy)/Re
			C(4, 1) = Pe*T
			C(4, 2) = -1.0d0*n_Miux/Re
			C(4, 3) = -1.0d0*n_Miuy/Re
			C(4, 4) = rho*W - d4d3*n_Miuz/Re
			C(4, 5) = Pe*rho - n_MiuT*(-d2d3*Ux-d2d3*Vy+d4d3*Wz)/Re
			C(5, 1) = -1.0d0*Pe*W*T
			C(5, 2) = -2.0d0*n_Miu*(Wx+Uz)/Re
			C(5, 3) = -2.0d0*n_Miu*(Vz+Wy)/Re
			C(5, 4) = -2.0d0*n_Miu*(-d2d3*Ux-d2d3*Vy+d4d3*Wz)/Re
			C(5, 5) = rho*W*g2-rho*W*Pe-2.0d0*Tz*n_MiuT*g2/Re/Pr
			! C(5, 1) = g1*T*W
			! C(5, 2) = -2.0d0*g2*n_Miu*(Uz+Wx)/Re
			! C(5, 3) = -2.0d0*g2*n_Miu*(Vz+Wy)/Re
			! C(5, 4) = -1.0d0*d4d3*g2*n_Miu*(2.0d0*Wz-Ux-Vy)/Re
			! C(5, 5) = gf*rho*W - 2.0d0*n_Miuz/Re/Pr

			C_c(1, 1) = W
			C_c(1, 4) = rho
			C_c(2, 2) = rho*W
			C_c(3, 3) = rho*W
			C_c(4, 1) = Pe*T
			C_c(4, 4) = rho*W
			C_c(4, 5) = Pe*rho
			C_c(5, 1) = -1.0d0*Pe*W*T
			C_c(5, 5) = rho*W*g2-rho*W*Pe

			C_v(2, 2) = -1.0d0*n_Miuz/Re
			C_v(2, 4) = d2d3*n_Miux/Re
			C_v(2, 5) = -1.0d0*n_MiuT*(Wx+Uz)/Re
			C_v(3, 3) = -1.0d0*n_Miuz/Re
			C_v(3, 4) = d2d3*n_Miuy/Re
			C_v(3, 5) = -1.0d0*n_MiuT*(Vz+Wy)/Re
			C_v(4, 2) = -1.0d0*n_Miux/Re
			C_v(4, 3) = -1.0d0*n_Miuy/Re
			C_v(4, 4) = -1.0d0*d4d3*n_Miuz/Re
			C_v(4, 5) = -1.0d0*n_MiuT*(-d2d3*Ux-d2d3*Vy+d4d3*Wz)/Re
			C_v(5, 2) = -2.0d0*n_Miu*(Wx+Uz)/Re
			C_v(5, 3) = -2.0d0*n_Miu*(Vz+Wy)/Re
			C_v(5, 4) = -2.0d0*n_Miu*(-d2d3*Ux-d2d3*Vy+d4d3*Wz)/Re
			C_v(5, 5) = -2.0d0*Tz*n_MiuT*g2/Re/Pr
			
			D(1, 1) = Ux+Vy+Wz
			D(1, 2) = rhox
			D(1, 3) = rhoy 
			D(1, 4) = rhoz
			D(2, 1) = U*Ux+V*Uy+W*Uz+Pe*Tx
			D(2, 2) = rho*Ux
			D(2, 3) = rho*Uy
			D(2, 4) = rho*Uz
			D(2, 5) = Pe*rhox - ( n_MiuTT*Tx*d2d3*(2.0d0*Ux-Vy-Wz) + n_MiuTT*Ty*(Uy+Vx) &
			+ n_MiuTT*Tz*(Wx+Uz) + n_MiuT*(d4d3*Uxx+Uyy+Uzz+d1d3*nVxy+d1d3*Wxz) )/Re
			D(3, 1) = U*Vx+V*Vy+W*Vz+Pe*Ty
			D(3, 2) = rho*Vx
			D(3, 3) = rho*Vy
			D(3, 4) = rho*Vz
			D(3, 5) = Pe*rhoy - ( n_MiuTT*Tx*(Uy+Vx) + n_MiuTT*Ty*d2d3*(-1.0d0*Ux+2.0d0*Vy-Wz) &
			+ n_MiuTT*Tz*(Vz+Wy) + n_MiuT*(nVxx+d4d3*nVyy+nVzz+d1d3*Uxy+d1d3*Wyz) )/Re
			D(4, 1) = U*Wx+V*Wy+W*Wz+Pe*Tz
			D(4, 2) = rho*Wx
			D(4, 3) = rho*Wy
			D(4, 4) = rho*Wz
			D(4, 5) = Pe*rhoz - ( n_MiuTT*Tx*(Wx+Uz) + n_MiuTT*Ty*(Vz+Wy) + &
			n_MiuTT*Tz*d2d3*(-1.0d0*Ux-Vy+2.0d0*Wz) + n_MiuT*(Wxx+Wyy+d4d3*Wzz+d1d3*Uxz+d1d3*nVyz) )/Re
			D(5, 1) = (1.0d0*g2-Pe)*(U*Tx+V*Ty+W*Tz)
			D(5, 2) = rho*Tx*g2-Pe*(rho*Tx+T*rhox)
			D(5, 3) = rho*Ty*g2-Pe*(rho*Ty+T*rhoy)
			D(5, 4) = rho*Tz*g2-Pe*(rho*Tz+T*rhoz)
			D(5, 5) = -1.0d0*Pe*(U*rhox+V*rhoy+W*rhoz)-(Txx+Tyy+Tzz)*n_MiuT*g2/Re/Pr&
			-(Tx*Tx+Ty*Ty+Tz*Tz)*n_MiuTT*g2/Re/Pr&
			-n_MiuT*(Uy*Uy+Vx*Vx+2.0d0*Uy*Vx+Uz*Uz+Wx*Wx+2.0d0*Uz*Wx+Vz*Vz+Wy*Wy+2.0d0*Vz*Wy)/Re&
			-n_MiuT*d4d3*(Ux*Ux+Vy*Vy+Wz*Wz-Ux*Wz-Ux*Vy-Vy*Wz)/Re
			! D(5, 1) = gf*U*Tx+gf*V*Ty+gf*W*Tz
			! D(5, 2) = gf*rho*Tx+g1*T*rhox
			! D(5, 3) = gf*rho*Ty+g1*T*rhoy
			! D(5, 4) = gf*rho*Tz+g1*T*rhoz
			! D(5, 5) = g1*U*rhox+g1*V*rhoy+g1*W*rhoz - ( n_MiuTT*(Tx*Tx+Ty*Ty+Tz*Tz) &
			! + n_MiuT*(Txx+Tyy+Tzz) )/Re/Pr - g2*( d4d3*n_MiuT*(Ux*Ux+Vy*Vy+Wz*Wz-Ux*Vy-Ux*Wz-Vy*Wz) &
			! + n_MiuT*(Uy*Uy+Vx*Vx+2.0d0*Uy*Vx+Uz*Uz+Wx*Wx+2.0d0*Uz*Wx+Vz*Vz+Wy*Wy+2.0d0*Vz*Wy) )/Re
			
			Vxx(2, 2) = d4d3*n_Miu/Re
			Vxx(3, 3) = n_Miu/Re
			Vxx(4, 4) = n_Miu/Re
			Vxx(5, 5) = n_Miu*g2/Re/Pr
			! Vxx(5, 5) = n_Miu/Re/Pr
			
			Vyy(2, 2) = n_Miu/Re
			Vyy(3, 3) = d4d3*n_Miu/Re
			Vyy(4, 4) = n_Miu/Re
			Vyy(5, 5) = n_Miu*g2/Re/Pr
			! Vyy(5, 5) = n_Miu/Re/Pr 
			
			Vzz(2, 2) = n_Miu/Re
			Vzz(3, 3) = n_Miu/Re
			Vzz(4, 4) = d4d3*n_Miu/Re
			Vzz(5, 5) = n_Miu*g2/Re/Pr
			! Vzz(5, 5) = n_Miu/Re/Pr
			
			Vxy(2, 3) = d1d3*n_Miu/Re
			Vxy(3, 2) = d1d3*n_Miu/Re
		
			Vxz(2, 4) = d1d3*n_Miu/Re
			Vxz(4, 2) = d1d3*n_Miu/Re
		
			Vyz(3, 4) = d1d3*n_Miu/Re
			Vyz(4, 3) = d1d3*n_Miu/Re
		end associate

		call split(A_c,G,A_p,A_m)
		call split(B_c,G,B_p,B_m)
		call split(C_c,G,C_p,C_m)

		this%G=G;     this%D=D
		this%A=A;     this%B=B;     this%C=C
		this%A_p=A_p; this%A_m=A_m; this%A_v=A_v
		this%B_p=B_p; this%B_m=B_m; this%B_v=B_v
		this%C_p=C_p; this%C_m=C_m; this%C_v=C_v
		this%Vxx=Vxx; this%Vyy=Vyy; this%Vzz=Vzz
		this%Vxy=Vxy; this%Vxz=Vxz; this%Vyz=Vyz
	end subroutine get_unadorned_cubes

	subroutine get_adorned_cubes(this,i,j,k)
		use global_parameters
		implicit none
		class(lns_OP_point_type),intent(inout) :: this
		type(lns_OP_point_type) :: Jor
		integer,intent(in) :: i,j,k
		complex(R_P) :: G(5, 5), D(5, 5)
		complex(R_P) :: A(5,5), B(5,5), C(5,5)
		complex(R_P) :: A_p(5, 5), A_m(5, 5), A_v(5, 5)
		complex(R_P) :: B_p(5, 5), B_m(5, 5), B_v(5, 5)
		complex(R_P) :: C_p(5, 5), C_m(5, 5), C_v(5, 5)
		complex(R_P) :: Vxx(5, 5), Vyy(5, 5), Vzz(5, 5), Vxy(5, 5), Vxz(5, 5), Vyz(5, 5)
		! 初始化矩阵
		G=0.0d0;D=0.0d0
		A=0.0d0;B=0.0d0;C=0.0d0
		A_p=0.0d0;A_m=0.0d0;A_v=0.0d0
		B_p=0.0d0;B_m=0.0d0;B_v=0.0d0
		C_p=0.0d0;C_m=0.0d0;C_v=0.0d0
		Vxx=0.0d0;Vyy=0.0d0;Vzz=0.0d0
		Vxy=0.0d0;Vxz=0.0d0;Vyz=0.0d0
		call Jor%get_unadorned_cubes(i,j,k)
		associate( &
			xi_x => xi_x(i,j,k), &
			xi_y => xi_y(i,j,k), &
			xi_z => xi_z(i,j,k), &
			eta_x => eta_x(i,j,k), &
			eta_y => eta_y(i,j,k), &
			eta_z => eta_z(i,j,k), &
			phi_x => phi_x(i,j,k), &
			phi_y => phi_y(i,j,k), &
			phi_z => phi_z(i,j,k), &
			xi_xx => xi_xx(i,j,k), &
			xi_yy => xi_yy(i,j,k), &
			xi_zz => xi_zz(i,j,k), &
			eta_xx => eta_xx(i,j,k), &
			eta_yy => eta_yy(i,j,k), &
			eta_zz => eta_zz(i,j,k), &
			phi_xx => phi_xx(i,j,k), &
			phi_yy => phi_yy(i,j,k), &
			phi_zz => phi_zz(i,j,k), &
			xi_xy => xi_xy(i,j,k), &
			xi_xz => xi_xz(i,j,k), &
			xi_yz => xi_yz(i,j,k), &
			eta_xy => eta_xy(i,j,k), &
			eta_xz => eta_xz(i,j,k), &
			eta_yz => eta_yz(i,j,k), &
			phi_xy => phi_xy(i,j,k), &
			phi_xz => phi_xz(i,j,k), &
			phi_yz => phi_yz(i,j,k) )
			G   = Jor%G
			A   = xi_x*Jor%A+xi_y*Jor%B+xi_z*Jor%C-xi_xx*Jor%Vxx-xi_yy*Jor%Vyy-xi_zz*Jor%Vzz &
			-xi_xy*Jor%Vxy-xi_xz*Jor%Vxz-xi_yz*Jor%Vyz
			A_p = xi_x*Jor%A_p+xi_y*Jor%B_p+xi_z*Jor%C_p
			A_m = xi_x*Jor%A_m+xi_y*Jor%B_m+xi_z*Jor%C_m
			A_v = xi_x*Jor%A_v+xi_y*Jor%B_v+xi_z*Jor%C_v-xi_xx*Jor%Vxx-xi_yy*Jor%Vyy-xi_zz*Jor%Vzz &
			-xi_xy*Jor%Vxy-xi_xz*Jor%Vxz-xi_yz*Jor%Vyz
			B   = eta_x*Jor%A+eta_y*Jor%B+eta_z*Jor%C-eta_xx*Jor%Vxx-eta_yy*Jor%Vyy-eta_zz*Jor%Vzz &
			-eta_xy*Jor%Vxy-eta_xz*Jor%Vxz-eta_yz*Jor%Vyz
			B_p = eta_x*Jor%A_p+eta_y*Jor%B_p+eta_z*Jor%C_p
			B_m = eta_x*Jor%A_m+eta_y*Jor%B_m+eta_z*Jor%C_m
			B_v = eta_x*Jor%A_v+eta_y*Jor%B_v+eta_z*Jor%C_v-eta_xx*Jor%Vxx-eta_yy*Jor%Vyy-eta_zz*Jor%Vzz &
			-eta_xy*Jor%Vxy-eta_xz*Jor%Vxz-eta_yz*Jor%Vyz
			C   = phi_x*Jor%A+phi_y*Jor%B+phi_z*Jor%C-phi_xx*Jor%Vxx-phi_yy*Jor%Vyy-phi_zz*Jor%Vzz &
			-phi_xy*Jor%Vxy-phi_xz*Jor%Vxz-phi_yz*Jor%Vyz
			C_p = phi_x*Jor%A_p+phi_y*Jor%B_p+phi_z*Jor%C_p
			C_m = phi_x*Jor%A_m+phi_y*Jor%B_m+phi_z*Jor%C_m
			C_v = phi_x*Jor%A_v+phi_y*Jor%B_v+phi_z*Jor%C_v-phi_xx*Jor%Vxx-phi_yy*Jor%Vyy-phi_zz*Jor%Vzz &
			-phi_xy*Jor%Vxy-phi_xz*Jor%Vxz-phi_yz*Jor%Vyz
			D   = Jor%D
			Vxx = xi_x*xi_x*Jor%Vxx+xi_y*xi_y*Jor%Vyy+xi_z*xi_z*Jor%Vzz+xi_x*xi_y*Jor%Vxy &
			+xi_x*xi_z*Jor%Vxz+xi_y*xi_z*Jor%Vyz
			Vyy = eta_x*eta_x*Jor%Vxx+eta_y*eta_y*Jor%Vyy+eta_z*eta_z*Jor%Vzz+eta_x*eta_y*Jor%Vxy &
			+eta_x*eta_z*Jor%Vxz+eta_y*eta_z*Jor%Vyz
			Vzz = phi_x*phi_x*Jor%Vxx+phi_y*phi_y*Jor%Vyy+phi_z*phi_z*Jor%Vzz+phi_x*phi_y*Jor%Vxy &
			+phi_x*phi_z*Jor%Vxz+phi_y*phi_z*Jor%Vyz
			Vxy = 2.0d0*xi_x*eta_x*Jor%Vxx+2.0d0*xi_y*eta_y*Jor%Vyy+2.0d0*xi_z*eta_z*Jor%Vzz &
			+(xi_x*eta_y+eta_x*xi_y)*Jor%Vxy+(xi_x*eta_z+eta_x*xi_z)*Jor%Vxz &
			+(xi_y*eta_z+eta_y*xi_z)*Jor%Vyz
			Vxz = 2.0d0*xi_x*phi_x*Jor%Vxx+2.0d0*xi_y*phi_y*Jor%Vyy+2.0d0*xi_z*phi_z*Jor%Vzz &
			+(xi_x*phi_y+phi_x*xi_y)*Jor%Vxy+(xi_x*phi_z+phi_x*xi_z)*Jor%Vxz &
			+(xi_y*phi_z+phi_y*xi_z)*Jor%Vyz
			Vyz = 2.0d0*eta_x*phi_x*Jor%Vxx+2.0d0*eta_y*phi_y*Jor%Vyy+2.0d0*eta_z*phi_z*Jor%Vzz &
			+(eta_x*phi_y+phi_x*eta_y)*Jor%Vxy+(eta_x*phi_z+phi_x*eta_z)*Jor%Vxz &
			+(eta_y*phi_z+phi_y*eta_z)*Jor%Vyz
		end associate

		this%G=G;	  this%D=D
		this%A=A;     this%B=B;     this%C=C
		this%A_p=A_p; this%A_m=A_m; this%A_v=A_v
		this%B_p=B_p; this%B_m=B_m; this%B_v=B_v
		this%C_p=C_p; this%C_m=C_m; this%C_v=C_v
		this%Vxx=Vxx; this%Vyy=Vyy; this%Vzz=Vzz 
		this%Vxy=Vxy; this%Vxz=Vxz; this%Vyz=Vyz
	end subroutine get_adorned_cubes

	subroutine colored_cubes(this,i,j,k)
		use global_parameters,only:lns_mode
		implicit none
		class(lns_OP_point_type),intent(inout) :: this
		integer,intent(in) :: i,j,k
		select case (lns_mode)
			case(0)
				call this%mint_cubes(i,j,k)
			case(1)
				call this%lilac_cubes(i,j,k)
		end select
	end subroutine colored_cubes

	subroutine teal_cubes(this,i,j,k)
		use global_parameters,only:Beta,Omega 
		implicit none
		class(lns_OP_point_type),intent(inout) :: this
		type(lns_OP_point_type) :: Jor
		integer,intent(in) :: i,j,k
		call Jor%get_adorned_cubes(i,j,k)
		this%G=Jor%G
		this%A=Jor%A;this%B=Jor%B;this%C=Jor%C
		! Notice that the cubes above has been merged into below ones. Do not touch those durning assembing.
		this%A_p=Jor%A_p; this%A_m=Jor%A_m
		this%A_v=Jor%A_v-Li*Beta*Jor%Vxz
		this%B_p=Jor%B_p; this%B_m=Jor%B_m
		this%B_v=Jor%B_v-Li*Beta*Jor%Vyz 
		this%C_p=0.0d0;   this%C_m=0.0d0;   this%C_v=0.0d0
		this%D=Jor%D-Li*Omega*Jor%G+Li*Beta*Jor%C+Beta*Beta*Jor%Vzz
		this%Vxx=Jor%Vxx; this%Vyy=Jor%Vyy; this%Vzz=0.0d0
		this%Vxy=Jor%Vxy; this%Vxz=0.0d0;   this%Vyz=0.0d0 
	end subroutine teal_cubes

	subroutine mint_cubes(this,i,j,k)
		use global_parameters,only:Alpha,Beta,Omega 
		implicit none 
		class(lns_OP_point_type),intent(inout) :: this
		type(lns_OP_point_type) :: Jor
		integer,intent(in) :: i,j,k
		call Jor%get_adorned_cubes(i,j,k)
		this%G=Jor%G
		this%A=Jor%A;this%B=Jor%B;this%C=Jor%C
		! Notice that the cubes above has been merged into below ones. Do not touch those durning assembing.
		this%A_p=Jor%A_p; this%A_m=Jor%A_m
		this%A_v=Jor%A_v-2.0d0*Li*Alpha*Jor%Vxx-Li*Beta*Jor%Vxz
		this%B_p=Jor%B_p; this%B_m=Jor%B_m
		this%B_v=Jor%B_v-Li*Alpha*Jor%Vxy-Li*Beta*Jor%Vyz 
		this%C_p=0.0d0;   this%C_m=0.0d0;   this%C_v=0.0d0
		this%D=Jor%D-Li*Omega*Jor%G+Li*Alpha*Jor%A &
		+Li*Beta*Jor%C+Alpha*Alpha*Jor%Vxx+Beta*Beta*Jor%Vzz+Alpha*Beta*Jor%Vxz
		this%Vxx=Jor%Vxx; this%Vyy=Jor%Vyy; this%Vzz=0.0d0
		this%Vxy=Jor%Vxy; this%Vxz=0.0d0;   this%Vyz=0.0d0 
	end subroutine mint_cubes

	subroutine skyblue_cubes(this,i,j,k)
		use global_parameters,only:Omega
		implicit none 
		class(lns_OP_point_type),intent(inout) :: this
		type(lns_OP_point_type) :: Jor
		integer,intent(in) :: i,j,k
		call Jor%get_adorned_cubes(i,j,k)
		this%G=Jor%G 
		this%A=Jor%A;this%B=Jor%B;this%C=Jor%C
		! Notice that the cubes above has been merged into below ones. Do not touch those durning assembing.
		this%A_p=Jor%A_p; this%A_m=Jor%A_m; this%A_v=Jor%A_v
		this%B_p=Jor%B_p; this%B_m=Jor%B_m; this%B_v=Jor%B_v
		this%C_p=Jor%C_p; this%C_m=Jor%C_m; this%C_v=Jor%C_v
		this%D=Jor%D-Li*Omega*Jor%G
		this%Vxx=Jor%Vxx; this%Vyy=Jor%Vyy; this%Vzz=Jor%Vzz 
		this%Vxy=Jor%Vxy; this%Vxz=Jor%Vxz; this%Vyz=Jor%Vyz
	end subroutine skyblue_cubes

	subroutine lilac_cubes(this,i,j,k)
		use global_parameters,only:Alpha,Omega
		implicit none 
		class(lns_OP_point_type),intent(inout) :: this
		type(lns_OP_point_type) :: Jor
		integer,intent(in) :: i,j,k
		call Jor%get_adorned_cubes(i,j,k)
		this%G=Jor%G 
		this%A=Jor%A;this%B=Jor%B;this%C=Jor%C
		! Notice that the cubes above has been merged into below ones. Do not touch those durning assembing.
		this%A_p=Jor%A_p; this%A_m=Jor%A_m; this%A_v=Jor%A_v-2*Li*Alpha*Jor%Vxx
		this%B_p=Jor%B_p; this%B_m=Jor%B_m; this%B_v=Jor%B_v-Li*Alpha*Jor%Vxy
		this%C_p=Jor%C_p; this%C_m=Jor%C_m; this%C_v=Jor%C_v-Li*Alpha*Jor%Vxz
		this%D=Jor%D-Li*Omega*Jor%G+Li*Alpha*Jor%A+Alpha*Alpha*Jor%Vxx
		this%Vxx=Jor%Vxx; this%Vyy=Jor%Vyy; this%Vzz=Jor%Vzz 
		this%Vxy=Jor%Vxy; this%Vxz=Jor%Vxz; this%Vyz=Jor%Vyz
	end subroutine lilac_cubes

	subroutine split(A,G,Aplus,Aminus)
		implicit none
		real(R_P), dimension(5, 5), intent(out) :: Aplus,Aminus
		real(R_P), dimension(5, 5), intent(in) :: A, G
		real(R_P) :: diag_plus(5, 5),diag_minus(5, 5)
		real(R_P), dimension(5, 5) :: At, Gt, Ab
		real(R_P) :: alfr(5), alfi(5), beta(5)
		real(R_P) :: vl(5, 5), vr(5, 5)
		complex(R_P) :: lambda_(5)
		integer :: ipiv(5), info
		real(R_P) :: work(100)
		complex(R_P) :: ZI=(0.0d0,1.0d0)
		integer :: i
		Aplus=0.0d0; Aminus=0.0d0; work=0.0d0
		alfr=0.0d0; alfi=0.0d0; beta=0.0d0
		vl=0.0d0; vr=0.0d0; lambda_=0.0d0
		diag_plus=0.0d0; diag_minus=0.0d0
		At=A; Gt=G
		call dggev('N', 'V', 5, At, 5, Gt, 5, alfr, alfi, &
		beta, vl, 5, vr, 5, work, 100, info)
		lambda_=(alfr+ZI*alfi)/beta
		do i=1, 5
			diag_plus(i, i)=max(0.0d0, lambda_(i)%re)
			diag_minus(i, i)=min(0.0d0, lambda_(i)%re)
		enddo
		vl=vr !! vl is inv(vr)
		call dgetrf(5, 5, vl, 5, ipiv, info)
		call dgetri(5, vl, 5, ipiv, work, 5, info)
		At=matmul(G, vr)
		At=matmul(At, diag_plus)
		At=matmul(At, vl)
		Ab=matmul(G, vr)
		Ab=matmul(Ab, diag_minus)
		Ab=matmul(Ab, vl)
		Aplus=At
		Aminus=Ab
	end subroutine split
end module matrix_used_as_cofficient
