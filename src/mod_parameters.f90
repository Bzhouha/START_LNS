! #include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_parameters
    use penf, only: R_P
    use mod_flowtype
    use petsc
    implicit none
    real(R_P), dimension(:, :, :), allocatable :: xi_xx,xi_yy,xi_zz,eta_xx,eta_yy,eta_zz,phi_xx,phi_yy,phi_zz
    real(R_P), dimension(:, :, :), allocatable :: xi_xy,xi_xz,xi_yz,eta_xy,eta_yz,eta_xz,phi_xy,phi_yz,phi_xz 
    real(R_P), dimension(:, :, :), allocatable :: xi_x,xi_y,xi_z,eta_x,eta_y,eta_z,phi_x,phi_y,phi_z
    complex(R_P), dimension(:, :, :), allocatable :: disturb
    real(R_P), dimension(:, :, :), allocatable :: xx,yy,zz 
    type(flowtype), dimension(:,:,:), allocatable :: bf 
    real(R_P), dimension(:, :, :, :), allocatable :: qq
    real(R_P) :: GAMMA=1.4d0, MA, Pr=0.72d0, Te, Re
    integer :: igs,jgs,kgs,igl,jgl,kgl,ige,jge,kge
    integer :: is,js,ks,il,jl,kl,ie,je,ke
    DM :: disturbDA,coordDA,meshDA,DA
    character(len=256) :: gridfile 
    character(len=256) :: flowfile  
    character(len=256) :: turbfile
    character(len=256) :: initfile
    logical :: ksp_mat_free_flg
    logical :: init_guess_flg
    character(len=256) :: dir  
    integer :: solver_mode=0
    integer :: split_mode=0
    PetscInt :: in,jn,kn,ln
    PetscScalar :: fk=1.0d0 ! 松弛系数
    complex(R_P) :: Alpha
    complex(R_P) :: Omega
    complex(R_P) :: Beta
    integer :: lns_mode
    integer :: BC_type 
    Vec :: tinkle_bell ! the local vec variable using in matrix_free implementation
    integer :: rank   
    integer :: sink 
    Mat :: Dolphin ! KSP: 免矩阵形式的线性系统
    Vec :: Turtle
    Mat :: Whale ! KSP: 显式矩阵形式的线性系统
    Mat :: Shark ! SNES: 牛顿法框架下的雅各比矩阵
    Vec :: RHS ! KSP: 右端项 
end module mod_parameters
