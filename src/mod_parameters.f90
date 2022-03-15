! #include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_parameters
    use penf, only: R_P
    use mod_flowtype
    use petsc
    implicit none
    public
    real(R_P), dimension(:, :, :), allocatable :: xi_xx,xi_yy,xi_zz,eta_xx,eta_yy,eta_zz,phi_xx,phi_yy,phi_zz ! 数组：度量系数
    real(R_P), dimension(:, :, :), allocatable :: xi_xy,xi_xz,xi_yz,eta_xy,eta_yz,eta_xz,phi_xy,phi_yz,phi_xz ! 数组：度量系数
    real(R_P), dimension(:, :, :), allocatable :: xi_x,xi_y,xi_z,eta_x,eta_y,eta_z,phi_x,phi_y,phi_z ! 数组：度量系数
    complex(R_P), dimension(:, :, :), allocatable :: disturb
    real(R_P), dimension(:, :, :), allocatable :: xx,yy,zz ! 数组：坐标
    type(flowtype), dimension(:,:,:), allocatable :: bf ! 数组：基本流信息
    real(R_P), dimension(:, :, :, :), allocatable :: qq ! 数组：基本流
    real(R_P) :: GAMMA=1.4d0, MA, Pr=0.72d0, Te, Re ! 流场参数
    integer :: igs,jgs,kgs,igl,jgl,kgl,ige,jge,kge ! MPI网格分块位置
    integer :: is,js,ks,il,jl,kl,ie,je,ke ! MPI网格分块位置
    DM :: disturbDA,coordDA,meshDA,DA ! DM.Object
    character(len=256) :: gridfile ! 文件名：网格文件
    character(len=256) :: flowfile ! 文件名：流场文件
    character(len=256) :: turbfile ! 文件名：边界文件
    character(len=256) :: initfile ! 文件名：初值文件
    logical :: ksp_mat_free_flg ! KSP中是否使用免矩阵生成矩阵
    logical :: init_guess_flg ! 是否赋初值
    character(len=256) :: dir ! 文件目录变量
    integer :: solver_mode=0 ! KSP.or.SNES
    integer :: split_mode=0 ! 对流系数矩阵拆分方式选择
    PetscInt :: in,jn,kn,ln ! 流场网格数、自由度数
    PetscScalar :: fk=1.0d0 ! 松弛系数
    complex(R_P) :: Alpha ! 波数
    complex(R_P) :: Omega ! 频率
    complex(R_P) :: Beta ! 波数
    integer :: lns_mode ! 2D-HLNS.or.3D-HLNS
    integer :: BC_type ! 边界条件类型：Dirichlet.or.Neumann.or.Robbin
    Vec :: tinkle_bell ! the local vec variable using in matrix_free implementation
    integer :: rank ! 进程编号
    integer :: sink ! 进程数
    Mat :: dolphin ! KSP: 免矩阵形式的左端矩阵
    Vec :: turtle ! 解
    Mat :: whale ! KSP: 显式矩阵形式的左端矩阵
    Mat :: shark ! SNES: 雅各比矩阵
    Vec :: RHS ! KSP: 右端项 
end module mod_parameters
