! #include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_parameters
    use penf, only: R_P
    use mod_flowtype
    use petsc
    implicit none
    public

    real(R_P), parameter, dimension(-2:2,-2:2) :: FDM_1nd_4ORD_CENTER=reshape( [&
        0.0d0       ,        0.0d0, -3.0d0/2.0d0, 2.0d0      ,  -1.0d0/2.0d0, &
        0.0d0       , -1.0d0/3.0d0, -1.0d0/2.0d0, 1.0d0      ,  -1.0d0/6.0d0, &
        1.0d0/12.0d0, -2.0d0/3.0d0,        0.0d0, 2.0d0/3.0d0, -1.0d0/12.0d0, &
        1.0d0/6.0d0 ,       -1.0d0,  1.0d0/2.0d0, 1.0d0/3.0d0,         0.0d0, &
        1.0d0/2.0d0 ,       -2.0d0,  3.0d0/2.0d0, 0.0d0      ,         0.0d0  &
        ], [5,5])
    real(R_P), parameter, dimension(-2:2,-2:2) :: FDM_2nd_4ORD_CENTER=reshape( [&
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0, &
        0.0d0        ,       1.0d0,        -2.0d0,       1.0d0,         0.0d0, &
        -1.0d0/12.0d0, 4.0d0/3.0d0, -15.0d0/6.0d0, 4.0d0/3.0d0, -1.0d0/12.0d0, &
        0.0d0        ,       1.0d0,        -2.0d0,       1.0d0,         0.0d0, &
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0  &
        ], [5,5])
    real(R_P), parameter, dimension(-2:2,-2:2) :: FDM_1nd_4ORD_Backward=reshape( [&
        0.0d0        ,       0.0d0,        0.0d0,       0.0d0,         0.0d0, &
        0.0d0        ,      -1.0d0,        1.0d0,       0.0d0,         0.0d0, &
        1.0d0/6.0d0  ,      -1.0d0,  1.0d0/2.0d0, 1.0d0/3.0d0,         0.0d0, &
        1.0d0/6.0d0  ,      -1.0d0,  1.0d0/2.0d0, 1.0d0/3.0d0,         0.0d0, &
        0.0d0        ,       0.0d0,        0.0d0,       0.0d0,         0.0d0  &
        ],[5,5])
    real(R_P), parameter, dimension(-2:2,-2:2) :: FDM_1nd_4ORD_Forward=reshape( [&
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0, &
        0.0d0        ,-1.0d0/3.0d0,  -1.0d0/2.0d0,       1.0d0,  -1.0d0/6.0d0, &
        0.0d0        ,-1.0d0/3.0d0,  -1.0d0/2.0d0,       1.0d0,  -1.0d0/6.0d0, &
        0.0d0        ,       0.0d0,        -1.0d0,       1.0d0,         0.0d0, &
        0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0  &
        ],[5,5])
    real(R_P), parameter, dimension(-1:1,-1:1) :: FDM_1nd_1ORD_Backward=reshape( [&
           0.0d0,           -1.0d0,                1.0d0, &
          -1.0d0,            1.0d0,                0.0d0, &
          -1.0d0,            1.0d0,                0.0d0  &
        ],[3,3])
    real(R_P), parameter, dimension(-1:1,-1:1) :: FDM_1nd_1ORD_Forward=reshape( [&
           0.0d0,           -1.0d0,                1.0d0, &
           0.0d0,           -1.0d0,                1.0d0, &
          -1.0d0,            1.0d0,                0.0d0  &
        ],[3,3])
    real(R_P), parameter, dimension(-1:1,-1:1) :: FDM_2nd_2ORD_CENTER=reshape( [&
           0.0d0,            0.0d0,                0.0d0, &
           1.0d0,           -2.0d0,                1.0d0, &
           0.0d0,            0.0d0,                0.0d0  &
        ], [3,3])
    real(R_P), parameter, dimension(-1:1,-1:1) :: FDM_1nd_2ORD_CENTER=reshape( [&
           0.0d0,          -1.0d0,                 1.0d0, &
    -1.0d0/2.0d0,           0.0d0,           1.0d0/2.0d0, &
          -1.0d0,           1.0d0,                 0.0d0  &
        ], [3,3])

    real(R_P), parameter :: delta_i(-2:2)=[0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]
    real(R_P), parameter :: delta_j(-2:2)=[0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]
    real(R_P), parameter :: delta_k(-2:2)=[0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]

    real(R_P), dimension(:,:,:), allocatable :: xi_xx,xi_yy,xi_zz,eta_xx,eta_yy,eta_zz,phi_xx,phi_yy,phi_zz ! 数组：度量系数
    real(R_P), dimension(:,:,:), allocatable :: xi_xy,xi_xz,xi_yz,eta_xy,eta_yz,eta_xz,phi_xy,phi_yz,phi_xz ! 数组：度量系数
    real(R_P), dimension(:, :,:), allocatable :: xi_x,xi_y,xi_z,eta_x,eta_y,eta_z,phi_x,phi_y,phi_z ! 数组：度量系数
    complex(R_P), dimension(:,:,:), allocatable :: disturb
    real(R_P), dimension(:,:,:), allocatable :: xx,yy,zz ! 数组：坐标
    type(flowtype), dimension(:,:,:), allocatable :: bf ! 数组：基本流信息
    real(R_P), dimension(:,:,:,:), allocatable :: qq ! 数组：基本流
    real(R_P) :: GAMMA=1.4d0, MA, Pr=0.72d0, Te, Re ! 流场参数
    integer :: igs,jgs,kgs,igl,jgl,kgl,ige,jge,kge ! MPI网格分块位置
    integer :: is,js,ks,il,jl,kl,ie,je,ke ! MPI网格分块位置
    real(R_P),parameter :: TAG = -1.0d0
    character(len=10) :: solver_mode ! KSP.or.SNES.or.newt.or.newt_sub
    character(len=256) :: bigridfile ! 文件名：流场文件
    character(len=256) :: biflowfile ! 文件名：流场文件
    character(len=256) :: turbfiles ! 文件名：边界文件
    PetscBool :: usedt,calculate_dt
    character(len=256) :: gridfile ! 文件名：网格文件
    character(len=256) :: flowfile ! 文件名：流场文件
    character(len=256) :: hdf5file ! 文件名：HDF5文件
    character(len=256) :: turbfile ! 文件名：边界文件
    character(len=256) :: initfile ! 文件名：初值文件
    character(len=256) :: pltfile ! 文件名：plot3d文件
    integer :: nx = PETSC_DECIDE
    integer :: ny = PETSC_DECIDE
    integer :: nz = PETSC_DECIDE
    character(len=7) :: io_type ! file I/O type
    logical :: init_guess_flg ! 是否赋初值
    PetscInt :: in,jn,kn,ln=5 ! 流场网格数、自由度数
    DM :: coordDA,meshDA,DA ! DM.Object
    integer :: split_mode=0 ! 对流系数矩阵拆分方式选择
    real(R_P) :: lm=1.0d0 ! 松弛系数
    complex(R_P) :: Alpha ! 波数
    complex(R_P) :: Omega ! 频率
    complex(R_P) :: Beta ! 波数
    real(R_P) :: dt=TAG ! 时间步长
    integer :: lns_mode ! 2D-HLNS.or.3D-HLNS
    real(R_P) :: cfl=30
    integer :: levels
    integer :: rank ! 进程编号
    integer :: sink ! 进程数
    Vec :: turtle ! 解向量
    Vec :: localx
    Mat :: whale ! 矩阵
    DM :: med1DA
    DM :: subDA
    Vec :: subx
    Vec :: RHS ! 右端项

end module mod_parameters
