!------------------------------------------------------------------------------
!
! Copyright (c) 2019-2024 Bzhouha 
! All rights reserved.
! 
! This file is part of START_LNS.
! START_LNS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! START_LNS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>. 
!
!------------------------------------------------------------------------------
!
! This work used the PETSc library, which is developed by the PETSc Development Team. See https://petsc.org/ for more information.
! Copyright (c) 2023, PETSc Development Team
! All rights reserved.
! This file is subject to the terms and conditions of the BSD 2-Clause License. See the file LICENSE in the top-level directory for more information.
!
!------------------------------------------------------------------------------
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

    integer, parameter :: wall_bc = 0 ! wall_bc = {'0';'797':'heat source'}
    real(kind=8) :: time,time0,time1

    real(R_P), dimension(:,:,:), allocatable :: xi_xx,xi_yy,xi_zz,eta_xx,eta_yy,eta_zz,phi_xx,phi_yy,phi_zz ! 数组：度量系数
    real(R_P), dimension(:,:,:), allocatable :: xi_xy,xi_xz,xi_yz,eta_xy,eta_yz,eta_xz,phi_xy,phi_yz,phi_xz ! 数组：度量系数
    real(R_P), dimension(:, :,:), allocatable :: xi_x,xi_y,xi_z,eta_x,eta_y,eta_z,phi_x,phi_y,phi_z ! 数组：度量系数
    complex(R_P), dimension(:,:,:), allocatable :: inlet
    complex(R_P), dimension(:,:,:), allocatable :: wall
    real(R_P), dimension(:,:,:), allocatable :: xx,yy,zz ! 数组：坐标
    type(flowtype), dimension(:,:,:), allocatable :: bf ! 数组：基本流信息
    real(R_P), dimension(:,:,:,:), allocatable :: qq ! 数组：基本流
    real(R_P) :: GAMMA=1.4d0, MA, Pr=0.713d0, Te, Re ! 流场参数
    integer :: igs,jgs,kgs,igl,jgl,kgl,ige,jge,kge ! MPI网格分块位置
    character(len=10) :: solver_mode='nasf' ! asf or nasf or lnasf
    integer :: is,js,ks,il,jl,kl,ie,je,ke ! MPI网格分块位置
    real(R_P),parameter :: TAG = -1.0d0
    character(len=256) :: output_prefix
    character(len=7) :: io_type='hdf5' ! file I/O type
    logical :: ex_ini_gus_flg=.False. ! 是否赋初值
    logical :: inlet_file_flg=.False. ! 是否读入入口
    character(len=256) :: bigridfile ! 文件名：流场文件
    character(len=256) :: biflowfile ! 文件名：流场文件
    character(len=256) :: biresfile ! 文件名：结果文件
    PetscBool :: usedt,calculate_dt
    character(len=256) :: inletfile ! 文件名：边界文件
    character(len=256) :: gridfile ! 文件名：网格文件
    character(len=256) :: flowfile ! 文件名：流场文件
    character(len=256) :: hdf5file ! 文件名：HDF5文件
    character(len=256) :: initfile ! 文件名：初值文件
    logical :: ini_gus_flg=.False. ! 是否赋初值
    character(len=256) :: pltfile ! 文件名：plot3d文件
    integer :: nx = PETSC_DECIDE
    integer :: ny = PETSC_DECIDE
    integer :: nz = PETSC_DECIDE
    PetscInt :: in,jn,kn,ln=5 ! 流场网格数、自由度数
    integer :: split_mode=0 ! 对流系数矩阵拆分方式选择
    real(R_P) :: lm=1.0d0 ! 松弛系数
    complex(R_P) :: Alpha ! 波数
    complex(R_P) :: Omega ! 频率
    complex(R_P) :: Beta ! 波数
    real(R_P) :: dt=TAG ! 时间步长
    integer :: lns_mode ! 2D-HLNS.or.3D-HLNS
    real(R_P) :: cfl=30
    integer :: levels
    DM :: meshDA,DA ! DM.Object
    integer :: rank ! 进程编号
    integer :: sink ! 进程数
    Vec :: turtle ! 解向量
    Mat :: whale ! 矩阵
    DM :: med1DA
    DM :: subDA
    Vec :: subx
    Vec :: RHS ! 右端项

end module mod_parameters
