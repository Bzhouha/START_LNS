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

module mod_flowtype 
    use penf, only: R_P
    implicit none
    
    type, public :: basetype
        real(R_P) :: rho=0.0d0
        real(R_P) :: x=0.0d0
        real(R_P) :: y=0.0d0
        real(R_P) :: z=0.0d0
        real(R_P) :: T=0.0d0 
    end type basetype

    type(basetype),parameter::BF_FLUX_NULL= &
    &basetype(0.0d0,0.0d0,0.0d0,0.0d0,0.0d0)

    type, public :: flowtype
        type(basetype) :: BF   = BF_FLUX_NULL   !< 基本流通量\private
        type(basetype) :: BFDx = BF_FLUX_NULL   !< 基本流通量流向一阶导数\private
        type(basetype) :: BFDy = BF_FLUX_NULL   !< 基本流通量法向一阶导数\private
        type(basetype) :: BFDz = BF_FLUX_NULL   !< 基本流通量展向一阶导数\private
        type(basetype) :: BFDxx= BF_FLUX_NULL   !< 基本流通量流向二阶导数\private
        type(basetype) :: BFDyy= BF_FLUX_NULL   !< 基本流通量法向二阶导数\private
        type(basetype) :: BFDzz= BF_FLUX_NULL   !< 基本流通量展向二阶导数\private
        type(basetype) :: BFDxy= BF_FLUX_NULL   !< 基本流通量流向法向二阶导数\private
        type(basetype) :: BFDyz= BF_FLUX_NULL   !< 基本流通量法向展向二阶导数\private
        type(basetype) :: BFDxz= BF_FLUX_NULL   !< 基本流通量流向展向二阶导数\private
    end type flowtype
    
end module mod_flowtype