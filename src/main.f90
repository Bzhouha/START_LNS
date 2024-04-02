!------------------------------------------------------------------------------
!
! Copyright (C) 2019-2024 Bzhouha
! 
! This file is part of START_LNS.
!
! START_LNS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! START_LNS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>. 
!
!------------------------------------------------------------------------------

!#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

program main
    use mod_solving
    use mod_files
    use petsc
    implicit none
    PetscErrorCode :: ierr

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

    call istream(PETSC_COMM_WORLD)
    call dstream(PETSC_COMM_WORLD)
    call ostream(PETSC_COMM_WORLD)

    call PetscFinalize(ierr)
end program main
