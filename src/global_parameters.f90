! #include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module global_parameters
    use penf, only: R_P
    use petsc
    use bf_point_org
    implicit none
    real(R_P), dimension(:, :, :), allocatable :: xi_xx,xi_yy,xi_zz,eta_xx,eta_yy,eta_zz,phi_xx,phi_yy,phi_zz
    real(R_P), dimension(:, :, :), allocatable :: xi_xy,xi_xz,xi_yz,eta_xy,eta_yz,eta_xz,phi_xy,phi_yz,phi_xz 
    real(R_P), dimension(:, :, :), allocatable :: xi_x,xi_y,xi_z,eta_x,eta_y,eta_z,phi_x,phi_y,phi_z
    real(R_P), dimension(:, :, :), allocatable :: xx, yy, zz 
    type(bf_point_type),dimension(:,:,:),allocatable :: bf 
    real(R_P), dimension(:, :, :, :), allocatable :: qq
    complex(R_P),dimension(:,:,:),allocatable :: inflow
    real(R_P) :: GAMMA=1.4d0, MA, Pr=0.72d0, Te, Re
    integer :: igs,jgs,kgs,igl,jgl,kgl,ige,jge,kge
    integer :: is,js,ks,il,jl,kl,ie,je,ke
    character(len=256) :: FileLocation     
    character(len=256) :: gridfile 
    character(len=256) :: flowfile     
    DM :: meshDA, coordDA, DA
    PetscInt :: in,jn,kn,ln
    complex(R_P) :: Alpha
    complex(R_P) :: Omega
    complex(R_P) :: Beta
    integer :: BC_type 
    integer :: mode=0
    integer :: rank   
    integer :: size 
    Vec :: turtle
    Mat :: Dolphin
    Mat :: Whale      
end module global_parameters
