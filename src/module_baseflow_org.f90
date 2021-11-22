module mod_baseflow_org 
    use penf, only: R_P
    implicit none
    type, public :: bf_flux_org_type
        real(R_P) :: rho=0.0d0
        real(R_P) :: x=0.0d0
        real(R_P) :: y=0.0d0
        real(R_P) :: z=0.0d0
        real(R_P) :: T=0.0d0 
    end type bf_flux_org_type
    type(bf_flux_org_type),parameter::BF_FLUX_NULL= &
    &bf_flux_org_type(0.0d0,0.0d0,0.0d0,0.0d0,0.0d0)
end module mod_baseflow_org