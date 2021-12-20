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