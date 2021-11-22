
! ===============================================
!                BaseFlow_Point
! -----------------------------------------------
! 这个模块创建一个结构，在结构中保存一个点的基本留信息和相应的导数信息。(rho,u,v,w,T)
!
!  所有子项使用bf_flux_org_type格式，即：
!      type, public :: bf_flux_org_type
!          real(R_P) :: rho     !< density
!          type(vector) :: vel  !< velocity
!          real(R_P) :: T       !< temperature
!      end type bf_flux_org_type
!
! ===============================================

module bf_point_org
    use mod_baseflow_org
    implicit none
    type, public :: bf_point_type
        type(bf_flux_org_type) :: BF   = BF_FLUX_NULL   !< 基本流通量\private
        type(bf_flux_org_type) :: BFDx = BF_FLUX_NULL   !< 基本流通量流向一阶导数\private
        type(bf_flux_org_type) :: BFDy = BF_FLUX_NULL   !< 基本流通量法向一阶导数\private
        type(bf_flux_org_type) :: BFDz = BF_FLUX_NULL   !< 基本流通量展向一阶导数\private
        type(bf_flux_org_type) :: BFDxx= BF_FLUX_NULL   !< 基本流通量流向二阶导数\private
        type(bf_flux_org_type) :: BFDyy= BF_FLUX_NULL   !< 基本流通量法向二阶导数\private
        type(bf_flux_org_type) :: BFDzz= BF_FLUX_NULL   !< 基本流通量展向二阶导数\private
        type(bf_flux_org_type) :: BFDxy= BF_FLUX_NULL   !< 基本流通量流向法向二阶导数\private
        type(bf_flux_org_type) :: BFDyz= BF_FLUX_NULL   !< 基本流通量法向展向二阶导数\private
        type(bf_flux_org_type) :: BFDxz= BF_FLUX_NULL   !< 基本流通量流向展向二阶导数\private
        contains
    end type bf_point_type
end module bf_point_org

