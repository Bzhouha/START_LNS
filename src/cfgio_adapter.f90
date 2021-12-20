!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: cfgio_adapter.f90
!> @file
!> @breif cfgio代码适配器文件.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
module mod_cfgio_adapter

    use cfgio_mod
    use mod_parameters

    implicit none

    public :: cfg_loader, cfg_writer

    private

    type(cfg_t), save :: cfg                            !<cfg类

    contains

    !> read infomation from cfg_file
    !! @param[in] cfgfn filename of cfgio
    subroutine cfg_loader(cfgfn)
        implicit none
        character(len=*), intent(in) :: cfgfn
        if(rank==0)then
        cfg=parse_cfg(trim(cfgfn))
        call cfg_filename(cfg)
        call cfg_freestream(cfg)
        !call cfg_bc(cfg)
        call cfg_hlns(cfg)
        endif
    end subroutine cfg_loader

    subroutine cfg_writer(cfgfn)
        implicit none
        character(len=*), intent(in) :: cfgfn
        cfg=parse_cfg(trim(cfgfn))
        call cfg%set("domain", "in", in)
        call cfg%set("domain", "jn", jn)
        call cfg%set("domain", "kn", kn)
        call cfg%write(cfgfn)
    end subroutine cfg_writer

    !> read filenames
    !> @param[in] cfg cfg_type
    subroutine cfg_filename(cfg)
        implicit none
        character(len=256) :: grid
        character(len=256) :: flow
        character(len=256) :: turb
        character(len=256) :: dir
        type(cfg_t) :: cfg

        if(cfg%has_key("filenames", "dir")) then
          call cfg%get("filenames", "dir", dir)
        else
          dir="./"
        endif

        call cfg%get("filenames", "grid", grid)
        call cfg%get("filenames", "flow", flow)
        call cfg%get("filenames", "disturb", turb)

        gridfile=trim(dir)//trim(grid)
        flowfile=trim(dir)//trim(flow)
        turbfile=trim(dir)//trim(turb)

    end subroutine cfg_filename

    !> 来流信息提取函数
    !> @param[in] cfg cfg类
    !> @return 返回配置文件freestream字段相关信息
    subroutine cfg_freestream(cfg)
        implicit none
        type(cfg_t) :: cfg

        if(cfg%has_key("freestream", "ma"))then
            call cfg%get("freestream", "ma", Ma)
        else
            print*, 'No Mach number is input,'
            Ma=0.001d0
            print*, 'The Ma is set to 0.001.'
        endif

        if(cfg%has_key("freestream", "re"))then
            call cfg%get("freestream", "re", Re)
        else
            print*, 'No Reynolds number is input,'
            Re=1.0d20
            print*, 'The Re is set to 1.0d20.'
        endif

        if(cfg%has_key("freestream", "te"))then
            call cfg%get("freestream", "te", Te)
        else
            print*, 'No freestream temperature is input,'
            Te=300.0d0
            print*, 'The Te is set to 300K.'
        endif
    end subroutine cfg_freestream

    subroutine cfg_bc(cfg)
        implicit none
        type(cfg_t) :: cfg

        if(cfg%has_key("Boundary Conditions", "BC_type"))then
            call cfg%get("Boundary Conditions", "BC_type", BC_type)
        else
            print*, 'No Boundary_Condition type is input,'
            BC_type=1
            print*, 'The Boundary_Condition type is set to Dirichlet.'
        endif
    end subroutine cfg_bc

    subroutine cfg_hlns(cfg)
        use penf,only:R_P
        implicit none
        real(R_P), allocatable, dimension(:) :: Ber
        character(len=256) :: init
        type(cfg_t) :: cfg
        integer :: npar

        npar=2

        if(cfg%has_key("hlns", "lns mode"))then
          call cfg%get("hlns", "lns mode", lns_mode)
        endif

        if(cfg%has_key("hlns", "Alpha"))then
          call cfg%get("hlns", "Alpha", Ber, npar)
          Alpha=cmplx(Ber(1),Ber(2),R_P)
        else
          print*, "No Alpha is input."
        endif

        if(cfg%has_key("hlns", "Beta"))then
          call cfg%get("hlns", "Beta", Ber, npar)
          Beta=cmplx(Ber(1),Ber(2),R_P)
        else
          print*, "No Beta is input."
        endif

        if(cfg%has_key("hlns", "Omega"))then
          call cfg%get("hlns", "Omega", Ber, npar)
          Omega=cmplx(Ber(1),Ber(2),R_P)
        else
          print*, "No Omega is input."
        endif

        call cfg%get("hlns", "initial guess", init)
        initfile=trim(init)

    end subroutine cfg_hlns
end module mod_cfgio_adapter
