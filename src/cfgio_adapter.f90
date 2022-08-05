!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: cfgio_adapter.f90
!> @file
!> @breif cfgio代码适配器文件.
!  DESCRIPTION:
!     call cfg%set("domain", "in", in)
!>
!------------------------------------------------------------------------------
module mod_cfgio_adapter
    use cfgio_mod
    use mod_parameters
    implicit none
    public :: cfg_loader
    private
    type(cfg_t), save :: cfg
    character(len=256) :: dir

    contains

    subroutine cfg_loader(cfgfn)
        implicit none
        character(len=*), intent(in) :: cfgfn
        if(rank==0)then
        cfg=parse_cfg(trim(cfgfn))
        call cfg_filename(cfg)
        call cfg_domain(cfg)
        call cfg_freestream(cfg)
        call cfg_hlns(cfg)
        endif
    end subroutine cfg_loader

    subroutine cfg_filename(cfg)
        implicit none
        character(len=256) :: grid
        character(len=256) :: flow
        character(len=256) :: hdf5
        type(cfg_t) :: cfg

        if(cfg%has_key("filenames", "dir")) then
            call cfg%get("filenames", "dir", dir)
        else
            dir="./"
        endif

        call cfg%get("filenames", "grid", grid)
        call cfg%get("filenames", "flow", flow)
        call cfg%get("filenames", "hdf5", hdf5)
        call cfg%get("filenames", "output prefix", output_prefix)

        gridfile=trim(dir)//trim(grid)
        flowfile=trim(dir)//trim(flow)
        hdf5file=trim(dir)//trim(hdf5)
    end subroutine cfg_filename

    subroutine cfg_domain(cfg)
        implicit none
        type(cfg_t) :: cfg
        call cfg%get("domain", "in", in)
        call cfg%get("domain", "jn", jn)
        call cfg%get("domain", "kn", kn)
        if(kn==1)then
            lns_mode = 2
        else
            lns_mode = 3
        endif
    end subroutine cfg_domain

    subroutine cfg_freestream(cfg)
        implicit none
        type(cfg_t) :: cfg

        if(cfg%has_key("freestream", "ma"))then
            call cfg%get("freestream", "ma", Ma)
        else
            write(*,*) 'No Mach number is input,'
            Ma=0.001d0
            write(*,*) 'The Ma is set to 0.001.'
        endif

        if(cfg%has_key("freestream", "re"))then
            call cfg%get("freestream", "re", Re)
        else
            write(*,*) 'No Reynolds number is input,'
            Re=1.0d20
            write(*,*) 'The Re is set to 1.0d20.'
        endif

        if(cfg%has_key("freestream", "te"))then
            call cfg%get("freestream", "te", Te)
        else
            write(*,*) 'No freestream temperature is input,'
            Te=300.0d0
            write(*,*) 'The Te is set to 300K.'
        endif
    end subroutine cfg_freestream

    subroutine cfg_hlns(cfg)
        use penf,only:R_P
        implicit none
        real(R_P), allocatable, dimension(:) :: Ber
        character(len=256) :: hfile
        character(len=256) :: pfile
        character(len=256) :: ffile
        character(len=256) :: gfile
        character(len=256) :: ifile
        character(len=256) :: dfile
        character(len=32)  :: mode
        type(cfg_t) :: cfg
        integer :: npar
        npar=2

        call system("mkdir -p data")

        if(cfg%has_key("hlns", "Alpha"))then
            call cfg%get("hlns", "Alpha", Ber, npar)
            Alpha=cmplx(Ber(1),Ber(2),R_P)
        else
            write(*,*) "No Alpha is input."
            stop
        endif
        if(cfg%has_key("hlns", "Beta"))then
            call cfg%get("hlns", "Beta", Ber, npar)
            Beta=cmplx(Ber(1),Ber(2),R_P)
        else
            write(*,*) "No Beta is input."
            stop
        endif
        if(cfg%has_key("hlns", "Omega"))then
            call cfg%get("hlns", "Omega", Ber, npar)
            Omega=cmplx(Ber(1),Ber(2),R_P)
        else
            write(*,*) "No Omega is input."
            stop
        endif

        ! call cfg%get("hlns", "flow file", ffile)
        ! biflowfile=trim(ffile)

        ! call cfg%get("hlns", "grid file", gfile)
        ! bigridfile=trim(gfile)

        ! call cfg%get("hlns", "hdf5 file", hfile)
        ! hdf5file=trim(hfile)

        ! call cfg%get("hlns", "plot3d file", pfile)
        ! pltfile=trim(pfile)

        ! call cfg%get("hlns", "initial guess", ifile)
        ! initfile=trim(ifile)
        ! inquire(file=trim(initfile),exist=ex_ini_gus_flg)

        ! call cfg%get("hlns", "inlet", dfile)
        ! inletfile=trim(dfile)
        ! inquire(file=trim(inletfile),exist=inlet_file_flg)

    end subroutine cfg_hlns
end module mod_cfgio_adapter
