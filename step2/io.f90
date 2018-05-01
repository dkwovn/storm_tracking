!Module:
!	subroutines for input/output, variables that read data or are output.
!History:
!11/10/2014	cliu	1st version

Module io
implicit none
private
public init_io, open_dfile, open_dfile2, &
	   open_ofile_txt, open_ofile, &
		nml_unit, nml_file
logical :: module_is_initialized = .false.
integer :: nml_unit, dfile_unit, dfile2_unit, ofile_unit, ofile_txt_unit
character(len=100) :: nml_file = 'cyclone_tracking_namelist.nml'
!=== define namelist variables ===
character(len=100) :: dfile_path = '/home/cliu/data/merra/6hourly/rvor_cycl/'
character(len=100) :: dfile_var	= 'rvor_cycl_list'
character(len=100) :: dfile2_path = '/home/cliu/data/merra/6hourly/rvor/'
character(len=100) :: dfile2_var = 'rvor850'
character(len=100) :: ofile_path = '/home/cliu/data/merra/6hourly/cycloTRACK/'
character(len=100) :: ofile_var	= 'cyclo_track'
character(len=100) :: ofile_txt_var = 'cyclo_track'
!=== define namelist ===
namelist /io_nml/	dfile_path, dfile_var, &
					dfile2_path, dfile2_var, &
					ofile_path, ofile_var, &
					ofile_txt_var
contains

subroutine init_io
nml_unit = 20
dfile_unit = 10
ofile_unit = 11
ofile_txt_unit = 12
dfile2_unit = 13
open(unit=nml_unit, file=nml_file, status='old')
read(unit=nml_unit, nml=io_nml)
close(nml_unit)
end subroutine init_io

function open_dfile(yr)
integer :: open_dfile, yr
character(len=150) :: dfile
character(len=5) :: yr_str
write(yr_str, "('.',i4)")yr
dfile=trim(dfile_path)//trim(dfile_var)//yr_str//'.bin'
open(unit=dfile_unit, file=dfile, form='unformatted', access='stream', status='old')
open_dfile=dfile_unit
end function open_dfile

function open_dfile2(yr)
integer :: open_dfile2, yr
character(len=150) :: dfile2
character(len=5) :: yr_str
write(yr_str, "('.',i4)")yr
dfile2=trim(dfile2_path)//trim(dfile2_var)//yr_str//'.bin'
open(unit=dfile2_unit, file=dfile2, form='unformatted', access='stream', status='old')
open_dfile2=dfile2_unit
end function open_dfile2

function open_ofile(yr)
integer :: open_ofile, yr
character(len=150) :: ofile
character(len=5) :: yr_str
write(yr_str, "('.',i4)")yr
ofile=trim(ofile_path)//trim(ofile_var)//yr_str//'.bin'
open(unit=ofile_unit, file=ofile, form='unformatted', access='stream')
open_ofile=ofile_unit
end function open_ofile

function open_ofile_txt(yr)
integer :: open_ofile_txt, yr
character(len=150) :: ofile_txt
character(len=5) :: yr_str
write(yr_str, "('.',i4)")yr
ofile_txt=trim(ofile_path)//trim(ofile_txt_var)//yr_str//'.txt'
open(unit=ofile_txt_unit, file=ofile_txt)
open_ofile_txt=ofile_txt_unit
end function open_ofile_txt

end 
