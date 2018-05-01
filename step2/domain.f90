
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Module domain
use io, only: nml_unit, nml_file
implicit none
real,parameter:: pi = 3.1415927
!=== define namelist variables ===
integer :: nx = 288, &
		   ny = 144, &
		   yr_st = 1979, &
		   yr_end = 2010, &
		   t_st_norm = 1, &
		   t_end_norm = 1460, &
		   t_st_leap = 1, &
		   t_end_leap = 1464
real ::	lon0 = 0.625, &
		lat0 = -89.375, &
		dlon = 1.25, &
		dlat = 1.25, &
		lon_min = 0.0, &
		lon_max = 360.0, &
		lat_min = 20.0, &
		lat_max = 90.0
!=== define namelist ===
namelist /metrics_nml/ nx, ny, yr_st, yr_end, &
					   lon0, lat0, dlon, dlat, &
					   t_st_norm, t_end_norm, t_st_leap, t_end_leap, &
					   lon_min, lon_max, &
					   lat_min, lat_max
!===============================
real :: dlon_rad, dlat_rad
real, allocatable :: lon(:), lat(:), lon_rad(:), lat_rad(:)
integer :: nt

contains

subroutine domain_init
open(unit=nml_unit,file=nml_file,status='old')
read(unit=nml_unit, nml=metrics_nml)
close(nml_unit)
call calculate_lon_lat
end subroutine domain_init

subroutine calculate_lon_lat
integer :: i, j
dlon_rad = dlon*pi/180.0
dlat_rad = dlat*pi/180.0
allocate(lon(nx))
allocate(lat(ny))
allocate(lon_rad(nx))
allocate(lat_rad(ny))
do i = 1, nx
	lon(i) = lon0 + (i-1)*dlon
	lon_rad(i)=lon(i)*pi/180.0
end do
do j = 1, ny
	lat(j) = lat0 + (j-1)*dlat
	lat_rad(j)=lat(j)*pi/180.0
end do
end subroutine calculate_lon_lat


end Module domain
