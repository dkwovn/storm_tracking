!Module:
!	functions and subroutines to manage tracks
!History:
!11/10/2014	cliu	1st version

Module track_manager
use io,		only: nml_unit, nml_file
use tools,	only: pi, distance_g, &
				  vector_cart, vector_product
implicit none
private
public	track, &
		track_init_modified, &
		track_optim

type point_chain
	real :: lon, lat, value
	real :: lon_rad, lat_rad
	type(point_chain),pointer :: next
end type
	
type point
	real :: lon, lat, value
	real :: lon_rad, lat_rad
	real :: cost
end type

type track
	type(point), allocatable :: track_point(:)
	type(track), pointer :: prev, next
	real :: track_length, track_cost
	integer :: track_code, points_num
end type

real,parameter :: eps = 1e-6
!=== define namelist variables ===
real :: phi_max = 0.8, &
		d_max_deg	= 5.0, &
		w1, w2
!=== define namelist ===
namelist /track_nml/ 	phi_max, d_max_deg, &
						w1, w2
!=======================
real :: d_max_rad 	
logical :: MODULE_IS_INITIALIZED = .false.
contains

subroutine track_manager_init
real :: lon1 = 140.0*pi/180.0
real :: lat1 = 50.0*pi/180.0
real :: lon2 = 141.5*pi/180.0
real :: lat2 = 52.0*pi/180.0
real :: lon3 = 143.5*pi/180.0
real :: lat3 = 53.1*pi/180.0
open(unit=nml_unit, file=nml_file, status='old')
read(unit=nml_unit, nml=track_nml)
close(nml_unit)
d_max_rad = d_max_deg * pi / 180.0
MODULE_IS_INITIALIZED = .true.
end subroutine track_manager_init

subroutine track_init_modified(dfile_unit, track_head, p_end, t_st, t_end)
integer,intent(in) :: dfile_unit, t_st, t_end
type(track),target,intent(inout) :: track_head
type(track),pointer,intent(inout) :: p_end

integer :: flag
integer :: t, last_t, error, track_cnt, track_code_tmp
real :: lon, lat, value, lon_rad, lat_rad
real :: d_rad_tmp_min, d_rad_tmp
type(point_chain),target :: init_frames(t_st:t_end)
type(point_chain),pointer :: p_init, p2_init, pmin_init
type(track), pointer :: p, p2

if(.not.MODULE_IS_INITIALIZED)then
	call track_manager_init
end if

!=================================
init_frames(:)%value = -999.0
track_cnt = 0
last_t = 0
p_init => init_frames(t_st)

do 
	read(unit=dfile_unit, iostat=error)t, lon, lat, value
	lon_rad = lon*pi/180.0
	lat_rad = lat*pi/180.0
	if(error/=0) then
		exit
	end if
	if(t<t_st) then
		cycle
	else if(t > t_end) then
		exit
	end if
	if(t == last_t) then
		allocate(p_init%next)
		p_init => p_init%next
		p_init%lon = lon
		p_init%lat = lat
		p_init%lon_rad = lon_rad
		p_init%lat_rad = lat_rad
		p_init%value   = value
		last_t = t
	else
		init_frames(t)%lon = lon
		init_frames(t)%lat = lat
		init_frames(t)%lon_rad = lon_rad
		init_frames(t)%lat_rad = lat_rad
		init_frames(t)%value   = value
		p_init%next => null()
		p_init => init_frames(t)
		last_t = t
	end if
end do
p_init%next => null()
!=== initialize the track_head ===
allocate(track_head%track_point(t_st:t_end))
track_head%track_point(:)%value = -998.0
track_head%track_point(:)%cost = phi_max
p_end => track_head
!=================================
p_init => init_frames(t_st)
do 
	call new_track(p_end, t_st, t_end, track_cnt)
	call assign_track_point_value(p_end, t_st, p_init%lon, p_init%lon_rad, &
								p_init%lat, p_init%lat_rad, value)
	!=== create a new phantom track for each real track ===
	if(associated(p_init%next)) then
		p_init => p_init%next
	else
		exit
	end if
end do

do t = t_st, t_end-1
	print*,'t=', t
	p => track_head
	do
		if(p%track_point(t)%value > -100) then
			!=== find the closest point for an existing track point ===
			d_rad_tmp_min = pi
			p2_init => init_frames(t+1)
			do
				d_rad_tmp = distance_g(p%track_point(t)%lon_rad, p%track_point(t)%lat_rad, &
									   p2_init%lon_rad, p2_init%lat_rad)
				if(d_rad_tmp < d_rad_tmp_min.and.p2_init%value>-100) then
					d_rad_tmp_min = d_rad_tmp
					pmin_init => p2_init
				end if
				if(associated(p2_init%next)) then
					p2_init => p2_init%next
				else 
					exit
				end if
			end do
			if(d_rad_tmp_min < d_max_rad) then
				call assign_track_point_value(p, t+1, pmin_init%lon, pmin_init%lon_rad, &
									  pmin_init%lat, pmin_init%lat_rad, pmin_init%value)
				pmin_init%value = -999.0
			else if(p%points_num == 1.and.t>t_st)then
				flag = 0
				p2 => track_head%next
				do
					if(p2%track_point(t-1)%value > -100)then
						d_rad_tmp = &
						distance_g(p%track_point(t)%lon_rad,p%track_point(t)%lat_rad,&
								p2%track_point(t-1)%lon_rad,p2%track_point(t-1)%lat_rad)
						if(d_rad_tmp < d_max_rad)then
							flag = 1
							exit
						end if
					end if
					if(associated(p2%next)) then
						p2 => p2%next
					else
						exit
					end if
				end do
				
				if(flag == 0)then
					!call delete_track(p, p_end, track_cnt)
				end if
			end if
	    end if
		if(associated(p%next)) then
			p => p%next
		else
			exit
		end if
	end do

	p2_init => init_frames(t+1)
	do
		if(p2_init%value > -100) then
			call new_track(p_end, t_st, t_end, track_cnt)
			call assign_track_point_value(p_end, t+1, p2_init%lon, p2_init%lon_rad, &
									  p2_init%lat, p2_init%lat_rad, p2_init%value)
		end if
		if(associated(p2_init%next)) then
			p2_init => p2_init%next
		else 
			exit
		end if
	end do
end do
!=====================================

print*,'track_cnt=', track_cnt
end subroutine

!=== use greedy exchange algorithm to optimize ===
function track_optim(track_head, p_end, t_st, t_end, option)
type(track),target,intent(in) :: track_head
type(track),pointer,intent(in) :: p_end
integer,intent(in) :: t_st, t_end, option
integer :: track_optim
!option=1 => forward
!option=-1 => backward
integer :: tt, t, exchange_flag
type(track), pointer :: p1, p2, pp1, pp2
real :: cost_change, cost_change_min 
!cost_change_min is supposed to be the negative change with the largest magnitude
!pp1, pp2 point to the pair to be swapped

if(.not.MODULE_IS_INITIALIZED)then
	call track_manager_init
end if
track_optim = 0
do tt = t_st+1, t_end-1
	if(option == 1)then
		t = tt
	else if(option == -1)then
		t = t_end - tt + t_st
	else
		print*,'error: wrong option for track_optim'
		stop
	end if
	cost_change_min = 0
	exchange_flag = 0
	p1 => track_head
	do
		p2 => p1%next
		do
			if(distance_g_point(p1%track_point(t),p2%track_point(t+option))<=d_max_rad.and. &
		       distance_g_point(p2%track_point(t),p1%track_point(t+option))<=d_max_rad) then
				cost_change = (cost(p1, p1, p2, t-option, t, t+option) + cost(p2, p2, p1, t-option, t, t+option)) &
							- (cost(p1, p1, p1, t-option, t, t+option) + cost(p2, p2, p2, t-option, t, t+option))
				if(cost_change < cost_change_min) then
					cost_change_min = cost_change
					pp1 => p1
					pp2 => p2
					exchange_flag = 1
				end if
			end if

			if(associated(p2%next)) then
				p2 => p2%next
			else
				exit
			end if
		end do
		if(associated(p1%next,p_end)) then
			exit
		else if(associated(p1%next)) then
			p1 => p1%next
		else
			print*,'Unexpected: zero track count'
		end if
	end do
	if(exchange_flag ==1)then
		call swap(pp1, pp2, t+option)
		track_optim = track_optim + 1
	end if
end do
end function 

subroutine swap(p1, p2, t)
type(track),pointer,intent(inout) :: p1, p2
integer,intent(in) :: t

type(point) :: temp

if(.not.MODULE_IS_INITIALIZED)then
	call track_manager_init
end if

temp = p1%track_point(t)
p1%track_point(t) = p2%track_point(t)
p2%track_point(t) = temp
end subroutine swap

!=== cost function is positive definite ===
function cost(p1, p2, p3, t1, t2, t3)
type(track),pointer,intent(in) :: p1, p2, p3
integer,intent(in) :: t1, t2, t3
real :: cost

real :: direction_dev, speed_dev
real :: distance1, distance2
real :: lon1, lon2, lon3
real :: lat1, lat2, lat3
real :: val1, val2, val3
!========================================
if(.not.MODULE_IS_INITIALIZED)then
	call track_manager_init
end if
!=== if at least one point is phantom point ===
val1 = p1%track_point(t1)%value
val2 = p2%track_point(t2)%value
val3 = p3%track_point(t3)%value

if(val1 < -100) then
	cost = 0 
	return
else if (.not. (val1 > -100 .and. val2 > -100 .and. val3 > -100)) then
	cost = phi_max
	return
end if

!=== if all the three points are real ===
lon1 = p1%track_point(t1)%lon_rad
lat1 = p1%track_point(t1)%lat_rad
lon2 = p2%track_point(t2)%lon_rad
lat2 = p2%track_point(t2)%lat_rad
lon3 = p3%track_point(t3)%lon_rad
lat3 = p3%track_point(t3)%lat_rad

if(lon1 == lon2 .and. lon2 == lon3 .and. &
    lat1 == lat2 .and. lat2 == lat3)then
    cost = 0
    return
else if(lon1 == lon2 .and. lat1 == lat2 .or. &
        lon2 == lon3 .and. lat2 == lat3)then
    cost = w2
    return
else
	direction_dev = 1 + vector_product(t_vector(lon2, lat2, lon1, lat1), &
								       t_vector(lon2, lat2, lon3, lat3))

	distance1 = distance_g(lon1, lat1, lon2, lat2)
	distance2 = distance_g(lon2, lat2, lon3, lat3)
	speed_dev = 1 - (2*sqrt(distance1*distance2)) / (distance1 + distance2)
	if(distance2 > d_max_rad) then
		cost = 2*phi_max
		return
	end if
end if
cost = 0.5*w1*direction_dev + w2*speed_dev
if(cost < -1*eps) then
	print*, 'ERROR: cost function being negative!'
	print*,distance1,distance2
	print*,cost
	stop
end if
end function


subroutine new_track(p, t_st, t_end, track_cnt)
type(track),pointer,intent(inout) :: p
integer,intent(in) :: t_st, t_end
integer,intent(inout) :: track_cnt

if(.not.MODULE_IS_INITIALIZED)then
	call track_manager_init
end if

track_cnt = track_cnt + 1
if(.not.associated(p%next))then
	allocate(p%next)
	p%next%prev => p
	p => p%next

	p%track_code = track_cnt
	allocate(p%track_point(t_st:t_end))
	p%track_point(:)%value = -999.0
	p%points_num = 0
	p%track_cost = 0
	p%next => null()
else
	print*,'Dangerous: p%next is already associated.'
	stop
end if
end subroutine

subroutine delete_track(p,p_end,track_cnt)
type(track),pointer,intent(inout) :: p, p_end
integer,intent(inout) :: track_cnt
type(track),pointer :: p_tmp
if(.not.associated(p%prev))then
	return
else if(.not.associated(p%next))then
	p%prev%next => null()
	p_tmp => p
	p => p%prev
	p_end => p
	deallocate(p_tmp)
	track_cnt = track_cnt - 1
else
	p%next%prev => p%prev
	p%prev%next => p%next
	p_tmp => p
	p => p%next
	deallocate(p_tmp)
	track_cnt = track_cnt - 1
end if
end subroutine

subroutine assign_track_point_value(p, t, lon, lon_rad, lat, lat_rad, value)
integer,intent(in) :: t
real,intent(in) :: lon, lon_rad, lat, lat_rad, value
type(track),pointer,intent(inout) :: p

if(.not.MODULE_IS_INITIALIZED)then
	call track_manager_init
end if

p%track_point(t)%lon = lon
p%track_point(t)%lat = lat
p%track_point(t)%value = value
p%track_point(t)%lon_rad = lon_rad
p%track_point(t)%lat_rad = lat_rad
p%points_num = p%points_num + 1
end subroutine assign_track_point_value

!=== t_vector is a unit vector ===
function t_vector(lon1, lat1, lon2, lat2)
real :: lon1, lon2, lat1, lat2
type(vector_cart) :: t_vector

real :: prod, x, y, z, norm

if(.not.MODULE_IS_INITIALIZED)then
	call track_manager_init
end if

prod = cos(lat1)*cos(lat2)*cos(lon1-lon2) + sin(lat1)*sin(lat2)
x = cos(lat2)*cos(lon2) - prod*cos(lat1)*cos(lon1)
y = cos(lat2)*sin(lon2) - prod*cos(lat1)*sin(lon1)
z = sin(lat2) - prod*sin(lat1)
norm = sqrt(x**2 + y**2 + z**2)

t_vector%x = x/norm
t_vector%y = y/norm
t_vector%z = z/norm
end function

function distance_g_point(point1, point2)
type(point),intent(in) :: point1, point2
real :: distance_g_point

real :: lon1, lon2, lat1, lat2
real :: val1, val2

if(.not.MODULE_IS_INITIALIZED)then
	call track_manager_init
end if

val1 = point1%value
val2 = point2%value

if(val1 < -100 .or. val2 < -100)then
	distance_g_point = d_max_rad
	return
end if

lon1 = point1%lon_rad
lat1 = point1%lat_rad
lon2 = point2%lon_rad
lat2 = point2%lat_rad

distance_g_point = acos(cos(lat1)*cos(lat2)*cos(lon1-lon2) + sin(lat1)*sin(lat2))
end function distance_g_point

end Module track_manager

