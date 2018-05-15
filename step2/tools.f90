!Program:
!    Subroutines and functions that are used in cyclone tracking
!History:
!11/07/2014 cliu    1st version

module tools
implicit none
public    neighbouring
public  tools_init
public  point
public  frame
public ifleap
public distance_g
public vector_cart
public pi
public vector_product

!=== customized type def ===
type point
integer :: x, y
real :: lon, lat
real :: value
type(point),pointer :: next
type(point),pointer :: prev
end type

type frame
type(point),pointer :: first
integer :: num
end type

type vector_cart
real::x, y, z
end type
!=== constants ===
real,parameter::pi=3.14159265359, eps=1e-6, a=6371000
character(len=100),parameter::nml_file='cyclone_tracking_namelist.nml'
!=== namelist variables ===
integer    :: nx        =    144, &
           ny        =     73 

real    :: lon0        =    0.0, &
           lat0        =    -90, &
           lon_min    =    0.0, &
           lon_max    =    360.0, &
           lat_min    =    20.0, &
           lat_max    =    90.0, &
           dlon        =    2.5, &
           dlat        =    2.5, &
           threshold    =    3e-5, &
           cntr_thre    =    5e-5, &
           area_thre    =    2e10

integer :: feature_point_option = 1
!=== allocable variables ===
real, allocatable :: field_orig(:,:)
real, allocatable :: lon(:), lat(:), lon_rad(:),lat_rad(:)
integer, allocatable :: label_field(:,:)
!===========================
logical :: MODULE_IS_INITIALIZED = .false.
!=== define namelist ===
namelist /metrics_nml/     nx, ny, &
                        lon_min, lon_max, &
                        lat_min, lat_max, &
                        dlon, dlat, &
                        lon0, lat0, &
                        threshold, cntr_thre, &
                        feature_point_option, &
                        area_thre
!=======================

contains

subroutine tools_init
open(unit=10,file=nml_file,status='old')
read(10,nml=metrics_nml)
close(10)
!=== allocate fields ===
allocate(field_orig(nx,ny))
allocate(label_field(nx,ny))
!=======================
call calculate_lon_lat

MODULE_IS_INITIALIZED = .true. 
print*,'Module tools is initialized.'
end subroutine tools_init

subroutine calculate_lon_lat
integer :: i, j
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

recursive subroutine neighbouring(field,x0,y0,threshold,label_field,area,label)
real,intent(in)::field(:,:),threshold
integer,intent(in)::x0,y0
integer,intent(inout)::label_field(:,:)
integer,intent(in) :: label
real,intent(inout) :: area

integer::en_x(0:7),en_y(0:7)
integer::nx,ny,i,j
nx=size(field,1)
ny=size(field,2)

label_field(x0,y0) = label
area = area + (dlon*pi/180.0*a*cos(lat_rad(y0))) * (dlat*pi/180.0*a)
call get_cord(x0,y0,en_x,en_y,nx,ny)
do i=0,7,2
    if(field(en_x(i),en_y(i))>=threshold.and.label_field(en_x(i),en_y(i))==0&
    .and.lat(y0)>=lat_min.and.lat(y0)<=lat_max)then
        call neighbouring(field,en_x(i),en_y(i),threshold,label_field,area,label)
    end if
end do

end subroutine neighbouring

subroutine get_feature_points(field_orig,label_field,option,feature_points_frame, area)
real,intent(in) :: field_orig(:,:)
integer,intent(inout) :: label_field(:,:)
type(frame),intent(inout) :: feature_points_frame
integer,intent(in) :: option
real,intent(in) :: area(:)

type(point),pointer :: p
integer :: i, j, k, en_x(0:7), en_y(0:7)
integer :: cnt
real :: mx_tmp, mn_tmp
if(.not.MODULE_IS_INITIALIZED)then
    print*,'Tools_mod is not initialized. get_feature_points'
    stop
end if
if(option == 1)then
    allocate(feature_points_frame%first)
    p=>feature_points_frame%first
    do i = 1, nx
        do j = 1, ny
            if(label_field(i,j) >= 1.and.area(label_field(i,j)) > area_thre)then
                call get_cord(i,j,en_x,en_y,nx,ny)
                mx_tmp = 0
                do k = 0, 7
                    if(field_orig(en_x(k),en_y(k)) < -999)then
                        mx_tmp = 999
                        exit
                    end if
                    if(field_orig(en_x(k),en_y(k)) >= mx_tmp)then
                        mx_tmp = field_orig(en_x(k),en_y(k))
                    end if
                end do
                if(field_orig(i,j) > mx_tmp.and.field_orig(i,j) > cntr_thre)then
                    mn_tmp = 0
                    cnt = 0
                    do k = 0, 7
                        if(field_orig(en_x(k),en_y(k))>0)then
                            mn_tmp = mn_tmp + field_orig(en_x(k),en_y(k))
                            cnt = cnt + 1
                        end if
                    end do
                    mn_tmp = mn_tmp/cnt
                    p%x = i
                    p%y = j
                    p%lon = lon(i)
                    p%lat = lat(j)
                    p%value = mn_tmp
                    feature_points_frame%num=feature_points_frame%num+1
                    allocate(p%next)
                    p%next%prev=>p
                    p=>p%next
                end if
            end if
        end do
    end do

    if(associated(p%prev))then
        p=>p%prev
        if(associated(p%next)) deallocate(p%next)
        p%next=>null()
    end if
end if
end subroutine get_feature_points

function ifleap(yr)
integer :: yr
logical :: ifleap
if(mod(yr,4) == 0 .and. mod(yr,100)/=0 .or. mod(yr,400) == 0)then
    ifleap = .TRUE.
end if
end function ifleap

subroutine init_type_frame(points_frame)
type(frame) :: points_frame(:)
integer :: nt, t
nt = size(points_frame,1)
do t = 1, nt
    points_frame(t)%num = 0
    points_frame(t)%first => null()
end do
end subroutine

subroutine get_cord(x,y,en_x,en_y,nx,ny)
implicit none
integer,intent(in)::nx,ny
integer,intent(in)::x,y
integer,intent(out)::en_x(0:7),en_y(0:7)
integer::k
real::dx,dy
do k=0,7
    if(abs(cos(pi/4*k))>eps)then
        dx=cos(pi/4*k)/abs(cos(pi/4*k))
    else
        dx=0
    end if
    if(abs(sin(pi/4*k))>eps)then
        dy=sin(pi/4*k)/abs(sin(pi/4*k))
    else
        dy=0
    end if
    en_x(k)=x+dx
    en_y(k)=y+dy
    if(en_x(k)>nx)then
        en_x(k)=en_x(k)-nx
    else if(en_x(k)<1)then
        en_x(k)=en_x(k)+nx
    end if
end do
end subroutine

function distance_g(lon1, lat1, lon2, lat2)
!=== Note that the lon, lat here are in radians ===
real :: lon1, lon2, lat1, lat2
real :: distance_g
distance_g = acos(cos(lat1)*cos(lat2)*cos(lon1-lon2) + sin(lat1)*sin(lat2))
end function distance_g


function vector_product(vector1, vector2)
type(vector_cart) :: vector1, vector2
real :: vector_product

vector_product = vector1%x*vector2%x + vector1%y*vector2%y + vector1%z*vector2%z
end function vector_product

end module tools
