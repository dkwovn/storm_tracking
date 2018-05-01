!Program:
!    the main program to find feature points
!History:
!11/7/2014    cliu    1st version

Program find_feature_points
use tools,    only: neighbouring, &
                    tools_init, &
                    get_feature_points, &
                    init_type_frame, &
                    nml_file, field_orig, label_field, &
                    lon, lat, nx, ny, threshold, lon_rad, lat_rad, &
                    lat_min, lat_max, &
                    feature_point_option, frame, point
implicit none
character(len=150) :: dfile, ofile_field, ofile_list
character(len=5) :: yr_str
integer :: yr, t, nt, i, j, cnt
type(frame),allocatable :: feature_points_frame(:)
type(point),pointer :: p
real :: area(1000)
!=== define namelist variables ===
character(len=100) :: dfile_path = '/home/cliu/data/merra/6hourly/rvor/', &
                      ofile_field_path = '/home/cliu/data/merra/6hourly/rvor_cycl/',&
                      ofile_list_path = '/home/cliu/data/merra/6hourly/rvor_cycl/'
character(len=100) :: dfile_var = 'rvor850', &
                      ofile_field_var = 'rvor_cycl', &
                      ofile_list_var = 'rvor_cycl_list'
integer :: yr_st = 1979, &
            yr_end = 2010
!=== define namelist ===
namelist /files_nml/ dfile_path, ofile_field_path, ofile_list_path, &
                     dfile_var, ofile_field_var, ofile_list_var, &
                     yr_st, yr_end
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

call tools_init

open(unit=11,file=nml_file,status='old')
read(11,nml=files_nml)
close(11)
do yr = yr_st, yr_end
    print*,'yr loop begin'
    if(mod(yr,4)==0.and.mod(yr,100)/=0.or.mod(yr,400)==0)then
        nt=364
    else
        nt=364
    end if
     allocate(feature_points_frame(nt))
    call init_type_frame(feature_points_frame)
    write(yr_str,"('.',i4)")yr
    dfile=trim(dfile_path)//trim(dfile_var)//yr_str//'.bin'
    ofile_field=trim(ofile_field_path)//trim(ofile_field_var)//yr_str//'.bin'
    ofile_list=trim(ofile_list_path)//trim(ofile_list_var)//yr_str//'.bin'
    open(unit=10,file=dfile,form='unformatted',access='stream',status='old')
    open(unit=11,file=ofile_field,form='unformatted',access='stream')
    open(unit=12,file=ofile_list,form='unformatted',access='stream')
    
    call init_type_frame(feature_points_frame)
    print*,nx,ny
    do t = 1, nt
        print*,'t=',t
        cnt=0
        area=0
        label_field=0
        read(10)field_orig
        do i = 1, nx
            do j = 1, ny
                if(field_orig(i,j) > threshold .and. label_field(i,j) == 0 .and. lat(j)>lat_min)then
                    cnt=cnt+1
                    call neighbouring(field_orig,i,j,threshold,label_field,area(cnt),cnt)
                end if
            end do
        end do
        !=== filtering the candidates ===
        !=== find local maxima ===
        call get_feature_points(field_orig,label_field,feature_point_option, feature_points_frame(t), area)
        write(11)label_field*1.0
        if(feature_points_frame(t)%num == 0) cycle
        p=>feature_points_frame(t)%first
        do 
            if(associated(p))then
                write(12)t,p%lon,p%lat,p%value
            else
                exit
            end if
            p=>p%next
        end do
        !================================
    end do
    close(10)
    close(11)
    close(12)
    print*,yr
    deallocate(feature_points_frame)
end do
print*,yr_st,yr_end
stop
end 


