!Program:
!	the main program to do the cyclone tracking
!History:
!11/10/2014	cliu	1st version

Program cyclone_tracking_main
use io, 		only: init_io,	open_dfile, &
					  open_ofile, open_ofile_txt
use track_manager, only: track, &
						 track_init_modified, &
						 track_optim
use domain, 	only: yr_st, yr_end, nt, &
					  domain_init, &
					  t_st_leap, t_end_leap, t_st_norm, t_end_norm,&
					  lon_min, lon_max, &
						lat_min, lat_max
use tools, 		only: ifleap
implicit none

integer, parameter :: TERM = 6
integer :: dfile_unit, ofile_unit, ofile_txt_unit
integer :: yr, t, i, leap, t_st, t_end
integer :: swap_cnt, tt_cnt
integer :: foreward, backward, option
type(track), target :: track_head
type(track), pointer :: p, track_end, p_tmp
character(len=100) :: something

call domain_init
call init_io
do yr = yr_st, yr_end
	if(ifleap(yr))then
		nt = 1464
		t_st = t_st_leap
		t_end = t_end_leap
	else 
		nt = 1460
		t_st = t_st_norm
		t_end = t_end_norm
	end if
	print*,t_st,'->',t_end
	!=== read and initialize the tracks ===
	dfile_unit = open_dfile(yr)
	call track_init_modified(dfile_unit, track_head, track_end, t_st, t_end,&
	lon_min, lon_max, lat_min, lat_max)
	print*,'track_init ends'
	close(dfile_unit)
	!=== swap and optiminize the tracks ===
	option = 1
	foreward = 1
	backward = 1
	tt_cnt = 0
	do
		do
			swap_cnt = track_optim(track_head, track_end, t_st, t_end, option)
			if(option == 1) then
				if(swap_cnt == 0) then
					foreward = 0
					exit
				else
					backward = 1
					print*,'foreward: swap = ', swap_cnt
				end if
			else if(option == -1) then
				if(swap_cnt == 0) then
					backward = 0
					exit
				else
					foreward = 1
					print*,'backward: swap = ', swap_cnt
				end if
			end if
		end do
		if(foreward == 0.and.backward == 0)then
			print*,'exit loop'
			exit
		end if
		tt_cnt = tt_cnt + 1
		print*,'TERM=', TERM, 'tt_cnt=', tt_cnt
		if(tt_cnt > TERM) exit
		option = option * (-1)
!		read(*,*) something
	end do
	!============================================================
	!=== write tracks ===
100	ofile_unit = open_ofile(yr)
	ofile_txt_unit = open_ofile_txt(yr)
	print*,'before writing'
	p => track_head
	do 
		if(associated(p%next)) then
			p_tmp => p
			p => p%next
			if(.not.associated(p_tmp,track_head)) deallocate(p_tmp)
		else
			exit
		end if
		do t = t_st, t_end
			write(ofile_unit) t, p%track_point(t)%lon, &
							  p%track_point(t)%lat, p%track_point(t)%value
			write(ofile_txt_unit, *) t, p%track_point(t)%lon, &
							  p%track_point(t)%lat, p%track_point(t)%value
		end do
		write(ofile_unit) 1, 0.0, 0.0, -999.0
		write(ofile_txt_unit, *) 1, 0.0, 0.0, -999.0
	end do
	print*,'after writing'
	deallocate(track_head%track_point)
	track_head%next => null()
	close(ofile_unit)
	close(ofile_txt_unit)
end do
end
