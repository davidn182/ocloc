program ARRIVAL_TIMES

!#############################################################################################
!#
!#
!#         Read time-averaged crosscorrelations computed from a large-N seismic array
!#      and determine the difference in arrival time between the direct wave at positive
!#     time and the direct wave at negative time. These measurements allow one to set up
!#     a system of equations which can be solved for potential timing errors of (some) of
!#      the stations. The details are given in the GJI publication ...
!#
!#      Modified by David Naranjo on Nov-13-2020
!#############################################################################################

use strings
use file_operations
use array_operations
use return_taper
use snr_routines
use ls_solutions
use, intrinsic :: iso_c_binding
use mpi

implicit none

include "sacf.h"

type cpl_vars
real(kind=4), dimension(2)                       :: center
real(kind=4), allocatable,dimension(:)           :: prop_online
real(kind=8)                                     :: azi, dist
integer(kind=4)                                  :: nwindows, no_xcorrwindows
integer(kind=8), allocatable,dimension(:)        :: cepsecs,mepsecs
character(len=3)                                 :: st1, st2
character(len=4), dimension(2)                   :: nw
character(len=7)                                 :: dir
character(len=500)                               :: path
end type


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

logical            :: file_spec,station_lt_trh

character(len=1)   :: slash
character(len=4)   :: taskid_char, wl_trh_char
character(len=5)   :: freq_char, lfreq_char, hfreq_char
character(len=7)   :: cpl_dir_tmp
character(len=8)   :: snr_trh_char
character(len=50)  :: buffst1ns, buffst1ew, buffst2ns, buffst2ew
character(len=50)  :: epsec_char, epsec_char1, epsec_char2, prop_online_char, no_xcorrwindowschar
character(len=250) :: crossfile_regex, crossfile_txt, crossfile_pos, crossfile_neg
character(len=250) :: cmd, dir_1, dir_2, dir_3, dir_4, dir_5, dir_6, dir_7, dir_8, dir_9, dir_10, dir_11, dir_20, dir_21, dir_22
character(len=250) :: dummy_char, xcorr_path, fname

character(len=350) :: st_info

real(kind=4)       :: lat1, lat2, lon1, lon2, xx, NaN, wl_trh, fwidth, slope, intercept
real(kind=4)       :: dvdf_beg, dvdf_end, v_beg, v_end, min_wl, BEGIN, DELTA
real(kind=4)       :: p_prem, ampint, c_maxampint, ac_maxampint
real(kind=8)       :: faz, edist, maxdist, dist_to_refcpl, lf, hf, t, sign_xcorr, temp_sign_xcorr, add_shift
real(kind=8)       :: v1, v2, c_snr, ac_snr, snr_trh, Ppvalu

integer(kind=4)    :: i, j, k, l, m, n, dum_ncpls, ncpls, ncpls_obs, ref_ncpls, nxcorrs, ill
integer(kind=4)    :: ns, ns_obs, ndt, nt, nn, nf, nftr, nf_all, nwl_trh, nsnr_trh, cpl_count
integer(kind=4)    :: apr_dt_shift, ind_cma, ind_acma, c_at, ac_at, peak_ind1, peak_ind2, trough_ind1, trough_ind2
integer(kind=4)    :: f1_mod_ind, f2_mod_ind, f1_data_ind, f2_data_ind, f_beg_ind, f_end_ind, ncepsecs, window, cepsec_ind
integer(kind=4)    :: rank, tag1, tag2, dummyint
integer(kind=4)    :: no_pzs, no_nzs, no_zs, nrl, nlen, nofts, nerr, no_peaks, no_troughs
integer(kind=4)    :: csl_length, int_length
integer(kind=4)    :: ios, loc, nargs, lda, ldb, info, ac1, c1, ac2, c2, ac2_env, c2_env, nsmpls_per_hp
integer(kind=4)    :: sig, sig_ext, sig_int, sig_int_ext, send, send_int
integer(kind=4)    :: numtasks,taskid,rc,ierr
integer(kind=4)    :: stat(MPI_STATUS_SIZE)

character(len=3)   , allocatable, dimension(:)      :: station, station_obs, sel_stations, station_obs_dum, cpl
character(len=4)   , allocatable, dimension(:)      :: nw
character(len=7)   , allocatable, dimension(:)      :: ref_cpl_dir, cpl_dir
character(len=100) , allocatable, dimension(:)      :: args, sacfiles
character(len=500) , allocatable, dimension(:)      :: ref_cpl_path

integer(kind=4)    , allocatable, dimension(:)      :: ipiv, refcpl_indices, cpl_indices
integer(kind=8)    , allocatable, dimension(:)      :: all_cepsecs,dum_all_cepsecs

real(kind=4)       , allocatable, dimension(:)      :: lat, lon, f_curve, c_curve, U_curve, Q_curve, v_curve, xcorr, xcorr_filt, wl_trhs, snr_trhs
real(kind=4)       , allocatable, dimension(:)      :: cfs, lfs, hfs, left_taper, right_taper, cfs_all
real(kind=4)       , allocatable, dimension(:)      :: sel_freqs, dummy_sel_freqs, sel_zeros, dummy_sel_zeros, dum_ref_curve
real(kind=4)       , allocatable, dimension(:,:)    :: ref_cpl_center,ref_curves
real(kind=8)       , allocatable, dimension(:)      :: ref_cpl_dist, weights
real(kind=8)       , allocatable, dimension(:)      :: peak_times, trough_times, peak_env, trough_env, c_sign_times, ac_sign_times
real(kind=8)       , allocatable, dimension(:)      :: int_c_sign, int_ac_sign
real(kind=8)       , allocatable, dimension(:)      :: Tobs, Tins ,Tobs_corr, apr_dt_init
real(kind=8)       , allocatable, dimension(:,:)    :: Tobs_lsqr, Tins_lsqr, Tins_var
real(kind=8)       , allocatable, dimension(:,:)    :: peaks, troughs, c_sign, ac_sign, A, B, result_shift, apr_dt_cpl
real(kind=8)       , allocatable, dimension(:,:,:)  :: cpl_shift
real(kind=8)       , allocatable, dimension(:,:,:,:):: apr_dt
real(kind=8)       , allocatable, dimension(:,:,:,:,:) :: R0

type(cpl_vars),      allocatable, dimension(:)      :: cpl_info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! PARAMETERS AND DEPENDENT VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer           , parameter :: int64 = selected_int_kind( 16 )
real(8)           , parameter :: pi=atan2(0.,-1.)
complex(8)        , parameter :: zero=cmplx(0.,0.),one=cmplx(1.),zi=cmplx(0.,1.)
real(8)           , parameter :: dt=0.04                      !!!!
integer(kind=4)   , parameter :: nsmpls=90000

!!! Write to file the filtered responses, window boundaries, etc
logical           , parameter :: wr_resp_details=.false.

!!! Filtering
real(kind=4)      , parameter :: fincr=0.01     ! dispersion curve (and eventual inversion) resolution
real(kind=4)      , parameter :: ftrans=0.12    !! width of frequency band used for transition from disp to reference curves (should be .ge. 5*fincr)
real(kind=4)      , parameter :: fmin=0.13      !! This value should be above 0.1 Hz, because of the computation of the reference velocity
real(kind=4)      , parameter :: fmax=0.75
!!! Minimum number of station couples a station must be part of
integer(kind=4)   , parameter :: ncpls_trh=3
!!! Minimum station separation to qualify at all for this analysis
real(kind=8)      , parameter :: min_sep=4000
!!! Range of minimum station separations in wavelengths
!! Uniform
!real(kind=4)      , parameter :: lwl_trh=1.0
!real(kind=4)      , parameter :: hwl_trh=1.0
!real(kind=4)      , parameter :: wl_trh_incr=0.2
!! Non-Uniform
!real(kind=4)      , parameter :: lwl_trh=1.0
!real(kind=4)      , parameter :: hwl_trh=6.0
!real(kind=4)      , parameter :: wl_trh_incr=0.25
!! Field data
real(kind=4)      , parameter :: lwl_trh=1.0
real(kind=4)      , parameter :: hwl_trh=4.0
real(kind=4)      , parameter :: wl_trh_incr=0.25

real(kind=4)      , parameter :: fwidth1=1.2    ! width (in octaves) at 0.15 Hz
real(kind=4)      , parameter :: fwidth2=0.8    ! width (in octaves) at 0.75 Hz
real(4)           , parameter :: twidth=0.04       ! width of taper (one side)
real(8)           , parameter :: dt_err=0.004      ! Resolution in time for which we determine the error

!!! Stations that whose timing is known to be correct (zero timing error), and who can therefore be eliminated from the system of equations.
character(len=3)  , parameter,dimension(30) :: stations_corr = (/'BER','EIN','GEV','HAH','HAS','HOS','KEF','KUG','LFE','ONG','PAT','PRE','RAH','RAR','RET','SDV','SKG','SKH','STA','SUH','NEW','HOP','KHR','KRV','SKF','STF','STK','VSV','ARN','MER'/)

character(len=3)  , parameter,dimension(21) :: stations_disc = (/'O08','O10','O21','O22','KEF','NYL','O07','O11','O06','O12','GEI','VOS','HDV','ELB','SEA','O01','O02','O16','O15','O17','O18'/)

integer           , parameter :: passes=2,order=4
double precision  , parameter :: transition_bandwidth=0.0, attenuation=0.0

real(8)           , parameter :: noise_st=240.     ! start of the noise window (in seconds) used for the computation of the SNR
!! Synthetic uniform
!real(4)           , parameter :: lsnr_trh=10.   ! Required signal-to-noise ratio of both causal and anti-causal part to be included
!real(4)           , parameter :: hsnr_trh=10.   ! Required signal-to-noise ratio of both causal and anti-causal part to be included
!real(4)           , parameter :: snr_trh_incr=10.   ! Required signal-to-noise ratio of both causal and anti-causal part to be included
!! Synthetic non-Uniform
!real(4)           , parameter :: lsnr_trh=5.   ! Required signal-to-noise ratio of both causal and anti-causal part to be included
!real(4)           , parameter :: hsnr_trh=60.   ! Required signal-to-noise ratio of both causal and anti-causal part to be included
!real(4)           , parameter :: snr_trh_incr=5.   ! Required signal-to-noise ratio of both causal and anti-causal part to be included
!! Field data
real(4)           , parameter :: lsnr_trh=5.   ! Required signal-to-noise ratio of both causal and anti-causal part to be included
real(4)           , parameter :: hsnr_trh=30.   ! Required signal-to-noise ratio of both causal and anti-causal part to be included
real(4)           , parameter :: snr_trh_incr=5.   ! Required signal-to-noise ratio of both causal and anti-causal part to be included

integer(kind=4)   , parameter :: ninv=3    ! Number of inversion strategies (currently 3: least squares, weighted least squares, and non-zero mean weighted least squares)

!!!! File containing station information (station code, lat, lon, etc.)
character(len=150), parameter  :: station_info='/vardim/home/weemstra/IMAGE/station_info'
!!!! Reference dispersion curves
!character(len=150), parameter  :: disp_info='/vardim/home/weemstra/IMAGE/disp_PREM300.data'
character(len=150), parameter  :: disp_info='/vardim/home/weemstra/IMAGE/Q_U_c_vs_T_n0_3km_1.5-20s.txt'
character(len=150), parameter  :: ref_path='/vardim/home/weemstra/IMAGE/FIELD_RESULTS/DISP_CURVES/SEL_REF_CURVES/'
!!! The initial errors are estimated using long period global phases and can be found in the following file:
character(len=150), parameter  :: init_instr_err='/vardim/home/weemstra/IMAGE/FIELD_RESULTS/INSTR_CORR/apriori_dt_estimates'
!!! This code is suitable for both synthetic and field data.
!character(len=150), parameter  :: td_path='/vardim/home/weemstra/IMAGE/SYNTH_RESULTS/CORRS/INSTR_CORR/CIRCLE_UNI/SRATE_25_NO_STALTA_TD_WHITE/'
!character(len=150), parameter  :: td_path='/vardim/home/weemstra/IMAGE/SYNTH_RESULTS/CORRS/INSTR_CORR/CIRCLE_NUNI/SRATE_25_NO_STALTA_TD_WHITE/'
character(len=150), parameter  :: td_path='/vardim/home/weemstra/IMAGE/FIELD_RESULTS/CORRS/INSTR_CORR/SRATE_25_NO_STALTA_TD_WHITE/'

!!! Location to which results are written
!character(len=150), parameter  :: wpath='/vardim/home/weemstra/IMAGE/SYNTH_RESULTS/CIRCLE_UNI/'
!character(len=150), parameter  :: wpath='/vardim/home/weemstra/IMAGE/SYNTH_RESULTS/CIRCLE_NUNI/'
!character(len=150), parameter  :: wpath='/vardim/home/weemstra/IMAGE/SYNTH_RESULTS/CIRCLE_NUNI_ALL/'
character(len=150), parameter  :: wpath='/vardim/home/weemstra/IMAGE/FIELD_RESULTS/INSTR_CORR/SINGLE_INVERSION/'

integer, parameter :: master=0                      ! The master task

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MPI initialization !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call MPI_INIT(ierr)
if(ierr .ne. MPI_SUCCESS)then
write(0,'(A)')'Error starting MPI program. Terminating.'
call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
end if
call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierr)
if (taskid .eq. master) then
print '(A16,I5.5,A9,I5.5)', 'Number of tasks=',numtasks,' My rank=',taskid
endif

if(numtasks.le.1)then
write(0,'(A)')'ERROR: This code should be run in parallel (at least two processors)'
stop
endif

!call MPI_GET_PROCESSOR_NAME(hostname,length,ierr)
write(6,*)'MPI task',taskid,'has started...'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SETTING SOME VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nt  = nsmpls
nn  = nsmpls*2
allocate(xcorr(nn), xcorr_filt(nn))

xx=-1.
NaN=sqrt(xx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REINITIALIZE DIRECTORIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

loc=len_trim(ref_path)
slash=ref_path(loc:loc)
if(slash.ne.'/')then
write(0,'(A)')'ERROR: The directory in which the reference dispersion curves are placed  does not end with a slash'
stop
endif
loc=len_trim(wpath)
slash=wpath(loc:loc)
if(slash.ne.'/')then
write(0,'(A)')'ERROR: The directory in which results will be written does not end with a slash'
stop
endif
loc=len_trim(td_path)
slash=td_path(loc:loc)
if(slash.ne.'/')then
write(0,'(A)')'ERROR: The directory from which the cross-correlations will be read does not end with a slash'
stop
endif

dir_1=trim(adjustl(wpath))//'XCORR_TSHIFTS/'
dir_2=trim(adjustl(wpath))//'COMPUTED_TSHIFTS/'
dir_3=trim(adjustl(wpath))//'FILTERED_RESPS_AND_ENVELOPES/'
dir_20=trim(adjustl(wpath))//'MISFITS/'

if(taskid .eq. master)then
write(cmd,'(A,A)')'rm -r ',trim(adjustl(dir_1))
call system(cmd)
write(cmd,'(A,A)')'mkdir ',trim(adjustl(dir_1))
call system(cmd)
write(cmd,'(A,A)')'rm -r ',trim(adjustl(dir_2))
call system(cmd)
write(cmd,'(A,A)')'mkdir ',trim(adjustl(dir_2))
call system(cmd)
write(cmd,'(A,A)')'rm -r ',trim(adjustl(dir_3))
call system(cmd)
write(cmd,'(A,A)')'mkdir ',trim(adjustl(dir_3))
call system(cmd)
write(cmd,'(A,A)')'rm -r ',trim(adjustl(dir_20))
call system(cmd)
write(cmd,'(A,A)')'mkdir ',trim(adjustl(dir_20))
call system(cmd)
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! OBTAIN NECESSARY STATION INFORMATION !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if(taskid .eq. master)then
write(6,'(A)')' '
write(6,'(A,A,A)')'Reading station information file'
write(6,'(A)')' '
endif
call nlines(station_info,ns,taskid,ios)
if(ios.ne.0) stop 'ERROR: Number of lines of the station-info file cannot be determined (variable: station_info).'

nargs=20;allocate(args(nargs))
allocate(nw(ns),station(ns),lon(ns),lat(ns))
open (unit=11,file=station_info)
do i=1,ns
read(11,'(A350)')st_info
call compact(st_info)
call parse(st_info,' ',args,nargs)
nw(i)=trim(adjustl(args(2)))
station(i)=trim(adjustl(args(3)))
if(trim(adjustl(args(5))).ne.'-')then
call value(trim(adjustl(args(5))),lat(i),ios)
if(ios.ne.0)stop 'Error converting latitude string'
call value(trim(adjustl(args(6))),lon(i),ios)
if(ios.ne.0)stop 'Error converting longitude string'
endif
enddo
close(11)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! SYNTHETIC OR REAL DATA !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! This piece determined wheter we are dealing with synthetic data (1 or 3), or with field data (2)
!!! Depending on the type of data, different parts of the code will be exceuted
call parse(trim(td_path),'/',args,nargs)
do i=1,nargs
if(trim(args(i)) .eq. 'CIRCLE_UNI') then
ill=1
exit
elseif(trim(args(i)) .eq. 'CIRCLE_NUNI') then
ill=3
exit
else
ill=2
endif
enddo
deallocate(args)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! RETRIEVE REFERENCE STATION COUPLES' DISPERSION CURVES !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate( cpl(2) )

select case ( ill )
case ( 1,3 )

ref_ncpls=1
allocate( ref_cpl_dist(ref_ncpls), ref_cpl_center(ref_ncpls,2))
ref_cpl_center(1,1)=sum(lat(:))/ns
ref_cpl_center(1,2)=sum(lon(:))/ns
ref_cpl_dist(1)=80000.  !!! Arbitrary

case ( 2 )

if(taskid .eq. master)then
write(6,'(A,A,A)')'Computation of the centers of the reference station couples (for real data only)'
write(6,'(A)')' '
endif

call nfiles(ref_path,ref_ncpls,taskid,ios)
if(ios.ne.0) stop 'ERROR: The reference dispersion curve directory does not exist (variable: ref_path).'
allocate( ref_cpl_dir(ref_ncpls), ref_cpl_path(ref_ncpls), ref_cpl_dist(ref_ncpls), ref_cpl_center(ref_ncpls,2))
call retr_files(ref_path,ref_ncpls,ref_cpl_dir,taskid,ios)
do i=1,ref_ncpls
ref_cpl_path(i)=trim(ref_path)//ref_cpl_dir(i)
call parse(ref_cpl_dir(i), '_', cpl, dummyint )
do j=1,ns
if(station(j)==cpl(1))then
lat1=lat(j)
lon1=lon(j)
endif
if(station(j)==cpl(2))then
lat2=lat(j)
lon2=lon(j)
endif
enddo
ref_cpl_center(i,1)=(lat1+lat2)/2.
ref_cpl_center(i,2)=(lon1+lon2)/2.
call writenum(lat1,buffst1ns,'f16.8')
call writenum(lon1,buffst1ew,'f16.8')
call writenum(lat2,buffst2ns,'f16.8')
call writenum(lon2,buffst2ew,'f16.8')
call Inverse(buffst1ns,buffst1ew,buffst2ns,buffst2ew,edist,faz)
ref_cpl_dist(i)=edist
enddo

end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! COMPUTE (NUMBER OF) FREQUENCY BANDS !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Number of frequency bands for which we will evaluate the instrument response
nf=nint((fmax-fmin)/fincr)+1
allocate(cfs(nf), lfs(nf), hfs(nf), ref_curves(nf,ref_ncpls))
ref_curves=0.

!!! Mind you: This part is hard code for an fmin of 0.13 and fmax of 0.75
slope=(fwidth2-fwidth1)/0.62
intercept=fwidth1-(slope*0.13)
nfdo1: do i=1,nf
cfs(i)=fmin+(i-1)*fincr
fwidth=intercept+slope*cfs(i)
lfs(i)=(cfs(i)*2)/((2**fwidth)+1)
hfs(i)=cfs(i)*2-lfs(i)
!cfs(i)=fmin+(i-1)*fincr
!lfs(i)=cfs(i)-0.075
!hfs(i)=cfs(i)+0.075
enddo nfdo1

!!! Number of frequency bands for which we compute an inital reference dispersion curve
nf_all=10001
allocate(cfs_all(nf_all))
do i=1,nf_all
cfs_all(i)=(i-1)*0.001
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! READ SORT AND INTERPOLATE PREM VELOCITIES !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if(taskid .eq. master)then
write(6,'(A,A,A)')'Retrieving and interpolating dispersion curve'
write(6,'(A)')' '
endif

call nlines(disp_info,nrl,taskid,ios)
if(ios.ne.0) stop 'ERROR: Number of lines of the dispersion curve file cannot be determined (variable: disp_info).'

allocate(f_curve(nrl),Q_curve(nrl),U_curve(nrl),c_curve(nrl))
open (unit=11,file=disp_info)
do i=1,nrl
read(11,*)f_curve(i),Q_curve(i),U_curve(i),c_curve(i)
!f_curve(i)=1./period
!U_curve(i)=U_curve(i)*1000.
!c_curve(i)=c_curve(i)*1000.
enddo
close(11)
deallocate(Q_curve,U_curve)


j=0; f1_mod_ind=0; f2_mod_ind=0
allocate( v_curve(nf_all) )
do i=1,nf_all
if( cfs_all(i) .ge. f_curve(nrl) )then
if(f2_mod_ind .eq. 0) then
f2_mod_ind = i-1
endif
v_curve(i)=c_curve(nrl)
elseif( cfs_all(i) .gt. f_curve(1) )then
do while (f_curve(j+1) .lt. cfs_all(i))
j=j+1
enddo
v_curve(i)=c_curve(j)+( (c_curve(j+1)-c_curve(j))/(f_curve(j+1)-f_curve(j))) * (cfs_all(i)-f_curve(j))
else
f1_mod_ind = i+1
v_curve(i)=c_curve(1)
endif
enddo
deallocate(f_curve, c_curve)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! RETRIEVE REFERENCE STATION COUPLES' DISPERSION CURVES !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


select case ( ill )
case ( 1,3 )

k=0
do j=1,nf

do while ( cfs_all(k+1) .lt. cfs(j) )
k=k+1
enddo

do_ref_cpls_1: do i=1,ref_ncpls

ref_curves(j,i)=v_curve(k)+ ( (v_curve(k+1)-v_curve(k) ) / 0.001 ) * ( cfs(j)-cfs_all(k) )
!write(6,*)cfs(j),lfs(j),hfs(j),ref_curves(j,i)

enddo do_ref_cpls_1
enddo

case( 2 )

if(taskid .eq. master)then
write(6,'(A,A,A)')'Retrieving and interpolating reference dispersion curves'
write(6,'(A)')' '
endif

allocate(dum_ref_curve(nf_all))

nftr=nint(ftrans/fincr)

allocate(sel_freqs(10000), dummy_sel_freqs(10000), sel_zeros(10000), dummy_sel_zeros(10000))
allocate(ipiv(4), A(4,4), B(4,1));lda=4;ldb=4
do_ref_cpls_2: do i=1,ref_ncpls

call writenum(taskid,taskid_char,'I4.4')

crossfile_regex=trim(ref_cpl_path(i))//'/'//'total_*_spectrum_pos_cross'
crossfile_txt='/tmp/crossfile_'//trim(adjustl(taskid_char))//'.txt'
call system('ls '//trim(crossfile_regex)//' > '//trim(crossfile_txt),ios)
open(unit=20,file=trim(crossfile_txt))
read(20,'(A)')crossfile_pos
close(20)
call nlines(trim(crossfile_pos),no_pzs,taskid,ios)

crossfile_regex=trim(ref_cpl_path(i))//'/'//'total_*_spectrum_neg_cross'
crossfile_txt='/tmp/crossfile_'//trim(adjustl(taskid_char))//'.txt'
call system('ls '//trim(crossfile_regex)//' > '//trim(crossfile_txt),ios)
open(unit=20,file=trim(crossfile_txt))
read(20,'(A)')crossfile_neg
close(20)
call nlines(trim(crossfile_neg),no_nzs,taskid,ios)

no_zs=no_pzs+no_nzs

open(15,file=crossfile_pos,access='sequential',form='formatted',status='old')
do j=1,no_pzs
read(15,'(F9.6,x,F10.1)')dummy_sel_freqs(j), dummy_sel_zeros(j)
enddo
close(15)
open(15,file=crossfile_neg,access='sequential',form='formatted',status='old')
do j=1,no_nzs
read(15,'(F9.6,x,F10.1)')dummy_sel_freqs(no_pzs+j), dummy_sel_zeros(no_pzs+j)
enddo
close(15)

sel_freqs=dummy_sel_freqs
call qsortC(sel_freqs(1:no_zs))
do j=1,no_zs
do k=1,no_zs
if(dummy_sel_freqs(k).eq.sel_freqs(j))then
sel_zeros(j)=dummy_sel_zeros(k)
cycle
endif
enddo
enddo

k=0; f1_data_ind=0; f2_data_ind=0
do j=1,nf_all
if( cfs_all(j) .ge. sel_freqs(no_zs) )then
if(f2_data_ind .eq. 0) then
f2_data_ind = j-1
endif
dum_ref_curve(j)=sel_zeros(no_zs)
elseif( cfs_all(j) .gt. sel_freqs(1) )then
do while (sel_freqs(k+1) .lt. cfs_all(j))
k=k+1
enddo
dum_ref_curve(j)=sel_zeros(k)+( (sel_zeros(k+1)-sel_zeros(k))/(sel_freqs(k+1)-sel_freqs(k)))*(cfs_all(j)-sel_freqs(k))
else
f1_data_ind = j+1
dum_ref_curve(j)=sel_zeros(1)
endif
enddo
!write(6,*)f1_data_ind,cfs_all(f1_data_ind),f2_data_ind,cfs_all(f2_data_ind)
!write(6,*)f1_mod_ind,cfs_all(f1_mod_ind),f2_mod_ind,cfs_all(f2_mod_ind)

if(f1_mod_ind .lt. f1_data_ind)then

if( cfs_all(f1_data_ind) .lt. (ftrans/4.) )then
f_beg_ind  = 1
f_end_ind  = (nint( (ftrans/4.) / 0.001 ) )
elseif( cfs_all(f1_data_ind) .le. ftrans)then
f_beg_ind  = 1
f_end_ind  = f1_data_ind
elseif( ( cfs_all(f1_data_ind)-cfs_all(f1_mod_ind) ) .lt. (ftrans/4.))then
f_beg_ind  = f1_data_ind - nint( (ftrans/4.) / 0.001 )
f_end_ind  = f1_data_ind
elseif( ( cfs_all(f1_data_ind)-cfs_all(f1_mod_ind) ) .lt. ftrans)then
f_beg_ind  = f1_mod_ind
f_end_ind  = f1_data_ind
else
f_beg_ind  = f1_data_ind - nint( ftrans / 0.001 )
f_end_ind  = f1_data_ind
endif
!write(6,*)f_beg_ind,cfs_all(f_beg_ind),f_end_ind,cfs_all(f_end_ind)

elseif(f1_data_ind .le. f1_mod_ind)then

if( cfs_all(f1_data_ind) .lt. (ftrans/4.) )then
f_beg_ind  = 1
f_end_ind  = (nint( (ftrans/4.) / 0.001 ) )
else
f_beg_ind  = 1
f_end_ind  = f1_data_ind
endif

endif

v_beg=dum_ref_curve(f_beg_ind)
dvdf_beg=( dum_ref_curve(f_beg_ind+1)-dum_ref_curve(f_beg_ind) )/0.001
v_end=dum_ref_curve(f_end_ind)
dvdf_end=( dum_ref_curve(f_end_ind)-dum_ref_curve(f_end_ind-1) )/0.001

B(1,1)=v_beg
B(2,1)=v_end
B(3,1)=dvdf_beg
B(4,1)=dvdf_end
A(1,1)=cfs_all(f_beg_ind)**3
A(1,2)=cfs_all(f_beg_ind)**2
A(1,3)=cfs_all(f_beg_ind)
A(1,4)=1.
A(2,1)=cfs_all(f_end_ind)**3
A(2,2)=cfs_all(f_end_ind)**2
A(2,3)=cfs_all(f_end_ind)
A(2,4)=1.
A(3,1)=3*cfs_all(f_beg_ind)**2
A(3,2)=2*cfs_all(f_beg_ind)
A(3,3)=1.
A(3,4)=0.
A(4,1)=3*cfs_all(f_end_ind)**2
A(4,2)=2*cfs_all(f_end_ind)
A(4,3)=1.
A(4,4)=0.
ipiv=0
call DGESV(4,1,A,lda,ipiv,B,ldb,info)

!!! Fill the transition gap with the third-degree polynomial
do j=f_beg_ind,f_end_ind
dum_ref_curve(j)=real( B(1,1)*cfs_all(j)**3 + B(2,1)*cfs_all(j)**2 + B(3,1)*cfs_all(j) + B(4,1) ,kind=4)
enddo

k=0
do j=1,nf

do while ( cfs_all(k+1) .lt. cfs(j) )
k=k+1
enddo

ref_curves(j,i)=dum_ref_curve(k)+ ( (dum_ref_curve(k+1)-dum_ref_curve(k) ) / 0.001 ) * ( cfs(j)-cfs_all(k) )
!write(6,*)cfs(j),lfs(j),hfs(j),ref_curves(j,i)

enddo

enddo do_ref_cpls_2

deallocate(sel_freqs, dummy_sel_freqs, sel_zeros, dummy_sel_zeros)
deallocate(A,B,ipiv)
deallocate(ref_cpl_dir, ref_cpl_path)

deallocate(dum_ref_curve)

end select

deallocate(v_curve,cfs_all)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! READ STATION PAIRS AND COMPUTE DISTANCES !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if(taskid .eq. master)then
write(6,'(A,A,A)')'Compute distances between the stations'
write(6,'(A)')' '
endif

!!! First discard autocorrelations and station couples located too close together and station couples with two correct stations
call nfiles(td_path,dum_ncpls,taskid,ios)
allocate( cpl_dir(dum_ncpls), cpl_indices(dum_ncpls) )
call retr_files(td_path,dum_ncpls,cpl_dir,taskid,ios)
if(ios.ne.0) stop 'ERROR: The data directory does not exist (variable: td_path).'
ncpls=0
!allocate(sel_stations(15))   !!! TEMPORARY CREATED
do i=1,dum_ncpls
if(taskid .eq. master)then
write(6,'(A,x,A)')'1',cpl_dir(i)
endif
call parse(cpl_dir(i), '_', cpl, dummyint )

if(cpl(1).eq.cpl(2)) cycle

if( (any(cpl(1) .eq. stations_corr,1) ) .and. (any(cpl(2) .eq. stations_corr,1) ) ) cycle

!!! As such, a subset of stations is selected
!if( (any(cpl(1) .eq. stations_disc,1) ) .or. (any(cpl(2) .eq. stations_disc,1) ) ) cycle


!!! THESE FEW LINES SPEED UP THE TESTING !!!
!sel_stations(1)='RAR'
!sel_stations(2)='STH'
!sel_stations(3)='KEF'
!sel_stations(4)='O22'
!sel_stations(5)='O23'
!sel_stations(6)='KUG'
!sel_stations(7)='YRD'
!sel_stations(8)='YRA'
!sel_stations(9)='YRN'
!sel_stations(10)='HDV'
!sel_stations(11)='O18'
!sel_stations(12)='PAT'
!sel_stations(13)='VOG'
!sel_stations(14)='LSF'
!sel_stations(15)='O10'
!if((.not. any(cpl(1) .eq. sel_stations,1) ) .or. (.not. any(cpl(2) .eq. sel_stations,1) )) cycle
!!! THESE FEW LINES SPEED UP THE TESTING !!!


do j=1,ns
if(station(j)==cpl(1))then
lat1=lat(j)
lon1=lon(j)
endif
if(station(j)==cpl(2))then
lat2=lat(j)
lon2=lon(j)
endif
enddo

dummy_char=trim(td_path)//cpl_dir(i)
call nfiles(dummy_char,nxcorrs,taskid,ios)
if(nxcorrs.eq.0) then
write(6,'(A)')'No crosscorrelations in directory '//dummy_char
cycle
endif

call writenum(lat1,buffst1ns,'f16.8')
call writenum(lon1,buffst1ew,'f16.8')
call writenum(lat2,buffst2ns,'f16.8')
call writenum(lon2,buffst2ew,'f16.8')
call Inverse(buffst1ns,buffst1ew,buffst2ns,buffst2ew,edist,faz)

if(edist.lt.min_sep) cycle


ncpls=ncpls+1
cpl_indices(ncpls)=i
enddo
!deallocate(sel_stations)   !!! TEMPORARY CREATED

maxdist=0.
allocate(cpl_info(ncpls))
do i=1,ncpls
if(taskid .eq. master)then
write(6,'(A,x,A)')'2',cpl_dir(cpl_indices(i))
endif
cpl_info(i)%dir=cpl_dir(cpl_indices(i))
cpl_info(i)%path=trim(td_path)//cpl_info(i)%dir
call parse(cpl_info(i)%dir, '_', cpl, dummyint )
cpl_info(i)%st1=cpl(1)
cpl_info(i)%st2=cpl(2)
do j=1,ns
if(station(j)==cpl(1))then
lat1=lat(j)
lon1=lon(j)
cpl_info(i)%nw(1)=nw(j)
endif
if(station(j)==cpl(2))then
lat2=lat(j)
lon2=lon(j)
cpl_info(i)%nw(2)=nw(j)
endif
enddo
cpl_info(i)%center(1)=(lat1+lat2)/2.
cpl_info(i)%center(2)=(lon1+lon2)/2.
call writenum(lat1,buffst1ns,'f16.8')
call writenum(lon1,buffst1ew,'f16.8')
call writenum(lat2,buffst2ns,'f16.8')
call writenum(lon2,buffst2ew,'f16.8')
call Inverse(buffst1ns,buffst1ew,buffst2ns,buffst2ew,edist,faz)
cpl_info(i)%azi=faz
cpl_info(i)%dist=edist
if(edist.gt.maxdist) maxdist=edist
enddo
deallocate(nw,lat,lon)
deallocate(cpl, cpl_dir, cpl_indices)

fname=trim(dir_1)//'distances_and_azimuths'
open(unit=23,file=fname)
do i=1,ncpls
write(23,'(A7,x,F8.2,x,F7.0)')cpl_info(i)%dir,cpl_info(i)%azi,cpl_info(i)%dist
enddo
close(23)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! READ, PER STATION PAIR, THE FILES !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Determine the number of time-lapse time-averaged crosscorrelations per station couple
nargs=5;allocate(args(nargs))
do i=1,ncpls
if(taskid .eq. master)then
write(6,'(A,x,A,x,A)')'3',cpl_info(i)%st1,cpl_info(i)%st2
endif
call nfiles(cpl_info(i)%path,nxcorrs,taskid,ios)
cpl_info(i)%nwindows=nxcorrs-1
allocate(sacfiles(nxcorrs))
call retr_files(cpl_info(i)%path,nxcorrs,sacfiles,taskid,ios)
!!!!!!!!!!!!!!! COMMENTED OUT !!!!!!!!!!!!!!!!!!!
!  if(nxcorrs .gt. 1)then  ! which is true if not only total_*.SAC exists
!    allocate(cpl_info(i)%cepsecs(nxcorrs-1))
!    allocate(cpl_info(i)%mepsecs(nxcorrs-1))
!    allocate(cpl_info(i)%prop_online(nxcorrs-1))
!    do j=1,nxcorrs-1
!      call compact(sacfiles(j))
!      call parse(sacfiles(j),'_',args,nargs)
!      call value(args(1)(1:10),cpl_info(i)%cepsecs(j),ios)
!      call value(args(2)(1:10),cpl_info(i)%mepsecs(j),ios)
!      call value(args(3)(1:4),cpl_info(i)%prop_online(j),ios)
!      !if(taskid .eq. master)then
!      !  write(6,*)cpl_info(i)%cepsecs(j), cpl_info(i)%mepsecs(j), cpl_info(i)%prop_online(j)
!      !endif
!    enddo
!  endif
!!!!!!!!!!!!!!! END COMMENTED OUT !!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!! ADDED !!!!!!!!!!!!!!!!!!!
cpl_info(i)%nwindows=0

!!!!!!!!!!!!!!! END ADDED !!!!!!!!!!!!!!!!!!!

call parse(sacfiles(nxcorrs),'_',args,nargs)
dummy_char=trim(adjustl(args(2)))
call parse(dummy_char,'.',args,nargs)
call value(args(1)(1:6),cpl_info(i)%no_xcorrwindows,ios)
write(6,*)cpl_info(i)%path,sacfiles(nxcorrs),cpl_info(i)%no_xcorrwindows
deallocate(sacfiles)
enddo

!!! Determine the number of unique time-lapse time-averaged crosscorrelations
allocate(dum_all_cepsecs(100));dum_all_cepsecs=9000000000000_int64
ncepsecs=0
do i=1,ncpls
if(taskid .eq. master)then
write(6,'(A,x,A,x,A)')'4',cpl_info(i)%st1,cpl_info(i)%st2
endif
do j=1,cpl_info(i)%nwindows
if(.not. any(cpl_info(i)%cepsecs(j) .eq. dum_all_cepsecs,1) )then
ncepsecs=ncepsecs+1
dum_all_cepsecs(ncepsecs)=cpl_info(i)%cepsecs(j)
endif
enddo
enddo
allocate(all_cepsecs(ncepsecs))
call int_qsortC(dum_all_cepsecs(:))
all_cepsecs=dum_all_cepsecs(1:ncepsecs)
deallocate(dum_all_cepsecs)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! COMPUTE (NUMBER OF) WAVELENGTH THRESHOLDS AND SNR THRESHOLDS !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nwl_trh=nint((hwl_trh-lwl_trh)/wl_trh_incr)+1
allocate(wl_trhs(nwl_trh))
do i=1,nwl_trh
wl_trhs(i)=lwl_trh+(i-1)*wl_trh_incr
enddo

nsnr_trh=nint((hsnr_trh-lsnr_trh)/snr_trh_incr)+1
allocate(snr_trhs(nsnr_trh))
do i=1,nsnr_trh
snr_trhs(i)=lsnr_trh+(i-1)*snr_trh_incr
enddo

allocate(result_shift(nsnr_trh,nwl_trh))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! READ APRIORI DT ESTIMATES !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(apr_dt_init(ns));apr_dt_init=0.

if(taskid .eq. master)then

allocate(apr_dt(ns,nsnr_trh,nwl_trh,ncepsecs+1))
allocate(R0(nsnr_trh,nwl_trh,ncepsecs+1,nf,ninv));R0=0.

write(6,'(A)')' '
write(6,'(A,A,A)')'Reading apriori_dt_estimates'
write(6,'(A)')' '

select case (ill)

case ( 1,3 )

apr_dt=0.

case ( 2 )

apr_dt=0.

call nlines(init_instr_err,ndt,taskid,ios)
if(ios.ne.0) stop 'ERROR: Number of lines of the aprior_dt_estimates file cannot be determined (variable: init_instr_err).'

open (unit=11,file=init_instr_err)
do i=1,ndt

read(11,'(A350)')st_info
call compact(st_info)
call parse(st_info,' ',args,nargs)

do j=1,ns
if(trim(adjustl(args(1))).eq.station(j))then
call value(trim(adjustl(args(2))),apr_dt_init(j),ios)
if(ios.ne.0)stop 'Error apriori dt estimate string'
endif
enddo

enddo
close(11)

do j=1,ncepsecs+1
do k=1,nwl_trh
do l=1,nsnr_trh
apr_dt(:,l,k,j)=apr_dt_init
enddo
enddo
enddo

end select

endif
deallocate(args)
deallocate(apr_dt_init)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! DETERMINE RELEVANT REFERENCE DISPERSION CURVE !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! For each station couple, determine reference station couple whose dispersion curve is used to compute wavelength
allocate(refcpl_indices(ncpls))
do i=1,ncpls
dist_to_refcpl=100000000.
lat1=cpl_info(i)%center(1)
lon1=cpl_info(i)%center(2)

do j=1,ref_ncpls
lat2=ref_cpl_center(j,1)
lon2=ref_cpl_center(j,2)

if(lat1.eq.lat2 .and. lon1.eq.lon2)then
refcpl_indices(i)=j
exit
else
call writenum(lat1,buffst1ns,'f16.8')
call writenum(lon1,buffst1ew,'f16.8')
call writenum(lat2,buffst2ns,'f16.8')
call writenum(lon2,buffst2ew,'f16.8')
call Inverse(buffst1ns,buffst1ew,buffst2ns,buffst2ew,edist,faz)
if(edist.lt.dist_to_refcpl)then
dist_to_refcpl=edist
refcpl_indices(i)=j
endif
endif
enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! MAIN LOOP OVER FREQUENCIES !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Time domain tapers
nofts=int(nint(twidth*nt),4)
allocate(cpl_indices(ns))
allocate(left_taper(nofts),right_taper(nofts))
call init_leftcos_taper(nofts,left_taper)
call init_rightcos_taper(nofts,right_taper)

tag1=1; tag2=2
allocate(apr_dt_cpl(nsnr_trh,nwl_trh))
nfdo5:do j=1,nf

!!! Directory to store envelopes and filtered xcorrelations
lf=real(lfs(j),kind=8)
hf=real(hfs(j),kind=8)
call writenum(cfs(j),freq_char,'F5.3')
call writenum(lf,lfreq_char,'F5.3')
call writenum(hf,hfreq_char,'F5.3')

dir_4=trim(dir_3)//trim(lfreq_char)//'_'//trim(hfreq_char)//'/'

if(taskid .eq. master)then
write(cmd,'(A,A)')'mkdir ',trim(adjustl(dir_4))
call system(cmd)

endif

nsmpls_per_hp=nint(1/cfs(j)/2/dt)
sig=nint(1/cfs(j)/3/dt)
sig_ext=2*sig
sig_int=nint(sig*(dt/dt_err))
sig_int_ext=2*sig_int

do window=1,ncepsecs+1

if(window .eq. ncepsecs+1)then
epsec_char='total'
else
cycle
call writenum(all_cepsecs(window),epsec_char,'I10.10')
endif
dir_5=trim(dir_1)//trim(epsec_char)//'/'
dir_6=trim(dir_2)//trim(epsec_char)//'/'
dir_7=trim(dir_4)//trim(epsec_char)//'/'
dir_21=trim(dir_20)//trim(epsec_char)//'/'

if(taskid .eq. master)then
write(cmd,'(A,A)')'mkdir ',trim(adjustl(dir_7))
call system(cmd)
if(j.eq.1)then
write(cmd,'(A,A)')'mkdir ',trim(adjustl(dir_5))
call system(cmd)
write(cmd,'(A,A)')'mkdir ',trim(adjustl(dir_6))
call system(cmd)
write(cmd,'(A,A)')'mkdir ',trim(adjustl(dir_21))
call system(cmd)
endif

!!!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Iteration over station couples and windows to determine shift per couple !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(cpl_shift(nsnr_trh,nwl_trh,ncpls))

write(6,'(A)')' '
write(6,'(A)')'ITERATION OVER STATION COUPLES AND WINDOWS FOR COMPUTATION OF THE REQUIRED PHASE SHIFT AT'
write(6,'(A)')'FREQUENCY: '//trim(freq_char)
write(6,'(A)')'WINDOW: '//trim(epsec_char)
write(6,'(A)')' '

i=0
docpls1: do while ( i .lt. ncpls )

i=i+1

apr_dt_cpl=0.
do k=1,ns
if(cpl_info(i)%st1 .eq. station(k))then
apr_dt_cpl(:,:)=apr_dt_cpl(:,:) + 2*apr_dt(k,:,:,window)
elseif(cpl_info(i)%st2 .eq. station(k))then
apr_dt_cpl(:,:)=apr_dt_cpl(:,:) - 2*apr_dt(k,:,:,window)
endif
enddo

if(i .lt. numtasks)then
rank=i
call MPI_SEND(i,1,MPI_INTEGER4,rank,tag1,MPI_COMM_WORLD,ierr)
call MPI_SEND(apr_dt_cpl(:,:),nsnr_trh*nwl_trh,MPI_REAL8,rank,tag1,MPI_COMM_WORLD,ierr)
else
call MPI_RECV(k,1,MPI_INTEGER4,MPI_ANY_SOURCE,tag2,MPI_COMM_WORLD,stat,ierr)
rank=stat(MPI_SOURCE)
call MPI_RECV(result_shift,nsnr_trh*nwl_trh,MPI_REAL8,rank,tag2,MPI_COMM_WORLD,stat,ierr)
cpl_shift(:,:,k)=result_shift
write(6,*)'Process with the following rank is finished',rank
call MPI_SEND(i,1,MPI_INTEGER4,rank,tag1,MPI_COMM_WORLD,ierr)
call MPI_SEND(apr_dt_cpl(:,:),nsnr_trh*nwl_trh,MPI_REAL8,rank,tag1,MPI_COMM_WORLD,ierr)
endif
enddo docpls1
!!! Wrap up by receiving last messages
do rank=1,numtasks-1
i=i+1
call MPI_RECV(k,1,MPI_INTEGER4,rank,tag2,MPI_COMM_WORLD,stat,ierr)
call MPI_RECV(result_shift,nsnr_trh*nwl_trh,MPI_REAL8,rank,tag2,MPI_COMM_WORLD,stat,ierr)
cpl_shift(:,:,k)=result_shift
call MPI_SEND(i,1,MPI_INTEGER4,rank,tag1,MPI_COMM_WORLD,ierr)
enddo

!!!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Least-squares inversion for all separation and snr thresholds !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

wl_do_master:do m=1,nwl_trh

!!! wl_trh_char
wl_trh=wl_trhs(m)
call writenum(wl_trh,wl_trh_char,'F4.2')

snr_do_master:do n=1,nsnr_trh

!!! snr_trh_char
snr_trh=real(snr_trhs(n),kind=8)
call writenum(snr_trh,snr_trh_char,'F8.1')

!!! Move on if none of the couples qualified
ncpls_obs=count(cpl_shift(n,m,:) .lt. 50.)
if(ncpls_obs.lt.3)then
cycle snr_do_master
endif

!!! Eliminate station couples of which one station forms less than ncpls_trh couples with other stations.
station_lt_trh=.true.
do while (station_lt_trh .eqv. .true.)
station_lt_trh = .false.
do k=1,ns

if( .not. any(station(k) .eq. stations_corr,1) )then
cpl_count=0

do_ncpls1:do i=1,ncpls
if(cpl_shift(n,m,i) .lt. 50.)then
if(cpl_info(i)%st1 .eq. station(k) .or. cpl_info(i)%st2 .eq. station(k))then
cpl_count=cpl_count+1
cpl_indices(cpl_count)=i
endif
endif
enddo do_ncpls1

if(cpl_count .gt. 0 .and. cpl_count .lt. ncpls_trh)then
station_lt_trh = .true.
do i=1,cpl_count
cpl_shift(n,m,cpl_indices(i)) = 100.
enddo
endif

endif

enddo
enddo

!!! (Again) move on if none of the couples qualified
ncpls_obs=count(cpl_shift(n,m,:) .lt. 50.)
if(ncpls_obs.lt.3)then
cycle snr_do_master
endif

!!! Determine which stations are present
ns_obs=0
allocate(station_obs_dum(ns))
do k=1,ns
if( .not. any(station(k) .eq. stations_corr,1) )then
do_ncpls2:do i=1,ncpls
if(cpl_shift(n,m,i) .lt. 50.)then
if(cpl_info(i)%st1 .eq. station(k) .or. cpl_info(i)%st2 .eq. station(k))then
ns_obs=ns_obs+1
station_obs_dum(ns_obs)=station(k)
exit do_ncpls2
endif
endif
enddo do_ncpls2
endif
enddo

!!! Store those stations in an array
allocate(station_obs(ns_obs))
do k=1,ns_obs
station_obs(k)=station_obs_dum(k)
enddo
deallocate(station_obs_dum)

!!! Initialize matrix A
allocate(A(max(1,ncpls_obs),ns_obs));A=0.
do k=1,ns_obs
ncpls_obs=0
do i=1,ncpls
if(cpl_shift(n,m,i) .lt. 50.)then
ncpls_obs=ncpls_obs+1
if(cpl_info(i)%st1 .eq. station_obs(k))then
A(ncpls_obs,k)=2.
elseif(cpl_info(i)%st2 .eq. station_obs(k))then
A(ncpls_obs,k)=-2.
endif
endif
enddo
enddo

!!! Fill the right-hand side
allocate(Tobs(max(ncpls_obs,ns_obs)));Tobs=0.
allocate(weights(ncpls_obs));weights=0.
ncpls_obs=0
do i=1,ncpls
if(cpl_shift(n,m,i) .lt. 50.)then
ncpls_obs=ncpls_obs+1
Tobs(ncpls_obs)=cpl_shift(n,m,i)
weights(ncpls_obs)=cpl_info(i)%dist/maxdist   !!! normalization by maxdist is to get values that are not too large; doesn't change the set of equations
endif
enddo


allocate(Tins_lsqr(ns_obs,ninv), Tins_var(ns_obs,ninv), Tobs_lsqr(ncpls_obs,ninv) )

!!! The system of equations is required to be overdetermined (also if we add an additional unknown for the non-zero mean)
if(ncpls_obs .le. (ns_obs+1) )then
write(6,'(7A)')'At a frequency of ',freq_char,' Hz, and for SNR threshold ',snr_trh_char,' and minimum station separation threshold ',wl_trh_char,' the system of equation is not overdetermined: no proper inversion possible.'
deallocate(station_obs,A,Tobs,Tobs_lsqr,Tins_lsqr,Tins_var,weights)
cycle snr_do_master
endif

!!! Simple least squares
call simple_ls(ncpls_obs, ns_obs, A, Tobs, Tins_lsqr(:,1), Tobs_lsqr(:,1), Tins_var(:,1), R0(n,m,window,j,1) )
if(R0(n,m,window,j,1) .eq. real(0.,kind=8))then         !!! This needs to be checked only once
write(6,'(7A)')'At a frequency of ',freq_char,' Hz, and for SNR threshold ',snr_trh_char,' and minimum station separation threshold ',wl_trh_char,' there is an isolated station: no inversion possible.'
deallocate(station_obs,A,Tobs,Tobs_lsqr,Tins_lsqr,Tins_var,weights)
cycle snr_do_master
endif

!!! weighted least squares
call weighted_ls(ncpls_obs, ns_obs, A, Tobs, weights, Tins_lsqr(:,2), Tobs_lsqr(:,2), Tins_var(:,2), R0(n,m,window,j,2) )

!!! weighted least squares, but accounting for the deviation of zero of the mean of the noise
call nzm_weighted_ls(ncpls_obs, ns_obs, A, Tobs, weights, Tins_lsqr(:,3), Tobs_lsqr(:,3), Tins_var(:,3), R0(n,m,window,j,3) )

!!!! weighted least squares, but accounting for the deviation of zero of the mean of the noise (but slightly different)
!call alt_nzm_weighted_ls(ncpls_obs, ns_obs, A, Tobs, weights, Tins_lsqr(:,4), Tobs_lsqr(:,4), Tins_var(:,4), R0(n,m,window,j,4) )

!!! Update the apriori time shift estimates using the results obtained with the best inversion strategy for that illumination pattern
select case(ill)

case(1)

do i=1,ns
do k=1,ns_obs
if(station_obs(k).eq.station(i))then
apr_dt(i,n,m,window)=Tins_lsqr(k,1)
endif
enddo
enddo

case(2,3)

do i=1,ns
do k=1,ns_obs
if(station_obs(k).eq.station(i))then
apr_dt(i,n,m,window)=Tins_lsqr(k,3)
endif
enddo
enddo

end select


!!! If synthetic data, read actual time shifts in order to compute correct couple time shifts
allocate(Tins(ns_obs),Tobs_corr(ncpls_obs))
Tins=0.;Tobs_corr=0.

if(ill .eq. 1 .or. ill .eq. 3)then

fname=trim(wpath)//'ACTUAL_TSHIFTS/time_shifts_0.1500'
call nlines(fname,ndt,taskid,ios)
open(unit=23,file=fname)
do k=1,ndt
read(23,'(A3,x,F6.3)')cpl_dir_tmp,add_shift
do i=1,ns_obs
if( station_obs(i) .eq. trim(cpl_dir_tmp) )then
Tins(i)=add_shift
exit
endif
enddo
enddo
close(23)

ncpls_obs=0
do i=1,ncpls
if(cpl_shift(n,m,i) .lt. 50.)then
ncpls_obs=ncpls_obs+1
Tobs_corr(ncpls_obs)=sum(A(ncpls_obs,:)*Tins(:))
endif
enddo

endif


!!! Write results per couple to file
dir_10=trim(dir_5)//trim(snr_trh_char)//'_'//trim(wl_trh_char)//'/'
inquire(file=dir_10, exist=file_spec)
if(.not. file_spec)then
cmd='mkdir '//trim(adjustl(dir_10))
call system(cmd)
endif

!!! inversion strategy 1
fname=trim(dir_10)//'simple_ls_cpl_shift_'//trim(freq_char); open(unit=23,file=fname)
ncpls_obs=0
do i=1,ncpls
if(cpl_shift(n,m,i) .lt. 50.)then
ncpls_obs=ncpls_obs+1
write(23,'(A7,x,F8.3,x,F8.3,x,F8.3)')cpl_info(i)%dir,cpl_shift(n,m,i),Tobs_lsqr(ncpls_obs,1),Tobs_corr(ncpls_obs)
endif
enddo
close(23)

!!! inversion strategy 2
fname=trim(dir_10)//'weighted_ls_cpl_shift_'//trim(freq_char); open(unit=23,file=fname)
ncpls_obs=0
do i=1,ncpls
if(cpl_shift(n,m,i) .lt. 50.)then
ncpls_obs=ncpls_obs+1
write(23,'(A7,x,F8.3,x,F8.3,x,F8.3)')cpl_info(i)%dir,cpl_shift(n,m,i),Tobs_lsqr(ncpls_obs,2),Tobs_corr(ncpls_obs)
endif
enddo
close(23)

!!! inversion strategy 3
fname=trim(dir_10)//'weighted_ls_nzm_noise_cpl_shift_'//trim(freq_char); open(unit=23,file=fname)
ncpls_obs=0
do i=1,ncpls
if(cpl_shift(n,m,i) .lt. 50.)then
ncpls_obs=ncpls_obs+1
write(23,'(A7,x,F8.3,x,F8.3,x,F8.3)')cpl_info(i)%dir,cpl_shift(n,m,i),Tobs_lsqr(ncpls_obs,3),Tobs_corr(ncpls_obs)
endif
enddo
close(23)

!!!! inversion strategy 4
!fname=trim(dir_10)//'weighted_ls_alt_nzm_noise_cpl_shift_'//trim(freq_char); open(unit=23,file=fname)
!ncpls_obs=0
!do i=1,ncpls
!  if(cpl_shift(n,m,i) .lt. 50.)then
!    ncpls_obs=ncpls_obs+1
!    write(23,'(A7,x,F8.3,x,F8.3,x,F8.3)')cpl_info(i)%dir,cpl_shift(n,m,i),Tobs_lsqr(ncpls_obs,4),Tobs_corr(ncpls_obs)
!  endif
!enddo
!close(23)


!!! Write results to file
dir_11=trim(dir_6)//trim(snr_trh_char)//'_'//trim(wl_trh_char)//'/'
inquire(file=dir_11, exist=file_spec)
if(.not. file_spec)then
cmd='mkdir '//trim(adjustl(dir_11))
call system(cmd)
endif

fname=trim(dir_11)//'time_shifts_'//trim(freq_char)
open(unit=23,file=fname)
do i=1,ns_obs
write(23,'(A3,x,F8.4,x,F12.8,x,F8.4,x,F12.8,x,F8.4,x,F12.8,x,F8.4)')station_obs(i),Tins_lsqr(i,1),Tins_var(i,1),Tins_lsqr(i,2),Tins_var(i,2),Tins_lsqr(i,3),Tins_var(i,3),Tins(i)
!write(23,'(A3,x,F8.4,x,F12.8,x,F8.4,x,F12.8,x,F8.4,x,F12.8,x,F8.4,x,F12.8,x,F8.4)')station_obs(i),Tins_lsqr(i,1),Tins_var(i,1),Tins_lsqr(i,2),Tins_var(i,2),Tins_lsqr(i,3),Tins_var(i,3),Tins_lsqr(i,4),Tins_var(i,4),Tins(i)
enddo
close(23)


dir_22=trim(dir_21)//trim(snr_trh_char)//'_'//trim(wl_trh_char)//'/'
inquire(file=dir_22, exist=file_spec)
if(.not. file_spec)then
cmd='mkdir '//trim(adjustl(dir_22))
call system(cmd)
endif
fname=trim(dir_22)//'misfits'
open(unit=23,file=fname,access='sequential',position='append')
write(23,'(A5,x,F20.5,x,F20.5,x,F20.5)')trim(freq_char),R0(n,m,window,j,1) ,R0(n,m,window,j,2) ,R0(n,m,window,j,3)
!write(23,'(A5,x,F20.5,x,F20.5,x,F20.5,x,F20.5)')trim(freq_char),R0(n,m,window,j,1) ,R0(n,m,window,j,2) , R0(n,m,window,j,3), R0(n,m,window,j,4)
close(23)

deallocate(station_obs,A,Tobs,Tobs_lsqr,Tins_lsqr,Tins_var,weights)
deallocate(Tins,Tobs_corr)
enddo snr_do_master
enddo wl_do_master

deallocate(cpl_shift)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! End of least-squares inversion !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

else
!!! Receive index of station couple from master (i)

slavedo: do while (tag1 .eq. 1)

call MPI_RECV(i,1,MPI_INTEGER4,master,tag1,MPI_COMM_WORLD,stat,ierr)

if(i .gt. ncpls) exit slavedo

call MPI_RECV(apr_dt_cpl,nsnr_trh*nwl_trh,MPI_REAL8,master,tag1,MPI_COMM_WORLD,stat,ierr)
result_shift=100.

!!! For testing purposes
!if(.not.  ((cpl_info(i)%dir.eq.'HOS_O17' .or. cpl_info(i)%dir.eq.'O17_HOS') &
!    .or.  (cpl_info(i)%dir.eq.'O17_O23' .or. cpl_info(i)%dir.eq.'O23_O17') &
!    .or.  (cpl_info(i)%dir.eq.'O23_O22' .or. cpl_info(i)%dir.eq.'O22_O23') &
!    .or.  (cpl_info(i)%dir.eq.'O22_O04' .or. cpl_info(i)%dir.eq.'O04_O22') &
!    .or.  (cpl_info(i)%dir.eq.'O04_O21' .or. cpl_info(i)%dir.eq.'O21_O04') &
!    .or.  (cpl_info(i)%dir.eq.'O21_O08' .or. cpl_info(i)%dir.eq.'O08_O21') &
!    .or.  (cpl_info(i)%dir.eq.'O08_O06' .or. cpl_info(i)%dir.eq.'O06_O08') &
!    .or.  (cpl_info(i)%dir.eq.'O06_HOS' .or. cpl_info(i)%dir.eq.'HOS_O06')))then
!  write(6,'(A)')'Following couple is discarded: '//trim(cpl_info(i)%dir)
!  result_shift=100.
!  call MPI_SEND(i,1,MPI_INTEGER4,master,tag2,MPI_COMM_WORLD,ierr)
!  call MPI_SEND(result_shift,1,MPI_REAL8,master,tag2,MPI_COMM_WORLD,ierr)
!  cycle slavedo
!endif

!if( (cpl_info(i)%st1.eq.'KEF') .and. (cpl_info(i)%st2.eq.'YRD') )then
!  write(6,'(A)')'Arrival-time difference of station couple '//trim(cpl_info(i)%dir)//' is evaluated'
!else
! write(6,'(A)')'Following couple is discarded: '//trim(cpl_info(i)%dir)
!  call MPI_SEND(i,1,MPI_INTEGER4,master,tag2,MPI_COMM_WORLD,ierr)
!  call MPI_SEND(result_shift,nsnr_trh*nwl_trh,MPI_REAL8,master,tag2,MPI_COMM_WORLD,ierr)
!  cycle slavedo
!endif


if(window.gt.ncepsecs)then
write(6,'(A)')'Arrival-time difference of station couple '//trim(cpl_info(i)%dir)//' is evaluated'
call writenum(cpl_info(i)%no_xcorrwindows,no_xcorrwindowschar,'I6.6')
fname='total_'//trim(adjustl(no_xcorrwindowschar))//'.SAC'
elseif(cpl_info(i)%nwindows .ge. 1 .and. any( all_cepsecs(window) .eq. cpl_info(i)%cepsecs ,1) )then
do k=1,cpl_info(i)%nwindows
if( cpl_info(i)%cepsecs(k) .eq. all_cepsecs(window) )then
cepsec_ind=k
endif
enddo
call writenum(cpl_info(i)%cepsecs(cepsec_ind),epsec_char1,'I10.10')
call writenum(cpl_info(i)%mepsecs(cepsec_ind),epsec_char2,'I10.10')
call writenum(cpl_info(i)%prop_online(cepsec_ind),prop_online_char,'F4.2')
fname=trim(epsec_char1)//'_'//trim(epsec_char2)//'_'//trim(adjustl(prop_online_char))//'.SAC'
else
write(6,'(A,x,I10)')'No crossscorrelation file exists for couple '//trim(cpl_info(i)%dir)//' at',all_cepsecs(window)
call MPI_SEND(i,1,MPI_INTEGER4,master,tag2,MPI_COMM_WORLD,ierr)
call MPI_SEND(result_shift,nsnr_trh*nwl_trh,MPI_REAL8,master,tag2,MPI_COMM_WORLD,ierr)
cycle slavedo
endif
xcorr_path=trim(cpl_info(i)%path)//'/'//trim(fname)

call rsac1(trim(xcorr_path),xcorr(:),nlen,BEGIN,DELTA,size(xcorr,1),nerr)
if(nerr .gt. 0) then    ! negative only a warning (e.g. nerr=-803 means number of samples is higher than size of array)
write(0,'(A6,I2,A18,A)') 'Error ',nerr,' reading in file: ',trim(xcorr_path)
stop
endif

dir_8=trim(dir_7)//trim(cpl_info(i)%dir)//'/'
cmd='mkdir '//trim(adjustl(dir_8))
call system(cmd)

!!! Velocities determining start and end of Signal window
v1=real(2*ref_curves(j,refcpl_indices(i)), kind=8)
v2=real(ref_curves(j,refcpl_indices(i))/2., kind=8)

!!! Filter xcorrelation (after tapering)
xcorr_filt=xcorr
xcorr_filt(1:nofts)=xcorr_filt(1:nofts)*left_taper
xcorr_filt(nn-nofts+1:nn)=xcorr_filt(nn-nofts+1:nn)*right_taper
call xapiir(xcorr_filt,nn,SAC_BUTTERWORTH,transition_bandwidth,attenuation,order,SAC_BANDPASS,lf,hf,dt,passes)

wl_do:do m=1,nwl_trh

!!! wl_trh_char
wl_trh=wl_trhs(m)
call writenum(wl_trh,wl_trh_char,'F4.2')

!!! Cycle if stations are not sufficiently separated at this frequncy
min_wl=ref_curves(j,refcpl_indices(i))/cfs(j)
if(cpl_info(i)%dist/min_wl .lt. wl_trh)then
write(6,'(A,F4.1,A)')'Station couple '//trim(cpl_info(i)%dir)//' does not exceed separation of ',wl_trh,' at this frequency'
cycle wl_do
endif

snr_do:do n=1,nsnr_trh

!!! snr_trh_char
snr_trh=real(snr_trhs(n),kind=8)
call writenum(snr_trh,snr_trh_char,'F8.1')

!!! Compute SNR accounting for the expected time shift
call SNR_with_shift(xcorr_filt,int(nt+1,kind=8),cpl_info(i)%dist,apr_dt_cpl(n,m)/2.,dt,noise_st,v1,v2,c_snr,ac_snr)

!!! Cycle if SNR's are not high enough at this frequncy
if(c_snr .lt. snr_trh .or. ac_snr .lt. snr_trh)then
write(6,'(A,F6.1,A)')'Station couple '//trim(cpl_info(i)%dir)//' does not exceed SNR of ',snr_trh,' at this frequency'
cycle snr_do
endif

!!! Compute envelope
apr_dt_shift=nint(apr_dt_cpl(n,m)/2./dt)
c2_env=int(nt+1)+nint(3*cpl_info(i)%dist/v2/dt)+apr_dt_shift
ac2_env=int(nt+1)-nint(3*cpl_info(i)%dist/v2/dt)+apr_dt_shift

!!! Determine number of peaks and troughs
no_peaks=0;no_troughs=0
do k=ac2_env,c2_env
if( xcorr_filt(k)-xcorr_filt(k-1) .gt. 0. .and. xcorr_filt(k+1)-xcorr_filt(k) .lt. 0. )then
no_peaks=no_peaks+1
endif
if( xcorr_filt(k)-xcorr_filt(k-1) .lt. 0. .and. xcorr_filt(k+1)-xcorr_filt(k) .gt. 0. )then
no_troughs=no_troughs+1
endif
enddo

!!! Store peaks and troughs and their times
allocate(peaks(4,no_peaks), troughs(4,no_troughs))
allocate(peak_times(no_peaks), trough_times(no_troughs))
no_peaks=0;no_troughs=0;peaks=0.;troughs=0.
do k=ac2_env,c2_env
if( xcorr_filt(k)-xcorr_filt(k-1) .gt. 0 .and. xcorr_filt(k+1)-xcorr_filt(k) .lt. 0 )then
no_peaks=no_peaks+1
peaks(1,no_peaks)=real(xcorr_filt(k),kind=8)
peak_times(no_peaks)=real(BEGIN,kind=8)+(k-1)*dt
endif
if( xcorr_filt(k)-xcorr_filt(k-1) .lt. 0 .and. xcorr_filt(k+1)-xcorr_filt(k) .gt. 0 )then
no_troughs=no_troughs+1
troughs(1,no_troughs)=real(xcorr_filt(k),kind=8)
trough_times(no_troughs)=real(BEGIN,kind=8)+(k-1)*dt
endif
enddo
!!! Construct envelope by fitting splines to peaks and troughs (separately)
allocate(peak_env(ac2_env:c2_env), trough_env(ac2_env:c2_env))
peak_env=0.;trough_env=0.
call cubspl ( peak_times(:), peaks(:,:), no_peaks, 0, 0 )
call cubspl ( trough_times(:), troughs(:,:), no_troughs, 0, 0 )
peak_ind1=0;peak_ind2=0
trough_ind1=0;trough_ind2=0
do k=ac2_env,c2_env

t= real(BEGIN,kind=8)+(k-1)*dt

if(t.gt.peak_times(1) .and. t.lt.peak_times(no_peaks))then
if(peak_ind1.eq.0)then
peak_ind1=k
endif
peak_ind2=k

peak_env(k)=Ppvalu(peak_times(:), peaks(:,:), no_peaks-1, 4, t, 0)
endif

if(t.gt.trough_times(1) .and. t.lt.trough_times(no_troughs))then
if(trough_ind1.eq.0)then
trough_ind1=k
endif
trough_ind2=k

trough_env(k)=Ppvalu(trough_times(:), troughs(:,:), no_troughs-1, 4, t, 0)
endif
enddo
deallocate(peak_times,trough_times,peaks,troughs)

!!! Determine start and ending indices of the signal windows
c1=int(nt+1)+nint(cpl_info(i)%dist/v1/dt)+apr_dt_shift   ! Starting index of causal signal window
c2=int(nt+1)+nint(cpl_info(i)%dist/v2/dt)+apr_dt_shift   ! Ending index causal signal window
ac1=int(nt+1)-nint(cpl_info(i)%dist/v1/dt)+apr_dt_shift  ! Ending index of acausal signal window
ac2=int(nt+1)-nint(cpl_info(i)%dist/v2/dt)+apr_dt_shift  ! Starting index acausal signal window

!!! Make sure that the signal windows are at least a period long
if( (c2-c1)*dt .lt. 1/cfs(j) )then
c2=c1+nint(1/cfs(j)/dt)
ac2=ac1-nint(1/cfs(j)/dt)
endif

if(wr_resp_details)then

dir_9=trim(dir_8)//trim(snr_trh_char)//'_'//trim(wl_trh_char)//'/'
cmd='mkdir '//trim(adjustl(dir_9))
call system(cmd)

!!! Write SNR bounds
fname=trim(dir_9)//'SNRs'
open(unit=30,file=fname)
write(30,'(x,F8.3)')c_snr
write(30,'(x,F8.3)')ac_snr
close(30)

!!! Write peaks
fname=trim(dir_9)//'top_envelope'
open(unit=30,file=fname)
do k=peak_ind1, peak_ind2
t=real(BEGIN,kind=8)+(k-1)*dt
write(30,'(x,F8.3,x,F10.3)')t,peak_env(k)
enddo
close(30)

!!! Write troughs
fname=trim(dir_9)//'bottom_envelope'
open(unit=30,file=fname)
do k=trough_ind1, trough_ind2
t=real(BEGIN,kind=8)+(k-1)*dt
write(30,'(x,F8.3,x,F10.3)')t,trough_env(k)
enddo
close(30)

!!! Write raw resp
fname=trim(dir_9)//'raw_resp'
open(unit=30,file=fname)
do k=ac2_env,c2_env
t=real(BEGIN,kind=8)+(k-1)*dt
write(30,'(x,F8.3,x,F10.3)')t,xcorr(k)
enddo
close(30)

!!! Write filtered_resp
fname=trim(dir_9)//'filtered_resp'
open(unit=30,file=fname)
do k=ac2_env,c2_env
t=real(BEGIN,kind=8)+(k-1)*dt
write(30,'(x,F8.3,x,F10.3)')t,xcorr_filt(k)
enddo
close(30)

!!! Write SNR bounds
fname=trim(dir_9)//'SNR_bounds'
open(unit=30,file=fname)
write(30,'(x,F8.3)')real(BEGIN,kind=8)+(ac2-1)*dt      !  start of acausal signal window
write(30,'(x,F8.3)')real(BEGIN,kind=8)+(ac1-1)*dt      !  end of acausal signal window
write(30,'(x,F8.3)')real(BEGIN,kind=8)+(c1-1)*dt       !  start of causal signal window
write(30,'(x,F8.3)')real(BEGIN,kind=8)+(c2-1)*dt       !  end of causal signal window
close(30)

endif

c_maxampint=0
do k=c1,c2
ampint=sum(peak_env(k-nsmpls_per_hp:k+nsmpls_per_hp)-trough_env(k-nsmpls_per_hp:k+nsmpls_per_hp))
if(ampint.gt.c_maxampint)then
c_maxampint=ampint
ind_cma=k
endif
enddo
ac_maxampint=0
do k=ac2,ac1
ampint=sum(peak_env(k-nsmpls_per_hp:k+nsmpls_per_hp)-trough_env(k-nsmpls_per_hp:k+nsmpls_per_hp))
if(ampint.gt.ac_maxampint)then
ac_maxampint=ampint
ind_acma=k
endif
enddo
deallocate(peak_env, trough_env)

!!! Recompute c1,c2,ac1,ac2
if(c_maxampint .gt. ac_maxampint)then
c_at = ind_cma
ac_at= 2*int(nt+1) - c_at + 2*apr_dt_shift
else
ac_at= ind_acma
c_at = 2*int(nt+1) - ac_at + 2*apr_dt_shift
endif
c1  =  c_at - nsmpls_per_hp   ! Starting index of causal signal window
c2  =  c_at + nsmpls_per_hp   ! Ending index causal signal window
ac1 = ac_at + nsmpls_per_hp   ! Ending index of acausal signal window
ac2 = ac_at - nsmpls_per_hp   ! Starting index acausal signal window

send=(c2-c1)+1
send_int=nint((send-1)*(dt/dt_err))+1

!!! Initialize causal signal (+ time vector)
csl_length=send+2*sig_ext
allocate( c_sign(4,csl_length), ac_sign(4,csl_length), c_sign_times(csl_length) ,ac_sign_times(csl_length))
do k=c1-sig_ext,c2+sig_ext
c_sign(1,k-(c1-sig_ext)+1)=real(xcorr_filt(k),kind=8)
c_sign_times(k-(c1-sig_ext)+1)=real(BEGIN,kind=8)+(k-1)*dt
enddo
do k=ac1+sig_ext,ac2-sig_ext,-1
ac_sign(1,ac1+sig_ext-k+1)=real(xcorr_filt(k),kind=8)
ac_sign_times(ac1+sig_ext-k+1)=-1.*(real(BEGIN,kind=8)+(k-1)*dt)
enddo

!!! Interpolate signal windows
call cubspl ( c_sign_times(:), c_sign(:,:), csl_length, 0, 0 )
call cubspl ( ac_sign_times(:), ac_sign(:,:), csl_length, 0, 0 )
int_length=send_int+2*sig_int_ext
allocate(int_c_sign(int_length), int_ac_sign(int_length))
int_c_sign=0.;int_ac_sign=0.
do k=1,int_length
t=c_sign_times(1)+(k-1)*dt_err
int_c_sign(k)=Ppvalu(c_sign_times(:), c_sign(:,:), csl_length-1, 4, t, 0)
t=ac_sign_times(1)+(k-1)*dt_err
int_ac_sign(k)=Ppvalu(ac_sign_times(:), ac_sign(:,:), csl_length-1, 4, t, 0)
enddo
!fname='/vardim/home/weemstra/tmp/c_sign.SAC'
!call wsac1(fname,real(int_c_sign(1:int_length),kind=4),int_length,real(c_sign_times(1)),real(dt_err),nerr)
!fname='/vardim/home/weemstra/tmp/ac_sign.SAC'
!call wsac1(fname,real(int_ac_sign(1:int_length),kind=4),int_length,real(ac_sign_times(1)),real(dt_err),nerr)

sign_xcorr=0.
do l=-1*sig_int,sig_int
temp_sign_xcorr=sum(int_c_sign(sig_int_ext+1+l:sig_int_ext+send_int+l)*int_ac_sign(sig_int_ext+1-l-1:sig_int_ext+send_int-l-1))
if(temp_sign_xcorr .gt. sign_xcorr)then
sign_xcorr=temp_sign_xcorr
add_shift=2*l*dt_err+dt_err
endif
temp_sign_xcorr=sum(int_c_sign(sig_int_ext+1+l:sig_int_ext+send_int+l)*int_ac_sign(sig_int_ext+1-l:sig_int_ext+send_int-l))
if(temp_sign_xcorr .gt. sign_xcorr)then
sign_xcorr=temp_sign_xcorr
add_shift=2*l*dt_err
endif
enddo
result_shift(n,m)=apr_dt_shift*2*dt+add_shift
!write(6,'(A,x,F8.3)')trim(cpl_info(i)%dir)//': Result (i.e., left-hand side of forward problem) is',result_shift(n,m)

if(wr_resp_details)then

!!! Write clock error (CE) window bounds
fname=trim(dir_9)//'CE_bounds_1'
open(unit=30,file=fname)
if(c_maxampint .gt. ac_maxampint)then
write(30,'(x,F8.3)')real(BEGIN,kind=8)+(c1-1)*dt
write(30,'(x,F8.3)')real(BEGIN,kind=8)+(c2-1)*dt
else
write(30,'(x,F8.3)')real(BEGIN,kind=8)+(ac2-1)*dt
write(30,'(x,F8.3)')real(BEGIN,kind=8)+(ac1-1)*dt
endif
close(30)

!!! Write clock error (CE) window bounds of the lower amplitude peak
fname=trim(dir_9)//'CE_bounds_2'
open(unit=30,file=fname)
if(c_maxampint .gt. ac_maxampint)then
write(30,'(x,F8.3)')real(BEGIN,kind=8)+(ac2-1)*dt
write(30,'(x,F8.3)')real(BEGIN,kind=8)+(ac1-1)*dt
else
write(30,'(x,F8.3)')real(BEGIN,kind=8)+(c1-1)*dt
write(30,'(x,F8.3)')real(BEGIN,kind=8)+(c2-1)*dt
endif
close(30)

!!! Write interpolated responses before adjustment
fname=trim(dir_9)//'interp_resp_before_corr'
open(unit=30,file=fname)
do k=1+sig_int_ext,1+send_int+sig_int_ext
t=c_sign_times(1)+(k-1)*dt_err
write(30,'(x,F8.3,x,F10.3,x,F10.3)')t,int_c_sign(k),int_ac_sign(k)
enddo
close(30)

!!! Write interpolated responses after adjustment
fname=trim(dir_9)//'interp_resp_after_corr'
open(unit=30,file=fname)
do k=1+sig_int_ext,1+send_int+sig_int_ext
t=c_sign_times(1)+(k-1)*dt_err+(add_shift/2.)
write(30,'(x,F8.3,x,F10.3,x,F10.3)')t,int_c_sign(k+nint(add_shift/2./dt_err)),int_ac_sign(k-nint(add_shift/2./dt_err))
enddo
close(30)

endif

deallocate(c_sign, ac_sign, c_sign_times, ac_sign_times)
deallocate(int_c_sign, int_ac_sign)

enddo snr_do

enddo wl_do

!!! Let master know that the work has been done and communicate result
call MPI_SEND(i,1,MPI_INTEGER4,master,tag2,MPI_COMM_WORLD,ierr)
call MPI_SEND(result_shift,nsnr_trh*nwl_trh,MPI_REAL8,master,tag2,MPI_COMM_WORLD,ierr)

enddo slavedo

endif

enddo

enddo nfdo5

do i=1,ncpls
if(allocated(cpl_info(i)%cepsecs))then
deallocate(cpl_info(i)%cepsecs)
deallocate(cpl_info(i)%mepsecs)
deallocate(cpl_info(i)%prop_online)
endif
enddo

if(taskid .eq. master)then
deallocate(apr_dt,R0)
endif

deallocate(station,apr_dt_cpl,wl_trhs,snr_trhs,result_shift)
deallocate(refcpl_indices,cpl_indices)
deallocate(left_taper,right_taper,xcorr,xcorr_filt)
deallocate(ref_cpl_dist, ref_cpl_center)
deallocate(cfs, lfs, hfs, ref_curves)
deallocate(cpl_info)


call MPI_BARRIER(MPI_COMM_WORLD,ierr)
if(ierr .ne. MPI_SUCCESS)then
write(0,'(A)')'Error in MPI_BARRIER. Terminating.'
call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
endif

call MPI_FINALIZE(ierr)
if(ierr .ne. MPI_SUCCESS)then
write(0,'(A)')'Error finalizing MPI program. Terminating.'
call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(0,'(A)') ' '
write(0,'(A)') '###############     PROGRAM COMPLETE      ################'
write(0,'(A)') ' '

end program TIMING_ERR_INV
