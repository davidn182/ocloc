program TIMING_ERR_INV

!#############################################################################################
!#                                              
!#
!#         Read time-averaged crosscorrelations computed from a large-N seismic array
!#      and determine the difference in arrival time between the direct wave at positive 
!#     time and the direct wave at negative time. These measurements allow one to set up 
!#     a system of equations which can be solved for potential timing errors of (some) of 
!#      the stations. The details are given in the GJI publication ... 
!#
!#
!#############################################################################################
use strings
use file_operations
use array_operations
use return_taper
use snr_routines
use, intrinsic :: iso_c_binding

implicit none

include "sacf.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(len=3)   :: station1, station2
character(len=4)   :: wl_trh_char
character(len=8)   :: snr_trh_char
character(len=250) :: cmd, dir_8, dir_9, dir_name
character(len=250) :: xcorr_path, fname, resp_details

real(kind=4)       :: ref_vel, dist_trh
real(kind=4)       :: BEGIN, DELTA
real(kind=8)       :: apr_dt_cpl, noise_st, lf, hf, t, sign_xcorr, temp_sign_xcorr, add_shift, cfs, min_wl, ampint, c_maxampint
real(kind=8)       :: result_shift, dt_err, dt, v1, v2, c_snr, ac_snr, snr_trh, Ppvalu, apr_dt_st1, apr_dt_st2, cpl_dist, ac_maxampint

integer(kind=4)    :: k, l
integer(kind=4)    :: nt, nn
integer(kind=4)    :: apr_dt_shift, ind_cma, ind_acma, c_at, ac_at, peak_ind1, peak_ind2, trough_ind1, trough_ind2
integer(kind=4)    :: nlen, nofts, nerr, no_peaks, no_troughs
integer(kind=4)    :: csl_length, int_length
integer(kind=4)    :: ac1, c1, ac2, c2, ac2_env, c2_env, nsmpls_per_hp
integer(kind=4)    :: sig, sig_ext, sig_int, sig_int_ext, send, send_int

real(kind=4)       , allocatable, dimension(:)      :: xcorr, xcorr_filt
real(kind=4)       , allocatable, dimension(:)      :: left_taper, right_taper
real(kind=8)       , allocatable, dimension(:)      :: peak_times, trough_times, peak_env, trough_env, c_sign_times, ac_sign_times
real(kind=8)       , allocatable, dimension(:)      :: int_c_sign, int_ac_sign
real(kind=8)       , allocatable, dimension(:,:)    :: peaks, troughs, c_sign, ac_sign


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! PARAMETERS AND DEPENDENT VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer           , parameter :: int64 = selected_int_kind( 16 )
real(8)           , parameter :: pi=atan2(0.,-1.)
complex(8)        , parameter :: zero=cmplx(0.,0.),one=cmplx(1.),zi=cmplx(0.,1.)                    !!!!
integer           , parameter :: passes=2,order=4
double precision  , parameter :: transition_bandwidth=0.0, attenuation=0.0
!logical           , parameter :: wr_resp_details=.false. !get details.
logical wr_resp_details
!########################################################
! variables needed
! xcorrpath : path to correlations
! xcorr_path = '/Users/davidnaranjohernandez/Dropbox/David/KAUST/Jupyter-notebooks/reykjanes/data/KEF_O15/total_013672.SAC'
! nn : number of samples (180000).

open(unit=10, file='params.txt')
read(10, *)
!params text file with:
!lf : low frq, hf: high freq. : filter, cpl_dist : couple distance in meters, ref_vel : reference velocity, dist_trh : Minimum separation between stations in terms of wavelength, snr_trh : signal-to-noise ratio threshold, noise_st=240 : start of the noise window (in seconds) used for the computation of the SNR, apr_dt_st1 and apr_dt_st2: estimated timing error of both stations in seconds. apr_dt_st1 is the first station of the directory name, dt_err : 0.004
read(10,'(A)')xcorr_path, station1, station2, dir_name
read(10, *)nn, lf, hf, cpl_dist, ref_vel, dist_trh, snr_trh, noise_st, apr_dt_st1, apr_dt_st2, dt_err, resp_details
close(10)
if (resp_details .eq. 'True') then
    wr_resp_details=.true.
else
    wr_resp_details=.false.
end if
nt = nn/2


! variables taken automatically
! BEGIN: time before the 0 in the trace
! DELTA: dt
! nlen: number of points in the trace
! size(xcorr,1) depends on nn
! nerr: error given from SAC.
allocate(xcorr(nn), xcorr_filt(nn))
!!! Time domain tapers
nofts=int(nint(0.1*nn),4) ! 0.1 is the fraction of total length used for tapering.
allocate(left_taper(nofts),right_taper(nofts))
call init_leftcos_taper(nofts,left_taper)
call init_rightcos_taper(nofts,right_taper)

call rsac1(trim(xcorr_path),xcorr(:),nlen,BEGIN,DELTA,size(xcorr,1),nerr)
dt = real(DELTA,8)
if(nerr .gt. 0) then    ! negative only a warning (e.g. nerr=-803 means number of samples is higher than size of array)
  write(0,'(A6,I2,A18,A)') 'Error ',nerr,' reading in file: ',trim(xcorr_path)
  stop
endif
 
dir_8='temp/'
cmd='mkdir '//trim(adjustl(dir_8))
call system(cmd)

!!! Velocities determining start and end of Signal window
v1=real(2*ref_vel, kind=8)
v2=real(ref_vel/2., kind=8)
cfs=(hf-lf)/2. ! central frequency.
nsmpls_per_hp=nint(1/cfs/2/dt) ! Number of samples per half-period
sig=nint(1/cfs/3/dt) ! Signal
sig_ext=2*sig
sig_int=nint(sig*(dt/dt_err))
sig_int_ext=2*sig_int

!!! Filter xcorrelation (after tapering)
xcorr_filt=xcorr
xcorr_filt(1:nofts)=xcorr_filt(1:nofts)*left_taper
xcorr_filt(nn-nofts+1:nn)=xcorr_filt(nn-nofts+1:nn)*right_taper
call xapiir(xcorr_filt,nn,SAC_BUTTERWORTH,transition_bandwidth,attenuation,order,SAC_BANDPASS,lf,hf,dt,passes)

call writenum(dist_trh,wl_trh_char,'F4.2')

!!! Cycle if stations are not sufficiently separated at this frequncy
min_wl=ref_vel/hf
if(cpl_dist/min_wl .lt. dist_trh)then
  stop 'Station couple does not exceed minimum separation.'
endif

call writenum(snr_trh,snr_trh_char,'F8.1')

apr_dt_cpl = 0
apr_dt_cpl=apr_dt_cpl + 2*apr_dt_st1
apr_dt_cpl=apr_dt_cpl - 2*apr_dt_st2

!!! Compute SNR accounting for the expected time shift
call SNR_with_shift(xcorr_filt,int(nt+1,kind=8),cpl_dist,apr_dt_cpl/2.,dt,noise_st,v1,v2,c_snr,ac_snr)

!!! Cycle if SNR's are not high enough at this frequncy
if(c_snr .lt. snr_trh .or. ac_snr .lt. snr_trh)then
  stop 'Station couple does not exceed SNR'
endif
    
!!! Compute envelope
apr_dt_shift=nint(apr_dt_cpl/2./dt)
c2_env=int(nt+1)+nint(3*cpl_dist/v2/dt)+apr_dt_shift !causal part
ac2_env=int(nt+1)-nint(3*cpl_dist/v2/dt)+apr_dt_shift ! acausal part
             
!!! Determine number of peaks and troughs
no_peaks=0;no_troughs=0
!write(6,*)ac2_env,c2_env
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
c1=int(nt+1)+nint(cpl_dist/v1/dt)+apr_dt_shift   ! Starting index of causal signal window
c2=int(nt+1)+nint(cpl_dist/v2/dt)+apr_dt_shift   ! Ending index causal signal window
ac1=int(nt+1)-nint(cpl_dist/v1/dt)+apr_dt_shift  ! Ending index of acausal signal window
ac2=int(nt+1)-nint(cpl_dist/v2/dt)+apr_dt_shift  ! Starting index acausal signal window

!!! Make sure that the signal windows are at least a period long
if( (c2-c1)*dt .lt. 1/cfs)then
  c2=c1+nint(1/cfs/dt)
  ac2=ac1-nint(1/cfs/dt)
endif

if(wr_resp_details)then

  !dir_9=trim(dir_8)//trim(snr_trh_char)//'_'//trim(wl_trh_char)//'/'
  dir_9=trim(dir_8)//trim(station1)//'_'//trim(station2)//'_'//trim(dir_name)//'/'
  cmd='mkdir '//trim(adjustl(dir_9))
  call system(cmd)

  !!! Write the begining of signal window
  fname=trim(dir_9)//'signal_window_indices'
  open(unit=30,file=fname)
  write(30,'(x,i8)')c1
  write(30,'(x,i8)')c2
  write(30,'(x,i8)')ac1
  write(30,'(x,i8)')ac2
  close(30)

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
allocate(c_sign(4,csl_length), ac_sign(4,csl_length), c_sign_times(csl_length) ,ac_sign_times(csl_length))
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

! For writting interpolated signal to SAC files.
!fname='/vardim/home/weemstra/tmp/c_sign.SAC' causal
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
result_shift=apr_dt_shift*2*dt+add_shift

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

deallocate(left_taper,right_taper,xcorr,xcorr_filt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(0,'(A)') '###############     PROGRAM COMPLETE      ################'
write(6,*)'Results saved in folder:', trim(dir_8)//trim(station1)//'_'//trim(station2)//'_'//trim(dir_name)//'/'
write(6,*)'Result shift: ',result_shift
write(6,*)'acausal signal from index: ',ac1
write(6,*)'acausal signal until index: ',ac2
write(6,*)'causal signal from index: ',c1
write(6,*)'causal signal until index: ',c2
write(6,*)'SNR causal wave',c_snr
write(6,*)'SNR acausal wave',ac_snr
end program TIMING_ERR_INV
