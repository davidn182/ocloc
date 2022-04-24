module snr_routines

implicit none
public

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SNR(cc,n3,dist,d1,noise_st,v1,v2,c_snr,ac_snr)
 
  real, intent(in)  :: d1,v1,v2,dist,noise_st
  integer, intent(in)  :: n3
  real, intent(in),dimension(:) :: cc 
  real, intent(out) :: c_snr,ac_snr
  
  integer :: i,n,c_a1i,c_a2i,c_b1i,c_b2i,ac_a1i,ac_a2i,ac_b1i,ac_b2i
  real :: a1r,a2r,b1r,b2r,c_amp,ac_amp,c_sum,ac_sum,central_sum,tsum,mean,var,nvar,std
  a1r=dist/v1                ! Starting time of signal window 
  a2r=dist/v2                ! Ending time of signal window 
  b1r=noise_st               ! Starting time of noise window
  b2r=2.0*b1r                ! Ending time of noise window 
  !!!!! causal snr
  c_a1i=n3+nint(a1r/d1)        ! Starting index of causal signal window 
  c_a2i=n3+nint(a2r/d1)        ! Ending index of causal signal window
  c_b1i=n3+nint(b1r/d1)        ! Starting index of causal noise window 
  c_b2i=n3+nint(b2r/d1)        ! Ending index of acausal noise window
  ac_a1i=n3-nint(a1r/d1)        ! Ending index of acausal signal window 
  ac_a2i=n3-nint(a2r/d1)        ! Starting index of acausal signal window
  ac_b1i=n3-nint(b1r/d1)        ! Ending index of acausal noise window 
  ac_b2i=n3-nint(b2r/d1)        ! Starting index of acausal noise window

  c_amp=maxval(abs(cc(c_a1i:c_a2i))) 
  ac_amp=maxval(abs(cc(ac_a2i:ac_a1i))) 
  c_sum=sum(cc(c_b1i:c_b2i))
  ac_sum=sum(cc(ac_b2i:ac_b1i))
  central_sum=sum(cc(ac_a1i:c_a1i))
  tsum=c_sum+ac_sum+central_sum
  n=(c_a1i-ac_a1i+1)+(c_b2i-c_b1i+1)+(ac_b1i-ac_b2i+1)
  mean=tsum/real(n)
  var=0.;
  do i=c_b1i,c_b2i
    var=var+((cc(i)-mean)**2)
  enddo
  do i=ac_b2i,ac_b1i
    var=var+((cc(i)-mean)**2)
  enddo
  do i=ac_a1i,c_a1i
    var=var+((cc(i)-mean)**2)
  enddo
  nvar=var/real(n) 
  std=sqrt(nvar)
  c_snr=c_amp/std
  ac_snr=ac_amp/std
  !write(6,*)cc(1),cc(45000),cc(90000),cc(90001),cc(90002),c_amp,ac_amp,c_snr,ac_snr,std 
  end subroutine SNR



  subroutine SNR_with_shift(cc,n3,dist,apr_dt,d1,noise_st,v1,v2,c_snr,ac_snr)
 
  real(kind=8),    intent(in)                :: d1,apr_dt,dist,v1,v2,noise_st
  integer(kind=8), intent(in)                :: n3
  real,            intent(in),dimension(:)   :: cc 
  real(kind=8),    intent(out)               :: c_snr,ac_snr
  
  integer(kind=8) :: i,n,c_a1i,c_a2i,c_b1i,c_b2i,ac_a1i,ac_a2i,ac_b1i,ac_b2i
  real(kind=8)    :: a1r,a2r,b1r,b2r,c_amp,ac_amp,c_sum,ac_sum,tsum,mean,var,nvar,std
   
  a1r=dist/v1+apr_dt
  a2r=dist/v2+apr_dt               ! Ending time of causal signal window 
  b1r=noise_st+apr_dt              ! Starting time of causal noise window
  b2r=2.0*noise_st+apr_dt               ! Ending time of causal noise window 
  
  c_a1i=n3+nint(a1r/d1)        ! Starting index of causal signal window 
  c_a2i=n3+nint(a2r/d1)        ! Ending index of causal signal window
  c_b1i=n3+nint(b1r/d1)        ! Starting index of causal noise window 
  c_b2i=n3+nint(b2r/d1)        ! Ending index of acausal noise window

  a1r=dist/v1-apr_dt                ! Starting time of acausal signal window 
  a2r=dist/v2-apr_dt                ! Ending time of acausal signal window 
  b1r=noise_st-apr_dt               ! Starting time of acausal noise window
  b2r=2.0*noise_st-apr_dt                ! Ending time of acausal noise window 

  ac_a1i=n3-nint(a1r/d1)        ! Ending index of acausal signal window 
  ac_a2i=n3-nint(a2r/d1)        ! Starting index of acausal signal window
  ac_b1i=n3-nint(b1r/d1)        ! Ending index of acausal noise window 
  ac_b2i=n3-nint(b2r/d1)        ! Starting index of acausal noise window

  c_amp=real(maxval(abs(cc(c_a1i:c_a2i))),kind=8)
  ac_amp=real(maxval(abs(cc(ac_a2i:ac_a1i))),kind=8)
  c_sum=real( sum( cc(c_b1i:c_b2i) ),kind=8)
  ac_sum=real( sum( cc(ac_b2i:ac_b1i) ),kind=8)
  tsum=c_sum+ac_sum
  n=(c_b2i-c_b1i+1)+(ac_b1i-ac_b2i+1)
  mean=tsum/real(n)
  var=0.
  do i=c_b1i,c_b2i
    var=var+((real(cc(i),kind=8)-mean)**2)
  enddo
  do i=ac_b2i,ac_b1i
    var=var+((real(cc(i),kind=8)-mean)**2)
  enddo
  nvar=var/real(n) 
  std=sqrt(nvar)
  c_snr=c_amp/std
  ac_snr=ac_amp/std
  end subroutine SNR_with_shift

end module snr_routines
