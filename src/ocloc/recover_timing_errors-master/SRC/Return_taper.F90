!APPLY A COSINE TAPER
module return_taper
  implicit none
  contains

!----------------------------------------------------------------
  subroutine init_cos_taper(n1,d1,minf0,minf1,maxf0,maxf1,tapered_window)
    integer :: i1
    integer :: i1_x,i0_x,i1_n,i0_n
    real    :: pi2
    integer, intent(in) :: n1
    real, intent(in) :: d1,minf0,minf1,maxf0,maxf1
    real, intent(out),dimension(:) :: tapered_window 

    pi2=atan(1.)*2. !pi/2.
    i0_n = min(n1,max(1,nint((minf0)/d1+1.)))   ! 0
    i1_n = min(n1,max(1,nint((minf1)/d1+1.)))
    i0_x = min(n1,max(1,nint((maxf0)/d1+1.)))
    i1_x = min(n1,max(1,nint((maxf1)/d1+1.)))
    !write(0,'(A,I5)')'i0_n = ',i0_n
    !write(0,'(A,I5)')'i1_n = ',i1_n
    !write(0,'(A,I5)')'i0_x = ',i0_x
    !write(0,'(A,I5)')'i1_x = ',i1_x
    if ( i0_n<1 .or. i0_n>n1 ) then
      write(0,'(A)') 'invalid min0 to init_cos_taper'
      stop
    end if
    if ( i1_n<1 .or. i1_n>n1 ) then
      write(0,'(A)') 'invalid min1 to init_cos_taper'
      stop
    end if
    if ( i0_x<1 .or. i0_x>n1 ) then
      write(0,'(A)') 'invalid max0 to init_cos_taper'
      stop
    end if
    if ( i1_x<1 .or. i1_x>n1 ) then
      write(0,'(A)') 'invalid max1 to init_cos_taper'
      stop
    end if
    if ( i0_n > i1_n) then
      write(0,'(A)') 'in0 > min1'
      stop
    end if
    if ( i0_x < i1_x) then
      write(0,'(A)') 'max0 < max1'
      stop
    end if
    do i1=i0_n,i1_n 
      tapered_window(i1)=(cos(pi2-(real(i1-i0_n)/real(i1_n-i0_n))*pi2))
    end do 
    do i1=i1_x,i0_x 
      tapered_window(i1)=(cos((real(i1-i1_x)/real(i0_x-i1_x))*pi2) )
    end do
    tapered_window(1:i0_n)=0.
    tapered_window(i0_x:n1)=0.
    tapered_window(i1_n:i1_x)=1.
  end subroutine init_cos_taper

!----------------------------------------------------------------
  subroutine init_dcos_taper(n1,d1,minf0,minf1,maxf0,maxf1,tapered_window)
    integer :: i1
    integer :: i1_x,i0_x,i1_n,i0_n
    real    :: pi2,pi
    integer, intent(in) :: n1
    real, intent(in) :: d1,minf0,minf1,maxf0,maxf1
    real, intent(out),dimension(:) :: tapered_window 
   
    pi2=atan(1.)*2. !pi/2.
    pi=atan2(0.,-1.)  !pi
    i0_n = min(n1,max(1,nint((minf0)/d1+1.)))   ! 0
    i1_n = min(n1,max(1,nint((minf1)/d1+1.)))
    i0_x = min(n1,max(1,nint((maxf0)/d1+1.)))
    i1_x = min(n1,max(1,nint((maxf1)/d1+1.)))
    !write(0,'(A,I5)')'i0_n = ',i0_n
    !write(0,'(A,I5)')'i1_n = ',i1_n
    !write(0,'(A,I5)')'i0_x = ',i0_x
    !write(0,'(A,I5)')'i1_x = ',i1_x
    if ( i0_n<1 .or. i0_n>n1 ) then
      write(0,'(A)') 'invalid min0 to init_cos_taper'
      stop
    end if
    if ( i1_n<1 .or. i1_n>n1 ) then
      write(0,'(A)') 'invalid min1 to init_cos_taper'
      stop
    end if
    if ( i0_x<1 .or. i0_x>n1 ) then
      write(0,'(A)') 'invalid max0 to init_cos_taper'
      stop
    end if
    if ( i1_x<1 .or. i1_x>n1 ) then
      write(0,'(A)') 'invalid max1 to init_cos_taper'
      stop
    end if
    if ( i0_n > i1_n) then
      write(0,'(A)') 'in0 > min1'
      stop
    end if
    if ( i0_x < i1_x) then
      write(0,'(A)') 'max0 < max1'
      stop
    end if
    do i1=i0_n,i1_n 
      tapered_window(i1)=(cos(pi-(real(i1-i0_n)/real(i1_n-i0_n))*pi)+1.)/2.
    end do 
    do i1=i1_x,i0_x 
      tapered_window(i1)=(cos((real(i1-i1_x)/real(i0_x-i1_x))*pi)+1.)/2
    end do
    tapered_window(1:i0_n)=0.
    tapered_window(i0_x:n1)=0.
    tapered_window(i1_n:i1_x)=1.
  end subroutine init_dcos_taper
!----------------------------------------------------------------
  subroutine init_leftcos_taper(n,tapered_window)
    integer :: i
    real    :: pi2,pi
    integer, intent(in) :: n
    real, intent(out),dimension(:) :: tapered_window 
   
    pi2=atan(1.)*2. !pi/2.
    pi=atan2(0.,-1.)  !pi
    do i=1,n
      tapered_window(i)=(cos(pi-(real(i-1)/real(n-1))*pi)+1.)/2.
    end do 
  end subroutine init_leftcos_taper
!----------------------------------------------------------------
  subroutine init_rightcos_taper(n,tapered_window)
    integer :: i
    real    :: pi2,pi
    integer, intent(in) :: n
    real, intent(out),dimension(:) :: tapered_window 
   
    pi2=atan(1.)*2. !pi/2.
    pi=atan2(0.,-1.)  !pi
    do i=1,n 
      tapered_window(i)=(cos((real(i-1)/real(n-1))*pi)+1.)/2.
    end do
  end subroutine init_rightcos_taper
end module 
