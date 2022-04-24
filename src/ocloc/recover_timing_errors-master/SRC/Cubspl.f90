subroutine cubspl ( tau, c, n, ibcbeg, ibcend )

!*****************************************************************************80
!
!! CUBSPL defines an interpolatory cubic spline.
!
!  Discussion:
!
!    A tridiagonal linear system for the unknown slopes S(I) of
!    F at TAU(I), I=1,..., N, is generated and then solved by Gauss
!    elimination, with S(I) ending up in C(2,I), for all I.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl de Boor
!
!  Reference:
!
!    Carl de Boor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TAU(N), the abscissas or X values of
!    the data points.  The entries of TAU are assumed to be
!    strictly increasing.
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N is
!    assumed to be at least 2.
!
!    Input/output, real ( kind = 8 ) C(4,N).
!    On input, if IBCBEG or IBCBEG is 1 or 2, then C(2,1)
!    or C(2,N) should have been set to the desired derivative
!    values, as described further under IBCBEG and IBCEND.
!    On output, C contains the polynomial coefficients of
!    the cubic interpolating spline with interior knots
!    TAU(2) through TAU(N-1).
!    In the interval interval (TAU(I), TAU(I+1)), the spline
!    F is given by
!      F(X) = 
!        C(1,I) + 
!        C(2,I) * H +
!        C(3,I) * H^2 / 2 + 
!        C(4,I) * H^3 / 6.
!    where H=X-TAU(I).  The routine PPVALU may be used to
!    evaluate F or its derivatives from TAU, C, L=N-1,
!    and K=4.
!
!    Input, integer ( kind = 4 ) IBCBEG, IBCEND, boundary condition indicators.
!    IBCBEG = 0 means no boundary condition at TAU(1) is given.
!    In this case, the "not-a-knot condition" is used.  That
!    is, the jump in the third derivative across TAU(2) is
!    forced to zero.  Thus the first and the second cubic
!    polynomial pieces are made to coincide.
!    IBCBEG = 1 means the slope at TAU(1) is to equal the
!    input value C(2,1).
!    IBCBEG = 2 means the second derivative at TAU(1) is
!    to equal C(2,1).
!    IBCEND = 0, 1, or 2 has analogous meaning concerning the
!    boundary condition at TAU(N), with the additional
!    information taken from C(2,N).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(4,n)
  real ( kind = 8 ) divdf1
  real ( kind = 8 ) divdf3
  real ( kind = 8 ) dtau
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  real ( kind = 8 ) tau(n)
!
!  C(3,*) and C(4,*) are used initially for temporary storage.
!
!  Store first differences of the TAU sequence in C(3,*).
!
!  Store first divided difference of data in C(4,*).
!
  do i = 2, n
    c(3,i) = tau(i) - tau(i-1)
  end do

  do i = 2, n 
    c(4,i) = ( c(1,i) - c(1,i-1) ) / ( tau(i) - tau(i-1) )
  end do
!
!  Construct the first equation from the boundary condition
!  at the left endpoint, of the form:
!
!    C(4,1) * S(1) + C(3,1) * S(2) = C(2,1)
!
!  IBCBEG = 0: Not-a-knot
!
  if ( ibcbeg == 0 ) then

    if ( n <= 2 ) then
      c(4,1) = 1.0D+00
      c(3,1) = 1.0D+00
      c(2,1) = 2.0D+00 * c(4,2)
      go to 120
    end if

    c(4,1) = c(3,3)
    c(3,1) = c(3,2) + c(3,3)
    c(2,1) = ( ( c(3,2) + 2.0D+00 * c(3,1) ) * c(4,2) * c(3,3) &
      + c(3,2)**2 * c(4,3) ) / c(3,1)
!
!  IBCBEG = 1: derivative specified.
!
  else if ( ibcbeg == 1 ) then

    c(4,1) = 1.0D+00
    c(3,1) = 0.0D+00

    if ( n == 2 ) then
      go to 120
    end if
!
!  Second derivative prescribed at left end.
!
  else

    c(4,1) = 2.0D+00
    c(3,1) = 1.0D+00
    c(2,1) = 3.0D+00 * c(4,2) - c(3,2) / 2.0D+00 * c(2,1)

    if ( n == 2 ) then
      go to 120
    end if

  end if
!
!  If there are interior knots, generate the corresponding
!  equations and carry out the forward pass of Gauss elimination,
!  after which the I-th equation reads:
!
!    C(4,I) * S(I) + C(3,I) * S(I+1) = C(2,I).
!
  do i = 2, n-1
    g = -c(3,i+1) / c(4,i-1)
    c(2,i) = g * c(2,i-1) + 3.0D+00 * ( c(3,i) * c(4,i+1) + c(3,i+1) * c(4,i) )
    c(4,i) = g * c(3,i-1) + 2.0D+00 * ( c(3,i) + c(3,i+1))
  end do
!
!  Construct the last equation from the second boundary condition, of
!  the form
!
!    -G * C(4,N-1) * S(N-1) + C(4,N) * S(N) = C(2,N)
!
!  If slope is prescribed at right end, one can go directly to
!  back-substitution, since the C array happens to be set up just
!  right for it at this point.
!
  if ( ibcend == 1 ) then
    go to 160
  end if

  if ( 1 < ibcend ) then
    go to 110
  end if
 
90    continue
!
!  Not-a-knot and 3 <= N, and either 3 < N or also not-a-knot
!  at left end point.
!
  if ( n /= 3 .or. ibcbeg /= 0 ) then
    g = c(3,n-1) + c(3,n)
    c(2,n) = ( ( c(3,n) + 2.0D+00 * g ) * c(4,n) * c(3,n-1) + c(3,n)**2 &
      * ( c(1,n-1) - c(1,n-2) ) / c(3,n-1) ) / g
    g = - g / c(4,n-1)
    c(4,n) = c(3,n-1)
    c(4,n) = c(4,n) + g * c(3,n-1)
    c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
    go to 160
  end if
!
!  N = 3 and not-a-knot also at left.
!
100   continue
 
  c(2,n) = 2.0D+00 * c(4,n)
  c(4,n) = 1.0D+00
  g = -1.0D+00 / c(4,n-1)
  c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
  c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
  go to 160
!
!  IBCEND = 2: Second derivative prescribed at right endpoint.
!
110   continue
 
  c(2,n) = 3.0D+00 * c(4,n) + c(3,n) / 2.0D+00 * c(2,n)
  c(4,n) = 2.0D+00
  g = -1.0D+00 / c(4,n-1)
  c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
  c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
  go to 160
!
!  N = 2.
!
120   continue
  
  if ( ibcend == 2  ) then

    c(2,n) = 3.0D+00 * c(4,n) + c(3,n) / 2.0D+00 * c(2,n)
    c(4,n) = 2.0D+00
    g = -1.0D+00 / c(4,n-1)
    c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
    c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
 
  else if ( ibcend == 0 .and. ibcbeg /= 0 ) then

    c(2,n) = 2.0D+00 * c(4,n)
    c(4,n) = 1.0D+00
    g = -1.0D+00 / c(4,n-1)
    c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
    c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)

  else if ( ibcend == 0 .and. ibcbeg == 0 ) then

    c(2,n) = c(4,n)

  end if
!
!  Back solve the upper triangular system 
!
!    C(4,I) * S(I) + C(3,I) * S(I+1) = B(I)
!
!  for the slopes C(2,I), given that S(N) is already known.
!
160   continue
 
  do i = n-1, 1, -1
    c(2,i) = ( c(2,i) - c(3,i) * c(2,i+1) ) / c(4,i)
  end do
!
!  Generate cubic coefficients in each interval, that is, the
!  derivatives at its left endpoint, from value and slope at its
!  endpoints.
!
  do i = 2, n
    dtau = c(3,i)
    divdf1 = ( c(1,i) - c(1,i-1) ) / dtau
    divdf3 = c(2,i-1) + c(2,i) - 2.0D+00 * divdf1
    c(3,i-1) = 2.0D+00 * ( divdf1 - c(2,i-1) - divdf3 ) / dtau
    c(4,i-1) = 6.0D+00 * divdf3 / dtau**2
  end do
 
  return
end
