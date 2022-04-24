module ls_solutions

 
implicit none   ! If you do implicit none here, you don't need to do it in the functions and subroutines 

integer(kind=4)                                 :: i, j, k, smlsiz, lwork, liwork, nlvl, ldb, info 
integer(kind=4)                                 :: rank
real(kind=8)                                    :: rcond
real(kind=8)      , allocatable, dimension(:,:) :: A_dum, AtA
real(kind=8)      , allocatable, dimension(:)   :: Tobs_dum, Tobs_lsqr_dum, S, work, Tins_lsqr_dum
integer(kind=4)   , allocatable, dimension(:)   :: iwork
real(kind=8)      , allocatable, dimension(:,:) :: Tobs_lsqr, Tins_lsqr, Tins_var


private

public   ::  simple_ls, weighted_ls, nzm_weighted_ls, alt_nzm_weighted_ls


contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine simple_ls(nrows, ncols, A, Tobs, Tins_lsqr, Tobs_lsqr, Tins_var, R0)

    integer(kind=4), intent(in)  :: nrows, ncols
    real(kind=8),    intent(in), dimension(nrows,ncols) :: A 
    real(kind=8),    intent(in), dimension(nrows)       :: Tobs 

    real(kind=8),    intent(inout) :: R0
    real(kind=8),    intent(out), dimension(ncols)     :: Tins_lsqr, Tins_var
    real(kind=8),    intent(out), dimension(nrows)     :: Tobs_lsqr

    !!! allocate local arrays
    call init_subr(nrows,ncols)
    
    !!! Initialize design matrix and data vector
    Tobs_dum(1:nrows) = Tobs
    A_dum(1:nrows,1:ncols) = A

    !!! Solve the system of equation
    S=0. ; work=0.; iwork=0; rcond=0.001; info=0
    call dgelsd( nrows, ncols, 1, A_dum, max(1,nrows), Tobs_dum, ldb, S, rcond, rank, work, lwork, iwork, info )

    if(rank .lt. ncols)then

      write(6,*)'nrows = ',nrows,'; number of cols = ',ncols
      call deinit_subr()

      return
    endif

    !!! Save the solutions
    Tins_lsqr(:)=Tobs_dum(1:ncols)
    Tins_lsqr_dum=Tins_lsqr

    !!! Compute least square data vector (how would the data look with the best solution)
    do i=1,nrows
      Tobs_lsqr(i)=sum(A(i,:)*Tins_lsqr(:))
    enddo

    !!! Reinitialize design matrix and data vector
    Tobs_dum = Tobs
    A_dum = A

    !!! Computation of covariance matrix for purpose of estimation of model variance
    call cov_matrix(nrows, ncols, Tins_var, R0)

    call deinit_subr()
 
    return

  end subroutine simple_ls

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine weighted_ls(nrows, ncols, A, Tobs, weights, Tins_lsqr, Tobs_lsqr, Tins_var, R0)

    integer(kind=4), intent(in)  :: nrows, ncols
    real(kind=8),    intent(in), dimension(nrows,ncols) :: A 
    real(kind=8),    intent(in), dimension(nrows)       :: Tobs, weights

    real(kind=8),    intent(inout) :: R0
    real(kind=8),    intent(out), dimension(ncols)        :: Tins_lsqr, Tins_var
    real(kind=8),    intent(out), dimension(nrows)        :: Tobs_lsqr

    !!! allocate local arrays
    call init_subr(nrows,ncols)

    !!! Initialize design matrix and data vector
    do i=1,nrows
      Tobs_dum(i) =Tobs(i)*weights(i)
      A_dum(i,:)  =A(i,:)*weights(i)
    enddo

    !!! Solve the system of equation
    S=0. ; work=0.; iwork=0; rcond=0.001; info=0
    call dgelsd( nrows, ncols, 1, A_dum, max(1,nrows), Tobs_dum, ldb, S, rcond, rank, work, lwork, iwork, info )
    write(6,*)'rank = ',rank,'; number of unknows = ',ncols

    if(rank .lt. ncols)then

      call deinit_subr()

      return
    endif

    !!! Save the solutions
    Tins_lsqr(:)=Tobs_dum(1:ncols)
    Tins_lsqr_dum=Tins_lsqr

    !!! Compute least square data vector (how would the data look with the best solution)
    do i=1,nrows
      Tobs_lsqr(i)=sum(A(i,:)*Tins_lsqr(:))
    enddo

    !!! Reinitialize design matrix and data vector
    do i=1,nrows
      A_dum(i,:)  =A(i,:)*weights(i)
      Tobs_dum(i) =Tobs(i)*weights(i)
    enddo

    call cov_matrix(nrows, ncols, Tins_var, R0)

    call deinit_subr()
 
    return

  end subroutine weighted_ls

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nzm_weighted_ls(nrows, ncols, A, Tobs, weights, Tins_lsqr, Tobs_lsqr, Tins_var, R0)

    integer(kind=4), intent(in)  :: nrows, ncols
    real(kind=8),    intent(in), dimension(nrows,ncols) :: A 
    real(kind=8),    intent(in), dimension(nrows)       :: Tobs, weights

    real(kind=8),    intent(inout) :: R0
    real(kind=8),    intent(out), dimension(ncols)     :: Tins_lsqr, Tins_var
    real(kind=8),    intent(out), dimension(nrows)     :: Tobs_lsqr
    real(kind=8),    allocatable, dimension(:)         :: Tins_var_dum

    !!! allocate local arrays
    call init_subr(nrows,ncols+1)

    !!! Initialize design matrix and data vector
    A_dum(1:nrows,1:ncols) = A
    do i=1,nrows
      A_dum(i,ncols+1) = 1./weights(i)
    enddo
    do i=1,nrows
      Tobs_dum(i) =Tobs(i)*weights(i)
      A_dum(i,:)  =A_dum(i,:)*weights(i)
    enddo

    !!! Solve the system of equation
    S=0. ; work=0.; iwork=0; rcond=0.001; info=0
    call dgelsd( nrows, ncols+1, 1, A_dum, max(1,nrows), Tobs_dum, ldb, S, rcond, rank, work, lwork, iwork, info )
    write(6,*)'rank = ',rank,'; number of unknows = ',ncols+1

    if(rank .lt. (ncols+1))then

      call deinit_subr()

      return
    endif

    !!! Save the solutions
    Tins_lsqr(:)=Tobs_dum(1:ncols)  !!! the element containing the mean of the noise is not send back
    Tins_lsqr_dum(:)=Tobs_dum(1:ncols+1)
    write(6,*)'The mean of value of Nazi = ',Tobs_dum(ncols+1)

    !!! Compute least square data vector (how would the data look with the best solution)
    do i=1,nrows
      Tobs_lsqr(i)=sum(A(i,1:ncols)*Tins_lsqr(:))
    enddo

    !!! Reinitialize design matrix and data vector
    A_dum(1:nrows,1:ncols) = A
    do i=1,nrows
      A_dum(i,ncols+1) = 1./weights(i)
    enddo
    do i=1,nrows
      Tobs_dum(i) =Tobs(i)*weights(i)
      A_dum(i,:)  =A_dum(i,:)*weights(i)
    enddo

    allocate(Tins_var_dum(ncols+1))
    call cov_matrix(nrows, ncols+1, Tins_var_dum, R0)
    Tins_var(:)=Tins_var_dum(1:ncols)
 
    deallocate(Tins_var_dum)
    call deinit_subr()
   
    return

  end subroutine nzm_weighted_ls

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alt_nzm_weighted_ls(nrows, ncols, A, Tobs, weights, Tins_lsqr, Tobs_lsqr, Tins_var, R0)

    integer(kind=4), intent(in)  :: nrows, ncols
    real(kind=8),    intent(in), dimension(nrows,ncols) :: A 
    real(kind=8),    intent(in), dimension(nrows)       :: Tobs, weights

    real(kind=8),    intent(inout) :: R0
    real(kind=8),    intent(out), dimension(ncols)     :: Tins_lsqr, Tins_var
    real(kind=8),    intent(out), dimension(nrows)     :: Tobs_lsqr
    real(kind=8),    allocatable, dimension(:)         :: Tins_var_dum

    !!! allocate local arrays
    call init_subr(nrows,ncols+1)

    !!! Initialize design matrix and data vector
    A_dum(1:nrows,1:ncols) = A
    do i=1,nrows
      A_dum(i,ncols+1) = 1.              !!! Note, THIS IS DIFFERENT WITH RESPECT TO PREVIOUS SUBROUTINE
    enddo
    do i=1,nrows
      Tobs_dum(i) =Tobs(i)*weights(i)
      A_dum(i,:)  =A_dum(i,:)*weights(i)
    enddo

    !!! Solve the system of equation
    S=0. ; work=0.; iwork=0; rcond=0.001; info=0
    call dgelsd( nrows, ncols+1, 1, A_dum, max(1,nrows), Tobs_dum, ldb, S, rcond, rank, work, lwork, iwork, info )
    write(6,*)'rank = ',rank,'; number of unknows = ',ncols+1

    if(rank .lt. (ncols+1))then

      call deinit_subr()

      return
    endif

    !!! Save the solutions
    Tins_lsqr(:)=Tobs_dum(1:ncols)  !!! the element containing the mean of the noise is not send back
    Tins_lsqr_dum(:)=Tobs_dum(1:ncols+1)
    write(6,*)'The mean of value of Nazi = ',Tobs_dum(ncols+1)

    !!! Compute least square data vector (how would the data look with the best solution)
    do i=1,nrows
      Tobs_lsqr(i)=sum(A(i,1:ncols)*Tins_lsqr(:))
    enddo

    !!! Reinitialize design matrix and data vector
    A_dum(1:nrows,1:ncols) = A
    do i=1,nrows
      A_dum(i,ncols+1) = 1.              !!! Note, THIS IS DIFFERENT WITH RESPECT TO PREVIOUS SUBROUTINE
    enddo
    do i=1,nrows
      Tobs_dum(i) =Tobs(i)*weights(i)
      A_dum(i,:)  =A_dum(i,:)*weights(i)
    enddo

    allocate(Tins_var_dum(ncols+1))
    call cov_matrix(nrows, ncols+1, Tins_var_dum, R0)
    Tins_var(:)=Tins_var_dum(1:ncols)
 
    deallocate(Tins_var_dum)
    call deinit_subr()
   
    return

  end subroutine alt_nzm_weighted_ls
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cov_matrix(nrows, ncols, Tins_var, R0)

    integer(kind=4), intent(in)    :: nrows, ncols
    real(kind=8),    intent(inout) :: R0
    real(kind=8),    intent(out), dimension(ncols)        :: Tins_var

    do i=1,ncols
      do j=1,ncols
        AtA(i,j)=sum(A_dum(:,i)*A_dum(:,j))
      enddo
    enddo

    call dpotrf('L',ncols,AtA,ncols,info)
    call dpotri('L',ncols,AtA,ncols,info)

    !!! Compute least square data vecter
    do i=1,nrows
      Tobs_lsqr_dum(i)=sum(A_dum(i,:)*Tins_lsqr_dum(:))
    enddo

    !!! Computation of minimum error for purpose of estimation of model variance
    R0=sum((Tobs_dum-Tobs_lsqr_dum)**2)

    !!! Computation of variance estimates
    do i=1,ncols
      Tins_var(i)=AtA(i,i)* ( R0 / (nrows-ncols) )
    enddo

    return

  end subroutine cov_matrix

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_subr(nrows, ncols)

    integer(kind=4), intent(in)  :: nrows, ncols

    allocate(A_dum(max(1,nrows),ncols))
    allocate(Tobs_dum(max(nrows,ncols)))
    allocate(S(min(ncols,nrows)) )
    allocate(AtA(ncols,ncols))
    allocate(Tobs_lsqr_dum(nrows));Tobs_lsqr_dum=0.
    allocate(Tins_lsqr_dum(ncols))

    smlsiz=50     !! 25 is probably the approximate value. 50 is on the safe side.
    nlvl = max(0, nint(log10( real(min( nrows, ncols ),kind=8) / ( smlsiz+1 ) )/log10(2.)) + 1 )
    lwork=12*ncols + 2*ncols*smlsiz + 8*ncols*nlvl + ncols*1 + (smlsiz+1)**2
    liwork=( 3 * min(nrows,ncols) * nlvl ) + ( 11 * min(nrows,ncols) )
    allocate(work(lwork),iwork(liwork))
    ldb=max(nrows,ncols)
  
    return

  end subroutine init_subr

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine deinit_subr()

      deallocate(A_dum, Tobs_dum, S, AtA, work, iwork, Tobs_lsqr_dum, Tins_lsqr_dum)
      !if(allocated(A_dum)) write(6,*) 'A_dum is now being deallocated'
      !deallocate(A_dum)
      !if(allocated(Tobs_dum)) write(6,*) 'Tobs_dum is now being deallocated'
      !deallocate(Tobs_dum)
      !if(allocated(S)) write(6,*) 'S is now being deallocated'
      !deallocate(S)
      !if(allocated(AtA)) write(6,*) 'AtA is now being deallocated'
      !deallocate(AtA)
      !if(allocated(work)) write(6,*) 'work is now being deallocated'
      !deallocate(work)
      !if(allocated(iwork)) write(6,*) 'iwork is now being deallocated'
      !deallocate(iwork)
      !if(allocated(Tobs_lsqr_dum)) write(6,*) 'Tobs_lsqr_dum is now being deallocated'
      !deallocate(Tobs_lsqr_dum)
      !if(allocated(Tins_lsqr_dum)) write(6,*) 'Tins_lsqr_dum is now being deallocated'
      !deallocate(Tins_lsqr_dum)

      return

  end subroutine deinit_subr

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ls_solutions
