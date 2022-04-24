! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

module array_operations

implicit none
public :: qsortC,int_qsortC,remove_dups
private :: partition,int_partition

contains

recursive subroutine int_qsortC(A)
  integer(8), intent(in out), dimension(:) :: A
  integer :: iq
  if(size(A) > 1) then
     call int_partition(A, iq)
     call int_qsortC(A(:iq-1))
     call int_qsortC(A(iq:))
  endif
end subroutine int_qsortC

recursive subroutine qsortC(A)
  real, intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call partition(A, iq)
     call qsortC(A(:iq-1))
     call qsortC(A(iq:))
  endif
end subroutine qsortC

subroutine int_partition(A, marker)
  integer(8), intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  integer(8) :: temp
  integer(8) :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do
end subroutine int_partition

subroutine partition(A, marker)
  real, intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  real :: temp
  real :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do
end subroutine partition
subroutine remove_dups(A,B,k)
  integer, intent(in), dimension(:) :: A
  integer, intent(out), dimension(:) :: B
  integer, intent(out) :: k
  integer :: i, j

!call remove_dups(ws,ws_dummy,nw) 
  !example = [1, 2, 3, 2, 2, 4, 5, 5, 4, 6, 6, 5]
  k = 1
  B(1) = A(1)
  outer: do i=2,size(A)
     do j=1,k
        if (B(j) == A(i)) then
           ! Found a match so start looking again
           cycle outer
        end if
     end do
     ! No match found so add it to the output
     k = k + 1
     B(k) = A(i)
  end do outer
  !write(*,advance='no',fmt='(a,i0,a)') 'Unique list has ',k,' elements: '
  !write(*,*) res(1:k)
end subroutine remove_dups

end module array_operations
