module file_operations

use strings
implicit none

contains

subroutine nlines(path,nl,taskid,ios)
  character(len=*), intent(in) :: path
  integer, intent(in)  :: taskid
  integer, intent(out) :: nl,ios
  logical :: file_spec
  character(len=3)                 :: taskid_char
  !character(len=20)                :: fname
  character(len=50)                :: fname

  inquire(file=path, exist=file_spec)
  if(file_spec)then
    ios=0
    call writenum(taskid,taskid_char,'I3.3')
    fname='/vardim/home/weemstra/tmp/no_of_lines_'//trim(taskid_char)
    call system('wc -l '//trim(adjustl(path))//' > '//fname)
    open(unit=20,file=fname)
      read(20,*)nl
    close(20)
    call system('rm '//trim(fname))
  else
    ios=1
  endif
  return
end subroutine nlines

subroutine nfiles(path,nf,taskid,ios)
  character(len=*), intent(in) :: path
  integer, intent(in)  :: taskid
  integer, intent(out) :: nf,ios
  logical :: file_spec
  character(len=3)                 :: taskid_char
  !character(len=20)                :: fname
  character(len=50)                :: fname

  inquire(file=path, exist=file_spec)
  if(file_spec)then
    ios=0
    call writenum(taskid,taskid_char,'I3.3')
    fname='/vardim/home/weemstra/tmp/no_of_files_'//trim(taskid_char)
    call system('ls '//trim(adjustl(path))//'| wc -l > '//fname)
    open(unit=20,file=fname)
    read(20,*)nf
    close(20)
    call system('rm '//trim(fname))
  else
    ios=1
  endif
  return
end subroutine nfiles

subroutine retr_files(path,nf,farray,taskid,ios)
  character(len=*), intent(in) :: path
  integer, intent(in)  :: nf, taskid
  integer, intent(out) :: ios
  logical :: file_spec
  character(len=*), intent(out), dimension(:) :: farray
  integer :: i
  character(len=3)                 :: taskid_char
  !character(len=18)                :: fname
  character(len=50)                :: fname

  inquire(file=path, exist=file_spec)
  if(file_spec)then
    ios=0
    call writenum(taskid,taskid_char,'I3.3')
    fname='/vardim/home/weemstra/tmp/input_'//trim(taskid_char)//'.txt'
    call system('ls '//trim(adjustl(path))//' > '//fname)
    open(unit=20,file=fname)
    do i=1,nf
      read(20,*)farray(i)
    enddo
    close(20)
    call system('rm '//fname)
  else
    ios=1
  endif
  return
end subroutine retr_files

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

end module file_operations
