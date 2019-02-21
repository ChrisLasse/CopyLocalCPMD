program allocater_new
  use scratch_interface
  implicit none
  doubleprecision,pointer :: a(:),b(:)
  integer :: i,n,len(1),task
  character(255) :: tag
  task=-1
!  call segment_interface(task=task)
  write(*,*) associated(a)
  task=1
  len=10
  tag="test1"
!  call segment_interface(len=len,tag=tag,task=task,arrayout=a)
  call request_scratch(len,a,tag)
  do i=1,len(1)
     a(i)=i*i
  end do
  task=3
!  call segment_interface(len=len,tag=tag,task=task,arrayin=a) 
  call save_scratch(len,a,tag)
  task=1
  tag="test2"
!  call segment_interface(len=len,tag=tag,task=task,arrayout=b)
  call request_scratch(len,a,tag)
  task=4
  tag="test1"
!  call segment_interface(len=len,tag=tag,task=task,arrayout=a)
  call request_saved_scratch(len,a,tag)
  do i=1,len(1)
     write(*,*) a(i)
  end do
  
  write(*,*) associated(a)
end program allocater_new
