module scratch_interface
  use scratch_data, only : dp
  use, intrinsic :: iso_c_binding
  use scratch_segment_managment_interface
  implicit none
  private
  public :: request_scratch
  public :: save_scratch
  public :: free_scratch
  public :: request_saved_scratch

  interface request_scratch
     module procedure :: request_scratch_r1
     module procedure :: request_scratch_r2
     module procedure :: request_scratch_r3
     module procedure :: request_scratch_r4
     module procedure :: request_scratch_r5
     module procedure :: request_scratch_r6
     module procedure :: request_scratch_c1
     module procedure :: request_scratch_c2
     module procedure :: request_scratch_c3
     module procedure :: request_scratch_c4
     module procedure :: request_scratch_c5
     module procedure :: request_scratch_c6
  end interface request_scratch

  interface save_scratch
     module procedure :: save_scratch_r1
     module procedure :: save_scratch_r2
     module procedure :: save_scratch_r3
     module procedure :: save_scratch_r4
     module procedure :: save_scratch_r5
     module procedure :: save_scratch_r6
     module procedure :: save_scratch_c1
     module procedure :: save_scratch_c2
     module procedure :: save_scratch_c3
     module procedure :: save_scratch_c4
     module procedure :: save_scratch_c5
     module procedure :: save_scratch_c6
  end interface save_scratch

  interface free_scratch
     module procedure :: free_scratch_r1
     module procedure :: free_scratch_r2
     module procedure :: free_scratch_r3
     module procedure :: free_scratch_r4
     module procedure :: free_scratch_r5
     module procedure :: free_scratch_r6
     module procedure :: free_scratch_c1
     module procedure :: free_scratch_c2
     module procedure :: free_scratch_c3
     module procedure :: free_scratch_c4
     module procedure :: free_scratch_c5
     module procedure :: free_scratch_c6
  end interface free_scratch

  interface request_saved_scratch
     module procedure :: request_saved_scratch_r1
     module procedure :: request_saved_scratch_r2
     module procedure :: request_saved_scratch_r3
     module procedure :: request_saved_scratch_r4
     module procedure :: request_saved_scratch_r5
     module procedure :: request_saved_scratch_r6
     module procedure :: request_saved_scratch_c1
     module procedure :: request_saved_scratch_c2
     module procedure :: request_saved_scratch_c3
     module procedure :: request_saved_scratch_c4
     module procedure :: request_saved_scratch_c5
     module procedure :: request_saved_scratch_c6
  end interface request_saved_scratch

  contains

  subroutine request_scratch_r1 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(1)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(out),contiguous :: ptr(:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    call segment_interface(len=len,tag=tag,task=1,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_r(1)/))
    return
  end subroutine request_scratch_r1

  subroutine request_saved_scratch_r1 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(1)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(out),contiguous :: ptr(:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    call segment_interface(len=len,tag=tag,task=4,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_r(1)/))
    return
  end subroutine request_saved_scratch_r1

  subroutine free_scratch_r1 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(1)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(in),contiguous :: ptr(:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=2,arrayin=tmp)
    return
  end subroutine free_scratch_r1

  subroutine save_scratch_r1 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(1)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(in),contiguous :: ptr(:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=3,arrayin=tmp)
    return
  end subroutine save_scratch_r1
  subroutine request_scratch_r2 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(2)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(out),contiguous :: ptr(:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    call segment_interface(len=len,tag=tag,task=1,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_r(1),il_r(2)/))
    return
  end subroutine request_scratch_r2

  subroutine request_saved_scratch_r2 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(2)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(out),contiguous :: ptr(:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    call segment_interface(len=len,tag=tag,task=4,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_r(1),il_r(2)/))
    return
  end subroutine request_saved_scratch_r2

  subroutine free_scratch_r2 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(2)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(in),contiguous :: ptr(:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=2,arrayin=tmp)
    return
  end subroutine free_scratch_r2

  subroutine save_scratch_r2 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(2)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(in),contiguous :: ptr(:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=3,arrayin=tmp)
    return
  end subroutine save_scratch_r2
  subroutine request_scratch_r3 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(3)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(out),contiguous :: ptr(:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    call segment_interface(len=len,tag=tag,task=1,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_r(1),il_r(2),il_r(3)/))
    return
  end subroutine request_scratch_r3

  subroutine request_saved_scratch_r3 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(3)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(out),contiguous :: ptr(:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    call segment_interface(len=len,tag=tag,task=4,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_r(1),il_r(2),il_r(3)/))
    return
  end subroutine request_saved_scratch_r3

  subroutine free_scratch_r3 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(3)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(in),contiguous :: ptr(:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=2,arrayin=tmp)
    return
  end subroutine free_scratch_r3

  subroutine save_scratch_r3 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(3)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(in),contiguous :: ptr(:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=3,arrayin=tmp)
    return
  end subroutine save_scratch_r3
  subroutine request_scratch_r4 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(4)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(out),contiguous :: ptr(:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    call segment_interface(len=len,tag=tag,task=1,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_r(1),il_r(2),il_r(3),il_r(4)/))
    return
  end subroutine request_scratch_r4

  subroutine request_saved_scratch_r4 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(4)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(out),contiguous :: ptr(:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    call segment_interface(len=len,tag=tag,task=4,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_r(1),il_r(2),il_r(3),il_r(4)/))
    return
  end subroutine request_saved_scratch_r4

  subroutine free_scratch_r4 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(4)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(in),contiguous :: ptr(:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=2,arrayin=tmp)
    return
  end subroutine free_scratch_r4

  subroutine save_scratch_r4 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(4)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(in),contiguous :: ptr(:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=3,arrayin=tmp)
    return
  end subroutine save_scratch_r4
  subroutine request_scratch_r5 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(5)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(out),contiguous :: ptr(:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    call segment_interface(len=len,tag=tag,task=1,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_r(1),il_r(2),il_r(3),il_r(4),il_r(5)/))
    return
  end subroutine request_scratch_r5

  subroutine request_saved_scratch_r5 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(5)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(out),contiguous :: ptr(:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    call segment_interface(len=len,tag=tag,task=4,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_r(1),il_r(2),il_r(3),il_r(4),il_r(5)/))
    return
  end subroutine request_saved_scratch_r5

  subroutine free_scratch_r5 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(5)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(in),contiguous :: ptr(:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=2,arrayin=tmp)
    return
  end subroutine free_scratch_r5

  subroutine save_scratch_r5 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(5)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(in),contiguous :: ptr(:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=3,arrayin=tmp)
    return
  end subroutine save_scratch_r5
  subroutine request_scratch_r6 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(6)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(out),contiguous :: ptr(:,:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    call segment_interface(len=len,tag=tag,task=1,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_r(1),il_r(2),il_r(3),il_r(4),il_r(5),il_r(6)/))
    return
  end subroutine request_scratch_r6

  subroutine request_saved_scratch_r6 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(6)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(out),contiguous :: ptr(:,:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    call segment_interface(len=len,tag=tag,task=4,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_r(1),il_r(2),il_r(3),il_r(4),il_r(5),il_r(6)/))
    return
  end subroutine request_saved_scratch_r6

  subroutine free_scratch_r6 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(6)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(in),contiguous :: ptr(:,:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=2,arrayin=tmp)
    return
  end subroutine free_scratch_r6

  subroutine save_scratch_r6 (il_r,ptr,tagin)
    integer,intent(in) :: il_r(6)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    REAL(dp),pointer,intent(in),contiguous :: ptr(:,:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_r)
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=3,arrayin=tmp)
    return
  end subroutine save_scratch_r6
  subroutine request_scratch_c1 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(1)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(out),contiguous :: ptr(:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    call segment_interface(len=len,tag=tag,task=1,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_c/))
    return
  end subroutine request_scratch_c1

  subroutine request_saved_scratch_c1 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(1)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(out),contiguous :: ptr(:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    call segment_interface(len=len,tag=tag,task=4,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_c/))
    return
  end subroutine request_saved_scratch_c1

  subroutine free_scratch_c1 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(1)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(in),contiguous :: ptr(:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=2,arrayin=tmp)
    return
  end subroutine free_scratch_c1

  subroutine save_scratch_c1 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(1)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(in),contiguous :: ptr(:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=3,arrayin=tmp)
    return
  end subroutine save_scratch_c1
  subroutine request_scratch_c2 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(2)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(out),contiguous :: ptr(:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    call segment_interface(len=len,tag=tag,task=1,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_c(1),il_c(2)/))
    return
  end subroutine request_scratch_c2

  subroutine request_saved_scratch_c2 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(2)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(out),contiguous :: ptr(:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    call segment_interface(len=len,tag=tag,task=4,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_c(1),il_c(2)/))
    return
  end subroutine request_saved_scratch_c2

  subroutine free_scratch_c2 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(2)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(in),contiguous :: ptr(:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=2,arrayin=tmp)
    return
  end subroutine free_scratch_c2

  subroutine save_scratch_c2 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(2)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(in),contiguous :: ptr(:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=3,arrayin=tmp)
    return
  end subroutine save_scratch_c2
  subroutine request_scratch_c3 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(3)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(out),contiguous :: ptr(:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    call segment_interface(len=len,tag=tag,task=1,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_c(1),il_c(2),il_c(3)/))
    return
  end subroutine request_scratch_c3

  subroutine request_saved_scratch_c3 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(3)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(out),contiguous :: ptr(:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    call segment_interface(len=len,tag=tag,task=4,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_c(1),il_c(2),il_c(3)/))
    return
  end subroutine request_saved_scratch_c3

  subroutine free_scratch_c3 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(3)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(in),contiguous :: ptr(:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=2,arrayin=tmp)
    return
  end subroutine free_scratch_c3

  subroutine save_scratch_c3 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(3)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(in),contiguous :: ptr(:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=3,arrayin=tmp)
    return
  end subroutine save_scratch_c3
  subroutine request_scratch_c4 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(4)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(out),contiguous :: ptr(:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    call segment_interface(len=len,tag=tag,task=1,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_c(1),il_c(2),il_c(3),il_c(4)/))
    return
  end subroutine request_scratch_c4

  subroutine request_saved_scratch_c4 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(4)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(out),contiguous :: ptr(:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    call segment_interface(len=len,tag=tag,task=4,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_c(1),il_c(2),il_c(3),il_c(4)/))
    return
  end subroutine request_saved_scratch_c4

  subroutine free_scratch_c4 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(4)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(in),contiguous :: ptr(:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=2,arrayin=tmp)
    return
  end subroutine free_scratch_c4

  subroutine save_scratch_c4 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(4)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(in),contiguous :: ptr(:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=3,arrayin=tmp)
    return
  end subroutine save_scratch_c4
  subroutine request_scratch_c5 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(5)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(out),contiguous :: ptr(:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    call segment_interface(len=len,tag=tag,task=1,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_c(1),il_c(2),il_c(3),il_c(4),il_c(5)/))
    return
  end subroutine request_scratch_c5

  subroutine request_saved_scratch_c5 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(5)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(out),contiguous :: ptr(:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    call segment_interface(len=len,tag=tag,task=4,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_c(1),il_c(2),il_c(3),il_c(4),il_c(5)/))
    return
  end subroutine request_saved_scratch_c5

  subroutine free_scratch_c5 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(5)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(in),contiguous :: ptr(:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=2,arrayin=tmp)
    return
  end subroutine free_scratch_c5

  subroutine save_scratch_c5 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(5)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(in),contiguous :: ptr(:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=3,arrayin=tmp)
    return
  end subroutine save_scratch_c5
  subroutine request_scratch_c6 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(6)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(out),contiguous :: ptr(:,:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    call segment_interface(len=len,tag=tag,task=1,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_c(1),il_c(2),il_c(3),il_c(4),il_c(5),il_c(6)/))
    return
  end subroutine request_scratch_c6

  subroutine request_saved_scratch_c6 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(6)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(out),contiguous :: ptr(:,:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    call segment_interface(len=len,tag=tag,task=4,arrayout=tmp)
    c_addr=c_loc(tmp)
    call c_f_pointer(c_addr,ptr,(/il_c(1),il_c(2),il_c(3),il_c(4),il_c(5),il_c(6)/))
    return
  end subroutine request_saved_scratch_c6

  subroutine free_scratch_c6 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(6)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(in),contiguous :: ptr(:,:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=2,arrayin=tmp)
    return
  end subroutine free_scratch_c6

  subroutine save_scratch_c6 (il_c,ptr,tagin)
    integer,intent(in) :: il_c(6)
    character(*),intent(in) ::tagin
    type(c_ptr) :: c_addr
    character(255) :: tag
    real(dp),pointer,contiguous :: tmp(:)
    COMPLEX(dp),pointer,intent(in),contiguous :: ptr(:,:,:,:,:,:)
    integer :: len   
    tag=trim(adjustl(tagin))
    len=product(il_c)*2
    c_addr=c_loc(ptr)
    call c_f_pointer(c_addr,tmp,(/len/))
    call segment_interface(len=len,tag=tag,task=3,arrayin=tmp)
    return
  end subroutine save_scratch_c6

end module scratch_interface
