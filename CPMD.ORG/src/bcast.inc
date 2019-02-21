subroutine mp_bcast_char_r0(data,  root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) ::  root,comm
  character(*), intent(inout)  :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_char_r0'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_CHARACTER,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+len(data)*size
  tim1=m_walltime()      
  call mpi_bcast(data,len(data),MPI_CHARACTER,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_char_r0
subroutine mp_bcast_char_r1(data, n, root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root,comm
  character(*), intent(inout) , dimension(*) :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_char_r1'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_CHARACTER,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+len(data)*n*size
  tim1=m_walltime()      
  call mpi_bcast(data,len(data)*n,MPI_CHARACTER,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_char_r1
subroutine mp_bcast_int4_r0(data,  root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) ::  root,comm
  integer(int_4), intent(inout)  :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_int4_r0'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_INTEGER,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+1*size
  tim1=m_walltime()      
  call mpi_bcast(data,1,MPI_INTEGER,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_int4_r0
subroutine mp_bcast_int4_r1(data, n, root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root,comm
  integer(int_4), intent(inout) , dimension(*) :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_int4_r1'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_INTEGER,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+n*size
  tim1=m_walltime()      
  call mpi_bcast(data,n,MPI_INTEGER,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_int4_r1
subroutine mp_bcast_int4_r2(data, n, root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root,comm
  integer(int_4), intent(inout) , dimension(1,*) :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_int4_r2'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_INTEGER,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+n*size
  tim1=m_walltime()      
  call mpi_bcast(data,n,MPI_INTEGER,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_int4_r2
subroutine mp_bcast_int8_r0(data,  root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) ::  root,comm
  integer(int_8), intent(inout)  :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_int8_r0'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_INTEGER8,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+1*size
  tim1=m_walltime()      
  call mpi_bcast(data,1,MPI_INTEGER8,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_int8_r0
subroutine mp_bcast_int8_r1(data, n, root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root,comm
  integer(int_8), intent(inout) , dimension(*) :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_int8_r1'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_INTEGER8,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+n*size
  tim1=m_walltime()      
  call mpi_bcast(data,n,MPI_INTEGER8,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_int8_r1
subroutine mp_bcast_real8_r0(data,  root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) ::  root,comm
  real(real_8), intent(inout)  :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_real8_r0'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_DOUBLE_PRECISION,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+1*size
  tim1=m_walltime()      
  call mpi_bcast(data,1,MPI_DOUBLE_PRECISION,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_real8_r0
subroutine mp_bcast_real8_r1(data, n, root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root,comm
  real(real_8), intent(inout) , dimension(*) :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_real8_r1'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_DOUBLE_PRECISION,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+n*size
  tim1=m_walltime()      
  call mpi_bcast(data,n,MPI_DOUBLE_PRECISION,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_real8_r1
subroutine mp_bcast_real8_r2(data, n, root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root,comm
  real(real_8), intent(inout) , dimension(1,*) :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_real8_r2'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_DOUBLE_PRECISION,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+n*size
  tim1=m_walltime()      
  call mpi_bcast(data,n,MPI_DOUBLE_PRECISION,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_real8_r2
subroutine mp_bcast_real8_r3(data, n, root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root,comm
  real(real_8), intent(inout) , dimension(1,1,*) :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_real8_r3'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_DOUBLE_PRECISION,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+n*size
  tim1=m_walltime()      
  call mpi_bcast(data,n,MPI_DOUBLE_PRECISION,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_real8_r3
subroutine mp_bcast_real8_r4(data, n, root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root,comm
  real(real_8), intent(inout) , dimension(1,1,1,*) :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_real8_r4'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_DOUBLE_PRECISION,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+n*size
  tim1=m_walltime()      
  call mpi_bcast(data,n,MPI_DOUBLE_PRECISION,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_real8_r4
subroutine mp_bcast_real8_r5(data, n, root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root,comm
  real(real_8), intent(inout) , dimension(1,1,1,1,*) :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_real8_r5'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_DOUBLE_PRECISION,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+n*size
  tim1=m_walltime()      
  call mpi_bcast(data,n,MPI_DOUBLE_PRECISION,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_real8_r5
subroutine mp_bcast_complex8_r0(data,  root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) ::  root,comm
  complex(real_8), intent(inout)  :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_complex8_r0'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_DOUBLE_COMPLEX,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+1*size
  tim1=m_walltime()      
  call mpi_bcast(data,1,MPI_DOUBLE_COMPLEX,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_complex8_r0
subroutine mp_bcast_complex8_r1(data, n, root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root,comm
  complex(real_8), intent(inout) , dimension(*) :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_complex8_r1'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_DOUBLE_COMPLEX,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+n*size
  tim1=m_walltime()      
  call mpi_bcast(data,n,MPI_DOUBLE_COMPLEX,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_complex8_r1
subroutine mp_bcast_complex8_r2(data, n, root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root,comm
  complex(real_8), intent(inout) , dimension(1,*) :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_complex8_r2'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_DOUBLE_COMPLEX,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+n*size
  tim1=m_walltime()      
  call mpi_bcast(data,n,MPI_DOUBLE_COMPLEX,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_complex8_r2
subroutine mp_bcast_complex8_r3(data, n, root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root,comm
  complex(real_8), intent(inout) , dimension(1,1,*) :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_complex8_r3'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_DOUBLE_COMPLEX,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+n*size
  tim1=m_walltime()      
  call mpi_bcast(data,n,MPI_DOUBLE_COMPLEX,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_complex8_r3
subroutine mp_bcast_complex8_r4(data, n, root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root,comm
  complex(real_8), intent(inout) , dimension(1,1,1,*) :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_complex8_r4'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_DOUBLE_COMPLEX,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+n*size
  tim1=m_walltime()      
  call mpi_bcast(data,n,MPI_DOUBLE_COMPLEX,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_complex8_r4
subroutine mp_bcast_logical_r0(data,  root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) ::  root,comm
  logical, intent(inout)  :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_logical_r0'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_LOGICAL,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+1*size
  tim1=m_walltime()      
  call mpi_bcast(data,1,MPI_LOGICAL,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_logical_r0
subroutine mp_bcast_logical_r1(data, n, root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root,comm
  logical, intent(inout) , dimension(*) :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_logical_r1'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_LOGICAL,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+n*size
  tim1=m_walltime()      
  call mpi_bcast(data,n,MPI_LOGICAL,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_logical_r1
subroutine mp_bcast_logical_r2(data, n, root,comm)
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, root,comm
  logical, intent(inout) , dimension(1,*) :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_logical_r2'
#ifdef __PARALLEL
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  if (first) then
    call mpi_type_size(MPI_LOGICAL,size,ierr)
    first=.false.
  end if
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+n*size
  tim1=m_walltime()      
  call mpi_bcast(data,n,MPI_LOGICAL,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_logical_r2
