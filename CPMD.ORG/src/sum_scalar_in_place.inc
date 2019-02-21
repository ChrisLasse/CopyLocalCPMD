subroutine mp_sum_in_place_int4_r0(data,comm)
#ifdef __PARALLEL
  use mpi
#endif
  implicit none
  ! Arguments
  integer, intent(in) :: comm
  integer(int_4), intent(inout) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_int4_r0'
#ifdef __PARALLEL
  integer :: ierr
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  integer(int_4) :: stack_buff
  if (first) then
    call mpi_type_size(MPI_INTEGER,size,ierr)
    first=.false.
  end if
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+1*size
  tim1=m_walltime()      
  call mpi_allreduce(data,stack_buff,1,MPI_INTEGER,&
         MPI_SUM,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  data=stack_buff
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_in_place_int4_r0
subroutine mp_sum_in_place_int8_r0(data,comm)
#ifdef __PARALLEL
  use mpi
#endif
  implicit none
  ! Arguments
  integer, intent(in) :: comm
  integer(int_8), intent(inout) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_int8_r0'
#ifdef __PARALLEL
  integer :: ierr
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  integer(int_8) :: stack_buff
  if (first) then
    call mpi_type_size(MPI_INTEGER8,size,ierr)
    first=.false.
  end if
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+1*size
  tim1=m_walltime()      
  call mpi_allreduce(data,stack_buff,1,MPI_INTEGER8,&
         MPI_SUM,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  data=stack_buff
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_in_place_int8_r0
subroutine mp_sum_in_place_real8_r0(data,comm)
#ifdef __PARALLEL
  use mpi
#endif
  implicit none
  ! Arguments
  integer, intent(in) :: comm
  real(real_8), intent(inout) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_real8_r0'
#ifdef __PARALLEL
  integer :: ierr
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  real(real_8) :: stack_buff
  if (first) then
    call mpi_type_size(MPI_DOUBLE_PRECISION,size,ierr)
    first=.false.
  end if
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+1*size
  tim1=m_walltime()      
  call mpi_allreduce(data,stack_buff,1,MPI_DOUBLE_PRECISION,&
         MPI_SUM,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  data=stack_buff
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_in_place_real8_r0
subroutine mp_sum_in_place_complex8_r0(data,comm)
#ifdef __PARALLEL
  use mpi
#endif
  implicit none
  ! Arguments
  integer, intent(in) :: comm
  complex(real_8), intent(inout) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_sum_in_place_complex8_r0'
#ifdef __PARALLEL
  integer :: ierr
  integer,save :: size
  logical :: first=.true.
  real(real_8) :: tim1,tim2
  complex(real_8) :: stack_buff
  if (first) then
    call mpi_type_size(MPI_DOUBLE_COMPLEX,size,ierr)
    first=.false.
  end if
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+1*size
  tim1=m_walltime()      
  call mpi_allreduce(data,stack_buff,1,MPI_DOUBLE_COMPLEX,&
         MPI_SUM,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  data=stack_buff
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_in_place_complex8_r0
