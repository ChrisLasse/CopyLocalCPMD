subroutine mp_all2all_int4_r1(data_in, data_out, n, comm )
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, comm
  integer(int_4), intent(in) , dimension(*) ::  data_in
  integer(int_4), intent(out) , dimension(*) :: data_out
  ! Variables
  integer :: ierr,nnodes
  character(*),parameter :: procedureN='mp_all2all_int4_r1'
  real(real_8) :: tim1,tim2
  logical :: first=.true.
  integer,save :: size
#ifdef __PARALLEL
  if (first) then
    call mpi_type_size(MPI_INTEGER,size,ierr)
    first=.false.
  end if
  call mpi_comm_size ( comm, nnodes, ierr )
  cmcal(ipar_aall)=cmcal(ipar_aall)+1.0d0
  cmlen(ipar_aall)=cmlen(ipar_aall)+n*nnodes*size
  tim1=m_walltime()     
  CALL mpi_alltoall(data_in,n,MPI_INTEGER,data_out,n,MPI_INTEGER,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_aall)=cmtim(ipar_aall)+tim2-tim1
#else
  data_out(1:n)=data_in(1:n)
#endif
end subroutine mp_all2all_int4_r1
subroutine mp_all2all_int8_r1(data_in, data_out, n, comm )
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, comm
  integer(int_8), intent(in) , dimension(*) ::  data_in
  integer(int_8), intent(out) , dimension(*) :: data_out
  ! Variables
  integer :: ierr,nnodes
  character(*),parameter :: procedureN='mp_all2all_int8_r1'
  real(real_8) :: tim1,tim2
  logical :: first=.true.
  integer,save :: size
#ifdef __PARALLEL
  if (first) then
    call mpi_type_size(MPI_INTEGER8,size,ierr)
    first=.false.
  end if
  call mpi_comm_size ( comm, nnodes, ierr )
  cmcal(ipar_aall)=cmcal(ipar_aall)+1.0d0
  cmlen(ipar_aall)=cmlen(ipar_aall)+n*nnodes*size
  tim1=m_walltime()     
  CALL mpi_alltoall(data_in,n,MPI_INTEGER8,data_out,n,MPI_INTEGER8,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_aall)=cmtim(ipar_aall)+tim2-tim1
#else
  data_out(1:n)=data_in(1:n)
#endif
end subroutine mp_all2all_int8_r1
subroutine mp_all2all_real8_r1(data_in, data_out, n, comm )
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, comm
  real(real_8), intent(in) , dimension(*) ::  data_in
  real(real_8), intent(out) , dimension(*) :: data_out
  ! Variables
  integer :: ierr,nnodes
  character(*),parameter :: procedureN='mp_all2all_real8_r1'
  real(real_8) :: tim1,tim2
  logical :: first=.true.
  integer,save :: size
#ifdef __PARALLEL
  if (first) then
    call mpi_type_size(MPI_DOUBLE_PRECISION,size,ierr)
    first=.false.
  end if
  call mpi_comm_size ( comm, nnodes, ierr )
  cmcal(ipar_aall)=cmcal(ipar_aall)+1.0d0
  cmlen(ipar_aall)=cmlen(ipar_aall)+n*nnodes*size
  tim1=m_walltime()     
  CALL mpi_alltoall(data_in,n,MPI_DOUBLE_PRECISION,data_out,n,MPI_DOUBLE_PRECISION,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_aall)=cmtim(ipar_aall)+tim2-tim1
#else
  data_out(1:n)=data_in(1:n)
#endif
end subroutine mp_all2all_real8_r1
subroutine mp_all2all_complex4_r1(data_in, data_out, n, comm )
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, comm
  complex(real_4), intent(in) , dimension(*) ::  data_in
  complex(real_4), intent(out) , dimension(*) :: data_out
  ! Variables
  integer :: ierr,nnodes
  character(*),parameter :: procedureN='mp_all2all_complex4_r1'
  real(real_8) :: tim1,tim2
  logical :: first=.true.
  integer,save :: size
#ifdef __PARALLEL
  if (first) then
    call mpi_type_size(MPI_COMPLEX,size,ierr)
    first=.false.
  end if
  call mpi_comm_size ( comm, nnodes, ierr )
  cmcal(ipar_aall)=cmcal(ipar_aall)+1.0d0
  cmlen(ipar_aall)=cmlen(ipar_aall)+n*nnodes*size
  tim1=m_walltime()     
  CALL mpi_alltoall(data_in,n,MPI_COMPLEX,data_out,n,MPI_COMPLEX,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_aall)=cmtim(ipar_aall)+tim2-tim1
#else
  data_out(1:n)=data_in(1:n)
#endif
end subroutine mp_all2all_complex4_r1
subroutine mp_all2all_complex8_r1(data_in, data_out, n, comm )
#ifdef __PARALLEL
  use mpi
#endif
  ! Arguments
  integer, intent(in) :: n, comm
  complex(real_8), intent(in) , dimension(*) ::  data_in
  complex(real_8), intent(out) , dimension(*) :: data_out
  ! Variables
  integer :: ierr,nnodes
  character(*),parameter :: procedureN='mp_all2all_complex8_r1'
  real(real_8) :: tim1,tim2
  logical :: first=.true.
  integer,save :: size
#ifdef __PARALLEL
  if (first) then
    call mpi_type_size(MPI_DOUBLE_COMPLEX,size,ierr)
    first=.false.
  end if
  call mpi_comm_size ( comm, nnodes, ierr )
  cmcal(ipar_aall)=cmcal(ipar_aall)+1.0d0
  cmlen(ipar_aall)=cmlen(ipar_aall)+n*nnodes*size
  tim1=m_walltime()     
  CALL mpi_alltoall(data_in,n,MPI_DOUBLE_COMPLEX,data_out,n,MPI_DOUBLE_COMPLEX,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_aall)=cmtim(ipar_aall)+tim2-tim1
#else
  data_out(1:n)=data_in(1:n)
#endif
end subroutine mp_all2all_complex8_r1
