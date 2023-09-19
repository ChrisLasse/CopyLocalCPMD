!=----------------------------------------------------------------------=
MODULE fftpw_batching
!=----------------------------------------------------------------------=

  USE fftpw_param
  USE fftpw_types,                              ONLY: PW_fft_type_descriptor
  USE kinds,                                    ONLY: real_8
  USE mltfft_utils,                             ONLY: mltfft_fftw,&
                                                      mltfft_fftw_t
  USE mp_interface,                             ONLY: mp_startall,&
                                                      mp_waitall
  USE system,                                   ONLY: parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE iso_fortran_env
  IMPLICIT NONE
  PRIVATE

  LOGICAL, ALLOCATABLE, SAVE :: locks_omp(:,:,:)

  PUBLIC :: locks_omp

  PUBLIC :: Prepare_Psi
  PUBLIC :: fft_com
  PUBLIC :: invfft_z_section
  PUBLIC :: invfft_y_section
  PUBLIC :: invfft_x_section
  PUBLIC :: fwfft_x_section
  PUBLIC :: fwfft_y_section
  PUBLIC :: fwfft_z_section

CONTAINS

SUBROUTINE Prepare_Psi( dfft, psi, aux, remswitch, mythread )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: remswitch, mythread
  COMPLEX(DP), INTENT(IN)  :: psi ( : , : )
  COMPLEX(DP), INTENT(OUT)  :: aux ( dfft%nr3 , * ) !ns(dfft%mype+1)*z_group_size )

  INTEGER :: j, i, iter
  INTEGER :: offset, offset2
!  INTEGER :: isub, isub4
!  CHARACTER(*), PARAMETER :: procedureN = 'Prepare_Psi'

  INTEGER(INT64) :: time(2), cr

!  IF( dfft%fft_tuning ) THEN
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
!  ELSE
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
!  END IF

  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!----------Prepare_Psi Start---------------------------

  ! dfft%mype+1 in dfft%thread_z_start eigentlich (dfft%my_node-1)*dfft%node_task_size+dfft%my_node_rank+1

  DO i = dfft%thread_z_start( mythread+1, remswitch, dfft%mype+1 ), dfft%thread_z_end( mythread+1, remswitch, dfft%mype+1 )
     iter = mod( i-1, dfft%nsw(dfft%mype+1) ) + 1
     offset  = ( iter - 1 ) * dfft%nr3
     offset2 = 2 * ( ( (i-1) / dfft%nsw(dfft%mype+1) ) + 1 )

     DO j = 1, dfft%prep_map(1,iter)-1
        aux( j, i ) = conjg( psi( dfft%nlm_r( offset + j ), offset2 - 1 ) - (0.0d0,1.0d0) * psi( dfft%nlm_r( offset + j ), offset2 ) )
     ENDDO
     DO j = 1, dfft%prep_map(2,iter)-1
        aux( j, i ) = psi( dfft%nl_r( offset + j ), offset2 - 1 ) + (0.0d0,1.0d0) * psi( dfft%nl_r( offset + j ), offset2 )
     ENDDO

     DO j = dfft%prep_map(3,iter), dfft%prep_map(4,iter)
        aux( j, i ) = (0.d0, 0.d0)
     ENDDO

     DO j = dfft%prep_map(5,iter)+1, dfft%nr3
        aux( j, i ) = psi( dfft%nl_r( offset + j ), offset2 - 1 ) + (0.0d0,1.0d0) * psi( dfft%nl_r( offset + j ), offset2 )
     ENDDO
     DO j = dfft%prep_map(6,iter)+1, dfft%nr3
        aux( j, i ) = conjg( psi( dfft%nlm_r( offset + j ), offset2 - 1 ) - (0.0d0,1.0d0) * psi( dfft%nlm_r( offset + j ), offset2 ) )
     ENDDO

  ENDDO

!----------Prepare_Psi End-----------------------------
!------------------------------------------------------
  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) dfft%time_adding( 1 ) = dfft%time_adding( 1 ) + ( time(2) - time(1) )

!  CALL SYSTEM_CLOCK( count_rate = cr )
!  write(6,*) "ACTUAL", REAL( dfft%time_adding( 1 ) / REAL( cr ) )

!  IF( dfft%fft_tuning ) THEN
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
!  ELSE
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
!  END IF

END SUBROUTINE Prepare_Psi

SUBROUTINE fft_com( dfft, remswitch, work_buffer )
  IMPLICIT NONE

  INTEGER, INTENT(IN)                            :: remswitch, work_buffer
  TYPE(PW_fft_type_descriptor), INTENT(INOUT)    :: dfft

  CHARACTER(*), PARAMETER :: procedureN = 'fft_com'

  INTEGER :: ierr, isub, isub4

  IF( dfft%fft_tuning ) THEN
     CALL tiset(procedureN//'_tuning',isub4)
  ELSE
     CALL tiset(procedureN,isub)
  END IF

  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 1 ), ierr )
  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )

  IF( remswitch .eq. 1 ) THEN
  
     CALL MP_STARTALL( dfft%comm_sendrecv(1), dfft%send_handle(:,work_buffer) )
     CALL MP_STARTALL( dfft%comm_sendrecv(2), dfft%recv_handle(:,work_buffer) )
  
     CALL MP_WAITALL( dfft%comm_sendrecv(1), dfft%send_handle(:,work_buffer) )
     CALL MP_WAITALL( dfft%comm_sendrecv(2), dfft%recv_handle(:,work_buffer) )
  
  ELSE
  
     CALL MP_STARTALL( dfft%comm_sendrecv(1), dfft%send_handle_rem(:,work_buffer) )
     CALL MP_STARTALL( dfft%comm_sendrecv(2), dfft%recv_handle_rem(:,work_buffer) )
  
     CALL MP_WAITALL( dfft%comm_sendrecv(1), dfft%send_handle_rem(:,work_buffer) )
     CALL MP_WAITALL( dfft%comm_sendrecv(2), dfft%recv_handle_rem(:,work_buffer) )
  
  END IF

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )
  !CALL mpi_win_unlock_all( dfft%mpi_window( 1 ), ierr )

  IF( dfft%fft_tuning ) THEN
     CALL tihalt(procedureN//'_tuning',isub4)
  ELSE
     CALL tihalt(procedureN,isub)
  END IF

END SUBROUTINE fft_com

SUBROUTINE invfft_z_section( dfft, aux, comm_mem_send, comm_mem_recv, batch_size, remswitch, mythread )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: batch_size, remswitch, mythread
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(INOUT) :: comm_mem_send( * ), comm_mem_recv( * )
  COMPLEX(DP), INTENT(INOUT)  :: aux ( dfft%nr3 , * ) !ns(dfft%mype+1)*z_group_size )

  INTEGER :: l, m, j, k, i
  INTEGER :: offset, kdest, ierr
!  INTEGER :: isub, isub4
!  CHARACTER(*), PARAMETER :: procedureN = 'invfft_z_section'

  INTEGER(INT64) :: time(3)

!  IF( dfft%fft_tuning ) THEN
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
!  ELSE
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
!  END IF

  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!------------z-FFT Start-------------------------------

  CALL mltfft_fftw_t('n','n',aux( : , dfft%thread_z_start( mythread+1, remswitch, dfft%mype+1 ) : dfft%thread_z_end( mythread+1, remswitch, dfft%mype+1 ) ), &
                   dfft%nr3, dfft%thread_z_sticks(mythread+1,remswitch,dfft%mype+1), &
                   aux( : , dfft%thread_z_start( mythread+1, remswitch, dfft%mype+1 ) : dfft%thread_z_end( mythread+1, remswitch, dfft%mype+1 ) ), &
                   dfft%nr3, dfft%thread_z_sticks(mythread+1,remswitch,dfft%mype+1), &
                   dfft%nr3, dfft%thread_z_sticks(mythread+1,remswitch,dfft%mype+1),-1,scal,.FALSE.,mythread,dfft%nthreads)

!-------------z-FFT End--------------------------------
!------------------------------------------------------
  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
!------------------------------------------------------
!---------Pre-Com-Copy Start---------------------------

  IF( .not. dfft%single_node ) THEN

     !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 1 ), ierr )

     DO l = 1, dfft%nodes_numb
        IF( dfft%non_blocking .and. l .eq. dfft%my_node+1 ) CYCLE 
        DO m = 1, dfft%node_task_size
           j = (l-1)*dfft%node_task_size + m
           offset = ( dfft%my_node_rank + (j-1)*dfft%node_task_size ) * dfft%small_chunks + (l-1)*(batch_size-1) * dfft%big_chunks
           DO k = dfft%thread_z_start( mythread+1, remswitch, dfft%mype+1 ), dfft%thread_z_end( mythread+1, remswitch, dfft%mype+1 )
              kdest = offset + dfft%nr3px * mod( (k-1), dfft%nsw(dfft%mype+1) ) + ( (k-1) / dfft%nsw(dfft%mype+1) ) * dfft%big_chunks
              DO i = 1, dfft%nr3p( j )
                 comm_mem_send( kdest + i ) = aux( i + dfft%nr3p_offset( j ), k )
              ENDDO
           ENDDO
        ENDDO
     ENDDO

     !CALL mpi_win_unlock_all( dfft%mpi_window( 1 ), ierr )

  END IF

  IF( dfft%single_node .or. dfft%non_blocking ) THEN

     !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )

     DO m = 1, dfft%node_task_size
        j = dfft%my_node*dfft%node_task_size + m
        offset = ( dfft%my_node_rank + (j-1)*dfft%node_task_size ) * dfft%small_chunks + dfft%my_node*(batch_size-1) * dfft%big_chunks
        DO k = dfft%thread_z_start( mythread+1, remswitch, dfft%mype+1 ), dfft%thread_z_end( mythread+1, remswitch, dfft%mype+1 )
           kdest = offset + dfft%nr3px * mod( k-1, dfft%nsw(dfft%mype+1) ) + ( (k-1) / dfft%nsw(dfft%mype+1) ) * dfft%big_chunks
           DO i = 1, dfft%nr3p( j )
              comm_mem_recv( kdest + i ) = aux( i + dfft%nr3p_offset( j ), k )
           ENDDO
        ENDDO
     ENDDO

     !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

  END IF

!----------Pre-Com-Copy End----------------------------
!------------------------------------------------------
  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(3) )
  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) dfft%time_adding( 2 ) = dfft%time_adding( 2 ) + ( time(2) - time(1) )
  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) dfft%time_adding( 3 ) = dfft%time_adding( 3 ) + ( time(3) - time(2) )

!  IF( dfft%fft_tuning ) THEN
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
!  ELSE
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
!  END IF

END SUBROUTINE invfft_z_section 

SUBROUTINE invfft_y_section( dfft, aux, comm_mem_recv, aux2_r, map_acinv, map_acinv_rem, counter, remswitch, mythread )
  !$ USE omp_lib
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: remswitch, mythread, counter
  COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( * )
  COMPLEX(DP), INTENT(INOUT) :: aux ( dfft%my_nr3p * dfft%nr2 * dfft%nr1 , * ) !y_group_size
  COMPLEX(DP), INTENT(INOUT) :: aux2_r( : , : )
  INTEGER, INTENT(IN) :: map_acinv( * ), map_acinv_rem( * )

  INTEGER :: i, k, offset, iter
!  INTEGER :: isub, isub4
!  CHARACTER(*), PARAMETER :: procedureN = 'invfft_y_section'

  INTEGER(INT64) :: time(5)

!  IF( dfft%fft_tuning ) THEN
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
!  ELSE
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
!  END IF

  IF( remswitch .eq. 1 ) THEN 
     Call First_Part_y_section( aux2_r, map_acinv )
  ELSE
     Call First_Part_y_section( aux2_r, map_acinv_rem )
  END IF

  !$  locks_omp( mythread+1, counter, 7 ) = .false.
  !$omp flush( locks_omp )
  !$  DO WHILE( ANY( locks_omp( :, counter, 7 ) ) )
  !$omp flush( locks_omp )
  !$  END DO

  Call Second_Part_y_section( aux2_r )

!  IF( dfft%fft_tuning ) THEN
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
!  ELSE
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
!  END IF

  CONTAINS

    SUBROUTINE First_Part_y_section( aux2, map )
      
      Implicit NONE
      COMPLEX(DP), INTENT(INOUT) :: aux2( dfft%nr2 , * ) !nr1s(dfft%mype2+1) * dfft%my_nr3p *  y_group_size )
      INTEGER, INTENT(IN)        :: map( * )


        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
      !------------------------------------------------------
      !--------After-Com-Copy Start--------------------------
      
        !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )
        
        DO i = dfft%thread_y_start( mythread+1, remswitch ), dfft%thread_y_end( mythread+1, remswitch )
           iter = mod( i-1, dfft%my_nr1p ) + 1
           offset = ( mod( i-1, dfft%my_nr1p * dfft%my_nr3p ) + ( (i-1) / ( dfft%my_nr1p * dfft%my_nr3p ) ) * dfft%my_nr3p * dfft%my_nr1p ) * dfft%nr2
           DO k = 1, dfft%zero_acinv_start( iter ) - 1
              aux2( k, i ) = comm_mem_recv( map( offset + k ) )
           END DO
           DO k = dfft%zero_acinv_start( iter ), dfft%zero_acinv_end( iter )
              aux2( k, i ) = (0.0_DP,0.0_DP)
           END DO
           DO k = dfft%zero_acinv_end( iter ) + 1, dfft%nr2
              aux2( k, i ) = comm_mem_recv( map( offset + k ) )
           END DO
        END DO
      
        !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )
      
      !---------After-Com-Copy End---------------------------
      !------------------------------------------------------
        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
      !------------------------------------------------------
      !------------y-FFT Start-------------------------------

        CALL mltfft_fftw_t('n','n',aux2( : , dfft%thread_y_start( mythread+1, remswitch ) : dfft%thread_y_end( mythread+1, remswitch ) ), &
                         dfft%nr2, dfft%thread_y_sticks(mythread+1,remswitch), &
                         aux2( : , dfft%thread_y_start( mythread+1, remswitch ) : dfft%thread_y_end( mythread+1, remswitch ) ), &
                         dfft%nr2, dfft%thread_y_sticks(mythread+1,remswitch), &
                         dfft%nr2, dfft%thread_y_sticks(mythread+1,remswitch),-1,scal,.FALSE.,mythread,dfft%nthreads)
      
      !-------------y-FFT End--------------------------------
      !------------------------------------------------------
        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(3) )
        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) dfft%time_adding( 4 ) = dfft%time_adding( 4 ) + ( time(2) - time(1) )
        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) dfft%time_adding( 5 ) = dfft%time_adding( 5 ) + ( time(3) - time(2) )

    END SUBROUTINE First_Part_y_section

    SUBROUTINE Second_Part_y_section( aux2 )

      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT) :: aux2  ( dfft%my_nr1p * dfft%my_nr3p * dfft%nr2, * ) !y_group_size
      INTEGER :: ibatch

        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(4) )
      !------------------------------------------------------
      !-------------yx-scatter-------------------------------
      
        DO i = dfft%thread_x_start( mythread+1, remswitch ), dfft%thread_x_end( mythread+1, remswitch )
           offset = mod( i-1, dfft%my_nr3p * dfft%my_nr2p ) * dfft%nr1
           ibatch = ( ( (i-1) / ( dfft%my_nr3p * dfft%my_nr2p ) ) + 1 )
           DO k = 1, dfft%zero_scatter_start - 1
              aux( offset + k, ibatch ) = aux2( dfft%map_scatter_inv( offset + k ), ibatch )
           END DO
           DO k = dfft%zero_scatter_start, dfft%zero_scatter_end
              aux( offset + k, ibatch ) = (0.0_DP, 0.0_DP)
           END DO
           DO k = dfft%zero_scatter_end + 1, dfft%nr1
              aux( offset + k, ibatch ) = aux2( dfft%map_scatter_inv( offset + k ), ibatch )
           END DO
        END DO
      
      !-------------yx-scatter-------------------------------
      !------------------------------------------------------
        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(5) )
        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) dfft%time_adding( 6 ) = dfft%time_adding( 6 ) + ( time(5) - time(4) )

    END SUBROUTINE Second_Part_y_section

END SUBROUTINE invfft_y_section

SUBROUTINE invfft_x_section( dfft, aux, remswitch, mythread )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: remswitch, mythread
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(INOUT) :: aux( dfft%nr1, * ) !dfft%my_nr3p * dfft%nr2 * x_group_size )

!  INTEGER :: isub, isub4
!  CHARACTER(*), PARAMETER :: procedureN = 'invfft_x_section'

  INTEGER(INT64) :: time(2)

!  IF( dfft%fft_tuning ) THEN
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
!  ELSE
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
!  END IF

  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!------------x-FFT Start-------------------------------

  CALL mltfft_fftw_t('n','n',aux( : , dfft%thread_x_start( mythread+1, remswitch ) : dfft%thread_x_end( mythread+1, remswitch ) ), &
                   dfft%nr1, dfft%thread_x_sticks(mythread+1,remswitch), &
                   aux( : , dfft%thread_x_start( mythread+1, remswitch ) : dfft%thread_x_end( mythread+1, remswitch ) ), &
                   dfft%nr1, dfft%thread_x_sticks(mythread+1,remswitch), &
                   dfft%nr1, dfft%thread_x_sticks(mythread+1,remswitch),-1,scal,.FALSE.,mythread,dfft%nthreads)

!-------------x-FFT End--------------------------------
!------------------------------------------------------
  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) dfft%time_adding( 7 ) = dfft%time_adding( 7 ) + ( time(2) - time(1) )

!  IF( dfft%fft_tuning ) THEN
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
!  ELSE
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
!  END IF

END SUBROUTINE invfft_x_section

SUBROUTINE fwfft_x_section( dfft, aux_r, aux2, counter, remswitch, mythread )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: counter, remswitch, mythread
  COMPLEX(DP), INTENT(INOUT)  :: aux_r( : , : ) !dfft%my_nr3p * dfft%nr2 * dfft%nr1 , * ) !x_group_size )
  COMPLEX(DP), INTENT(INOUT) :: aux2( dfft%nr1s * dfft%my_nr3p * dfft%nr2 , * ) !x_group_size )

!  INTEGER :: isub, isub4
!  CHARACTER(*), PARAMETER :: procedureN = 'fwfft_x_section'

  INTEGER(INT64) :: time(4)

!  IF( dfft%fft_tuning ) THEN
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
!  ELSE
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
!  END IF

  Call First_Part_x_section( aux_r )

  !$  locks_omp( mythread+1, counter, 8 ) = .false.
  !$omp flush( locks_omp )
  !$  DO WHILE( ANY( locks_omp( :, counter, 8 ) ) )
  !$omp flush( locks_omp )
  !$  END DO
 
  Call Second_Part_x_section( aux_r )

!  IF( dfft%fft_tuning ) THEN
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
!  ELSE
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
!  END IF

  CONTAINS

    SUBROUTINE First_Part_x_section( aux )

      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT) :: aux( dfft%nr1 , * ) !dfft%my_nr3p * dfft%nr2 * x_group_size )

      !------------------------------------------------------
      !------------x-FFT Start-------------------------------
      
        !$  locks_omp( mythread+1, counter, 9 ) = .false.
        !$omp flush( locks_omp )
        !$  DO WHILE( ANY( locks_omp( :, counter, 9 ) ) )
        !$omp flush( locks_omp )
        !$  END DO

        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )

        CALL mltfft_fftw_t('n','n',aux( : , dfft%thread_x_start( mythread+1, remswitch ) : dfft%thread_x_end( mythread+1, remswitch ) ), &
                         dfft%nr1, dfft%thread_x_sticks(mythread+1,remswitch), &
                         aux( : , dfft%thread_x_start( mythread+1, remswitch ) : dfft%thread_x_end( mythread+1, remswitch ) ), &
                         dfft%nr1, dfft%thread_x_sticks(mythread+1,remswitch), &
                         dfft%nr1, dfft%thread_x_sticks(mythread+1,remswitch),1,scal,.FALSE.,mythread,dfft%nthreads)
      
      !-------------x-FFT End--------------------------------
      !------------------------------------------------------
        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) dfft%time_adding( 9 ) = dfft%time_adding( 9 ) + ( time(2) - time(1) )

    END SUBROUTINE First_Part_x_section

    SUBROUTINE Second_Part_x_section( aux )

      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT)  :: aux( dfft%my_nr3p * dfft%nr2 * dfft%nr1 , * ) !x_group_size )
      INTEGER :: i, j, offset, ibatch

        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(3) )
      !------------------------------------------------------
      !------Forward xy-scatter Start------------------------
      
        DO i = dfft%thread_y_start( mythread+1, remswitch ), dfft%thread_y_end( mythread+1, remswitch )
           ibatch = ( ( (i-1) / ( dfft%my_nr3p * dfft%nr1s ) ) + 1 )
           offset = dfft%my_nr2p * mod( i-1, dfft%nr1s * dfft%my_nr3p )
           DO j = 1, dfft%my_nr2p
              aux2( j + offset , ibatch ) = aux( dfft%map_scatter_fw( j + offset ), ibatch )
           ENDDO
        ENDDO
      
      !-------Forward xy-scatter End-------------------------
      !------------------------------------------------------
        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(4) )
        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) dfft%time_adding( 10 ) = dfft%time_adding( 10 ) + ( time(4) - time(3) )

    END SUBROUTINE Second_Part_x_section

END SUBROUTINE fwfft_x_section

SUBROUTINE fwfft_y_section( dfft, aux, comm_mem_send, comm_mem_recv, map_pcfw, batch_size, counter, remswitch, mythread )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: counter, batch_size, remswitch, mythread
  COMPLEX(DP), INTENT(INOUT)  :: comm_mem_send( * )
  COMPLEX(DP), INTENT(INOUT)  :: comm_mem_recv( * )
  COMPLEX(DP), INTENT(INOUT) :: aux( : , : )
  INTEGER, INTENT(IN) :: map_pcfw( * )

  INTEGER :: l, m, i, offset, j, k, ibatch, jter, offset2
!  INTEGER :: isub, isub4
  CHARACTER(*), PARAMETER :: procedureN = 'fwfft_y_section'

  INTEGER(INT64) :: time(4)

!  IF( dfft%fft_tuning ) THEN
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
!  ELSE
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
!  END IF

  Call First_Part_y_section( aux )

  !$  locks_omp( mythread+1, counter, 10 ) = .false.
  !$omp flush( locks_omp )
  !$  DO WHILE( ANY( locks_omp( :, counter, 10 ) ) )
  !$omp flush( locks_omp )
  !$  END DO
 
  Call Second_Part_y_section( aux )

!  IF( dfft%fft_tuning ) THEN
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
!  ELSE
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
!  END IF

  CONTAINS

    SUBROUTINE First_Part_y_section( aux2 )
      
      Implicit NONE
      COMPLEX(DP), INTENT(INOUT) :: aux2( dfft%nr2 , * ) !nr1s(dfft%mype2+1) * dfft%my_nr3p *  y_group_size )

        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
      !------------------------------------------------------
      !------------y-FFT Start-------------------------------

        CALL mltfft_fftw_t('n','n',aux2( : , dfft%thread_y_start( mythread+1, remswitch ) : dfft%thread_y_end( mythread+1, remswitch ) ), &
                         dfft%nr2, dfft%thread_y_sticks(mythread+1,remswitch), &
                         aux2( : , dfft%thread_y_start( mythread+1, remswitch ) : dfft%thread_y_end( mythread+1, remswitch ) ), &
                         dfft%nr2, dfft%thread_y_sticks(mythread+1,remswitch), &
                         dfft%nr2, dfft%thread_y_sticks(mythread+1,remswitch),1,scal,.FALSE.,mythread,dfft%nthreads)
      
      !-------------y-FFT End--------------------------------
      !------------------------------------------------------
        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) dfft%time_adding( 11 ) = dfft%time_adding( 11 ) + ( time(2) - time(1) )

    END SUBROUTINE First_Part_y_section

    SUBROUTINE Second_Part_y_section( aux2 )

      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT) :: aux2  ( dfft%my_nr1p * dfft%my_nr3p * dfft%nr2, * ) !y_group_size

        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(3) )

      !------------------------------------------------------
      !---------Pre-Com-Copy Start---------------------------
      
        IF( .not. dfft%single_node ) THEN
      
           !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 1 ), ierr )
      
           DO l = 1, dfft%nodes_numb
              IF( dfft%non_blocking .and. l .eq. dfft%my_node+1 ) CYCLE
              DO m = 1, dfft%node_task_size
                 i = (l-1)*dfft%node_task_size + m
                 offset = ( dfft%my_node_rank + (i-1)*dfft%node_task_size ) * dfft%small_chunks + (l-1)*(batch_size-1) * dfft%big_chunks
                 DO j = dfft%thread_z_start( mythread+1, remswitch, i ), dfft%thread_z_end( mythread+1, remswitch, i )
                    jter = mod( j-1, dfft%nsw( i ) ) 
                    ibatch = ( (j-1) / dfft%nsw( i ) )
                    offset2 =  ibatch * dfft%big_chunks
                    DO k = 1, dfft%my_nr3p
                       comm_mem_send( offset + offset2 + jter*dfft%nr3px + k ) = &
                       aux2( dfft%map_pcfw( (i-1)*dfft%small_chunks + jter*dfft%nr3px + k ), ibatch+1 )
                    END DO
                 END DO
              END DO
           END DO
      
           !CALL mpi_win_unlock_all( dfft%mpi_window( 1 ), ierr )
      
        END IF
      
        !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )
      
        IF( .not. dfft%single_node .and. dfft%non_blocking ) THEN
      
           DO m = 1, dfft%node_task_size
              i = dfft%my_node*dfft%node_task_size + m
              offset = ( dfft%my_node_rank + (i-1)*dfft%node_task_size ) * dfft%small_chunks + dfft%my_node*(batch_size-1) * dfft%big_chunks
              DO j = dfft%thread_z_start( mythread+1, remswitch, i ), dfft%thread_z_end( mythread+1, remswitch, i )
                 jter = mod( j-1, dfft%nsw( i ) ) 
                 ibatch = ( (j-1) / dfft%nsw( i ) )
                 offset2 =  ibatch * dfft%big_chunks
                 DO k = 1, dfft%my_nr3p
                    comm_mem_recv( offset + offset2 + jter*dfft%nr3px + k ) = &
                    aux2( dfft%map_pcfw( (i-1)*dfft%small_chunks + jter*dfft%nr3px + k ), ibatch+1 )
                 END DO
              END DO
           END DO
      
        END IF
      
        IF( dfft%single_node ) THEN
      
           DO m = 1, dfft%node_task_size
              i = dfft%my_node*dfft%node_task_size + m
              offset = ( dfft%my_node_rank + (i-1)*dfft%node_task_size ) * dfft%small_chunks + dfft%my_node*(batch_size-1) * dfft%big_chunks
              DO j = dfft%thread_z_start( mythread+1, remswitch, i ), dfft%thread_z_end( mythread+1, remswitch, i )
                 jter = mod( j-1, dfft%nsw( i ) ) 
                 ibatch = ( (j-1) / dfft%nsw( i ) )
                 offset2 =  ibatch * dfft%big_chunks
                 DO k = 1, dfft%my_nr3p
                    comm_mem_recv( offset + offset2 + jter*dfft%nr3px + k ) = &
                    aux2( dfft%map_pcfw( (i-1)*dfft%small_chunks + jter*dfft%nr3px + k ), ibatch+1 )
                 END DO
              END DO
           END DO
      
        END IF
      
        !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )
      
      !----------Pre-Com-Copy End----------------------------
      !------------------------------------------------------
        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(4) )
        IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) dfft%time_adding( 12 ) = dfft%time_adding( 12 ) + ( time(4) - time(3) )

    END SUBROUTINE Second_Part_y_section
      
END SUBROUTINE fwfft_y_section

SUBROUTINE fwfft_z_section( dfft, comm_mem_recv, aux, counter, batch_size, remswitch, mythread )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: counter, batch_size, remswitch, mythread
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( * )
  COMPLEX(DP), INTENT(INOUT)  :: aux ( dfft%nr3 , * ) !ns(dfft%mype+1)*z_group_size )

  INTEGER :: j, l, k, i
  INTEGER :: offset, kfrom, ierr
!  INTEGER :: isub, isub4
!  CHARACTER(*), PARAMETER :: procedureN = 'fwfft_z_section'

  INTEGER(INT64) :: time(3)

!  IF( dfft%fft_tuning ) THEN
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
!  ELSE
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
!  END IF

  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!--------After-Com-Copy Start--------------------------

  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )
 
  DO j = 1, dfft%nodes_numb
     DO l = 1, dfft%node_task_size
        offset = ( dfft%my_node_rank*dfft%node_task_size + (l-1) ) * dfft%small_chunks + (j-1)*batch_size * dfft%big_chunks
        DO k = dfft%thread_z_start( mythread+1, remswitch, dfft%mype+1 ), dfft%thread_z_end( mythread+1, remswitch, dfft%mype+1 )
           kfrom = offset + dfft%nr3px * mod( k-1, dfft%nsw(dfft%mype+1) ) + ( (k-1) / dfft%nsw(dfft%mype+1) ) * dfft%big_chunks
           DO i = 1, dfft%nr3p( (j-1)*dfft%node_task_size + l )
              aux( dfft%nr3p_offset( (j-1)*dfft%node_task_size + l ) + i, k ) = comm_mem_recv( kfrom + i )
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

!---------After-Com-Copy End---------------------------
!------------------------------------------------------
  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
!------------------------------------------------------
!------------z-FFT Start-------------------------------

  !$  locks_omp( mythread+1, counter, 11 ) = .false.
  !$omp flush( locks_omp )
  !$  DO WHILE( ANY( locks_omp( :, counter, 11 ) ) )
  !$omp flush( locks_omp )
  !$  END DO

  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(3) )


  CALL mltfft_fftw_t('n','n',aux( : , dfft%thread_z_start( mythread+1, remswitch, dfft%mype+1 ) : dfft%thread_z_end( mythread+1, remswitch, dfft%mype+1 ) ), &
                   dfft%nr3, dfft%thread_z_sticks(mythread+1,remswitch,dfft%mype+1), &
                   aux( : , dfft%thread_z_start( mythread+1, remswitch, dfft%mype+1 ) : dfft%thread_z_end( mythread+1, remswitch, dfft%mype+1 ) ), &
                   dfft%nr3, dfft%thread_z_sticks(mythread+1,remswitch,dfft%mype+1), &
                   dfft%nr3, dfft%thread_z_sticks(mythread+1,remswitch,dfft%mype+1),1,scal,.FALSE.,mythread,dfft%nthreads)

!-------------z-FFT End--------------------------------
!------------------------------------------------------
  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) CALL SYSTEM_CLOCK( time(4) )
  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) dfft%time_adding( 13 ) = dfft%time_adding( 13 ) + ( time(2) - time(1) )
  IF( mythread .eq. 1 .or. dfft%nthreads .eq. 1 ) dfft%time_adding( 14 ) = dfft%time_adding( 14 ) + ( time(4) - time(3) )

!  IF( dfft%fft_tuning ) THEN
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
!  ELSE
!     IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
!  END IF

END SUBROUTINE fwfft_z_section

!=----------------------------------------------------------------------=
END MODULE fftpw_batching
!=----------------------------------------------------------------------=
