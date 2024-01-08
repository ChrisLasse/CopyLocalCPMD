!=----------------------------------------------------------------------=
MODULE fftpw_batching
!=----------------------------------------------------------------------=

  USE cppt,                                     ONLY: nzh_r,&
                                                      indz_r
  USE fft,                                      ONLY: FFT_TYPE_DESCRIPTOR
  USE fftpw_types,                              ONLY: PW_fft_type_descriptor
  USE kinds,                                    ONLY: real_8
  USE mltfft_utils,                             ONLY: mltfft_fftw,&
                                                      mltfft_fftw_t
  USE mp_interface,                             ONLY: mp_startall,&
                                                      mp_waitall
  USE parac,                                    ONLY: parai
  USE system,                                   ONLY: parm,&
                                                      cntl
  USE timer,                                    ONLY: tihalt,&
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

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  REAL(DP), PARAMETER :: scal = 1.0

CONTAINS

SUBROUTINE Prepare_Psi( plac, psi, aux, remswitch, mythread )
  IMPLICIT NONE

  TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: plac
  INTEGER, INTENT(IN) :: remswitch, mythread
  COMPLEX(DP), INTENT(IN)  :: psi ( : , : )
  COMPLEX(DP), INTENT(OUT)  :: aux ( plac%nr3 , * ) !ns(parai%me+1)*z_group_size )

  INTEGER :: j, i, iter
  INTEGER :: offset, offset2
!  INTEGER :: isub, isub4
!  CHARACTER(*), PARAMETER :: procedureN = 'Prepare_Psi'

  INTEGER(INT64) :: time(2), cr

!  IF( cntl%fft_tune_batchsize ) THEN
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
!  ELSE
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
!  END IF

  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!----------Prepare_Psi Start---------------------------

  ! parai%me+1 in plac%thread_z_start eigentlich (parai%my_node-1)*parai%node_nproc+parai%node_me+1

  DO i = plac%thread_z_start( mythread+1, remswitch, parai%me+1, plac%which ), plac%thread_z_end( mythread+1, remswitch, parai%me+1, plac%which )
     iter = mod( i-1, plac%nsw(parai%me+1) ) + 1
     offset  = ( iter - 1 ) * plac%nr3
     offset2 = 2 * ( ( (i-1) / plac%nsw(parai%me+1) ) + 1 )

     DO j = 1, plac%prep_map(1,iter)-1
        aux( j, i ) = conjg( psi( indz_r( offset + j ), offset2 - 1 ) - (0.0d0,1.0d0) * psi( indz_r( offset + j ), offset2 ) )
     ENDDO
     DO j = 1, plac%prep_map(2,iter)-1
        aux( j, i ) = psi( nzh_r( offset + j ), offset2 - 1 ) + (0.0d0,1.0d0) * psi( nzh_r( offset + j ), offset2 )
     ENDDO

     DO j = plac%prep_map(3,iter), plac%prep_map(4,iter)
        aux( j, i ) = (0.d0, 0.d0)
     ENDDO

     DO j = plac%prep_map(5,iter)+1, plac%nr3
        aux( j, i ) = psi( nzh_r( offset + j ), offset2 - 1 ) + (0.0d0,1.0d0) * psi( nzh_r( offset + j ), offset2 )
     ENDDO
     DO j = plac%prep_map(6,iter)+1, plac%nr3
        aux( j, i ) = conjg( psi( indz_r( offset + j ), offset2 - 1 ) - (0.0d0,1.0d0) * psi( indz_r( offset + j ), offset2 ) )
     ENDDO

  ENDDO

!----------Prepare_Psi End-----------------------------
!------------------------------------------------------
  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) plac%time_adding( 1 ) = plac%time_adding( 1 ) + ( time(2) - time(1) )

!  CALL SYSTEM_CLOCK( count_rate = cr )
!  write(6,*) "ACTUAL", REAL( plac%time_adding( 1 ) / REAL( cr ) )

!  IF( cntl%fft_tune_batchsize ) THEN
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
!  ELSE
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
!  END IF

END SUBROUTINE Prepare_Psi

SUBROUTINE fft_com( plac, remswitch, work_buffer, which )
  IMPLICIT NONE

  INTEGER, INTENT(IN)                            :: remswitch, work_buffer, which
  TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT)    :: plac 

  CHARACTER(*), PARAMETER :: procedureN = 'fft_com'

  INTEGER :: ierr, isub, isub4

  IF( cntl%fft_tune_batchsize ) THEN
     CALL tiset(procedureN//'_tuning',isub4)
  ELSE
     CALL tiset(procedureN,isub)
  END IF

  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 1 ), ierr )
  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )
  
  CALL MP_STARTALL( plac%comm_sendrecv(1,plac%which), parai%send_handle(:,work_buffer,remswitch,which) )
  CALL MP_STARTALL( plac%comm_sendrecv(2,plac%which), parai%recv_handle(:,work_buffer,remswitch,which) )
  
  CALL MP_WAITALL( plac%comm_sendrecv(1,plac%which), parai%send_handle(:,work_buffer,remswitch,which) )
  CALL MP_WAITALL( plac%comm_sendrecv(2,plac%which), parai%recv_handle(:,work_buffer,remswitch,which) )

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )
  !CALL mpi_win_unlock_all( dfft%mpi_window( 1 ), ierr )

  IF( cntl%fft_tune_batchsize ) THEN
     CALL tihalt(procedureN//'_tuning',isub4)
  ELSE
     CALL tihalt(procedureN,isub)
  END IF

END SUBROUTINE fft_com

SUBROUTINE invfft_z_section( plac, aux, comm_mem_send, comm_mem_recv, batch_size, remswitch, mythread, nss )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: batch_size, remswitch, mythread
  TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: plac
  COMPLEX(DP), INTENT(INOUT) :: comm_mem_send( * ), comm_mem_recv( * )
  COMPLEX(DP), INTENT(INOUT)  :: aux ( plac%nr3 , * ) !ns(parai%me+1)*z_group_size )
  INTEGER, INTENT(IN) :: nss(*)

  INTEGER :: l, m, j, k, i
  INTEGER :: offset, kdest, ierr
!  INTEGER :: isub, isub4
!  CHARACTER(*), PARAMETER :: procedureN = 'invfft_z_section'

  INTEGER(INT64) :: time(3)

!  IF( cntl%fft_tune_batchsize ) THEN
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
!  ELSE
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
!  END IF

  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!------------z-FFT Start-------------------------------

  CALL mltfft_fftw_t('n','n',aux( : , plac%thread_z_start( mythread+1, remswitch, parai%me+1, plac%which ) : plac%thread_z_end( mythread+1, remswitch, parai%me+1, plac%which ) ), &
                   plac%nr3, plac%thread_z_sticks(mythread+1,remswitch,parai%me+1,plac%which), &
                   aux( : , plac%thread_z_start( mythread+1, remswitch, parai%me+1, plac%which ) : plac%thread_z_end( mythread+1, remswitch, parai%me+1, plac%which ) ), &
                   plac%nr3, plac%thread_z_sticks(mythread+1,remswitch,parai%me+1,plac%which), &
                   plac%nr3, plac%thread_z_sticks(mythread+1,remswitch,parai%me+1,plac%which),-1,scal,.FALSE.,mythread,parai%ncpus)

!-------------z-FFT End--------------------------------
!------------------------------------------------------
  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
!------------------------------------------------------
!---------Pre-Com-Copy Start---------------------------

  IF( parai%nnode .ne. 1 ) THEN

     !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 1 ), ierr )

     DO l = 1, parai%nnode
        IF( l .eq. parai%my_node+1 ) CYCLE 
        DO m = 1, parai%node_nproc
           j = (l-1)*parai%node_nproc + m
           offset = ( parai%node_me + (j-1)*parai%node_nproc ) * plac%small_chunks(plac%which) + (l-1)*(batch_size-1) * plac%big_chunks(plac%which)
           DO k = plac%thread_z_start( mythread+1, remswitch, parai%me+1, plac%which ), plac%thread_z_end( mythread+1, remswitch, parai%me+1, plac%which )
              kdest = offset + plac%nr3px * mod( (k-1), nss(parai%me+1) ) + ( (k-1) / nss(parai%me+1) ) * plac%big_chunks(plac%which)
              DO i = 1, plac%nr3p( j )
                 comm_mem_send( kdest + i ) = aux( i + plac%nr3p_offset( j ), k )
              ENDDO
           ENDDO
        ENDDO
     ENDDO

     !CALL mpi_win_unlock_all( dfft%mpi_window( 1 ), ierr )

  END IF

  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )

  DO m = 1, parai%node_nproc
     j = parai%my_node*parai%node_nproc + m
     offset = ( parai%node_me + (j-1)*parai%node_nproc ) * plac%small_chunks(plac%which) + parai%my_node*(batch_size-1) * plac%big_chunks(plac%which)
     DO k = plac%thread_z_start( mythread+1, remswitch, parai%me+1, plac%which ), plac%thread_z_end( mythread+1, remswitch, parai%me+1, plac%which )
        kdest = offset + plac%nr3px * mod( k-1, nss(parai%me+1) ) + ( (k-1) / nss(parai%me+1) ) * plac%big_chunks(plac%which)
        DO i = 1, plac%nr3p( j )
           comm_mem_recv( kdest + i ) = aux( i + plac%nr3p_offset( j ), k )
        ENDDO
     ENDDO
  ENDDO

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

!----------Pre-Com-Copy End----------------------------
!------------------------------------------------------
  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(3) )
  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) plac%time_adding( 2 ) = plac%time_adding( 2 ) + ( time(2) - time(1) )
  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) plac%time_adding( 3 ) = plac%time_adding( 3 ) + ( time(3) - time(2) )

!  IF( cntl%fft_tune_batchsize ) THEN
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
!  ELSE
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
!  END IF

END SUBROUTINE invfft_z_section 

SUBROUTINE invfft_y_section( plac, aux, comm_mem_recv, aux2_r, map_acinv, map_acinv_rem, counter, remswitch, mythread, my_nr1s )
  !$ USE omp_lib
  IMPLICIT NONE

  TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: plac
  INTEGER, INTENT(IN) :: remswitch, mythread, counter, my_nr1s
  COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( * )
  COMPLEX(DP), INTENT(INOUT) :: aux ( plac%my_nr3p * plac%nr2 * plac%nr1 , * ) !y_group_size
  COMPLEX(DP), INTENT(INOUT) :: aux2_r( : , : )
  INTEGER, INTENT(IN) :: map_acinv( * ), map_acinv_rem( * )

  INTEGER :: i, k, offset, iter
!  INTEGER :: isub, isub4
!  CHARACTER(*), PARAMETER :: procedureN = 'invfft_y_section'

  INTEGER(INT64) :: time(5)

!  IF( cntl%fft_tune_batchsize ) THEN
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
!  ELSE
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
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

!  IF( cntl%fft_tune_batchsize ) THEN
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
!  ELSE
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
!  END IF

  CONTAINS

    SUBROUTINE First_Part_y_section( aux2, map )
      
      Implicit NONE
      COMPLEX(DP), INTENT(INOUT) :: aux2( plac%nr2 , * ) !nr1s(parai%me+1) * plac%my_nr3p *  y_group_size )
      INTEGER, INTENT(IN)        :: map( * )


        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
      !------------------------------------------------------
      !--------After-Com-Copy Start--------------------------
      
        !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )
        
        DO i = plac%thread_y_start( mythread+1, remswitch, plac%which ), plac%thread_y_end( mythread+1, remswitch, plac%which )
           iter = mod( i-1, my_nr1s ) + 1
           offset = ( mod( i-1, my_nr1s * plac%my_nr3p ) + ( (i-1) / ( my_nr1s * plac%my_nr3p ) ) * plac%my_nr3p * my_nr1s ) * plac%nr2
           DO k = 1, plac%zero_acinv_start( iter, plac%which ) - 1
              aux2( k, i ) = comm_mem_recv( map( offset + k ) )
           END DO
           DO k = plac%zero_acinv_start( iter, plac%which ), plac%zero_acinv_end( iter, plac%which )
              aux2( k, i ) = (0.0_DP,0.0_DP)
           END DO
           DO k = plac%zero_acinv_end( iter, plac%which ) + 1, plac%nr2
              aux2( k, i ) = comm_mem_recv( map( offset + k ) )
           END DO
        END DO
      
        !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )
      
      !---------After-Com-Copy End---------------------------
      !------------------------------------------------------
        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
      !------------------------------------------------------
      !------------y-FFT Start-------------------------------

        CALL mltfft_fftw_t('n','n',aux2( : , plac%thread_y_start( mythread+1, remswitch, plac%which ) : plac%thread_y_end( mythread+1, remswitch, plac%which ) ), &
                         plac%nr2, plac%thread_y_sticks(mythread+1,remswitch,plac%which), &
                         aux2( : , plac%thread_y_start( mythread+1, remswitch, plac%which ) : plac%thread_y_end( mythread+1, remswitch, plac%which ) ), &
                         plac%nr2, plac%thread_y_sticks(mythread+1,remswitch,plac%which), &
                         plac%nr2, plac%thread_y_sticks(mythread+1,remswitch,plac%which),-1,scal,.FALSE.,mythread,parai%ncpus)
      
      !-------------y-FFT End--------------------------------
      !------------------------------------------------------
        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(3) )
        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) plac%time_adding( 4 ) = plac%time_adding( 4 ) + ( time(2) - time(1) )
        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) plac%time_adding( 5 ) = plac%time_adding( 5 ) + ( time(3) - time(2) )

    END SUBROUTINE First_Part_y_section

    SUBROUTINE Second_Part_y_section( aux2 )

      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT) :: aux2  ( my_nr1s * plac%my_nr3p * plac%nr2, * ) !y_group_size
      INTEGER :: ibatch

        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(4) )
      !------------------------------------------------------
      !-------------yx-scatter-------------------------------
      
        DO i = plac%thread_x_start( mythread+1, remswitch, plac%which ), plac%thread_x_end( mythread+1, remswitch, plac%which )
           offset = mod( i-1, plac%my_nr3p * plac%nr2 ) * plac%nr1
           ibatch = ( ( (i-1) / ( plac%my_nr3p * plac%nr2 ) ) + 1 )
           DO k = 1, plac%zero_scatter_start( plac%which ) - 1
              aux( offset + k, ibatch ) = aux2( plac%map_scatter_inv( offset + k, plac%which ), ibatch )
           END DO
           DO k = plac%zero_scatter_start( plac%which ), plac%zero_scatter_end( plac%which )
              aux( offset + k, ibatch ) = (0.0_DP, 0.0_DP)
           END DO
           DO k = plac%zero_scatter_end( plac%which ) + 1, plac%nr1
              aux( offset + k, ibatch ) = aux2( plac%map_scatter_inv( offset + k, plac%which ), ibatch )
           END DO
        END DO
      
      !-------------yx-scatter-------------------------------
      !------------------------------------------------------
        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(5) )
        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) plac%time_adding( 6 ) = plac%time_adding( 6 ) + ( time(5) - time(4) )

    END SUBROUTINE Second_Part_y_section

END SUBROUTINE invfft_y_section

SUBROUTINE invfft_x_section( plac, aux, remswitch, mythread )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: remswitch, mythread
  TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: plac
  COMPLEX(DP), INTENT(INOUT) :: aux( plac%nr1, * ) !plac%my_nr3p * plac%nr2 * x_group_size )

!  INTEGER :: isub, isub4
!  CHARACTER(*), PARAMETER :: procedureN = 'invfft_x_section'

  INTEGER(INT64) :: time(2)

!  IF( cntl%fft_tune_batchsize ) THEN
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
!  ELSE
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
!  END IF

  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!------------x-FFT Start-------------------------------

  CALL mltfft_fftw_t('n','n',aux( : , plac%thread_x_start( mythread+1, remswitch, plac%which ) : plac%thread_x_end( mythread+1, remswitch, plac%which ) ), &
                   plac%nr1, plac%thread_x_sticks(mythread+1,remswitch, plac%which), &
                   aux( : , plac%thread_x_start( mythread+1, remswitch, plac%which ) : plac%thread_x_end( mythread+1, remswitch, plac%which ) ), &
                   plac%nr1, plac%thread_x_sticks(mythread+1,remswitch, plac%which), &
                   plac%nr1, plac%thread_x_sticks(mythread+1,remswitch, plac%which),-1,scal,.FALSE.,mythread,parai%ncpus)

!-------------x-FFT End--------------------------------
!------------------------------------------------------
  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) plac%time_adding( 7 ) = plac%time_adding( 7 ) + ( time(2) - time(1) )

!  IF( cntl%fft_tune_batchsize ) THEN
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
!  ELSE
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
!  END IF

END SUBROUTINE invfft_x_section

SUBROUTINE fwfft_x_section( plac, aux_r, aux2, counter, remswitch, mythread, my_nr1s )
  IMPLICIT NONE

  TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: plac
  INTEGER, INTENT(IN) :: counter, remswitch, mythread, my_nr1s
  COMPLEX(DP), INTENT(INOUT)  :: aux_r( : ) !plac%my_nr3p * plac%nr2 * plac%nr1 , * ) !x_group_size )
  COMPLEX(DP), INTENT(INOUT) :: aux2( my_nr1s * plac%my_nr3p * plac%nr2 , * ) !x_group_size )

!  INTEGER :: isub, isub4
!  CHARACTER(*), PARAMETER :: procedureN = 'fwfft_x_section'

  INTEGER(INT64) :: time(4)

!  IF( cntl%fft_tune_batchsize ) THEN
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
!  ELSE
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
!  END IF

  Call First_Part_x_section( aux_r )

  !$  locks_omp( mythread+1, counter, 8 ) = .false.
  !$omp flush( locks_omp )
  !$  DO WHILE( ANY( locks_omp( :, counter, 8 ) ) )
  !$omp flush( locks_omp )
  !$  END DO
 
  Call Second_Part_x_section( aux_r )

!  IF( cntl%fft_tune_batchsize ) THEN
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
!  ELSE
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
!  END IF

  CONTAINS

    SUBROUTINE First_Part_x_section( aux )

      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT) :: aux( plac%nr1 , * ) !plac%my_nr3p * plac%nr2 * x_group_size )

      !------------------------------------------------------
      !------------x-FFT Start-------------------------------
      
        !$  locks_omp( mythread+1, counter, 9 ) = .false.
        !$omp flush( locks_omp )
        !$  DO WHILE( ANY( locks_omp( :, counter, 9 ) ) )
        !$omp flush( locks_omp )
        !$  END DO

        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )

        CALL mltfft_fftw_t('n','n',aux( : , plac%thread_x_start( mythread+1, remswitch, plac%which ) : plac%thread_x_end( mythread+1, remswitch, plac%which ) ), &
                         plac%nr1, plac%thread_x_sticks(mythread+1,remswitch, plac%which), &
                         aux( : , plac%thread_x_start( mythread+1, remswitch, plac%which ) : plac%thread_x_end( mythread+1, remswitch, plac%which ) ), &
                         plac%nr1, plac%thread_x_sticks(mythread+1,remswitch, plac%which), &
                         plac%nr1, plac%thread_x_sticks(mythread+1,remswitch, plac%which),1,scal,.FALSE.,mythread,parai%ncpus)
      
      !-------------x-FFT End--------------------------------
      !------------------------------------------------------
        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) plac%time_adding( 9 ) = plac%time_adding( 9 ) + ( time(2) - time(1) )

    END SUBROUTINE First_Part_x_section

    SUBROUTINE Second_Part_x_section( aux )

      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT)  :: aux( plac%my_nr3p * plac%nr2 * plac%nr1 , * ) !x_group_size )
      INTEGER :: i, j, offset, ibatch

        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(3) )
      !------------------------------------------------------
      !------Forward xy-scatter Start------------------------
      
        DO i = plac%thread_y_start( mythread+1, remswitch, plac%which ), plac%thread_y_end( mythread+1, remswitch, plac%which )
           ibatch = ( ( (i-1) / ( plac%my_nr3p * my_nr1s ) ) + 1 )
           offset = plac%nr2 * mod( i-1, my_nr1s * plac%my_nr3p )
           DO j = 1, plac%nr2
              aux2( j + offset , ibatch ) = aux( plac%map_scatter_fw( j + offset, plac%which ), ibatch )
           ENDDO
        ENDDO
      
      !-------Forward xy-scatter End-------------------------
      !------------------------------------------------------
        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(4) )
        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) plac%time_adding( 10 ) = plac%time_adding( 10 ) + ( time(4) - time(3) )

    END SUBROUTINE Second_Part_x_section

END SUBROUTINE fwfft_x_section

SUBROUTINE fwfft_y_section( plac, aux, comm_mem_send, comm_mem_recv, map_pcfw, batch_size, counter, remswitch, mythread, my_nr1s, nss )
  IMPLICIT NONE

  TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: plac
  INTEGER, INTENT(IN) :: counter, batch_size, remswitch, mythread, my_nr1s
  COMPLEX(DP), INTENT(INOUT)  :: comm_mem_send( * )
  COMPLEX(DP), INTENT(INOUT)  :: comm_mem_recv( * )
  COMPLEX(DP), INTENT(INOUT) :: aux( : , : )
  INTEGER, INTENT(IN) :: map_pcfw( * ), nss( * )

  INTEGER :: l, m, i, offset, j, k, ibatch, jter, offset2
!  INTEGER :: isub, isub4
  CHARACTER(*), PARAMETER :: procedureN = 'fwfft_y_section'

  INTEGER(INT64) :: time(4)

!  IF( cntl%fft_tune_batchsize ) THEN
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
!  ELSE
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
!  END IF

  Call First_Part_y_section( aux )

  !$  locks_omp( mythread+1, counter, 10 ) = .false.
  !$omp flush( locks_omp )
  !$  DO WHILE( ANY( locks_omp( :, counter, 10 ) ) )
  !$omp flush( locks_omp )
  !$  END DO
 
  Call Second_Part_y_section( aux )

!  IF( cntl%fft_tune_batchsize ) THEN
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
!  ELSE
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
!  END IF

  CONTAINS

    SUBROUTINE First_Part_y_section( aux2 )
      
      Implicit NONE
      COMPLEX(DP), INTENT(INOUT) :: aux2( plac%nr2 , * ) !nr1s(parai%me+1) * plac%my_nr3p *  y_group_size )

        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
      !------------------------------------------------------
      !------------y-FFT Start-------------------------------

        CALL mltfft_fftw_t('n','n',aux2( : , plac%thread_y_start( mythread+1, remswitch, plac%which ) : plac%thread_y_end( mythread+1, remswitch, plac%which ) ), &
                         plac%nr2, plac%thread_y_sticks(mythread+1,remswitch,plac%which), &
                         aux2( : , plac%thread_y_start( mythread+1, remswitch, plac%which ) : plac%thread_y_end( mythread+1, remswitch, plac%which ) ), &
                         plac%nr2, plac%thread_y_sticks(mythread+1,remswitch,plac%which), &
                         plac%nr2, plac%thread_y_sticks(mythread+1,remswitch,plac%which),1,scal,.FALSE.,mythread,parai%ncpus)
      
      !-------------y-FFT End--------------------------------
      !------------------------------------------------------
        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) plac%time_adding( 11 ) = plac%time_adding( 11 ) + ( time(2) - time(1) )

    END SUBROUTINE First_Part_y_section

    SUBROUTINE Second_Part_y_section( aux2 )

      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT) :: aux2  ( my_nr1s * plac%my_nr3p * plac%nr2, * ) !y_group_size

        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(3) )

      !------------------------------------------------------
      !---------Pre-Com-Copy Start---------------------------
      
        IF( parai%nnode .ne. 1 ) THEN
      
           !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 1 ), ierr )
      
           DO l = 1, parai%nnode
              IF( l .eq. parai%my_node+1 ) CYCLE
              DO m = 1, parai%node_nproc
                 i = (l-1)*parai%node_nproc + m
                 offset = ( parai%node_me + (i-1)*parai%node_nproc ) * plac%small_chunks(plac%which) + (l-1)*(batch_size-1) * plac%big_chunks(plac%which)
                 DO j = plac%thread_z_start( mythread+1, remswitch, i, plac%which ), plac%thread_z_end( mythread+1, remswitch, i, plac%which )
                    jter = mod( j-1, nss( i ) ) 
                    ibatch = ( (j-1) / nss( i ) )
                    offset2 =  ibatch * plac%big_chunks(plac%which)
                    DO k = 1, plac%my_nr3p
                       comm_mem_send( offset + offset2 + jter*plac%nr3px + k ) = &
                       aux2( map_pcfw( (i-1)*plac%small_chunks(plac%which) + jter*plac%nr3px + k ), ibatch+1 )
                    END DO
                 END DO
              END DO
           END DO
      
           !CALL mpi_win_unlock_all( dfft%mpi_window( 1 ), ierr )
      
        END IF
      
        !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )
      
        DO m = 1, parai%node_nproc
           i = parai%my_node*parai%node_nproc + m
           offset = ( parai%node_me + (i-1)*parai%node_nproc ) * plac%small_chunks(plac%which) + parai%my_node*(batch_size-1) * plac%big_chunks(plac%which)
           DO j = plac%thread_z_start( mythread+1, remswitch, i, plac%which ), plac%thread_z_end( mythread+1, remswitch, i, plac%which )
              jter = mod( j-1, nss( i ) ) 
              ibatch = ( (j-1) / nss( i ) )
              offset2 =  ibatch * plac%big_chunks(plac%which)
              DO k = 1, plac%my_nr3p
                 comm_mem_recv( offset + offset2 + jter*plac%nr3px + k ) = &
                 aux2( map_pcfw( (i-1)*plac%small_chunks(plac%which) + jter*plac%nr3px + k ), ibatch+1 )
              END DO
           END DO
        END DO
      
        !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )
      
      !----------Pre-Com-Copy End----------------------------
      !------------------------------------------------------
        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(4) )
        IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) plac%time_adding( 12 ) = plac%time_adding( 12 ) + ( time(4) - time(3) )

    END SUBROUTINE Second_Part_y_section
      
END SUBROUTINE fwfft_y_section

SUBROUTINE fwfft_z_section( plac, comm_mem_recv, aux, counter, batch_size, remswitch, mythread, nss, factor_in )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: counter, batch_size, remswitch, mythread
  TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: plac
  COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( * )
  COMPLEX(DP), INTENT(INOUT)  :: aux ( plac%nr3 , * ) !ns(parai%me+1)*z_group_size )
  INTEGER, INTENT(IN) :: nss( * )
  DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: factor_in

  DOUBLE PRECISION :: factor
  INTEGER :: j, l, k, i
  INTEGER :: offset, kfrom, ierr
!  INTEGER :: isub, isub4
!  CHARACTER(*), PARAMETER :: procedureN = 'fwfft_z_section'

  INTEGER(INT64) :: time(3)

!  IF( cntl%fft_tune_batchsize ) THEN
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
!  ELSE
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
!  END IF

  IF( present( factor_in ) ) THEN
     factor = factor_in
  ELSE
     factor = 1
  END IF

  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!--------After-Com-Copy Start--------------------------

  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )
 
  DO j = 1, parai%nnode
     DO l = 1, parai%node_nproc
        offset = ( parai%node_me*parai%node_nproc + (l-1) ) * plac%small_chunks(plac%which) + (j-1)*batch_size * plac%big_chunks(plac%which)
        DO k = plac%thread_z_start( mythread+1, remswitch, parai%me+1, plac%which ), plac%thread_z_end( mythread+1, remswitch, parai%me+1, plac%which )
           kfrom = offset + plac%nr3px * mod( k-1, nss(parai%me+1) ) + ( (k-1) / nss(parai%me+1) ) * plac%big_chunks(plac%which)
           DO i = 1, plac%nr3p( (j-1)*parai%node_nproc + l )
              aux( plac%nr3p_offset( (j-1)*parai%node_nproc + l ) + i, k ) = comm_mem_recv( kfrom + i ) * factor
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

!---------After-Com-Copy End---------------------------
!------------------------------------------------------
  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
!------------------------------------------------------
!------------z-FFT Start-------------------------------

  !$  locks_omp( mythread+1, counter, 11 ) = .false.
  !$omp flush( locks_omp )
  !$  DO WHILE( ANY( locks_omp( :, counter, 11 ) ) )
  !$omp flush( locks_omp )
  !$  END DO

  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(3) )


  CALL mltfft_fftw_t('n','n',aux( : , plac%thread_z_start( mythread+1, remswitch, parai%me+1, plac%which ) : plac%thread_z_end( mythread+1, remswitch, parai%me+1, plac%which ) ), &
                   plac%nr3, plac%thread_z_sticks(mythread+1,remswitch,parai%me+1,plac%which), &
                   aux( : , plac%thread_z_start( mythread+1, remswitch, parai%me+1, plac%which ) : plac%thread_z_end( mythread+1, remswitch, parai%me+1, plac%which ) ), &
                   plac%nr3, plac%thread_z_sticks(mythread+1,remswitch,parai%me+1,plac%which), &
                   plac%nr3, plac%thread_z_sticks(mythread+1,remswitch,parai%me+1,plac%which),1,scal,.FALSE.,mythread,parai%ncpus)

!-------------z-FFT End--------------------------------
!------------------------------------------------------
  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(4) )
  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) plac%time_adding( 13 ) = plac%time_adding( 13 ) + ( time(2) - time(1) )
  IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) plac%time_adding( 14 ) = plac%time_adding( 14 ) + ( time(4) - time(3) )

!  IF( cntl%fft_tune_batchsize ) THEN
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
!  ELSE
!     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
!  END IF

END SUBROUTINE fwfft_z_section

!=----------------------------------------------------------------------=
END MODULE fftpw_batching
!=----------------------------------------------------------------------=
