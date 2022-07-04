!=----------------------------------------------------------------------=
MODULE fftpw_batching
!=----------------------------------------------------------------------=

  USE fftpw_legacy_routines,                    ONLY: fft_1D,&
                                                      fft_scatter_xy
  USE fftpw_param
  USE fftpw_types,                              ONLY: PW_fft_type_descriptor
  USE iso_fortran_env
  IMPLICIT NONE
  PRIVATE


  PUBLIC :: Prepare_Psi_overlapp
  PUBLIC :: invfft_pre_com
  PUBLIC :: fft_com
  PUBLIC :: invfft_after_com
  PUBLIC :: Apply_V
  PUBLIC :: fwfft_pre_com
  PUBLIC :: fwfft_after_com
  PUBLIC :: Accumulate_Psi_overlapp

CONTAINS

SUBROUTINE Prepare_Psi_overlapp( dfft, psi, ibatch, ngms, batch_size, howmany )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ibatch, ngms, batch_size, howmany
  COMPLEX(DP), INTENT(IN)  :: psi ( ngms, howmany )
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft

  INTEGER :: j, i
  INTEGER :: offset

  INTEGER(INT64) :: time(2)

  dfft%counter = dfft%counter + 1

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!----------Prepare_Psi Start---------------------------

  IF( howmany .eq. 2 ) THEN
     !$omp parallel do private( i, j, offset )
     DO i = 1, dfft%nsw( dfft%mype+1 )
        offset = (i-1)*dfft%nr3
        !$omp simd
        DO j = 1, dfft%zero_prep_start(i,1)-1
           dfft%aux( offset + j ) = psi( dfft%nl_r( offset + j ), 1 ) + (0.0d0,1.0d0) * psi( dfft%nl_r( offset + j ), 2 )
        ENDDO
        !$omp simd
        DO j = dfft%zero_prep_start(i,2), dfft%zero_prep_end(i,2)
           dfft%aux( offset + j ) = (0.d0, 0.d0)
        ENDDO
        !$omp simd
        DO j = dfft%zero_prep_end(i,1)+1, dfft%nr3
           dfft%aux( offset + j ) = psi( dfft%nl_r( offset + j ), 1 ) + (0.0d0,1.0d0) * psi( dfft%nl_r( offset + j ), 2 )
        ENDDO
     ENDDO
     !$omp end parallel do
  
     !$omp parallel do private( i, j, offset )
     DO i = 1, dfft%nsw( dfft%mype+1 )
        offset = (i-1)*dfft%nr3
        !$omp simd
        DO j = 1, dfft%zero_prep_start(i,3)-1
           dfft%aux( offset + j ) = conjg( psi( dfft%nlm_r( offset + j ), 1 ) - (0.0d0,1.0d0) * psi( dfft%nlm_r( offset + j ), 2 ) )
        ENDDO
        !$omp simd
        DO j = dfft%zero_prep_end(i,3)+1, dfft%nr3
           dfft%aux( offset + j ) = conjg( psi( dfft%nlm_r( offset + j ), 1 ) - (0.0d0,1.0d0) * psi( dfft%nlm_r( offset + j ), 2 ) )
        ENDDO
     ENDDO
     !$omp end parallel do
  ELSE
     !$omp parallel do private( i, j, offset )
     DO i = 1, dfft%nsw( dfft%mype+1 )
        offset = (i-1)*dfft%nr3
        !$omp simd
        DO j = 1, dfft%zero_prep_start(i,1)-1
           dfft%aux( offset + j ) = psi( dfft%nl_r( offset + j ), 1 )
        ENDDO
        !$omp simd
        DO j = dfft%zero_prep_start(i,2), dfft%zero_prep_end(i,2)
           dfft%aux( offset + j ) = (0.d0, 0.d0)
        ENDDO
        !$omp simd
        DO j = dfft%zero_prep_end(i,1)+1, dfft%nr3
           dfft%aux( offset + j ) = psi( dfft%nl_r( offset + j ), 1 )
        ENDDO
     ENDDO
     !$omp end parallel do
  
     !$omp parallel do private( i, j, offset )
     DO i = 1, dfft%nsw( dfft%mype+1 )
        offset = (i-1)*dfft%nr3
        !$omp simd
        DO j = 1, dfft%zero_prep_start(i,3)-1
           dfft%aux( offset + j ) = conjg( psi( dfft%nlm_r( offset + j ), 1 ) )
        ENDDO
        !$omp simd
        DO j = dfft%zero_prep_end(i,3)+1, dfft%nr3
           dfft%aux( offset + j ) = conjg( psi( dfft%nlm_r( offset + j ), 1 ) )
        ENDDO
     ENDDO
     !$omp end parallel do
  END IF

!----------Prepare_Psi End-----------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )

  dfft%time_adding( 1 ) = dfft%time_adding( 1 ) + ( time(2) - time(1) )

END SUBROUTINE Prepare_Psi_overlapp

SUBROUTINE invfft_pre_com( dfft, comm_mem_send, comm_mem_recv, ibatch, batch_size )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ibatch, batch_size
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(INOUT) :: comm_mem_send( : ), comm_mem_recv( : )

  INTEGER :: l, m, j, k, i
  INTEGER :: offset, kfrom, kdest, ierr

  INTEGER(INT64) :: time(3)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!------------z-FFT Start-------------------------------

  CALL fft_1D( dfft%aux, dfft%nsw(dfft%mype+1), dfft%nr3, 2 )

!-------------z-FFT End--------------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )
!------------------------------------------------------
!---------Pre-Com-Copy Start---------------------------

  IF( .not. dfft%single_node ) THEN

     !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 1 ), ierr )

     !$omp parallel private( l, m, j, offset, k, kdest, kfrom, i )
     DO l = 1, dfft%nodes_numb
        IF( dfft%non_blocking .and. l .eq. dfft%my_node+1 ) CYCLE 
        DO m = 1, dfft%node_task_size
           j = (l-1)*dfft%node_task_size + m
           offset = ( dfft%my_node_rank + (j-1)*dfft%node_task_size ) * dfft%small_chunks + ( (l-1)*(batch_size-1) + (ibatch-1) ) * dfft%big_chunks
           !$omp do
           DO k = 0, dfft%nsw(dfft%mype+1)-1
              kdest = offset + dfft%nr3px * k 
              kfrom = dfft%nr3p_offset( j ) + dfft%nr3 * k
              !omp simd
              DO i = 1, dfft%nr3p( j )
                 comm_mem_send( kdest + i ) = dfft%aux( kfrom + i )
              ENDDO
           ENDDO
           !$omp end do nowait
        ENDDO
     ENDDO
     !$omp end parallel

     !CALL mpi_win_unlock_all( dfft%mpi_window( 1 ), ierr )

  END IF

  IF( dfft%single_node .or. dfft%non_blocking ) THEN

     !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )

     !$omp parallel private( m, j, offset, k, kdest, kfrom, i )
     DO m = 1, dfft%node_task_size
        j = dfft%my_node*dfft%node_task_size + m
        offset = ( dfft%my_node_rank + (j-1)*dfft%node_task_size ) * dfft%small_chunks + ( dfft%my_node*(batch_size-1) + (ibatch-1) ) * dfft%big_chunks
        !$omp do
        DO k = 0, dfft%nsw(dfft%mype+1)-1
           kdest = offset + dfft%nr3px * k 
           kfrom = dfft%nr3p_offset( j ) + dfft%nr3 * k
           !omp simd
           DO i = 1, dfft%nr3p( j )
              comm_mem_recv( kdest + i ) = dfft%aux( kfrom + i )
           ENDDO
        ENDDO
        !$omp end do nowait
     ENDDO
     !$omp end parallel

     !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

  END IF

!----------Pre-Com-Copy End----------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(3) )

  dfft%time_adding( 2 ) = dfft%time_adding( 2 ) + ( time(2) - time(1) )
  dfft%time_adding( 3 ) = dfft%time_adding( 3 ) + ( time(3) - time(2) )

END SUBROUTINE invfft_pre_com 

SUBROUTINE fft_com( dfft, comm_mem_send, comm_mem_recv, sendsize, intra_me, inter_node_comm, nodes_numb, inter_me, non_blocking )
  IMPLICIT NONE

  INTEGER, INTENT(IN)                            :: sendsize, intra_me, nodes_numb, inter_me
  TYPE(MPI_COMM), INTENT(IN)                     :: inter_node_comm
  LOGICAL, INTENT(IN)                            :: non_blocking
  TYPE(PW_fft_type_descriptor), INTENT(IN)       :: dfft
  COMPLEX(DP), INTENT(IN)                        :: comm_mem_send( * )
  COMPLEX(DP), INTENT(INOUT)                     :: comm_mem_recv( * )

  INTEGER :: ierr, i, f
  TYPE( MPI_REQUEST ) :: handle( (nodes_numb-1)*2 )

  IF( intra_me .eq. 0 ) THEN

     !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 1 ), ierr )
     !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )

     IF( .not. non_blocking ) THEN
   
        CALL mpi_alltoall( comm_mem_send, sendsize, MPI_DOUBLE_COMPLEX, comm_mem_recv, &
                           sendsize, MPI_DOUBLE_COMPLEX, inter_node_comm, ierr)
 
     ELSE

        f = 0
        DO i = 1, nodes_numb !SENDING

           IF( i-1 .eq. inter_me ) CYCLE
           f = f + 1
           CALL MPI_ISEND( comm_mem_send( 1 + (i-1)*sendsize ), sendsize, MPI_DOUBLE_COMPLEX, i-1, &
                           inter_me, inter_node_comm, handle( f ), ierr )

        END DO

        DO i = 1, nodes_numb !RECEIVING

           IF( i-1 .eq. inter_me ) CYCLE
           f = f + 1
           CALL MPI_IRECV( comm_mem_recv( 1 + (i-1)*sendsize ), sendsize, MPI_DOUBLE_COMPLEX, i-1, &
                           MPI_ANY_TAG, inter_node_comm, handle( f ), ierr )

        END DO

        CALL MPI_WAITALL( (nodes_numb-1)*2, handle, MPI_STATUSES_IGNORE, ierr )

     END IF

     !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )
     !CALL mpi_win_unlock_all( dfft%mpi_window( 1 ), ierr )

  END IF

END SUBROUTINE fft_com


SUBROUTINE invfft_after_com( dfft, f, comm_mem_recv, ibatch )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ibatch
  COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( : )
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(OUT) :: f( : )

  INTEGER :: i, j, k
  INTEGER :: offset1, offset2, ierr

  INTEGER(INT64) :: time(6)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!--------After-Com-Copy Start--------------------------

  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )

  IF( dfft%rem ) THEN
     !$omp parallel do private( i, j, k, offset1, offset2 )
     DO i = 1, dfft%my_nr3p
        DO j = 1, dfft%my_nr1p
           offset1 = (i-1)*dfft%my_nr1p*dfft%nr2 + (j-1)*dfft%nr2
           offset2 = offset1 + (ibatch-1)*dfft%my_nr3p*dfft%my_nr1p*dfft%nr2
           !omp simd
           DO k = 1, dfft%zero_acinv_start( j )-1
              dfft%aux( offset1 + k ) = &
              comm_mem_recv( dfft%map_acinv_rem( offset2 + k ) )
           END DO
           !omp simd
           DO k = dfft%zero_acinv_start( j ), dfft%zero_acinv_end( j )
              dfft%aux( offset1 + k ) = (0.0_DP,0.0_DP)
           END DO
           !omp simd
           DO k = dfft%zero_acinv_end( j )+1, dfft%nr2
              dfft%aux( offset1 + k ) = &
              comm_mem_recv( dfft%map_acinv_rem( offset2 + k ) )
           END DO
        END DO
     END DO
     !$omp end parallel do
  ELSE
     !$omp parallel do private( i, j, k, offset1, offset2 )
     DO i = 1, dfft%my_nr3p
        DO j = 1, dfft%my_nr1p
           offset1 = (i-1)*dfft%my_nr1p*dfft%nr2 + (j-1)*dfft%nr2
           offset2 = offset1 + (ibatch-1)*dfft%my_nr3p*dfft%my_nr1p*dfft%nr2
           !omp simd
           DO k = 1, dfft%zero_acinv_start( j )-1
              dfft%aux( offset1 + k ) = &
              comm_mem_recv( dfft%map_acinv( offset2 + k ) )
           END DO
           !omp simd
           DO k = dfft%zero_acinv_start( j ), dfft%zero_acinv_end( j )
              dfft%aux( offset1 + k ) = (0.0_DP,0.0_DP)
           END DO
           !omp simd
           DO k = dfft%zero_acinv_end( j )+1, dfft%nr2
              dfft%aux( offset1 + k ) = &
              comm_mem_recv( dfft%map_acinv( offset2 + k ) )
           END DO
        END DO
     END DO
     !$omp end parallel do
  END IF

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

!---------After-Com-Copy End---------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )
!------------------------------------------------------
!------------y-FFT Start-------------------------------

  CALL fft_1D( dfft%aux, dfft%nr1w(dfft%mype+1) * dfft%my_nr3p, dfft%nr2, 2 )

!-------------y-FFT End--------------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(3) )
!------------------------------------------------------
!------Forward xy-scatter Start------------------------
    
  CALL fft_scatter_xy( dfft, dfft%aux, f, dfft%nnr, 2 )

!-------Forward xy-scatter End-------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(4) )
!------------------------------------------------------
!------------z-FFT Start-------------------------------

  CALL fft_1D( f, dfft%my_nr2p * dfft%my_nr3p, dfft%nr1, 2 )

!-------------z-FFT End--------------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(5) )

  dfft%time_adding( 4 ) = dfft%time_adding( 4 ) + ( time(2) - time(1) )
  dfft%time_adding( 5 ) = dfft%time_adding( 5 ) + ( time(3) - time(2) )
  dfft%time_adding( 6 ) = dfft%time_adding( 6 ) + ( time(4) - time(3) )
  dfft%time_adding( 7 ) = dfft%time_adding( 7 ) + ( time(5) - time(4) )

END SUBROUTINE invfft_after_com

SUBROUTINE Apply_V( dfft, f, v )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(INOUT) :: f( : )
  REAL(DP), INTENT(IN) :: v( : )

  INTEGER :: j

  INTEGER(INT64) :: time(2)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!-----------Apply V Start------------------------------

  !$omp parallel do private ( j )
  DO j = 1, dfft%nnr
     f( j )= f( j ) * v(j)
  END DO
  !$omp end parallel do

!------------Apply V End-------------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )

  dfft%time_adding( 8 ) = dfft%time_adding( 8 ) + ( time(2) - time(1) )

END SUBROUTINE Apply_V

SUBROUTINE fwfft_pre_com( dfft, f, comm_mem_send, comm_mem_recv, ibatch, batch_size )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ibatch, batch_size
  COMPLEX(DP), INTENT(INOUT)  :: f( : )
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(OUT) :: comm_mem_send( : )
  COMPLEX(DP), INTENT(INOUT) :: comm_mem_recv( : )

  INTEGER :: i, j, k, l, m
  INTEGER :: offset, ierr

  INTEGER(INT64) :: time(5)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!------------z-FFT Start-------------------------------

  CALL fft_1D( f, dfft%my_nr2p * dfft%my_nr3p, dfft%nr1, -2 )

!-------------z-FFT End--------------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )
!------------------------------------------------------
!------Forward xy-scatter Start------------------------
    
  CALL fft_scatter_xy( dfft, dfft%aux, f, dfft%nnr, -2 )

!-------Forward xy-scatter End-------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(3) )
!------------------------------------------------------
!------------y-FFT Start-------------------------------

  CALL fft_1D( dfft%aux, dfft%nr1w(dfft%mype+1) * dfft%my_nr3p, dfft%nr2, -2 )

!-------------y-FFT End--------------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(4) )
!------------------------------------------------------
!---------Pre-Com-Copy Start---------------------------

  IF( .not. dfft%single_node ) THEN

     !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 1 ), ierr )

     !$omp parallel private( i, j, k, offset, l, m )
     DO l = 1, dfft%nodes_numb
        IF( dfft%non_blocking .and. l .eq. dfft%my_node+1 ) CYCLE
        DO m = 1, dfft%node_task_size
           i = (l-1)*dfft%node_task_size + m
           offset = ( dfft%my_node_rank + (i-1)*dfft%node_task_size ) * dfft%small_chunks + ( (l-1)*(batch_size-1) + (ibatch-1) ) * dfft%big_chunks
           !$omp do
           DO j = 1, dfft%nsw( i )
              !$omp simd
              DO k = 1, dfft%my_nr3p
                 comm_mem_send( offset + (j-1)*dfft%nr3px + k ) = &
                 dfft%aux( dfft%map_pcfw( (i-1)*dfft%small_chunks + (j-1)*dfft%nr3px + k ) )
              END DO
           END DO
           !$omp end do nowait
        END DO
     END DO
     !$omp end parallel

     !CALL mpi_win_unlock_all( dfft%mpi_window( 1 ), ierr )

  END IF

  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )

  IF( .not. dfft%single_node .and. dfft%non_blocking ) THEN

     !$omp parallel private( i, j, k, offset, m )
     DO m = 1, dfft%node_task_size
        i = dfft%my_node*dfft%node_task_size + m
        offset = ( dfft%my_node_rank + (i-1)*dfft%node_task_size ) * dfft%small_chunks + ( dfft%my_node*(batch_size-1) + (ibatch-1) ) * dfft%big_chunks
        !$omp do
        DO j = 1, dfft%nsw( i )
           !$omp simd
           DO k = 1, dfft%my_nr3p
              comm_mem_recv( offset + (j-1)*dfft%nr3px + k ) = &
              dfft%aux( dfft%map_pcfw( (i-1)*dfft%small_chunks + (j-1)*dfft%nr3px + k ) )
           END DO
        END DO
        !$omp end do nowait
     END DO
     !$omp end parallel

  END IF

  IF( dfft%single_node ) THEN

     !$omp parallel private( i, j, k, offset, m )
     DO m = 1, dfft%node_task_size
        i = dfft%my_node*dfft%node_task_size + m
        offset = ( dfft%my_node_rank+ (i-1)*dfft%node_task_size ) * dfft%small_chunks + ( dfft%my_node*(batch_size-1) + (ibatch-1) ) * dfft%big_chunks
        !$omp do
        DO j = 1, dfft%nsw( i )
           !$omp simd
           DO k = 1, dfft%my_nr3p
              comm_mem_recv( offset + (j-1)*dfft%nr3px + k ) = &
              dfft%aux( dfft%map_pcfw( (i-1)*dfft%small_chunks + (j-1)*dfft%nr3px + k ) )
           END DO
        END DO
        !$omp end do nowait
     END DO
     !$omp end parallel

  END IF

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

!----------Pre-Com-Copy End----------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(5) )

  dfft%time_adding( 9 ) = dfft%time_adding( 9 ) + ( time(2) - time(1) )
  dfft%time_adding( 10 ) = dfft%time_adding( 10 ) + ( time(3) - time(2) )
  dfft%time_adding( 11 ) = dfft%time_adding( 11 ) + ( time(4) - time(3) )
  dfft%time_adding( 12 ) = dfft%time_adding( 12 ) + ( time(5) - time(4) )

END SUBROUTINE fwfft_pre_com

SUBROUTINE fwfft_after_com( dfft, comm_mem_recv, ibatch, batch_size )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ibatch, batch_size
  COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( : )
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft

  INTEGER :: j, l, k, i
  INTEGER :: offset, kfrom, kdest, ierr

  INTEGER(INT64) :: time(3)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!--------After-Com-Copy Start--------------------------

  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )
 
  !$omp parallel private( kdest, kfrom, j, i, k, offset, l )
  DO j = 1, dfft%nodes_numb
     DO l = 1, dfft%node_task_size
        offset = ( dfft%my_node_rank*dfft%node_task_size + (l-1) ) * dfft%small_chunks + ( (j-1)*batch_size + (ibatch-1) ) * dfft%big_chunks
        !$omp do
        DO k = 0, dfft%nsw(dfft%mype+1)-1
           kfrom = offset + dfft%nr3px * k
           kdest = dfft%nr3p_offset( (j-1)*dfft%node_task_size + l ) + dfft%nr3 * k
           !$omp simd
           DO i = 1, dfft%nr3p( (j-1)*dfft%node_task_size + l )
              dfft%aux( kdest + i ) = &
              comm_mem_recv( kfrom + i )
           ENDDO
        ENDDO
        !$omp end do nowait
     ENDDO
  ENDDO
  !$omp end parallel

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

!---------After-Com-Copy End---------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )
!------------------------------------------------------
!------------z-FFT Start-------------------------------

  CALL fft_1D( dfft%aux, dfft%nsw(dfft%mype+1), dfft%nr3, -2 )

!-------------z-FFT End--------------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(3) )

  dfft%time_adding( 13 ) = dfft%time_adding( 13 ) + ( time(2) - time(1) )
  dfft%time_adding( 14 ) = dfft%time_adding( 14 ) + ( time(3) - time(2) )

END SUBROUTINE fwfft_after_com

SUBROUTINE Accumulate_Psi_overlapp( dfft, hpsi, ibatch, ngms, batch_size, howmany )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ibatch, batch_size, ngms, howmany
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(INOUT) :: hpsi( ngms, howmany )

  COMPLEX(DP) :: fp, fm
  INTEGER :: j, l, k, i
  INTEGER :: offset, kfrom, kdest

  INTEGER(INT64) :: time(2)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!--------Accumulate_Psi Start--------------------------

  IF( howmany .eq. 2 ) THEN
  
     !$omp parallel do private( j, fp, fm )
     DO j = 1, ngms
        fp = ( dfft%aux( dfft%nl(j) ) + dfft%aux( dfft%nlm(j) ) ) * dfft%tscale_gamma
        fm = ( dfft%aux( dfft%nl(j) ) - dfft%aux( dfft%nlm(j) ) ) * dfft%tscale_gamma
        hpsi ( j, 1 ) = hpsi ( j, 1 ) + cmplx(  dble(fp),  aimag(fm), KIND=DP )
        hpsi ( j, 2 ) = hpsi ( j, 2 ) + cmplx( aimag(fp), - dble(fm), KIND=DP )        
     END DO
     !$omp end parallel do
  
  ELSE
  
     !$omp parallel do private( j )
     DO j = 1, ngms
        hpsi( j, 1 ) = hpsi ( j, 1 ) + dfft%aux( dfft%nl( j ) ) * dfft%tscale
     END DO
     !$omp end parallel do
    
  END IF

!---------Accumulate_Psi End---------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )

  dfft%time_adding( 15 ) = dfft%time_adding( 15 ) + ( time(2) - time(1) )

END SUBROUTINE Accumulate_Psi_overlapp

!=----------------------------------------------------------------------=
END MODULE fftpw_batching
!=----------------------------------------------------------------------=
