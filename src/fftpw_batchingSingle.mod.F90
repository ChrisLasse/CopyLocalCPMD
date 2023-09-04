!=----------------------------------------------------------------------=
MODULE fftpw_batchingSingle
!=----------------------------------------------------------------------=

  USE cppt,                                     ONLY: hg
  USE density_utils,                            ONLY: build_density_sum
  USE elct,                                     ONLY: crge
  USE fftpw_legacy_routines,                    ONLY: fft_1D
  USE fftpw_param
  USE fftpw_types,                              ONLY: PW_fft_type_descriptor
  USE kinds,                                    ONLY: real_8
  USE system,                                   ONLY: parm
  USE timer,                                    ONLY: tihalt,&
                                                      tiset
  USE iso_fortran_env
  IMPLICIT NONE
  PRIVATE


  PUBLIC :: invfft_pre_com
  PUBLIC :: fft_com_single
  PUBLIC :: invfft_after_com
  PUBLIC :: fwfft_pre_com
  PUBLIC :: fwfft_after_com

CONTAINS

SUBROUTINE invfft_pre_com( dfft, f, comm_mem_send, comm_mem_recv, ns )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: ns( * )
  COMPLEX(DP), INTENT(INOUT) :: comm_mem_send( * ), comm_mem_recv( * )
  COMPLEX(DP), INTENT(INOUT)  :: f( dfft%nr3 , * ) !ns(dfft%mype+1)*z_group_size )

  INTEGER :: l, m, j, k, i
  INTEGER :: offset, kdest, ierr

!------------------------------------------------------
!------------z-FFT Start-------------------------------

  CALL fft_1D( f, ns(dfft%mype+1), dfft%nr3, 2 )

!-------------z-FFT End--------------------------------
!------------------------------------------------------
!------------------------------------------------------
!---------Pre-Com-Copy Start---------------------------

  IF( .not. dfft%single_node ) THEN

     !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 1 ), ierr )

     !$omp parallel private( l, m, j, offset, k, kdest, i )
     DO l = 1, dfft%nodes_numb
        IF( dfft%non_blocking .and. l .eq. dfft%my_node+1 ) CYCLE 
        DO m = 1, dfft%node_task_size
           j = (l-1)*dfft%node_task_size + m
           offset = ( dfft%my_node_rank + (j-1)*dfft%node_task_size ) * dfft%small_chunks
           !$omp do
           DO k = 1, ns(dfft%mype+1)
              kdest = offset + dfft%nr3px * (k-1)
              !omp simd
              DO i = 1, dfft%nr3p( j )
                 comm_mem_send( kdest + i ) = f( i + dfft%nr3p_offset( j ), k )
              ENDDO
           ENDDO
           !$omp end do nowait
        ENDDO
     ENDDO
     !$omp end parallel

     !CALL mpi_win_unlock_all( dfft%mpi_window( 1 ), ierr )

  END IF

  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )

  !$omp parallel private( m, j, offset, k, kdest, i )
  DO m = 1, dfft%node_task_size
     j = dfft%my_node*dfft%node_task_size + m
     offset = ( dfft%my_node_rank + (j-1)*dfft%node_task_size ) * dfft%small_chunks
     !$omp do
     DO k = 1, ns(dfft%mype+1)
        kdest = offset + dfft%nr3px * (k-1)
        !omp simd
        DO i = 1, dfft%nr3p( j )
           comm_mem_recv( kdest + i ) = f( i + dfft%nr3p_offset( j ), k )
        ENDDO
     ENDDO
     !$omp end do nowait
  ENDDO
  !$omp end parallel

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

!----------Pre-Com-Copy End----------------------------
!------------------------------------------------------

END SUBROUTINE invfft_pre_com 

SUBROUTINE fft_com_single( dfft, sendsize, nodes_numb )
  IMPLICIT NONE

  INTEGER, INTENT(IN)                            :: sendsize, nodes_numb
  TYPE(PW_fft_type_descriptor), INTENT(INOUT)    :: dfft

  CHARACTER(*), PARAMETER :: procedureN = 'fft_com_single'

  INTEGER :: ierr, isub, isub4

  IF( dfft%fft_tuning ) THEN
     CALL tiset(procedureN//'_tuning',isub4)
  ELSE
     CALL tiset(procedureN,isub)
  END IF

  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 1 ), ierr )
  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )

  CALL MP_STARTALL( dfft%comm_sendrecv(1), dfft%send_handle(:,1) )
  CALL MP_STARTALL( dfft%comm_sendrecv(2), dfft%recv_handle(:,1) )
  
  CALL MP_WAITALL( dfft%comm_sendrecv(1), dfft%send_handle(:,1) )
  CALL MP_WAITALL( dfft%comm_sendrecv(2), dfft%recv_handle(:,1) )

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )
  !CALL mpi_win_unlock_all( dfft%mpi_window( 1 ), ierr )

  IF( dfft%fft_tuning ) THEN
     CALL tihalt(procedureN//'_tuning',isub4)
  ELSE
     CALL tihalt(procedureN,isub)
  END IF

END SUBROUTINE fft_com_single

SUBROUTINE invfft_after_com( dfft, f, comm_mem_recv, aux, map_acinv, nr1s )
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( * )
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(OUT) :: f( * ) !dfft%my_nr3p * dfft%nr2 * dfft%nr1, * ) !batch_size )
  INTEGER, INTENT(IN) :: nr1s( * )
  COMPLEX(DP), INTENT(INOUT) :: aux( * ) !nr1s(dfft%mype2+1) * dfft%my_nr3p * dfft%nr2 * batch_size )
  INTEGER, INTENT(IN) :: map_acinv( * )

  INTEGER(INT64) :: time(6)
  INTEGER(INT64) :: auto_time(6)

  CALL invfft_y_portion( dfft, comm_mem_recv, aux( 1:nr1s(dfft%mype2+1)*dfft%my_nr3p*dfft%nr2 ), map_acinv, nr1s )

!------------------------------------------------------
!------Forward xy-scatter Start------------------------

  CALL fft_scatter_xy( dfft, aux( 1:nr1s(dfft%mype2+1)*dfft%my_nr3p*dfft%nr2 ), f, 2 )

!-------Forward xy-scatter End-------------------------
!------------------------------------------------------
!------------------------------------------------------
!------------x-FFT Start-------------------------------

  CALL fft_1D( f, dfft%my_nr2p * dfft%my_nr3p, dfft%nr1, 2 )

!-------------x-FFT End--------------------------------
!------------------------------------------------------

END SUBROUTINE invfft_after_com

SUBROUTINE invfft_y_portion( dfft, comm_mem_recv, aux, map_acinv, nr1s )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( * )
  INTEGER, INTENT(IN) :: nr1s( * )
  COMPLEX(DP), INTENT(INOUT) :: aux( dfft%nr2 , * ) !nr1s(dfft%mype2+1) * dfft%my_nr3p *  y_group_size )
  INTEGER, INTENT(IN) :: map_acinv( * )

  INTEGER :: i, jter, k, offset

!------------------------------------------------------
!--------After-Com-Copy Start--------------------------

  !CALL mpi_win_unlock_all( dfft%mpi_window( 1 ), ierr )

  !$omp parallel do private( i, jter, k, offset )
  DO i = 1, dfft%my_nr1p * dfft%my_nr3p
     jter = mod( i-1, dfft%my_nr1p ) + 1
     offset = (i-1) * dfft%nr2
     !omp simd
     DO k = 1, dfft%zero_acinv_start( jter ) - 1
        aux( k, i ) = comm_mem_recv( map_acinv( offset + k ) )
     END DO
     !omp simd
     DO k = dfft%zero_acinv_start( jter ), dfft%zero_acinv_end( jter )
        aux( k, i ) = (0.0_DP,0.0_DP)
     END DO
     !omp simd
     DO k = dfft%zero_acinv_end( jter ) + 1, dfft%nr2
        aux( k, i ) = comm_mem_recv( map_acinv( offset + k ) )
     END DO
  END DO
  !$omp end parallel do

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

!---------After-Com-Copy End---------------------------
!------------------------------------------------------
!------------------------------------------------------
!------------y-FFT Start-------------------------------

  CALL fft_1D( aux, nr1s(dfft%mype2+1) * dfft%my_nr3p, dfft%nr2, 2 )

!-------------y-FFT End--------------------------------
!------------------------------------------------------

END SUBROUTINE invfft_y_portion

SUBROUTINE fft_scatter_xy ( dfft, f_in, f_aux, isgn )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: isgn
  COMPLEX(DP), INTENT(INOUT) :: f_in  ( * ) !dfft%my_nr1p * dfft%my_nr3p * dfft%nr2, * ) !scatter_group_size )
  COMPLEX(DP), INTENT(INOUT) :: f_aux ( * ) !dfft%my_nr3p * dfft%nr2 * dfft%nr1 , * ) !dfft%my_nr3p * dfft%nr2 * dfft%nr1 )

  INTEGER :: i, k, j, offset

  IF( isgn .gt. 0 ) THEN

     !$omp parallel do private( i, k )
     DO i = 1, dfft%my_nr3p * dfft%my_nr2p
        DO k = 1, dfft%zero_scatter_start - 1
           f_aux( (i-1)*dfft%nr1 + k ) = f_in( dfft%map_scatter_inv( (i-1)*dfft%nr1 + k ) )
        END DO
        DO k = dfft%zero_scatter_start, dfft%zero_scatter_end
           f_aux( (i-1)*dfft%nr1 + k ) = (0.0_DP, 0.0_DP)
        END DO
        DO k = dfft%zero_scatter_end + 1, dfft%nr1
           f_aux( (i-1)*dfft%nr1 + k ) = f_in( dfft%map_scatter_inv( (i-1)*dfft%nr1 + k ) )
        END DO
     END DO
     !$omp end parallel do

  ELSE

     !$omp parallel do private( i, offset, j )
     DO i = 1, dfft%nr1s * dfft%my_nr3p
        offset = dfft%my_nr2p * (i-1)
        DO j = 1, dfft%my_nr2p
           f_in( j + offset ) = f_aux( dfft%map_scatter_fw( j + offset ) )
        ENDDO
     ENDDO
     !$omp end parallel do

  END IF

END SUBROUTINE fft_scatter_xy

SUBROUTINE fwfft_pre_com( dfft, f, aux, comm_mem_send, comm_mem_recv, nr1s, ns )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(INOUT)  :: f( * ) !dfft%my_nr3p * dfft%nr2 * dfft%nr1 , * ) !batch_size )
  COMPLEX(DP), INTENT(OUT) :: comm_mem_send( * )
  COMPLEX(DP), INTENT(INOUT) :: comm_mem_recv( * )
  INTEGER, INTENT(IN) :: nr1s( * ), ns( * )
  COMPLEX(DP), INTENT(INOUT) :: aux( * ) !nr1s(dfft%mype2+1) * dfft%my_nr3p * dfft%nr2 * batch_size )

  INTEGER :: i, j, k, l, m
  INTEGER :: offset, ierr, ibatch

!------------------------------------------------------
!------------x-FFT Start-------------------------------

  CALL fft_1D( f, dfft%my_nr2p * dfft%my_nr3p, dfft%nr1, -2 )

!-------------x-FFT End--------------------------------
!------------------------------------------------------
!------------------------------------------------------
!------Forward xy-scatter Start------------------------

  CALL fft_scatter_xy( dfft, aux( 1:nr1s(dfft%mype2+1)*dfft%my_nr3p*dfft%nr2 ), f, -2 )

!-------Forward xy-scatter End-------------------------
!------------------------------------------------------

  CALL fwfft_y_portion( dfft, comm_mem_send, comm_mem_recv, aux( 1: nr1s(dfft%mype2+1)*dfft%my_nr3p*dfft%nr2 ), dfft%map_pcfw, nr1s, ns )

END SUBROUTINE fwfft_pre_com

SUBROUTINE fwfft_y_portion( dfft, comm_mem_send, comm_mem_recv, aux, map_pcfw, nr1s, ns )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(INOUT)  :: comm_mem_send( * )
  COMPLEX(DP), INTENT(INOUT)  :: comm_mem_recv( * )
  INTEGER, INTENT(IN) :: nr1s( * ), ns( * )
  COMPLEX(DP), INTENT(INOUT) :: aux( dfft%nr2 * nr1s(dfft%mype2+1) * dfft%my_nr3p , * ) !y_group_size )
  INTEGER, INTENT(IN) :: map_pcfw( * )

  INTEGER :: l, m, i, offset, j, k

!------------------------------------------------------
!------------y-FFT Start-------------------------------

  CALL fft_1D( aux( 1 : dfft%nr2 * nr1s(dfft%mype2+1) * dfft%my_nr3p , 1 : 1 ), nr1s(dfft%mype2+1) * dfft%my_nr3p, dfft%nr2, -2 )

!-------------y-FFT End--------------------------------
!------------------------------------------------------
!------------------------------------------------------
!---------Pre-Com-Copy Start---------------------------

  IF( .not. dfft%single_node ) THEN

     !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 1 ), ierr )

     !$omp parallel private( i, j, k, offset, l, m )
     DO l = 1, dfft%nodes_numb
        IF( dfft%non_blocking .and. l .eq. dfft%my_node+1 ) CYCLE
        DO m = 1, dfft%node_task_size
           i = (l-1)*dfft%node_task_size + m
           offset = ( dfft%my_node_rank + (i-1)*dfft%node_task_size ) * dfft%small_chunks
           !$omp do
           DO j = 1, ns( i )
              !$omp simd
              DO k = 1, dfft%my_nr3p
                 comm_mem_send( offset + (j-1)*dfft%nr3px + k ) = &
                 aux( dfft%map_pcfw( (i-1)*dfft%small_chunks + (j-1)*dfft%nr3px + k ), 1 )
              END DO
           END DO
           !$omp end do nowait
        END DO
     END DO
     !$omp end parallel

  END IF

  !$omp parallel private( i, j, k, offset, m )
  DO m = 1, dfft%node_task_size
     i = dfft%my_node*dfft%node_task_size + m
     offset = ( dfft%my_node_rank + (i-1)*dfft%node_task_size ) * dfft%small_chunks
     !$omp do
     DO j = 1, ns( i )
        !$omp simd
        DO k = 1, dfft%my_nr3p
           comm_mem_recv( offset + (j-1)*dfft%nr3px + k ) = &
           aux( dfft%map_pcfw( (i-1)*dfft%small_chunks + (j-1)*dfft%nr3px + k ), 1 )
        END DO
     END DO
     !$omp end do nowait
  END DO
  !$omp end parallel

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

!----------Pre-Com-Copy End----------------------------
!------------------------------------------------------

END SUBROUTINE fwfft_y_portion

SUBROUTINE fwfft_after_com( dfft, comm_mem_recv, f, ns )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ns( * )
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( * )
  COMPLEX(DP), INTENT(INOUT)  :: f ( dfft%nr3 , * ) !ns(dfft%mype+1)*z_group_size )

  INTEGER :: j, l, k, i
  INTEGER :: offset, kfrom, ierr

!------------------------------------------------------
!--------After-Com-Copy Start--------------------------

  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )
 
  !$omp parallel private( kfrom, j, i, k, offset, l )
  DO j = 1, dfft%nodes_numb
     DO l = 1, dfft%node_task_size
        offset = ( dfft%my_node_rank*dfft%node_task_size + (l-1) ) * dfft%small_chunks + (j-1) * dfft%big_chunks
        !$omp do
        DO k = 1, ns(dfft%mype3+1)
           kfrom = offset + dfft%nr3px * (k-1)
           !$omp simd
           DO i = 1, dfft%nr3p( (j-1)*dfft%node_task_size + l )
              f( dfft%nr3p_offset( (j-1)*dfft%node_task_size + l ) + i, k ) = comm_mem_recv( kfrom + i ) * dfft%tscale
           ENDDO
        ENDDO
        !$omp end do
     ENDDO
  ENDDO
  !$omp end parallel

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

!---------After-Com-Copy End---------------------------
!------------------------------------------------------
!------------------------------------------------------
!------------z-FFT Start-------------------------------

  CALL fft_1D( f, ns(dfft%mype+1), dfft%nr3, -2 )

!-------------z-FFT End--------------------------------
!------------------------------------------------------

END SUBROUTINE fwfft_after_com

!=----------------------------------------------------------------------=
END MODULE fftpw_batchingSingle
!=----------------------------------------------------------------------=
