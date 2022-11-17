!=----------------------------------------------------------------------=
MODULE fftpw_batching
!=----------------------------------------------------------------------=

  USE cppt,                                     ONLY: hg
  USE density_utils,                            ONLY: build_density_sum
  USE elct,                                     ONLY: crge
  USE fftpw_legacy_routines,                    ONLY: fft_1D,&
                                                      fft_scatter_xy,&
                                                      fft_scatter_xy_t,&
                                                      fft_1d_t
  USE fftpw_param
  USE fftpw_types,                              ONLY: PW_fft_type_descriptor
  USE kinds,                                    ONLY: real_8
  USE system,                                   ONLY: parm
  USE iso_fortran_env
  IMPLICIT NONE
  PRIVATE


  PUBLIC :: Prepare_Psi_overlapp
  PUBLIC :: invfft_pre_com
  PUBLIC :: fft_com
  PUBLIC :: invfft_after_com
  PUBLIC :: Apply_V
  PUBLIC :: Build_CD
  PUBLIC :: fwfft_pre_com
  PUBLIC :: fwfft_after_com
  PUBLIC :: Accumulate_Psi_overlapp

CONTAINS

SUBROUTINE Prepare_Psi_overlapp( dfft, psi, aux, ngms, z_group_size, ns, last )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: ns(:)
  INTEGER, INTENT(IN) :: ngms, z_group_size
  LOGICAL, INTENT(IN) :: last
  COMPLEX(DP), INTENT(IN)  :: psi (:,:)
  COMPLEX(DP), INTENT(OUT)  :: aux ( dfft%nr3 , ns(dfft%mype+1)*z_group_size )

  INTEGER :: j, i, iter
  INTEGER :: offset, offset2

  INTEGER(INT64) :: time(2), cr

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!----------Prepare_Psi Start---------------------------

  !$omp parallel do private( i, j, iter, offset, offset2 )
  DO i = 1, ns( dfft%mype+1 ) * z_group_size
     iter = mod( i-1, ns(dfft%mype+1) ) + 1
     offset  = ( iter - 1 ) * dfft%nr3
     offset2 = 2 * ( ( (i-1) / ns(dfft%mype+1) ) + 1 )

     !$omp simd
     DO j = 1, dfft%prep_map(1,iter)-1
        aux( j, i ) = conjg( psi( dfft%nlm_r( offset + j ), offset2 - 1 ) - (0.0d0,1.0d0) * psi( dfft%nlm_r( offset + j ), offset2 ) )
     ENDDO
     !$omp simd
     DO j = 1, dfft%prep_map(2,iter)-1
        aux( j, i ) = psi( dfft%nl_r( offset + j ), offset2 - 1 ) + (0.0d0,1.0d0) * psi( dfft%nl_r( offset + j ), offset2 )
     ENDDO

     !$omp simd
     DO j = dfft%prep_map(3,iter), dfft%prep_map(4,iter)
        aux( j, i ) = (0.d0, 0.d0)
     ENDDO

     !$omp simd
     DO j = dfft%prep_map(5,iter)+1, dfft%nr3
        aux( j, i ) = psi( dfft%nl_r( offset + j ), offset2 - 1 ) + (0.0d0,1.0d0) * psi( dfft%nl_r( offset + j ), offset2 )
     ENDDO
     !$omp simd
     DO j = dfft%prep_map(6,iter)+1, dfft%nr3
        aux( j, i ) = conjg( psi( dfft%nlm_r( offset + j ), offset2 - 1 ) - (0.0d0,1.0d0) * psi( dfft%nlm_r( offset + j ), offset2 ) )
     ENDDO

  ENDDO
  !$omp end parallel do

! Waiting for case where the number of states is uneven
!  !$omp parallel do private( i, j, iter, offset, offset2 )
!  DO i = 1, ns( dfft%mype+1 ) * z_group_size
!     iter = mod( i-1, ns(dfft%mype+1) ) + 1
!     offset  = ( iter - 1 ) * dfft%nr3
!     offset2 = 2 * ( ( (i-1) / ns(dfft%mype+1) ) + 1 )
!     !$omp simd
!     DO j = 1, dfft%zero_prep_start(iter,1)-1
!        aux( j, i ) = psi( dfft%nl_r( offset + j ), offset2 - 1 )
!     ENDDO
!     !$omp simd
!     DO j = dfft%zero_prep_start(iter,2), dfft%zero_prep_end(iter,2)
!        aux( j, i ) = (0.d0, 0.d0)
!     ENDDO
!     !$omp simd
!     DO j = dfft%zero_prep_end(iter,1)+1, dfft%nr3
!        aux( j, i ) = psi( dfft%nl_r( offset + j ), offset2 - 1 )
!     ENDDO
!  ENDDO
!  !$omp end parallel do
!  
!  !$omp parallel do private( i, j, iter, offset, offset2 )
!  DO i = 1, ns( dfft%mype+1 ) * z_group_size
!     iter = mod( i-1, ns(dfft%mype+1) ) + 1
!     offset  = ( iter - 1 ) * dfft%nr3
!     offset2 = 2 * ( ( (i-1) / ns(dfft%mype+1) ) + 1 )
!     !$omp simd
!     DO j = 1, dfft%zero_prep_start(iter,3)-1
!        aux( j, i ) = conjg( psi( dfft%nlm_r( offset + j ), offset2 - 1 ) )
!     ENDDO
!     !$omp simd
!     DO j = dfft%zero_prep_end(iter,3)+1, dfft%nr3
!        aux( j, i ) = conjg( psi( dfft%nlm_r( offset + j ), offset2 - 1 ) )
!     ENDDO
!  ENDDO
!  !$omp end parallel do

!----------Prepare_Psi End-----------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )

!  dfft%time_adding( 1 ) = dfft%time_adding( 1 ) + ( time(2) - time(1) )
!  CALL SYSTEM_CLOCK( count_rate = cr )
!  write(6,*) "ACTUAL", REAL( dfft%time_adding( 1 ) / REAL( cr ) )

END SUBROUTINE Prepare_Psi_overlapp

SUBROUTINE invfft_pre_com( dfft, aux, comm_mem_send, comm_mem_recv, iset, batch_size, z_group_size, ns )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iset, batch_size, z_group_size
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: ns(:)
  COMPLEX(DP), INTENT(INOUT) :: comm_mem_send( : ), comm_mem_recv( : )
  COMPLEX(DP), INTENT(INOUT)  :: aux ( dfft%nr3 , ns(dfft%mype+1)*z_group_size )

  INTEGER :: l, m, j, k, i
  INTEGER :: offset, kdest, ierr

  INTEGER(INT64) :: time(3)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!------------z-FFT Start-------------------------------

  CALL fft_1D_t( aux, ns(dfft%mype+1)*z_group_size, dfft%nr3, 2 )

!-------------z-FFT End--------------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )
!------------------------------------------------------
!---------Pre-Com-Copy Start---------------------------

  IF( .not. dfft%single_node ) THEN

     !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 1 ), ierr )

     !$omp parallel private( l, m, j, offset, k, kdest, i )
     DO l = 1, dfft%nodes_numb
        IF( dfft%non_blocking .and. l .eq. dfft%my_node+1 ) CYCLE 
        DO m = 1, dfft%node_task_size
           j = (l-1)*dfft%node_task_size + m
           offset = ( dfft%my_node_rank + (j-1)*dfft%node_task_size ) * dfft%small_chunks + ( (l-1)*(batch_size-1) + (iset-1)*dfft%z_group_size_save ) * dfft%big_chunks
           !$omp do
           DO k = 1, ns(dfft%mype+1) * z_group_size
              kdest = offset + dfft%nr3px * mod( (k-1), ns(dfft%mype+1) ) + ( (k-1) / ns(dfft%mype+1) ) * dfft%big_chunks
              !omp simd
              DO i = 1, dfft%nr3p( j )
                 comm_mem_send( kdest + i ) = aux( i + dfft%nr3p_offset( j ), k )
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

     !$omp parallel private( m, j, offset, k, kdest, i )
     DO m = 1, dfft%node_task_size
        j = dfft%my_node*dfft%node_task_size + m
        offset = ( dfft%my_node_rank + (j-1)*dfft%node_task_size ) * dfft%small_chunks + ( dfft%my_node*(batch_size-1) + (iset-1)*dfft%z_group_size_save ) * dfft%big_chunks
        !$omp do
        DO k = 1, ns(dfft%mype+1) * z_group_size 
           kdest = offset + dfft%nr3px * mod( k-1, ns(dfft%mype+1) ) + ( (k-1) / ns(dfft%mype+1) ) * dfft%big_chunks
           !omp simd
           DO i = 1, dfft%nr3p( j )
              comm_mem_recv( kdest + i ) = aux( i + dfft%nr3p_offset( j ), k )
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
  dfft%counter(2) = dfft%counter(2) + 1

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

SUBROUTINE invfft_after_com( dfft, f, comm_mem_recv, aux, map_acinv, map_acinv_rem, y_set_size, batch_size, nr1s )
  USE mpi_f08
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: y_set_size, batch_size
  COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( : )
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(OUT) :: f( dfft%nnr, batch_size )
  INTEGER, INTENT(IN) :: nr1s(:)
  COMPLEX(DP), INTENT(INOUT) :: aux( nr1s(dfft%mype2+1) * dfft%my_nr3p * dfft%nr2 * batch_size )
  INTEGER, INTENT(IN) :: map_acinv(:), map_acinv_rem(:)

  INTEGER :: iset, y_group_size, ibatch

  INTEGER(INT64) :: time(6)

  DO iset = 1, y_set_size

     !Too slow? lets wait and see!
     IF( iset .ne. y_set_size .or. iset .eq. 1 .or. y_group_size * y_set_size .eq. batch_size ) THEN
        IF( mod( batch_size, y_set_size ) .eq. 0 ) THEN
           y_group_size = batch_size / y_set_size
        ELSE
           y_group_size = batch_size / y_set_size + 1
        ENDIF
        dfft%y_group_size_save = y_group_size
     ELSE
        y_group_size = batch_size - (y_set_size-1) * y_group_size
     END IF

     CALL invfft_y_portion( dfft, comm_mem_recv, aux( 1+(iset-1)*nr1s(dfft%mype2+1)*dfft%my_nr3p*dfft%nr2*dfft%y_group_size_save : &
                            nr1s(dfft%mype2+1)*dfft%my_nr3p*dfft%nr2*y_group_size+(iset-1)*nr1s(dfft%mype2+1)*dfft%my_nr3p*dfft%nr2*dfft%y_group_size_save ), &
                            map_acinv, map_acinv_rem, y_group_size, nr1s, iset )

  ENDDO

  CALL SYSTEM_CLOCK( time(3) )
!------------------------------------------------------
!------Forward xy-scatter Start------------------------

  DO ibatch = 1, batch_size

     CALL fft_scatter_xy_t( dfft, aux( 1+nr1s(dfft%mype2+1)*dfft%my_nr3p*dfft%nr2*(ibatch-1) : nr1s(dfft%mype2+1)*dfft%my_nr3p*dfft%nr2*ibatch ), f(:,ibatch), dfft%nnr, 2 )

  ENDDO

!-------Forward xy-scatter End-------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(4) )
!------------------------------------------------------
!------------x-FFT Start-------------------------------

  CALL fft_1D_t( f( 1 : dfft%my_nr3p * dfft%nr2 * dfft%nr1 , : ), dfft%my_nr2p * dfft%my_nr3p * batch_size, dfft%nr1, 2 )

!-------------x-FFT End--------------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(5) )

  dfft%time_adding( 6 ) = dfft%time_adding( 6 ) + ( time(4) - time(3) )
  dfft%time_adding( 7 ) = dfft%time_adding( 7 ) + ( time(5) - time(4) )
  dfft%counter(3) = dfft%counter(3) + 1

END SUBROUTINE invfft_after_com

SUBROUTINE invfft_y_portion( dfft, comm_mem_recv, aux, map_acinv, map_acinv_rem, y_group_size, nr1s, iset )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: y_group_size, iset
  COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( : )
  INTEGER, INTENT(IN) :: nr1s(:)
  COMPLEX(DP), INTENT(INOUT) :: aux( dfft%nr2 , nr1s(dfft%mype2+1) * dfft%my_nr3p *  y_group_size )
  INTEGER, INTENT(IN) :: map_acinv(:), map_acinv_rem(:)

  INTEGER :: i, k, offset, jter

  INTEGER(INT64) :: time(6)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!--------After-Com-Copy Start--------------------------

  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )
  
  IF( dfft%rem ) THEN

     !$omp parallel do private( i, jter, k, offset )
     DO i = 1, dfft%my_nr1p * dfft%my_nr3p * y_group_size
        jter = mod( i-1, dfft%my_nr1p ) + 1
        offset = ( mod( i-1, dfft%my_nr1p * dfft%my_nr3p ) + ( ( (i-1) / ( dfft%my_nr1p * dfft%my_nr3p ) ) + (iset-1)*dfft%y_group_size_save ) * dfft%my_nr3p * dfft%my_nr1p ) * dfft%nr2
        !omp simd
        DO k = 1, dfft%zero_acinv_start( jter ) - 1
           aux( k, i ) = comm_mem_recv( map_acinv_rem( offset + k ) )
        END DO
        !omp simd
        DO k = dfft%zero_acinv_start( jter ), dfft%zero_acinv_end( jter )
           aux( k, i ) = (0.0_DP,0.0_DP)
        END DO
        !omp simd
        DO k = dfft%zero_acinv_end( jter ) + 1, dfft%nr2
           aux( k, i ) = comm_mem_recv( map_acinv_rem( offset + k ) )
        END DO
     END DO
     !$omp end parallel do

  ELSE

     !$omp parallel do private( i, jter, k, offset )
     DO i = 1, dfft%my_nr1p * dfft%my_nr3p * y_group_size
        jter = mod( i-1, dfft%my_nr1p ) + 1
        offset = ( mod( i-1, dfft%my_nr1p * dfft%my_nr3p ) + ( ( (i-1) / ( dfft%my_nr1p * dfft%my_nr3p ) ) + (iset-1)*dfft%y_group_size_save ) * dfft%my_nr3p * dfft%my_nr1p ) * dfft%nr2
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

  END IF

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

!---------After-Com-Copy End---------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )
!------------------------------------------------------
!------------y-FFT Start-------------------------------

  CALL fft_1D_t( aux, nr1s(dfft%mype2+1) * dfft%my_nr3p * y_group_size, dfft%nr2, 2 )

!-------------y-FFT End--------------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(3) )

  dfft%time_adding( 4 ) = dfft%time_adding( 4 ) + ( time(2) - time(1) )
  dfft%time_adding( 5 ) = dfft%time_adding( 5 ) + ( time(3) - time(2) )

END SUBROUTINE invfft_y_portion

SUBROUTINE Apply_V( dfft, f, v, ibatch ) 
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(INOUT) :: f( : )
  REAL(DP), INTENT(IN) :: v( : )
  INTEGER, INTENT(IN)  :: ibatch

  INTEGER :: j

  INTEGER(INT64) :: time(2)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!-----------Apply V Start------------------------------

  IF( .not. dfft%rsactive ) THEN
     !$omp parallel do private ( j )
     DO j = 1, dfft%nnr
        f( j )= - f( j ) * v( j )
     END DO
     !$omp end parallel do
  ELSE
     !$omp parallel do private ( j )
     DO j = 1, dfft%nnr
        f( j )= dfft%wfn_keep( j, ibatch ) * v(j)
     END DO
     !$omp end parallel do
  END IF

!------------Apply V End-------------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )

  dfft%time_adding( 8 ) = dfft%time_adding( 8 ) + ( time(2) - time(1) )
  dfft%counter(4) = dfft%counter(4) + 1

END SUBROUTINE Apply_V

SUBROUTINE Build_CD( dfft, f, rhoe, num ) 
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(IN) :: f( : )
  REAL(real_8), INTENT(OUT) :: rhoe( : )
  INTEGER, INTENT(IN)  :: num

  REAL(real_8) :: coef3, coef4

  INTEGER(INT64) :: time(2)

!------------------------------------------------------
!-----------Build CD Start-----------------------------

  ! Compute the charge density from the wave functions
  ! in real space
  coef3=crge%f(num,1)/parm%omega
!  IF (is2.GT.nstate) THEN
!     coef4=0.0_real_8
!  ELSE
     coef4=crge%f(num+1,1)/parm%omega
!  ENDIF
  CALL build_density_sum(coef3,coef4,f,rhoe,dfft%nnr)

!------------Build CD End------------------------------
!------------------------------------------------------

END SUBROUTINE Build_CD

SUBROUTINE fwfft_pre_com( dfft, f, comm_mem_send, comm_mem_recv, ibatch, batch_size, nr1s, ns )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ibatch, batch_size
  COMPLEX(DP), INTENT(INOUT)  :: f( : )
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(OUT) :: comm_mem_send( : )
  COMPLEX(DP), INTENT(INOUT) :: comm_mem_recv( : )
  INTEGER, INTENT(IN) :: nr1s(:), ns(:)

  INTEGER :: i, j, k, l, m
  INTEGER :: offset, ierr

  INTEGER(INT64) :: time(5)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!------------x-FFT Start-------------------------------

  CALL fft_1D( f, dfft%my_nr2p * dfft%my_nr3p, dfft%nr1, -2 )

!-------------x-FFT End--------------------------------
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

  CALL fft_1D( dfft%aux, nr1s(dfft%mype2+1) * dfft%my_nr3p, dfft%nr2, -2 )

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
           DO j = 1, ns( i )
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
        DO j = 1, ns( i )
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
        DO j = 1, ns( i )
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
  dfft%counter(5) = dfft%counter(5) + 1

END SUBROUTINE fwfft_pre_com

SUBROUTINE fwfft_after_com( dfft, comm_mem_recv, aux, iset, batch_size, z_group_size, ns )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iset, batch_size, z_group_size
  INTEGER, INTENT(IN) :: ns(:)
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( : )
  COMPLEX(DP), INTENT(INOUT)  :: aux ( dfft%nr3 , ns(dfft%mype+1)*z_group_size )

  INTEGER :: j, l, k, i
  INTEGER :: offset, kfrom, ierr

  INTEGER(INT64) :: time(3)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!--------After-Com-Copy Start--------------------------

  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )
 
  !$omp parallel private( kfrom, j, i, k, offset, l )
  DO j = 1, dfft%nodes_numb
     DO l = 1, dfft%node_task_size
        offset = ( dfft%my_node_rank*dfft%node_task_size + (l-1) ) * dfft%small_chunks + ( (j-1)*batch_size + (iset-1)*dfft%z_group_size_save ) * dfft%big_chunks
        !$omp do
        DO k = 1, ns(dfft%mype3+1) * z_group_size
           kfrom = offset + dfft%nr3px * mod( k-1, ns(dfft%mype+1) ) + ( (k-1) / ns(dfft%mype+1) ) * dfft%big_chunks
           !$omp simd
           DO i = 1, dfft%nr3p( (j-1)*dfft%node_task_size + l )
              aux( dfft%nr3p_offset( (j-1)*dfft%node_task_size + l ) + i, k ) = comm_mem_recv( kfrom + i )
           ENDDO
        ENDDO
        !$omp end do
     ENDDO
  ENDDO
  !$omp end parallel

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

!---------After-Com-Copy End---------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )
!------------------------------------------------------
!------------z-FFT Start-------------------------------

  CALL fft_1D_t( aux, ns(dfft%mype+1)*z_group_size, dfft%nr3, -2 )

!-------------z-FFT End--------------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(3) )

  dfft%time_adding( 13 ) = dfft%time_adding( 13 ) + ( time(2) - time(1) )
  dfft%time_adding( 14 ) = dfft%time_adding( 14 ) + ( time(3) - time(2) )
  dfft%counter(6) = dfft%counter(6) + 1

END SUBROUTINE fwfft_after_com

SUBROUTINE Accumulate_Psi_overlapp( dfft, aux, hpsi, ngms, z_group_size, last, ns, psi )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: z_group_size, ngms
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: ns(:)
  COMPLEX(DP), INTENT(IN), OPTIONAL :: psi(:,:)
  COMPLEX(DP), INTENT(INOUT) :: hpsi(:,:)
  COMPLEX(DP), INTENT(IN)  :: aux ( dfft%nr3 * ns(dfft%mype+1), z_group_size )
  LOGICAL, INTENT(IN) :: last

  COMPLEX(DP) :: fp, fm
  INTEGER :: j, l, k, i, igroup
  INTEGER :: offset

  INTEGER(INT64) :: time(2)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!--------Accumulate_Psi Start--------------------------
   
  !Somehow slower in the beginning than previous  implementation... needs investigation
  !$omp parallel private( j, igroup, fp, fm )
  DO igroup = 1, z_group_size
     !$omp do
     DO j = 1, ngms
        fp = ( aux( dfft%nl(j), igroup ) + aux( dfft%nlm(j), igroup ) ) * (- dfft%tscale )
        fm = ( aux( dfft%nl(j), igroup ) - aux( dfft%nlm(j), igroup ) ) * (- dfft%tscale )
        hpsi ( j, (2*igroup)-1 ) = -1 * ((parm%tpiba2*dfft%gg_pw(j))*psi( j, (2*igroup)-1 ) + cmplx(  dble(fp) , aimag(fm), KIND=DP ) )
        hpsi ( j, (2*igroup)   ) = -1 * ((parm%tpiba2*dfft%gg_pw(j))*psi( j, (2*igroup)   ) + cmplx(  aimag(fp), -dble(fm), KIND=DP ) )
     END DO
     !$omp end do nowait
  ENDDO
  !$omp end parallel
     
! Waiting for case where the number of states is uneven (Probably faulty even then)
!        !$omp parallel do private( j )
!        DO j = 1, ngms
!           hpsi( j, (2*igroup)-1 ) = dfft%aux2( dfft%nl( j ) + offset ) * dfft%tscale
!        END DO
!        !$omp end parallel do

!---------Accumulate_Psi End---------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )

!  dfft%time_adding( 15 ) = dfft%time_adding( 15 ) + ( time(2) - time(1) )
  dfft%counter(7) = dfft%counter(7) + 1

END SUBROUTINE Accumulate_Psi_overlapp

!=----------------------------------------------------------------------=
END MODULE fftpw_batching
!=----------------------------------------------------------------------=
