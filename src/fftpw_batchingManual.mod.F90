!=----------------------------------------------------------------------=
MODULE fftpw_batchingManual
!=----------------------------------------------------------------------=

  USE fftpw_legacy_routines
  USE fftpw_param
  USE fftpw_types,                              ONLY: PW_fft_type_descriptor
  USE kinds,                                    ONLY: real_8
  USE system,                                   ONLY: parm
  USE iso_fortran_env
  IMPLICIT NONE
  PRIVATE


  PUBLIC :: Prepare_Psi_Man
  PUBLIC :: fft_com_Man
  PUBLIC :: Accumulate_Psi_Man
  PUBLIC :: invfft_z_section_Man
  PUBLIC :: invfft_y_section_Man
  PUBLIC :: invfft_x_section_Man
  PUBLIC :: fwfft_x_section_Man
  PUBLIC :: fwfft_y_section_Man
  PUBLIC :: fwfft_z_section_Man

CONTAINS

SUBROUTINE Prepare_Psi_Man( dfft, psi, aux, ngms, remswitch, mythread, ns, last, priv, ip, jp )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: ns( * )
  INTEGER, INTENT(IN) :: ngms, remswitch, mythread
  INTEGER, INTENT(INOUT) :: ip, jp
  LOGICAL, INTENT(IN) :: last
  INTEGER, INTENT(INOUT) :: priv(:)
  COMPLEX(DP), INTENT(IN)  :: psi ( : , : )
  COMPLEX(DP), INTENT(OUT)  :: aux ( dfft%nr3 , * ) !ns(dfft%mype+1)*z_group_size )

  INTEGER :: j, i, iter
  INTEGER :: offset, offset2

  INTEGER(INT64) :: time(2), cr

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!----------Prepare_Psi Start---------------------------

  ! dfft%mype+1 in dfft%thread_z_start eigentlich (dfft%my_node-1)*dfft%node_task_size+dfft%my_node_rank+1

  DO i = dfft%thread_z_start( mythread+1, remswitch, dfft%mype+1 ), dfft%thread_z_end( mythread+1, remswitch, dfft%mype+1 )
!  DO i = dfft%thread_z_start( mythread+1, remswitch, dfft%my_node_rank+1 ), dfft%thread_z_end( mythread+1, remswitch, dfft%my_node_rank+1 )
     iter = mod( i-1, ns(dfft%mype+1) ) + 1
     offset  = ( iter - 1 ) * dfft%nr3
     offset2 = 2 * ( ( (i-1) / ns(dfft%mype+1) ) + 1 )

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
  CALL SYSTEM_CLOCK( time(2) )

!  dfft%time_adding( 1 ) = dfft%time_adding( 1 ) + ( time(2) - time(1) )
!  CALL SYSTEM_CLOCK( count_rate = cr )
!  write(6,*) "ACTUAL", REAL( dfft%time_adding( 1 ) / REAL( cr ) )

END SUBROUTINE Prepare_Psi_Man

SUBROUTINE fft_com_Man( dfft, comm_mem_send, comm_mem_recv, sendsize, intra_me, inter_node_comm, nodes_numb, inter_me, non_blocking, work_buffer )
  IMPLICIT NONE

  INTEGER, INTENT(IN)                            :: sendsize, intra_me, nodes_numb, inter_me, work_buffer
  TYPE(MPI_COMM), INTENT(IN)                     :: inter_node_comm
  LOGICAL, INTENT(IN)                            :: non_blocking
  TYPE(PW_fft_type_descriptor), INTENT(INOUT)       :: dfft
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

        IF( work_buffer .lt. 0 .or. .false. ) THEN

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

        ELSE

           IF( sendsize .eq. dfft%sendsize_save ) THEN
   
              CALL MPI_STARTALL( nodes_numb-1, dfft%send_handle(:,work_buffer) )
              CALL MPI_STARTALL( nodes_numb-1, dfft%recv_handle(:,work_buffer) )
   
              CALL MPI_WAITALL( nodes_numb-1, dfft%send_handle(:,work_buffer), MPI_STATUSES_IGNORE, ierr )
              CALL MPI_WAITALL( nodes_numb-1, dfft%recv_handle(:,work_buffer), MPI_STATUSES_IGNORE, ierr )
     
           ELSE
   
              CALL MPI_STARTALL( nodes_numb-1, dfft%send_handle_rem(:,work_buffer) )
              CALL MPI_STARTALL( nodes_numb-1, dfft%recv_handle_rem(:,work_buffer) )
   
              CALL MPI_WAITALL( nodes_numb-1, dfft%send_handle_rem(:,work_buffer), MPI_STATUSES_IGNORE, ierr )
              CALL MPI_WAITALL( nodes_numb-1, dfft%recv_handle_rem(:,work_buffer), MPI_STATUSES_IGNORE, ierr )
   
           END IF

        END IF

     END IF

     !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )
     !CALL mpi_win_unlock_all( dfft%mpi_window( 1 ), ierr )

  END IF

END SUBROUTINE fft_com_Man

SUBROUTINE Accumulate_Psi_Man( dfft, aux, hpsi, ngms, z_group_size, last, ns, psi, mythread )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: z_group_size, ngms, mythread
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: ns( * )
  COMPLEX(DP), INTENT(IN), OPTIONAL :: psi( : , : )
  COMPLEX(DP), INTENT(INOUT) :: hpsi( : , : )
  COMPLEX(DP), INTENT(IN)  :: aux ( dfft%nr3 * ns(dfft%mype+1), * ) !z_group_size )
  LOGICAL, INTENT(IN) :: last

  COMPLEX(DP) :: fp, fm
  INTEGER :: j, l, k, i, igroup
  INTEGER :: offset

  INTEGER(INT64) :: time(2)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!--------Accumulate_Psi Start--------------------------
   
  !Somehow slower in the beginning than previous  implementation... needs investigation
  DO igroup = 1, z_group_size
     DO j = dfft%thread_ngms_start( mythread+1 ), dfft%thread_ngms_end( mythread+1 )
        fp = ( aux( dfft%nl(j), igroup ) + aux( dfft%nlm(j), igroup ) ) * (- dfft%tscale )
        fm = ( aux( dfft%nl(j), igroup ) - aux( dfft%nlm(j), igroup ) ) * (- dfft%tscale )
        hpsi ( j, (2*igroup)-1 ) = -1 * ((parm%tpiba2*dfft%gg_pw(j))*psi( j, (2*igroup)-1 ) + cmplx(  dble(fp) , aimag(fm), KIND=DP ) )
        hpsi ( j, (2*igroup)   ) = -1 * ((parm%tpiba2*dfft%gg_pw(j))*psi( j, (2*igroup)   ) + cmplx(  aimag(fp), -dble(fm), KIND=DP ) )
     END DO
  ENDDO
     
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

END SUBROUTINE Accumulate_Psi_Man

SUBROUTINE invfft_z_section_Man( dfft, aux, comm_mem_send, comm_mem_recv, iset, batch_size, remswitch, mythread, ns )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iset, batch_size, remswitch, mythread
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: ns( * )
  COMPLEX(DP), INTENT(INOUT) :: comm_mem_send( * ), comm_mem_recv( * )
  COMPLEX(DP), INTENT(INOUT)  :: aux ( dfft%nr3 , * ) !ns(dfft%mype+1)*z_group_size )

  INTEGER :: l, m, j, k, i
  INTEGER :: offset, kdest, ierr

  INTEGER(INT64) :: time(3)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!------------z-FFT Start-------------------------------

  CALL fft_noOMP_1D( aux( : , dfft%thread_z_start( mythread+1, remswitch, dfft%mype+1 ) : dfft%thread_z_end( mythread+1, remswitch, dfft%mype+1 ) ) , &
                     dfft%thread_z_sticks(:,:,dfft%mype+1), dfft%nr3, remswitch, mythread, dfft%nthreads, 2 )

!-------------z-FFT End--------------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )
!------------------------------------------------------
!---------Pre-Com-Copy Start---------------------------

!  DO i = dfft%thread_z_start( mythread+1 ), dfft%thread_z_end( mythread+1 )

  IF( .not. dfft%single_node ) THEN

     !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 1 ), ierr )

     DO l = 1, dfft%nodes_numb
        IF( dfft%non_blocking .and. l .eq. dfft%my_node+1 ) CYCLE 
        DO m = 1, dfft%node_task_size
           j = (l-1)*dfft%node_task_size + m
           offset = ( dfft%my_node_rank + (j-1)*dfft%node_task_size ) * dfft%small_chunks + ( (l-1)*(batch_size-1) + (iset-1)*batch_size ) * dfft%big_chunks
           DO k = dfft%thread_z_start( mythread+1, remswitch, dfft%mype+1 ), dfft%thread_z_end( mythread+1, remswitch, dfft%mype+1 )
              kdest = offset + dfft%nr3px * mod( (k-1), ns(dfft%mype+1) ) + ( (k-1) / ns(dfft%mype+1) ) * dfft%big_chunks
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
        offset = ( dfft%my_node_rank + (j-1)*dfft%node_task_size ) * dfft%small_chunks + ( dfft%my_node*(batch_size-1) + (iset-1)*batch_size ) * dfft%big_chunks
        DO k = dfft%thread_z_start( mythread+1, remswitch, dfft%mype+1 ), dfft%thread_z_end( mythread+1, remswitch, dfft%mype+1 )
           kdest = offset + dfft%nr3px * mod( k-1, ns(dfft%mype+1) ) + ( (k-1) / ns(dfft%mype+1) ) * dfft%big_chunks
           DO i = 1, dfft%nr3p( j )
              comm_mem_recv( kdest + i ) = aux( i + dfft%nr3p_offset( j ), k )
           ENDDO
        ENDDO
     ENDDO

     !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

  END IF

!----------Pre-Com-Copy End----------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(3) )

  dfft%time_adding( 2 ) = dfft%time_adding( 2 ) + ( time(2) - time(1) )
  dfft%time_adding( 3 ) = dfft%time_adding( 3 ) + ( time(3) - time(2) )
  dfft%counter(2) = dfft%counter(2) + 1

END SUBROUTINE invfft_z_section_Man 

SUBROUTINE invfft_y_section_Man( dfft, aux, comm_mem_recv, aux2_r, map_acinv, map_acinv_rem, nr1s, y_group_size, yset, remswitch, mythread )
  !$ USE omp_lib
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: y_group_size, yset, remswitch, mythread
  COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( * )
  INTEGER, INTENT(IN) :: nr1s( * )
  COMPLEX(DP), INTENT(INOUT) :: aux ( dfft%my_nr3p * dfft%nr2 * dfft%nr1 , * ) !y_group_size
  COMPLEX(DP), INTENT(INOUT) :: aux2_r( : , : )
  INTEGER, INTENT(IN) :: map_acinv( * ), map_acinv_rem( * )

  INTEGER :: i, k, offset, iter, igroup

  INTEGER(INT64) :: time(5)

  Call First_Part_y_section( aux2_r )

  !$OMP barrier
 
  Call Second_Part_y_section( aux2_r )

  CONTAINS

    SUBROUTINE First_Part_y_section( aux2 )
      
      Implicit NONE
      COMPLEX(DP), INTENT(INOUT) :: aux2( dfft%nr2 , * ) !nr1s(dfft%mype2+1) * dfft%my_nr3p *  y_group_size )


        CALL SYSTEM_CLOCK( time(1) )
      !------------------------------------------------------
      !--------After-Com-Copy Start--------------------------
      
        !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )
        
        IF( dfft%rem ) THEN
      
           DO i = dfft%thread_y_start( mythread+1, remswitch ), dfft%thread_y_end( mythread+1, remswitch )
              iter = mod( i-1, dfft%my_nr1p ) + 1
              offset = ( mod( i-1, dfft%my_nr1p * dfft%my_nr3p ) + ( ( (i-1) / ( dfft%my_nr1p * dfft%my_nr3p ) ) + (yset-1)*y_group_size ) * dfft%my_nr3p * dfft%my_nr1p ) * dfft%nr2
              DO k = 1, dfft%zero_acinv_start( iter ) - 1
                 aux2( k, i ) = comm_mem_recv( map_acinv_rem( offset + k ) )
              END DO
              DO k = dfft%zero_acinv_start( iter ), dfft%zero_acinv_end( iter )
                 aux2( k, i ) = (0.0_DP,0.0_DP)
              END DO
              DO k = dfft%zero_acinv_end( iter ) + 1, dfft%nr2
                 aux2( k, i ) = comm_mem_recv( map_acinv_rem( offset + k ) )
              END DO
           END DO
      
        ELSE
      
           DO i = dfft%thread_y_start( mythread+1, remswitch ), dfft%thread_y_end( mythread+1, remswitch )
              iter = mod( i-1, dfft%my_nr1p ) + 1
              offset = ( mod( i-1, dfft%my_nr1p * dfft%my_nr3p ) + ( ( (i-1) / ( dfft%my_nr1p * dfft%my_nr3p ) ) + (yset-1)*y_group_size ) * dfft%my_nr3p * dfft%my_nr1p ) * dfft%nr2
              DO k = 1, dfft%zero_acinv_start( iter ) - 1
                 aux2( k, i ) = comm_mem_recv( map_acinv( offset + k ) )
!                 IF( mythread .eq. 0 ) THEN
!                 write(6,*) mythread, k, offset, map_acinv( offset + k )
!                 write(6,*) mythread, comm_mem_recv( map_acinv( offset + k ) )
!                 END IF
              END DO
              DO k = dfft%zero_acinv_start( iter ), dfft%zero_acinv_end( iter )
                 aux2( k, i ) = (0.0_DP,0.0_DP)
              END DO
              DO k = dfft%zero_acinv_end( iter ) + 1, dfft%nr2
                 aux2( k, i ) = comm_mem_recv( map_acinv( offset + k ) )
              END DO
           END DO
      
        END IF
      
        !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )
      
      !---------After-Com-Copy End---------------------------
      !------------------------------------------------------
        CALL SYSTEM_CLOCK( time(2) )
      !------------------------------------------------------
      !------------y-FFT Start-------------------------------
      
        CALL fft_noOMP_1D( aux2( : , dfft%thread_y_start( mythread+1, remswitch ) : dfft%thread_y_end( mythread+1, remswitch ) ), dfft%thread_y_sticks, dfft%nr2, remswitch, mythread, dfft%nthreads, 2 )
      
      !-------------y-FFT End--------------------------------
      !------------------------------------------------------
        CALL SYSTEM_CLOCK( time(3) )

        dfft%time_adding( 4 ) = dfft%time_adding( 4 ) + ( time(2) - time(1) )
        dfft%time_adding( 5 ) = dfft%time_adding( 5 ) + ( time(3) - time(2) )

    END SUBROUTINE First_Part_y_section

    SUBROUTINE Second_Part_y_section( aux2 )

      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT) :: aux2  ( dfft%my_nr1p * dfft%my_nr3p * dfft%nr2, * ) !y_group_size
        CALL SYSTEM_CLOCK( time(4) )
      !------------------------------------------------------
      !-------------yx-scatter-------------------------------
      
        DO i = dfft%thread_x_start( mythread+1, remswitch ), dfft%thread_x_end( mythread+1, remswitch )
           offset = mod( i-1, dfft%my_nr3p * dfft%my_nr2p ) * dfft%nr1
           igroup = ( ( (i-1) / ( dfft%my_nr3p * dfft%my_nr2p ) ) + 1 )
           DO k = 1, dfft%zero_scatter_start - 1
              aux( offset + k, igroup ) = aux2( dfft%map_scatter_inv( offset + k ), igroup )
           END DO
           DO k = dfft%zero_scatter_start, dfft%zero_scatter_end
              aux( offset + k, igroup ) = (0.0_DP, 0.0_DP)
           END DO
           DO k = dfft%zero_scatter_end + 1, dfft%nr1
              aux( offset + k, igroup ) = aux2( dfft%map_scatter_inv( offset + k ), igroup )
           END DO
        END DO
      
      !-------------yx-scatter-------------------------------
      !------------------------------------------------------
        CALL SYSTEM_CLOCK( time(5) )
      
        dfft%time_adding( 6 ) = dfft%time_adding( 6 ) + ( time(5) - time(4) )

    END SUBROUTINE Second_Part_y_section

END SUBROUTINE invfft_y_section_Man

SUBROUTINE invfft_x_section_Man( dfft, aux, remswitch, mythread )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: remswitch, mythread
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(INOUT) :: aux( dfft%nr1, * ) !dfft%my_nr3p * dfft%nr2 * x_group_size )

  INTEGER(INT64) :: time(2)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!------------x-FFT Start-------------------------------

  CALL fft_noOMP_1D( aux( : , dfft%thread_x_start( mythread+1, remswitch ) : dfft%thread_x_end( mythread+1, remswitch ) ), dfft%thread_x_sticks, dfft%nr1, remswitch, mythread, dfft%nthreads, 2 )

!-------------x-FFT End--------------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )

  dfft%time_adding( 7 ) = dfft%time_adding( 7 ) + ( time(2) - time(1) )

END SUBROUTINE invfft_x_section_Man

SUBROUTINE fwfft_x_section_Man( dfft, aux_r, aux2, nr1s, remswitch, mythread )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: remswitch, mythread
  INTEGER, INTENT(IN) :: nr1s( * )
  COMPLEX(DP), INTENT(INOUT)  :: aux_r( : , : ) !dfft%my_nr3p * dfft%nr2 * dfft%nr1 , * ) !x_group_size )
  COMPLEX(DP), INTENT(INOUT) :: aux2( nr1s(dfft%mype2+1) * dfft%my_nr3p * dfft%nr2 , * ) !x_group_size )

  INTEGER :: i, j, igroup, offset

  INTEGER(INT64) :: time(4)

  Call First_Part_x_section( aux_r )

  !$OMP barrier
 
  Call Second_Part_x_section( aux_r )

  CONTAINS

    SUBROUTINE First_Part_x_section( aux )

      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT) :: aux( dfft%nr1 , * ) !dfft%my_nr3p * dfft%nr2 * x_group_size )


        CALL SYSTEM_CLOCK( time(1) )
      !------------------------------------------------------
      !------------x-FFT Start-------------------------------
      
        CALL fft_noOMP_1D( aux( : , dfft%thread_x_start( mythread+1, remswitch ) : dfft%thread_x_end( mythread+1, remswitch ) ), dfft%thread_x_sticks, dfft%nr1, remswitch, mythread, dfft%nthreads, -2 )
      
      !-------------x-FFT End--------------------------------
      !------------------------------------------------------
        CALL SYSTEM_CLOCK( time(2) )

        dfft%time_adding( 9 ) = dfft%time_adding( 9 ) + ( time(2) - time(1) )

    END SUBROUTINE First_Part_x_section

    SUBROUTINE Second_Part_x_section( aux )

      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT)  :: aux( dfft%my_nr3p * dfft%nr2 * dfft%nr1 , * ) !x_group_size )

        CALL SYSTEM_CLOCK( time(3) )
      !------------------------------------------------------
      !------Forward xy-scatter Start------------------------
      
        DO i = dfft%thread_y_start( mythread+1, remswitch ), dfft%thread_y_end( mythread+1, remswitch )
           igroup = ( ( (i-1) / ( dfft%my_nr3p * dfft%nr1s ) ) + 1 )
           offset = dfft%my_nr2p * mod( i-1, dfft%nr1s * dfft%my_nr3p )
           DO j = 1, dfft%my_nr2p
              aux2( j + offset , igroup ) = aux( dfft%map_scatter_fw( j + offset ), igroup )
           ENDDO
        ENDDO
      
      !-------Forward xy-scatter End-------------------------
      !------------------------------------------------------
        CALL SYSTEM_CLOCK( time(4) )
      
        dfft%time_adding( 10 ) = dfft%time_adding( 10 ) + ( time(4) - time(3) )

    END SUBROUTINE Second_Part_x_section

END SUBROUTINE fwfft_x_section_Man

SUBROUTINE fwfft_y_section_Man( dfft, aux, comm_mem_send, comm_mem_recv, map_pcfw, nr1s, ns, batch_size, y_group_size, yset, remswitch, mythread )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: y_group_size, yset, batch_size, remswitch, mythread
  COMPLEX(DP), INTENT(INOUT)  :: comm_mem_send( * )
  COMPLEX(DP), INTENT(INOUT)  :: comm_mem_recv( * )
  INTEGER, INTENT(IN) :: nr1s( * ), ns( * )
!  COMPLEX(DP), INTENT(INOUT) :: aux( dfft%nr2 * nr1s(dfft%mype2+1) * dfft%my_nr3p , * ) !y_group_size )
  COMPLEX(DP), INTENT(INOUT) :: aux( : , : )
  INTEGER, INTENT(IN) :: map_pcfw( * )

  INTEGER :: l, m, i, offset, j, k, igroup, jter, offset2, offset3

  INTEGER(INT64) :: time(3)

  Call First_Part_y_section( aux )

  !$OMP barrier
 
  Call Second_Part_y_section( aux )

  CONTAINS

    SUBROUTINE First_Part_y_section( aux2 )
      
      Implicit NONE
      COMPLEX(DP), INTENT(INOUT) :: aux2( dfft%nr2 , * ) !nr1s(dfft%mype2+1) * dfft%my_nr3p *  y_group_size )

        CALL SYSTEM_CLOCK( time(1) )
      !------------------------------------------------------
      !------------y-FFT Start-------------------------------
      
        CALL fft_noOMP_1D( aux2( : , dfft%thread_y_start( mythread+1, remswitch ) : dfft%thread_y_end( mythread+1, remswitch ) ), dfft%thread_y_sticks, dfft%nr2, remswitch, mythread, dfft%nthreads, -2 )
      
      !-------------y-FFT End--------------------------------
      !------------------------------------------------------
        CALL SYSTEM_CLOCK( time(2) )

        dfft%time_adding( 11 ) = dfft%time_adding( 11 ) + ( time(2) - time(1) )

    END SUBROUTINE First_Part_y_section

    SUBROUTINE Second_Part_y_section( aux2 )

      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT) :: aux2  ( dfft%my_nr1p * dfft%my_nr3p * dfft%nr2, * ) !y_group_size

        CALL SYSTEM_CLOCK( time(2) )

      !------------------------------------------------------
      !---------Pre-Com-Copy Start---------------------------
      
        IF( .not. dfft%single_node ) THEN
      
           !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 1 ), ierr )
      
           DO l = 1, dfft%nodes_numb
              IF( dfft%non_blocking .and. l .eq. dfft%my_node+1 ) CYCLE
              DO m = 1, dfft%node_task_size
                 i = (l-1)*dfft%node_task_size + m
                 offset = ( dfft%my_node_rank + (i-1)*dfft%node_task_size ) * dfft%small_chunks + ( (l-1)*(batch_size-1) + (yset-1)*batch_size ) * dfft%big_chunks
                 DO j = dfft%thread_z_start( mythread+1, remswitch, i ), dfft%thread_z_end( mythread+1, remswitch, i )
                    jter = mod( j-1, ns( i ) ) 
                    igroup = ( (j-1) / ns( i ) )
                    offset2 =  igroup * dfft%big_chunks
                    DO k = 1, dfft%my_nr3p
                       comm_mem_send( offset + offset2 + jter*dfft%nr3px + k ) = &
                       aux2( dfft%map_pcfw( (i-1)*dfft%small_chunks + jter*dfft%nr3px + k ), igroup+1 )
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
              offset = ( dfft%my_node_rank + (i-1)*dfft%node_task_size ) * dfft%small_chunks + ( dfft%my_node*(batch_size-1) + (yset-1)*batch_size ) * dfft%big_chunks
              DO j = dfft%thread_z_start( mythread+1, remswitch, i ), dfft%thread_z_end( mythread+1, remswitch, i )
                 jter = mod( j-1, ns( i ) ) 
                 igroup = ( (j-1) / ns( i ) )
                 offset2 =  igroup * dfft%big_chunks
                 DO k = 1, dfft%my_nr3p
                    comm_mem_recv( offset + offset2 + jter*dfft%nr3px + k ) = &
                    aux2( dfft%map_pcfw( (i-1)*dfft%small_chunks + jter*dfft%nr3px + k ), igroup+1 )
                 END DO
              END DO
           END DO
      
        END IF
      
        IF( dfft%single_node ) THEN
      
           DO m = 1, dfft%node_task_size
              i = dfft%my_node*dfft%node_task_size + m
              offset = ( dfft%my_node_rank + (i-1)*dfft%node_task_size ) * dfft%small_chunks + ( dfft%my_node*(batch_size-1) + (yset-1)*batch_size ) * dfft%big_chunks
              DO j = dfft%thread_z_start( mythread+1, remswitch, i ), dfft%thread_z_end( mythread+1, remswitch, i )
                 jter = mod( j-1, ns( i ) ) 
                 igroup = ( (j-1) / ns( i ) )
                 offset2 =  igroup * dfft%big_chunks
                 DO k = 1, dfft%my_nr3p
                    comm_mem_recv( offset + offset2 + jter*dfft%nr3px + k ) = &
                    aux2( dfft%map_pcfw( (i-1)*dfft%small_chunks + jter*dfft%nr3px + k ), igroup+1 )
                 END DO
              END DO
           END DO
      
        END IF
      
        !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )
      
      !----------Pre-Com-Copy End----------------------------
      !------------------------------------------------------
        CALL SYSTEM_CLOCK( time(3) )
      
        dfft%time_adding( 12 ) = dfft%time_adding( 12 ) + ( time(3) - time(2) )

    END SUBROUTINE Second_Part_y_section
      
END SUBROUTINE fwfft_y_section_Man

SUBROUTINE fwfft_z_section_Man( dfft, comm_mem_recv, aux, iset, batch_size, remswitch, mythread, ns )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iset, batch_size, remswitch, mythread
  INTEGER, INTENT(IN) :: ns( * )
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( * )
  COMPLEX(DP), INTENT(INOUT)  :: aux ( dfft%nr3 , * ) !ns(dfft%mype+1)*z_group_size )

  INTEGER :: j, l, k, i
  INTEGER :: offset, kfrom, ierr

  INTEGER(INT64) :: time(3)

  CALL SYSTEM_CLOCK( time(1) )
!------------------------------------------------------
!--------After-Com-Copy Start--------------------------

  !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, dfft%mpi_window( 2 ), ierr )
 
  DO j = 1, dfft%nodes_numb
     DO l = 1, dfft%node_task_size
        offset = ( dfft%my_node_rank*dfft%node_task_size + (l-1) ) * dfft%small_chunks + ( (j-1)*batch_size + (iset-1)*batch_size ) * dfft%big_chunks
        DO k = dfft%thread_z_start( mythread+1, remswitch, dfft%mype+1 ), dfft%thread_z_end( mythread+1, remswitch, dfft%mype+1 )
           kfrom = offset + dfft%nr3px * mod( k-1, ns(dfft%mype+1) ) + ( (k-1) / ns(dfft%mype+1) ) * dfft%big_chunks
           DO i = 1, dfft%nr3p( (j-1)*dfft%node_task_size + l )
              aux( dfft%nr3p_offset( (j-1)*dfft%node_task_size + l ) + i, k ) = comm_mem_recv( kfrom + i )
!              write(6,*) dfft%my_node_rank, mythread, dfft%nr3p_offset( (j-1)*dfft%node_task_size + l ) + i, k,  kfrom + i
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  !CALL mpi_win_unlock_all( dfft%mpi_window( 2 ), ierr )

!---------After-Com-Copy End---------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(2) )
!------------------------------------------------------
!------------z-FFT Start-------------------------------

  CALL fft_noOMP_1D( aux( : , dfft%thread_z_start( mythread+1, remswitch, dfft%mype+1 ) : dfft%thread_z_end( mythread+1, remswitch, dfft%mype+1 ) ), &
                     dfft%thread_z_sticks(:,:,dfft%mype+1), dfft%nr3, remswitch, mythread, dfft%nthreads, -2 )

!-------------z-FFT End--------------------------------
!------------------------------------------------------
  CALL SYSTEM_CLOCK( time(3) )

  dfft%time_adding( 13 ) = dfft%time_adding( 13 ) + ( time(2) - time(1) )
  dfft%time_adding( 14 ) = dfft%time_adding( 14 ) + ( time(3) - time(2) )
  dfft%counter(6) = dfft%counter(6) + 1

END SUBROUTINE fwfft_z_section_Man

!=----------------------------------------------------------------------=
END MODULE fftpw_batchingManual
!=----------------------------------------------------------------------=
