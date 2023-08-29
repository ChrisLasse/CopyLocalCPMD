!=----------------------------------------------------------------------=
MODULE fftpw_make_maps
!=----------------------------------------------------------------------=

  USE fft,                            ONLY: wfn_r,&
                                            fft_batchsize,&
                                            fft_numbatches
  USE fftpw_param
  USE fftpw_types,                    ONLY: PW_fft_type_descriptor
  USE timer,                          ONLY: tihalt,&
                                            tiset
  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: Prep_single_fft_com
  PUBLIC :: Prep_fft_com
  PUBLIC :: Make_Manual_Maps

CONTAINS

SUBROUTINE Prep_single_fft_com( comm_send, comm_recv, sendsize, comm, nodes_numb, mype, my_node, my_node_rank, node_task_size, send_handle, recv_handle, comm_sendrecv, do_comm )
  IMPLICIT NONE

  INTEGER, INTENT(IN)                                 :: sendsize, nodes_numb, mype, my_node, my_node_rank, node_task_size
  TYPE(MPI_COMM), INTENT(IN)                          :: comm
  COMPLEX(DP), INTENT(IN)                             :: comm_send( : )
  COMPLEX(DP), INTENT(INOUT)                          :: comm_recv( : )
  TYPE( MPI_REQUEST ), ALLOCATABLE, INTENT(INOUT)     :: send_handle( : , : )
  TYPE( MPI_REQUEST ), ALLOCATABLE, INTENT(INOUT)     :: recv_handle( : , : )
  INTEGER, ALLOCATABLE, INTENT(OUT)                   :: comm_sendrecv( : )
  LOGICAL, INTENT(OUT)                                :: do_comm

  INTEGER, ALLOCATABLE :: comm_info_send(:,:), comm_info_recv(:,:)

  INTEGER :: i, j, k, l, m, n, p
  INTEGER :: ierr, eff_nodes, rem, origin_node, target_node
  INTEGER :: howmany_sending( node_task_size, nodes_numb ) , howmany_receiving( node_task_size, nodes_numb )
  LOGICAL :: com_done( nodes_numb )

  eff_nodes = nodes_numb - 1

  howmany_sending   = eff_nodes / node_task_size
  howmany_receiving = eff_nodes / node_task_size
  rem = mod( eff_nodes, node_task_size )
  DO i = 1, rem * 2
     k = mod( i-1, node_task_size ) + 1
     IF( i .le. rem ) howmany_sending( k , : )   = howmany_sending( k , : )   + 1
     IF( i .gt. rem ) howmany_receiving( k , : ) = howmany_receiving( k , : ) + 1
  ENDDO

  ALLOCATE( comm_info_send( howmany_sending( 1 , 1 ), node_task_size*nodes_numb ) )
  ALLOCATE( comm_info_recv( howmany_sending( 1 , 1 ), node_task_size*nodes_numb ) )
  ALLOCATE( comm_sendrecv( 2 ) )
  comm_info_send = 0
  comm_info_recv = 0
  comm_sendrecv(1) = howmany_sending  ( my_node_rank+1 , my_node+1 )
  comm_sendrecv(2) = howmany_receiving( my_node_rank+1 , my_node+1 )

  do_comm = .true.
  IF( comm_sendrecv(1) .eq. 0 .and. comm_sendrecv(2) .eq. 0 ) do_comm = .false.

  DO i = 1, nodes_numb

     com_done = .false.

     DO k = 1, node_task_size 

s_loop: DO l = 1, howmany_sending( k, i )

           DO m = 1, nodes_numb
              IF( m .eq. i ) CYCLE

              DO n = 1, node_task_size 

                 IF( .not. com_done( m ) .and. howmany_receiving( n , m ) .ne. 0 ) THEN
                    com_done( m ) = .true.
                    howmany_receiving( n , m ) = howmany_receiving( n , m ) - 1
                    comm_info_send( l , k + (i-1)*node_task_size ) = n + (m-1)*node_task_size
                    DO p = 1, howmany_sending( 1, 1 )
                       IF( comm_info_recv( p , n + (m-1)*node_task_size ) .eq. 0 ) THEN
                          comm_info_recv( p , n + (m-1)*node_task_size ) = k + (i-1)*node_task_size
                          EXIT
                       END IF
                    ENDDO
                    CYCLE s_loop
                 END IF 

              ENDDO

           ENDDO

        ENDDO s_loop

     ENDDO

  ENDDO
  
  IF( ALLOCATED( send_handle ) )       DEALLOCATE( send_handle )
  IF( ALLOCATED( recv_handle ) )       DEALLOCATE( recv_handle )
  
  ALLOCATE( send_handle    ( comm_sendrecv(1), 1 ) )
  ALLOCATE( recv_handle    ( comm_sendrecv(2), 1 ) )
  
  DO j = 1, comm_sendrecv( 1 ) !INITIALIZE SENDING AND RECEIVING

     target_node = (comm_info_send(j,mype+1)-1) / node_task_size

     CALL MPI_SEND_INIT( comm_send( 1 + target_node*sendsize ), sendsize, MPI_DOUBLE_COMPLEX, comm_info_send( j , mype+1 ) - 1, &
                         mype, comm, send_handle( j , 1 ), ierr )

  ENDDO        

  DO j = 1, comm_sendrecv( 2 )

     origin_node = (comm_info_recv(j,mype+1)-1) / node_task_size

     CALL MPI_RECV_INIT( comm_recv( 1 + origin_node*sendsize ), sendsize, MPI_DOUBLE_COMPLEX, comm_info_recv( j , mype+1 ) - 1, &
                         MPI_ANY_TAG, comm, recv_handle( j , 1 ), ierr )

  ENDDO        

END SUBROUTINE Prep_single_fft_com

SUBROUTINE Prep_fft_com( comm_send, comm_recv, sendsize, sendsize_rem, comm, nodes_numb, mype, my_node, my_node_rank, node_task_size, buffer_size, send_handle, recv_handle, send_handle_rem, recv_handle_rem, comm_sendrecv, do_comm )
  IMPLICIT NONE

  INTEGER, INTENT(IN)                                 :: sendsize, sendsize_rem, nodes_numb, mype, my_node, my_node_rank, node_task_size, buffer_size
  TYPE(MPI_COMM), INTENT(IN)                          :: comm
  COMPLEX(DP), INTENT(IN)                             :: comm_send( : , : )
  COMPLEX(DP), INTENT(INOUT)                          :: comm_recv( : , : )
  TYPE( MPI_REQUEST ), ALLOCATABLE, INTENT(INOUT)     :: send_handle( : , : )
  TYPE( MPI_REQUEST ), ALLOCATABLE, INTENT(INOUT)     :: recv_handle( : , : )
  TYPE( MPI_REQUEST ), ALLOCATABLE, INTENT(INOUT)     :: send_handle_rem( : , : )
  TYPE( MPI_REQUEST ), ALLOCATABLE, INTENT(INOUT)     :: recv_handle_rem( : , : )
  INTEGER, ALLOCATABLE, INTENT(OUT)                   :: comm_sendrecv( : )
  LOGICAL, INTENT(OUT)                                :: do_comm

  LOGICAL, SAVE :: first = .true.
  INTEGER, SAVE :: buffer_size_save

  INTEGER, ALLOCATABLE :: comm_info_send(:,:), comm_info_recv(:,:)

  INTEGER :: i, j, k, l, m, n, p
  INTEGER :: ierr, eff_nodes, rem, origin_node, target_node
  INTEGER :: howmany_sending( node_task_size, nodes_numb ) , howmany_receiving( node_task_size, nodes_numb )
  LOGICAL :: com_done( nodes_numb )

  eff_nodes = nodes_numb - 1

  howmany_sending   = eff_nodes / node_task_size
  howmany_receiving = eff_nodes / node_task_size
  rem = mod( eff_nodes, node_task_size )
  DO i = 1, rem * 2
     k = mod( i-1, node_task_size ) + 1
     IF( i .le. rem ) howmany_sending( k , : )   = howmany_sending( k , : )   + 1
     IF( i .gt. rem ) howmany_receiving( k , : ) = howmany_receiving( k , : ) + 1
  ENDDO

  ALLOCATE( comm_info_send( howmany_sending( 1 , 1 ), node_task_size*nodes_numb ) )
  ALLOCATE( comm_info_recv( howmany_sending( 1 , 1 ), node_task_size*nodes_numb ) )
  ALLOCATE( comm_sendrecv( 2 ) )
  comm_info_send = 0
  comm_info_recv = 0
  comm_sendrecv(1) = howmany_sending  ( my_node_rank+1 , my_node+1 )
  comm_sendrecv(2) = howmany_receiving( my_node_rank+1 , my_node+1 )

  do_comm = .true.
  IF( comm_sendrecv(1) .eq. 0 .and. comm_sendrecv(2) .eq. 0 ) do_comm = .false.

  DO i = 1, nodes_numb

     com_done = .false.

     DO k = 1, node_task_size 

s_loop: DO l = 1, howmany_sending( k, i )

           DO m = 1, nodes_numb
              IF( m .eq. i ) CYCLE

              DO n = 1, node_task_size 

                 IF( .not. com_done( m ) .and. howmany_receiving( n , m ) .ne. 0 ) THEN
                    com_done( m ) = .true.
                    howmany_receiving( n , m ) = howmany_receiving( n , m ) - 1
                    comm_info_send( l , k + (i-1)*node_task_size ) = n + (m-1)*node_task_size
                    DO p = 1, howmany_sending( 1, 1 )
                       IF( comm_info_recv( p , n + (m-1)*node_task_size ) .eq. 0 ) THEN
                          comm_info_recv( p , n + (m-1)*node_task_size ) = k + (i-1)*node_task_size
                          EXIT
                       END IF
                    ENDDO
                    CYCLE s_loop
                 END IF 

              ENDDO

           ENDDO

        ENDDO s_loop

     ENDDO

  ENDDO

  IF( first ) THEN
     first = .false.
     buffer_size_save = buffer_size
  ELSE
     DO i = 1, buffer_size_save
        DO j = 1, comm_sendrecv(1)
           CALL MPI_REQUEST_FREE( send_handle( j , i ) )
           CALL MPI_REQUEST_FREE( send_handle_rem( j , i ) )
        ENDDO
        DO j = 1, comm_sendrecv(2)
           CALL MPI_REQUEST_FREE( recv_handle( j , i ) )
           CALL MPI_REQUEST_FREE( recv_handle_rem( j , i ) )
        ENDDO
     ENDDO
     buffer_size_save = buffer_size
  END IF
  
  IF( ALLOCATED( send_handle ) )       DEALLOCATE( send_handle )
  IF( ALLOCATED( send_handle_rem ) )   DEALLOCATE( send_handle_rem )
  IF( ALLOCATED( recv_handle ) )       DEALLOCATE( recv_handle )
  IF( ALLOCATED( recv_handle_rem ) )   DEALLOCATE( recv_handle_rem )
  
  ALLOCATE( send_handle    ( comm_sendrecv(1), buffer_size ) )
  ALLOCATE( send_handle_rem( comm_sendrecv(1), buffer_size ) )
  ALLOCATE( recv_handle    ( comm_sendrecv(2), buffer_size ) )
  ALLOCATE( recv_handle_rem( comm_sendrecv(2), buffer_size ) )
  
  DO i = 1, buffer_size !INITIALIZE SENDING AND RECEIVING

     DO j = 1, comm_sendrecv( 1 )

        target_node = (comm_info_send(j,mype+1)-1) / node_task_size

        CALL MPI_SEND_INIT( comm_send( 1 + target_node*sendsize, i ), sendsize, MPI_DOUBLE_COMPLEX, comm_info_send( j , mype+1 ) - 1, &
                            mype, comm, send_handle( j , i ), ierr )

     ENDDO        

     DO j = 1, comm_sendrecv( 2 )

        origin_node = (comm_info_recv(j,mype+1)-1) / node_task_size

        CALL MPI_RECV_INIT( comm_recv( 1 + origin_node*sendsize, i ), sendsize, MPI_DOUBLE_COMPLEX, comm_info_recv( j , mype+1 ) - 1, &
                            MPI_ANY_TAG, comm, recv_handle( j , i ), ierr )

     ENDDO        

  ENDDO

  IF( sendsize_rem .ne. 0 ) THEN

     DO i = 1, buffer_size !INITIALIZE SENDING AND RECEIVING
  
        DO j = 1, comm_sendrecv( 1 )

           target_node = (comm_info_send(j,mype+1)-1) / node_task_size
  
           CALL MPI_SEND_INIT( comm_send( 1 + target_node*sendsize_rem, i ), sendsize_rem, MPI_DOUBLE_COMPLEX, comm_info_send( j , mype+1 ) - 1, &
                               mype, comm, send_handle_rem( j , i ), ierr )
  
        ENDDO        
  
        DO j = 1, comm_sendrecv( 2 )
  
           origin_node = (comm_info_recv(j,mype+1)-1) / node_task_size
  
           CALL MPI_RECV_INIT( comm_recv( 1 + origin_node*sendsize_rem, i ), sendsize_rem, MPI_DOUBLE_COMPLEX, comm_info_recv( j , mype+1 ) - 1, &
                               MPI_ANY_TAG, comm, recv_handle_rem( j , i ), ierr )
  
        ENDDO        
  
     ENDDO

  END IF

END SUBROUTINE Prep_fft_com

SUBROUTINE Make_Manual_Maps( dfft, batch_size, rem_size )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: batch_size, rem_size

  INTEGER :: i, j, overlap_cor

  overlap_cor = dfft%nthreads - dfft%eff_nthreads

! z things

  IF( ALLOCATED( dfft%thread_z_sticks ) ) DEALLOCATE( dfft%thread_z_sticks )
  IF( ALLOCATED( dfft%thread_z_start ) ) DEALLOCATE( dfft%thread_z_start )
  IF( ALLOCATED( dfft%thread_z_end ) ) DEALLOCATE( dfft%thread_z_end )
  
  ALLOCATE( dfft%thread_z_sticks( dfft%nthreads, 2, dfft%node_task_size * dfft%nodes_numb ) )
  ALLOCATE( dfft%thread_z_start( dfft%nthreads, 2, dfft%node_task_size * dfft%nodes_numb ) )
  ALLOCATE( dfft%thread_z_end( dfft%nthreads, 2, dfft%node_task_size * dfft%nodes_numb ) )

  DO j = 1, dfft%node_task_size * dfft%nodes_numb

     DO i = 1+overlap_cor, dfft%nthreads
        dfft%thread_z_sticks( i, 1, j ) = ( dfft%nsw( j ) * batch_size ) / dfft%eff_nthreads
     ENDDO
     DO i = 1+overlap_cor, mod( dfft%nsw( j ) * batch_size, dfft%eff_nthreads ) + overlap_cor
        dfft%thread_z_sticks( i, 1, j ) = dfft%thread_z_sticks( i, 1, j ) + 1
     ENDDO
   
     dfft%thread_z_start( 1+overlap_cor, 1, j ) = 1
     DO i = 2+overlap_cor, dfft%nthreads
        dfft%thread_z_start( i, 1, j ) = dfft%thread_z_start( i-1, 1, j ) + dfft%thread_z_sticks( i-1, 1, j )
     ENDDO
     DO i = 1+overlap_cor, dfft%nthreads
        dfft%thread_z_end( i, 1, j ) = dfft%thread_z_start( i, 1, j ) + dfft%thread_z_sticks( i, 1, j ) - 1
     ENDDO
   
     IF( rem_size .ne. 0 ) THEN
        DO i = 1+overlap_cor, dfft%nthreads
           dfft%thread_z_sticks( i, 2, j ) = ( dfft%nsw( j ) * rem_size ) / dfft%eff_nthreads
        ENDDO
        DO i = 1+overlap_cor, mod( dfft%nsw( j ) * rem_size, dfft%eff_nthreads ) + overlap_cor
           dfft%thread_z_sticks( i, 2, j ) = dfft%thread_z_sticks( i, 2, j ) + 1
        ENDDO
      
        dfft%thread_z_start( 1+overlap_cor, 2, j ) = 1
        DO i = 2+overlap_cor, dfft%nthreads
           dfft%thread_z_start( i, 2, j ) = dfft%thread_z_start( i-1, 2, j ) + dfft%thread_z_sticks( i-1, 2, j )
        ENDDO
        DO i = 1+overlap_cor, dfft%nthreads
           dfft%thread_z_end( i, 2, j ) = dfft%thread_z_start( i, 2, j ) + dfft%thread_z_sticks( i, 2, j ) - 1
        ENDDO
     END IF

  ENDDO

! y things

  IF( ALLOCATED( dfft%thread_y_sticks ) ) DEALLOCATE( dfft%thread_y_sticks )
  IF( ALLOCATED( dfft%thread_y_start ) ) DEALLOCATE( dfft%thread_y_start )
  IF( ALLOCATED( dfft%thread_y_end ) ) DEALLOCATE( dfft%thread_y_end )
  
  ALLOCATE( dfft%thread_y_sticks( dfft%nthreads, 2 ) )
  ALLOCATE( dfft%thread_y_start( dfft%nthreads, 2 ) )
  ALLOCATE( dfft%thread_y_end( dfft%nthreads, 2 ) )

  DO i = 1+overlap_cor, dfft%nthreads
     dfft%thread_y_sticks( i, 1 ) = ( dfft%my_nr1p * dfft%my_nr3p * batch_size ) / dfft%eff_nthreads
  ENDDO
  DO i = 1+overlap_cor, mod( dfft%my_nr1p * dfft%my_nr3p * batch_size, dfft%eff_nthreads ) + overlap_cor
     dfft%thread_y_sticks( i, 1 ) = dfft%thread_y_sticks( i, 1 ) + 1
  ENDDO

  dfft%thread_y_start( 1+overlap_cor, 1 ) = 1
  DO i = 2+overlap_cor, dfft%nthreads
     dfft%thread_y_start( i, 1 ) = dfft%thread_y_start( i-1, 1 ) + dfft%thread_y_sticks( i-1, 1 )
  ENDDO
  DO i = 1+overlap_cor, dfft%nthreads
     dfft%thread_y_end( i, 1 ) = dfft%thread_y_start( i, 1 ) + dfft%thread_y_sticks( i, 1 ) - 1
  ENDDO

  IF( rem_size .ne. 0 ) THEN
     DO i = 1+overlap_cor, dfft%nthreads
        dfft%thread_y_sticks( i, 2 ) = ( dfft%my_nr1p * dfft%my_nr3p * rem_size ) / dfft%eff_nthreads
     ENDDO
     DO i = 1+overlap_cor, mod( dfft%my_nr1p * dfft%my_nr3p * rem_size, dfft%eff_nthreads ) + overlap_cor
        dfft%thread_y_sticks( i, 2 ) = dfft%thread_y_sticks( i, 2 ) + 1
     ENDDO
   
     dfft%thread_y_start( 1+overlap_cor, 2 ) = 1
     DO i = 2+overlap_cor, dfft%nthreads
        dfft%thread_y_start( i, 2 ) = dfft%thread_y_start( i-1, 2 ) + dfft%thread_y_sticks( i-1, 2 )
     ENDDO
     DO i = 1+overlap_cor, dfft%nthreads
        dfft%thread_y_end( i, 2 ) = dfft%thread_y_start( i, 2 ) + dfft%thread_y_sticks( i, 2 ) - 1
     ENDDO
  END IF

! x things

  IF( ALLOCATED( dfft%thread_x_sticks ) ) DEALLOCATE( dfft%thread_x_sticks )
  IF( ALLOCATED( dfft%thread_x_start ) ) DEALLOCATE( dfft%thread_x_start )
  IF( ALLOCATED( dfft%thread_x_end ) ) DEALLOCATE( dfft%thread_x_end )
  
  ALLOCATE( dfft%thread_x_sticks( dfft%nthreads, 2 ) )
  ALLOCATE( dfft%thread_x_start( dfft%nthreads, 2 ) )
  ALLOCATE( dfft%thread_x_end( dfft%nthreads, 2 ) )

  DO i = 1+overlap_cor, dfft%nthreads
     dfft%thread_x_sticks( i, 1 ) = ( dfft%my_nr3p * dfft%my_nr2p * batch_size ) / dfft%eff_nthreads
  ENDDO
  DO i = 1+overlap_cor, mod( dfft%my_nr3p * dfft%my_nr2p * batch_size, dfft%eff_nthreads ) + overlap_cor
     dfft%thread_x_sticks( i, 1 ) = dfft%thread_x_sticks( i, 1 ) + 1
  ENDDO

  dfft%thread_x_start( 1+overlap_cor, 1 ) = 1
  DO i = 2+overlap_cor, dfft%nthreads
     dfft%thread_x_start( i, 1 ) = dfft%thread_x_start( i-1, 1 ) + dfft%thread_x_sticks( i-1, 1 )
  ENDDO
  DO i = 1+overlap_cor, dfft%nthreads
     dfft%thread_x_end( i, 1 ) = dfft%thread_x_start( i, 1 ) + dfft%thread_x_sticks( i, 1 ) - 1
  ENDDO

  IF( rem_size .ne. 0 ) THEN
     DO i = 1+overlap_cor, dfft%nthreads
        dfft%thread_x_sticks( i, 2 ) = ( dfft%my_nr3p * dfft%my_nr2p * rem_size ) / dfft%eff_nthreads
     ENDDO
     DO i = 1+overlap_cor, mod( dfft%my_nr3p * dfft%my_nr2p * rem_size, dfft%eff_nthreads ) + overlap_cor
        dfft%thread_x_sticks( i, 2 ) = dfft%thread_x_sticks( i, 2 ) + 1
     ENDDO
   
     dfft%thread_x_start( 1+overlap_cor, 2 ) = 1
     DO i = 2+overlap_cor, dfft%nthreads
        dfft%thread_x_start( i, 2 ) = dfft%thread_x_start( i-1, 2 ) + dfft%thread_x_sticks( i-1, 2 )
     ENDDO
     DO i = 1+overlap_cor, dfft%nthreads
        dfft%thread_x_end( i, 2 ) = dfft%thread_x_start( i, 2 ) + dfft%thread_x_sticks( i, 2 ) - 1
     ENDDO
  END IF

! gspace things

  IF( ALLOCATED( dfft%thread_ngms ) ) DEALLOCATE( dfft%thread_ngms )
  IF( ALLOCATED( dfft%thread_ngms_start ) ) DEALLOCATE( dfft%thread_ngms_start )
  IF( ALLOCATED( dfft%thread_ngms_end ) ) DEALLOCATE( dfft%thread_ngms_end )
  
  ALLOCATE( dfft%thread_ngms( dfft%nthreads ) )
  ALLOCATE( dfft%thread_ngms_start( dfft%nthreads ) )
  ALLOCATE( dfft%thread_ngms_end( dfft%nthreads ) )

  DO i = 1+overlap_cor, dfft%nthreads
     dfft%thread_ngms( i ) = ( dfft%ngms ) / dfft%eff_nthreads
  ENDDO
  DO i = 1+overlap_cor, mod( dfft%ngms, dfft%eff_nthreads ) + overlap_cor
     dfft%thread_ngms( i ) = dfft%thread_ngms( i ) + 1
  ENDDO

  dfft%thread_ngms_start( 1+overlap_cor ) = 1
  DO i = 2+overlap_cor, dfft%nthreads
     dfft%thread_ngms_start( i ) = dfft%thread_ngms_start( i-1 ) + dfft%thread_ngms( i-1 )
  ENDDO
  DO i = 1+overlap_cor, dfft%nthreads
     dfft%thread_ngms_end( i ) = dfft%thread_ngms_start( i ) + dfft%thread_ngms( i ) - 1
  ENDDO

! rspace things

  IF( ALLOCATED( dfft%thread_rspace ) ) DEALLOCATE( dfft%thread_rspace )
  IF( ALLOCATED( dfft%thread_rspace_start ) ) DEALLOCATE( dfft%thread_rspace_start )
  IF( ALLOCATED( dfft%thread_rspace_end ) ) DEALLOCATE( dfft%thread_rspace_end )
  
  ALLOCATE( dfft%thread_rspace( dfft%nthreads ) )
  ALLOCATE( dfft%thread_rspace_start( dfft%nthreads ) )
  ALLOCATE( dfft%thread_rspace_end( dfft%nthreads ) )

  DO i = 1+overlap_cor, dfft%nthreads
     dfft%thread_rspace( i ) = ( dfft%my_nr3p * dfft%nr2 * dfft%nr1 ) / dfft%eff_nthreads
  ENDDO
  DO i = 1+overlap_cor, mod( dfft%my_nr3p * dfft%nr2 * dfft%nr1, dfft%eff_nthreads ) + overlap_cor
     dfft%thread_rspace( i ) = dfft%thread_rspace( i ) + 1
  ENDDO

  dfft%thread_rspace_start( 1+overlap_cor ) = 1
  DO i = 2+overlap_cor, dfft%nthreads
     dfft%thread_rspace_start( i ) = dfft%thread_rspace_start( i-1 ) + dfft%thread_rspace( i-1 )
  ENDDO
  DO i = 1+overlap_cor, dfft%nthreads
     dfft%thread_rspace_end( i ) = dfft%thread_rspace_start( i ) + dfft%thread_rspace( i ) - 1
  ENDDO
  

END SUBROUTINE Make_Manual_Maps

!=----------------------------------------------------------------------=
END MODULE fftpw_make_maps
!=----------------------------------------------------------------------=
