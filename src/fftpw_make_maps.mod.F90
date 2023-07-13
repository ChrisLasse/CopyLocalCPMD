!=----------------------------------------------------------------------=
MODULE fftpw_make_maps
!=----------------------------------------------------------------------=

  USE fft,                            ONLY: wfn_r,&
                                            fft_batchsize,&
                                            fft_numbatches
  USE fftpw_param
  USE fftpw_types,                    ONLY: PW_fft_type_descriptor
  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: Prep_Copy_Maps
  PUBLIC :: Set_Req_Vals
  PUBLIC :: MapVals_CleanUp
  PUBLIC :: Prep_fft_com
  PUBLIC :: GIMME_GROUP_SIZES
  PUBLIC :: Make_Manual_Maps

CONTAINS

SUBROUTINE Set_Req_Vals( dfft, nbnd, batch_size, rem_size, num_buff, ir1, ns )
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: batch_size, nbnd, num_buff
  INTEGER, INTENT(OUT) :: rem_size
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft

  INTEGER, INTENT(IN)  :: ir1(:), ns(:)

  INTEGER, SAVE :: last_batch  = 0
  INTEGER, SAVE :: last_buffer = 0
  INTEGER :: i, j, k, l, f

  IF( dfft%make_first ) THEN

     dfft%make_first = .false.
     dfft%nr3px = MAXVAL ( dfft%nr3p )
     dfft%my_nr1p = count ( ir1 > 0 )
     dfft%small_chunks = dfft%nr3px * MAXVAL( ns )
     dfft%big_chunks = dfft%small_chunks * dfft%node_task_size * dfft%node_task_size
   
     dfft%tscale = 1.0d0 / dble( dfft%nr1 * dfft%nr2 * dfft%nr3 )
     dfft%tscale_gamma = 0.5d0 / dble( dfft%nr1 * dfft%nr2 * dfft%nr3 )

     ALLOCATE( dfft%aux( dfft%nnr ) )
     ALLOCATE( dfft%aux2( dfft%nr1w(dfft%mype2+1) * dfft%my_nr3p * dfft%nr2 * dfft%max_batch_size ) ) 
     ALLOCATE( dfft%aux_array( dfft%nr1w(dfft%mype2+1) * dfft%my_nr3p * dfft%nr2 , dfft%max_batch_size ) ) 

  END IF

  IF( last_buffer .ne. num_buff .and. dfft%wave ) THEN

     last_buffer = num_buff
  
     IF( ALLOCATED( dfft%first_step ) )       DEALLOCATE( dfft%first_step )
     IF( ALLOCATED( dfft%first_loading ) )    DEALLOCATE( dfft%first_loading )
     IF( ALLOCATED( dfft%buffer_sequence ) )  DEALLOCATE( dfft%buffer_sequence )

     ALLOCATE( dfft%first_step( num_buff ) )
     ALLOCATE( dfft%first_loading( num_buff ) )
     ALLOCATE( dfft%buffer_sequence( num_buff ) )
     DO i = 1, num_buff
        dfft%buffer_sequence( i ) = i
     END DO
     IF( num_buff .gt. 1 ) THEN
        dfft%buffer_sequence( 1 ) = 2
        dfft%buffer_sequence( 2 ) = 1
     END IF

  END IF

  IF( last_batch .ne. batch_size .and. dfft%wave ) THEN

     last_batch = batch_size

     IF( ALLOCATED( dfft%bench_aux ) )        DEALLOCATE( dfft%bench_aux )
     ALLOCATE( dfft%bench_aux( dfft%my_nr3p * dfft%nr2 * dfft%nr1, batch_size ) )

  END IF

  dfft%first_step = .true.
  rem_size = mod( (nbnd+1)/2, batch_size )

END SUBROUTINE Set_Req_Vals

SUBROUTINE Prep_Copy_Maps( dfft, ngms, batch_size, rem_size, ir1, ns )
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ngms, batch_size, rem_size
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft

  INTEGER, INTENT(IN)  :: ir1(:), ns(:)

  INTEGER, SAVE :: last_batch = 0
  LOGICAL, SAVE :: one_done = .false.

  IF( .not. allocated( dfft%zero_prep_start ) ) THEN

     !Prepare_Psi
     ALLOCATE( dfft%zero_prep_start( ns( dfft%mype+1), 4 ) )
     ALLOCATE( dfft%zero_prep_end( ns( dfft%mype+1), 4 ) )
     ALLOCATE( dfft%prep_map( 6, ns( dfft%mype+1 ) ) )
     CALL Make_PrepPsi_Maps( dfft%zero_prep_start, dfft%zero_prep_end, dfft%prep_map, dfft%nr3, ns(dfft%mype+1), ngms, dfft%nl, dfft%nlm )
   
     !INV_After_Com
     ALLOCATE( dfft%zero_acinv_start( dfft%my_nr1p ) ) 
     ALLOCATE( dfft%zero_acinv_end( dfft%my_nr1p ) ) 
   
     !Scatter_xy
     ALLOCATE( dfft%map_scatter_inv( dfft%nnr ) )
     ALLOCATE( dfft%map_scatter_fw ( dfft%nnr ) )
     IF( dfft%what .eq. 1 ) THEN 
        CALL Make_scatter_Map( dfft, dfft%nr1p, dfft%indp )
     ELSE
        CALL Make_scatter_Map( dfft, dfft%nr1w, dfft%indw )
     END IF 

     !FW_Pre_Com
     ALLOCATE( dfft%map_pcfw( dfft%nr3px * dfft%nproc3 * MAXVAL( ns ) ) )
     CALL Make_fw_yzCOM_Map( dfft, ir1, ns )

  END IF

  IF( .not. one_done .and. batch_size .eq. 1 ) THEN

     ALLOCATE( dfft%map_acinv_one( dfft%my_nr3p * dfft%my_nr1p * dfft%nr2 * batch_size ) )
     ALLOCATE( dfft%map_acinv_rem_one( dfft%my_nr3p * dfft%my_nr1p * dfft%nr2 * rem_size ) )
     CALL Make_inv_yzCOM_Maps( dfft, dfft%map_acinv_one, dfft%map_acinv_rem_one, batch_size, rem_size, ir1, ns )
     one_done = .true.

  END IF

  IF( last_batch .ne. batch_size .and. dfft%wave ) THEN

     last_batch = batch_size

     IF( ALLOCATED( dfft%map_acinv ) )        DEALLOCATE( dfft%map_acinv )
     IF( ALLOCATED( dfft%map_acinv_rem ) )    DEALLOCATE( dfft%map_acinv_rem )
     ALLOCATE( dfft%map_acinv( dfft%my_nr3p * dfft%my_nr1p * dfft%nr2 * batch_size ) )
     ALLOCATE( dfft%map_acinv_rem( dfft%my_nr3p * dfft%my_nr1p * dfft%nr2 * rem_size ) )
     CALL Make_inv_yzCOM_Maps( dfft, dfft%map_acinv, dfft%map_acinv_rem, batch_size, rem_size, ir1, ns )

  END IF

END SUBROUTINE Prep_Copy_Maps


SUBROUTINE Make_PrepPsi_Maps( zero_start, zero_end, prep_map, nr3, my_nsw, ngms, nl, nlm )
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nr3, my_nsw, ngms
  INTEGER, INTENT(IN)  :: nl( : )
  INTEGER, OPTIONAL, INTENT(IN)  :: nlm( : )
  INTEGER, INTENT(OUT) :: zero_start( :, : ), zero_end( :, : ), prep_map( :, : )

  LOGICAL :: l_map( nr3 * my_nsw )
  LOGICAL :: l_map_m( nr3 * my_nsw )
  LOGICAL :: l_map_z( nr3 * my_nsw )
  INTEGER :: i, j
  LOGICAL :: first

  !$omp parallel private( i )
  !$omp do
  DO i = 1, nr3 * my_nsw
     l_map(i) = .false.
     l_map_m(i) = .false.
     l_map_z(i) = .false.
  ENDDO
  !$omp end do
  !$omp do
  DO i = 1, ngms
     l_map( nl(i) ) = .true.
     IF( present( nlm ) ) THEN
        l_map_m( nlm(i) ) = .true.
        l_map_z( nl(i) ) = .true.
        l_map_z( nlm(i) ) = .true.
     END IF
  ENDDO
  !$omp end do nowait
  !$omp end parallel

  zero_start = 1
  zero_end   = nr3
  first = .true.

  !$omp parallel do private( i, first, j )
  DO i = 1, my_nsw
     first = .true.
     DO j = 1, nr3

        IF( l_map( (i-1)*nr3 + j ) .eqv. .true. ) THEN
           IF( first .eqv. .false. ) THEN
              zero_end(i,1) = j-1
              first = .true.
           END IF
        ELSE
           IF( first .eqv. .true. ) THEN
              zero_start(i,1) = j
              first = .false.
           END IF
        END IF

     ENDDO 
  ENDDO
  !$omp end parallel do

  IF( present( nlm ) ) THEN
   
     !$omp parallel do private( i, first, j )
     DO i = 1, my_nsw
        first = .true.
        DO j = 1, nr3
   
           IF( l_map_z( (i-1)*nr3 + j ) .eqv. .true. ) THEN
              IF( first .eqv. .false. ) THEN
                 zero_end(i,2) = j-1
                 first = .true.
              END IF
           ELSE
              IF( first .eqv. .true. ) THEN
                 zero_start(i,2) = j
                 first = .false.
              END IF
           END IF
   
        ENDDO 
     ENDDO
     !$omp end parallel do
   
     !$omp parallel do private( i, first, j )
     DO i = 1, my_nsw
        first = .true.
        DO j = 1, nr3
   
           IF( l_map_m( (i-1)*nr3 + j ) .eqv. .true. ) THEN
              IF( first .eqv. .false. ) THEN
                 zero_end(i,3) = j-1
                 first = .true.
              END IF
           ELSE
              IF( first .eqv. .true. ) THEN
                 zero_start(i,3) = j
                 first = .false.
              END IF
           END IF
   
        ENDDO 
     ENDDO
     !$omp end parallel do

  END IF

  DO i = 1, my_nsw
     zero_start( i, 4 ) = MAXVAL( zero_start( i, 1:3 ), 1 )
     zero_end  ( i, 4 ) = MINVAL( zero_end  ( i, 1:3 ), 1 )
  ENDDO

  DO i = 1, my_nsw

     prep_map( 1 , i ) = zero_start( i , 3 )
     prep_map( 2 , i ) = zero_start( i , 1 )
     prep_map( 3 , i ) = zero_start( i , 2 )
     prep_map( 4 , i ) = zero_end  ( i , 2 )
     prep_map( 5 , i ) = zero_end  ( i , 1 )
     prep_map( 6 , i ) = zero_end  ( i , 3 )

  ENDDO

END SUBROUTINE Make_PrepPsi_Maps

SUBROUTINE Make_inv_yzCOM_Maps( dfft, map_acinv, map_acinv_rem, batch_size, rem_size, ir1, ns )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN)  :: batch_size, rem_size
  INTEGER, INTENT(OUT) :: map_acinv(:), map_acinv_rem(:)

  INTEGER, INTENT(IN)  :: ir1(:), ns(:)

  LOGICAL :: l_map( dfft%my_nr3p * dfft%my_nr1p * dfft%nr2 )
  LOGICAL :: first
  INTEGER :: ibatch, j, l, i, k
  INTEGER :: iproc3, offset, it, mc, m1, m2, i1
  INTEGER :: ierr

  !FOR BATCH_SIZE
  
  l_map = .false.
  dfft%map_acinv = 0
  
  !$omp parallel private( ibatch, j, l, iproc3, offset, i, it, mc, m1, m2, i1, k)
  DO ibatch = 1, batch_size
     DO j = 1, dfft%nodes_numb
        DO l = 1, dfft%node_task_size
           iproc3 = (j-1)*dfft%node_task_size + l
           offset = ( dfft%my_node_rank*dfft%node_task_size + (l-1) ) * dfft%small_chunks + ( (j-1)*batch_size + (ibatch-1) ) * dfft%big_chunks
           !$omp do    
           DO i = 1, ns( iproc3 )
              it = offset + dfft%nr3px * (i-1) 
              mc = dfft%ismap( i + dfft%iss(iproc3) ) ! this is  m1+(m2-1)*nr1x  of the  current pencil
              m1 = mod ( mc-1, dfft%nr1 ) + 1
              m2 = (mc-1)/dfft%nr1 + 1
              i1 = m2 + ( ir1(m1) - 1 ) * dfft%nr2 + (ibatch-1)*dfft%my_nr3p*dfft%my_nr1p*dfft%nr2
              DO k = 1, dfft%my_nr3p
                 IF( ibatch .eq. 1 ) l_map( i1 ) = .true.
                 map_acinv( i1 ) = k + it
                 i1 = i1 + dfft%nr2*dfft%my_nr1p
              ENDDO
           ENDDO
           !$omp end do
        ENDDO
     ENDDO
  ENDDO
  !$omp end parallel
  
  
  !FOR REM_SIZE
  
  dfft%map_acinv_rem = 0
  
  !$omp parallel private( ibatch, j, l, iproc3, offset, i, it, mc, m1, m2, i1, k)
  DO ibatch = 1, rem_size
     DO j = 1, dfft%nodes_numb
        DO l = 1, dfft%node_task_size
           iproc3 = (j-1)*dfft%node_task_size + l
           offset = ( dfft%my_node_rank*dfft%node_task_size + (l-1) ) * dfft%small_chunks + ( (j-1)*rem_size + (ibatch-1) ) * dfft%big_chunks
           !$omp do    
           DO i = 1, ns( iproc3 )
              it = offset + dfft%nr3px * (i-1) 
              mc = dfft%ismap( i + dfft%iss(iproc3) ) ! this is  m1+(m2-1)*nr1x  of the  current pencil
              m1 = mod ( mc-1, dfft%nr1 ) + 1
              m2 = (mc-1)/dfft%nr1 + 1
              i1 = m2 + ( ir1(m1) - 1 ) * dfft%nr2 + (ibatch-1)*dfft%my_nr3p*dfft%my_nr1p*dfft%nr2
              DO k = 1, dfft%my_nr3p
                 map_acinv_rem( i1 ) = k + it
                 i1 = i1 + dfft%nr2*dfft%my_nr1p
              ENDDO
           ENDDO
           !$omp end do
        ENDDO
     ENDDO
  ENDDO
  !$omp end parallel
  
  dfft%zero_acinv_start = 0
  dfft%zero_acinv_end = dfft%nr2
  first = .true.
  
  IF( dfft%mype3 .eq. 0 ) THEN 
     DO j = 1, dfft%my_nr1p
        first = .true.
        DO l = 1, dfft%nr2
  
           IF( l_map( (j-1)*dfft%nr2 + l ) .eqv. .true. ) THEN
              IF( first .eqv. .false. ) THEN
                 dfft%zero_acinv_end( j ) = l-1
                 first = .true.
              END IF
           ELSE
              IF( first .eqv. .true. ) THEN
                 dfft%zero_acinv_start( j ) = l
                 first = .false.
              END IF
           END IF
  
        ENDDO
     ENDDO
  END IF
  IF( dfft%my_node_rank .eq. 0 .and. .not. dfft%single_node ) THEN
     CALL MPI_BCAST( dfft%zero_acinv_start, dfft%my_nr1p, MPI_INTEGER, 0, dfft%inter_node_comm, ierr)
     CALL MPI_BCAST( dfft%zero_acinv_end  , dfft%my_nr1p, MPI_INTEGER, 0, dfft%inter_node_comm, ierr)
  END IF
  CALL MPI_BCAST( dfft%zero_acinv_start, dfft%my_nr1p, MPI_INTEGER, 0, dfft%node_comm, ierr)
  CALL MPI_BCAST( dfft%zero_acinv_end  , dfft%my_nr1p, MPI_INTEGER, 0, dfft%node_comm, ierr)
  
END SUBROUTINE Make_inv_yzCOM_Maps

SUBROUTINE Make_scatter_Map( dfft, nr1s, inds )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN)  :: nr1s( dfft%nproc2 )
  INTEGER, INTENT(IN)  :: inds( dfft%nr1 , dfft%nproc2 )

  INTEGER :: ncpx, sendsize, it, m3, i1, m1, icompact, iproc2, i, j, l
  LOGICAL :: first

  dfft%nr1s = nr1s( 1 )
  ncpx = MAXVAL( nr1s ) * dfft%my_nr3p       ! maximum number of Y columns to be disributed
  sendsize = ncpx * MAXVAL( dfft%nr2p )       ! dimension of the scattered chunks (safe value)

  dfft%map_scatter_inv = 0

  !$omp parallel private(it,m3,i1,m1,icompact)
  DO iproc2 = 1, dfft%nproc2
     !$omp do 
     DO i = 0, ncpx-1
        IF(i>=nr1s(iproc2)*dfft%my_nr3p) CYCLE ! control i from 0 to nr1p_(iproc2)*desc%my_nr3p-1
        it = ( iproc2 - 1 ) * sendsize + MAXVAL( dfft%nr2p ) * i
        m3 = i/nr1s(iproc2)+1
        i1 = mod(i,nr1s(iproc2))+1
        m1 = inds(i1,iproc2)
        icompact = m1 + (m3-1)*dfft%nr1*dfft%my_nr2p
        DO j = 1, dfft%my_nr2p
           dfft%map_scatter_inv( icompact ) = j + it
           icompact = icompact + dfft%nr1
        ENDDO
     ENDDO
     !$omp end do nowait
  ENDDO
  !$omp end parallel

  first = .true.
  
  do l = 1, dfft%nr1
  
     IF( dfft%map_scatter_inv( l ) .eqv. .true. ) THEN
        IF( first .eqv. .false. ) THEN
           dfft%zero_scatter_end = l-1
           EXIT
        END IF
     ELSE
        IF( first .eqv. .true. ) THEN
           dfft%zero_scatter_start = l
           first = .false.
        END IF
     END IF
  
  end do

  !$omp parallel do collapse(2) private(it,m3,i1,m1,icompact)
  DO iproc2 = 1, dfft%nproc2
     DO i = 0, ncpx-1
        IF(i>=nr1s(iproc2)*dfft%my_nr3p) CYCLE ! control i from 0 to nr1p_(iproc2)*desc%my_nr3p-1
        it = ( iproc2 - 1 ) * sendsize + MAXVAL( dfft%nr2p ) * i
        m3 = i/nr1s(iproc2)+1
        i1 = mod(i,nr1s(iproc2))+1
        m1 = inds(i1,iproc2)
        icompact = m1 + (m3-1)*dfft%nr1*dfft%my_nr2p
        DO j = 1, dfft%my_nr2p
           dfft%map_scatter_fw( j + it ) = icompact 
           icompact = icompact + dfft%nr1
        ENDDO
     ENDDO
  ENDDO
  !$omp end parallel do

END SUBROUTINE Make_scatter_Map

SUBROUTINE Make_fw_yzCOM_Map( dfft, ir1, ns )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft

  INTEGER, INTENT(IN)  :: ir1(:), ns(:)

  INTEGER :: iproc3, i, k, it, mc, m1, m2, i1

  dfft%map_pcfw = 0

  !$omp parallel private( iproc3, i, it, mc, m1, m2, i1, k )
  DO iproc3 = 1, dfft%nproc3
     !$omp do
!     DO i = 1, MAXVAL( ns ) !( iproc3 )
     DO i = 1, ns( iproc3 )
!     DO i = 1, ns( dfft%my_node_rank+1 )
        it = ( iproc3 - 1 ) * dfft%small_chunks + dfft%nr3px * (i-1)
        mc = dfft%ismap( i + dfft%iss( iproc3 ) ) ! this is  m1+(m2-1)*nr1x  of the  current pencil
        m1 = mod (mc-1,dfft%nr1) + 1
        m2 = (mc-1)/dfft%nr1 + 1
        i1 = m2 + ( ir1(m1) - 1 ) * dfft%nr2
        DO k = 1, dfft%my_nr3p
           dfft%map_pcfw( k + it ) = i1
           i1 = i1 + dfft%nr2 * dfft%my_nr1p
!           write(6,*) dfft%my_node_rank, k+it, i1
        ENDDO
     ENDDO
     !$omp end do
  ENDDO
  !$omp end parallel

END SUBROUTINE Make_fw_yzCOM_Map

SUBROUTINE MapVals_CleanUp( dfft )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft

  dfft%make_first = .true.
  IF( ALLOCATED( dfft%aux ) )              DEALLOCATE( dfft%aux )
  IF( ALLOCATED( dfft%aux2 ) )             DEALLOCATE( dfft%aux2 )
  IF( ALLOCATED( dfft%bench_aux ) )        DEALLOCATE( dfft%bench_aux )
  IF( ALLOCATED( dfft%first_step ) )       DEALLOCATE( dfft%first_step )
  IF( ALLOCATED( dfft%first_loading ) )    DEALLOCATE( dfft%first_loading )
  IF( ALLOCATED( dfft%buffer_sequence ) )  DEALLOCATE( dfft%buffer_sequence )
  IF( ALLOCATED( dfft%zero_prep_start ) )  DEALLOCATE( dfft%zero_prep_start )
  IF( ALLOCATED( dfft%zero_prep_end ) )    DEALLOCATE( dfft%zero_prep_end )
  IF( ALLOCATED( dfft%map_acinv ) )        DEALLOCATE( dfft%map_acinv )
  IF( ALLOCATED( dfft%map_acinv_rem ) )    DEALLOCATE( dfft%map_acinv_rem )
  IF( ALLOCATED( dfft%zero_acinv_start ) ) DEALLOCATE( dfft%zero_acinv_start )
  IF( ALLOCATED( dfft%zero_acinv_end ) )   DEALLOCATE( dfft%zero_acinv_end )
  IF( ALLOCATED( dfft%map_pcfw ) )         DEALLOCATE( dfft%map_pcfw )
  IF( ALLOCATED( dfft%map_scatter_inv ) )  DEALLOCATE( dfft%map_scatter_inv )
  IF( ALLOCATED( dfft%map_scatter_fw ) )   DEALLOCATE( dfft%map_scatter_fw )

END SUBROUTINE MapVals_CleanUp

SUBROUTINE Prep_fft_com( comm_send, comm_recv, sendsize, sendsize_rem, inter_node_comm, nodes_numb, inter_me, buffer_size, send_handle, recv_handle, send_handle_rem, recv_handle_rem )
  IMPLICIT NONE

  INTEGER, INTENT(IN)                            :: sendsize, sendsize_rem, nodes_numb, inter_me, buffer_size
  TYPE(MPI_COMM), INTENT(IN)                     :: inter_node_comm
  COMPLEX(DP), INTENT(IN)                        :: comm_send( : , : )
  COMPLEX(DP), INTENT(INOUT)                     :: comm_recv( : , : )
  TYPE( MPI_REQUEST ), ALLOCATABLE, INTENT(INOUT)               :: send_handle( : , : )
  TYPE( MPI_REQUEST ), ALLOCATABLE, INTENT(INOUT)               :: recv_handle( : , : )
  TYPE( MPI_REQUEST ), ALLOCATABLE, INTENT(INOUT)               :: send_handle_rem( : , : )
  TYPE( MPI_REQUEST ), ALLOCATABLE, INTENT(INOUT)               :: recv_handle_rem( : , : )

  INTEGER, SAVE :: sendsize_save = 0
  LOGICAL, SAVE :: first = .true.
  INTEGER :: ierr, i, f, j, k
  INTEGER, SAVE :: buffer_size_save

  IF( first ) THEN
     first = .false.
     buffer_size_save = buffer_size
  ELSE
     DO i = 1, buffer_size_save
        DO j = 1, nodes_numb-1
!           IF( send_handle( j , i ) .ne. 0 )     CALL MPI_REQUEST_FREE( send_handle( j , i ) )
!           IF( recv_handle( j , i ) .ne. 0 )     CALL MPI_REQUEST_FREE( recv_handle( j , i ) )
!           IF( send_handle_rem( j , i ) .ne. 0 ) CALL MPI_REQUEST_FREE( send_handle_rem( j , i ) )
!           IF( recv_handle_rem( j , i ) .ne. 0 ) CALL MPI_REQUEST_FREE( recv_handle_rem( j , i ) )
           CALL MPI_REQUEST_FREE( send_handle( j , i ) )
           CALL MPI_REQUEST_FREE( recv_handle( j , i ) )
           CALL MPI_REQUEST_FREE( send_handle_rem( j , i ) )
           CALL MPI_REQUEST_FREE( recv_handle_rem( j , i ) )
        ENDDO
     ENDDO
     buffer_size_save = buffer_size
  END IF

  IF( ALLOCATED( send_handle ) )       DEALLOCATE( send_handle )
  IF( ALLOCATED( recv_handle ) )       DEALLOCATE( recv_handle )
  IF( ALLOCATED( send_handle_rem ) )   DEALLOCATE( send_handle_rem )
  IF( ALLOCATED( recv_handle_rem ) )   DEALLOCATE( recv_handle_rem )

  ALLOCATE( send_handle    ( nodes_numb-1, buffer_size ) )
  ALLOCATE( recv_handle    ( nodes_numb-1, buffer_size ) )
  ALLOCATE( send_handle_rem( nodes_numb-1, buffer_size ) )
  ALLOCATE( recv_handle_rem( nodes_numb-1, buffer_size ) )

  DO k = 1, buffer_size !INITIALIZE SENDING AND RECEIVING

     f = 0
     DO i = 1, nodes_numb
     
        IF( i-1 .eq. inter_me ) CYCLE
        f = f + 1
        CALL MPI_SEND_INIT( comm_send( 1 + (i-1)*sendsize, k ), sendsize, MPI_DOUBLE_COMPLEX, i-1, &
                        inter_me, inter_node_comm, send_handle( f , k ), ierr )
        CALL MPI_RECV_INIT( comm_recv( 1 + (i-1)*sendsize, k ), sendsize, MPI_DOUBLE_COMPLEX, i-1, &
                        MPI_ANY_TAG, inter_node_comm, recv_handle( f , k ), ierr )
     
     ENDDO

  ENDDO

  IF( sendsize_rem .ne. 0 ) THEN

     DO k = 1, buffer_size !INITIALIZE SENDING AND RECEIVING
   
        f = 0
        DO i = 1, nodes_numb
        
           IF( i-1 .eq. inter_me ) CYCLE
           f = f + 1
           CALL MPI_SEND_INIT( comm_send( 1 + (i-1)*sendsize_rem, k ), sendsize_rem, MPI_DOUBLE_COMPLEX, i-1, &
                           inter_me, inter_node_comm, send_handle_rem( f , k ), ierr )
           CALL MPI_RECV_INIT( comm_recv( 1 + (i-1)*sendsize_rem, k ), sendsize_rem, MPI_DOUBLE_COMPLEX, i-1, &
                           MPI_ANY_TAG, inter_node_comm, recv_handle_rem( f , k ), ierr )
        
        ENDDO
   
     ENDDO

  END IF

END SUBROUTINE Prep_fft_com

SUBROUTINE GIMME_GROUP_SIZES( dfft, last_time )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  LOGICAL, INTENT(OUT) :: last_time

  INTEGER :: i

  INTEGER :: x_last_size, y_last_size, z_last_size
  CHARACTER(len=12) :: z_format, y_format, x_format, z_remformat, y_remformat, x_remformat

  last_time = .false.

  IF( ( .not. dfft%autotune_finished .and. dfft%tunning ) .or. .not. last_time ) THEN

     IF( dfft%autotune_finished .or. .not. dfft%tunning ) last_time = .true.

     IF( dfft%y_z_optimal ) THEN
        DO i = 1, dfft%max_batch_size
           IF( dfft%optimal_groups( i, dfft%buffer_size_save, 2 ) .le. dfft%batch_size_save ) THEN
              dfft%y_group_size = dfft%optimal_groups( i, dfft%buffer_size_save, 2 )
              EXIT
           ENDIF
        ENDDO
        DO i = 1, dfft%max_batch_size
           IF( dfft%optimal_groups( i, dfft%buffer_size_save, 3 ) .le. dfft%batch_size_save ) THEN
              dfft%z_group_size = dfft%optimal_groups( i, dfft%buffer_size_save, 3 )
              EXIT
           ENDIF
        ENDDO
     ELSE
        dfft%y_group_size = dfft%y_group_autosize
        dfft%z_group_size = dfft%z_group_autosize
     END IF

     IF( mod( dfft%batch_size_save, dfft%y_group_size ) .eq. 0 ) THEN
         dfft%y_loop_size(1) = dfft%batch_size_save / dfft%y_group_size
         dfft%y_groups(:,1) = dfft%y_group_size
     ELSE
         dfft%y_rem1 = .true.
         dfft%y_loop_size(1) = dfft%batch_size_save / dfft%y_group_size + 1
         y_last_size = dfft%batch_size_save - (dfft%y_loop_size(1)-1) * dfft%y_group_size
         dfft%y_groups(:,1) = dfft%y_group_size
         dfft%y_groups(dfft%y_loop_size(1),1) = y_last_size
     ENDIF
     IF( mod( dfft%batch_size_save, dfft%z_group_size ) .eq. 0 ) THEN
         dfft%z_loop_size(1) = dfft%batch_size_save / dfft%z_group_size
         dfft%z_groups(:,1) = dfft%z_group_size
     ELSE
         dfft%z_loop_size(1) = dfft%batch_size_save / dfft%z_group_size + 1
         z_last_size = dfft%batch_size_save - (dfft%z_loop_size(1)-1) * dfft%z_group_size
         dfft%z_groups(:,1) = dfft%z_group_size
         dfft%z_groups(dfft%z_loop_size(1),1) = z_last_size
     ENDIF

     IF( dfft%rem_size .ne. 0 ) THEN

        IF( dfft%y_z_optimal ) THEN
           DO i = 1, dfft%max_batch_size
              IF( dfft%optimal_groups( i, dfft%buffer_size_save, 2 ) .le. dfft%rem_size ) THEN
                 dfft%y_group_remsize = dfft%optimal_groups( i, dfft%buffer_size_save, 2 )
                 EXIT
              ENDIF
           ENDDO
           DO i = 1, dfft%max_batch_size
              IF( dfft%optimal_groups( i, dfft%buffer_size_save, 3 ) .le. dfft%rem_size ) THEN
                 dfft%z_group_remsize = dfft%optimal_groups( i, dfft%buffer_size_save, 3 )
                 EXIT
              ENDIF
           ENDDO
        ELSE
           dfft%y_group_remsize = MIN( dfft%rem_size, dfft%y_group_autosize )
           dfft%z_group_remsize = MIN( dfft%rem_size, dfft%z_group_autosize )
        END IF

        IF( mod( dfft%rem_size, dfft%y_group_remsize ) .eq. 0 ) THEN
            dfft%y_loop_size(2) = dfft%rem_size / dfft%y_group_remsize
            dfft%y_groups(:,2) = dfft%y_group_remsize
        ELSE
            dfft%y_rem2 = .true.
            dfft%y_loop_size(2) = dfft%rem_size / dfft%y_group_remsize + 1
            y_last_size = dfft%rem_size - (dfft%y_loop_size(2)-1) * dfft%y_group_remsize
            dfft%y_groups(:,2) = dfft%y_group_remsize
            dfft%y_groups(dfft%y_loop_size(2),2) = y_last_size
        ENDIF
        IF( mod( dfft%rem_size, dfft%z_group_remsize ) .eq. 0 ) THEN
            dfft%z_loop_size(2) = dfft%rem_size / dfft%z_group_remsize
            dfft%z_groups(:,2) = dfft%z_group_remsize
        ELSE
            dfft%z_loop_size(2) = dfft%rem_size / dfft%z_group_remsize + 1
            z_last_size = dfft%rem_size - (dfft%z_loop_size(2)-1) * dfft%z_group_remsize
            dfft%z_groups(:,2) = dfft%z_group_remsize
            dfft%z_groups(dfft%z_loop_size(2),2) = z_last_size
        ENDIF

     END IF

!========== x_size -> 4 arrays possible ================

     IF( dfft%x_optimal ) THEN
        DO i = 1, dfft%max_batch_size
           IF( dfft%optimal_groups( i, dfft%buffer_size_save, 1 ) .le. dfft%y_groups(1,1) ) THEN
              dfft%x_group_size = dfft%optimal_groups( i, dfft%buffer_size_save, 1 )
              EXIT
           ENDIF
        ENDDO
     ELSE
        dfft%x_group_size = dfft%x_group_autosize
     END IF
     IF( mod( dfft%y_groups(1,1), dfft%x_group_size ) .eq. 0 ) THEN
         dfft%x_loop_size(1) = dfft%y_groups(1,1) / dfft%x_group_size
         dfft%x_groups(:,1) = dfft%x_group_size
     ELSE
         dfft%x_loop_size(1) = dfft%y_groups(1,1) / dfft%x_group_size + 1
         x_last_size = dfft%y_groups(1,1) - (dfft%x_loop_size(1)-1) * dfft%x_group_size
         dfft%x_groups(:,1) = dfft%x_group_size
         dfft%x_groups(dfft%x_loop_size(1),1) = x_last_size
     ENDIF

     IF( dfft%y_rem1 ) THEN

        IF( dfft%x_optimal ) THEN
           DO i = 1, dfft%max_batch_size
              IF( dfft%optimal_groups( i, dfft%buffer_size_save, 1 ) .le. dfft%y_groups(dfft%y_loop_size(1),1) ) THEN
                 dfft%x_group_size = dfft%optimal_groups( i, dfft%buffer_size_save, 1 )
                 EXIT
              ENDIF
           ENDDO
        ELSE
           dfft%x_group_size = MIN( dfft%y_groups(dfft%y_loop_size(1),1), dfft%x_group_autosize )
        END IF
        IF( mod( dfft%y_groups(dfft%y_loop_size(1),1), dfft%x_group_size ) .eq. 0 ) THEN
            dfft%x_loop_size(2) = dfft%y_groups(dfft%y_loop_size(1),1) / dfft%x_group_size
            dfft%x_groups(:,2) = dfft%x_group_size
        ELSE
            dfft%x_loop_size(2) = dfft%y_groups(dfft%y_loop_size(1),1) / dfft%x_group_size + 1
            x_last_size = dfft%y_groups(dfft%y_loop_size(1),1) - (dfft%x_loop_size(2)-1) * dfft%x_group_size
            dfft%x_groups(:,2) = dfft%x_group_size
            dfft%x_groups(dfft%x_loop_size(2),2) = x_last_size
        ENDIF

     END IF

     IF( dfft%rem_size .ne. 0 ) THEN

        IF( dfft%x_optimal ) THEN
           DO i = 1, dfft%max_batch_size
              IF( dfft%optimal_groups( i, dfft%buffer_size_save, 1 ) .le. dfft%y_groups(1,2) ) THEN
                 dfft%x_group_remsize = dfft%optimal_groups( i, dfft%buffer_size_save, 1 )
                 EXIT
              ENDIF
           ENDDO
        ELSE
           dfft%x_group_remsize = MIN( dfft%y_groups(1,2), dfft%x_group_autosize )
        END IF
        IF( mod( dfft%y_groups(1,2), dfft%x_group_remsize ) .eq. 0 ) THEN
            dfft%x_loop_size(3) = dfft%y_groups(1,2) / dfft%x_group_remsize
            dfft%x_groups(:,3) = dfft%x_group_remsize
        ELSE
            dfft%x_loop_size(3) = dfft%y_groups(1,2) / dfft%x_group_remsize + 1
            x_last_size = dfft%y_groups(1,2) - (dfft%x_loop_size(3)-1) * dfft%x_group_remsize
            dfft%x_groups(:,3) = dfft%x_group_remsize
            dfft%x_groups(dfft%x_loop_size(3),3) = x_last_size
        ENDIF

        IF( dfft%y_rem2 ) THEN

           IF( dfft%x_optimal ) THEN
              DO i = 1, dfft%max_batch_size
                 IF( dfft%optimal_groups( i, dfft%buffer_size_save, 1 ) .le. dfft%y_groups(dfft%y_loop_size(2),2) ) THEN
                    dfft%x_group_remsize = dfft%optimal_groups( i, dfft%buffer_size_save, 1 )
                    EXIT
                 ENDIF
              ENDDO
           ELSE
              dfft%x_group_remsize = MIN( dfft%y_groups(dfft%y_loop_size(2),2), dfft%x_group_autosize )
           END IF
           IF( mod( dfft%y_groups(dfft%y_loop_size(2),2), dfft%x_group_remsize ) .eq. 0 ) THEN
               dfft%x_loop_size(4) = dfft%y_groups(dfft%y_loop_size(2),2) / dfft%x_group_remsize
               dfft%x_groups(:,4) = dfft%x_group_remsize
           ELSE
               dfft%x_loop_size(4) = dfft%y_groups(dfft%y_loop_size(2),2) / dfft%x_group_remsize + 1
               x_last_size = dfft%y_groups(dfft%y_loop_size(2),2) - (dfft%x_loop_size(4)-1) * dfft%x_group_remsize
               dfft%x_groups(:,4) = dfft%x_group_remsize
               dfft%x_groups(dfft%x_loop_size(4),4) = x_last_size
           ENDIF

        END IF

     END IF

  END IF

  IF( dfft%mype .eq. 0 .and. ( .false. .or. last_time ) ) THEN
     write(6,*) " "
     write(6,*) "buffer, batch sizes", dfft%buffer_size_save, dfft%batch_size_save
     IF( .not. last_time ) write(6,*) "x,y,z autogroup sizes", dfft%x_group_autosize, dfft%y_group_autosize, dfft%z_group_autosize
     write(6,*) " "
     write(z_format,'(A5,I2,A5)') "(A10,", dfft%z_loop_size(1), "(I3))"
     write(6,z_format) "z_groups:", ( dfft%z_groups(i,1), i = 1, dfft%z_loop_size(1) )
     write(z_remformat,'(A5,I2,A5)') "(A14,", dfft%z_loop_size(2), "(I3))"
     IF( dfft%rem_size .ne. 0 ) write(6,z_remformat) "z_rem_groups:", ( dfft%z_groups(i,2), i = 1, dfft%z_loop_size(2) )
     write(6,*) " "
     write(y_format,'(A5,I2,A5)') "(A10,", dfft%y_loop_size(1), "(I3))"
     write(6,y_format) "y_groups:", ( dfft%y_groups(i,1), i = 1, dfft%y_loop_size(1) )
     write(y_remformat,'(A5,I2,A5)') "(A14,", dfft%y_loop_size(2), "(I3))"
     IF( dfft%rem_size .ne. 0 ) write(6,y_remformat) "y_rem_groups:", ( dfft%y_groups(i,2), i = 1, dfft%y_loop_size(2) )
     write(6,*) " "
     write(x_format,'(A5,I2,A5)') "(A10,", dfft%x_loop_size(1), "(I3))"
     write(6,x_format) "x_groups:", ( dfft%x_groups(i,1), i = 1, dfft%x_loop_size(1) )
     write(x_remformat,'(A5,I2,A5)') "(A15,", dfft%x_loop_size(2), "(I3))"
     IF( dfft%y_rem1 ) write(6,x_remformat) "x_yrem_groups:", ( dfft%x_groups(i,2), i = 1, dfft%x_loop_size(2) )
     write(x_format,'(A5,I2,A5)') "(A14,", dfft%x_loop_size(3), "(I3))"
     IF( dfft%rem_size .ne. 0 ) write(6,x_format) "x_rem_groups:", ( dfft%x_groups(i,3), i = 1, dfft%x_loop_size(3) )
     write(x_remformat,'(A5,I2,A5)') "(A19,", dfft%x_loop_size(4), "(I3))"
     IF( dfft%y_rem2 ) write(6,x_remformat) "x_rem_yrem_groups:", ( dfft%x_groups(i,4), i = 1, dfft%x_loop_size(4) )
     write(6,*) " "
  END IF

  dfft%y_rem1 = .false.
  dfft%y_rem2 = .false.

END SUBROUTINE

SUBROUTINE Make_Manual_Maps( dfft, batch_size, rem_size )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  INTEGER, INTENT(IN) :: batch_size, rem_size

  INTEGER :: i, j

! z things

  IF( ALLOCATED( dfft%thread_z_sticks ) ) DEALLOCATE( dfft%thread_z_sticks )
  IF( ALLOCATED( dfft%thread_z_start ) ) DEALLOCATE( dfft%thread_z_start )
  IF( ALLOCATED( dfft%thread_z_end ) ) DEALLOCATE( dfft%thread_z_end )
  
  ALLOCATE( dfft%thread_z_sticks( dfft%nthreads, 2, dfft%node_task_size * dfft%nodes_numb ) )
  ALLOCATE( dfft%thread_z_start( dfft%nthreads, 2, dfft%node_task_size * dfft%nodes_numb ) )
  ALLOCATE( dfft%thread_z_end( dfft%nthreads, 2, dfft%node_task_size * dfft%nodes_numb ) )

  DO j = 1, dfft%node_task_size * dfft%nodes_numb

     DO i = 1, dfft%nthreads
        dfft%thread_z_sticks( i, 1, j ) = ( dfft%nsw( j ) * batch_size ) / dfft%nthreads
     ENDDO
     DO i = 1, mod( dfft%nsw( j ) * batch_size, dfft%nthreads )
        dfft%thread_z_sticks( i, 1, j ) = dfft%thread_z_sticks( i, 1, j ) + 1
     ENDDO
   
     dfft%thread_z_start( 1, 1, j ) = 1
     DO i = 2, dfft%nthreads
        dfft%thread_z_start( i, 1, j ) = dfft%thread_z_start( i-1, 1, j ) + dfft%thread_z_sticks( i-1, 1, j )
     ENDDO
     DO i = 1, dfft%nthreads
        dfft%thread_z_end( i, 1, j ) = dfft%thread_z_start( i, 1, j ) + dfft%thread_z_sticks( i, 1, j ) - 1
     ENDDO
   
     IF( rem_size .ne. 0 ) THEN
        DO i = 1, dfft%nthreads
           dfft%thread_z_sticks( i, 2, j ) = ( dfft%nsw( j ) * rem_size ) / dfft%nthreads
        ENDDO
        DO i = 1, mod( dfft%nsw( j ) * rem_size, dfft%nthreads )
           dfft%thread_z_sticks( i, 2, j ) = dfft%thread_z_sticks( i, 2, j ) + 1
        ENDDO
      
        dfft%thread_z_start( 1, 2, j ) = 1
        DO i = 2, dfft%nthreads
           dfft%thread_z_start( i, 2, j ) = dfft%thread_z_start( i-1, 2, j ) + dfft%thread_z_sticks( i-1, 2, j )
        ENDDO
        DO i = 1, dfft%nthreads
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

  DO i = 1, dfft%nthreads
     dfft%thread_y_sticks( i, 1 ) = ( dfft%my_nr1p * dfft%my_nr3p * batch_size ) / dfft%nthreads
  ENDDO
  DO i = 1, mod( dfft%my_nr1p * dfft%my_nr3p * batch_size, dfft%nthreads )
     dfft%thread_y_sticks( i, 1 ) = dfft%thread_y_sticks( i, 1 ) + 1
  ENDDO

  dfft%thread_y_start( 1, 1 ) = 1
  DO i = 2, dfft%nthreads
     dfft%thread_y_start( i, 1 ) = dfft%thread_y_start( i-1, 1 ) + dfft%thread_y_sticks( i-1, 1 )
  ENDDO
  DO i = 1, dfft%nthreads
     dfft%thread_y_end( i, 1 ) = dfft%thread_y_start( i, 1 ) + dfft%thread_y_sticks( i, 1 ) - 1
  ENDDO

  IF( rem_size .ne. 0 ) THEN
     DO i = 1, dfft%nthreads
        dfft%thread_y_sticks( i, 2 ) = ( dfft%my_nr1p * dfft%my_nr3p * rem_size ) / dfft%nthreads
     ENDDO
     DO i = 1, mod( dfft%my_nr1p * dfft%my_nr3p * rem_size, dfft%nthreads )
        dfft%thread_y_sticks( i, 2 ) = dfft%thread_y_sticks( i, 2 ) + 1
     ENDDO
   
     dfft%thread_y_start( 1, 2 ) = 1
     DO i = 2, dfft%nthreads
        dfft%thread_y_start( i, 2 ) = dfft%thread_y_start( i-1, 2 ) + dfft%thread_y_sticks( i-1, 2 )
     ENDDO
     DO i = 1, dfft%nthreads
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

  DO i = 1, dfft%nthreads
     dfft%thread_x_sticks( i, 1 ) = ( dfft%my_nr3p * dfft%my_nr2p * batch_size ) / dfft%nthreads
  ENDDO
  DO i = 1, mod( dfft%my_nr3p * dfft%my_nr2p * batch_size, dfft%nthreads )
     dfft%thread_x_sticks( i, 1 ) = dfft%thread_x_sticks( i, 1 ) + 1
  ENDDO

  dfft%thread_x_start( 1, 1 ) = 1
  DO i = 2, dfft%nthreads
     dfft%thread_x_start( i, 1 ) = dfft%thread_x_start( i-1, 1 ) + dfft%thread_x_sticks( i-1, 1 )
  ENDDO
  DO i = 1, dfft%nthreads
     dfft%thread_x_end( i, 1 ) = dfft%thread_x_start( i, 1 ) + dfft%thread_x_sticks( i, 1 ) - 1
  ENDDO

  IF( rem_size .ne. 0 ) THEN
     DO i = 1, dfft%nthreads
        dfft%thread_x_sticks( i, 2 ) = ( dfft%my_nr3p * dfft%my_nr2p * rem_size ) / dfft%nthreads
     ENDDO
     DO i = 1, mod( dfft%my_nr3p * dfft%my_nr2p * rem_size, dfft%nthreads )
        dfft%thread_x_sticks( i, 2 ) = dfft%thread_x_sticks( i, 2 ) + 1
     ENDDO
   
     dfft%thread_x_start( 1, 2 ) = 1
     DO i = 2, dfft%nthreads
        dfft%thread_x_start( i, 2 ) = dfft%thread_x_start( i-1, 2 ) + dfft%thread_x_sticks( i-1, 2 )
     ENDDO
     DO i = 1, dfft%nthreads
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

  DO i = 1, dfft%nthreads
     dfft%thread_ngms( i ) = ( dfft%ngms ) / dfft%nthreads
  ENDDO
  DO i = 1, mod( dfft%ngms, dfft%nthreads )
     dfft%thread_ngms( i ) = dfft%thread_ngms( i ) + 1
  ENDDO

  dfft%thread_ngms_start( 1 ) = 1
  DO i = 2, dfft%nthreads
     dfft%thread_ngms_start( i ) = dfft%thread_ngms_start( i-1 ) + dfft%thread_ngms( i-1 )
  ENDDO
  DO i = 1, dfft%nthreads
     dfft%thread_ngms_end( i ) = dfft%thread_ngms_start( i ) + dfft%thread_ngms( i ) - 1
  ENDDO

! rspace things

  IF( ALLOCATED( dfft%thread_rspace ) ) DEALLOCATE( dfft%thread_rspace )
  IF( ALLOCATED( dfft%thread_rspace_start ) ) DEALLOCATE( dfft%thread_rspace_start )
  IF( ALLOCATED( dfft%thread_rspace_end ) ) DEALLOCATE( dfft%thread_rspace_end )
  
  ALLOCATE( dfft%thread_rspace( dfft%nthreads ) )
  ALLOCATE( dfft%thread_rspace_start( dfft%nthreads ) )
  ALLOCATE( dfft%thread_rspace_end( dfft%nthreads ) )

  DO i = 1, dfft%nthreads
     dfft%thread_rspace( i ) = ( dfft%my_nr3p * dfft%nr2 * dfft%nr1 ) / dfft%nthreads
  ENDDO
  DO i = 1, mod( dfft%my_nr3p * dfft%nr2 * dfft%nr1, dfft%nthreads )
     dfft%thread_rspace( i ) = dfft%thread_rspace( i ) + 1
  ENDDO

  dfft%thread_rspace_start( 1 ) = 1
  DO i = 2, dfft%nthreads
     dfft%thread_rspace_start( i ) = dfft%thread_rspace_start( i-1 ) + dfft%thread_rspace( i-1 )
  ENDDO
  DO i = 1, dfft%nthreads
     dfft%thread_rspace_end( i ) = dfft%thread_rspace_start( i ) + dfft%thread_rspace( i ) - 1
  ENDDO
  

END SUBROUTINE Make_Manual_Maps

!=----------------------------------------------------------------------=
END MODULE fftpw_make_maps
!=----------------------------------------------------------------------=
