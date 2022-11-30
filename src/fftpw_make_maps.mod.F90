!=----------------------------------------------------------------------=
MODULE fftpw_make_maps
!=----------------------------------------------------------------------=

  USE fft,                            ONLY: wfn_r,&
                                            fft_batchsize,&
                                            fft_numbatches
  USE fftpw_param
  USE rswfmod,                        ONLY: rsactive
  USE fftpw_types,                    ONLY: PW_fft_type_descriptor
  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: Prep_Copy_Maps, Set_Req_Vals, MapVals_CleanUp

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
     dfft%rsactive = rsactive
     dfft%nr3px = MAXVAL ( dfft%nr3p )
     dfft%my_nr1p = count ( ir1 > 0 )
     dfft%nsx = MAXVAL( ns )
     dfft%small_chunks = dfft%nr3px * MAXVAL( ns )
     dfft%big_chunks = dfft%small_chunks * dfft%node_task_size * dfft%node_task_size
   
     dfft%tscale = 1.0d0 / dble( dfft%nr1 * dfft%nr2 * dfft%nr3 )
     dfft%tscale_gamma = 0.5d0 / dble( dfft%nr1 * dfft%nr2 * dfft%nr3 )

     ALLOCATE( dfft%aux( dfft%nnr ) )

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
     IF( ALLOCATED( dfft%bench_aux_rem ) )    DEALLOCATE( dfft%bench_aux_rem )
   
     ALLOCATE( dfft%bench_aux( dfft%nnr, batch_size ) )
     ALLOCATE( dfft%bench_aux_rem( dfft%nnr, rem_size ) )

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
     DO i = 1, ns( iproc3 )
        it = ( iproc3 - 1 ) * dfft%small_chunks + dfft%nr3px * (i-1)
        mc = dfft%ismap( i + dfft%iss( iproc3 ) ) ! this is  m1+(m2-1)*nr1x  of the  current pencil
        m1 = mod (mc-1,dfft%nr1) + 1
        m2 = (mc-1)/dfft%nr1 + 1
        i1 = m2 + ( ir1(m1) - 1 ) * dfft%nr2
        DO k = 1, dfft%my_nr3p
           dfft%map_pcfw( k + it ) = i1
           i1 = i1 + dfft%nr2 * dfft%my_nr1p
        ENDDO
     ENDDO
     !$omp end do
  ENDDO
  !$omp end parallel

  dfft%fw_off = 0
  DO i = 1, dfft%node_task_size
     dfft%fw_off = dfft%fw_off + ns( dfft%my_node * dfft%node_task_size + i )
  ENDDO

END SUBROUTINE Make_fw_yzCOM_Map

SUBROUTINE MapVals_CleanUp( dfft )
  IMPLICIT NONE

  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft

  IF( ALLOCATED( dfft%aux ) )              DEALLOCATE( dfft%aux )
  IF( ALLOCATED( dfft%bench_aux ) )        DEALLOCATE( dfft%bench_aux )
  IF( ALLOCATED( dfft%bench_aux_rem ) )    DEALLOCATE( dfft%bench_aux_rem )
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


END SUBROUTINE MapVals_CleanUp

!=----------------------------------------------------------------------=
END MODULE fftpw_make_maps
!=----------------------------------------------------------------------=
