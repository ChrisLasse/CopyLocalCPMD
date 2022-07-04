!#INCLUDE "/apps/intel/ComposerXE2019/compilers_and_libraries_2019.5.281/linux/mkl/include/mkl_service.f90"

!=----------------------------------------------------------------------=
MODULE fftpw_overlapp
!=----------------------------------------------------------------------=


  USE fftpw_batching
  USE fftpw_make_maps,                          ONLY: Set_Req_Vals,&
                                                      Prep_Copy_Maps,&
                                                      MapVals_CleanUp
  USE fftpw_param
  USE fftpw_single_calls,                       ONLY: Prepare_Psi_single,&
                                                      invfft_single
  USE fftpw_types,                              ONLY: PW_fft_type_descriptor,&
                                                      create_shared_memory_window_2d 
  USE mkl_service
  USE iso_fortran_env

  !$ USE omp_lib, ONLY: omp_get_thread_num, &
  !$                    omp_get_num_threads,&
  !$                    omp_set_num_threads,&
  !$                    omp_in_parallel

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Overlapp_FFT
  PUBLIC :: Compare_invffts

CONTAINS

SUBROUTINE Overlapp_FFT( dfft, psi, hpsi, v, ngms, nbnd )

  IMPLICIT NONE
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(IN)  :: psi (ngms, nbnd)
  COMPLEX(DP), INTENT(INOUT) :: hpsi(ngms, nbnd)
  INTEGER, INTENT(IN)    :: ngms, nbnd
  REAL(DP), INTENT(IN) :: v(dfft%nnr)

  COMPLEX(DP), POINTER, SAVE, CONTIGUOUS :: comm_send(:,:), comm_recv(:,:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_calc_inv(:,:), locks_com_inv(:), locks_calc_fw(:,:), locks_com_fw(:), locks_calc_1(:,:), locks_calc_2(:,:), locks_calc_finished(:,:)

  INTEGER, SAVE :: rem_size
  LOGICAL, SAVE :: first = .true.

  INTEGER :: sendsize, sendsize_rem, nthreads, nested_threads
  INTEGER :: i, j, iter, ierr, buffer_size, batch_size

!  IF( dfft%autotune ) THEN
!     CALL Autotune_Batchsize( dfft, psi, hpsi, nthreads, nested_threads, ngms, nbnd, max_batch_size, buffer_size, batch_size, max_buffer_size )
!  END IF
  !$omp parallel
  dfft%cpus_per_task = omp_get_num_threads()
  !$omp end parallel
  nthreads = MIN( 2, dfft%cpus_per_task )
  nested_threads = MAX( 1, dfft%cpus_per_task - 1 )

  buffer_size = 3
  batch_size  = 1
  sendsize = MAXVAL ( dfft%nr3p ) * MAXVAL( dfft%nsw ) * dfft%node_task_size * dfft%node_task_size * batch_size
  sendsize_rem = MAXVAL ( dfft%nr3p ) * MAXVAL( dfft%nsw ) * dfft%node_task_size * dfft%node_task_size * rem_size

  IF( first ) THEN

     first = .false.
     CALL Set_Req_Vals( dfft, nbnd, batch_size, rem_size, buffer_size )
     CALL Prep_Copy_Maps( dfft, ngms, batch_size, rem_size )
   
   
     CALL create_shared_memory_window_2d( comm_send, i, dfft, sendsize*dfft%nodes_numb, buffer_size ) 
     CALL create_shared_memory_window_2d( comm_recv, i, dfft, sendsize*dfft%nodes_numb, buffer_size ) 
   
!     CALL create_shared_locks_2d( locks_calc_inv, i, dfft, dfft%node_task_size, ( nbnd / batch_size ) + 1 )
!     CALL create_shared_locks_2d( locks_calc_fw,  i, dfft, dfft%node_task_size, ( nbnd / batch_size ) + 1 )
!     CALL create_shared_locks_1d( locks_com_inv,  i, dfft, ( nbnd / batch_size ) + 1 )
!     CALL create_shared_locks_1d( locks_com_fw,   i, dfft, ( nbnd / batch_size ) + 1 )
!     
!     CALL create_shared_locks_2d( locks_calc_1  , i, dfft, dfft%node_task_size, nbnd + batch_size + (buffer_size-1)*batch_size )
!     CALL create_shared_locks_2d( locks_calc_2  , i, dfft, dfft%node_task_size, nbnd + batch_size + (buffer_size-1)*batch_size )

  END IF

  locks_calc_inv = .true.
  locks_calc_fw  = .true.
  locks_com_inv  = .true.
  locks_com_fw   = .true.

  locks_calc_1   = .true.
  DO i = 1, batch_size*buffer_size
     locks_calc_1( : , i ) = .false.
  ENDDO
  locks_calc_2   = .true.

  dfft%first_loading = .true.
  dfft%rem = .false.

  CALL MPI_BARRIER(dfft%comm, ierr)

  CALL Placeholder_Name( dfft, psi, hpsi, v, comm_send, comm_recv, locks_calc_inv, locks_com_inv, locks_calc_fw, locks_com_fw, locks_calc_1, locks_calc_2, dfft%first_step, nthreads, nested_threads, batch_size, rem_size, ngms, nbnd, sendsize, sendsize_rem, buffer_size )

!  CALL clean_up_shared( dfft )
!  CALL MapVals_CleanUp( dfft )

END SUBROUTINE Overlapp_FFT

SUBROUTINE Autotune_Batchsize( dfft, psi, hpsi, nthreads, nested_threads, ngms, nbnd, max_batch_size, lowest_buffer_size, lowest_batch_size, max_buffer_size )

  IMPLICIT NONE
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(IN)  :: psi (ngms, nbnd)
  COMPLEX(DP), INTENT(OUT) :: hpsi(ngms, nbnd)
  INTEGER, INTENT(IN)    :: nested_threads, ngms, nbnd, max_batch_size, max_buffer_size
  INTEGER, INTENT(INOUT) :: nthreads
  INTEGER, INTENT(OUT)   :: lowest_batch_size, lowest_buffer_size
  COMPLEX(DP), POINTER, CONTIGUOUS :: comm_send(:,:), comm_recv(:,:)
  LOGICAL, POINTER, CONTIGUOUS :: locks_calc_inv(:,:), locks_com_inv(:), locks_calc_fw(:,:), locks_com_fw(:), locks_calc_1(:,:), locks_calc_2(:,:), locks_calc_finished(:,:)

  INTEGER :: i, k, ierr, sendsize, tr, p, m, rem_size, sendsize_rem, num_buff, lowest_batch_size_ar( max_buffer_size )
  INTEGER(INT64) :: time(2), total_time( max_batch_size, max_buffer_size ), cr, lowest_batch_size_val( max_buffer_size )
  LOGICAL        :: total_time_MASK( max_batch_size, max_buffer_size )
  REAL(DP) :: total_time_av, total_time_acp

  !Try the variable-buffer way

  IF( dfft%my_node .eq. 0 .and. dfft%my_node_rank .eq. 0 ) write(6,*) "Trying the variable way:"
  total_time = 0
  total_time_MASK = .true.
!  CALL clean_up_shared( dfft )
!  CALL MapVals_CleanUp( dfft )
  CALL SYSTEM_CLOCK( count_rate = cr )

  DO num_buff = max_buffer_size, max_buffer_size

     DO i = max_batch_size, 1, -1
        IF(i*num_buff .gt. (nbnd+1)/2 ) THEN
           total_time_MASK(i,num_buff) = .false.
           CYCLE
        END IF

        IF( dfft%my_node .eq. 0 ) dfft%writ = .true.
        hpsi = (0.0d0,0.0d0)
        CALL Set_Req_Vals( dfft, nbnd, i, rem_size, num_buff )
        CALL Prep_Copy_Maps( dfft, ngms, i, rem_size )

        sendsize     = MAXVAL ( dfft%nr3p ) * MAXVAL( dfft%nsw ) * dfft%node_task_size * dfft%node_task_size * i
        CALL create_shared_memory_window_2d( comm_send, tr, dfft, sendsize*dfft%nodes_numb, num_buff ) 
        CALL create_shared_memory_window_2d( comm_recv, tr, dfft, sendsize*dfft%nodes_numb, num_buff ) 
        sendsize_rem = MAXVAL ( dfft%nr3p ) * MAXVAL( dfft%nsw ) * dfft%node_task_size * dfft%node_task_size * rem_size

!        CALL create_shared_locks_2d( locks_calc_inv, tr, dfft, dfft%node_task_size, ( nbnd / i ) + 1 )
!        CALL create_shared_locks_2d( locks_calc_fw , tr, dfft, dfft%node_task_size, ( nbnd / i ) + 1 )
!        CALL create_shared_locks_1d( locks_com_inv , tr, dfft, ( nbnd / i ) + 1 )
!        CALL create_shared_locks_1d( locks_com_fw  , tr, dfft, ( nbnd / i ) + 1 )
!
!        CALL create_shared_locks_2d( locks_calc_1  , tr, dfft, dfft%node_task_size, nbnd + i + (num_buff-1)*i )
!        CALL create_shared_locks_2d( locks_calc_2  , tr, dfft, dfft%node_task_size, nbnd + i + (num_buff-1)*i )

        DO k = 1, 4

           locks_calc_inv = .true.
           locks_calc_fw  = .true.
           locks_com_inv  = .true.
           locks_com_fw   = .true.

           locks_calc_1   = .true.
           DO m = 1, num_buff*i
              locks_calc_1( : , m ) = .false.
           ENDDO
           locks_calc_2   = .true.

           dfft%first_loading = .true.
           dfft%rem = .false.

           CALL MPI_BARRIER(dfft%comm, ierr)

           CALL SYSTEM_CLOCK( time( 1 ) )

!           CALL Placeholder_Name( dfft, psi, hpsi, comm_send, comm_recv, locks_calc_inv, locks_com_inv, locks_calc_fw, locks_com_fw, locks_calc_1, locks_calc_2, dfft%first_step, nthreads, nested_threads, i , rem_size, ngms, nbnd, sendsize, sendsize_rem, num_buff )

           CALL SYSTEM_CLOCK( time( 2 ) )

           CALL MPI_BARRIER(dfft%comm, ierr)

           dfft%writ = .false.

           IF( k .eq. 1 ) CYCLE

           total_time(i,num_buff) = total_time(i,num_buff) + ( time(2) - time(1) )

        END DO

!        DO m = 1, nbnd
!           DO p = 1, ngms
!              IF( ABS(REAL(psi(p,m)*2.d0-hpsi(p,m),KIND = real64)) > eps8 ) THEN
!                 write(*,*)  "FFT ERROR (Autotune)", dfft%mype, 2.d0*psi(p,m),hpsi(p,m), p,m, i
!                 !STOP
!              end IF
!              IF( hpsi(p,m) .ne. hpsi(p,m) ) write(*,*)  "FFT ERROR (Autotune)", dfft%mype, hpsi(p,m)
!           END DO
!        END do

!        CALL clean_up_shared( dfft )
        CALL MapVals_CleanUp( dfft )

        total_time_av = REAL( total_time(i,num_buff) / 3, KIND = REAL64 ) / REAL ( cr , KIND = REAL64 )

        CALL MPI_ALLREDUCE( total_time_av, total_time_acp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dfft%comm, ierr)

        total_time_acp = total_time_acp / ( dfft%nodes_numb * dfft%node_task_size )

        IF( dfft%my_node .eq. 0 .and. dfft%my_node_rank .eq. 0 ) write(6,*) "Buffer_size", num_buff, "Batch_size", i, "Time:", total_time_acp

     END DO

     IF( dfft%my_node .eq. 0 .and. dfft%my_node_rank .eq. 0 ) write(6,*) " "

     lowest_batch_size_ar( num_buff ) = MINLOC( total_time(:,num_buff), 1, total_time_MASK(:,num_buff) )
     lowest_batch_size_val( num_buff ) = MINVAL( total_time(:,num_buff), 1, total_time_MASK(:,num_buff) )

  END DO

  ! Choose Buffer- and Batchsize

  IF( dfft%my_node .eq. 0 .and. dfft%my_node_rank .eq. 0 ) THEN

!     lowest_buffer_size = MINLOC( lowest_batch_size_val, 1 )
     lowest_buffer_size = max_buffer_size
     lowest_batch_size  = lowest_batch_size_ar( lowest_buffer_size )

  END IF
  CALL MPI_BCAST( lowest_batch_size, 1, MPI_INTEGER, 0, dfft%comm, ierr )
  CALL MPI_BCAST( lowest_buffer_size   , 1, MPI_INTEGER, 0, dfft%comm, ierr )

END SUBROUTINE Autotune_Batchsize


SUBROUTINE Placeholder_Name( dfft, psi, hpsi, v, comm_send, comm_recv, locks_calc_inv, locks_com_inv, locks_calc_fw, locks_com_fw, locks_calc_1, locks_calc_2, first_step, nthreads, nested_threads, batch_size_save, rem_size, ngms, nbnd_source, sendsize_save, sendsize_rem, num_buff )
  IMPLICIT NONE


  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(IN)  :: psi (ngms, nbnd_source)
  COMPLEX(DP), INTENT(INOUT) :: hpsi(ngms, nbnd_source)
  INTEGER, INTENT(IN)    :: nested_threads, ngms, nbnd_source, sendsize_save, rem_size, sendsize_rem, batch_size_save, num_buff
  INTEGER, INTENT(INOUT) :: nthreads
  COMPLEX(DP), INTENT(INOUT) :: comm_send(:,:), comm_recv(:,:)
  LOGICAL, INTENT(INOUT)     :: locks_calc_inv(:,:), locks_com_inv(:), locks_calc_fw(:,:), locks_com_fw(:), locks_calc_1(:,:), locks_calc_2(:,:)
  LOGICAL, INTENT(INOUT)  :: first_step( num_buff )
  REAL(DP), INTENT(IN) :: v(:)

  INTEGER :: my_thread_num, ibatch, batch_size, sendsize, nbnd, ierr
  INTEGER :: next, max_nbnd, work_buffer, last_buffer
  INTEGER :: counter( 2 , 3)
  INTEGER(INT64):: b_tim_n(100)
  LOGICAL :: do_calc, do_com, finished_all, last_start, last_start_triggered, uneven
  LOGICAL :: finished( 2 , 3 )

  INTEGER(INT64) :: time(100) 

  !$  IF( dfft%my_node_rank .ne. 0 .or. dfft%single_node ) THEN
  !$     nthreads = 1
!  !$     CALL omp_set_num_threads( 5 )
  !$  END IF
  !$omp parallel IF( nthreads .eq. 2 ) num_threads( nthreads ) &
  !$omp private( my_thread_num, ibatch, batch_size, sendsize, do_calc, do_com, next, counter, last_buffer, finished_all, work_buffer, first_step, finished, uneven, time ) &
  !$omp proc_bind( close )
  !$  my_thread_num = omp_get_thread_num()
  !$  IF( my_thread_num .eq. 1 ) THEN
  !$     CALL omp_set_max_active_levels( 2 )
  !$     CALL omp_set_num_threads( nested_threads )
  !$     CALL mkl_set_dynamic(0) 
  !$     ierr = mkl_set_num_threads_local( nested_threads )
  !$  END IF

  IF( .not. omp_in_parallel() ) my_thread_num = 1

  batch_size = batch_size_save
  sendsize = sendsize_save
  next = -1
  counter = 0
  finished = .false.
  last_buffer = 0
  do_calc = .false.
  do_com  = .false.
  first_step = .true.
  finished_all = .false.
  last_start = .true.
  last_start_triggered = .false.
  uneven  = .false.
  IF( my_thread_num .eq. 1 ) do_calc = .true.
  IF( .not. dfft%single_node .and. ( my_thread_num .eq. 0 .or. ( nthreads .eq. 1 .and. dfft%my_node_rank .eq. 0 ) ) ) do_com = .true.

!  IF( .not. dfft%lgamma ) THEN
!     nbnd = nbnd_source
!  ELSE
     nbnd = ( nbnd_source + 1 ) / 2
     IF( mod( nbnd_source, 2 ) .ne. 0 ) uneven = .true.
!  END IF

  IF( rem_size .ne. 0 ) THEN
     max_nbnd = ( nbnd / batch_size ) + 1
  ELSE
     max_nbnd = ( nbnd / batch_size )
  END IF

  DO WHILE( .not. finished_all )

     next = mod( next, num_buff ) + 1
     IF( next .eq. 0 ) THEN
        work_buffer = 1
     ELSE
        work_buffer = dfft%buffer_sequence( next )
     END IF

     IF( .false. ) THEN
        IF( dfft%writ .and. dfft%my_node_rank .eq. 0 .and. do_calc .and. counter(1,1) .lt. 10 .and. batch_size_save .eq. 1 ) THEN
           IF( first_step( work_buffer ) ) THEN
              IF( counter(1,1) .ge. num_buff ) write(6,*) "FFT     3", " Buffer ", work_buffer, " base-state ", counter(1,3)*batch_size_save+1
              write(6,*) "FFT/COM 1", " Buffer ", work_buffer, " base-state ", counter(1,1)*batch_size_save+1
           ELSE
              write(6,*) "FFT/COM 2", " Buffer ", work_buffer, " base-state ", counter(1,2)*batch_size_save+1
           END IF
        END IF
     END IF

     IF( rem_size .ne. 0 .and. ( ( first_step( work_buffer ) .and. last_buffer .eq. 0 .and. ( counter( 1, 1 ) .eq. max_nbnd - 1 .or. counter( 2, 1 ) .eq. max_nbnd - 1 ) ) .or. last_buffer .eq. work_buffer ) ) THEN
        last_buffer = work_buffer
        batch_size = rem_size
        sendsize = sendsize_rem
        IF( do_calc ) dfft%rem = .true.
     END IF

     IF( first_step( work_buffer ) ) THEN

        ! Last and First Step

        CALL SYSTEM_CLOCK( time(1) )

        IF( do_calc ) THEN

           ! Last Step

           IF( dfft%first_loading( work_buffer ) ) THEN
              dfft%first_loading( work_buffer ) = .false.
           ELSE

              counter( 1, 3 ) = counter( 1, 3 ) + 1
        
              CALL SYSTEM_CLOCK( time(10) )

              !$omp flush( locks_com_fw )
              !$  DO WHILE( locks_com_fw( counter(1,3) ) .and. .not. dfft%single_node )
              !$omp flush( locks_com_fw )
              !$  END DO

              CALL SYSTEM_CLOCK( time(11) )

              IF( batch_size .ne. batch_size_save .and. last_start ) THEN
                 last_start = .false.
                 last_start_triggered = .true.
                 batch_size = batch_size_save
              END IF

              DO ibatch = 1, batch_size
                 CALL fwfft_after_com( dfft, comm_recv(:,work_buffer), ibatch, batch_size )
!                 IF( .not. dfft%lgamma ) THEN
!                    CALL Accumulate_Psi_overlapp( dfft, hpsi(:,ibatch+((counter(1,3)-1)*batch_size_save)), ibatch, ngms, batch_size, 1 )
!                 ELSE
                    IF( counter( 1, 3 ) .eq. max_nbnd .and. ibatch .eq. batch_size .and. uneven ) THEN
                       CALL Accumulate_Psi_overlapp( dfft, hpsi(:,1+((ibatch-1)+((counter(1,3)-1)*batch_size_save))*2), ibatch, ngms, batch_size, 1 )
                    ELSE
                       CALL Accumulate_Psi_overlapp( dfft, hpsi(:,1+((ibatch-1)+((counter(1,3)-1)*batch_size_save))*2:2+((ibatch-1)+((counter(1,3)-1)*batch_size_save))*2), ibatch, ngms, batch_size, 2 )
                    END IF
!                 END IF

                 !$  locks_calc_1( dfft%my_node_rank+1, ibatch+((counter(1,3)+num_buff-1)*batch_size_save) ) = .false.
                 !$omp flush( locks_calc_1 )

              ENDDO

              IF( last_start_triggered ) batch_size = rem_size

              IF( counter( 1, 3 ) .eq. max_nbnd ) finished_all = .true.

           END IF

           ! First Step

           IF( finished( 1, 1 ) ) THEN
              CONTINUE
           ELSE

              counter( 1, 1 ) = counter( 1, 1 ) + 1
   
              DO ibatch = 1, batch_size

                 CALL SYSTEM_CLOCK( time(30+(ibatch-1)*2) )

                 IF( batch_size .eq. rem_size ) THEN
                    !$omp flush( locks_calc_1 )
                    !$  DO WHILE( ANY(locks_calc_1( :, 1+((counter(1,1)-1)*batch_size_save):batch_size_save+((counter(1,1)-1)*batch_size_save) ) ) )
                    !$omp flush( locks_calc_1 )
                    !$  END DO
                 ELSE
                    !$omp flush( locks_calc_1 )
                    !$  DO WHILE( ANY(locks_calc_1( :, ibatch+((counter(1,1)-1)*batch_size_save) ) ) )
                    !$omp flush( locks_calc_1 )
                    !$  END DO
                 END IF

                 CALL SYSTEM_CLOCK( time(31+(ibatch-1)*2) )
 
!                 IF( .not. dfft%lgamma ) THEN
!                    CALL Prepare_Psi_overlapp( dfft, psi(:,ibatch+((counter(1,1)-1)*batch_size_save)), ibatch, ngms, batch_size, 1 )
!                 ELSE
                    IF( counter( 1, 1 ) .eq. max_nbnd .and. ibatch .eq. batch_size .and. uneven ) THEN
                       CALL Prepare_Psi_overlapp( dfft, psi(:,1+((ibatch-1)+((counter(1,1)-1)*batch_size_save))*2), ibatch, ngms, batch_size, 1 )
                    ELSE
                       CALL Prepare_Psi_overlapp( dfft, psi(:,1+((ibatch-1)+((counter(1,1)-1)*batch_size_save))*2:2+((ibatch-1)+((counter(1,1)-1)*batch_size_save))*2), ibatch, ngms, batch_size, 2 )
                    END IF
!                 END IF
                 CALL invfft_pre_com( dfft, comm_send(:,work_buffer), comm_recv(:,work_buffer), ibatch, batch_size )
              ENDDO

              !$  locks_calc_inv( dfft%my_node_rank+1, counter(1,1) ) = .false.
              !$omp flush( locks_calc_inv )

              IF( counter( 1, 1 ) .eq. max_nbnd ) finished( 1, 1 ) = .true.

           END IF
   
        END IF

        CALL SYSTEM_CLOCK( time(2) )
   
        IF( do_com ) THEN
   
           IF( finished( 2, 1 ) ) THEN
              CONTINUE
           ELSE

              counter( 2, 1 ) = counter( 2, 1 ) + 1

              CALL SYSTEM_CLOCK( time(20) )

              !$omp flush( locks_calc_inv )
              !$  DO WHILE( ANY( locks_calc_inv( :, counter(2,1) ) ) )
              !$omp flush( locks_calc_inv )
              !$  END DO

              CALL SYSTEM_CLOCK( time(7) )

              CALL fft_com( dfft, comm_send(:,work_buffer), comm_recv(:,work_buffer), sendsize, dfft%my_node_rank, &
                            dfft%inter_node_comm, dfft%nodes_numb, dfft%my_inter_node_rank, dfft%non_blocking )

              CALL SYSTEM_CLOCK( time(8) )

              !$  locks_com_inv( counter(2,1) ) = .false.
              !$omp flush( locks_com_inv )


              IF( counter( 2, 1 ) .eq. max_nbnd ) finished( 2, 1 ) = .true.

           END IF
   
        END IF

        IF( dfft%single_node ) CALL MPI_BARRIER(dfft%comm, ierr) 

        first_step( work_buffer ) = .false.


        IF( do_calc ) THEN
           dfft%time_adding(30) = dfft%time_adding(30) + ( time(2) - time(1) )
           dfft%time_adding(31) = dfft%time_adding(31) + ( time(11) - time(10) )
           DO ibatch = 1, batch_size
              dfft%time_adding(32) = dfft%time_adding(32) + ( time(31+(ibatch-1)*2) - time(30+(ibatch-1)*2) )
           ENDDO
        END IF
        IF( do_com ) THEN
           dfft%time_adding(100) = dfft%time_adding(100) + ( time(8) - time(7) )
           dfft%time_adding(90) = dfft%time_adding(90) + ( time(7) - time(20) )
        END IF

     ELSE


        ! Middle Step

        CALL SYSTEM_CLOCK( time(4) )

        IF( do_calc ) THEN
   
           IF( finished( 1, 2 ) ) THEN
              CONTINUE
           ELSE

              counter( 1, 2 ) = counter( 1, 2 ) + 1

              CALL SYSTEM_CLOCK( time(12) )

              !$omp flush( locks_com_inv )
              !$  DO WHILE( locks_com_inv( counter(1,2) ) .and. .not. dfft%single_node )
              !$omp flush( locks_com_inv )
              !$  END DO

              CALL SYSTEM_CLOCK( time(13) )

              DO ibatch = 1, batch_size
                 CALL invfft_after_com( dfft, dfft%bench_aux(:,ibatch), comm_recv(:,work_buffer), ibatch )
                 CALL Apply_V( dfft, dfft%bench_aux(:,ibatch), v )

                 !$  locks_calc_2( dfft%my_node_rank+1, ibatch+((counter(1,2)-1)*batch_size_save) ) = .false.
                 !$omp flush( locks_calc_2 )

              ENDDO
      
              DO ibatch = 1, batch_size

                 CALL SYSTEM_CLOCK( time(50+(ibatch-1)*2) )

                 IF( batch_size .eq. rem_size ) THEN
                    !$omp flush( locks_calc_2 )
                    !$  DO WHILE( ANY(locks_calc_2( :, 1+((counter(1,2)-1)*batch_size_save):batch_size+((counter(1,2)-1)*batch_size_save) ) ) )
                    !$omp flush( locks_calc_2 )
                    !$  END DO
                 ELSE
                    !$omp flush( locks_calc_2 )
                    !$  DO WHILE( ANY(locks_calc_2( :, ibatch+((counter(1,2)-1)*batch_size_save) ) ) )
                    !$omp flush( locks_calc_2 )
                    !$  END DO
                 END IF

                 CALL SYSTEM_CLOCK( time(51+(ibatch-1)*2) )

                 CALL fwfft_pre_com( dfft, dfft%bench_aux(:,ibatch), comm_send(:,work_buffer), comm_recv(:,work_buffer), ibatch, batch_size )
              ENDDO

              !$  locks_calc_fw( dfft%my_node_rank+1, counter(1,2) ) = .false.
              !$omp flush( locks_calc_fw )

              IF( counter( 1, 2 ) .eq. max_nbnd ) finished( 1, 2 ) = .true.

           END IF

        END IF

        CALL SYSTEM_CLOCK( time(5) )
   
        IF( do_com ) THEN

           IF( finished( 2, 2 ) ) THEN
              CONTINUE
           ELSE

              counter( 2, 2 ) = counter( 2, 2 ) + 1

              CALL SYSTEM_CLOCK( time(21) )

              !$omp flush( locks_calc_fw )
              !$  DO WHILE( ANY( locks_calc_fw( :, counter(2,2) ) ) )
              !$omp flush( locks_calc_fw )
              !$  END DO

              CALL SYSTEM_CLOCK( time(9) )
   
              CALL fft_com( dfft, comm_send(:,work_buffer), comm_recv(:,work_buffer), sendsize, dfft%my_node_rank, &
                            dfft%inter_node_comm, dfft%nodes_numb, dfft%my_inter_node_rank, dfft%non_blocking )

              CALL SYSTEM_CLOCK( time(10) )

              !$  locks_com_fw( counter(2,2) ) = .false.
              !$omp flush( locks_com_fw )

              IF( counter( 2, 2 ) .eq. max_nbnd ) THEN
                 finished( 2, 2 ) = .true.
                 IF( .not. do_calc ) THEN
                    finished_all = .true.
                 END IF
              END IF

           END IF
   
        END IF

        IF( dfft%single_node ) CALL MPI_BARRIER(dfft%comm, ierr) 

        first_step( work_buffer ) = .true.

        CALL SYSTEM_CLOCK( time(6) )

        IF( do_calc ) THEN
           dfft%time_adding(30) = dfft%time_adding(30) + ( time(5) - time(4) )
           dfft%time_adding(33) = dfft%time_adding(33) + ( time(13) - time(12) )
           DO ibatch = 1, batch_size
              dfft%time_adding(34) = dfft%time_adding(34) + ( time(51+(ibatch-1)*2) - time(50+(ibatch-1)*2) )
           ENDDO
        END IF
        IF( do_com ) THEN
           dfft%time_adding(100) = dfft%time_adding(100) + ( time(10) - time(9) )
           dfft%time_adding(91) = dfft%time_adding(91) + ( time(9) - time(21) )
        END IF

     END IF

     IF( batch_size .ne. batch_size_save ) THEN
        batch_size = batch_size_save
        sendsize = sendsize_save
        IF( do_calc ) dfft%rem = .false.
     END IF
 

  ENDDO 

  !$  IF( my_thread_num .eq. 1 ) THEN
  !$     CALL omp_set_max_active_levels( 1 )
  !$     CALL omp_set_num_threads( dfft%cpus_per_task )
  !$     CALL mkl_set_dynamic(1) 
  !$     ierr = mkl_set_num_threads_local( dfft%cpus_per_task )
  !$  END IF
  !$omp barrier
  !$omp end parallel

END SUBROUTINE Placeholder_Name

SUBROUTINE Compare_invffts( dfft, c0, psi, ngms, nbnd)
  IMPLICIT NONE
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(IN)  :: c0 ( ngms, nbnd )
  INTEGER, INTENT(IN)    :: ngms, nbnd

  COMPLEX(DP), POINTER, SAVE, CONTIGUOUS :: comm_send(:,:), comm_recv(:,:)

  COMPLEX(DP) :: psi( dfft%nnr, 2 )
  INTEGER, SAVE :: rem_size
  INTEGER :: sendsize, sendsize_rem, nthreads, nested_threads
  INTEGER :: i, j, iter, ierr, buffer_size, batch_size
  



  buffer_size = 1
  batch_size  = 1
  CALL Set_Req_Vals( dfft, nbnd, batch_size, rem_size, buffer_size )
  CALL Prep_Copy_Maps( dfft, ngms, batch_size, rem_size )


  sendsize = MAXVAL ( dfft%nr3p ) * MAXVAL( dfft%nsw ) * dfft%node_task_size * dfft%node_task_size * batch_size
  CALL create_shared_memory_window_2d( comm_send, i, dfft, sendsize*dfft%nodes_numb, buffer_size ) 
  CALL create_shared_memory_window_2d( comm_recv, i, dfft, sendsize*dfft%nodes_numb, buffer_size ) 

  CALL Prepare_Psi_single( dfft, c0( : ,1:2 ) )
  
  CALL invfft_single( dfft, psi( :, 1 ), comm_send(:,1), comm_recv(:,1), sendsize)

  CALL Prepare_Psi_single( dfft, c0( : ,3:4 ) )
  
  CALL invfft_single( dfft, psi( :, 2 ), comm_send(:,1), comm_recv(:,1), sendsize)

END SUBROUTINE Compare_invffts

!=----------------------------------------------------------------------=
END MODULE fftpw_overlapp
!=----------------------------------------------------------------------=
