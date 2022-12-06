#include "cpmd_global.h"

MODULE autotune_utils
  USE elct,                            ONLY: crge
  USE spin,                            ONLY: clsd
  USE rhoofr_utils,                    ONLY: rhoofr_batchfft
  USE vpsi_utils,                      ONLY: vpsi_batchfft
  USE rnlsm_utils,                     ONLY: rnlsm
  USE fft,                             ONLY: batch_fft,&
                                             fft_tune_max_it
  USE fftprp_utils,                    ONLY: autotune_fftbatchsize
  USE fftpw_base,                      ONLY: dfft
  USE fftpw_param,                     ONLY: DP
  USE mp_interface,                    ONLY: mp_bcast
  USE system,                          ONLY: cnti, cntl
  USE rswfmod,                         ONLY: rsactive
  USE kinds,                           ONLY: real_8
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE iso_fortran_env

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: autotune
  PUBLIC :: autotune_pw_vpsi
  
CONTAINS
  ! ==================================================================
  SUBROUTINE autotune(c0,c2,rhoe,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! == Run all iterations for autotuning rnlsm/vpsi/rhoofr          ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(IN) __CONTIGUOUS  :: c0(:,:)
    REAL(real_8),INTENT(OUT) __CONTIGUOUS    :: rhoe(:,:)
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: psi(:,:), c2(:,:)
    INTEGER, INTENT(IN)                      :: nstate

    CHARACTER(*), PARAMETER                  :: procedureN = 'autotuning'

    INTEGER                                  :: isub, it, num_it

    CALL tiset(procedureN,isub)
    
    IF(cntl%fft_tune_batchsize.OR.cntl%rnlsm_autotune)THEN
       num_it=MAX(fft_tune_max_it,cnti%rnlsm_autotune_maxit)
       DO it=1,num_it
          IF(it.LE.cnti%rnlsm_autotune_maxit.AND.cntl%overlapp_comm_comp)THEN
             CALL rnlsm(c0,nstate,1,1,.FALSE.,unpack_dfnl_fnl=.FALSE.)
          END IF
          IF(it.LE.fft_tune_max_it.AND.batch_fft)THEN
             rsactive = cntl%krwfn
             CALL autotune_fftbatchsize()
             CALL rhoofr_batchfft(c0,rhoe,psi(:,1),nstate)
             CALL vpsi_batchfft(c0,c2,crge%f(:,1),rhoe,psi(:,1),nstate,1,clsd%nlsd,.TRUE.)
             rsactive = .FALSE.
          END IF
       END DO
    END IF
    CALL autotune_fftbatchsize()
    cntl%rnlsm_autotune=.FALSE.
    CALL tihalt(procedureN,isub)
  END SUBROUTINE autotune

  SUBROUTINE autotune_pw_vpsi( time, finished )

  REAL(DP), INTENT(IN) :: time
  LOGICAL, INTENT(OUT) :: finished

  INTEGER, SAVE :: bb_counter = 1
  INTEGER, SAVE  :: mbuff, mbatch, repeats, omit, set_repeats
  REAL(DP), SAVE, ALLOCATABLE :: final_time(:), set_time(:,:), save_time(:)
  LOGICAL, SAVE, ALLOCATABLE :: set_time_mask(:)
  INTEGER, SAVE, ALLOCATABLE :: sets(:,:)
  LOGICAL, SAVE :: first = .true.
  LOGICAL, SAVE :: alloc = .true.
  LOGICAL, SAVE :: set_timing = .true.
  INTEGER, SAVE :: timed = 0
  INTEGER, SAVE :: waiting = 0
  LOGICAL :: next_batch

  INTEGER, SAVE :: set_called = 0
  INTEGER(INT64), SAVE :: cr

  finished = .false.
  next_batch = .false.

  IF( first ) THEN 
     mbuff = 3
     mbatch = dfft%max_batch_size
     repeats = 5
     set_repeats = 3
     omit = 3
     ALLOCATE( final_time( mbuff*mbatch ) )
     ALLOCATE( save_time( mbatch ) )
     ALLOCATE( set_time( mbatch, 4 ) )
     ALLOCATE( set_time_mask( mbatch ) )
     ALLOCATE( sets( mbuff*mbatch, 4 ) )
     final_time = 0.d0
     set_time = 0.d0
     save_time = 0.d0
     set_time_mask = .false.
     first = .false.
     CALL SYSTEM_CLOCK( count_rate = cr )
  END IF

  IF( alloc .and. waiting .le. omit-1 ) THEN
     waiting = waiting + 1
     dfft%auto_timings = 0
  ELSE
     alloc = .false.
     IF( set_timing ) THEN
        set_called = set_called + 1
        save_time( dfft%z_set_size_save ) = save_time( dfft%z_set_size_save ) + time 
        set_time( dfft%z_set_size_save, 1 ) = set_time( dfft%z_set_size_save, 1 ) + REAL( dfft%auto_timings(1) / REAL( cr ) )
        set_time( dfft%y_set_size_save, 2 ) = set_time( dfft%y_set_size_save, 2 ) + REAL( dfft%auto_timings(2) / REAL( cr ) )
        set_time( dfft%scatter_set_size_save, 3 ) = set_time( dfft%scatter_set_size_save, 3 ) + REAL( dfft%auto_timings(3) / REAL( cr ) )
        dfft%auto_timings = 0
        set_time_mask( dfft%z_set_size_save ) = .true.
        IF( set_called .gt. set_repeats-1 ) THEN
           IF( dfft%mype .eq. 0 .and. .true. ) THEN
              write(6,*) "vpsi z set tunning:", dfft%z_set_size_save, "TIME:", set_time( dfft%z_set_size_save, 1 ) / set_repeats
              write(6,*) "vpsi y set tunning:", dfft%y_set_size_save, "TIME:", set_time( dfft%y_set_size_save, 2 ) / set_repeats
              write(6,*) "vpsi scatter set tunning:", dfft%scatter_set_size_save, "TIME:", set_time( dfft%scatter_set_size_save, 3 ) / set_repeats
           END IF
           IF( dfft%z_set_size_save .ge. (dfft%batch_size_save+1) / 2 ) THEN
              IF( dfft%z_set_size_save .eq. dfft%batch_size_save ) THEN
                 IF( dfft%mype .eq. 0 ) THEN
                    sets( bb_counter, 1 ) = MINLOC( set_time(:,1), 1, set_time_mask )
                    sets( bb_counter, 2 ) = MINLOC( set_time(:,2), 1, set_time_mask )
                    sets( bb_counter, 3 ) = MINLOC( set_time(:,3), 1, set_time_mask )
                    dfft%z_set_size_save = MINLOC( set_time(:,1), 1, set_time_mask )
                    dfft%y_set_size_save = MINLOC( set_time(:,2), 1, set_time_mask )
                    dfft%scatter_set_size_save = MINLOC( set_time(:,3), 1, set_time_mask )
                    IF( .true. ) THEN
                       write(6,*) "vpsi z set chosen:", dfft%z_set_size_save, "TIME:", set_time( dfft%z_set_size_save, 1 ) / set_repeats
                       write(6,*) "vpsi y set chosen:", dfft%y_set_size_save, "TIME:", set_time( dfft%y_set_size_save, 2 ) / set_repeats
                       write(6,*) "vpsi scatter set chosen:", dfft%scatter_set_size_save, "TIME:", set_time( dfft%scatter_set_size_save, 3 ) / set_repeats
                    END IF
                 END IF
                 CALL MP_BCAST( dfft%z_set_size_save, 0, dfft%comm )
                 CALL MP_BCAST( dfft%y_set_size_save, 0, dfft%comm )
                 CALL MP_BCAST( dfft%scatter_set_size_save, 0, dfft%comm )
                 set_timing = .false.
              ELSE
                 set_called = 0
                 dfft%z_set_size_save = dfft%batch_size_save
                 dfft%y_set_size_save = dfft%batch_size_save
                 dfft%scatter_set_size_save = dfft%batch_size_save
              END IF
           ELSE
              set_called = 0
              dfft%z_set_size_save = dfft%z_set_size_save + 1
              dfft%y_set_size_save = dfft%y_set_size_save + 1
              dfft%scatter_set_size_save = dfft%scatter_set_size_save + 1
           END IF
        END IF
     ELSE
        final_time( bb_counter ) = final_time( bb_counter ) + time
        timed = timed + 1
        IF( timed .eq. repeats ) next_batch = .true.
     END IF
  END IF

  IF( next_batch ) THEN
 
     final_time( bb_counter ) = final_time( bb_counter ) / ( repeats ) !+ set_repeats )
     IF( dfft%mype .eq. 0 ) write(6,*) "vpsi tunning buffer:", dfft%buffer_size_save, "batch:", dfft%batch_size_save, "TIME:", final_time( bb_counter )

     dfft%buffer_size_save = ( bb_counter / mbatch ) + 1
     dfft%batch_size_save  = mod( bb_counter, mbatch ) + 1
     bb_counter = bb_counter + 1
  
     IF( dfft%buffer_size_save .eq. mbuff+1 ) THEN
        IF( dfft%mype .eq. 0 ) THEN
           dfft%buffer_size_save = ( ( MINLOC( final_time, 1 ) - 1 ) / mbatch ) + 1
           dfft%batch_size_save  = mod( MINLOC( final_time, 1 ) - 1, mbatch ) + 1
           dfft%z_set_size_save = sets( dfft%batch_size_save + (dfft%buffer_size_save-1) * mbatch, 1 )
           dfft%y_set_size_save = sets( dfft%batch_size_save + (dfft%buffer_size_save-1) * mbatch, 2 )
           dfft%scatter_set_size_save = sets( dfft%batch_size_save + (dfft%buffer_size_save-1) * mbatch, 3 )
!           dfft%x_set_size_save = sets( dfft%batch_size_save + (dfft%buffer_size_save-1) * mbatch, 4 )
           write(6,*) " "
           write(6,*) "VPSI TUNNING FINISHED"
           write(6,*) "USED BUFFERSIZE:", dfft%buffer_size_save, "USED BATCHSIZE:", dfft%batch_size_save
           write(6,*) "USED z_set_size:", dfft%z_set_size_save
           write(6,*) "USED y_set_size:", dfft%y_set_size_save
           write(6,*) "USED scatter_set_size:", dfft%scatter_set_size_save
!           write(6,*) "USED x_set_size:", dfft%x_set_size_save
           write(6,*) " "
           write(6,*) final_time
           write(6,*) " "
        END IF
        CALL MP_BCAST( dfft%buffer_size_save, 0, dfft%comm )
        CALL MP_BCAST(  dfft%batch_size_save, 0, dfft%comm )
        CALL MP_BCAST(  dfft%z_set_size_save, 0, dfft%comm )
        CALL MP_BCAST(  dfft%y_set_size_save, 0, dfft%comm )
        CALL MP_BCAST(  dfft%scatter_set_size_save, 0, dfft%comm )
!        CALL MP_BCAST(  dfft%x_set_size_save, 0, dfft%comm )
        finished = .true.
     END IF
     alloc = .true. 
     dfft%z_set_size_save = 1
     dfft%y_set_size_save = 1
     dfft%scatter_set_size_save = 1
!     dfft%x_set_size_save = 1
     set_time = 0.d0
     save_time = 0.d0
     set_time_mask = .false.
     set_called = 0
     set_timing = .true.
     timed = 0
     waiting = 0
  
  END IF

  END SUBROUTINE autotune_pw_vpsi

END MODULE autotune_utils
