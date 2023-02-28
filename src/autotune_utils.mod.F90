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
  USE fftpw_make_maps,                 ONLY: GIMME_GROUP_SIZES
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
  PUBLIC :: autotune_pw_4Svpsi
  
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
  LOGICAL, SAVE :: auto_write = .true.

  INTEGER, SAVE :: set_called = 0
  INTEGER(INT64), SAVE :: cr

  finished = .false.
  next_batch = .false.

  IF( first ) THEN 
     mbuff = dfft%max_buffer_size
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
        set_time( dfft%apply_set_size_save, 4 ) = set_time( dfft%apply_set_size_save, 4 ) + REAL( dfft%auto_timings(4) / REAL( cr ) )
        dfft%auto_timings = 0
        set_time_mask( dfft%z_set_size_save ) = .true.
        IF( set_called .gt. set_repeats-1 ) THEN
           IF( dfft%mype .eq. 0 .and. auto_write ) THEN
              write(6,*) "vpsi z set tunning:", dfft%z_set_size_save, "TIME:", set_time( dfft%z_set_size_save, 1 ) / set_repeats
              write(6,*) "vpsi y set tunning:", dfft%y_set_size_save, "TIME:", set_time( dfft%y_set_size_save, 2 ) / set_repeats
              write(6,*) "vpsi scatter set tunning:", dfft%scatter_set_size_save, "TIME:", set_time( dfft%scatter_set_size_save, 3 ) / set_repeats
              write(6,*) "vpsi apply set tunning:", dfft%apply_set_size_save, "TIME:", set_time( dfft%apply_set_size_save, 4 ) / set_repeats
           END IF
           IF( dfft%z_set_size_save .ge. (dfft%batch_size_save+1) / 2 ) THEN
              IF( dfft%z_set_size_save .eq. dfft%batch_size_save ) THEN
                 IF( dfft%mype .eq. 0 ) THEN
                    sets( bb_counter, 1 ) = MINLOC( set_time(:,1), 1, set_time_mask )
                    sets( bb_counter, 2 ) = MINLOC( set_time(:,2), 1, set_time_mask )
                    sets( bb_counter, 3 ) = MINLOC( set_time(:,3), 1, set_time_mask )
                    sets( bb_counter, 4 ) = MINLOC( set_time(:,4), 1, set_time_mask )
                    dfft%z_set_size_save = MINLOC( set_time(:,1), 1, set_time_mask )
                    dfft%y_set_size_save = MINLOC( set_time(:,2), 1, set_time_mask )
                    dfft%scatter_set_size_save = MINLOC( set_time(:,3), 1, set_time_mask )
                    dfft%apply_set_size_save = MINLOC( set_time(:,4), 1, set_time_mask )
                    IF( auto_write ) THEN
                       write(6,*) "vpsi z set chosen:", dfft%z_set_size_save, "TIME:", set_time( dfft%z_set_size_save, 1 ) / set_repeats
                       write(6,*) "vpsi y set chosen:", dfft%y_set_size_save, "TIME:", set_time( dfft%y_set_size_save, 2 ) / set_repeats
                       write(6,*) "vpsi scatter set chosen:", dfft%scatter_set_size_save, "TIME:", set_time( dfft%scatter_set_size_save, 3 ) / set_repeats
                       write(6,*) "vpsi apply set chosen:", dfft%apply_set_size_save, "TIME:", set_time( dfft%apply_set_size_save, 4 ) / set_repeats
                    END IF
                 END IF
                 CALL MP_BCAST( dfft%z_set_size_save, 0, dfft%comm )
                 CALL MP_BCAST( dfft%y_set_size_save, 0, dfft%comm )
                 CALL MP_BCAST( dfft%scatter_set_size_save, 0, dfft%comm )
                 CALL MP_BCAST( dfft%apply_set_size_save, 0, dfft%comm )
                 set_timing = .false.
              ELSE
                 set_called = 0
                 dfft%z_set_size_save = dfft%batch_size_save
                 dfft%y_set_size_save = dfft%batch_size_save
                 dfft%scatter_set_size_save = dfft%batch_size_save
                 dfft%apply_set_size_save = dfft%batch_size_save
              END IF
           ELSE
              set_called = 0
              dfft%z_set_size_save = dfft%z_set_size_save + 1
              dfft%y_set_size_save = dfft%y_set_size_save + 1
              dfft%scatter_set_size_save = dfft%scatter_set_size_save + 1
              dfft%apply_set_size_save = dfft%apply_set_size_save + 1
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
           dfft%apply_set_size_save = sets( dfft%batch_size_save + (dfft%buffer_size_save-1) * mbatch, 4 )
           write(6,*) " "
           write(6,*) "VPSI TUNNING FINISHED"
           write(6,*) "USED BUFFERSIZE:", dfft%buffer_size_save, "USED BATCHSIZE:", dfft%batch_size_save
           write(6,*) "USED z_set_size:", dfft%z_set_size_save
           write(6,*) "USED y_set_size:", dfft%y_set_size_save
           write(6,*) "USED scatter_set_size:", dfft%scatter_set_size_save
!           write(6,*) "USED x_set_size:", dfft%x_set_size_save
           write(6,*) "USED apply_set_size:", dfft%apply_set_size_save
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
        CALL MP_BCAST(  dfft%apply_set_size_save, 0, dfft%comm )
        finished = .true.
     END IF
     alloc = .true. 
     dfft%z_set_size_save = 1
     dfft%y_set_size_save = 1
     dfft%scatter_set_size_save = 1
!     dfft%x_set_size_save = 1
     dfft%apply_set_size_save = 1
     set_time = 0.d0
     save_time = 0.d0
     set_time_mask = .false.
     set_called = 0
     set_timing = .true.
     timed = 0
     waiting = 0
  
  END IF

  END SUBROUTINE autotune_pw_vpsi

  SUBROUTINE autotune_pw_4Svpsi_old( time, finished )

  REAL(DP), INTENT(IN) :: time
  LOGICAL, INTENT(OUT) :: finished

  INTEGER, SAVE :: bb_counter = 1
  INTEGER, SAVE  :: mbuff, mbatch, repeats, omit, set_repeats
  REAL(DP), SAVE, ALLOCATABLE :: final_time(:), set_z_time(:), set_yx_time(:), save_time(:)
  INTEGER, SAVE, ALLOCATABLE :: sets(:,:), set_z_save(:), set_yx_save(:,:)
  LOGICAL, SAVE, ALLOCATABLE :: set_time_mask(:,:)
  LOGICAL, SAVE :: first = .true.
  LOGICAL, SAVE :: alloc = .true.
  LOGICAL, SAVE :: set_timing = .true.
  INTEGER, SAVE :: timed = 0
  INTEGER, SAVE :: waiting = 0
  LOGICAL :: next_batch
  LOGICAL, SAVE :: auto_write = .true.
  INTEGER, SAVE :: which_z  = 1
  INTEGER, SAVE :: which_yx = 1
  INTEGER, SAVE :: z_repeats = 0

  INTEGER :: yx_loc, z_loc, y_group_size

  INTEGER, SAVE :: set_called = 0
  INTEGER(INT64), SAVE :: cr

  finished = .false.
  next_batch = .false.

  IF( first ) THEN 
     mbuff = dfft%max_buffer_size
     mbatch = dfft%max_batch_size
     repeats = 5
     set_repeats = 3
     omit = 3
     ALLOCATE( final_time( mbuff*mbatch ) )
     ALLOCATE( save_time( mbatch ) )
     ALLOCATE( set_z_time( mbatch*mbatch ) )
     ALLOCATE( set_z_save( mbatch ) )
     ALLOCATE( set_yx_time( mbatch*mbatch ) )
     ALLOCATE( set_yx_save( mbatch*mbatch, 2 ) )
     ALLOCATE( sets( mbuff*mbatch, 3 ) )
     ALLOCATE( set_time_mask( mbatch*mbatch, 2 ) )
     set_time_mask = .false.
     final_time = 0.d0
     set_z_time = 0.d0
     set_yx_time = 0.d0
     save_time = 0.d0
     set_z_save = 0
     set_yx_save = 0
     first = .false.
     CALL SYSTEM_CLOCK( count_rate = cr )
  END IF

  IF( alloc .and. waiting .le. omit-1 ) THEN
     waiting = waiting + 1
     dfft%auto_4Stimings = 0
  ELSE
     alloc = .false.
     IF( set_timing ) THEN
        set_called = set_called + 1
        z_repeats = z_repeats + 1
        set_z_time( which_z ) = set_z_time( which_z ) + REAL( dfft%auto_4Stimings(1) / REAL( cr ) )
        set_z_save( which_z ) = dfft%z_set_size_save
        set_yx_time( which_yx ) = set_yx_time( which_yx ) + REAL( dfft%auto_4Stimings(2) / REAL( cr ) )
        set_yx_save( which_yx, 1 ) = dfft%y_set_size_save
        set_yx_save( which_yx, 2 ) = dfft%x_set_size_save
        set_time_mask( which_z, 1 ) = .true.
        set_time_mask( which_yx, 2 ) = .true.
        dfft%auto_4Stimings = 0
        IF( mod( dfft%batch_size_save, dfft%y_set_size_save ) .eq. 0 ) THEN 
           y_group_size = dfft%batch_size_save / dfft%y_set_size_save
        ELSE
           y_group_size = dfft%batch_size_save / dfft%y_set_size_save + 1
        END IF
        IF( set_called .gt. set_repeats-1 ) THEN
           IF( dfft%mype .eq. 0 .and. auto_write ) THEN
              write(6,'(A7,I2,A23,I6,A6,E15.8)') "batch:", dfft%batch_size_save, "; z set tunning:", dfft%z_set_size_save, "TIME:", set_z_time( which_z ) / z_repeats
              write(6,'(A7,I2,A23,I3,I3,A6,E15.8)') "batch:", dfft%batch_size_save, "; y and x set tunning:", dfft%y_set_size_save, dfft%x_set_size_save, "TIME:", set_yx_time( which_yx ) / set_repeats
           END IF
           IF( dfft%x_set_size_save .ge. (y_group_size+1) / 2 ) THEN
              IF( dfft%x_set_size_save .eq. y_group_size ) THEN 
                 IF( dfft%y_set_size_save .ge. (dfft%batch_size_save+1) / 2 ) THEN
                    IF( dfft%y_set_size_save .eq. dfft%batch_size_save ) THEN
                       IF( dfft%mype .eq. 0 ) THEN
                          z_loc = MINLOC( set_z_time, 1, set_time_mask( : , 1 ) )
                          dfft%z_set_size_save = set_z_save( z_loc )
                          sets( bb_counter, 1 ) = set_z_save( z_loc )
                          yx_loc = MINLOC( set_yx_time, 1, set_time_mask( : , 2 ) )
                          dfft%y_set_size_save = set_yx_save( yx_loc, 1 )
                          dfft%x_set_size_save = set_yx_save( yx_loc, 2 )
                          sets( bb_counter, 2 ) = set_yx_save( yx_loc, 1 )
                          sets( bb_counter, 3 ) = set_yx_save( yx_loc, 2 )
                          IF( auto_write ) THEN
                             write(6,'(A7,I2,A23,I6,A6,E15.8)') "batch:", dfft%batch_size_save, "; z set chosen:", dfft%z_set_size_save, "TIME:", set_z_time( z_loc ) / set_repeats
                             write(6,'(A7,I2,A23,I3,I3,A6,E15.8)') "batch:", dfft%batch_size_save, "; y and x set chosen:", dfft%y_set_size_save, dfft%x_set_size_save, "TIME:", set_yx_time( yx_loc ) / set_repeats
                          END IF
                       END IF
                       CALL MP_BCAST( dfft%z_set_size_save, 0, dfft%comm )
                       CALL MP_BCAST( dfft%y_set_size_save, 0, dfft%comm )
                       CALL MP_BCAST( dfft%x_set_size_save, 0, dfft%comm )
                       set_timing = .false.
                    ELSE
                       set_z_time( which_z ) = set_z_time( which_z ) / ( z_repeats / set_repeats )
                       set_called = 0
                       z_repeats = 0
                       which_z = which_z + 1
                       which_yx = which_yx + 1
                       dfft%z_set_size_save = dfft%batch_size_save
                       dfft%y_set_size_save = dfft%batch_size_save
                       dfft%x_set_size_save = 1
                    END IF
                 ELSE
                    set_z_time( which_z ) = set_z_time( which_z ) / ( z_repeats / set_repeats )
                    set_called = 0
                    z_repeats = 0
                    which_z = which_z + 1
                    which_yx = which_yx + 1
                    dfft%z_set_size_save = dfft%z_set_size_save + 1
                    dfft%y_set_size_save = dfft%y_set_size_save + 1
                    dfft%x_set_size_save = 1
                 END IF
              ELSE
                 set_called = 0
                 which_yx = which_yx + 1
                 dfft%x_set_size_save = y_group_size
              END IF
           ELSE
              set_called = 0
              which_yx = which_yx + 1
              dfft%x_set_size_save = dfft%x_set_size_save + 1
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
     IF( dfft%mype .eq. 0 ) write(6,'(A21,I3,A7,I3,A6,E15.8)') "vpsi tunning buffer:", dfft%buffer_size_save, "batch:", dfft%batch_size_save, "TIME:", final_time( bb_counter )

     dfft%buffer_size_save = ( bb_counter / mbatch ) + 1
     dfft%batch_size_save  = mod( bb_counter, mbatch ) + 1
     bb_counter = bb_counter + 1
  
     IF( dfft%buffer_size_save .eq. mbuff+1 ) THEN
        IF( dfft%mype .eq. 0 ) THEN
           dfft%buffer_size_save = ( ( MINLOC( final_time, 1 ) - 1 ) / mbatch ) + 1
           dfft%batch_size_save  = mod( MINLOC( final_time, 1 ) - 1, mbatch ) + 1
           dfft%z_set_size_save = sets( dfft%batch_size_save + (dfft%buffer_size_save-1) * mbatch, 1 )
           dfft%y_set_size_save = sets( dfft%batch_size_save + (dfft%buffer_size_save-1) * mbatch, 2 )
           dfft%x_set_size_save = sets( dfft%batch_size_save + (dfft%buffer_size_save-1) * mbatch, 3 )
           write(6,*) " "
           write(6,*) "VPSI TUNNING FINISHED"
           write(6,*) "USED BUFFERSIZE:", dfft%buffer_size_save, "USED BATCHSIZE:", dfft%batch_size_save
           write(6,*) "USED z_set_size:", dfft%z_set_size_save
           write(6,*) "USED y_set_size:", dfft%y_set_size_save
           write(6,*) "USED x_set_size:", dfft%x_set_size_save
           write(6,*) " "
           write(6,*) final_time
           write(6,*) " "
        END IF
        CALL MP_BCAST( dfft%buffer_size_save, 0, dfft%comm )
        CALL MP_BCAST(  dfft%batch_size_save, 0, dfft%comm )
        CALL MP_BCAST(  dfft%z_set_size_save, 0, dfft%comm )
        CALL MP_BCAST(  dfft%y_set_size_save, 0, dfft%comm )
        CALL MP_BCAST(  dfft%x_set_size_save, 0, dfft%comm )
        finished = .true.
     ELSE
        alloc = .true. 
        dfft%z_set_size_save = 1
        dfft%y_set_size_save = 1
        dfft%x_set_size_save = 1
        set_z_time = 0.d0
        set_yx_time = 0.d0
        save_time = 0.d0
        set_called = 0
        z_repeats = 0
        set_timing = .true.
        timed = 0
        waiting = 0
        which_z = 1
        which_yx = 1
        set_time_mask = .false.
        set_z_save = 0
        set_yx_save = 0
     END IF
  
  END IF

  END SUBROUTINE autotune_pw_4Svpsi_old

  SUBROUTINE autotune_pw_4Svpsi( time )

  REAL(DP), INTENT(IN) :: time

  INTEGER, SAVE  :: mbuff, mbatch, repeats, omit, set_repeats
  REAL(DP), SAVE, ALLOCATABLE :: final_time(:), z_group_time(:), y_group_time(:), x_group_time(:)
  LOGICAL, SAVE, ALLOCATABLE :: group_mask(:,:)
  LOGICAL, SAVE :: first = .true.
  LOGICAL, SAVE :: alloc = .true.
  LOGICAL, SAVE :: group_timing = .true.
  INTEGER, SAVE :: timed = 0
  INTEGER, SAVE :: waiting = 0
  LOGICAL :: next_batch
  LOGICAL, SAVE :: auto_write = .true.
  INTEGER, SAVE :: x_repeats = 0
  INTEGER, SAVE :: y_z_repeats = 0
  INTEGER, SAVE :: group_repeats = 0

  INTEGER :: bb_counter, i, j

  INTEGER, SAVE :: group_called = 0
  INTEGER(INT64), SAVE :: cr

  next_batch = .false.

  IF( first ) THEN 
     mbuff = dfft%max_buffer_size
     mbatch = dfft%max_batch_size
     repeats = 5
     group_repeats = 3
     omit = 3
     ALLOCATE( final_time( mbuff*mbatch ) )
     ALLOCATE( z_group_time( mbatch ) )
     ALLOCATE( y_group_time( mbatch ) )
     ALLOCATE( x_group_time( mbatch ) )
     ALLOCATE( group_mask( mbatch, 3 ) )
     final_time = 0.d0
     z_group_time = 0.d0
     y_group_time = 0.d0
     x_group_time = 0.d0
     first = .false.
     CALL SYSTEM_CLOCK( count_rate = cr )
  END IF

  IF( alloc .and. waiting .le. omit-1 ) THEN
     waiting = waiting + 1
     dfft%auto_4Stimings = 0
  ELSE
     bb_counter = (dfft%buffer_size_save-1) * dfft%max_batch_size + dfft%batch_size_save
     alloc = .false.
     IF( group_timing .and. dfft%batch_size_save .eq. dfft%max_batch_size ) THEN
        z_group_time( dfft%z_group_autosize ) = z_group_time( dfft%z_group_autosize ) + REAL( REAL( dfft%auto_4Stimings(1) / ( dfft%z_divider * dfft%z_group_autosize ) ) / REAL( cr ) )
        y_group_time( dfft%y_group_autosize ) = y_group_time( dfft%y_group_autosize ) + REAL( REAL( dfft%auto_4Stimings(2) / ( dfft%y_divider * dfft%y_group_autosize ) ) / REAL( cr ) )
        y_z_repeats = y_z_repeats + 1
        IF( .not. dfft%x_optimal ) THEN
           x_group_time( dfft%x_group_autosize ) = x_group_time( dfft%x_group_autosize ) + REAL( REAL( dfft%auto_4Stimings(3) / ( dfft%x_divider * dfft%x_group_autosize ) ) / REAL( cr ) )
           x_repeats = x_repeats + 1
        END IF
        dfft%auto_4Stimings = 0
        group_called = group_called + 1
        IF( group_called .gt. group_repeats-1 ) THEN
           group_called = 0
           IF( dfft%y_group_autosize .eq. dfft%max_batch_size .and. dfft%x_group_autosize .ne. 1 ) THEN
              x_group_time( dfft%x_group_autosize ) = x_group_time( dfft%x_group_autosize ) / x_repeats
              IF( dfft%mype .eq. 0 ) write(6,'(A7,I3,A17,I3,A6,E15.8)') "batch:", dfft%batch_size_save ,"tunning x_group:", dfft%x_group_autosize, "TIME:", x_group_time( dfft%x_group_autosize )
              dfft%x_group_autosize = dfft%x_group_autosize - 1
              x_repeats = 0
           ELSE
              IF( .not. dfft%x_optimal ) THEN
                 x_group_time( dfft%x_group_autosize ) = x_group_time( dfft%x_group_autosize ) / x_repeats
                 IF( dfft%mype .eq. 0 ) write(6,'(A7,I3,A17,I3,A6,E15.8)') "batch:", dfft%batch_size_save ,"tunning x_group:", dfft%x_group_autosize, "TIME:", x_group_time( dfft%x_group_autosize )
                 group_mask = .true.
                 x_repeats = 0
                 IF( dfft%mype .eq. 0 ) THEN
                    DO i = 1, dfft%max_batch_size
                       dfft%optimal_groups( i , dfft%buffer_size_save , 1 ) = MINLOC( x_group_time, 1, group_mask(:,1) )
                       group_mask( dfft%optimal_groups( i, dfft%buffer_size_save, 1 ) , 1 ) = .false.
                    ENDDO
                 END IF
                 DO i = 1, dfft%max_batch_size
                    CALL MP_BCAST( dfft%optimal_groups(i,dfft%buffer_size_save,1), 0, dfft%comm )
                 ENDDO
                 dfft%x_optimal = .true.
              END IF
              y_group_time( dfft%y_group_autosize ) = y_group_time( dfft%y_group_autosize ) / y_z_repeats
              IF( dfft%mype .eq. 0 ) write(6,'(A7,I3,A17,I3,A6,E15.8)') "batch:", dfft%batch_size_save ,"tunning y_group:", dfft%y_group_autosize, "TIME:", y_group_time( dfft%y_group_autosize )
              z_group_time( dfft%z_group_autosize ) = z_group_time( dfft%z_group_autosize ) / y_z_repeats
              IF( dfft%mype .eq. 0 ) write(6,'(A7,I3,A17,I3,A6,E15.8)') "batch:", dfft%batch_size_save ,"tunning z_group:", dfft%z_group_autosize, "TIME:", z_group_time( dfft%z_group_autosize )
              y_z_repeats = 0
              IF( dfft%y_group_autosize .gt. 1 ) THEN
                 dfft%y_group_autosize = dfft%y_group_autosize - 1
                 dfft%z_group_autosize = dfft%z_group_autosize - 1
              ELSE
                 IF( dfft%mype .eq. 0 ) THEN
                    DO i = 1, dfft%max_batch_size
                       dfft%optimal_groups( i , dfft%buffer_size_save , 2 ) = MINLOC( y_group_time, 1, group_mask(:,2) )
                       group_mask( dfft%optimal_groups( i, dfft%buffer_size_save, 2 ) , 2 ) = .false.
                       dfft%optimal_groups( i , dfft%buffer_size_save , 3 ) = MINLOC( z_group_time, 1, group_mask(:,3) )
                       group_mask( dfft%optimal_groups( i, dfft%buffer_size_save, 3 ) , 3 ) = .false.
                    ENDDO
                 END IF
                 DO i = 1, dfft%max_batch_size
                    DO j = 2, 3
                       CALL MP_BCAST( dfft%optimal_groups(i,dfft%buffer_size_save,j), 0, dfft%comm )
                    ENDDO
                 ENDDO
                 dfft%y_z_optimal = .true.
                 group_timing = .false.
              END IF
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
     IF( dfft%mype .eq. 0 ) write(6,'(A21,I3,A7,I3,A6,E15.8)') "vpsi tunning buffer:", dfft%buffer_size_save, "batch:", dfft%batch_size_save, "TIME:", final_time( bb_counter )

     IF( dfft%batch_size_save .eq. 1 ) THEN
        IF( dfft%buffer_size_save .ne. mbuff ) THEN
           dfft%batch_size_save = dfft%max_batch_size
           dfft%buffer_size_save = dfft%buffer_size_save + 1
           dfft%x_group_autosize = dfft%max_batch_size
           dfft%y_group_autosize = dfft%max_batch_size
           dfft%z_group_autosize = dfft%max_batch_size
           dfft%x_optimal = .false.
           dfft%y_z_optimal = .false.
        ELSE
           dfft%buffer_size_save = dfft%buffer_size_save + 1
        END IF
     ELSE
        dfft%batch_size_save = dfft%batch_size_save - 1
     END IF
  
     IF( dfft%buffer_size_save .eq. mbuff+1 ) THEN
        IF( dfft%mype .eq. 0 ) THEN
           dfft%buffer_size_save = ( ( MINLOC( final_time, 1 ) - 1 ) / mbatch ) + 1
           dfft%batch_size_save  = mod( MINLOC( final_time, 1 ) - 1, mbatch ) + 1
           write(6,*) " "
           write(6,*) "VPSI TUNNING FINISHED"
           write(6,*) " "
           write(6,*) final_time(1:10)
           write(6,*) final_time(11:20)
           IF( mbuff .eq. 3 ) write(6,*) final_time(21:30)
           write(6,*) " "
        END IF
        CALL MP_BCAST( dfft%buffer_size_save, 0, dfft%comm )
        CALL MP_BCAST(  dfft%batch_size_save, 0, dfft%comm )
        dfft%autotune_finished = .true.
     ELSE
        alloc = .true. 
        x_group_time = 0.d0
        y_group_time = 0.d0
        z_group_time = 0.d0
        timed = 0
        waiting = 0
        group_timing = .true.
        y_z_repeats = 0
        x_repeats = 0
     END IF
  
  END IF

  END SUBROUTINE autotune_pw_4Svpsi

END MODULE autotune_utils
