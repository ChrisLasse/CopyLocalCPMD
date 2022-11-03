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

  INTEGER, SAVE :: times_called = 0
  INTEGER, SAVE :: bb_counter = 1
  INTEGER, SAVE  :: mbuff, mbatch, repeats, omit
  REAL(DP), SAVE, ALLOCATABLE :: final_time(:)

  times_called = times_called + 1
  finished = .false.

  IF( times_called .eq. 1 ) THEN 
     mbuff = 3
     mbatch = 10
     repeats = 10
     omit = 3
     ALLOCATE( final_time( mbuff*mbatch ) )
     final_time = 0.d0
  END IF
  IF( mod( times_called-1, repeats ) .gt. omit-1 ) final_time( bb_counter ) = final_time( bb_counter ) + time
  IF( mod( times_called, repeats ) .eq. 0 ) THEN
  
     final_time( bb_counter ) = final_time( bb_counter ) / ( repeats - omit )
     IF( dfft%mype .eq. 0 ) write(6,*) "vpsi tunning buffer:", dfft%buffer_size_save, "batch:", dfft%batch_size_save, "TIME:", final_time( bb_counter )

     dfft%buffer_size_save = ( bb_counter / mbatch ) + 1
     dfft%batch_size_save  = mod( bb_counter, mbatch ) + 1
     bb_counter = bb_counter + 1
  
     IF( dfft%buffer_size_save .eq. mbuff+1 ) THEN
        IF( dfft%mype .eq. 0 ) THEN
           dfft%buffer_size_save = ( ( MINLOC( final_time, 1 ) - 1 ) / mbatch ) + 1
           dfft%batch_size_save  = mod( MINLOC( final_time, 1 ) - 1, mbatch ) + 1
           write(6,*) " "
           write(6,*) "VPSI TUNNING FINISHED"
           write(6,*) "USED BUFFERSIZE:", dfft%buffer_size_save, "USED BATCHSIZE:", dfft%batch_size_save
           write(6,*) " "
           write(6,*) final_time
           write(6,*) " "
        END IF
        CALL MP_BCAST( dfft%buffer_size_save, 0, dfft%comm )
        CALL MP_BCAST(  dfft%batch_size_save, 0, dfft%comm )
        finished = .true.
     END IF
  
  END IF

  END SUBROUTINE autotune_pw_vpsi

END MODULE autotune_utils
