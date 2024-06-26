#if defined(__FFT_HAS_LOW_LEVEL_TIMERS)
#define HAS_LOW_LEVEL_TIMERS .TRUE.
#else
#define HAS_LOW_LEVEL_TIMERS .FALSE.
#endif

#include "cpmd_global.h"

MODULE fftmain_utils





  USE cppt,                            ONLY: indz,&
                                             nzh
  USE cp_cuda_types,                   ONLY: cp_cuda_env
  USE cp_cufft_types,                  ONLY: cp_cufft,&
                                             cp_cufft_device_get_ptrs,&
                                             cp_cufft_plans_t,&
                                             cp_cufft_stream_get_ptrs
  USE cp_cufft_utils,                  ONLY: cp_cufft_get_plan
  USE cublas_types,                    ONLY: cublas_handle_t
  USE cuda_types,                      ONLY: cuda_memory_t,&
                                             cuda_stream_t
  USE cuda_utils,                      ONLY: cuda_memcpy_async_device_to_host,&
                                             cuda_memcpy_async_host_to_device,&
                                             cuda_stream_synchronize
  USE cufft_types,                     ONLY: cufft_plan_t
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: &
       lfrm, lmsq, lr1, lr1m, lr1s, lr2s, lr3s, lrxpl, lsrm, mfrays, msqf, &
       msqs, msrays, qr1, qr1s, qr2s, qr3max, qr3min, qr3s, sp5, sp8, sp9, &
       xf, yf, fft_residual, fft_total, fft_numbatches, fft_batchsize, locks_inv, locks_fw, FFT_TYPE_DESCRIPTOR, fft_buffsize
  USE fft_maxfft,                      ONLY: maxfftn, maxfft
  USE fftcu_methods,                   ONLY: fftcu_frw_full_1,&
                                             fftcu_frw_full_2,&
                                             fftcu_frw_sprs_1,&
                                             fftcu_frw_sprs_2,&
                                             fftcu_inv_full_1,&
                                             fftcu_inv_full_2,&
                                             fftcu_inv_sprs_1,&
                                             fftcu_inv_sprs_2
  USE fftnew_utils,                    ONLY: Prep_fft_com,&
                                             Make_Manual_Maps
  USE fftutil_utils,                   ONLY: fft_comm,&
                                             getz,&
                                             pack_x2y,&
                                             pack_y2x,&
                                             phasen,&
                                             putz,&
                                             unpack_x2y,&
                                             unpack_y2x,&
                                             putz_n,&
                                             getz_n,&
                                             pack_y2x_n,&
                                             unpack_x2y_n,&
                                             pack_x2y_n,&
                                             unpack_y2x,&
                                             unpack_y2x_n,&
                                             locks_omp,&
                                             locks_calc_1,&
                                             locks_calc_2,&
                                             locks_omp_big,&
                                             Prepare_Psi,&
                                             fft_com,&
                                             invfft_z_section,&
                                             invfft_y_section,&
                                             invfft_x_section,&
                                             fwfft_z_section,&
                                             fwfft_y_section,&
                                             fwfft_x_section
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE machine,                         ONLY: m_walltime
  USE mltfft_utils,                    ONLY: mltfft_cuda,&
                                             mltfft_default,&
                                             mltfft_essl,&
                                             mltfft_fftw,&
                                             mltfft_hp
  USE mp_interface,                    ONLY: mp_win_alloc_shared_mem_central, &
                                             mp_win_alloc_shared_mem
  USE parac,                           ONLY: parai
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw
  USE thread_view_types,               ONLY: thread_view_t
  USE timer,                           ONLY: tihalt,&
                                             tiset

  !$ USE omp_lib, ONLY: omp_in_parallel, omp_get_thread_num, &
  !$ omp_set_num_threads, omp_set_nested

  USE, INTRINSIC :: ISO_C_BINDING,     ONLY: C_PTR,&
                                             C_NULL_PTR,&
                                             C_F_POINTER
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  USE iso_fortran_env

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mltfft
  PUBLIC :: invfftn
  PUBLIC :: fwfftn
  PUBLIC :: invfftn_batch
  PUBLIC :: fwfftn_batch

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

  COMPLEX(DP), POINTER, SAVE, CONTIGUOUS :: comm_send(:,:)
  COMPLEX(DP), POINTER, SAVE, CONTIGUOUS :: comm_recv(:,:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_calc_inv(:,:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_com_inv(:,:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_calc_fw(:,:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_com_fw(:,:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_all_com(:,:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_all_com2(:,:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_sing_1(:,:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_sing_2(:,:)

  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_flush(:,:)
  PUBLIC :: locks_flush


  PUBLIC :: comm_send
  PUBLIC :: comm_recv
  PUBLIC :: locks_calc_inv
  PUBLIC :: locks_com_inv
  PUBLIC :: locks_calc_fw
  PUBLIC :: locks_com_fw
  PUBLIC :: locks_all_com
  PUBLIC :: locks_all_com2
  PUBLIC :: locks_sing_1
  PUBLIC :: locks_sing_2

  PUBLIC :: invfft_batch
  PUBLIC :: fwfft_batch
  !public :: fftnew


CONTAINS


  ! ==================================================================
  SUBROUTINE fftnew(isign,f,sparse,comm)
    ! ==--------------------------------------------------------------==
#ifdef __PARALLEL
    USE mpi_f08
#endif
    INTEGER, INTENT(IN)                      :: isign
    COMPLEX(real_8), INTENT(INOUT) &
                               __CONTIGUOUS  :: f(:)
    LOGICAL, INTENT(IN)                      :: sparse
#ifdef __PARALLEL
    type(MPI_COMM), INTENT(IN)               :: comm
#else
    INTEGER, INTENT(IN)                      :: comm
#endif
    CHARACTER(*), PARAMETER                  :: procedureN = 'fftnew'

    COMPLEX(real_8), DIMENSION(:), &
      POINTER __CONTIGUOUS                   :: xf_ptr, yf_ptr
    INTEGER                                  :: lda, m, mm, n1o, n1u, ierr,isub
    INTEGER(int_8)                           :: il_xf(2)
    REAL(real_8)                             :: scale
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN//'get_scratch',isub)
#ifdef _USE_SCRATCHLIBRARY
    il_xf(1)=maxfft
    il_xf(2)=1
    CALL request_scratch(il_xf,xf,procedureN//'_xf',ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL request_scratch(il_xf,yf,procedureN//'_yf',ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
#endif
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN//'get_scratch',isub)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub)
    xf_ptr => xf(:,1)
    yf_ptr => yf(:,1)

    LDA=HUGE(0);MM=HUGE(0);N1U=HUGE(0);N1O=HUGE(0);M=HUGE(0)
    scale=HUGE(0.0_real_8)
    IF (isign.EQ.-1) THEN
       scale=1._real_8
       IF (sparse) THEN
          m=msrays
          CALL mltfft('N','T',f,qr1s,m,xf_ptr,m,qr1s,lr1s,m,isign,scale )
          lda=lsrm*lr1m
          mm=qr2s*(qr3max-qr3min+1)
          CALL pack_x2y(xf_ptr,yf_ptr,msrays,lda,lrxpl,sp5,maxfftn,parai%nproc,cntl%tr4a2a)
          CALL fft_comm(yf_ptr,xf_ptr,lda,cntl%tr4a2a,comm)
          CALL unpack_x2y(xf_ptr,yf_ptr,mm,lr1,lda,msqs,lmsq,sp9,maxfftn,parai%nproc,cntl%tr4a2a)
          m=(qr3max-qr3min+1)*qr1
          CALL mltfft('N','T',yf_ptr,qr2s,m,xf_ptr,m,qr2s,lr2s,m,isign,scale )
          m=qr1*qr2s
          CALL putz(xf_ptr,yf_ptr,qr3min,qr3max,qr3s,m)
          CALL mltfft('N','T',yf_ptr,qr3s,m,f,m,qr3s,lr3s,m,isign,scale )
       ELSE
          m=mfrays
          CALL mltfft('N','T',f,qr1s,m,xf_ptr,m,qr1s,lr1s,m,isign,scale )
          lda=lfrm*lr1m
          mm=qr2s*qr3s
          CALL pack_x2y(xf_ptr,yf_ptr,mfrays,lda,lrxpl,sp5,maxfftn,parai%nproc,cntl%tr4a2a)
          CALL fft_comm(yf_ptr,xf_ptr,lda,cntl%tr4a2a,comm)
          CALL unpack_x2y(xf_ptr,yf_ptr,mm,lr1,lda,msqf,lmsq,sp8,maxfftn,parai%nproc,cntl%tr4a2a)
          m=qr1*qr3s
          CALL mltfft('N','T',yf_ptr,qr2s,m,xf_ptr,m,qr2s,lr2s,m,isign,scale )
          m=qr1*qr2s
          CALL mltfft('N','T',xf_ptr,qr3s,m,f,m,qr3s,lr3s,m,isign,scale )
          n1u=lrxpl(parai%mepos,1)
          n1o=lrxpl(parai%mepos,2)
          CALL phasen(f,qr1,qr2s,qr3s,n1u,n1o,lr2s,lr3s)
       ENDIF
    ELSE
       IF (sparse) THEN
          scale=1._real_8
          m=qr1*qr2s
          CALL mltfft('T','N',f,m,qr3s,xf_ptr,qr3s,m,lr3s,m,isign,scale )
          CALL getz(xf_ptr,f,qr3min,qr3max,qr3s,m)
          m=(qr3max-qr3min+1)*qr1
          CALL mltfft('T','N',f,m,qr2s,yf_ptr,qr2s,m,lr2s,m,isign,scale )
          lda=lsrm*lr1m
          mm=qr2s*(qr3max-qr3min+1)
          CALL pack_y2x(xf_ptr,yf_ptr,mm,lr1,lda,msqs,lmsq,sp9,maxfftn,parai%nproc,cntl%tr4a2a)
          CALL fft_comm(xf_ptr,yf_ptr,lda,cntl%tr4a2a,comm)
          CALL unpack_y2x(xf_ptr,yf_ptr,mm,msrays,lda,lrxpl,sp5,maxfftn,parai%nproc,cntl%tr4a2a)
          scale=1._real_8/REAL(lr1s*lr2s*lr3s,kind=real_8)
          m=msrays
          CALL mltfft('T','N',xf_ptr,m,qr1s,f,qr1s,m,lr1s,m,isign,scale )
       ELSE
          scale=1._real_8
          n1u=lrxpl(parai%mepos,1)
          n1o=lrxpl(parai%mepos,2)
          CALL phasen(f,qr1,qr2s,qr3s,n1u,n1o,lr2s,lr3s)
          m=qr1*qr2s
          CALL mltfft('T','N',f,m,qr3s,xf_ptr,qr3s,m,lr3s,m,isign,scale )
          m=qr1*qr3s
          CALL mltfft('T','N',xf_ptr,m,qr2s,yf_ptr,qr2s,m,lr2s,m,isign,scale )
          lda=lfrm*lr1m
          mm=qr2s*qr3s
          CALL pack_y2x(xf_ptr,yf_ptr,mm,lr1,lda,msqf,lmsq,sp8,maxfftn,parai%nproc,cntl%tr4a2a)
          CALL fft_comm(xf_ptr,yf_ptr,lda,cntl%tr4a2a,comm)
          CALL unpack_y2x(xf_ptr,yf_ptr,mm,mfrays,lda,lrxpl,sp5,maxfftn,parai%nproc,cntl%tr4a2a)
          scale=1._real_8/REAL(lr1s*lr2s*lr3s,kind=real_8)
          m=mfrays
          CALL mltfft('T','N',xf_ptr,m,qr1s,f,qr1s,m,lr1s,m,isign,scale )
       ENDIF
    ENDIF
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN//'release_scratch',isub)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_xf,yf,procedureN//'_yf',ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL free_scratch(il_xf,xf,procedureN//'_xf',ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#endif
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN//'release_scratch',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fftnew

  ! ==================================================================
  SUBROUTINE fftnew_batch(isign,f,n,swap,step,comm,ibatch)
    ! ==--------------------------------------------------------------==
    ! Author:
    ! Tobias Kloeffel, CCC,FAU Erlangen-Nuernberg tobias.kloeffel@fau.de
    ! Gerald Mathias, LRZ, Garching Gerald.Mathias@lrz.de
    ! Bernd Meyer, CCC, FAU Erlangen-Nuernberg bernd.meyer@fau.de
    ! Date March 2019
    ! Special FFT driver that operates on batches of FFTs of multiple
    ! states
    ! This works fine because the 3D-FFT is anyway split into multiple
    ! 1D FFTs
    ! overlapping communication with computation possible and effective
    ! Full performance only with saved arrays or external scratch_library
    ! TODO:
    ! Try MPI_IALLTOALL
    ! Write version for FULL ffts, usefull for transforming batches of
    ! pair densities during HFX calculations!

#ifdef __PARALLEL
    USE mpi_f08
#endif
    IMPLICIT NONE
    INTEGER, INTENT(IN)                      :: isign, n, swap, step, ibatch
    COMPLEX(real_8),INTENT(INOUT) &
                __CONTIGUOUS                 :: f(:)
#ifdef __PARALLEL
    type(MPI_COMM), INTENT(IN)               :: comm
#else
    INTEGER, INTENT(IN)                      :: comm
#endif

    CHARACTER(*), PARAMETER                  :: procedureN = 'fftnew_batch'

    REAL(real_8)                             :: scale
    INTEGER                                  :: lda, m, mm, n1o, n1u
    REAL(real_8)                             :: temp

    LDA=HUGE(0);MM=HUGE(0);N1U=HUGE(0);N1O=HUGE(0);M=HUGE(0)
    scale=HUGE(0.0_real_8)
    IF (isign.EQ.-1) THEN
       scale=1._real_8
       IF(step.EQ.1)THEN
          IF(n.NE.0)THEN
             m=msrays*n
             CALL mltfft('N','T',f,qr1s,m,xf(:,swap),&
                  m,qr1s,lr1s,m,isign,scale )
             lda=lsrm*lr1m
             m=msrays
             CALL pack_x2y_n(xf(:,swap),yf(:,swap),m,lda,lrxpl,sp5,maxfftn,&
                  parai%nproc,cntl%tr4a2a,n)
             !$ locks_inv(ibatch,1) = .FALSE.
             !$omp flush(locks_inv)
          END IF
       ELSE IF(step.EQ.2)THEN
          IF (n.NE.0) THEN
             lda=lsrm*lr1m*n
             !$omp flush(locks_inv)
             !$ DO WHILE ( locks_inv(ibatch,1) )
             !$omp flush(locks_inv)
             !$ END DO
             !$ locks_inv(ibatch,1) = .TRUE.
             !$omp flush(locks_inv)
             CALL fft_comm(yf(:,swap),xf(:,swap),lda,cntl%tr4a2a,comm)
             !$ locks_inv(ibatch,1)=.FALSE.
             !$omp flush(locks_inv)
             !$ locks_inv(ibatch,2)=.FALSE.
             !$omp flush(locks_inv)
          END IF
       ELSE IF(step.EQ.3)THEN
          IF (n.NE.0) THEN
             lda=lsrm*lr1m
             mm=qr2s*(qr3max-qr3min+1)
             m=(qr3max-qr3min+1)*qr1*qr2s
             !$omp flush(locks_inv)
             !$ DO WHILE ( locks_inv(ibatch,2) )
             !$omp flush(locks_inv)
             !$ END DO
             !$omp flush(locks_inv)
             !$ DO WHILE ( locks_inv(ibatch,1) )
             !$omp flush(locks_inv)
             !$ END DO
             CALL unpack_x2y_n(xf(:,swap),yf(:,swap),mm,lr1,lda,msqs,lmsq,sp9,&
                  fpar%nnr1,parai%nproc,cntl%tr4a2a,n,m)
             m=(qr3max-qr3min+1)*qr1*n
             CALL mltfft('N','T',yf(:,swap),qr2s,m,xf(:,swap),m,qr2s,lr2s,m,isign,&
                  scale)
             m=qr1*qr2s
             CALL putz_n(xf(:,swap),yf(:,swap),qr3min,qr3max,qr3s,qr1,qr2s,n)
             m=qr1*qr2s*n
             !copy_out copies the wavefcuntion out in order
             !this is a very time consuming process
             !we better live with mixed order of the wavefunctions and
             !offload that logic into vpsi/rhoofr
             !IF (n.gt.1) THEN
             !   CALL mltfft('N','T',yf(:,swap),qr3s,m,xf(:,swap),m,qr3s,lr3s,m,isign,scale)
             !   lda=qr1*qr2s
             !   m=(loop-1)*fft_batchsize+1
             !   call copy_out(xf(:,swap),lda,qr3s,n,f_out,m)
             !ELSE
             CALL mltfft('N','T',yf(:,swap),qr3s,m,f,m,&
                  qr3s,lr3s,m,isign,scale)
             !END IF
          END IF
       END IF

    ELSE

       IF(step.EQ.1)THEN
          IF(n.NE.0)THEN
             scale=1._real_8
             lda=qr1*qr2s
             m=qr1*qr2s*n
             CALL mltfft('T','N',f,m,qr3s,yf(:,swap),&
                  qr3s,m,lr3s,m,isign,scale )
             CALL getz_n(yf(:,swap),xf(:,swap),qr3min,qr3max,qr3s,qr1,qr2s,n)
             m=(qr3max-qr3min+1)*qr1*n
             CALL mltfft('T','N',xf(:,swap),m,qr2s,yf(:,swap),qr2s,m,lr2s,m,isign,&
                  scale )
             lda=lsrm*lr1m
             mm=qr2s*(qr3max-qr3min+1)
             m=(qr3max-qr3min+1)*qr1*qr2s
             CALL pack_y2x_n(xf(:,swap),yf(:,swap),mm,lr1,lda,msqs,lmsq,sp9,&
                  maxfftn,parai%nproc,cntl%tr4a2a,n,m)
             !$ locks_fw(ibatch,1) = .FALSE.
             !$omp flush(locks_fw)
          END IF
       ELSE IF(step.EQ.2)THEN
          IF(n.NE.0)THEN
             lda=lsrm*lr1m*n
             !$omp flush(locks_fw)
             !$ DO WHILE ( locks_fw(ibatch,1) )
             !$omp flush(locks_fw)
             !$ END DO
             !$ locks_fw(ibatch,1) = .TRUE.
             !$omp flush(locks_fw)
             CALL fft_comm(xf(:,swap),yf(:,swap),lda,cntl%tr4a2a,comm)
             !$ locks_fw(ibatch,1)=.FALSE.
             !$omp flush(locks_fw)
             !$ locks_fw(ibatch,2)=.FALSE.
             !$omp flush(locks_fw)
          END IF
       ELSE IF(step.EQ.3)THEN
          IF(n.NE.0)THEN
             lda=lsrm*lr1m
             mm=qr2s*(qr3max-qr3min+1)
             !$omp flush(locks_fw)
             !$ DO WHILE ( locks_fw(ibatch,2) )
             !$omp flush(locks_fw)
             !$ END DO
             !$omp flush(locks_fw)
             !$ DO WHILE ( locks_fw(ibatch,1) )
             !$omp flush(locks_fw)
             !$ END DO
             CALL unpack_y2x_n(xf(:,swap),yf(:,swap),mm,msrays,lda,lrxpl,sp5,&
                  maxfftn,parai%nproc,cntl%tr4a2a,n)
             scale=1._real_8/REAL(lr1s*lr2s*lr3s,kind=real_8)
             m=msrays*n
             CALL mltfft('T','N',xf(:,swap),m,qr1s,&
                  f(:),qr1s,m,lr1s,m,isign,scale )
          END IF
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fftnew_batch
  ! ==================================================================
  SUBROUTINE fftnew_cuda(isign,f,sparse,comm,thread_view, copy_data_to_device, copy_data_to_host )
    ! ==--------------------------------------------------------------==
#ifdef __PARALLEL
    USE mpi_f08
#endif
    INTEGER, INTENT(IN)                      :: isign
    COMPLEX(real_8)                          :: f(:)
    LOGICAL, INTENT(IN)                      :: sparse
#ifdef __PARALLEL
    type(MPI_COMM), INTENT(IN)                      :: comm
#else
    INTEGER, INTENT(IN)                      :: comm
#endif
    TYPE(thread_view_t), INTENT(IN), &
      OPTIONAL                               :: thread_view
    LOGICAL, INTENT(IN), OPTIONAL            :: copy_data_to_device, &
                                                copy_data_to_host

    CHARACTER(*), PARAMETER                  :: procedureN = 'fftnew_cuda'

    COMPLEX(real_8), POINTER __CONTIGUOUS    :: xf_ptr(:), yf_ptr(:)
    INTEGER                                  :: device_idx, host_buff_ptr, &
                                                lda, stream_idx
    LOGICAL                                  :: copy_to_device, copy_to_host
    REAL(real_8)                             :: scale
    TYPE(cp_cufft_plans_t), POINTER          :: plans_d
    TYPE(cublas_handle_t), POINTER           :: blas_handle_p
    TYPE(cuda_memory_t), POINTER             :: lrxpl_d, msqf_d, msqs_d, &
                                                sp5_d, sp8_d, sp9_d, t1_d, &
                                                t2_d
    TYPE(cuda_stream_t), POINTER             :: stream_p

    copy_to_host = .TRUE.
    copy_to_device = .TRUE.
    IF( PRESENT( copy_data_to_host ) ) copy_to_host = copy_data_to_host
    IF( PRESENT( copy_data_to_device ) ) copy_to_device = copy_data_to_device

    !vw if thread view is available, the proben the device/stream
    !vw else NEED TO CLEAN THAT

    IF( PRESENT( thread_view ) ) THEN
       device_idx = thread_view%device_idx
       stream_idx = thread_view%stream_idx - 1 !vw FIX that -1
       host_buff_ptr = thread_view%host_buff_ptr
    ELSE
       device_idx = 1
       stream_idx = 0
       !$ stream_idx = omp_get_thread_num()
       host_buff_ptr = stream_idx + 1 !vw FIX that +1
    ENDIF

    xf_ptr => xf(:, host_buff_ptr )
    yf_ptr => yf(:, host_buff_ptr )

    CALL cp_cufft_stream_get_ptrs ( cp_cufft, device_idx, stream_idx, stream=stream_p, &
         t1_d=t1_d, t2_d=t2_d, plans=plans_d, blas_handle=blas_handle_p )
    CALL cp_cufft_device_get_ptrs ( cp_cufft, device_idx, sp5_d=sp5_d, sp8_d=sp8_d, &
         sp9_d=sp9_d, msqs_d=msqs_d, msqf_d=msqf_d, lrxpl_d=lrxpl_d ) 

    LDA=HUGE(0)
    scale=HUGE(0.0_real_8)
    IF (isign.EQ.-1) THEN
       scale=1._real_8
       IF (sparse) THEN
          CALL fftcu_inv_sprs_1 ( f, xf_ptr, yf_ptr, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp5_d, lrxpl_d, copy_to_device, plans_d,&
               blas_handle_p, stream_p )

          CALL cuda_stream_synchronize(stream_p)
          lda=lsrm*lr1m
          CALL fft_comm(yf_ptr,xf_ptr,lda,cntl%tr4a2a,comm)

          CALL fftcu_inv_sprs_2 ( f, xf_ptr, yf_ptr, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp9_d, msqs_d, copy_to_host, plans_d,&
               blas_handle_p, stream_p )

       ELSE
          CALL fftcu_inv_full_1 ( f, xf_ptr, yf_ptr, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp5_d, lrxpl_d, copy_to_device, plans_d,&
               blas_handle_p, stream_p )

          CALL cuda_stream_synchronize(stream_p)
          lda=lfrm*lr1m
          CALL fft_comm(yf_ptr,xf_ptr,lda,cntl%tr4a2a,comm)

          CALL fftcu_inv_full_2 ( f, xf_ptr, yf_ptr, parai%mepos, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp8_d, msqf_d, copy_to_host, plans_d,&
               blas_handle_p, stream_p )

       ENDIF
    ELSE
       IF (sparse) THEN
          CALL fftcu_frw_sprs_1 ( f, xf_ptr, yf_ptr, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp9_d, msqs_d, copy_to_device, plans_d,&
               blas_handle_p, stream_p )

          CALL cuda_stream_synchronize(stream_p)
          lda=lsrm*lr1m
          CALL fft_comm(xf_ptr,yf_ptr,lda,cntl%tr4a2a,comm)

          CALL fftcu_frw_sprs_2 ( f, xf_ptr, yf_ptr, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp5_d, lrxpl_d, copy_to_host, plans_d,&
               blas_handle_p, stream_p )

       ELSE
          CALL fftcu_frw_full_1 ( f, xf_ptr, yf_ptr, parai%mepos, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp8_d, msqf_d, copy_to_device, plans_d,&
               blas_handle_p, stream_p )

          CALL cuda_stream_synchronize(stream_p)
          lda=lfrm*lr1m
          CALL fft_comm(xf_ptr,yf_ptr,lda,cntl%tr4a2a,comm)

          CALL fftcu_frw_full_2 ( f, xf_ptr, yf_ptr, parai%nproc, cntl%tr4a2a,&
               t1_d, t2_d, sp5_d, lrxpl_d, copy_to_host, plans_d,&
               blas_handle_p, stream_p )

       ENDIF
    ENDIF

    IF( copy_to_host ) CALL cuda_stream_synchronize(stream_p)

    ! ==--------------------------------------------------------------==
  END SUBROUTINE fftnew_cuda

  ! ==================================================================
  SUBROUTINE mltfft(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,isign,scale,thread_view)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: transa, transb
    COMPLEX(real_8),INTENT(IN) __CONTIGUOUS  :: a(:)
    INTEGER                                  :: ldax, lday
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: b(:)
    INTEGER                                  :: ldbx, ldby, n, m, isign
    REAL(real_8)                             :: scale
    TYPE(thread_view_t), INTENT(IN), &
      OPTIONAL                               :: thread_view

    CHARACTER(*), PARAMETER                  :: procedureN = 'mltfft'

    INTEGER                                  :: device_idx, host_buff_ptr, &
                                                stream_idx, isub
    TYPE(cp_cufft_plans_t), POINTER          :: plans_d
    TYPE(cublas_handle_t), POINTER           :: blas_handle_p
    TYPE(cuda_memory_t), POINTER             :: t1_d, t2_d
    TYPE(cuda_stream_t), POINTER             :: stream_p
    TYPE(cufft_plan_t), POINTER              :: plan_p

! ==--------------------------------------------------------------==
!vw for the moment we use CPU FFT 
!IF( cp_cuda_env%use_fft ) THEN
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub)
    IF( .FALSE. ) THEN

       PRINT *,'START DEBUG ----------------------------------------'

       IF( PRESENT( thread_view ) ) THEN
          device_idx = thread_view%device_idx
          stream_idx = thread_view%stream_idx - 1 !vw FIX that -1
          host_buff_ptr = thread_view%host_buff_ptr
       ELSE
          device_idx = 1
          stream_idx = 0
          !$ IF( omp_in_parallel() ) THEN
          !$    stream_idx = omp_get_thread_num()
          !    !$    WRITE(6,*) 'FFTNEW_CUDA - Stream n.: ',stream_idx ! acm: this is just for debug
          !$ ENDIF
          host_buff_ptr = stream_idx + 1 !vw FIX that +1
       ENDIF

       t1_d => cp_cufft%devices(device_idx)%streams(stream_idx)%t1
       t2_d => cp_cufft%devices(device_idx)%streams(stream_idx)%t2
       plans_d => cp_cufft%devices(device_idx)%streams(stream_idx)%plans
       stream_p => cp_cufft%devices(device_idx)%streams(stream_idx)%stream
       blas_handle_p => cp_cufft%devices(device_idx)%streams(stream_idx)%blas_handle

       PRINT *,'cuda_memcpy_async_host_to_device'

       CALL cuda_memcpy_async_host_to_device ( a(1:maxfftn), t1_d, stream_p )!vw check if maxfftn is correct

       PRINT *,'cp_cufft_get_plan'
       CALL cp_cufft_get_plan ( transa,transb,     ldax,lday,               n,m, plan_p, plans_d, stream_p )

       PRINT *,'mltfft_cuda'
       CALL mltfft_cuda(        transa,transb, t1_d,ldax,lday,t2_d,ldbx,ldby,n,m,isign,scale, plan_p, blas_handle_p, stream_p )

       PRINT *,'cuda_memcpy_async_device_to_host'
       WRITE (6,*) 'yf_d%init',t2_d%init
       WRITE (6,*) 'yf_d%idevice',t2_d%device
       WRITE (6,*) 'yf_d%n_bytes',t2_d%n_bytes
       WRITE (6,*) 'stream_p%init',stream_p%init
       WRITE (6,*) 'stream_p%device',stream_p%device

       CALL cuda_memcpy_async_device_to_host ( t2_d, b(1:maxfftn), stream_p )!vw check if maxfftn is correct

       PRINT *,'cuda_stream_synchronize'
       CALL cuda_stream_synchronize(stream_p)


       PRINT *,'DONE DEBUG ----------------------------------------'
    ELSE

#if defined(__HAS_FFT_DEFAULT)
       CALL mltfft_default(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,isign,scale)
#elif defined(__HAS_FFT_ESSL)
       CALL mltfft_essl(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,isign,scale)
#elif defined(__HAS_FFT_HP)
       CALL mltfft_hp(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,isign,scale)
#elif defined(_HAS_FFT_FFTW3)
       CALL mltfft_fftw(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,isign,scale,.FALSE.)
#else
       CALL stopgm(procedureN,"MLTFFT ROUTINE NOT AVAILABLE",&
            __LINE__,__FILE__)
#endif
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE mltfft
  ! ==================================================================
  ! NEW FFT CODE
  ! ==================================================================
  SUBROUTINE invfftn(f,sparse,comm,thread_view, copy_data_to_device, copy_data_to_host )
    ! ==--------------------------------------------------------------==
    ! == COMPUTES THE INVERSE FOURIER TRANSFORM OF A COMPLEX          ==
    ! == FUNCTION F. THE FOURIER TRANSFORM IS                         ==
    ! == RETURNED IN F (THE INPUT F IS OVERWRITTEN).                  ==
    ! ==--------------------------------------------------------------==
    USE fft,                             ONLY : tfft
#ifdef __PARALLEL
    USE mpi_f08
#endif
    COMPLEX(real_8),INTENT(INOUT) &
                                __CONTIGUOUS :: f(:)
    LOGICAL                                  :: sparse
#ifdef __PARALLEL
    type(MPI_COMM), INTENT(IN)               :: comm
#else
    INTEGER, INTENT(IN)                      :: comm
#endif
    TYPE(thread_view_t), INTENT(IN), &
      OPTIONAL                               :: thread_view
    LOGICAL, INTENT(IN), OPTIONAL            :: copy_data_to_device, &
                                                copy_data_to_host

    CHARACTER(*), PARAMETER                  :: procedureN = 'invfftn'

    INTEGER                                  :: isign, isub

    CALL tiset(procedureN,isub)
    isign=-1
    IF( cp_cuda_env%use_fft ) THEN
       CALL fftnew_cuda(isign,f,sparse, comm, thread_view=thread_view, &
            & copy_data_to_device=copy_data_to_device, copy_data_to_host=copy_data_to_host )
    ELSE
       IF( .false. ) THEN
          CALL fftnew(isign,f,sparse, parai%allgrp )
       ELSE
          CALL fft_improved( isign, tfft, f, tfft%nhg, tfft%nr1p, tfft%ir1p, tfft%nsp )
       END IF
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE invfftn
  ! ==================================================================
  SUBROUTINE fwfftn(f,sparse,comm,thread_view, copy_data_to_device, copy_data_to_host )
    ! ==--------------------------------------------------------------==
    ! == COMPUTES THE FORWARD FOURIER TRANSFORM OF A COMPLEX          ==
    ! == FUNCTION F. THE FOURIER TRANSFORM IS                         ==
    ! == RETURNED IN F IN OUTPUT (THE INPUT F IS OVERWRITTEN).        ==
    ! ==--------------------------------------------------------------==
    USE fft,                             ONLY :tfft 
#ifdef __PARALLEL
    USE mpi_f08
#endif
    COMPLEX(real_8),INTENT(INOUT) &
     __CONTIGUOUS                            :: f(:)
    LOGICAL                                  :: sparse
#ifdef __PARALLEL
    type(MPI_COMM), INTENT(IN)                      :: comm
#else
    INTEGER, INTENT(IN)                      :: comm
#endif
    TYPE(thread_view_t), INTENT(IN), &
      OPTIONAL                               :: thread_view
    LOGICAL, INTENT(IN), OPTIONAL            :: copy_data_to_device, &
                                                copy_data_to_host

    CHARACTER(*), PARAMETER                  :: procedureN = 'fwfftn'

    INTEGER                                  :: isign, isub

    CALL tiset(procedureN,isub)
    isign=1
    IF( cp_cuda_env%use_fft ) THEN
       CALL fftnew_cuda(isign,f,sparse, comm, thread_view=thread_view, &
            & copy_data_to_device=copy_data_to_device, copy_data_to_host=copy_data_to_host )
    ELSE
       IF( .false. ) THEN
          CALL fftnew(isign,f,sparse, parai%allgrp )
       ELSE
          CALL fft_improved( isign, tfft, f, tfft%nhg, tfft%nr1p, tfft%ir1p, tfft%nsp )
       END IF
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fwfftn
  ! ==================================================================
  SUBROUTINE fwfftn_batch(f,len_f,n,swap,step,ibatch)
    ! ==--------------------------------------------------------------==
    ! == COMPUTES THE FORWARD FOURIER TRANSFORM OF A COMPLEX          ==
    ! == FUNCTION F_IN. THE FOURIER TRANSFORM IS                      ==
    ! == RETURNED IN F_OUT                                            ==
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: n, swap, step, ibatch, len_f
    COMPLEX(real_8), INTENT(INOUT)           :: f(len_f)
    CHARACTER(*), PARAMETER                  :: procedureN = 'fwfftn_batch'

    INTEGER                                  :: isign, isub, isub1

    IF(cntl%fft_tune_batchsize) THEN
       CALL tiset(procedureN//'tune',isub1)
    ELSE
       CALL tiset(procedureN,isub)
    END IF
    isign=1
    CALL fftnew_batch(isign,f,n,swap,step,parai%allgrp,ibatch)
    IF(cntl%fft_tune_batchsize) THEN
       CALL tihalt(procedureN//'tune',isub1)
    ELSE
       CALL tihalt(procedureN,isub)
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fwfftn_batch
  ! ==================================================================
  SUBROUTINE invfftn_batch(f,len_f,n,swap,step,ibatch)
    ! ==--------------------------------------------------------------==
    ! == COMPUTES THE INVERSE FOURIER TRANSFORM OF A COMPLEX          ==
    ! == FUNCTION F_IN. THE FOURIER TRANSFORM IS                      ==
    ! == RETURNED IN F_OUT                                            ==
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: n, swap, step, ibatch, len_f
    COMPLEX(real_8), INTENT(INOUT)           :: f(len_f)
    CHARACTER(*), PARAMETER                  :: procedureN = 'invfftn_batch'

    INTEGER                                  :: isign, isub, isub1

    
    IF(cntl%fft_tune_batchsize) THEN
       CALL tiset(procedureN//'tune',isub1)
    ELSE
       CALL tiset(procedureN,isub)
    END IF
    isign=-1
    CALL fftnew_batch(isign,f,n,swap,step,parai%allgrp,ibatch)
    IF(cntl%fft_tune_batchsize) THEN
       CALL tihalt(procedureN//'tune',isub1)
    ELSE
       CALL tihalt(procedureN,isub)
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE invfftn_batch

  SUBROUTINE invfft_batch( tfft, step, batch_size, ispec, remswitch, mythread, counter, work_buffer, f_inout1, f_inout2, f_inout3 )
    IMPLICIT NONE
  
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    INTEGER, INTENT(IN) :: step, batch_size, remswitch, mythread, counter, work_buffer, ispec
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout1(:,:)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout2(:,:)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout3(:,:)

    CHARACTER(*), PARAMETER :: procedureN = 'invfft_batch'

    INTEGER :: isub, isub4

    IF( cntl%fft_tune_batchsize ) THEN
       IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
    ELSE
!       IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
       CALL tiset(procedureN,isub)
    END IF
  
    IF( step .eq. 1 ) THEN
       CALL fft_improved_batch( tfft, -1, step, batch_size, ispec, remswitch, mythread, counter, work_buffer, &
                            f_inout1=f_inout1, f_inout2=f_inout2, f_inout3=f_inout3 )
    ELSE IF( step .eq. 2 ) THEN                                          
       CALL fft_improved_batch( tfft, -1, step, batch_size, ispec, remswitch, mythread, counter, work_buffer )
    ELSE IF( step .eq. 3 ) THEN                                       
       CALL fft_improved_batch( tfft, -1, step, batch_size, ispec, remswitch, mythread, counter, work_buffer, &
                            f_inout1=f_inout1, f_inout2=f_inout2, f_inout3=f_inout3 )
    ELSE IF( step .eq. 4 ) THEN                                       
       CALL fft_improved_batch( tfft, -1, step, batch_size, ispec, remswitch, mythread, counter, work_buffer, &
                            f_inout1=f_inout1, f_inout2=f_inout2, f_inout3=f_inout3 )
    END IF

    IF( cntl%fft_tune_batchsize ) THEN
       IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
    ELSE
!       IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
       CALL tihalt(procedureN,isub)
    END IF
  
  END SUBROUTINE invfft_batch
  
  SUBROUTINE fwfft_batch( tfft, step, batch_size, ispec, remswitch, mythread, counter, work_buffer, f_inout1, f_inout2, f_inout3 )
    IMPLICIT NONE
  
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft 
    INTEGER, INTENT(IN) :: step, batch_size, remswitch, mythread, counter, work_buffer, ispec
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout1(:,:)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout2(:,:)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout3(:,:)

    CHARACTER(*), PARAMETER :: procedureN = 'fwfft_batch'

    INTEGER :: isub, isub4

    IF( cntl%fft_tune_batchsize ) THEN
       IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
    ELSE
!       IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
       CALL tiset(procedureN,isub)
    END IF
  
    IF( step .eq. 1 ) THEN
       CALL fft_improved_batch( tfft, 1, step, batch_size, ispec, remswitch, mythread, counter, work_buffer, &
                            f_inout1=f_inout1, f_inout2=f_inout2, f_inout3=f_inout3 )
    ELSE IF( step .eq. 2 ) THEN
       CALL fft_improved_batch( tfft, 1, step, batch_size, ispec, remswitch, mythread, counter, work_buffer, &
                            f_inout1=f_inout1, f_inout2=f_inout2, f_inout3=f_inout3 )
    ELSE IF( step .eq. 3 ) THEN
       CALL fft_improved_batch( tfft, 1, step, batch_size, ispec, remswitch, mythread, counter, work_buffer )
    ELSE IF( step .eq. 4 ) THEN
       CALL fft_improved_batch( tfft, 1, step, batch_size, ispec, remswitch, mythread, counter, work_buffer, &
                            f_inout1=f_inout1, f_inout2=f_inout2, f_inout3=f_inout3 )
    END IF

    IF( cntl%fft_tune_batchsize ) THEN
       IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
    ELSE
!       IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
       CALL tihalt(procedureN,isub)
    END IF
  
  END SUBROUTINE fwfft_batch

  SUBROUTINE fft_improved_batch( tfft, isign, step, batch_size, ispec, remswitch, mythread, counter, work_buffer, f_inout1, f_inout2, f_inout3 )
    IMPLICIT NONE
  
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) ::tfft 
    INTEGER, INTENT(IN) :: isign, step, batch_size, counter, work_buffer, ispec
    INTEGER, INTENT(IN) :: remswitch, mythread
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout1(:,:)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout2(:,:)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout3(:,:)

    CHARACTER(*), PARAMETER :: procedureN = 'fft_improved_batch'
  
    INTEGER :: current, isub, isub4, ierr
    INTEGER(INT64) :: time(30)

    IF( cntl%fft_tune_batchsize ) THEN
       IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
    ELSE
       IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
    END IF
 
    current = (counter-1)*fft_batchsize
  
    IF( isign .eq. -1 ) THEN !!  invfft
  
       IF( step .eq. 1 ) THEN

          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) CALL SYSTEM_CLOCK( time(1) )

          !In theory, locks could be made faster by checking each iset individually for non-remainder cases 
!          !$omp flush( locks_calc_1 )
!          !$  DO WHILE( ANY(locks_calc_1( :, 1+current:fft_batchsize+current ) ) )
!          !$omp flush( locks_calc_1 )
!          !$  END DO

          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) CALL SYSTEM_CLOCK( time(2) )
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) tfft%time_adding( 17 ) = tfft%time_adding( 17 ) + ( time(2) - time(1) )

          CALL invfft_z_section( tfft, f_inout1, f_inout2(:,work_buffer), f_inout3(:,work_buffer), batch_size, remswitch, mythread, tfft%nsw, current )

          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) CALL SYSTEM_CLOCK( time(3) )

          !$  locks_omp( mythread+1, counter, 1 ) = .false.
          !$omp flush( locks_omp )
          !$  DO WHILE( ANY( locks_omp( :, counter, 1 ) ) )
          !$omp flush( locks_omp )
          !$  END DO

          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) CALL SYSTEM_CLOCK( time(4) )
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) tfft%time_adding( 18 ) = tfft%time_adding( 18 ) + ( time(4) - time(3) )

          !$  IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) THEN
          !$     locks_calc_inv( parai%node_me+1, counter ) = .false.
          !$omp flush( locks_calc_inv )
          !$  END IF
  
       ELSE IF( step .eq. 2 ) THEN
  
          IF( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) CALL SYSTEM_CLOCK( time(5) )

          !$omp flush( locks_calc_inv )
          !$  DO WHILE( ANY( locks_calc_inv( :, counter ) ) )
          !$omp flush( locks_calc_inv )
          !$  END DO

          IF( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) CALL SYSTEM_CLOCK( time(6) )
          IF( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) tfft%time_adding( 19 ) = tfft%time_adding( 19 ) + ( time(6) - time(5) )
  
          CALL fft_com( tfft, remswitch, work_buffer, 1 )

          IF( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) CALL SYSTEM_CLOCK( time(7) )
          IF( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) tfft%time_adding( 20 ) = tfft%time_adding( 20 ) + ( time(7) - time(6) )
  
          !$  locks_com_inv( parai%node_me+1, counter ) = .false.
          !$omp flush( locks_com_inv )
  
       ELSE IF( step .eq. 3 ) THEN
  
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) CALL SYSTEM_CLOCK( time(8) )

          !$omp flush( locks_com_inv )
          !$  DO WHILE( ANY( locks_com_inv( :, counter ) ) .and. parai%nnode .ne. 1 )
          !$omp flush( locks_com_inv )
          !$  END DO

          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) CALL SYSTEM_CLOCK( time(9) )
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) tfft%time_adding( 21 ) = tfft%time_adding( 21 ) + ( time(9) - time(8) )

          CALL invfft_y_section( tfft, f_inout1, f_inout2(:,work_buffer), f_inout3, &
                                 tfft%map_acinv_wave, tfft%map_acinv_wave_rem, counter, remswitch, mythread, tfft%nr1w, ispec )

       ELSE IF( step .eq. 4 ) THEN

          CALL invfft_x_section( tfft, f_inout1, remswitch, mythread )

          !Wouldnt all threads need to be finished before this lock can be lifted?
!          IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) THEN
!             IF( tfft%which_wave .eq. 1 ) THEN
!                !$  locks_calc_1( parai%node_me+1, 1+current+(fft_batchsize*fft_buffsize): &
!                !$                             ispec+current+(fft_batchsize*fft_buffsize) ) = .false.
!                !$omp flush( locks_calc_1 )
!             ELSE
!                !$  locks_calc_2( parai%node_me+1, 1+current:ispec+current ) = .false.
!                !$omp flush( locks_calc_2 )
!             END IF
!          END IF

          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) CALL SYSTEM_CLOCK( time(10) )

          !$  locks_omp_big( mythread+1, ispec, counter, 6 ) = .false.
          !$omp flush( locks_omp_big )
          !$  DO WHILE( ANY( locks_omp_big( :, ispec, counter, 6 ) ) )
          !$omp flush( locks_omp_big )
          !$  END DO

          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) CALL SYSTEM_CLOCK( time(11) )
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) tfft%time_adding( 22 ) = tfft%time_adding( 22 ) + ( time(11) - time(10) )

       END IF
  
    ELSE !! fwfft
  
       IF( step .eq. 1 ) THEN

!          !$OMP Barrier !Should be here?
          CALL fwfft_x_section( tfft, f_inout1(:,1), f_inout2, ispec, counter, remswitch, mythread, tfft%nr1w )

       ELSE IF( step .eq. 2 ) THEN

          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) CALL SYSTEM_CLOCK( time(12) )

          !$omp flush( locks_calc_2 )
          !$  DO WHILE( ANY(locks_calc_2(:,1+current:ispec+current ) ) )
          !$omp flush( locks_calc_2 )
          !$  END DO
  
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) CALL SYSTEM_CLOCK( time(13) )
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) tfft%time_adding( 23 ) = tfft%time_adding( 23 ) + ( time(13) - time(12) )
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) &
              tfft%time_adding_extra( 23, counter ) = tfft%time_adding_extra( 23, counter ) + ( time(13) - time(12) )
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) &
              tfft%time_adding_extra2( 23, ispec, counter ) = tfft%time_adding_extra2( 23, ispec, counter ) + ( time(13) - time(12) )
  
          CALL fwfft_y_section( tfft, f_inout1, f_inout2(:,work_buffer), f_inout3(:,work_buffer), &
                                    tfft%map_pcfw(:,1), batch_size, ispec, counter, remswitch, mythread, tfft%nr1w, tfft%nsw )

          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) CALL SYSTEM_CLOCK( time(14) )

          !$  locks_omp_big( mythread+1, ispec, counter, 4 ) = .false.
          !$omp flush( locks_omp_big )
          !$  DO WHILE( ANY( locks_omp_big( :, ispec, counter, 4 ) ) )
          !$omp flush( locks_omp_big )
          !$  END DO
  
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) CALL SYSTEM_CLOCK( time(15) )
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) tfft%time_adding( 24 ) = tfft%time_adding( 24 ) + ( time(15) - time(14) )
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) &
              tfft%time_adding_extra( 24, counter ) = tfft%time_adding_extra( 24, counter ) + ( time(15) - time(14) )
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) &
              tfft%time_adding_extra2( 24, ispec, counter ) = tfft%time_adding_extra2( 24, ispec, counter ) + ( time(15) - time(14) )

!          !$  IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) THEN
!          !$     locks_calc_fw( parai%node_me+1, counter ) = .false.
!          !$omp flush( locks_calc_fw )
!          !$  END IF

       ELSE IF( step .eq. 3 ) THEN

          IF( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) CALL SYSTEM_CLOCK( time(16) )

          !$omp flush( locks_calc_fw )
          !$  DO WHILE( ANY( locks_calc_fw( :, counter ) ) )
          !$omp flush( locks_calc_fw )
          !$  END DO

          IF( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) CALL SYSTEM_CLOCK( time(21) )
          IF( ( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) .and. remswitch .eq. 1 ) tfft%time_adding( 25 ) = tfft%time_adding( 25 ) + ( time(21) - time(16) )
          IF( ( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) .and. remswitch .eq. 1 ) &
              tfft%time_adding_extra( 25, counter ) = ( time(21) - time(16) )

          !$ locks_all_com2( parai%node_me+1, counter ) = .false.
          !$omp flush( locks_all_com2 )
          !$  DO WHILE( ANY( locks_all_com2( :, counter ) ) )
          !$omp flush( locks_all_com2 )
          !$  END DO
 
          CALL MPI_BARRIER(parai%allgrp, ierr)
 
          IF( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) CALL SYSTEM_CLOCK( time(17) )
          IF( ( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) .and. remswitch .eq. 1 ) tfft%time_adding( 33 ) = tfft%time_adding( 33 ) + ( time(17) - time(21) )
          IF( ( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) .and. remswitch .eq. 1 ) &
              tfft%time_adding_extra( 33, counter ) = ( time(17) - time(21) )
     
          CALL fft_com( tfft, remswitch, work_buffer, 1 )

          IF( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) CALL SYSTEM_CLOCK( time(18) )
          IF( ( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) .and. remswitch .eq. 1 ) tfft%time_adding( 26 ) = tfft%time_adding( 26 ) + ( time(18) - time(17) )
          IF( ( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) .and. remswitch .eq. 1 ) &
              tfft%time_adding_extra( 26, counter ) = ( time(18) - time(17) )
          IF( ( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) .and. remswitch .eq. 1 ) tfft%com_count = tfft%com_count + 1

          !$ locks_all_com( parai%node_me+1, counter ) = .false.
          !$omp flush( locks_all_com )
          !$  DO WHILE( ANY( locks_all_com( :, counter ) ) )
          !$omp flush( locks_all_com )
          !$  END DO

          CALL MPI_BARRIER(parai%allgrp, ierr)

          IF( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) CALL SYSTEM_CLOCK( time(19) )
          IF( ( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) .and. remswitch .eq. 1 ) tfft%time_adding( 31 ) = tfft%time_adding( 31 ) + ( time(19) - time(18) )
          IF( ( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) .and. remswitch .eq. 1 ) &
              tfft%time_adding_extra( 31, counter ) = ( time(19) - time(18) )
          IF( ( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) .and. remswitch .eq. 1 ) tfft%time_adding( 32 ) = tfft%time_adding( 32 ) + ( time(19) - time(17) )
          IF( ( mythread .eq. 0 .or. parai%ncpus_FFT .eq. 1 ) .and. remswitch .eq. 1 ) &
              tfft%time_adding_extra( 32, counter ) = ( time(19) - time(17) )

          !$  locks_com_fw( parai%node_me+1, counter ) = .false.
          !$omp flush( locks_com_fw )
  
       ELSE IF( step .eq. 4 ) THEN
  
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) CALL SYSTEM_CLOCK( time(19) )

          !$omp flush( locks_com_fw )
          !$  DO WHILE( ANY( locks_com_fw( :, counter ) ) .and. parai%nnode .ne. 1 )
          !$omp flush( locks_com_fw )
          !$  END DO

          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) CALL SYSTEM_CLOCK( time(20) )
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) tfft%time_adding( 27 ) = tfft%time_adding( 27 ) + ( time(20) - time(19) )
          IF( parai%ncpus_FFT .eq. 1 .or. ( mythread .eq. 1 .and. cntl%overlapp_comm_comp ) .or. ( mythread .eq. 0 .and. .not. cntl%overlapp_comm_comp ) ) &
              tfft%time_adding_extra( 27, counter ) = ( time(20) - time(19) )
  
          CALL fwfft_z_section( tfft, f_inout1(:,work_buffer), f_inout2, counter, batch_size, remswitch, mythread, tfft%nsw )

!          !$  locks_omp( mythread+1, counter, 11 ) = .false.
!          !$omp flush( locks_omp )
!          !$  DO WHILE( ANY( locks_omp( :, counter, 11 ) ) )
!          !$omp flush( locks_omp )
!          !$  END DO
        
          IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) THEN
             IF( cntl%krwfn ) THEN
             !$   locks_calc_2( parai%node_me+1, 1+(counter+fft_buffsize-1)*fft_batchsize:batch_size+(counter+fft_buffsize-1)*fft_batchsize ) = .false.
             !$omp flush( locks_calc_2 )
             ELSE 
             !$   locks_calc_1( parai%node_me+1, 1+(counter+fft_buffsize-1)*fft_batchsize:fft_batchsize+(counter+fft_buffsize-1)*fft_batchsize ) = .false.
             !$omp flush( locks_calc_1 )
             END IF
          END IF
 
       END IF
  
    END IF

    IF( cntl%fft_tune_batchsize ) THEN
       IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
    ELSE
       IF( parai%ncpus_FFT .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
    END IF
  
  END SUBROUTINE fft_improved_batch

  SUBROUTINE fft_improved( isign, tfft, f, ngs, nr1s, ir1s, nss )

    IMPLICIT NONE
  
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft 
    INTEGER, INTENT(IN) :: isign, ngs
    COMPLEX(DP), TARGET, INTENT(INOUT) :: f(:)
    INTEGER, INTENT(IN) :: ir1s(:), nss(:), nr1s

    COMPLEX(DP), POINTER, CONTIGUOUS :: f_2d(:,:)

#ifdef _USE_SCRATCHLIBRARY
    COMPLEX(DP), POINTER, SAVE __CONTIGUOUS, ASYNCHRONOUS :: aux(:,:)
    INTEGER(int_8) :: il_aux(2)
#else
    COMPLEX(DP), ALLOCATABLE, SAVE, TARGET, ASYNCHRONOUS  :: aux(:,:)
#endif

    INTEGER :: i, ierr, isub, mythread
    CHARACTER(*), PARAMETER                  :: procedureN = 'fftpw'
    LOGICAL, SAVE :: first = .true.
    INTEGER, SAVE :: sendsize
    TYPE(C_PTR) :: baseptr( 0:parai%node_nproc-1 )
    INTEGER :: arrayshape(3)
    COMPLEX(DP), SAVE, POINTER, CONTIGUOUS   :: Big_Pointer(:,:,:)

    CALL tiset(procedureN,isub)

#ifdef _USE_SCRATCHLIBRARY
    il_aux(1) = tfft%nr2 * nr1s * tfft%my_nr3p  
    il_aux(2) = 1
    CALL request_scratch(il_aux,aux,procedureN//'_aux',ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate aux', &
         __LINE__,__FILE__)
#else
    IF( .not. allocated( aux ) ) ALLOCATE( aux( tfft%nr2 * nr1s * tfft%my_nr3p, 1 ) )
    IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate aux', &
         __LINE__,__FILE__)
#endif

    tfft%which = 2

    IF( first ) THEN
  
       first = .false.

       sendsize = MAXVAL ( tfft%nr3p ) * MAXVAL( nss ) * parai%max_node_nproc * parai%max_node_nproc
   
       CALL mp_win_alloc_shared_mem( 'c', sendsize*parai%nnode*2, 1, baseptr, parai%node_nproc, parai%node_me, parai%node_grp )

       arrayshape(1) = sendsize*parai%nnode
       arrayshape(2) = 1
       arrayshape(3) = 2
       CALL C_F_POINTER( baseptr(0), Big_Pointer, arrayshape )
       comm_send => Big_Pointer(:,:,1) 
       comm_recv => Big_Pointer(:,:,2) 

       CALL Prep_fft_com( comm_send, comm_recv, sendsize, 0, parai%nnode, parai%me, parai%my_node, parai%node_me, &
                          parai%node_nproc, parai%max_node_nproc, parai%cp_overview, 1, tfft%comm_sendrecv(:,2), tfft%do_comm(2), 2 )

       CALL Make_Manual_Maps( tfft, 1, 0, nss, nr1s, ngs, tfft%which, 0 )

       IF( .not. allocated( locks_omp ) ) ALLOCATE( locks_omp( parai%ncpus_FFT, 1, 20 ) )
       !$ locks_omp = .true.
       IF( .not. allocated( locks_omp_big ) ) ALLOCATE( locks_omp_big( parai%ncpus_FFT, 1, 1, 20 ) )
       !$ locks_omp_big = .true.

    END IF

!    comm_send => Big_Pointer(:,:,1) 
!    comm_recv => Big_Pointer(:,:,2) 
    CALL MPI_BARRIER( parai%allgrp, ierr )
    !$ locks_omp = .true.
    !$ locks_omp_big = .true.

    !$OMP parallel num_threads( parai%ncpus_FFT ) &
    !$omp private(mythread) &
    !$omp proc_bind(close)
    !$ mythread = omp_get_thread_num()

    IF( isign .eq. -1 ) THEN !!  invfft


       CALL invfft_z_section( tfft, f, comm_send(:,1), comm_recv(:,1), 1, 1, mythread, nss, 1 )

       !$OMP barrier
       !$OMP master
          CALL MPI_BARRIER( parai%allgrp, ierr )
          IF( tfft%do_comm(2) ) CALL fft_com( tfft, 1, 1, 2 )
          CALL MPI_BARRIER( parai%allgrp, ierr )
       !$OMP end master
       !$OMP barrier

       CALL invfft_y_section( tfft, f, comm_recv(:,1), aux, tfft%map_acinv_pot, tfft%map_acinv_pot, 1, 1, mythread, tfft%nr1p, 1 )

       CALL invfft_x_section( tfft, f, 1, mythread )

    ELSE !! fw fft

       CALL fwfft_x_section( tfft, f, aux, 1, 1, 1, mythread, tfft%nr1p )

       CALL fwfft_y_section( tfft, aux, comm_send(:,1), comm_recv(:,1), tfft%map_pcfw(:,2), 1, 1, 1, 1, mythread, tfft%nr1p, tfft%nsp )
    
       !$OMP barrier
       !$OMP master
          CALL MPI_BARRIER( parai%allgrp, ierr )
          IF( tfft%do_comm(2) ) CALL fft_com( tfft, 1, 1, 2 )
          CALL MPI_BARRIER( parai%allgrp, ierr )
       !$OMP end master
       !$OMP barrier
    
       CALL fwfft_z_section( tfft, comm_recv(:,1), f, 1, 1, 1, mythread, tfft%nsp, tfft%tscale )
    
    END IF

    !$omp end parallel
 
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_aux,aux,procedureN//'_aux',ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate aux', &
         __LINE__,__FILE__)
#else
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate aux', &
         __LINE__,__FILE__)
#endif

    tfft%which = 1

    CALL tihalt(procedureN,isub)
 
  END SUBROUTINE fft_improved

END MODULE fftmain_utils
