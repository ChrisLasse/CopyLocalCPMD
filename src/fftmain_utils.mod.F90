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
       xf, yf, fft_residual, fft_total, fft_numbatches, fft_batchsize, locks_inv, locks_fw
  USE fft_maxfft,                      ONLY: maxfftn, maxfft
  USE fftcu_methods,                   ONLY: fftcu_frw_full_1,&
                                             fftcu_frw_full_2,&
                                             fftcu_frw_sprs_1,&
                                             fftcu_frw_sprs_2,&
                                             fftcu_inv_full_1,&
                                             fftcu_inv_full_2,&
                                             fftcu_inv_sprs_1,&
                                             fftcu_inv_sprs_2
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
                                             unpack_y2x_n
  USE fftpw_converting,                ONLY: ConvertFFT_Coeffs
  USE fftpw_make_maps,                 ONLY: Prep_copy_Maps,&
                                             Set_Req_Vals
  USE fftpw_param
  USE fftpw_types,                     ONLY: PW_fft_type_descriptor,&
                                             create_shared_memory_window_1d
  USE fftpw_batching
  USE fftpw_batchingManual
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE machine,                         ONLY: m_walltime
  USE mltfft_utils,                    ONLY: mltfft_cuda,&
                                             mltfft_default,&
                                             mltfft_essl,&
                                             mltfft_fftw,&
                                             mltfft_hp
  USE parac,                           ONLY: parai
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw
  USE thread_view_types,               ONLY: thread_view_t
  USE timer,                           ONLY: tihalt,&
                                             tiset

  !$ USE omp_lib, ONLY: omp_in_parallel, omp_get_thread_num, &
  !$ omp_set_num_threads, omp_set_nested

  USE, INTRINSIC :: ISO_C_BINDING,     ONLY: C_PTR, C_NULL_PTR
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mltfft
  PUBLIC :: invfftn
  PUBLIC :: fwfftn
  PUBLIC :: invfftn_batch
  PUBLIC :: invfftn_batch_com
  PUBLIC :: fwfftn_batch_com
  PUBLIC :: fwfftn_batch


  COMPLEX(DP), POINTER, SAVE, CONTIGUOUS :: comm_send(:,:)
  COMPLEX(DP), POINTER, SAVE, CONTIGUOUS :: comm_recv(:,:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_calc_inv(:,:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_com_inv(:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_calc_fw(:,:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_com_fw(:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_calc_1(:,:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_calc_2(:,:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_sing_1(:,:)
  LOGICAL, POINTER, SAVE, CONTIGUOUS :: locks_sing_2(:,:)

  PUBLIC :: comm_send
  PUBLIC :: comm_recv
  PUBLIC :: locks_calc_inv
  PUBLIC :: locks_com_inv
  PUBLIC :: locks_calc_fw
  PUBLIC :: locks_com_fw
  PUBLIC :: locks_calc_1
  PUBLIC :: locks_calc_2
  PUBLIC :: locks_sing_1
  PUBLIC :: locks_sing_2

  PUBLIC :: invfft_pwbatch
  PUBLIC :: fwfft_pwbatch
  PUBLIC :: invfftu
  PUBLIC :: fwfftu
  PUBLIC :: invfft_4S
  PUBLIC :: fwfft_4S
  !public :: fftnew


CONTAINS


  ! ==================================================================
  SUBROUTINE fftnew(isign,f,sparse,comm)
    ! ==--------------------------------------------------------------==
#ifdef __PARALLEL
    USE mpi_f08
#endif
    INTEGER, INTENT(IN)                      :: isign
    COMPLEX(real_8)                          :: f(:)
    LOGICAL, INTENT(IN)                      :: sparse
#ifdef __PARALLEL
    type(MPI_COMM), INTENT(IN)               :: comm
#else
    INTEGER, INTENT(IN)                      :: comm
#endif
    CHARACTER(*), PARAMETER                  :: procedureN = 'fftnew'

    COMPLEX(real_8), DIMENSION(:), &
      POINTER __CONTIGUOUS                   :: xf_ptr, yf_ptr
    INTEGER                                  :: lda, m, mm, n1o, n1u, ierr
    INTEGER(int_8)                           :: il_xf(2)
    REAL(real_8)                             :: scale

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
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_xf,yf,procedureN//'_yf',ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL free_scratch(il_xf,xf,procedureN//'_xf',ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#endif
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
    COMPLEX(real_8),INTENT(INOUT)            :: f(:)
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
    COMPLEX(real_8)                          :: a(:)
    INTEGER                                  :: ldax, lday
    COMPLEX(real_8)                          :: b(:)
    INTEGER                                  :: ldbx, ldby, n, m, isign
    REAL(real_8)                             :: scale
    TYPE(thread_view_t), INTENT(IN), &
      OPTIONAL                               :: thread_view

    CHARACTER(*), PARAMETER                  :: procedureN = 'mltfft'

    INTEGER                                  :: device_idx, host_buff_ptr, &
                                                stream_idx
    TYPE(cp_cufft_plans_t), POINTER          :: plans_d
    TYPE(cublas_handle_t), POINTER           :: blas_handle_p
    TYPE(cuda_memory_t), POINTER             :: t1_d, t2_d
    TYPE(cuda_stream_t), POINTER             :: stream_p
    TYPE(cufft_plan_t), POINTER              :: plan_p

! ==--------------------------------------------------------------==
!vw for the moment we use CPU FFT 
!IF( cp_cuda_env%use_fft ) THEN

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
#elif defined(__HAS_FFT_FFTW3)
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
    USE fftpw_base,                             ONLY: dfftp
#ifdef __PARALLEL
    USE mpi_f08
#endif
    COMPLEX(real_8)                          :: f(:)
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
          CALL fftpw( isign, dfftp, f, dfftp%ngm, dfftp%nr1p, dfftp%ir1p, dfftp%nsp )
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
    USE fftpw_base,                             ONLY: dfftp
#ifdef __PARALLEL
    USE mpi_f08
#endif
    COMPLEX(real_8)                          :: f(:)
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
          CALL fftpw( isign, dfftp, f, dfftp%ngm, dfftp%nr1p, dfftp%ir1p, dfftp%nsp )
       END IF
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fwfftn
  ! ==================================================================
  ! ==================================================================
  ! NEW-er FFT CODE -> includes PW things
  ! ==================================================================
  SUBROUTINE invfftu(f,sparse,comm,thread_view, copy_data_to_device, copy_data_to_host )
    ! ==--------------------------------------------------------------==
    ! == COMPUTES THE INVERSE FOURIER TRANSFORM OF A COMPLEX          ==
    ! == FUNCTION F. THE FOURIER TRANSFORM IS                         ==
    ! == RETURNED IN F (THE INPUT F IS OVERWRITTEN).                  ==
    ! ==--------------------------------------------------------------==
    USE fftpw_base,                             ONLY: dfftp
#ifdef __PARALLEL
    USE mpi_f08
#endif
    COMPLEX(real_8)                          :: f(:)
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
          CALL fftpw( isign, dfftp, f, dfftp%ngm, dfftp%nr1p, dfftp%ir1p, dfftp%nsp )
          write(6,*) "NEW INVFFT CALLED"
       END IF
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE invfftu
  ! ==================================================================
  SUBROUTINE fwfftu(f,sparse,comm,thread_view, copy_data_to_device, copy_data_to_host )
    ! ==--------------------------------------------------------------==
    ! == COMPUTES THE FORWARD FOURIER TRANSFORM OF A COMPLEX          ==
    ! == FUNCTION F. THE FOURIER TRANSFORM IS                         ==
    ! == RETURNED IN F IN OUTPUT (THE INPUT F IS OVERWRITTEN).        ==
    ! ==--------------------------------------------------------------==
    USE fftpw_base,                             ONLY: dfftp
#ifdef __PARALLEL
    USE mpi_f08
#endif
    COMPLEX(real_8)                          :: f(:)
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

    CHARACTER(*), PARAMETER                  :: procedureN = 'fwfftu'

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
          CALL fftpw( isign, dfftp, f, dfftp%ngm, dfftp%nr1p, dfftp%ir1p, dfftp%nsp )
          write(6,*) "NEW FWFFT CALLED"
       END IF
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fwfftu
  ! ==================================================================

  SUBROUTINE fwfftn_batch(f,n,swap,step,ibatch)
    ! ==--------------------------------------------------------------==
    ! == COMPUTES THE FORWARD FOURIER TRANSFORM OF A COMPLEX          ==
    ! == FUNCTION F_IN. THE FOURIER TRANSFORM IS                      ==
    ! == RETURNED IN F_OUT                                            ==
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: n, swap, step, ibatch
    COMPLEX(real_8), INTENT(INOUT)           :: f(fpar%kr1*fpar%kr2s*fpar%kr3s*fft_batchsize)
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
  SUBROUTINE invfftn_batch(f,n,swap,step,ibatch)
    ! ==--------------------------------------------------------------==
    ! == COMPUTES THE INVERSE FOURIER TRANSFORM OF A COMPLEX          ==
    ! == FUNCTION F_IN. THE FOURIER TRANSFORM IS                      ==
    ! == RETURNED IN F_OUT                                            ==
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: n, swap, step, ibatch
    COMPLEX(real_8), INTENT(INOUT)           :: f(fpar%kr1*fpar%kr2s*fpar%kr3s*fft_batchsize)
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
  ! ==================================================================
  SUBROUTINE invfftn_batch_com(int_mod)
    ! ==--------------------------------------------------------------==
    ! == COMPUTES THE INVERSE FOURIER TRANSFORM OF A COMPLEX          ==
    ! == FUNCTION F_IN. THE FOURIER TRANSFORM IS                      ==
    ! == RETURNED IN F_OUT                                            ==
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: int_mod
    CHARACTER(*), PARAMETER                  :: procedureN = 'invfftn_batch_com'

    INTEGER                                  :: ibatch,n,lda,isub, swap

    CALL tiset(procedureN,isub)
    DO ibatch=1,fft_numbatches+1
       IF(ibatch.LE.fft_numbatches)THEN
          n=fft_batchsize
       ELSE
          n=fft_residual
       END IF
       IF(n.NE.0)THEN
          lda=lsrm*lr1m*n
          swap=mod(ibatch,int_mod)+1
          !$omp flush(locks_inv)
          !$ DO WHILE ( locks_inv(ibatch,1) )
          !$omp flush(locks_inv)
          !$ END DO
          !$ locks_inv(ibatch,1) = .TRUE.
          !$omp flush(locks_inv)
          CALL fft_comm(yf(:,swap),xf(:,swap),lda,cntl%tr4a2a,parai%allgrp)
          !$ locks_inv(ibatch,1)=.FALSE.
          !$omp flush(locks_inv)
          !$ locks_inv(ibatch,2)=.FALSE.
          !$omp flush(locks_inv)
       END IF
    END DO

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE invfftn_batch_com
  ! ==================================================================
  SUBROUTINE fwfftn_batch_com(int_mod)
    ! ==--------------------------------------------------------------==
    ! == COMPUTES THE INVERSE FOURIER TRANSFORM OF A COMPLEX          ==
    ! == FUNCTION F_IN. THE FOURIER TRANSFORM IS                      ==
    ! == RETURNED IN F_OUT                                            ==
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: int_mod
    CHARACTER(*), PARAMETER                  :: procedureN = 'fwfftn_batch_com'

    INTEGER                                  :: ibatch,n,lda,isub, swap, isub1

    CALL tiset(procedureN,isub)
    DO ibatch=1,fft_numbatches+1
       IF(ibatch.LE.fft_numbatches)THEN
          n=fft_batchsize
       ELSE
          n=fft_residual
       END IF
       IF(n.NE.0)THEN
          lda=lsrm*lr1m*n
          swap=mod(ibatch,int_mod)+1
          !$omp flush(locks_fw)
          !$ DO WHILE ( locks_fw(ibatch,1) )
          !$omp flush(locks_fw)
          !$ END DO
          !$ locks_fw(ibatch,1) = .TRUE.
          !$omp flush(locks_fw)
          CALL fft_comm(xf(:,swap),yf(:,swap),lda,cntl%tr4a2a,parai%allgrp)
          !$ locks_fw(ibatch,1)=.FALSE.
          !$omp flush(locks_fw)
          !$ locks_fw(ibatch,2)=.FALSE.
          !$omp flush(locks_fw)
       END IF
    END DO
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fwfftn_batch_com

  SUBROUTINE invfft_pwbatch( dfft, step, batch_size, set_size_1, set_size_2, set_size_3, counter, work_buffer, f_in, f_inout1, f_inout2 )
    IMPLICIT NONE
  
    TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
    INTEGER, INTENT(IN) :: step, batch_size, counter, work_buffer
    INTEGER, INTENT(IN) :: set_size_1, set_size_2, set_size_3
    COMPLEX(DP), INTENT(INOUT) :: f_in(:,:)
    COMPLEX(DP), INTENT(INOUT) :: f_inout1(:,:)
  
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout2(:)
  
    IF( step .eq. 1 ) THEN
       CALL fftpw_batch( dfft, -1, step, batch_size, set_size_1, set_size_2, set_size_3, counter, work_buffer, f_in, f_inout1, f_inout2 )
    ELSE IF( step .eq. 2 ) THEN                                          
       CALL fftpw_batch( dfft, -1, step, batch_size, set_size_1, set_size_2, set_size_3, counter, work_buffer, f_in, f_inout1 )
    ELSE IF( step .eq. 3 ) THEN                                       
       CALL fftpw_batch( dfft, -1, step, batch_size, set_size_1, set_size_2, set_size_3, counter, work_buffer, f_in, f_inout1 )
    END IF
  
  END SUBROUTINE invfft_pwbatch
  
  SUBROUTINE fwfft_pwbatch( dfft, step, batch_size, set_size, counter, work_buffer, f_in, f_inout1, f_inout2 )
    IMPLICIT NONE
  
    TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
    INTEGER, INTENT(IN) :: step, batch_size, counter, work_buffer, set_size
    COMPLEX(DP), INTENT(INOUT) :: f_in(:,:)
    COMPLEX(DP), INTENT(INOUT) :: f_inout1(:,:)
  
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout2(:)
  
    IF( step .eq. 1 ) THEN
       CALL fftpw_batch( dfft, 1, step, batch_size, 1, 1, 1, counter, work_buffer, f_in, f_inout1, f_inout2 )
    ELSE IF( step .eq. 2 ) THEN
       CALL fftpw_batch( dfft, 1, step, batch_size, 1, 1, 1, counter, work_buffer, f_in, f_inout1 )
    ELSE IF( step .eq. 3 ) THEN
       CALL fftpw_batch( dfft, 1, step, batch_size, set_size, 1, 1, counter, work_buffer, f_in, f_inout1, f_inout2 )
    END IF
  
  END SUBROUTINE fwfft_pwbatch

  SUBROUTINE invfft_4S( dfft, step, batch_size, divparam_1, divparam_2, divparam_3, counter, work_buffer, first_dim, priv, f_inout1, f_inout2, f_inout3, f_inout4, ip, jp )
    IMPLICIT NONE
  
    TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
    INTEGER, INTENT(IN) :: step, batch_size, counter, work_buffer, first_dim
    INTEGER, INTENT(IN) :: divparam_1, divparam_2, divparam_3
    INTEGER, INTENT(INOUT) :: priv(:)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout1(first_dim,*)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout2(:,:)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout3(:,:)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout4(:,:)
    INTEGER,     OPTIONAL, INTENT(INOUT) :: ip, jp
  
    IF( step .eq. 1 ) THEN
       CALL fftpw_4S( dfft, -1, step, batch_size, divparam_1, divparam_2, divparam_3, counter, work_buffer, first_dim, priv, &
                      f_inout1=f_inout1, f_inout2=f_inout2, f_inout3=f_inout3, f_inout4=f_inout4, ip=ip, jp=jp )
    ELSE IF( step .eq. 2 ) THEN                                          
       CALL fftpw_4S( dfft, -1, step, batch_size, 0, 0, 0, counter, work_buffer, 0, priv, f_inout2=f_inout2, f_inout3=f_inout3 )
    ELSE IF( step .eq. 3 ) THEN                                       
       CALL fftpw_4S( dfft, -1, step, batch_size, divparam_1, divparam_2, divparam_3, counter, work_buffer, first_dim, priv, &
                      f_inout1=f_inout1, f_inout2=f_inout2, f_inout3=f_inout3, f_inout4=f_inout4 )
    ELSE IF( step .eq. 4 ) THEN                                       
       CALL fftpw_4S( dfft, -1, step, batch_size, divparam_1, divparam_2, divparam_3, counter, work_buffer, first_dim, priv, &
                      f_inout1=f_inout1, f_inout2=f_inout2, f_inout3=f_inout3, f_inout4=f_inout4 )
    END IF
  
  END SUBROUTINE invfft_4S
  
  SUBROUTINE fwfft_4S( dfft, step, batch_size, divparam_1, divparam_2, divparam_3, counter, work_buffer, first_dim, priv, f_inout1, f_inout2, f_inout3, f_inout4 )
    IMPLICIT NONE
  
    TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
    INTEGER, INTENT(IN) :: step, batch_size, counter, work_buffer, divparam_1, divparam_2, divparam_3, first_dim
    INTEGER, INTENT(INOUT) :: priv(:)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout1(first_dim,*)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout2(:,:)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout3(:,:)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout4(:,:)
  
    IF( step .eq. 1 ) THEN
       CALL fftpw_4S( dfft, 1, step, batch_size, divparam_1, divparam_2, divparam_3, counter, work_buffer, first_dim, priv, &
                      f_inout1=f_inout1, f_inout2=f_inout2, f_inout3=f_inout3, f_inout4=f_inout4 )
    ELSE IF( step .eq. 2 ) THEN
       CALL fftpw_4S( dfft, 1, step, batch_size, divparam_1, divparam_2, divparam_3, counter, work_buffer, first_dim, priv, &
                      f_inout1=f_inout1, f_inout2=f_inout2, f_inout3=f_inout3, f_inout4=f_inout4 )
    ELSE IF( step .eq. 3 ) THEN
       CALL fftpw_4S( dfft, 1, step, batch_size, 0, 0, 0, counter, work_buffer, 0, priv, f_inout2=f_inout2, f_inout3=f_inout3 )
    ELSE IF( step .eq. 4 ) THEN
       CALL fftpw_4S( dfft, 1, step, batch_size, divparam_1, divparam_2, divparam_3, counter, work_buffer, first_dim, priv, &
                      f_inout1=f_inout1, f_inout2=f_inout2, f_inout3=f_inout3, f_inout4=f_inout4 )
    END IF
  
  END SUBROUTINE fwfft_4S

  SUBROUTINE fftpw_batch( dfft, isign, step, batch_size, set_size_1, set_size_2, set_size_3, counter, work_buffer, f_in, f_inout1, f_inout2 )
  USE iso_fortran_env
    IMPLICIT NONE
  
    TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
    INTEGER, INTENT(IN) :: isign, step, batch_size, counter, work_buffer
    INTEGER, INTENT(IN) :: set_size_1, set_size_2, set_size_3
    COMPLEX(DP), INTENT(INOUT) :: f_in(:,:)
    COMPLEX(DP), INTENT(INOUT) :: f_inout1(:,:)
  
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout2(:)
  
    INTEGER :: ibatch, iset, current, howmany, group_size
    LOGICAL :: last
    INTEGER(INT64) :: auto_time(4)
    INTEGER(INT64) :: time(17), cr
 
    current = (counter-1)*dfft%batch_size_save
  
    IF( isign .eq. -1 ) THEN !!  invfft
  
       IF( step .eq. 1 ) THEN

          CALL SYSTEM_CLOCK( time(1) )
  
          !In theory, locks could be made faster by checking each iset individually for non-remainder cases 
          !$omp flush( locks_calc_1 )
          !$  DO WHILE( ANY(locks_calc_1( :, 1+current:dfft%batch_size_save+current ) ) )
          !$omp flush( locks_calc_1 )
          !$  END DO

          CALL SYSTEM_CLOCK( time(2) )
          dfft%time_adding( 22 ) = dfft%time_adding( 22 ) + ( time(2) - time(1) )

          CALL SYSTEM_CLOCK( auto_time(1) )
          DO iset = 1, set_size_1

             !Too slow? lets wait and see!
             IF( iset .ne. set_size_1 .or. iset .eq. 1 .or. group_size * set_size_1 .eq. batch_size ) THEN
                IF( mod( batch_size, set_size_1 ) .eq. 0 ) THEN
                   group_size = batch_size / set_size_1
                ELSE
                   group_size = batch_size / set_size_1 + 1
                ENDIF
                dfft%z_group_size_save = group_size
             ELSE
                group_size = batch_size - (set_size_1-1) * group_size
             END IF
 
!             IF( counter .eq. dfft%max_nbnd .and. iset .eq. set_size_1 .and. dfft%uneven ) THEN
!                last = .true.
!             ELSE
!                last = .false.
!             END IF

             CALL SYSTEM_CLOCK( time(2) )
     
             CALL Prepare_Psi_overlapp( dfft, f_in(:,1+(((iset-1)*dfft%z_group_size_save)+current)*2:2*group_size+(((iset-1)*dfft%z_group_size_save)+current)*2), &
                                        dfft%aux2( 1 : dfft%nr3 * dfft%nsw(dfft%mype+1) * group_size ), dfft%ngms, group_size, dfft%nsw, last )

             CALL SYSTEM_CLOCK( time(3) )
             dfft%time_adding( 1 ) = dfft%time_adding( 1 ) + ( time(3) - time(2) )

             CALL invfft_pre_com( dfft, dfft%aux2( 1 : dfft%nr3 * dfft%nsw(dfft%mype+1) * group_size ), f_inout1(:,work_buffer), f_inout2, iset, batch_size, group_size, dfft%nsw )

          ENDDO 
          CALL SYSTEM_CLOCK( auto_time(2) )
          dfft%auto_timings(1) = dfft%auto_timings(1) + ( auto_time(2) - auto_time(1) )
          dfft%auto_4Stimings(1) = dfft%auto_4Stimings(1) + ( auto_time(2) - auto_time(1) )
 
          !$  locks_calc_inv( dfft%my_node_rank+1, counter ) = .false.
          !$omp flush( locks_calc_inv )
  
       ELSE IF( step .eq. 2 ) THEN
  
          CALL SYSTEM_CLOCK( time(4) )

          !$omp flush( locks_calc_inv )
          !$  DO WHILE( ANY( locks_calc_inv( :, counter ) ) )
          !$omp flush( locks_calc_inv )
          !$  END DO

          CALL SYSTEM_CLOCK( time(5) )
          dfft%time_adding( 23 ) = dfft%time_adding( 23 ) + ( time(5) - time(4) )
  
          CALL fft_com( dfft, f_in(:,work_buffer), f_inout1(:,work_buffer), dfft%sendsize, dfft%my_node_rank, &
                        dfft%inter_node_comm, dfft%nodes_numb, dfft%my_inter_node_rank, dfft%non_blocking, work_buffer )
          CALL SYSTEM_CLOCK( time(16) )
          dfft%time_adding( 28 ) = dfft%time_adding( 28 ) + ( time(16) - time(5) )
  
          !$  locks_com_inv( counter ) = .false.
          !$omp flush( locks_com_inv )
  
       ELSE IF( step .eq. 3 ) THEN
  
          CALL SYSTEM_CLOCK( time(6) )

          !$omp flush( locks_com_inv )
          !$  DO WHILE( locks_com_inv( counter ) .and. .not. dfft%single_node )
          !$omp flush( locks_com_inv )
          !$  END DO

          CALL SYSTEM_CLOCK( time(7) )
          dfft%time_adding( 24 ) = dfft%time_adding( 24 ) + ( time(7) - time(6) )

          CALL invfft_after_com( dfft, f_inout1(:,1:batch_size), f_in(:,work_buffer), &
                                 dfft%aux2( 1 : dfft%my_nr1p * dfft%my_nr3p * dfft%nr2 * batch_size ), &
                                 dfft%map_acinv, dfft%map_acinv_rem, set_size_1, set_size_2, set_size_3, batch_size, dfft%nr1w )
  
          IF( dfft%vpsi ) THEN
             !$  locks_calc_2( dfft%my_node_rank+1, 1+current:batch_size+current ) = .false.
             !$omp flush( locks_calc_2 )
          ELSE
             !$  locks_calc_1( dfft%my_node_rank+1, 1+current+(dfft%batch_size_save*dfft%buffer_size_save):batch_size+current+(dfft%batch_size_save*dfft%buffer_size_save) ) = .false.
             !$omp flush( locks_calc_1 )
          END IF
  
       END IF
  
    ELSE !! fw fft
  
       IF( step .eq. 1 ) THEN

          CALL SYSTEM_CLOCK( time(8) )
  
          !$omp flush( locks_calc_2 )
          !$  DO WHILE( ANY(locks_calc_2( :, 1+current:batch_size+current ) ) )
          !$omp flush( locks_calc_2 )
          !$  END DO

          CALL SYSTEM_CLOCK( time(9) )
          dfft%time_adding( 25 ) = dfft%time_adding( 25 ) + ( time(9) - time(8) )
  
          CALL fwfft_pre_com( dfft, f_in, dfft%aux2( 1 : dfft%my_nr1p * dfft%my_nr3p * dfft%nr2 * batch_size ), &
                              f_inout1(:,work_buffer), f_inout2, set_size_1, set_size_2, batch_size, dfft%nr1w, dfft%nsw )
  
          !$  locks_calc_fw( dfft%my_node_rank+1, counter ) = .false.
          !$omp flush( locks_calc_fw )
  
       ELSE IF( step .eq. 2 ) THEN

          CALL SYSTEM_CLOCK( time(10) )
  
          !$omp flush( locks_calc_fw )
          !$  DO WHILE( ANY( locks_calc_fw( :, counter ) ) )
          !$omp flush( locks_calc_fw )
          !$  END DO

          CALL SYSTEM_CLOCK( time(11) )
          dfft%time_adding( 26 ) = dfft%time_adding( 26 ) + ( time(11) - time(10) )
     
          CALL fft_com( dfft, f_in(:,work_buffer), f_inout1(:,work_buffer), dfft%sendsize, dfft%my_node_rank, &
                        dfft%inter_node_comm, dfft%nodes_numb, dfft%my_inter_node_rank, dfft%non_blocking, work_buffer )
          CALL SYSTEM_CLOCK( time(17) )
          dfft%time_adding( 29 ) = dfft%time_adding( 29 ) + ( time(17) - time(11) )
  
          !$  locks_com_fw( counter ) = .false.
          !$omp flush( locks_com_fw )
  
       ELSE IF( step .eq. 3 ) THEN
  
          CALL SYSTEM_CLOCK( time(12) )

          !$omp flush( locks_com_fw )
          !$  DO WHILE( locks_com_fw( counter ) .and. .not. dfft%single_node )
          !$omp flush( locks_com_fw )
          !$  END DO

          CALL SYSTEM_CLOCK( time(13) )
          dfft%time_adding( 27 ) = dfft%time_adding( 27 ) + ( time(13) - time(12) )
  
          CALL SYSTEM_CLOCK( auto_time(3) )
          DO iset = 1, set_size_1

             !Too slow? lets wait and see!
             IF( iset .ne. set_size_1 .or. iset .eq. 1 .or. group_size * set_size_1 .eq. batch_size ) THEN
                IF( mod( batch_size, set_size_1 ) .eq. 0 ) THEN
                   group_size = batch_size / set_size_1
                ELSE
                   group_size = batch_size / set_size_1 + 1
                ENDIF
                dfft%z_group_size_save = group_size
             ELSE
                group_size = batch_size - (set_size_1-1) * group_size
             END IF
 
!             IF( counter .eq. dfft%max_nbnd .and. iset .eq. set_size_1 .and. dfft%uneven ) THEN
!                last = .true.
!             ELSE
!                last = .false.
!             END IF

             CALL fwfft_after_com( dfft, f_inout2, dfft%aux2( 1 : dfft%nr3 * dfft%nsw(dfft%mype+1) * group_size ), iset, batch_size, group_size, dfft%nsw )

             CALL SYSTEM_CLOCK( time(14) )

             CALL Accumulate_Psi_overlapp( dfft, dfft%aux2( 1 : dfft%nr3 * dfft%nsw(dfft%mype+1) * group_size ), &
                                           f_inout1(:,1+(((iset-1)*dfft%z_group_size_save)+current)*2:2*group_size+(((iset-1)*dfft%z_group_size_save)+current)*2), &
                                           dfft%ngms, group_size, last, dfft%nsw, &
                                           f_in    (:,1+(((iset-1)*dfft%z_group_size_save)+current)*2:2*group_size+(((iset-1)*dfft%z_group_size_save)+current)*2) )

             CALL SYSTEM_CLOCK( time(15) )
  
             dfft%time_adding( 15 ) = dfft%time_adding( 15 ) + ( time(15) - time(14) )
  
          ENDDO
          CALL SYSTEM_CLOCK( auto_time(4) )
          dfft%auto_timings(1) = dfft%auto_timings(1) + ( auto_time(4) - auto_time(3) )
          dfft%auto_4Stimings(1) = dfft%auto_4Stimings(1) + ( auto_time(4) - auto_time(3) )
 
          !$  IF( dfft%rsactive ) THEN
          !$     locks_calc_2( dfft%my_node_rank+1, 1+(counter+dfft%num_buff-1)*dfft%batch_size_save:batch_size+(counter+dfft%num_buff-1)*dfft%batch_size_save ) = .false.
          !$omp flush( locks_calc_2 )
          !$  ELSE 
          !$     locks_calc_1( dfft%my_node_rank+1, 1+(counter+dfft%num_buff-1)*dfft%batch_size_save:batch_size+(counter+dfft%num_buff-1)*dfft%batch_size_save ) = .false.
          !$omp flush( locks_calc_1 )
          !$  END IF
  
       END IF
  
    END IF
  
  
  END SUBROUTINE fftpw_batch

  SUBROUTINE fftpw_4S( dfft, isign, step, batch_size, divparam_1, remswitch, mythread, counter, work_buffer, first_dim, priv, f_inout1, f_inout2, f_inout3, f_inout4, ip, jp )
  USE iso_fortran_env
    IMPLICIT NONE
  
    TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
    INTEGER, INTENT(IN) :: isign, step, batch_size, counter, work_buffer, first_dim
    INTEGER, INTENT(IN) :: divparam_1, remswitch, mythread
    INTEGER, INTENT(INOUT) :: priv(:)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout1(first_dim,*)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout2(:,:)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout3(:,:)
    COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: f_inout4(:,:)
    INTEGER,     OPTIONAL, INTENT(INOUT) :: ip, jp
  
    INTEGER :: ibatch, current, howmany
    LOGICAL :: last
    INTEGER(INT64) :: auto_time(4)
    INTEGER(INT64) :: time(17), cr
 
    current = (counter-1)*dfft%batch_size_save
  
    IF( isign .eq. -1 ) THEN !!  invfft
  
       IF( step .eq. 1 ) THEN

          CALL SYSTEM_CLOCK( time(1) )

          !$OMP Barrier 

          CALL SYSTEM_CLOCK( time(2) )
          dfft%time_adding( 22 ) = dfft%time_adding( 22 ) + ( time(2) - time(1) )

          CALL SYSTEM_CLOCK( time(2) )

          CALL Prepare_Psi_Man( dfft, f_inout2, f_inout1, dfft%ngms, remswitch, mythread, dfft%nsw, last, priv, ip, jp )

          CALL SYSTEM_CLOCK( time(3) )
          dfft%time_adding( 1 ) = dfft%time_adding( 1 ) + ( time(3) - time(2) )

          !In theory, locks could be made faster by checking each iset individually for non-remainder cases 
          !$omp flush( locks_calc_1 )
          !$  DO WHILE( ANY(locks_calc_1( :, 1+current:dfft%batch_size_save+current ) ) )
          !$omp flush( locks_calc_1 )
          !$  END DO

          CALL invfft_z_section_Man( dfft, f_inout1, f_inout3(:,work_buffer), f_inout4(:,work_buffer), 1, batch_size, remswitch, mythread, dfft%nsw )

          !$OMP Barrier
          !$  IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) THEN
          !$     locks_calc_inv( dfft%my_node_rank+1, counter ) = .false.
          !$omp flush( locks_calc_fw )
          !$  END IF
  
       ELSE IF( step .eq. 2 ) THEN
  
          CALL SYSTEM_CLOCK( time(4) )

          !$omp flush( locks_calc_inv )
          !$  DO WHILE( ANY( locks_calc_inv( :, counter ) ) )
          !$omp flush( locks_calc_inv )
          !$  END DO

          CALL SYSTEM_CLOCK( time(5) )
          dfft%time_adding( 23 ) = dfft%time_adding( 23 ) + ( time(5) - time(4) )
  
          CALL fft_com_Man( dfft, f_inout2(:,work_buffer), f_inout3(:,work_buffer), dfft%sendsize, dfft%my_node_rank, &
                            dfft%inter_node_comm, dfft%nodes_numb, dfft%my_inter_node_rank, dfft%non_blocking, work_buffer )
          CALL SYSTEM_CLOCK( time(16) )
          dfft%time_adding( 28 ) = dfft%time_adding( 28 ) + ( time(16) - time(5) )
  
          !$  locks_com_inv( counter ) = .false.
          !$omp flush( locks_com_inv )
  
       ELSE IF( step .eq. 3 ) THEN
  
          CALL SYSTEM_CLOCK( time(6) )

          !$omp flush( locks_com_inv )
          !$  DO WHILE( locks_com_inv( counter ) .and. .not. dfft%single_node )
          !$omp flush( locks_com_inv )
          !$  END DO

          CALL SYSTEM_CLOCK( time(7) )
          dfft%time_adding( 24 ) = dfft%time_adding( 24 ) + ( time(7) - time(6) )

          CALL invfft_y_section_Man( dfft, f_inout1, f_inout2(:,work_buffer), f_inout3, &
                                     dfft%map_acinv, dfft%map_acinv_rem, dfft%nr1w, divparam_1, 1, remswitch, mythread )

       ELSE IF( step .eq. 4 ) THEN

          CALL invfft_x_section_Man( dfft, f_inout1, remswitch, mythread )

          IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) THEN
             IF( dfft%vpsi ) THEN
                !$  locks_calc_2( dfft%my_node_rank+1, 1+current:dfft%batch_size_save+current ) = .false.
                !$omp flush( locks_calc_2 )
             ELSE
                !$  locks_calc_1( dfft%my_node_rank+1, 1+current+(dfft%batch_size_save*dfft%buffer_size_save):dfft%batch_size_save+current+(dfft%batch_size_save*dfft%buffer_size_save) ) = .false.
                !$omp flush( locks_calc_1 )
             END IF
          END IF

          !$OMP Barrier 
       END IF
  
    ELSE !! fw fft
  
       IF( step .eq. 1 ) THEN

!          !$OMP Barrier !Should be here?
          CALL fwfft_x_section_Man( dfft, f_inout2, f_inout3, dfft%nr1w, remswitch, mythread )

       ELSE IF( step .eq. 2 ) THEN

          CALL SYSTEM_CLOCK( time(8) )

          !$omp flush( locks_calc_2 )
          !$  DO WHILE( ANY(locks_calc_2(:,1+current:dfft%batch_size_save+current ) ) )
          !$omp flush( locks_calc_2 )
          !$  END DO
  
          CALL SYSTEM_CLOCK( time(9) )
          dfft%time_adding( 25 ) = dfft%time_adding( 25 ) + ( time(9) - time(8) )
  
          CALL fwfft_y_section_Man( dfft, f_inout4, f_inout2(:,work_buffer), f_inout3(:,work_buffer), &
                                    dfft%map_pcfw, dfft%nr1w, dfft%nsw, batch_size, divparam_1, 1, remswitch, mythread )

          !$OMP Barrier
          !$  IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) THEN
          !$     locks_calc_fw( dfft%my_node_rank+1, counter ) = .false.
          !$omp flush( locks_calc_fw )
          !$  END IF

       ELSE IF( step .eq. 3 ) THEN

          CALL SYSTEM_CLOCK( time(10) )

          !$omp flush( locks_calc_fw )
          !$  DO WHILE( ANY( locks_calc_fw( :, counter ) ) )
          !$omp flush( locks_calc_fw )
          !$  END DO
  
          CALL SYSTEM_CLOCK( time(11) )
          dfft%time_adding( 26 ) = dfft%time_adding( 26 ) + ( time(11) - time(10) )
     
          CALL fft_com_Man( dfft, f_inout2(:,work_buffer), f_inout3(:,work_buffer), dfft%sendsize, dfft%my_node_rank, &
                            dfft%inter_node_comm, dfft%nodes_numb, dfft%my_inter_node_rank, dfft%non_blocking, work_buffer )
          CALL SYSTEM_CLOCK( time(17) )
          dfft%time_adding( 29 ) = dfft%time_adding( 29 ) + ( time(17) - time(11) )

          !$  locks_com_fw( counter ) = .false.
          !$omp flush( locks_com_fw )
  
       ELSE IF( step .eq. 4 ) THEN
  
          CALL SYSTEM_CLOCK( time(12) )

          !$omp flush( locks_com_fw )
          !$  DO WHILE( locks_com_fw( counter ) .and. .not. dfft%single_node )
          !$omp flush( locks_com_fw )
          !$  END DO

          CALL SYSTEM_CLOCK( time(13) )
          dfft%time_adding( 27 ) = dfft%time_adding( 27 ) + ( time(13) - time(12) )
  
          CALL fwfft_z_section_Man( dfft, f_inout4(:,work_buffer), f_inout1, 1, batch_size, remswitch, mythread, dfft%nsw )

          CALL SYSTEM_CLOCK( time(14) )

          IF( dfft%nthreads .eq. 1 .or. mythread .eq. 1 ) THEN
             !$  IF( dfft%rsactive ) THEN
             !$     locks_calc_2( dfft%my_node_rank+1, 1+(counter+dfft%num_buff-1)*dfft%batch_size_save:batch_size+(counter+dfft%num_buff-1)*dfft%batch_size_save ) = .false.
             !$omp flush( locks_calc_2 )
             !$  ELSE 
             !$     locks_calc_1( dfft%my_node_rank+1, 1+(counter+dfft%num_buff-1)*dfft%batch_size_save:dfft%batch_size_save+(counter+dfft%num_buff-1)*dfft%batch_size_save ) = .false.
             !$omp flush( locks_calc_1 )
             !$  END IF
          END IF

          !$OMP Barrier 

          CALL Accumulate_Psi_Man( dfft, f_inout1, f_inout3, dfft%ngms, divparam_1, last, dfft%nsw, f_inout2, mythread )

          CALL SYSTEM_CLOCK( time(15) )
  
          dfft%time_adding( 15 ) = dfft%time_adding( 15 ) + ( time(15) - time(14) )
 
       END IF
  
    END IF
  
  
  END SUBROUTINE fftpw_4S

  SUBROUTINE fftpw( isign, dfft, f, ng, nr1s, ir1, ns )

!    USE mpi_f08

    IMPLICIT NONE
  
    TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
    INTEGER, INTENT(IN) :: isign, ng
    COMPLEX(DP), INTENT(INOUT) :: f(:)
    INTEGER, INTENT(IN) :: ir1(:), ns(:), nr1s(:)

    COMPLEX(DP), ALLOCATABLE, SAVE :: temp(:)

    COMPLEX(DP), POINTER, SAVE, CONTIGUOUS :: shared1(:)
    COMPLEX(DP), POINTER, SAVE, CONTIGUOUS :: shared2(:)

    INTEGER :: i, ierr
    CHARACTER(*), PARAMETER                  :: procedureN = 'fftpw'

    dfft%singl = .true. 


    IF( .not. allocated( dfft%aux ) ) THEN
  
       ALLOCATE( temp( dfft%nnr ) )

       CALL Set_Req_Vals( dfft, dfft%nstate, 1, dfft%rem_size, 1, ir1, ns )
       CALL Prep_Copy_Maps( dfft, ng, 1, dfft%rem_size, ir1, ns )
  
       dfft%sendsize = MAXVAL ( dfft%nr3p ) * MAXVAL( dfft%nsp ) * dfft%node_task_size * dfft%node_task_size * 1
     
       CALL create_shared_memory_window_1d( shared1, 80, dfft, dfft%sendsize*dfft%nodes_numb ) 
       CALL create_shared_memory_window_1d( shared2, 81, dfft, dfft%sendsize*dfft%nodes_numb ) 
     
    END IF

    IF( isign .eq. -1 ) THEN !!  invfft

       dfft%aux2 = f

       IF( dfft%single_node ) CALL MPI_BARRIER( dfft%comm, ierr )
       IF( .not. dfft%single_node ) CALL MPI_BARRIER( dfft%node_comm, ierr )
CALL MPI_BARRIER( dfft%comm, ierr )

       CALL invfft_pre_com( dfft, dfft%aux2, shared1, shared2, 1, 1, 1, ns )

       IF( dfft%single_node ) THEN
          CALL MPI_BARRIER( dfft%comm, ierr )
       ELSE 
CALL MPI_BARRIER( dfft%comm, ierr )
          CALL MPI_BARRIER( dfft%node_comm, ierr )
          CALL fft_com( dfft, shared1, shared2, dfft%sendsize, dfft%my_node_rank, &
                        dfft%inter_node_comm, dfft%nodes_numb, dfft%my_inter_node_rank, dfft%non_blocking, -1 )
          CALL MPI_BARRIER( dfft%node_comm, ierr )
CALL MPI_BARRIER( dfft%comm, ierr )
       END IF
   
       CALL invfft_after_com( dfft, f, shared2, dfft%aux2, dfft%map_acinv_one, dfft%map_acinv_rem_one, 1, 1, 1, 1, nr1s )
CALL MPI_BARRIER( dfft%comm, ierr )
       IF( .not. dfft%single_node ) CALL MPI_BARRIER( dfft%node_comm, ierr )
       IF( dfft%single_node ) CALL MPI_BARRIER( dfft%comm, ierr )

    ELSE !! fw fft
    
       IF( dfft%single_node ) CALL MPI_BARRIER( dfft%comm, ierr )
       IF( .not. dfft%single_node ) CALL MPI_BARRIER( dfft%node_comm, ierr )
CALL MPI_BARRIER( dfft%comm, ierr )

       CALL fwfft_pre_com( dfft, f, dfft%aux2, shared1, shared2, 1, 1, 1, nr1s, ns )
    
CALL MPI_BARRIER( dfft%comm, ierr )
       IF( dfft%single_node ) THEN
          CALL MPI_BARRIER( dfft%comm, ierr )
       ELSE 
          CALL MPI_BARRIER( dfft%node_comm, ierr )
          CALL fft_com( dfft, shared1, shared2, dfft%sendsize, dfft%my_node_rank, &
                        dfft%inter_node_comm, dfft%nodes_numb, dfft%my_inter_node_rank, dfft%non_blocking, -1 )
          CALL MPI_BARRIER( dfft%node_comm, ierr )
       END IF
CALL MPI_BARRIER( dfft%comm, ierr ) 
    
       CALL fwfft_after_com( dfft, shared2, dfft%aux2, 1, 1, 1, ns )
CALL MPI_BARRIER( dfft%comm, ierr )

       IF( .not. dfft%single_node ) CALL MPI_BARRIER( dfft%node_comm, ierr )
       IF( dfft%single_node ) CALL MPI_BARRIER( dfft%comm, ierr )

       f = dfft%aux2 * dfft%tscale 
    
    END IF

    dfft%singl = .false. 
  
  END SUBROUTINE fftpw

END MODULE fftmain_utils
