#include "cpmd_global.h"

MODULE vpsi_utils
  USE cdftmod,                         ONLY: cdftcom,&
                                             cdftlog,&
                                             wdiff
  USE cnst,                            ONLY: uimag
  USE cp_cuda_types,                   ONLY: cp_cuda_devices_fft,&
                                             cp_cuda_env
  USE cp_cufft_types,                  ONLY: cp_cufft,&
                                             cp_cufft_device_get_ptrs,&
                                             cp_cufft_stream_get_ptrs
  USE cp_cuvpsi_types,                 ONLY: cp_cuvpsi_t
  USE cp_cuvpsi_utils,                 ONLY: cp_cuapply_potential,&
                                             cp_cuvpotdg_copy_to_device,&
                                             cp_cuvpsi_finalize_and_dealloc,&
                                             cp_cuvpsi_init_and_alloc
  USE cp_cuwfn_types,                  ONLY: cp_cuwfn,&
                                             cp_cuwfn_device_get_ptrs,&
                                             cp_cuwfn_get
  USE cp_grp_utils,                    ONLY: cp_grp_redist,&
                                             cp_grp_redist_array
  USE cppt,                            ONLY: gk,&
                                             hg,&
                                             indzs,&
                                             nzh,&
                                             nzhs
  USE cuda_types,                      ONLY: cuda_memory_t,&
                                             cuda_stream_t
  USE cuda_utils,                      ONLY: cuda_alloc_host,&
                                             cuda_dealloc_host,&
                                             cuda_mem_zero_bytes,&
                                             cuda_z_points_to
  USE cuuser_utils,                    ONLY: CuUser_setpsi_1_state_g,&
                                             CuUser_setpsi_2_states_g
  USE dg,                              ONLY: ipooldg,&
                                             tdgcomm
  USE efld,                            ONLY: extf
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: inzf,&
                                             jgw,&
                                             jhg,&
                                             llr1,&
                                             ngrm,&
                                             nzff,&
                                             msrays,&
                                             wfn_r,&
                                             wfn_g,&
                                             fft_batchsize,&
                                             fft_numbatches,&
                                             fft_residual,&
                                             fft_total,&
                                             fft_tune_num_it,&
                                             fft_time_total,&
                                             xf,&
                                             yf,&
                                             locks_inv,&
                                             locks_fw
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn,&
                                             fwfftn_batch,&
<<<<<<< HEAD
                                             invfftn_batch
=======
                                             fwfftn_batch_com,&
                                             invfftn_batch,&
                                             comm_send,&
                                             comm_recv,&
                                             locks_calc_inv,&
                                             locks_calc_fw,&
                                             locks_com_inv,&
                                             locks_com_fw,&
                                             locks_calc_1,&
                                             locks_calc_2,&
                                             invfft_pwbatch,&
                                             fwfft_pwbatch
>>>>>>> b8d086d55 (FFT works, application of v wrong/missing)
  USE fftnew_utils,                    ONLY: setfftn
  USE fftpw_base,                      ONLY: dfft,&
                                             wfn_real
  USE fftpw_batching,                  ONLY: Apply_v
  USE fftpw_make_maps,                 ONLY: Prep_copy_Maps,&
                                             Set_Req_Vals,&
                                             MapVals_CleanUp 
  USE fftpw_param,                     ONLY: DP
  USE fftpw_types,                     ONLY: create_shared_memory_window_2d,&
                                             create_shared_locks_2d,&
                                             create_shared_locks_1d,&
                                             Clean_up_shared
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE kpclean_utils,                   ONLY: c_clean
  USE kpnt,                            ONLY: hgkm,&
                                             hgkp,&
                                             rk
  USE kpts,                            ONLY: tkpts
  USE fft_maxfft,                      ONLY: maxfft
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_comm_dup,&
                                             mp_comm_free,&
                                             mp_comm_null,&
                                             mp_bcast
  USE parac,                           ONLY: parai
  USE part_1d,                         ONLY: part_1d_get_el_in_blk,&
                                             part_1d_nbr_el_in_blk
  USE prcp,                            ONLY: prcp_com
  USE reshaper,                        ONLY: reshape_inplace
  USE rswfmod,                         ONLY: maxstates,&
                                             rsactive,&
                                             rswf
  USE special_functions,               ONLY: cp_erf
  USE spin,                            ONLY: clsd,&
                                             lspin2,&
                                             spin_mod
  USE state_utils,                     ONLY: add_wfn,&
                                             set_psi_1_state_g,&
                                             set_psi_1_state_g_kpts,&
                                             set_psi_2_states_g,&
                                             set_psi_batch_g
  USE switch_functionals_utils,        ONLY: switch_functionals
  USE system,                          ONLY: &
       cntl, fpar, group, locpot2, ncpw, nkpt, parap, parm, spar
  USE td_input,                        ONLY: gampl,&
                                             gndir,&
                                             gpot,&
                                             gpotv,&
                                             td_prop
  USE thread_view_types,               ONLY: thread_view_get,&
                                             thread_view_t
  USE thread_view_utils,               ONLY: thread_views_finalize,&
                                             thread_views_init
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vpsi_lse_utils,                  ONLY: vpsi_lse
  USE vtaupsi_utils,                   ONLY: vtaupsi
  USE zeroing_utils,                   ONLY: zeroing

  !$ USE omp_lib, ONLY: omp_get_max_threads, omp_get_thread_num, &
  !$                    omp_get_num_threads, omp_set_num_threads
#ifdef _HASNT_OMP_SET_NESTED
  !$ USE omp_lib, ONLY: omp_get_max_active_levels, &
  !$                    omp_set_max_active_levels
#else
  !$ USE omp_lib, ONLY: omp_get_max_active_levels, omp_get_nested, &
  !$                    omp_set_max_active_levels, omp_set_nested
#endif

  USE nvtx_utils
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif
#ifdef __PARALLEL
  USE mpi_f08
#endif
#ifdef _INTEL_MKL
  use mkl_service
#endif

 USE iso_fortran_env

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vpsi
  PUBLIC :: vpsi_batchfft
  PUBLIC :: vpsimt
  PUBLIC :: do_the_vpsi_thing
  !public :: movepsid

CONTAINS

  ! if (.not.tdg) r*8 vpotdg(nnr1,*) => r*8 vpot(nnr1,*)
  ! else vpotdg = copy of psi(:: 2)
  !
  ! r*8 vpotx(*) => r*8 vpotdg(1,1)
  ! r*8 psix(*) => c*8 psi(*) for use in PARALLEL sections
  ! ==================================================================
  SUBROUTINE vpsi(c0,c2,f,vpot,psi,nstate,ikind,ispin,redist_c2)
    ! ==================================================================
    ! == K-POINT AND NOT K-POINT VERSION OF VPSI.                     ==
    ! ==--------------------------------------------------------------==
    ! == VPOT:   IN INPUT POTENTIAL                                   ==
    ! == ISPIN:  Need with LSD option for diagonalization scheme      ==
    ! ==         dimension of VPOT(NNR1,ISPIN)                        ==
    ! ==--------------------------------------------------------------==
    ! EHR[
    ! EHR]
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:)
    REAL(real_8)                             :: f(:)
    REAL(real_8), TARGET __CONTIGUOUS        :: vpot(:,:)
    COMPLEX(real_8), TARGET __CONTIGUOUS     :: psi(:)
    INTEGER                                  :: nstate, ikind, ispin
    LOGICAL                                  :: redist_c2

    CHARACTER(*), PARAMETER                  :: procedureN = 'vpsi'
    COMPLEX(real_8), PARAMETER               :: zone = (1.0_real_8,0.0_real_8)

    COMPLEX(real_8)                          :: fm, fp, psii, psin
    COMPLEX(real_8), ALLOCATABLE             :: C2_vpsi(:,:)
    COMPLEX(real_8), DIMENSION(:), &
      POINTER __CONTIGUOUS                   :: psi_p
    COMPLEX(real_8), DIMENSION(:, :), &
      POINTER __CONTIGUOUS                   :: psis
#ifdef __PARALLEL
    INTEGER :: device_idx, i, i_stream, iclpot = 0, id, ierr, ig, &
      ir, is1, is2, isub, isub2, isub3, isub4, iwf, ixx, ixxs, iyy, izz, jj, &
      leadx, lspin, n_max_threads, n_nested_threads, n_streams_per_task, &
      njump, nnrx, nostat, nrxyz1s, nrxyz2, stream_idx
    type(MPI_COMM)                           :: fft_comm
#else
    INTEGER :: device_idx, fft_comm, i, i_stream, iclpot = 0, id, ierr, ig, &
      ir, is1, is2, isub, isub2, isub3, isub4, iwf, ixx, ixxs, iyy, izz, jj, &
      leadx, lspin, n_max_threads, n_nested_threads, n_streams_per_task, &
      njump, nnrx, nostat, nrxyz1s, nrxyz2, stream_idx
#endif
    INTEGER, SAVE                            :: dgfirst = 0
    LOGICAL :: copy_data_to_device, copy_data_to_device_inv, &
      copy_data_to_host, lg_vpotx3a, lg_vpotx3b, ready, wfcopy
    REAL(real_8)                             :: chksum, csmult, fi, fip1, g2, &
                                                xskin
    REAL(real_8), ALLOCATABLE                :: vpotx3a(:,:,:), vpotx3b(:,:,:)
    REAL(real_8), POINTER __CONTIGUOUS       :: VPOTX(:)
    REAL(real_8), POINTER, SAVE __CONTIGUOUS :: vpotdg(:,:)
    TYPE(cp_cuvpsi_t)                        :: cp_cuvpsi
    TYPE(cuda_memory_t), POINTER             :: c0_d, inzs_d, nzfs_d, psi_d
    TYPE(cuda_stream_t), POINTER             :: stream_p
    TYPE(thread_view_t)                      :: thread_view
    TYPE(thread_view_t), ALLOCATABLE, &
      DIMENSION(:)                           :: thread_views

    !$ LOGICAL :: nested_orig
    !$ INTEGER :: max_level_orig
    ! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )
    __NVTX_TIMER_START ( procedureN//'_a' )

    IF (group%nogrp.GT.1)CALL stopgm(procedureN,&
         'OLD TASK GROUPS NOT SUPPORTED ANYMORE ',&
         __LINE__,__FILE__)
    IF (tdgcomm%tdg) THEN
       IF (lspin2%tlse)CALL stopgm(procedureN,&
            'TDG TOGETHER WITH  LSE IS NOT YET IMPLEMENTED ',&
            __LINE__,__FILE__)
       CALL setfftn(ipooldg)
       IF (dgfirst.EQ.0) THEN
          ALLOCATE(vpotdg(fpar%nnr1,2),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          dgfirst=1
       ENDIF
    ELSE
       CALL setfftn(0)
    ENDIF
    IF (lspin2%tlse.AND.lspin2%tlsets.AND.(lspin2%tross.OR.lspin2%tcas22.OR.lspin2%tpenal.OR.lspin2%troot)) THEN
       CALL stopgm(procedureN,&
            'NO SLATER TS WITH ROSS, CAS22, PENALTY, ROOTHAAN',&
            __LINE__,__FILE__)
    ENDIF

    !
    ! Need to allocate a buffer to accumulate the vpsi part of C2
    !
    ALLOCATE(C2_vpsi(nkpt%ngwk,nstate),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',&
         __LINE__,__FILE__)

    CALL zeroing(C2_vpsi)!, nkpt%ngwk*nstate)

    ! For the moment we synchronize the C0 and C2 if the cp groups are used
    ! That should be done else where
    ! IF(cp_nogrp.GT.1) THEN
    ! CALL MY_BCAST(C0,2*NGW*NSTATE*8/IRAT,0,CP_INTER_GRP)
    ! CALL MY_BCAST(C2,2*NGW*NSTATE*8/IRAT,0,CP_INTER_GRP)
    ! ENDIF

    leadx = fpar%nnr1
    nnrx  = llr1
    !
    IF (tdgcomm%tdg) THEN
       CALL zeroing(psi)!, maxfft)
       CALL dcopy (fpar%nnr1, vpot, 1, psi, 2)
       CALL fwfftn(psi,.FALSE.,parai%allgrp)
       CALL movepsid(psi)
       CALL invfftn(psi, .FALSE.,parai%allgrp)
       CALL dcopy(llr1, psi, 2, vpotdg, 1)
       IF (cntl%tlsd.AND.ispin.EQ.2) THEN
          CALL zeroing(psi)!, maxfft)
          CALL dcopy (fpar%nnr1, vpot(1,2), 1, psi, 2)
          CALL fwfftn(psi,.FALSE.,parai%allgrp)
          CALL movepsid(psi)
          CALL invfftn(psi, .FALSE.,parai%allgrp)
          CALL dcopy(llr1, psi, 2, vpotdg(1,2), 1)
       ENDIF
    ELSE
       vpotdg => vpot
    ENDIF
    !
    !vpotx(1:SIZE(vpotdg)) => vpotdg
    CALL reshape_inplace(vpotdg, (/fpar%nnr1*ispin/), vpotx)

    njump=2
    IF (tkpts%tkpnt) njump=1
    IF (lspin2%tlse) THEN
       nostat = clsd%ialpha-1
    ELSE
       nostat = nstate
    ENDIF

    IF (cntl%cdft)THEN
       IF (.NOT.cdftlog%tcall)THEN
          csmult=cdftcom%cdft_v(1)
       ELSE
          csmult=cdftcom%cdft_v(1)+cdftcom%cdft_v(2)
       ENDIF
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          vpotdg(ir,1)=vpotdg(ir,1)+csmult*wdiff(ir)
       ENDDO
       IF (.NOT.cdftlog%tcall)THEN
          IF (cdftlog%tspinc) THEN
             csmult=-cdftcom%cdft_v(1)
          ELSE
             csmult=cdftcom%cdft_v(1)
          ENDIF
       ELSE
          csmult=cdftcom%cdft_v(1)-cdftcom%cdft_v(2)
       ENDIF
       IF (cntl%tlsd.AND.ispin.EQ.2) THEN
          !$omp parallel do private(IR)
          DO ir=1,fpar%nnr1
             vpotdg(ir,2)=vpotdg(ir,2)+csmult*wdiff(ir)
          ENDDO
       ENDIF
    ENDIF

    ! ==--------------------------------------------------------------==

    !acm gpu version:
    !acm if possible do the accumulation of the potential on the gpu
    copy_data_to_host   = (cntl%tlsd.AND.ispin.EQ.2) .OR. ( .NOT. cp_cuda_env%use_fft )
    copy_data_to_device = copy_data_to_host
    CALL cp_cuwfn_get ( cp_cuwfn, ready=ready, chksum=chksum )
    copy_data_to_device_inv = .NOT. ( ready .AND. cp_cuda_env%use_fft .AND. &
         &                           .NOT. tkpts%tkpnt .AND. &
         &                            ABS( SUM( ABS( c0(:,1) ) ) - chksum ) < 1.0e-12_real_8 )
    !write(*,*) 'vpsi_utils.mod.F90: copy_data_to_device_inv=',copy_data_to_device_inv
    n_streams_per_task = MAX( cp_cuda_env%fft_n_streams_per_device * cp_cuda_env%fft_n_devices_per_task, 1 )!vw fix the min
    IF( .NOT. cp_cuda_env%use_fft ) n_streams_per_task = 1

    n_max_threads = 1
    !$ n_max_threads = omp_get_max_threads()
    IF( n_streams_per_task > n_max_threads ) CALL stopgm(procedureN,&
         'number of streams creater than max number of threads',&
         __LINE__,__FILE__)

    n_nested_threads = n_max_threads / n_streams_per_task

    fft_comm = mp_comm_null
    CALL thread_views_init ( cp_cuda_env%fft_n_devices_per_task, cp_cuda_env%fft_n_streams_per_device, thread_views )

    IF( cp_cuda_env%use_fft ) THEN
       CALL tiset(procedureN//'_gpu_init',isub2)
#ifndef _HASNT_OMP_SET_NESTED
       !$ nested_orig = omp_get_nested ( )
#endif
       !$ max_level_orig = omp_get_max_active_levels( )
#ifndef _HASNT_OMP_SET_NESTED
       !$ CALL omp_set_nested( .TRUE. )
#endif
       !$ CALL omp_set_max_active_levels( 2 )
       CALL cuda_alloc_host(psis,[SIZE(psi),n_streams_per_task], lbs=[1,0])

       IF( .NOT. copy_data_to_host ) THEN
          ! TODO acm: the init+allocation/deallocation should be done outside (once for entire code)
          CALL cp_cuvpsi_init_and_alloc( cp_cuda_env, cp_cuda_devices_fft, cp_cuvpsi, SIZE(vpotdg) )
          CALL cp_cuvpotdg_copy_to_device( cp_cuda_env, cp_cuvpsi, vpotdg )
       ENDIF

       IF( .NOT. copy_data_to_device_inv ) THEN
          DO device_idx = 1, cp_cuda_env%fft_n_devices_per_task
             !CALL cp_cuwfn_device_get_ptrs ( cp_cuwfn, device_idx, c0_d=c0_d )
             !CALL cuda_memcpy_host_to_device ( c0, c0_d )
          ENDDO
       ENDIF
       CALL tihalt(procedureN//'_gpu_init',isub2)
    ENDIF

    ! acm copy vpotdg on device
    __NVTX_TIMER_STOP

    !$omp parallel   if( cp_cuda_env%use_fft ) &
    !$omp            num_threads( n_streams_per_task ) &
    !$omp            default( none ) &
    !$omp            private( i, is1, is2, iwf, wfcopy, i_stream, psi_p, &
    !$omp                     xskin, fip1, fm, fp, fi, ixxs, izz, iyy, jj, ixx, &
    !$omp                     lg_vpotx3b, lg_vpotx3a, nrxyz2, nrxyz1s, iclpot, &
    !$omp                     ierr, lspin, vpotx3a, vpotx3b, g2, psii, psin, thread_view, &
    !$omp                     stream_idx, device_idx, stream_p, psi_d, c0_d, nzfs_d, &
    !$omp                     inzs_d, ig, id, ir ) &
    !$omp            firstprivate( fft_comm ) &
    !$omp            shared( rsactive, hg, jgw, prcp_com, geq0, gpotv, gampl, &
    !$omp                    gpot, gk, rk, hgkm, hgkp, indzs, nzhs, ncpw, gndir, &
    !$omp                    f, parap, parm, spar, locpot2, vpotdg, td_prop, spin_mod, &
    !$omp                    cntl, ispin, psi, psis, parai, nostat, njump, n_streams_per_task, &
    !$omp                    ikind, nstate, maxstates, tkpts, nnrx, rswf, c0, vpotx, &
    !$omp                    leadx, extf, c2_vpsi, cp_cufft, cp_cuvpsi, cp_cuda_env, n_nested_threads, &
    !$omp                    copy_data_to_host, copy_data_to_device, copy_data_to_device_inv, &
    !$omp                    thread_views, cp_cuwfn )

    !vw set number of children threads
    !$ CALL omp_set_num_threads ( n_nested_threads )

    i_stream = 0
    !$ i_stream = omp_get_thread_num()
    IF(cp_cuda_env%fft_n_devices_per_task>0) THEN
       thread_view = thread_views( i_stream ) !vw need to fix that
    ENDIF

    !vw create local comm for each stream/thread
    DO i = 0, n_streams_per_task - 1
       IF( i == i_stream ) THEN
          CALL mp_comm_dup ( parai%allgrp, fft_comm )
       ENDIF
       !$omp barrier
    ENDDO

    IF( cp_cuda_env%use_fft ) THEN
       psi_p => psis( :, i_stream )
    ELSE
       psi_p => psi( : )
    ENDIF

    IF( .NOT. copy_data_to_device_inv ) THEN
       CALL thread_view_get( thread_view, device_idx=device_idx, stream_idx=stream_idx )
       CALL cp_cuwfn_device_get_ptrs ( cp_cuwfn, device_idx, c0_d=c0_d )
       CALL cp_cufft_device_get_ptrs ( cp_cufft, device_idx, nzfs_d=nzfs_d, inzs_d=inzs_d )
       CALL cp_cufft_stream_get_ptrs ( cp_cufft, device_idx, stream_idx, stream=stream_p, t1_d=psi_d )
    ENDIF

    !$omp do schedule( static, 1 )
    DO i = 1,part_1d_nbr_el_in_blk(nostat,parai%cp_inter_me,parai%cp_nogrp),njump
       is1 = part_1d_get_el_in_blk(i,nostat,parai%cp_inter_me,parai%cp_nogrp)
       is2 = nostat+1
       IF (njump.EQ.2) THEN
          IF (i+1.LE.part_1d_nbr_el_in_blk(nostat,parai%cp_inter_me,parai%cp_nogrp))&
               is2 = part_1d_get_el_in_blk(i+1,nostat,parai%cp_inter_me,parai%cp_nogrp)
       ENDIF

       wfcopy=.FALSE.
       IF (rsactive) THEN
          iwf = is1+(ikind-1)*nstate
          IF (iwf.LE.maxstates) THEN
             wfcopy=.TRUE.
             copy_data_to_host=.TRUE.
             copy_data_to_device=.TRUE.
             copy_data_to_device_inv=.TRUE.
          ENDIF
       ENDIF
       IF (wfcopy) THEN
          iwf=is1/2+1
          IF (tkpts%tkpnt) iwf = is1+(ikind-1)*nstate
          CALL dcopy(2*nnrx,rswf(1,iwf),1,psi_p(1),1)
       ELSE
          __NVTX_TIMER_START ( procedureN//'_pre' )
          ! ==----------------------------------------------------------==
          ! == Store the wavefunctions in array PSI                     ==
          ! == in the following way:                                    ==
          ! == Since we use the Gamma point only, the wavefunctions in  ==
          ! == real space are real. This fact is used to save time by   ==
          ! == using complex arrays where the wavefunctions of          ==
          ! == 2 different states (i and i+1) are combined, and         ==
          ! == performing complex Fourier Transforms.                   ==
          ! == Warning: If you use K-points other than Gamma, or zone   ==
          ! == boundary, this symmetry is broken and the trick is not   ==
          ! == valid anymore. First step is to build the array PSI.     ==
          ! ==----------------------------------------------------------==
          IF( .NOT. copy_data_to_device_inv ) THEN
             CALL cuda_mem_zero_bytes ( psi_d, psi_d%n_bytes )
             IF (tkpts%tkpnt) THEN
                CALL stopgm(procedureN,'k-points with wfn on GPU NYI',&
                     __LINE__,__FILE__)
             ELSE
                IF (is2.GT.nostat) THEN
                   CALL CuUser_setpsi_1_state_g ( zone, cuda_z_points_to( c0_d, (is1-1)*SIZE(c0,1)+1 ), &
                        & psi_d, jgw, nzfs_d, inzs_d, geq0, stream_p )
                ELSE
                   CALL CuUser_setpsi_2_states_g ( cuda_z_points_to( c0_d, (is1-1)*SIZE(c0,1)+1 ), &
                        & cuda_z_points_to( c0_d, (is2-1)*SIZE(c0,1)+1 ), psi_d, jgw, nzfs_d, inzs_d, geq0, stream_p )
                ENDIF
             ENDIF
          ELSE
             CALL zeroing(psi_p)!, maxfft)
             IF (tkpts%tkpnt) THEN
                CALL set_psi_1_state_g_kpts(zone,c0(:,is1),psi_p)
             ELSE
                IF (is2.GT.nostat) THEN
                   CALL set_psi_1_state_g(zone,c0(:,is1),psi_p)
                ELSE
                   CALL set_psi_2_states_g(c0(:,is1),c0(:,is2),psi_p)
                ENDIF
             ENDIF
          ENDIF
          __NVTX_TIMER_STOP
          ! ==----------------------------------------------------------==
          ! == Transform the wavefunction to real space                 ==
          ! ==----------------------------------------------------------==
          CALL invfftn(psi_p,.TRUE., fft_comm, thread_view=thread_view, copy_data_to_host=copy_data_to_host, &
               & copy_data_to_device=copy_data_to_device_inv )
       ENDIF
       ! ==------------------------------------------------------------==
       ! == Apply the potential (V), which acts in real space.         ==
       ! ==------------------------------------------------------------==
       __NVTX_TIMER_START ( procedureN//'_apply' )
       IF (cntl%tlsd.AND.ispin.EQ.2) THEN
          IF (is1.EQ.spin_mod%nsup) THEN
             ! EHR[
             IF (td_prop%td_extpot) THEN
                !$omp parallel do private(IR)
                DO ir=1,nnrx
                   psi_p(ir)=(vpotx(ir)+extf(ir))*REAL(psi_p(ir))&
                        +uimag*(vpotx(leadx+ir)+extf(ir))*AIMAG(psi_p(ir))
                ENDDO
             ELSE
                !$omp parallel do private(IR)
                DO ir=1,nnrx
                   psi_p(ir)=vpotx(ir)*REAL(psi_p(ir))&
                        +uimag*vpotx(leadx+ir)*AIMAG(psi_p(ir))
                ENDDO
             ENDIF
             ! EHR]
          ELSE
             !lspin=0
             !IF (is1.GT.spin_mod%nsup) lspin=leadx
             lspin=1
             IF (is1.GT.spin_mod%nsup) lspin=2
             !$omp parallel do private(IR)
             DO ir=1,nnrx
                psi_p(ir)=vpotdg(ir,lspin)*psi_p(ir)
             ENDDO
          ENDIF
       ELSE
          IF( .NOT. copy_data_to_host ) THEN
             CALL cp_cuapply_potential ( nnrx, cp_cufft, cp_cuvpsi, thread_view )
          ELSE
             !$omp parallel do private(IR) shared(NNRX,PSI_p)
             DO ir=1,nnrx
                psi_p(ir)=vpotdg(ir,1)*psi_p(ir)
             ENDDO
          ENDIF
          ! kk-mb === print local potential (start) ===
          IF (locpot2%tlpot) THEN
             iclpot=iclpot+1
             nrxyz1s=spar%nr1s*spar%nr2s*spar%nr3s
             nrxyz2=0
             IF (iclpot .EQ. 1) THEN
                ALLOCATE(vpotx3a(spar%nr1s,spar%nr2s,spar%nr3s),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                CALL zeroing(vpotx3a)!,nrxyz1s)
                lg_vpotx3a=.TRUE.
                lg_vpotx3b=.FALSE.
             ELSE
                IF (lg_vpotx3a) THEN
                   ALLOCATE(vpotx3b(spar%nr1s,spar%nr2s,spar%nr3s),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   CALL zeroing(vpotx3b)!,nrxyz1s)
                   DEALLOCATE(vpotx3a,STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                        __LINE__,__FILE__)
                   lg_vpotx3a=.FALSE.
                   lg_vpotx3b=.TRUE.
                ELSE IF (lg_vpotx3b) THEN
                   ALLOCATE(vpotx3a(spar%nr1s,spar%nr2s,spar%nr3s),STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                        __LINE__,__FILE__)
                   CALL zeroing(vpotx3a)!,nrxyz1s)
                   DEALLOCATE(vpotx3b,STAT=ierr)
                   IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                        __LINE__,__FILE__)
                   lg_vpotx3a=.TRUE.
                   lg_vpotx3b=.FALSE.
                ENDIF
             ENDIF
             DO ir=1,nnrx
                IF (vpotx(ir) .NE. 0) THEN
                   nrxyz2=nrxyz2+1
                   ixx=MOD(nrxyz2-1,parm%nr1)+1
                   jj=(nrxyz2-1)/parm%nr1+1
                   iyy=MOD(jj-1,spar%nr2s)+1
                   izz=(nrxyz2-1)/(parm%nr1*spar%nr2s)+1
                   ixxs=ixx+parap%nrxpl(parai%mepos,1)-1
                   IF (lg_vpotx3a) THEN
                      vpotx3a(ixxs,iyy,izz)=vpotx(ir)
                   ELSE IF (lg_vpotx3b) THEN
                      vpotx3b(ixxs,iyy,izz)=vpotx(ir)
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
          ! kk-mb === print local potential (end)  ===
       ENDIF
       __NVTX_TIMER_STOP
       ! ==------------------------------------------------------------==
       ! == Back transform to reciprocal space the product V.PSI       ==
       ! ==------------------------------------------------------------==
       CALL fwfftn(psi_p,.TRUE., fft_comm, thread_view=thread_view, copy_data_to_device=copy_data_to_device )
       ! ==------------------------------------------------------------==
       ! == Decode the combination of wavefunctions (now multiplied by ==
       ! == the potential), back to those of states i and i+1. Multiply==
       ! == also by the occupation numbers f(i) and other constants in ==
       ! == order to obtain the force on the electronic degrees of     ==
       ! == freedom, stored in array FC.                               ==
       ! == Here we also add the force from the kinetic energy.        ==
       ! ==------------------------------------------------------------==
       __NVTX_TIMER_START ( procedureN//'_acc' )
       IF (tkpts%tkpnt) THEN
          fi=f(is1)
          IF (fi.EQ.0._real_8) fi=2._real_8
          ! EHR[
          IF (cntl%tgaugep.OR.cntl%tgaugef) THEN
             IF (gndir.EQ.0) THEN
                DO ig=1,ncpw%ngw
                   fp=psi_p(nzhs(ig))
                   fm=psi_p(indzs(ig))
                   C2_vpsi(ig,is1)=-fi*&
                        (0.5_real_8*parm%tpiba2*hgkp(ig,ikind)*c0(ig,is1)+fp)
                   C2_vpsi(ig+ncpw%ngw,is1)=-fi*&
                        (0.5_real_8*parm%tpiba2*hgkm(ig,ikind)*c0(ig+ncpw%ngw,is1)+fm)
                   C2_vpsi(ig,is1)=c2_vpsi(ig,is1)&
                        -fi*gpot*(rk(td_prop%pertdir,1)&
                        +gk(td_prop%pertdir,ig))*parm%tpiba*c0(ig ,is1)&
                        -fi*0.5_real_8*(gpot**2)*              c0(ig,is1)
                   C2_vpsi(ig+ncpw%ngw,is1)=c2_vpsi(ig+ncpw%ngw,is1)&
                        -fi*gpot*(rk(td_prop%pertdir,1)&
                        -gk(td_prop%pertdir,ig))*parm%tpiba*c0(ig+ncpw%ngw,is1)&
                        -fi*0.5_real_8*(gpot**2)*              c0(ig+ncpw%ngw,is1)
                ENDDO
             ELSE
                DO ig=1,ncpw%ngw
                   fp=psi_p(nzhs(ig))
                   fm=psi_p(indzs(ig))
                   C2_vpsi(ig,is1)=-fi*&
                        (0.5_real_8*parm%tpiba2*hgkp(ig,ikind)*c0(ig,is1)+fp)
                   C2_vpsi(ig+ncpw%ngw,is1)=-fi*&
                        (0.5_real_8*parm%tpiba2*hgkm(ig,ikind)*c0(ig+ncpw%ngw,is1)+fm)
                   ! elisa
                   DO id=1,3
                      IF (gampl(id).NE.0._real_8) THEN
                         C2_vpsi(ig,is1)=c2_vpsi(ig,is1)&
                              -fi*gpotv(id)*(rk(id,1)&
                              +gk(id,ig))*parm%tpiba*c0(ig ,is1)&
                              -fi*0.5_real_8*(gpotv(id)**2)*   c0(ig,is1)
                         C2_vpsi(ig+ncpw%ngw,is1)=c2_vpsi(ig+ncpw%ngw,is1)&
                              -fi*gpotv(id)*(rk(id,1)&
                              -gk(id,ig))*parm%tpiba*c0(ig+ncpw%ngw,is1)&
                              -fi*0.5_real_8*(gpotv(id)**2)*   c0(ig+ncpw%ngw,is1)
                      ENDIF
                   ENDDO
                   ! elisa
                ENDDO
             ENDIF

          ELSE
             !$omp parallel do private(IG,FP,FM) &
             !$omp   shared(ncpw,IS1,FI,IKIND,PSI_p,C0,C2_vpsi)
             DO ig=1,ncpw%ngw
                fp=psi_p(nzhs(ig))
                fm=psi_p(indzs(ig))
                C2_vpsi(ig,    is1)=-fi*&
                     (0.5_real_8*parm%tpiba2*hgkp(ig,ikind)*c0(ig,    is1)+fp)
                C2_vpsi(ig+ncpw%ngw,is1)=-fi*&
                     (0.5_real_8*parm%tpiba2*hgkm(ig,ikind)*c0(ig+ncpw%ngw,is1)+fm)
             ENDDO
          ENDIF
          ! EHR]
          IF (geq0)C2_vpsi(1+ncpw%ngw,is1)=CMPLX(0._real_8,0._real_8,kind=real_8)
       ELSE
          fi=f(is1)*0.5_real_8
          IF (fi.EQ.0._real_8.AND..NOT.cntl%tksham) fi=1._real_8
          IF (fi.EQ.0._real_8.AND.cntl%tksham) fi=0.5_real_8
          fip1=0._real_8
          IF (is2.LE.nostat) fip1=f(is2)*0.5_real_8
          IF (fip1.EQ.0._real_8.AND..NOT.cntl%tksham) fip1=1._real_8
          IF (fip1.EQ.0._real_8.AND.cntl%tksham) fip1=0.5_real_8
          IF (prcp_com%akin.GT.1.e-10_real_8) THEN
             xskin=1._real_8/prcp_com%gskin
             !$omp parallel do private(IG,FP,FM,G2,psin,psii) shared(JGW,PSI_p) &
             !$omp   shared(XSKIN,FI,C0,C2_vpsi,IS1) &
             !$omp   shared(NOSTAT,FIP1,IS2)
             !dir$ concurrent
             DO ig=1,jgw
                psin=psi_p(nzhs(ig))! access only once
                psii=psi_p(indzs(ig))! these mem locations
                fp=psin+psii
                fm=psin-psii
                g2=parm%tpiba2*(hg(ig)+&
                     prcp_com%gakin*(1._real_8+cp_erf((hg(ig)-prcp_com%gckin)*xskin)))
                C2_vpsi(ig,is1)=-fi*(g2*c0(ig,is1)+&
                     CMPLX(REAL(fp),AIMAG(fm),kind=real_8))
                IF (is2.LE.nostat) C2_vpsi(ig,is2)=-fip1*(g2*c0(ig,is2)+&
                     CMPLX(AIMAG(fp),-REAL(fm),kind=real_8))
             ENDDO
          ELSE
             !$omp parallel do private(IG,FP,FM,psin,psii) shared(JGW,PSI_p) &
             !$omp    shared(C2_vpsi,C0,FI,IS1,IS2,NOSTAT,FIP1)
             DO ig=1,jgw
                psin=psi_p(nzhs(ig))! access only once
                psii=psi_p(indzs(ig))! these mem locations
                fp=psin+psii
                fm=psin-psii
                C2_vpsi(ig,is1)=-fi*((parm%tpiba2*hg(ig))*c0(ig,is1)+&
                     CMPLX(REAL(fp),AIMAG(fm),kind=real_8))
                IF (is2.LE.nostat) C2_vpsi(ig,is2)=-fip1*&
                     ((parm%tpiba2*hg(ig))*&
                     c0(ig,is2)+CMPLX(AIMAG(fp),-REAL(fm),kind=real_8))
             ENDDO
          ENDIF
       ENDIF
       __NVTX_TIMER_STOP
    ENDDO                     ! Loop over IA
    !$omp end do

    !vw free local comm for each stream/thread
    CALL mp_comm_free ( fft_comm )

    !$omp end parallel

    __NVTX_TIMER_START ( procedureN//'_c' )

    CALL thread_views_finalize ( thread_views )

    IF( cp_cuda_env%use_fft ) THEN
       CALL tiset(procedureN//'_gpu_fin',isub4)
       CALL cuda_dealloc_host(psis)
#ifndef _HASNT_OMP_SET_NESTED
       !$ CALL omp_set_nested( nested_orig )
#endif
       !$ CALL omp_set_max_active_levels( max_level_orig )
       IF( .NOT. copy_data_to_host ) THEN
          ! TODO acm: the init+allocation/deallocation should be done outside (once for entire code)
          CALL cp_cuvpsi_finalize_and_dealloc( cp_cuvpsi )
       ENDIF
       CALL tihalt(procedureN//'_gpu_fin',isub4)
    ENDIF

    ! Special terms for MGGA. Parallelised over states, call needs
    ! to be before redist & sum to avoid useless allocation calls
    !
    IF (cntl%ttau) CALL vtaupsi(c0,c2_vpsi,f,psi,nstate,ispin)
    !
    ! redistribute C2_vpsi over the groups if needed
    !
    IF (redist_c2) THEN
       CALL tiset(procedureN//'_grps_b',isub3)
       CALL cp_grp_redist(C2_vpsi,nkpt%ngwk,nstate)
       CALL tihalt(procedureN//'_grps_b',isub3)
    ENDIF

    !
    ! add up the vpsi contribution to C2
    !
    CALL add_wfn(nkpt%ngwk,nstate,zone,C2_vpsi,nkpt%ngwk,c2,nkpt%ngwk)

    DEALLOCATE(C2_vpsi,stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',&
         __LINE__,__FILE__)

    IF (tkpts%tkpnt) CALL c_clean(c2,nstate,ikind)
    ! SPECIAL TERMS FOR LSE METHODS
    IF (lspin2%tlse) CALL vpsi_lse(c0,c2,f,vpot,psi,nstate,.TRUE.)
    !
    __NVTX_TIMER_STOP
    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vpsi

  ! ==================================================================
  SUBROUTINE vpsimt(c0,c2,f,vpot,psi,nstate,ispin,swfun)
    ! ==================================================================
    ! == No K-Points and no LSE  VERSION OF VPSI.                     ==
    ! == EXACTLY LIKE VPSI BUT NO CONTRIBUTION FROM KINETIC ENERGY    ==
    ! ==--------------------------------------------------------------==
    ! == VPOT:   IN INPUT POTENTIAL                                   ==
    ! == ISPIN:  Need with LSD option for diagonalization scheme      ==
    ! ==         dimension of VPOT(NNR1,ISPIN)                        ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:)
    REAL(real_8)                             :: f(:)
    REAL(real_8), TARGET                     :: vpot(:,:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate, ispin
    LOGICAL                                  :: swfun

    CHARACTER(*), PARAMETER                  :: procedureN = 'vpsimt'

    COMPLEX(real_8)                          :: c0tmp1, c0tmp2, fm, fp, psii, &
                                                psin
    INTEGER                                  :: i, ia, ib, ibb, ifft, ig, &
                                                is1, is2, isub, l, lead, &
                                                leadx, lspin, njump, nnrx, &
                                                nostat, nsta
    REAL(real_8)                             :: fi, fip1

!(nkpt%ngwk,nstate)
! Variables
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    ! Other illegal combinations already caught in VPSI
    IF (lspin2%tlse.AND.(lspin2%tcas22.OR.lspin2%tpenal)) THEN
       CALL stopgm('VPSIMT','LSE (CAS22, PENALTY) NOT IMPLEMENTED',&
            __LINE__,__FILE__)
    ENDIF
    IF (tkpts%tkpnt) CALL stopgm('VPSIMT','K-POINTS NOT IMPLEMENTED',&
         __LINE__,__FILE__)
    IF (group%nogrp.EQ.1) THEN
       lead  = fpar%kr1s*parai%ngrays
       leadx = fpar%nnr1
       ifft  = 1
       nnrx=fpar%nnr1
    ELSE
       CALL stopgm(procedureN,'this doesnt work anymore....',&
            __LINE__,__FILE__)
       lead  = fpar%kr1s*ngrm
       leadx = fpar%krx*fpar%kr2s*fpar%kr3s
       ifft=2
       nnrx=fpar%krx*fpar%kr2s*fpar%kr3s
    ENDIF
    njump=2*group%nogrp
    IF (lspin2%tlse) THEN
       nostat = clsd%ialpha-1
    ELSE
       nostat = nstate
    ENDIF
    DO ia=1,nostat,njump
       CALL zeroing(psi)!,maxfft*group%nogrp)
       ! ==----------------------------------------------------------==
       ! == Store the wavefunctions in array PSI                     ==
       ! == in the following way:                                    ==
       ! == Since we use the Gamma point only, the wavefunctions in  ==
       ! == real space are real. This fact is used to save time by   ==
       ! == using complex arrays where the wavefunctions of          ==
       ! == 2 different states (i and i+1) are combined, and         ==
       ! == performing complex Fourier Transforms.                   ==
       ! == Warning: If you use K-points other than Gamma, or zone   ==
       ! == boundary, this symmetry is broken and the trick is not   ==
       ! == valid anymore. First step is to build the array PSI.     ==
       ! ==----------------------------------------------------------==
       nsta=MIN((nostat-ia+2)/2,group%nogrp)
       DO ib=1,nsta
          i=ia+2*(ib-1)
          ibb=(ib-1)*lead
          is1 = i
          is2 = i+1
          IF (is2.GT.nostat) THEN
             !CDIR NODEP
#if defined(__VECTOR)
             !$omp parallel do private(IG)
             DO ig=1,ncpw%ngw
                psi(nzhs(ig)+ibb)=c0(ig,is1)
                psi(indzs(ig)+ibb)=CONJG(c0(ig,is1))
             ENDDO
#else
             !$omp parallel do private(IG,C0TMP1) schedule(static)
             DO ig=1,ncpw%ngw
                c0tmp1=c0(ig,is1)
                psi(nzhs(ig)+ibb)=c0tmp1
                psi(indzs(ig)+ibb)=CONJG(c0tmp1)
             ENDDO
#endif
             IF (geq0) psi(nzhs(1)+ibb)=c0(1,is1)
          ELSE
             !CDIR NODEP
#if defined(__VECTOR)
             !$omp parallel do private(IG)
             DO ig=1,ncpw%ngw
                psi(nzhs(ig)+ibb)=c0(ig,is1)+uimag*c0(ig,is2)
                psi(indzs(ig)+ibb)=CONJG(c0(ig,is1))+&
                     uimag*CONJG(c0(ig,is2))
             ENDDO
#else
             !$omp parallel do private(IG,C0TMP1,C0TMP2) schedule(static)
             DO ig=1,ncpw%ngw
                c0tmp1=c0(ig,is1)
                c0tmp2=c0(ig,is2)
                psi(nzhs(ig)+ibb)=c0tmp1+uimag*c0tmp2
                psi(indzs(ig)+ibb)=CONJG(c0tmp1)+uimag*CONJG(c0tmp2)
             ENDDO
#endif
             IF (geq0) psi(nzhs(1)+ibb)=c0(1,is1)+uimag*c0(1,is2)
          ENDIF
       ENDDO
       ! ==----------------------------------------------------------==
       ! == Transform the wavefunction to real space                 ==
       ! ==----------------------------------------------------------==
       IF (ifft.EQ.1) THEN
          CALL  invfftn(psi,.TRUE.,parai%allgrp)
       ELSE
          CALL stopgm("THIS","SG_INVFFT NOT AVAILABLE ANYMORE",&
               __LINE__,__FILE__)
       ENDIF
       ! ==------------------------------------------------------------==
       ! == Apply the potential (V), which acts in real space.         ==
       ! ==------------------------------------------------------------==
       IF (cntl%tlsd.AND.ispin.EQ.2) THEN
          i=ia
          IF (i.EQ.spin_mod%nsup) THEN
             !$omp parallel do private(L) schedule(static)
             DO l=1,nnrx
                psi(l)=vpot(l,1)*REAL(psi(l))&
                     +uimag*vpot(l,2)*AIMAG(psi(l))
             ENDDO
          ELSE
             lspin=1
             IF (i.GT.spin_mod%nsup) lspin=2
             !$omp parallel do private(L) schedule(static)
             DO l=1,nnrx
                psi(l)=vpot(l,lspin)*psi(l)
             ENDDO
          ENDIF
       ELSE
          !$omp parallel do private(L) schedule(static)
          DO l=1,nnrx
             psi(l)=vpot(l,1)*psi(l)
          ENDDO
       ENDIF
       ! ==------------------------------------------------------------==
       ! == Back transform to reciprocal space the product V.PSI       ==
       ! ==------------------------------------------------------------==
       IF (ifft.EQ.1) THEN
          CALL  fwfftn(psi,.TRUE.,parai%allgrp)
       ELSE
          CALL stopgm("THIS","SG_FWFFT NOT AVAILABLE ANYMORE",&
               __LINE__,__FILE__)
       ENDIF
       ! ==------------------------------------------------------------==
       ! == Decode the combination of wavefunctions (now multiplied by ==
       ! == the potential), back to those of states i and i+1. Multiply==
       ! == also by the occupation numbers f(i) and other constants in ==
       ! == order to obtain the force on the electronic degrees of     ==
       ! == freedom, stored in array FC.                               ==
       ! == Here we also add the force from the kinetic energy.        ==
       ! ==------------------------------------------------------------==
       DO ib=1,nsta
          i=ia+2*(ib-1)
          ibb=(ib-1)*lead
          is1 = i
          is2 = i+1
          fi=f(is1)*0.5_real_8
          IF (fi.EQ.0._real_8) fi=1._real_8
          fip1=0._real_8
          IF (is2.LE.nostat) fip1=f(is2)*0.5_real_8
          IF (fip1.EQ.0._real_8) fip1=1._real_8
          !$omp parallel do private(IG,psin,psii,FP,FM) schedule(static)
          DO ig=1,ncpw%ngw
             psin=psi(nzhs(ig)+ibb)
             psii=psi(indzs(ig)+ibb)
             fp=psin+psii
             fm=psin-psii
             ! mb            FP=PSI(NZHS(IG)+IBB)+PSI(INDZS(IG)+IBB)
             ! mb            FM=PSI(NZHS(IG)+IBB)-PSI(INDZS(IG)+IBB)
             c2(ig,is1)=c2(ig,is1)-fi*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
             IF (is2.LE.nostat)&
                  c2(ig,is2)=c2(ig,is2)-fip1*CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
          ENDDO
       ENDDO
    ENDDO                     ! Loop over IA
    ! Reestablish initial conditions for the potential
    ! SPECIAL TERMS FOR LSE METHODS
    IF (lspin2%tlse) THEN
       CALL vpsi_lse(c0,c2,f,vpot,psi,nstate,.FALSE.)
    ENDIF
    ! META FUNCTIONALS NEED SPECIAL TERM
    IF (cntl%ttau) THEN
       IF (swfun) CALL switch_functionals
       IF (cntl%ttau) CALL vtaupsi(c0,c2,f,psi,nstate,ispin)
       IF (swfun) CALL switch_functionals
    ENDIF
    !
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vpsimt
  ! ==================================================================
  SUBROUTINE movepsid(a)
    COMPLEX(real_8)                          :: a(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'movepsid'

    COMPLEX(real_8), ALLOCATABLE             :: b(:)
    INTEGER                                  :: ierr, ig

! ==--------------------------------------------------------------==

    ALLOCATE(b(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    !$omp parallel do private(IG) schedule(static)
    DO ig=1,ncpw%nhg
       b(ig)=a(nzh(ig))
    ENDDO
    CALL zeroing(a)!,maxfftn)
    !$omp parallel do private(IG) schedule(static)
    DO ig=1,jhg
       a(nzff(ig))=b(ig)
       a(inzf(ig))=CONJG(b(ig))
    ENDDO
    IF (geq0) a(nzff(1))=b(1)
    DEALLOCATE(b,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE movepsid
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE vpsi_batchfft(c0,c2,f,vpot,psi,nstate,ikind,ispin,redist_c2)
    ! ==================================================================
    ! == K-POINT AND NOT K-POINT VERSION OF VPSI.                     ==
    ! ==--------------------------------------------------------------==
    ! == VPOT:   IN INPUT POTENTIAL                                   ==
    ! == ISPIN:  Need with LSD option for diagonalization scheme      ==
    ! ==         dimension of VPOT(NNR1,ISPIN)                        ==
    ! ==--------------------------------------------------------------==
    ! EHR[
    ! EHR]
    ! Modified: Tobias Kloeffel, Erlangen
    ! Date May 2019
    ! special version of vpsi to use the batch fft driver
    ! TODO
    ! move communication phase into vpsi:
    ! benefits: reduces memory footprint as only two batches are needed
    ! in memory; expands the time for the communication phase as also
    ! the decobination phase of the wf's can take place during
    ! communication phse
    ! cons: code complexity will increase, e.g. calling alltoall from here?
    ! Full performance only with saved arrays or scratch_library

    COMPLEX(real_8) __CONTIGUOUS             :: c0(:,:), c2(:,:)
    REAL(real_8) __CONTIGUOUS                :: f(:)
    REAL(real_8), TARGET __CONTIGUOUS        :: vpot(:,:)
    COMPLEX(real_8), TARGET __CONTIGUOUS     :: psi(:)
    INTEGER                                  :: nstate, ikind, ispin
    LOGICAL                                  :: redist_c2
    LOGICAL                                  :: lg_vpotx3a, lg_vpotx3b
    CHARACTER(*), PARAMETER                  :: procedureN = 'vpsi_batchfft'
    COMPLEX(real_8), PARAMETER               :: zone = (1.0_real_8,0.0_real_8)

    COMPLEX(real_8)                          :: fm, fp, psii, psin
    COMPLEX(real_8), POINTER __CONTIGUOUS &
                           , ASYNCHRONOUS    :: wfn_r1(:)
    INTEGER :: i, iclpot = 0, id, ierr, ig, &
      ir, is1, is2, isub, isub2, isub3, isub4, iwf, ixx, ixxs, iyy, izz, jj, &
      leadx, njump, nnrx, nostat, nrxyz1s, nrxyz2, start_loop2, &
      ist,states_fft,  bsize, ibatch, istate, ir1, first_state, end_loop2, &
      i_start1, i_start2, i_start3, me_grp, n_grp, start_loop1, end_loop1,&
      offset_state, nthreads, nested_threads, methread, count, swap, int_mod
    INTEGER(int_8)                           :: il_wfng(2), il_wfnr(2), il_wfnr1(1), il_xf(2)
    REAL(real_8)                             :: chksum, csmult, fi, fip1,&
                                                xskin, temp_time
    REAL(real_8), ALLOCATABLE                :: vpotx3a(:,:,:), vpotx3b(:,:,:)
    REAL(real_8), POINTER __CONTIGUOUS       :: VPOTX(:),vpotdg(:,:,:),extf_p(:,:)
    INTEGER, ALLOCATABLE                     :: lspin(:,:)
 
    INTEGER                                  :: rem, i1, j1
    INTEGER(INT64) :: time(2)
    ! ==--------------------------------------------------------------==

    IF(cntl%fft_tune_batchsize) THEN
       CALL tiset(procedureN//'_tuning',isub4)
    ELSE
       CALL tiset(procedureN,isub)
    END IF

    IF (group%nogrp.GT.1)CALL stopgm(procedureN,&
         'OLD TASK GROUPS NOT SUPPORTED ANYMORE ',&
         __LINE__,__FILE__)
    IF (tdgcomm%tdg) &
         CALL stopgm(procedureN,'DOUBLE GRID NOT SUPPORTED', &
               __LINE__,__FILE__)
    CALL setfftn(0)

    IF (lspin2%tlse.AND.lspin2%tlsets.AND.(lspin2%tross.OR.lspin2%tcas22.OR.lspin2%tpenal&
         .OR.lspin2%troot)) THEN
       CALL stopgm(procedureN,&
            'NO SLATER TS WITH ROSS, CAS22, PENALTY, ROOTHAAN',&
            __LINE__,__FILE__)
    ENDIF

    IF(td_prop%td_extpot) CALL reshape_inplace(extf, (/fpar%kr1*fpar%kr2s,fpar%kr3s/), &
         extf_p)
    CALL reshape_inplace(vpot, (/fpar%kr1*fpar%kr2s,fpar%kr3s,ispin/), vpotdg)
    !
    !vpotx(1:SIZE(vpotdg)) => vpotdg
    CALL reshape_inplace(vpot, (/fpar%nnr1*ispin/), vpotx)

    ALLOCATE(lspin(2,fft_batchsize),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate lspin', &
         __LINE__,__FILE__)

    !$ ALLOCATE(locks_fw(fft_numbatches+1,2),STAT=ierr)
    !$ IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate locks_fw', &
    !$      __LINE__,__FILE__)
    !$ ALLOCATE(locks_inv(fft_numbatches+1,2),STAT=ierr)
    !$ IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate locks_inv', &
    !$      __LINE__,__FILE__)

    !if overlapping comm/comp active:
    ! if rsactive we operate on two batches
    ! else we operate on 3 batches
    !else just one batch
    int_mod=3
    il_wfng(1)=fpar%kr1s*msrays
    il_wfng(2)=fft_batchsize

    il_wfnr1(1)=fpar%kr1*fpar%kr2s*fpar%kr3s*fft_batchsize

    il_wfnr(1)=fpar%kr1*fpar%kr2s*fpar%kr3s*fft_batchsize
    IF(il_wfnr(1).EQ.0)il_wfnr(1)=fpar%kr2s*fpar%kr3s*fft_batchsize
    il_wfnr(1)=il_wfnr(1)+MOD(il_wfnr(1),4)
    il_wfnr(2)=1

    il_xf(1)=fpar%nnr1*fft_batchsize
    IF(il_xf(1).EQ.0) il_xf(1)=maxfft*fft_batchsize
    il_xf(1)=il_xf(1)+MOD(il_xf(1),4)
    il_xf(2)=3
    start_loop1=1
    end_loop1=fft_numbatches+2
    start_loop2=2
    end_loop2=fft_numbatches+3

    IF(rsactive)THEN
       il_wfnr(2)=fft_numbatches
       IF(fft_residual.GT.0) il_wfnr(2)=fft_numbatches+1
       il_xf(2)=2
       int_mod=2
       start_loop1=0
       end_loop1=fft_numbatches+1
       start_loop2=1
       end_loop2=fft_numbatches+2
    END IF

    IF(cntl%overlapp_comm_comp.AND.fft_numbatches.GT.1)THEN
       nthreads=MIN(2,parai%ncpus)
       nested_threads=(MAX(parai%ncpus-1,1))
#if !defined(_INTEL_MKL)
       CALL stopgm(procedureN, 'Overlapping communication and computation: Behavior of BLAS &
            routine inside parallel region not checked',&
            __LINE__,__FILE__)
#endif
    ELSE
       nthreads=1
       nested_threads=parai%ncpus
       int_mod=1
       il_xf(2)=1
       start_loop1=0
       start_loop2=0
       end_loop1=fft_numbatches+1
       end_loop2=fft_numbatches+1
    END IF
 
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_xf,xf,procedureN//'_xf',ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate xf', &
         __LINE__,__FILE__)
    CALL request_scratch(il_xf,yf,procedureN//'_yf',ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate yf', &
         __LINE__,__FILE__)
    CALL request_scratch(il_wfng,wfn_g,'wfn_g',ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate wfn_g', &
         __LINE__,__FILE__)
#endif
    IF(.NOT.rsactive)THEN
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_wfnr,wfn_r,'wfn_r',ierr)
#else
       ALLOCATE(wfn_r(il_wfnr(1),il_wfnr(2)),STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate wfn_r', &
            __LINE__,__FILE__)
#if !defined _USE_SCRATCHLIBRARY
       ALLOCATE(wfn_g(il_wfng(1),il_wfng(2)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate wfn_g', &
            __LINE__,__FILE__)
#endif
    END IF
    ! 

    leadx = fpar%nnr1
    nnrx  = llr1

    njump=2
    IF (tkpts%tkpnt) njump=1
    IF (lspin2%tlse) THEN
       nostat = clsd%ialpha-1
    ELSE
       nostat = nstate
    ENDIF
    me_grp=parai%cp_inter_me
    n_grp=parai%cp_nogrp
    i_start1=0
    i_start2=part_1d_get_el_in_blk(1,nostat,me_grp,n_grp)-1
    i_start3=part_1d_get_el_in_blk(1,nostat,me_grp,n_grp)-1


    IF (cntl%cdft)THEN
       IF (.NOT.cdftlog%tcall)THEN
          csmult=cdftcom%cdft_v(1)
       ELSE
          csmult=cdftcom%cdft_v(1)+cdftcom%cdft_v(2)
       ENDIF
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          vpotdg(ir,1,1)=vpotdg(ir,1,1)+csmult*wdiff(ir)
       ENDDO
       IF (.NOT.cdftlog%tcall)THEN
          IF (cdftlog%tspinc) THEN
             csmult=-cdftcom%cdft_v(1)
          ELSE
             csmult=cdftcom%cdft_v(1)
          ENDIF
       ELSE
          csmult=cdftcom%cdft_v(1)-cdftcom%cdft_v(2)
       ENDIF
       IF (cntl%tlsd.AND.ispin.EQ.2) THEN
          !$omp parallel do private(IR)
          DO ir=1,fpar%nnr1
             vpotdg(ir,2,1)=vpotdg(ir,2,1)+csmult*wdiff(ir)
          ENDDO
       ENDIF
    ENDIF
    CALL reshape_inplace(vpot, (/fpar%kr1*fpar%kr2s,fpar%kr3s,ispin/), vpotdg)

    methread=0


    IF(.NOT.rsactive) wfn_r1=>wfn_r(:,1)
    IF(cntl%fft_tune_batchsize) temp_time=m_walltime()
    CALL SYSTEM_CLOCK( time(1) )
    !$ locks_inv = .TRUE.
    !$ locks_fw = .TRUE.
    !$OMP parallel IF(nthreads.eq.2) num_threads(nthreads) &
    !$omp private(methread,ibatch,bsize,count,ist,is1,is2,ir,offset_state,swap) &
    !$omp proc_bind(close)
    !$ methread = omp_get_thread_num()
    !$ IF(methread.EQ.1)THEN
    !$    CALL omp_set_max_active_levels(2)
    !$    CALL omp_set_num_threads(nested_threads)
#if defined(_HAS_FFT_FFTW3)
    !$    CALL dfftw_plan_with_nthreads(nested_threads)
#endif
#ifdef _INTEL_MKL
    !$    CALL mkl_set_dynamic(0)
    !$    ierr = mkl_set_num_threads_local(nested_threads)
#endif
    !$ END IF
    !$OMP barrier
    
    !Loop over batches
    DO ibatch=1,fft_numbatches+3
       IF(.NOT.rsactive)THEN
          IF(methread.EQ.1.OR.nthreads.EQ.1)THEN
             !process batches starting from ibatch .eq. 1 until ibatch .eq. fft_numbatches+1
             IF(ibatch.LE.fft_numbatches+1)THEN
                IF(ibatch.LE.fft_numbatches)THEN
                   bsize=fft_batchsize
                ELSE
                   bsize=fft_residual
                END IF
                IF(bsize.NE.0)THEN
                   ! Loop over the electronic states of this batch
                   CALL set_psi_batch_g(c0,wfn_g,int(il_wfng(1)),i_start1,bsize,nstate,me_grp,n_grp)
                   ! ==--------------------------------------------------------------==
                   ! ==  Fourier transform the wave functions to real space.         ==
                   ! ==  In the array PSI was used also the fact that the wave       ==
                   ! ==  functions at Gamma are real, to form a complex array (PSI)  ==
                   ! ==  with the wave functions corresponding to two different      ==
                   ! ==  states (i and i+1) as the real and imaginary part. This     ==
                   ! ==  allows to call the FFT routine 1/2 of the times and save    ==
                   ! ==  time.                                                       ==
                   ! ==  Here we operate on a batch of states, containing njump*bsize==
                   ! ==  states. To achive better overlapping of the communication   ==
                   ! ==  and communication phase, we operate on two batches at once  ==
                   ! ==  ist revers to the current batch (rsactive) or is identical  ==
                   ! ==  to swap                                                     ==
                   ! ==--------------------------------------------------------------==
                   swap=mod(ibatch,int_mod)+1
                   CALL invfftn_batch(wfn_g,size(wfn_g),bsize,swap,1,ibatch)
                   i_start1=i_start1+bsize*njump
                END IF
             END IF
          END IF
          IF(methread.EQ.0.OR.nthreads.EQ.1)THEN
             !process batches starting from ibatch .eq. 1 until ibatch .eq. fft_numbatches+1
             !communication phase
             IF(ibatch.LE.fft_numbatches+1)THEN
                IF(ibatch.LE.fft_numbatches)THEN
                   bsize=fft_batchsize
                ELSE
                   bsize=fft_residual
                END IF
                IF(bsize.NE.0)THEN
                   swap=mod(ibatch,int_mod)+1
                   CALL invfftn_batch(wfn_r,int(il_wfnr(1)),bsize,swap,2,ibatch)
                END IF
             END IF
          END IF
          IF (methread.EQ.1.OR.nthreads.EQ.1)THEN
             !process batches starting from ibatch .eq. 2 until ibatch .eq. fft_numbatches+2
             !data related to ibatch-1!
!             IF(ibatch.GE.2.AND.ibatch.LE.fft_numbatches+2)THEN
             IF(ibatch.GT.start_loop1.AND.ibatch.LE.end_loop1)THEN
                IF (ibatch-start_loop1.LE.fft_numbatches)THEN
                   bsize=fft_batchsize
                ELSE
                   bsize=fft_residual
                END IF
                IF(bsize.NE.0)THEN
                   swap=mod(ibatch-start_loop1,int_mod)+1
                   CALL invfftn_batch(wfn_r1,int(il_wfnr(1)),bsize,swap,3,ibatch-start_loop1)
                END IF
             END IF
          END IF
          ! At this point locks_inv(:,ibatch-1) are .FALSE.
       END IF
       ! ==------------------------------------------------------------==
       ! == Apply the potential (V), which acts in real space.         ==

       IF(methread.EQ.1.OR.nthreads.EQ.1)THEN
!          IF(ibatch.GE.2.AND.ibatch.LE.fft_numbatches+2)THEN
          IF(ibatch.GT.start_loop1.AND.ibatch.LE.end_loop1)THEN
             IF(ibatch-start_loop1.LE.fft_numbatches)THEN
                bsize=fft_batchsize
             ELSE
                bsize=fft_residual
             END IF
             IF(bsize.NE.0)THEN
                swap=mod(ibatch-start_loop1,int_mod)+1
                lspin=1
                offset_state=i_start2
                DO count=1,bsize
                   is1=offset_state+1
                   offset_state=offset_state+1
                   IF (njump.EQ.2) THEN
                      is2=offset_state+1
                      offset_state=offset_state+1
                   ELSE
                      is2=0
                   END IF
                   IF (cntl%tlsd.AND.ispin.EQ.2) THEN
                      IF (is1.GT.spin_mod%nsup) lspin(1,count)=2
                      IF (is2.GT.spin_mod%nsup) lspin(2,count)=2
                   END IF
                   !njump states per single fft
                END DO
!                IF(rsactive)THEN
!                   CALL mult_vpot_psi_rsactive(wfn_r(:,ibatch-1),wfn_r1,vpotdg,&
!                        fpar%kr1*fpar%kr2s,bsize,fpar%kr3s,lspin,clsd%nlsd)
!                ELSE
                IF(rsactive) wfn_r1=>wfn_r(:,ibatch-start_loop1)
                   CALL mult_vpot_psi(wfn_r1,vpotdg,&
                        fpar%kr1*fpar%kr2s,bsize,fpar%kr3s,lspin,clsd%nlsd)
!                END IF
                IF (td_prop%td_extpot.AND.cntl%tlsd.AND.ispin.EQ.2) THEN
                   offset_state=i_start2
                   DO count=1,bsize
                      is1=offset_state+1
                      offset_state=offset_state+1
                      is2=0
                      IF(njump.EQ.2) THEN
                         is2=offset_state+1
                         offset_state=offset_state+1
                      END IF
                      IF (is1.EQ.spin_mod%nsup) EXIT
                   END DO
                   !count is now set to: count .eq. spin_mod%nsup or count .eq. bsize+1
                   IF (count .LE. bsize) THEN
                      CALL mult_extf_psi(wfn_r1,extf,fpar%kr1*fpar%kr2s,bsize,&
                           fpar%kr3s,count)
                   END IF
                END IF
             ! ==------------------------------------------------------------==
             ! == Back transform to reciprocal space the product V.PSI       ==
             ! ==------------------------------------------------------------==
                CALL fwfftn_batch(wfn_r1,int(il_wfnr(1)),bsize,swap,1,ibatch-start_loop1)
                i_start2=i_start2+bsize*njump
             END IF
          END IF
       END IF
       IF(methread.EQ.0.OR.nthreads.EQ.1)THEN
!          IF(ibatch.GE.2.AND.ibatch.LE.fft_numbatches+2)THEN
          IF(ibatch.GT.start_loop1.AND.ibatch.LE.end_loop1)THEN
             IF(ibatch-start_loop1.LE.fft_numbatches)THEN
                bsize=fft_batchsize
             ELSE
                bsize=fft_residual
             END IF
             IF(bsize.NE.0)THEN
                swap=mod(ibatch-start_loop1,int_mod)+1
                CALL fwfftn_batch(wfn_r,int(il_wfnr(1)),bsize,swap,2,ibatch-start_loop1)
             END IF
          END IF
       END IF
       IF(methread.EQ.1.OR.nthreads.EQ.1)THEN
          !data related to ibatch-2
!          IF(ibatch.GE.3.AND.ibatch.LE.fft_numbatches+3)THEN
          IF(ibatch.GT.start_loop2.AND.ibatch.LE.end_loop2)THEN
!             IF(ibatch-2.LE.fft_numbatches)THEN
             IF(ibatch-start_loop2.LE.fft_numbatches)THEN
                bsize=fft_batchsize
             ELSE
                bsize=fft_residual
             END IF
             IF(bsize.NE.0)THEN
                IF(rsactive) wfn_r1=>wfn_r(:,ibatch-start_loop2)
                swap=mod(ibatch-start_loop2,int_mod)+1
                CALL fwfftn_batch(wfn_g,size(wfn_g),bsize,swap,3,ibatch-start_loop2)
                IF (tkpts%tkpnt) THEN
                   CALL calc_c2_kpnt(c2,c0,wfn_g,f,bsize,i_start3,ikind)
                ELSE
                   CALL calc_c2(c2,c0,wfn_g,f,bsize,i_start3,njump,nostat)
                ENDIF
                i_start3=i_start3+bsize*njump
             END IF
          END IF
       END IF
    END DO

    !$ IF (methread.EQ.1) THEN
    !$    CALL omp_set_max_active_levels(1)
    !$    CALL omp_set_num_threads(parai%ncpus)
#ifdef _INTEL_MKL
    !$    CALL mkl_set_dynamic(1)
    !$    ierr = mkl_set_num_threads_local(0)
#endif
#if defined(_HAS_FFT_FFTW3)
    !$    CALL dfftw_plan_with_nthreads(parai%ncpus)
#endif
    !$ END IF
    !$OMP barrier

    !$omp end parallel

    CALL SYSTEM_CLOCK( time(2) )
    dfft%time_adding( 99 ) = time(2) - time(1)

    IF(cntl%fft_tune_batchsize) fft_time_total(fft_tune_num_it)=fft_time_total(fft_tune_num_it)+m_walltime()-temp_time

    DO i=1,2
       ! kk-mb === print local potential (start) ===
       IF (locpot2%tlpot) THEN
          iclpot=iclpot+1
          nrxyz1s=spar%nr1s*spar%nr2s*spar%nr3s
          nrxyz2=0
          IF (iclpot .EQ. 1) THEN
             ALLOCATE(vpotx3a(spar%nr1s,spar%nr2s,spar%nr3s),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             CALL zeroing(vpotx3a)!,nrxyz1s)
             lg_vpotx3a=.TRUE.
             lg_vpotx3b=.FALSE.
          ELSE
             IF (lg_vpotx3a) THEN
                ALLOCATE(vpotx3b(spar%nr1s,spar%nr2s,spar%nr3s),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                CALL zeroing(vpotx3b)!,nrxyz1s)
                DEALLOCATE(vpotx3a,STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                     __LINE__,__FILE__)
                lg_vpotx3a=.FALSE.
                lg_vpotx3b=.TRUE.
             ELSE IF (lg_vpotx3b) THEN
                ALLOCATE(vpotx3a(spar%nr1s,spar%nr2s,spar%nr3s),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                CALL zeroing(vpotx3a)!,nrxyz1s)
                DEALLOCATE(vpotx3b,STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                     __LINE__,__FILE__)
                lg_vpotx3a=.TRUE.
                lg_vpotx3b=.FALSE.
             ENDIF
          ENDIF
          DO ir=1,nnrx
             IF (vpotx(ir) .NE. 0) THEN
                nrxyz2=nrxyz2+1
                ixx=MOD(nrxyz2-1,parm%nr1)+1
                jj=(nrxyz2-1)/parm%nr1+1
                iyy=MOD(jj-1,spar%nr2s)+1
                izz=(nrxyz2-1)/(parm%nr1*spar%nr2s)+1
                ixxs=ixx+parap%nrxpl(1,parai%mepos)-1
                IF (lg_vpotx3a) THEN
                   vpotx3a(ixxs,iyy,izz)=vpotx(ir)
                ELSE IF (lg_vpotx3b) THEN
                   vpotx3b(ixxs,iyy,izz)=vpotx(ir)
                ENDIF
             ENDIF
          ENDDO
       ENDIF
       ! kk-mb === print local potential (end)  ===
    END DO
    !
    ! redistribute C2 over the groups if needed
    !
    IF (redist_c2) THEN
       CALL tiset(procedureN//'_grps_b',isub3)
       CALL cp_grp_redist_array(C2,nkpt%ngwk,nstate)
       CALL tihalt(procedureN//'_grps_b',isub3)
    ENDIF

    !free wfn_g,and wfn_r
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_xf,yf,procedureN//'_yf',ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate yf', &
         __LINE__,__FILE__)
    CALL free_scratch(il_xf,xf,procedureN//'_xf',ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate xf', &
         __LINE__,__FILE__)
    CALL free_scratch(il_wfnr,wfn_r,'wfn_r',ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate wfn_r', &
         __LINE__,__FILE__)
    CALL free_scratch(il_wfng,wfn_g,'wfn_g',ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate wfn_g', &
         __LINE__,__FILE__)
#else
    IF(.NOT.rsactive)THEN
       DEALLOCATE(wfn_g,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate wfn_g', &
            __LINE__,__FILE__)
       DEALLOCATE(wfn_r,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate wfn_r', &
            __LINE__,__FILE__)
    END IF
#endif
    !$ DEALLOCATE(locks_fw,STAT=ierr)
    !$ IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate locks_fw', &
    !$      __LINE__,__FILE__)
    !$ DEALLOCATE(locks_inv,STAT=ierr)
    !$ IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate locks_inv', &
    !$      __LINE__,__FILE__)
    DEALLOCATE(lspin,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate lspin', &
         __LINE__,__FILE__)

    IF (tkpts%tkpnt) CALL c_clean(c2,nstate,ikind)
    ! SPECIAL TERMS FOR LSE METHODS
    IF (lspin2%tlse) CALL vpsi_lse(c0,c2,f,vpot,psi,nstate,.TRUE.)
    ! META FUNCTIONALS NEED SPECIAL TERM
    IF (cntl%ttau) CALL vtaupsi(c0,c2,f,psi,nstate,ispin)
    !
    IF(cntl%fft_tune_batchsize) THEN
       CALL tihalt(procedureN//'_tuning',isub4)
    ELSE
       CALL tihalt(procedureN,isub)
    END IF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vpsi_batchfft
  ! ==================================================================
  SUBROUTINE mult_vpot_psi(psi,vpot,n1,n2,n3,spins,nspin)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: n1, n2, n3, nspin, spins(2,n2)
    COMPLEX(real_8), INTENT(INOUT)           :: psi(n1,n2,n3)
    REAL(real_8), INTENT(IN)                 :: vpot(n1,n3,nspin)

    INTEGER                                  :: l1,l2,l3
    
    !$omp parallel do private(l1,l2,l3)
    DO l3=1,n3
       DO l2=1,n2
          DO l1=1,n1
             psi(l1,l2,l3)=&
                  vpot(l1,l3,spins(1,l2))*REAL(psi(l1,l2,l3),KIND=real_8)&
                  +uimag*vpot(l1,l3,spins(2,l2))*AIMAG(psi(l1,l2,l3))
          END DO
       END DO
    END DO
  END SUBROUTINE mult_vpot_psi
  ! ==================================================================
  SUBROUTINE mult_vpot_psi_rsactive(psi_in,psi_out,vpot,n1,n2,n3,spins,nspin)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: n1, n2, n3, nspin, spins(2,n2)
    COMPLEX(real_8), INTENT(OUT)             :: psi_out(n1,n2,n3)
    COMPLEX(real_8), INTENT(IN)              :: psi_in(n1,n2,n3)
    REAL(real_8), INTENT(IN)                 :: vpot(n1,n3,nspin)

    INTEGER                                  :: l1,l2,l3
    
    !$omp parallel do private(l1,l2,l3)
    DO l3=1,n3
       DO l2=1,n2
          DO l1=1,n1
             psi_out(l1,l2,l3)=&
                  vpot(l1,l3,spins(1,l2))*REAL(psi_in(l1,l2,l3),KIND=real_8)&
                  +uimag*vpot(l1,l3,spins(2,l2))*AIMAG(psi_in(l1,l2,l3))
          END DO
       END DO
    END DO
  END SUBROUTINE mult_vpot_psi_rsactive
  ! ==================================================================
  SUBROUTINE mult_extf_psi(psi,vpot,n1,n2,n3,l2)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: n1, n2, n3, l2
    COMPLEX(real_8), INTENT(INOUT)           :: psi(n1,n2,n3)
    REAL(real_8), INTENT(IN)                 :: vpot(n1,n3,2)

    INTEGER                                  :: l1,l3
    
    !$omp parallel do private(l3,l1)
    DO l3=1,n3
       DO l1=1,n1
          psi(l1,l2,l3)=&
               vpot(l1,l3,1)*REAL(psi(l1,l2,l3),KIND=real_8)&
               +uimag*vpot(l1,l3,2)*AIMAG(psi(l1,l2,l3))
       END DO
    END DO
  END SUBROUTINE mult_extf_psi
  ! ==================================================================
  SUBROUTINE calc_c2_kpnt(c2,c0,psi,f,n,is_start,ikind)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: n, is_start, ikind
    COMPLEX(real_8), INTENT(IN)              :: psi(fpar%kr1s*msrays,n), c0(:,:)
    COMPLEX(real_8), INTENT(OUT)             :: c2(:,:)
    REAL(real_8), INTENT(IN)                 :: f(:)
    COMPLEX(real_8)                          :: fm, fp
    INTEGER                                  :: is, is1, ig, id
    REAL(real_8)                             :: fi

    !$omp parallel private(is,is1,fi,ig,fp,fm)
    DO is=1,n
       is1=is+is_start
       fi=f(is1)
       IF (fi.EQ.0._real_8) fi=2._real_8
       ! EHR[
       IF (cntl%tgaugep.OR.cntl%tgaugef) THEN
          IF (gndir.EQ.0) THEN
             !$omp do
             DO ig=1,ncpw%ngw
                fp=psi(nzhs(ig),is)
                fm=psi(indzs(ig),is)
                C2(ig,is1)=-fi*&
                     (0.5_real_8*parm%tpiba2*hgkp(ig,ikind)*c0(ig,is1)+fp)
                C2(ig+ncpw%ngw,is1)=-fi*&
                     (0.5_real_8*parm%tpiba2*hgkm(ig,ikind)*c0(ig+ncpw%ngw,is1)+fm)
                C2(ig,is1)=c2(ig,is1)&
                     -fi*gpot*(rk(td_prop%pertdir,1)&
                     +gk(td_prop%pertdir,ig))*parm%tpiba*c0(ig ,is1)&
                     -fi*0.5_real_8*(gpot**2)*              c0(ig,is1)
                C2(ig+ncpw%ngw,is1)=c2(ig+ncpw%ngw,is1)&
                     -fi*gpot*(rk(td_prop%pertdir,1)&
                     -gk(td_prop%pertdir,ig))*parm%tpiba*c0(ig+ncpw%ngw,is1)&
                     -fi*0.5_real_8*(gpot**2)*              c0(ig+ncpw%ngw,is1)
             ENDDO
             !$omp end do nowait
          ELSE
             !$omp do
             DO ig=1,ncpw%ngw
                fp=psi(nzhs(ig),is)
                fm=psi(indzs(ig),is)
                C2(ig,is1)=-fi*&
                     (0.5_real_8*parm%tpiba2*hgkp(ig,ikind)*c0(ig,is1)+fp)
                C2(ig+ncpw%ngw,is1)=-fi*&
                     (0.5_real_8*parm%tpiba2*hgkm(ig,ikind)*c0(ig+ncpw%ngw,is1)+fm)
                ! elisa
                DO id=1,3
                   IF (gampl(id).NE.0._real_8) THEN
                      C2(ig,is1)=c2(ig,is1)&
                           -fi*gpotv(id)*(rk(id,1)&
                           +gk(id,ig))*parm%tpiba*c0(ig ,is1)&
                           -fi*0.5_real_8*(gpotv(id)**2)*   c0(ig,is1)
                      C2(ig+ncpw%ngw,is1)=c2(ig+ncpw%ngw,is1)&
                           -fi*gpotv(id)*(rk(id,1)&
                           -gk(id,ig))*parm%tpiba*c0(ig+ncpw%ngw,is1)&
                           -fi*0.5_real_8*(gpotv(id)**2)*   c0(ig+ncpw%ngw,is1)
                   ENDIF
                ENDDO
                ! elisa
             ENDDO
             !$omp end do nowait
          ENDIF
       ELSE
          !$omp do
          DO ig=1,ncpw%ngw
             fp=psi(nzhs(ig),is)
             fm=psi(indzs(ig),is)
             C2(ig,    is1)=-fi*&
                  (0.5_real_8*parm%tpiba2*hgkp(ig,ikind)*c0(ig,    is1)+fp)
             C2(ig+ncpw%ngw,is1)=-fi*&
                  (0.5_real_8*parm%tpiba2*hgkm(ig,ikind)*c0(ig+ncpw%ngw,is1)+fm)
          END DO
          !$omp end do nowait
       ENDIF
       ! EHR]
    END DO
    !$omp end parallel
    IF(geq0)THEN
       DO is=1,n
          is1=is+is_start
          C2(1+ncpw%ngw,is1)=CMPLX(0._real_8,0._real_8,kind=real_8)
       END DO
    END IF
  END SUBROUTINE calc_c2_kpnt
  ! ==================================================================
  SUBROUTINE calc_c2(c2,c0,psi,f,n,is_start,njump,nostat)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: n, is_start, njump, nostat
    COMPLEX(real_8),INTENT(IN)               :: psi(fpar%kr1s*msrays,n)
    COMPLEX(real_8),INTENT(IN) __CONTIGUOUS  :: c0(:,:)
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: c2(:,:)
    REAL(real_8), INTENT(IN) __CONTIGUOUS    :: f(:)
    COMPLEX(real_8)                          :: fm, fp, psii, psin
    INTEGER                                  :: is, is1, is2, ig, offset
    REAL(real_8)                             :: fi, fip1, xskin, g2

    !$omp parallel private(is,is1,is2,ig,fp,fm,psin,psii,fi,fip1,xskin,g2)
    DO is=1,n
       IF(njump.EQ.1)THEN
          is1=is_start + is
          is2=nostat + 1
       ELSE
          is1=is_start + (is-1)*2 + 1
          is2=is1+1
       END IF
       fi=f(is1)*0.5_real_8
       IF (fi.EQ.0._real_8.AND..NOT.cntl%tksham) fi=1._real_8
       IF (fi.EQ.0._real_8.AND.cntl%tksham) fi=0.5_real_8
       fip1=0._real_8
       IF (is2.LE.nostat) fip1=f(is2)*0.5_real_8
       IF (fip1.EQ.0._real_8.AND..NOT.cntl%tksham) fip1=1._real_8
       IF (fip1.EQ.0._real_8.AND.cntl%tksham) fip1=0.5_real_8
       IF (prcp_com%akin.GT.1.e-10_real_8) THEN
          xskin=1._real_8/prcp_com%gskin
          IF(is2.LE.nostat)THEN
             !$omp do
             DO ig=1,jgw
                psin=psi(nzhs(ig),is)! access only once
                psii=psi(indzs(ig),is)! these mem locations
                fp=psin+psii
                fm=psin-psii
                g2=parm%tpiba2*(hg(ig)+&
                     prcp_com%gakin*(1._real_8+cp_erf((hg(ig)-prcp_com%gckin)*xskin)))
                C2(ig,is1)=-fi*(g2*c0(ig,is1)+&
                     CMPLX(REAL(fp),AIMAG(fm),kind=real_8))
                C2(ig,is2)=-fip1*(g2*c0(ig,is2)+&
                     CMPLX(AIMAG(fp),-REAL(fm),kind=real_8))
             END DO
             !$omp end do nowait
          ELSE
             !$omp do
             DO ig=1,jgw
                psin=psi(nzhs(ig),is)! access only once
                psii=psi(indzs(ig),is)! these mem locations
                fp=psin+psii
                fm=psin-psii
                g2=parm%tpiba2*(hg(ig)+&
                     prcp_com%gakin*(1._real_8+cp_erf((hg(ig)-prcp_com%gckin)*xskin)))
                C2(ig,is1)=-fi*(g2*c0(ig,is1)+&
                     CMPLX(REAL(fp),AIMAG(fm),kind=real_8))
             END DO
             !$omp end do nowait
          END IF
       ELSE
          IF(is2.LE.nostat)THEN
             !$omp do
             DO ig=1,jgw
                psin=psi(nzhs(ig),is)! access only once
                psii=psi(indzs(ig),is)! these mem locations
                fp=psin+psii
                fm=psin-psii
                C2(ig,is1)=-fi*((parm%tpiba2*hg(ig))*c0(ig,is1)+&
                     CMPLX(REAL(fp),AIMAG(fm),kind=real_8))
                C2(ig,is2)=-fip1*&
                     ((parm%tpiba2*hg(ig))*&
                     c0(ig,is2)+CMPLX(AIMAG(fp),-REAL(fm),kind=real_8))
             END DO
             !$omp end do nowait
          ELSE
             !$omp do
             DO ig=1,jgw
                psin=psi(nzhs(ig),is)! access only once
                psii=psi(indzs(ig),is)! these mem locations
                fp=psin+psii
                fm=psin-psii
                C2(ig,is1)=-fi*((parm%tpiba2*hg(ig))*c0(ig,is1)+&
                     CMPLX(REAL(fp),AIMAG(fm),kind=real_8))
             END DO
             !$omp end do nowait
          END IF
       END IF
    END DO
    !$omp end parallel
  END SUBROUTINE calc_c2

  SUBROUTINE do_the_vpsi_thing( psi, hpsi, v, nbnd_source )
  
    IMPLICIT NONE
    COMPLEX(DP), INTENT(INOUT)  :: psi(dfft%ngw, nbnd_source)
    COMPLEX(DP), INTENT(INOUT) :: hpsi(dfft%ngw, nbnd_source)
    INTEGER, INTENT(IN)    :: nbnd_source
    REAL(DP), INTENT(IN) :: v(dfft%nnr)

    COMPLEX(DP), SAVE, POINTER  :: rs_wave(:,:)
    COMPLEX(DP), ALLOCATABLE, SAVE, TARGET :: batch_aux(:,:)

    LOGICAL, ALLOCATABLE, SAVE :: first_step(:)
    LOGICAL, SAVE :: do_calc = .false.
    LOGICAL, SAVE :: do_com = .false.
    INTEGER, SAVE :: sendsize_rem, nthreads, nested_threads, nbnd, my_thread_num
    INTEGER, SAVE :: times_called = 0
    REAL(DP), SAVE, ALLOCATABLE :: final_time(:)
    INTEGER :: i, j, iter, ierr, buffer_size, batch_size, ng, ngms
  
    INTEGER :: counter( 2 , 3)
    LOGICAL :: finished( 2 , 3 )
    INTEGER :: sendsize_save, next, last_buffer, ibatch, work_buffer, z_set_size, y_set_size, scatter_set_size
    LOGICAL :: finished_all, last_start, last_start_triggered

    INTEGER(INT64) :: time(20)
    INTEGER(INT64), SAVE :: cr
    REAL(DP) :: timer(29)
    LOGICAL, SAVE :: first = .true.

    dfft%wave = .true.
    dfft%vpsi = .true.
    dfft%time_adding = 0
    dfft%counter = 0
    timer = 0.0d0

    batch_size =  dfft%batch_size_save
    buffer_size = dfft%buffer_size_save
    z_set_size = dfft%z_set_size_save
    y_set_size = dfft%y_set_size_save
    scatter_set_size = dfft%scatter_set_size_save
    ngms = dfft%ngw

    CALL SYSTEM_CLOCK( count_rate = cr )
    CALL SYSTEM_CLOCK( time(3) )

    IF( first ) THEN
 
       first = .false.
       dfft%use_maps = .true.
       dfft%ngms = ngms

       IF( dfft%overlapp ) THEN
          !$omp parallel
          dfft%cpus_per_task = omp_get_num_threads()
          !$omp end parallel
          nthreads = MIN( 2, dfft%cpus_per_task )
          nested_threads = MAX( 1, dfft%cpus_per_task - 1 )
       ELSE
          nthreads = 1
       END IF
       dfft%uneven  = .false.
  
       IF( dfft%my_node_rank .ne. 0 .or. dfft%single_node ) nthreads = 1
       
       nbnd = ( nbnd_source + 1 ) / 2
       IF( mod( nbnd_source, 2 ) .ne. 0 ) dfft%uneven = .true.
       IF( dfft%mype .eq. 0 .and. dfft%uneven ) THEN
          write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write(6,*) "WARNING: UNEVEN NUMBER OF STATES, NOT TESTED!"
          write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       END IF
  
    END IF

    IF( dfft%remember_batch_vpsi .ne. batch_size .or. dfft%remember_buffer_vpsi .ne. buffer_size ) THEN

       dfft%remember_batch_vpsi  = batch_size
       dfft%remember_buffer_vpsi = buffer_size
       CALL Clean_up_shared( dfft, 1 ) 

       CALL Set_Req_Vals( dfft, nbnd_source, batch_size, dfft%rem_size, buffer_size, dfft%ir1w, dfft%nsw )
       CALL Prep_Copy_Maps( dfft, ngms, batch_size, dfft%rem_size, dfft%ir1w, dfft%nsw )
  
       dfft%sendsize = MAXVAL ( dfft%nr3p ) * MAXVAL( dfft%nsw ) * dfft%node_task_size * dfft%node_task_size * batch_size
!       indi_size = MAXVAL ( dfft%nr3p ) * MAXVAL( dfft%nsw ) * dfft%node_task_size * dfft%node_task_size * batch_size
     
       CALL create_shared_memory_window_2d( comm_send, 1, dfft, dfft%sendsize*dfft%nodes_numb, buffer_size ) 
       CALL create_shared_memory_window_2d( comm_recv, 2, dfft, dfft%sendsize*dfft%nodes_numb, buffer_size ) 
     
       CALL create_shared_locks_2d( locks_calc_inv, 20, dfft, dfft%node_task_size, ( nbnd_source / batch_size ) + 1 )
       CALL create_shared_locks_2d( locks_calc_fw,  21, dfft, dfft%node_task_size, ( nbnd_source / batch_size ) + 1 )
       CALL create_shared_locks_1d( locks_com_inv,  50, dfft, ( nbnd_source / batch_size ) + 1 )
       CALL create_shared_locks_1d( locks_com_fw,   51, dfft, ( nbnd_source / batch_size ) + 1 )
       
       CALL create_shared_locks_2d( locks_calc_1  , 22, dfft, dfft%node_task_size, nbnd_source + batch_size + (buffer_size-1)*batch_size )
       CALL create_shared_locks_2d( locks_calc_2  , 23, dfft, dfft%node_task_size, nbnd_source + batch_size + (buffer_size-1)*batch_size )

       IF( ALLOCATED( batch_aux ) )        DEALLOCATE( batch_aux )
       ALLOCATE( batch_aux( dfft%my_nr3p * dfft%nr2 * dfft%nr1, batch_size ) )
  
       dfft%num_buff = buffer_size
       IF( dfft%rem_size .ne. 0 ) THEN
          dfft%max_nbnd = ( nbnd / batch_size ) + 1
       ELSE
          dfft%max_nbnd = ( nbnd / batch_size )
       END IF
       IF( ALLOCATED( first_step ) )   DEALLOCATE( first_step )
       ALLOCATE( first_step( dfft%buffer_size_save ) ) 
  
    END IF
  
    sendsize_rem = MAXVAL ( dfft%nr3p ) * MAXVAL( dfft%nsw ) * dfft%node_task_size * dfft%node_task_size * dfft%rem_size
    sendsize_save = dfft%sendsize
    last_start = .true.
    last_start_triggered = .false.
  
    locks_calc_inv = .true.
    locks_calc_fw  = .true.
    locks_com_inv  = .true.
    locks_com_fw   = .true.
  
    locks_calc_1   = .true.
    locks_calc_2   = .true.
    IF( .not. dfft%rsactive ) THEN
       DO i = 1, batch_size*buffer_size
          locks_calc_1( : , i ) = .false.
       ENDDO
    ELSE
       DO i = 1, batch_size*buffer_size
          locks_calc_2( : , i ) = .false.
       ENDDO
    END IF

    dfft%first_loading = .true.
    dfft%rem = .false.

    IF( .not. dfft%rsactive) THEN
       rs_wave => batch_aux( : , 1:batch_size )
    END IF
 
    CALL MPI_BARRIER(dfft%comm, ierr)

    CALL SYSTEM_CLOCK( time(1) )
  
    !!! Initializiation finished
  
    !$omp parallel IF( nthreads .eq. 2 ) num_threads( nthreads ) &
    !$omp private( my_thread_num, ibatch, batch_size, do_calc, do_com, next, counter, last_buffer, finished_all, first_step, work_buffer, finished, z_set_size, y_set_size, scatter_set_size ) &
    !$omp proc_bind( close )
    
    do_calc = .false.
    do_com  = .false.
    batch_size = dfft%batch_size_save
    z_set_size = dfft%z_set_size_save
    y_set_size = dfft%y_set_size_save
    scatter_set_size = dfft%scatter_set_size_save
    next = -1
    counter = 0
    finished = .false.
    last_buffer = 0
    finished_all = .false.
    first_step = dfft%first_step

    !$  IF( nthreads .eq. 2 ) THEN
    !$     my_thread_num = omp_get_thread_num()
    !$     IF( my_thread_num .eq. 1 ) THEN
    !$        CALL omp_set_max_active_levels( 2 )
    !$        CALL omp_set_num_threads( nested_threads )
    !$        CALL mkl_set_dynamic(0) 
    !$        ierr = mkl_set_num_threads_local( nested_threads )
    !$        do_calc = .true.
    !$     ELSE
    !$        do_com  = .true.
    !$     END IF
    !$  END IF
  
    IF( nthreads .eq. 1 ) THEN
       IF( dfft%my_node_rank .eq. 0 .and. .not. dfft%single_node ) do_com = .true.
!       CALL omp_set_num_threads( 5 )
       my_thread_num = 1
       do_calc = .true.
    ENDIF
   
    DO WHILE( .not. finished_all )
  
       next = mod( next, buffer_size ) + 1
       IF( next .eq. 0 ) THEN
          work_buffer = 1
       ELSE
          work_buffer = dfft%buffer_sequence( next )
       END IF
  
       IF( dfft%rem_size .ne. 0 .and. ( ( first_step( work_buffer ) .and. last_buffer .eq. 0 .and. ( counter( 1, 1 ) .eq. dfft%max_nbnd - 1 .or. counter( 2, 1 ) .eq. dfft%max_nbnd - 1 ) ) .or. last_buffer .eq. work_buffer ) ) THEN
          last_buffer = work_buffer
          batch_size = dfft%rem_size
          IF( z_set_size .gt. (batch_size+1)/2 ) z_set_size = batch_size
          IF( y_set_size .gt. (batch_size+1)/2 ) y_set_size = batch_size
          IF( scatter_set_size .gt. (batch_size+1)/2 ) scatter_set_size = batch_size
          IF( do_com  ) dfft%sendsize = sendsize_rem
          IF( do_calc ) dfft%rem = .true.
       END IF
    
       IF( first_step( work_buffer ) ) THEN
  
          ! Last and First Step
  
          IF( do_calc ) THEN
  
             ! Last Step
  
             IF( dfft%first_loading( work_buffer ) ) THEN
                dfft%first_loading( work_buffer ) = .false.
             ELSE
  
                counter( 1, 3 ) = counter( 1, 3 ) + 1
  
                IF( batch_size .ne. dfft%batch_size_save .and. last_start ) THEN
                   last_start = .false.
                   last_start_triggered = .true.
                   batch_size = dfft%batch_size_save
                   z_set_size = dfft%z_set_size_save
                END IF
  
                CALL SYSTEM_CLOCK( time(5) )
                CALL fwfft_pwbatch( dfft, 3, batch_size, z_set_size, counter( 1, 3 ), work_buffer, psi, hpsi, comm_recv(:,work_buffer) )
                CALL SYSTEM_CLOCK( time(6) )
                dfft%time_adding( 18 ) = dfft%time_adding( 18 ) + ( time(6) - time(5) )
  
                IF( last_start_triggered ) THEN
                   batch_size = dfft%rem_size
                   IF( z_set_size .gt. (batch_size+1)/2 ) z_set_size = batch_size
                END IF
                   
  
                IF( counter( 1, 3 ) .eq. dfft%max_nbnd ) finished_all = .true.
  
             END IF
  
             ! First Step
  
             IF( finished( 1, 1 ) .or. dfft%rsactive ) THEN
                CONTINUE
             ELSE
  
                counter( 1, 1 ) = counter( 1, 1 ) + 1

                CALL SYSTEM_CLOCK( time(10) )
                CALL invfft_pwbatch( dfft, 1, batch_size, z_set_size, 0, 0, counter( 1, 1 ), work_buffer, psi, comm_send, comm_recv(:,work_buffer) )
                CALL SYSTEM_CLOCK( time(11) )
                dfft%time_adding( 19 ) = dfft%time_adding( 19 ) + ( time(11) - time(10) )

                IF( counter( 1, 1 ) .eq. dfft%max_nbnd ) finished( 1, 1 ) = .true.
  
             END IF
     
          END IF
     
          IF( do_com ) THEN
     
             IF( finished( 2, 1 ) .or. dfft%rsactive) THEN
                CONTINUE
             ELSE
  
                counter( 2, 1 ) = counter( 2, 1 ) + 1
  
                CALL SYSTEM_CLOCK( time(12) )
                CALL invfft_pwbatch( dfft, 2, batch_size, 0, 0, 0, counter( 2, 1 ), work_buffer, comm_send, comm_recv )
                CALL SYSTEM_CLOCK( time(13) )
                dfft%time_adding( 16 ) = dfft%time_adding( 16 ) + ( time(13) - time(12) )
 
                IF( counter( 2, 1 ) .eq. dfft%max_nbnd ) finished( 2, 1 ) = .true.
  
             END IF
     
          END IF
  
          IF( dfft%single_node ) CALL MPI_BARRIER(dfft%comm, ierr) 
  
          first_step( work_buffer ) = .false.
  
       ELSE
  
          ! Middle Step
  
          IF( do_calc ) THEN
     
             IF( finished( 1, 2 ) ) THEN
                CONTINUE
             ELSE
  
                counter( 1, 2 ) = counter( 1, 2 ) + 1

                CALL SYSTEM_CLOCK( time(14) )
                IF( dfft%rsactive ) THEN
                   rs_wave => wfn_real( : , 1+(counter(1,2)-1)*batch_size : batch_size+(counter(1,2)-1)*batch_size )
                   CONTINUE
                ELSE
                   CALL invfft_pwbatch( dfft, 3, batch_size, y_set_size, scatter_set_size, 0, counter( 1, 2 ), work_buffer, comm_recv, rs_wave )
                END IF
 
                DO ibatch = 1, batch_size
                   CALL Apply_V( dfft, rs_wave(:,ibatch), v, ibatch )
                ENDDO
                CALL SYSTEM_CLOCK( time(15) )
                dfft%time_adding( 20 ) = dfft%time_adding( 20 ) + ( time(15) - time(14) )
  
                CALL SYSTEM_CLOCK( time(16) )
                CALL fwfft_pwbatch( dfft, 1, batch_size, 1, counter( 1, 2 ), work_buffer, rs_wave, comm_send, comm_recv(:,work_buffer) )
                CALL SYSTEM_CLOCK( time(17) )
                dfft%time_adding( 21 ) = dfft%time_adding( 21 ) + ( time(17) - time(16) )
  
                IF( counter( 1, 2 ) .eq. dfft%max_nbnd ) finished( 1, 2 ) = .true.
  
             END IF
  
          END IF
  
          IF( do_com ) THEN
  
             IF( finished( 2, 2 ) ) THEN
                CONTINUE
             ELSE
  
                counter( 2, 2 ) = counter( 2, 2 ) + 1
  
                CALL SYSTEM_CLOCK( time(18) )
                CALL fwfft_pwbatch( dfft, 2, batch_size, 1, counter( 2, 2 ), work_buffer, comm_send, comm_recv )
                CALL SYSTEM_CLOCK( time(19) )
                dfft%time_adding( 17 ) = dfft%time_adding( 17 ) + ( time(19) - time(18) )
  
                IF( counter( 2, 2 ) .eq. dfft%max_nbnd ) THEN
                   finished( 2, 2 ) = .true.
                   IF( .not. do_calc ) THEN
                      finished_all = .true.
                   END IF
                END IF
  
             END IF
     
          END IF
  
          IF( dfft%single_node ) CALL MPI_BARRIER(dfft%comm, ierr) 
  
          first_step( work_buffer ) = .true.
  
       END IF
  
       IF( batch_size .ne. dfft%batch_size_save ) THEN
          batch_size = dfft%batch_size_save
          z_set_size = dfft%z_set_size_save
          y_set_size = dfft%y_set_size_save
          scatter_set_size = dfft%scatter_set_size_save
          IF( do_com  ) dfft%sendsize = sendsize_save
          IF( do_calc ) dfft%rem = .false.
       END IF
   
  
    ENDDO 
  
    !$  IF( nthreads .eq. 2 .and. my_thread_num .eq. 1 ) THEN
    !$     CALL omp_set_max_active_levels( 1 )
    !$     CALL omp_set_num_threads( dfft%cpus_per_task )
    !$     CALL mkl_set_dynamic(1) 
    !$     ierr = mkl_set_num_threads_local( dfft%cpus_per_task )
    !$  END IF
    !$omp barrier
    !$omp end parallel
  
    CALL SYSTEM_CLOCK( time(4) )
    dfft%time_adding( 98 ) = time(4) - time(1)

    DO i = 1, 29
       timer( i ) = REAL( dfft%time_adding( i ), KIND = REAL64 ) / REAL ( cr , KIND = REAL64 )
    ENDDO

    IF( dfft%mype .eq. 0 .and. .false. ) THEN

          WRITE(6,*)" "
          WRITE(6,*)"Some extra VPSI times"
          write(6,*)"==================================="
          write(6,*)"INV FFT before Com"
          write(6,*)"CALC LOCK 1", timer(22)
          write(6,*)"Prepare_psi ",           timer(1)
          write(6,*)"INV z_fft ",             timer(2)
          write(6,*)"INV Pre_com_copy ",      timer(3)
  
          write(6,*)"INV FFT before Com sum  ",timer(1)+timer(2)+timer(3)+timer(22)
          write(6,*)"Control:", timer(19)
          write(6,*)"==================================="
          write(6,*)"INV FFT after Com"
          write(6,*)"CALC LOCK 2", timer(24)
          write(6,*)"INV After_com_copy ",    timer(4)
          write(6,*)"INV y_fft ",             timer(5)
          write(6,*)"INV xy_scatter ",        timer(6)
          write(6,*)"INV x_fft ",             timer(7)
          write(6,*)"Apply V ",               timer(8)
  
          write(6,*)"INV FFT after Com sum ",  timer(4)+timer(5)+timer(6)+timer(7)+timer(8)+timer(24)
          write(6,*)"Control:", timer(20)
          write(6,*)"==================================="
          write(6,*)"FW FFT before Com"
          write(6,*)"CALC LOCK 3", timer(25)
          write(6,*)"FW x_fft ",             timer(9)
          write(6,*)"FW xy_scatter ",        timer(10)
          write(6,*)"FW y_fft ",             timer(11)
          write(6,*)"FW Pre_com_copy ",      timer(12)
  
          write(6,*)"FW FFT before Com sum  ",timer(9)+timer(10)+timer(11)+timer(12)+timer(25)
          write(6,*)"Control:", timer(21)
          write(6,*)"==================================="
          write(6,*)"FW FFT after Com"
          write(6,*)"CALC LOCK 4", timer(27)
          write(6,*)"FW After_com_copy ",    timer(13)
          write(6,*)"FW z_fft ",             timer(14)
          write(6,*)"Accumulate_Psi ",       timer(15)
  
          write(6,*)"FW FFT after Com sum",  timer(13)+timer(14)+timer(15)+timer(27)
          write(6,*)"Control:", timer(18)
          write(6,*)"==================================="
          write(6,*)"COM LOCK 1", timer(23)
          write(6,*)"FIRST COMM TIMES:", timer(28)
          write(6,*)"Control:", timer(16)
          write(6,*)"COM LOCK 2", timer(26)
          write(6,*)"SECOND COMM TIMES:", timer(29)
          write(6,*)"Control:", timer(17)
          write(6,*)"==================================="
          write(6,*)"Adding up CALC:", timer(1)+timer(2)+timer(3)+timer(4)+timer(5)+timer(6)+&
                                  timer(7)+timer(8)+timer(9)+timer(10)+timer(11)+timer(12)+&
                                  timer(13)+timer(14)+timer(15)+timer(22)+timer(24)+&
                                  timer(25)+timer(27)
          write(6,*)"Adding up COMM:", timer(23)+timer(26)+timer(28)+timer(29)
          write(6,*)"Control:", REAL( REAL( dfft%time_adding( 98 ) ) / REAL( cr ), KIND = REAL64 )
          WRITE(6,*)" "

    END IF

    dfft%wave = .false.
    dfft%vpsi = .false.
  
  END SUBROUTINE do_the_vpsi_thing

END MODULE vpsi_utils
