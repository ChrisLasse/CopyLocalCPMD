MODULE fftnew_utils
  USE cell,                            ONLY: cell_com
  USE cnst,                            ONLY: pi
  USE cp_cuda_types,                   ONLY: cp_cuda_env
  USE cp_cufft_types,                  ONLY: cp_cufft,&
                                             cp_cufft_device_get_ptrs
  USE cppt,                            ONLY: hg,&
                                             indz,&
                                             indzs,&
                                             inyh,&
                                             nzh,&
                                             nzhs
  USE cuda_types,                      ONLY: cuda_memory_t
  USE cuda_utils,                      ONLY: cuda_memcpy_host_to_device
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: &
       fftpool, fftpoolsize, fpoolv, inzf, inzfp, inzh, inzhp, inzs, inzsp, &
       jgw, jgws, jhg, jhgs, kr2max, kr2min, kr3max, kr3min, lfrm, llr1, &
       lmsq, lmsqmax, lnzf, lnzs, lr1, lr1m, lr1s, lr2, lr2s, lr3, lr3s, &
       lrxpl, lrxpool, lsrm, mfrays, mg, msp, msqf, msqfpool, msqs, msqspool, &
       msrays, mz, ngrm, nhrm, nr1m, nzff, nzffp, nzfs, nzfsp, qr1, qr1s, &
       qr2, qr2max, qr2min, qr2s, qr3, qr3max, qr3min, qr3s, sp5, sp8, sp9, &
       spm, FFT_TYPE_DESCRIPTOR
  USE fft_maxfft,                      ONLY: maxfft,&
                                             maxfftn
  USE fftchk_utils,                    ONLY: fftchk
  USE kinds,                           ONLY: real_8
  USE loadpa_utils,                    ONLY: leadim
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum,&
                                             mp_send_init_COMPLEX,&
                                             mp_recv_init_COMPLEX
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             parap,&
                                             parm,&
                                             spar,&
                                             cntl
  USE utils,                           ONLY: icopy
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setfftn
  !public :: rmfftnset
  PUBLIC :: addfftnset
  !public :: setrays
  PUBLIC :: Prep_fft_com
  PUBLIC :: Make_Manual_Maps
  PUBLIC :: Make_inv_yzCOM_Maps

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

CONTAINS

  ! ==================================================================
  SUBROUTINE setfftn(ipool)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ipool

    INTEGER                                  :: i_device, ip, ipx
    TYPE(cuda_memory_t), POINTER             :: inzs_d, lrxpl_d, msqf_d, &
                                                msqs_d, nzfs_d, sp5_d, sp8_d, &
                                                sp9_d

    IF (ipool.EQ.0) THEN
       msrays   = parai%ngrays
       mfrays   = parai%nhrays
       llr1     = fpar%nnr1
       qr1s     = fpar%kr1s
       qr2s     = fpar%kr2s
       qr3s     = fpar%kr3s
       qr1      = fpar%kr1
       lr1s     = spar%nr1s
       lr2s     = spar%nr2s
       lr3s     = spar%nr3s
       lr1      = parm%nr1
       qr2max   = kr2max
       qr2min   = kr2min
       qr3max   = kr3max
       qr3min   = kr3min
       lsrm     = ngrm
       lfrm     = nhrm
       lr1m     = nr1m
       lmsq     = nhrm
       maxfftn  = maxfft
       jgw      = ncpw%ngw
       jgws     = spar%ngws
       jhg      = ncpw%nhg
       jhgs     = spar%nhgs
       nzff => nzh
       inzf => indz
       nzfs => nzhs
       inzs => indzs
       inzh => inyh

       DO ip=0,parai%nproc-1
          lrxpl(ip,1)    = parap%nrxpl(ip,1)
          lrxpl(ip,2)    = parap%nrxpl(ip,2)
          sp5(ip)        = parap%sparm(5,ip)
          sp8(ip)        = parap%sparm(8,ip)
          sp9(ip)        = parap%sparm(9,ip)
          ipx=lmsq*ip
          IF (nhrm.GT.0) THEN
             CALL icopy(nhrm,msp(1,1,ip+1),1,msqf(ipx+1),1)
             CALL icopy(nhrm,msp(1,2,ip+1),1,msqs(ipx+1),1)
          ENDIF
       ENDDO
    ELSEIF (ipool.GT.0 .AND. ipool.LE.fftpool) THEN
       ! LOAD FROM POOL
       msrays   = fpoolv( 1,ipool)
       mfrays   = fpoolv( 2,ipool)
       llr1     = fpoolv( 3,ipool)
       qr1s     = fpoolv( 4,ipool)
       qr2s     = fpoolv( 5,ipool)
       qr3s     = fpoolv( 6,ipool)
       qr1      = fpoolv( 7,ipool)
       qr2      = fpoolv( 8,ipool)
       qr3      = fpoolv( 9,ipool)
       lr1s     = fpoolv(10,ipool)
       lr2s     = fpoolv(11,ipool)
       lr3s     = fpoolv(12,ipool)
       lr1      = fpoolv(13,ipool)
       lr2      = fpoolv(14,ipool)
       lr3      = fpoolv(15,ipool)
       qr2max   = fpoolv(16,ipool)
       qr2min   = fpoolv(17,ipool)
       qr3max   = fpoolv(18,ipool)
       qr3min   = fpoolv(19,ipool)
       lsrm     = fpoolv(20,ipool)
       lfrm     = fpoolv(21,ipool)
       lr1m     = fpoolv(22,ipool)
       lmsq     = fpoolv(23,ipool)
       maxfftn  = fpoolv(24,ipool)
       jgw      = fpoolv(25,ipool)
       jgws     = fpoolv(26,ipool)
       jhg      = fpoolv(27,ipool)
       jhgs     = fpoolv(28,ipool)
       nzff => nzffp(:, ipool)
       inzf => inzfp(:, ipool)
       nzfs => nzfsp(:, ipool)
       inzs => inzsp(:, ipool)
       inzh => inzhp(:, :, ipool)

       DO ip=0,parai%nproc-1
          lrxpl(ip,1)    = lrxpool(ip,1,ipool)
          lrxpl(ip,2)    = lrxpool(ip,2,ipool)
          sp5(ip)        = spm(5,ip,ipool)
          sp8(ip)        = spm(8,ip,ipool)
          sp9(ip)        = spm(9,ip,ipool)
          ipx=lmsq*ip
          IF (lmsqmax.GT.0) THEN
             CALL icopy(lmsq,msqfpool(1,ip+1,ipool),1,msqf(ipx+1),1)
             CALL icopy(lmsq,msqspool(1,ip+1,ipool),1,msqs(ipx+1),1)
          ENDIF
       ENDDO
    ELSE
       CALL stopgm("SETFFTN","FFTPOOL NOT DEFINED",& 
            __LINE__,__FILE__)
    ENDIF

    !vw copy FFT arrays to GPU memory
    IF( cp_cuda_env%use_fft ) THEN
       DO i_device = 1, cp_cuda_env%fft_n_devices_per_task
          CALL cp_cufft_device_get_ptrs ( cp_cufft, i_device, sp5_d=sp5_d, sp8_d=sp8_d, sp9_d=sp9_d, &
               & msqs_d=msqs_d, msqf_d=msqf_d, lrxpl_d=lrxpl_d, nzfs_d=nzfs_d, inzs_d=inzs_d )
          CALL cuda_memcpy_host_to_device ( sp5, sp5_d )
          CALL cuda_memcpy_host_to_device ( sp8, sp8_d )
          CALL cuda_memcpy_host_to_device ( sp9, sp9_d )
          CALL cuda_memcpy_host_to_device ( lrxpl, lrxpl_d )
          CALL cuda_memcpy_host_to_device ( msqf, msqf_d )
          CALL cuda_memcpy_host_to_device ( msqs, msqs_d )
          CALL cuda_memcpy_host_to_device ( nzfs, nzfs_d )
          CALL cuda_memcpy_host_to_device ( inzs, inzs_d )
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE setfftn
  ! ==================================================================
  SUBROUTINE rmfftnset(ipool)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ipool

    INTEGER                                  :: i, ip, j

    IF (ipool.GT.0 .AND. ipool.LE.fftpool) THEN
       DO ip=ipool+1,fftpool
          i=ip
          j=ip-1
          CALL icopy(28,fpoolv(1,i),1,fpoolv(1,j),1)
          CALL icopy(lnzf,nzffp(1,i),1,nzffp(1,j),1)
          CALL icopy(lnzf,inzfp(1,i),1,inzfp(1,j),1)
          CALL icopy(lnzs,nzfsp(1,i),1,nzfsp(1,j),1)
          CALL icopy(lnzs,inzsp(1,i),1,inzsp(1,j),1)
          CALL icopy(3*lnzf,inzhp(1,1,i),1,inzhp(1,1,j),1)
          CALL icopy(2*SIZE(lrxpool,1),lrxpool(0,1,i),1,lrxpool(0,1,j),1)
          CALL icopy(9*SIZE(lrxpool,1),spm(1,0,i),1,spm(1,0,j),1)
          IF (lmsqmax.GT.0) THEN
             CALL icopy(lmsqmax*parai%nproc,msqfpool(1,1,i),1,&
                  msqfpool(1,1,j),1)
             CALL icopy(lmsqmax*parai%nproc,msqspool(1,1,i),1,&
                  msqspool(1,1,j),1)
          ENDIF
       ENDDO
       CALL zeroing(fpoolv(:,fftpool))!,28)
       CALL zeroing(nzffp(:,fftpool))!,lnzf)
       CALL zeroing(inzfp(:,fftpool))!,lnzf)
       CALL zeroing(nzfsp(:,fftpool))!,lnzf)
       CALL zeroing(inzsp(:,fftpool))!,lnzf)
       CALL zeroing(inzhp(:,:,fftpool))!,3*lnzf)
       CALL zeroing(lrxpool(:,:,fftpool))!,2*maxcpu+2)
       CALL zeroing(spm(:,:,fftpool))!,9*maxcpu+9)
       IF (lmsqmax.GT.0) THEN
          CALL zeroing(msqfpool(:,:,fftpool))!,lmsqmax*parai%nproc)
          CALL zeroing(msqspool(:,:,fftpool))!,lmsqmax*parai%nproc)
       ENDIF
       fftpool=fftpool-1
    ELSE
       CALL stopgm("RMFFTNSET","FFTPOOL NOT DEFINED",& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rmfftnset
  ! ==================================================================
  SUBROUTINE addfftnset(ecutf,ecuts,ipool)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: ecutf, ecuts
    INTEGER                                  :: ipool

    CHARACTER(*), PARAMETER                  :: procedureN = 'addfftnset'

    INTEGER                                  :: i, ierr, ig, j, k, l, lh1, &
                                                lh2, lh3, nh1, nh2, nh3
    INTEGER, SAVE                            :: icount = 0
    REAL(real_8)                             :: aa1, aa2, aa3, rr, xpaim, &
                                                xplanes, xpnow

! ==--------------------------------------------------------------==

    IF (icount.EQ.0) THEN
       lmsqmax=nhrm
       CALL zeroing(lrxpool)!,2*fftpoolsize*(maxcpu+1))
       CALL zeroing(spm)!,9*fftpoolsize*(maxcpu+1))
       l=(lmsqmax*parai%nproc)
       ALLOCATE(msqf(l),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(msqs(l),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       l=(lmsqmax*parai%nproc*fftpoolsize)
       ALLOCATE(msqfpool(lmsqmax,parai%nproc,fftpoolsize),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(msqspool(lmsqmax,parai%nproc,fftpoolsize),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       l=(ncpw%nhg*fftpoolsize)
       ALLOCATE(nzffp(lnzf,l/lnzf),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       IF (lnzs > 0) THEN
          ALLOCATE(nzfsp(lnzs,l/lnzs),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(nzfsp(1,l),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF

       ALLOCATE(inzfp(lnzf,l/lnzf),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       IF (lnzs > 0) THEN
          ALLOCATE(inzsp(lnzs,l/lnzs),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(inzsp(1,l),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF

       l=(3*ncpw%nhg*fftpoolsize)
       ALLOCATE(inzhp(3,lnzf,l/(3*lnzf)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    icount=1
    IF (ecutf.LT.0._real_8 .OR. ecuts.LT.0._real_8) RETURN
    ! 
    ! IF(TKPNT) CALL STOPGM("ADDFFTNSET","NOT IMPLEMENTED")
    ! 
    fftpool=fftpool+1
    IF (fftpool.GT.fftpoolsize) THEN
       CALL stopgm("ADDFFTNSET","TOO MANY ENTRIES IN POOL",& 
            __LINE__,__FILE__)
    ENDIF
    ipool=fftpool
    ! 
    jhg=ncpw%nhg
    DO ig=2,ncpw%nhg
       IF (parm%tpiba2*hg(ig).GT.ecutf) THEN
          jhg=ig-1
          GOTO 100
       ENDIF
    ENDDO
100 CONTINUE
    rr=REAL(jhg,kind=real_8)
    CALL mp_sum(rr,parai%allgrp)
    jhgs=NINT(rr)
    IF (jhgs.GT.spar%nhgs) CALL stopgm("ADDFFTNSET","JHGS TOO LARGE",& 
         __LINE__,__FILE__)
    jgw=ncpw%nhg
    DO ig=2,ncpw%nhg
       IF (parm%tpiba2*hg(ig).GT.ecuts) THEN
          jgw=ig-1
          GOTO 101
       ENDIF
    ENDDO
101 CONTINUE
    rr=REAL(jgw,kind=real_8)
    CALL mp_sum(rr,parai%allgrp)
    jgws=NINT(rr)
    IF (jgws.GT.spar%ngws) CALL stopgm("ADDFFTNSET","JGWS TOO LARGE",& 
         __LINE__,__FILE__)
    ! 
    aa1=parm%alat
    aa2=parm%alat*cell_com%celldm(2)
    aa3=parm%alat*cell_com%celldm(3)
    lr1s=NINT(aa1/pi*SQRT(ecutf)+0.5_real_8)
    lr2s=NINT(aa2/pi*SQRT(ecutf)+0.5_real_8)
    lr3s=NINT(aa3/pi*SQRT(ecutf)+0.5_real_8)
    ! 
    lr1s=fftchk(lr1s,2)
    lr2s=fftchk(lr2s,2)
    lr3s=fftchk(lr3s,2)
    CALL leadim(lr1s,lr2s,lr3s,qr1s,qr2s,qr3s)
    ! 
    nh1=spar%nr1s/2+1
    nh2=spar%nr2s/2+1
    nh3=spar%nr3s/2+1
    lh1=lr1s/2+1
    lh2=lr2s/2+1
    lh3=lr3s/2+1
    !$omp parallel do private(IG,I,J,K) shared(NH1,NH2,NH3,LH1,LH2,LH3)
    DO ig=1,jhg
       i=inyh(1,ig)-nh1
       j=inyh(2,ig)-nh2
       k=inyh(3,ig)-nh3
       inzhp(1,ig,ipool)=lh1+i
       inzhp(2,ig,ipool)=lh2+j
       inzhp(3,ig,ipool)=lh3+k
    ENDDO
    inzh => inzhp(:,:,ipool)
    ! 
    CALL zeroing(lrxpool(:,:,ipool))!,2*(maxcpu+1))
    xplanes=REAL(lr1s,kind=real_8)
    xpnow=0.0_real_8
    DO i=parai%nproc,1,-1
       xpaim = xpnow + xplanes/parai%nproc
       lrxpool(i-1,1,ipool)=NINT(xpnow)+1
       lrxpool(i-1,2,ipool)=NINT(xpaim)
       IF (NINT(xpaim).GT.lr1s) lrxpool(i-1,2,ipool)=lr1s
       IF (i.EQ.1) lrxpool(i-1,2,ipool)=lr1s
       xpnow = xpaim
    ENDDO
    lr1=lrxpool(parai%mepos,2,ipool)-lrxpool(parai%mepos,1,ipool)+1
    CALL leadim(lr1,lr2s,lr3s,qr1,qr2s,qr3s)
    lr2=lr2s
    lr3=lr3s
    qr2=qr2s
    qr3=qr3s
    llr1=qr1*qr2*qr3
    ! 
    CALL setrays(ipool)
    ! 
    maxfftn = maxfft
    ! 
    fpoolv( 1,ipool) = msrays
    fpoolv( 2,ipool) = mfrays
    fpoolv( 3,ipool) = llr1
    fpoolv( 4,ipool) = qr1s
    fpoolv( 5,ipool) = qr2s
    fpoolv( 6,ipool) = qr3s
    fpoolv( 7,ipool) = qr1
    fpoolv( 8,ipool) = qr2
    fpoolv( 9,ipool) = qr3
    fpoolv(10,ipool) = lr1s
    fpoolv(11,ipool) = lr2s
    fpoolv(12,ipool) = lr3s
    fpoolv(13,ipool) = lr1
    fpoolv(14,ipool) = lr2
    fpoolv(15,ipool) = lr3
    fpoolv(16,ipool) = qr2max
    fpoolv(17,ipool) = qr2min
    fpoolv(18,ipool) = qr3max
    fpoolv(19,ipool) = qr3min
    fpoolv(20,ipool) = lsrm
    fpoolv(21,ipool) = lfrm
    fpoolv(22,ipool) = lr1m
    fpoolv(23,ipool) = lmsq
    fpoolv(24,ipool) = maxfftn
    fpoolv(25,ipool) = jgw
    fpoolv(26,ipool) = jgws
    fpoolv(27,ipool) = jhg
    fpoolv(28,ipool) = jhgs
    ! 
    IF (paral%io_parent) THEN
       WRITE(6,*)
       WRITE(6,'(A,T50,A,I4)') ' ADD NEW FFT SET ',' SET NUMBER ',&
            ipooL
       WRITE(6,'(A,T51,3I5)') ' REAL SPACE GRID ',lr1s,lr2s,lr3S
       WRITE(6,'(A,T20,A,F6.0,T44,A,I10)') ' SPARSE FFT SETUP: ',&
            'CUTOFF [Ry]:',ecuts,'PLANE WAVES:',jgwS
       WRITE(6,'(A,T20,A,F6.0,T44,A,I10)') ' FULL FFT SETUP  : ',&
            'CUTOFF [Ry]:',ecutf,'PLANE WAVES:',jhgS
       WRITE(6,*)
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE addfftnset
  ! ==================================================================
  SUBROUTINE setrays(ipool)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ipool

    CHARACTER(*), PARAMETER                  :: procedureN = 'setrays'

    INTEGER :: i, ierr, ig, ij, img, iny1, iny2, iny3, ip, ipro, ixf, j, &
      jgwl, jhgl, jj, jmg, msglen, mxrp, nh1, nh2, nh3, ny1, ny2, ny3, qr1m
    INTEGER, ALLOCATABLE                     :: mq(:), my(:)

! Variables
! ==--------------------------------------------------------------==
! GATHER ARRAY FOR FFT ALONG X
! SPARSITY FOR FFT ALONG Y

    ALLOCATE(mg(qr2s,qr3s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(mz((2*qr3s)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(my((2*qr2s)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(mg)!,qr2s*qr3s)
    CALL zeroing(mz)!,2*qr3s)
    CALL zeroing(my)!,2*qr2s)
    nh1=lr1s/2+1
    nh2=lr2s/2+1
    nh3=lr3s/2+1
    DO ig=1,jgw
       ny2=inzh(2,ig)
       ny3=inzh(3,ig)
       iny2=-ny2+2*nh2
       iny3=-ny3+2*nh3
       mg(ny2,ny3)=mg(ny2,ny3)+1
       mg(iny2,iny3)=mg(iny2,iny3)+1
       my(ny2)=my(ny2)+1
       my(iny2)=my(iny2)+1
       mz(ny3)=mz(ny3)+1
       mz(iny3)=mz(iny3)+1
    ENDDO
    CALL mp_sum(my,my(lr2s+1:),lr2s,parai%allgrp)
    CALL icopy(lr2s,my(lr2s+1),1,my,1)
    CALL mp_sum(mz,mz(lr3s+1:),lr3s,parai%allgrp)
    CALL icopy(lr3s,mz(lr3s+1),1,mz,1)
    qr2min=1
    DO i=1,qr2s
       IF (my(i).NE.0) THEN
          qr2min=i
          GOTO 50
       ENDIF
    ENDDO
50  CONTINUE
    qr2max=qr2s
    DO i=qr2s,1,-1
       IF (my(i).NE.0) THEN
          qr2max=i
          GOTO 51
       ENDIF
    ENDDO
51  CONTINUE
    qr3min=1
    DO i=1,qr3
       IF (mz(i).NE.0) THEN
          qr3min=i
          GOTO 52
       ENDIF
    ENDDO
52  CONTINUE
    qr3max=1
    DO i=qr3,1,-1
       IF (mz(i).NE.0) THEN
          qr3max=i
          GOTO 53
       ENDIF
    ENDDO
53  CONTINUE
    ! ==--------------------------------------------------------------==
    img=0
    DO j=1,qr3s
       DO i=1,qr2s
          IF (mg(i,j).NE.0) THEN
             img=img+1
             mg(i,j)=img
          ENDIF
       ENDDO
    ENDDO
    msrays=img
    DO ig=jgw+1,jhg
       ny2=inzh(2,ig)
       ny3=inzh(3,ig)
       iny2=-ny2+2*nh2
       iny3=-ny3+2*nh3
       jmg=mg(ny2,ny3)
       IF (jmg.EQ.0) mg(ny2,ny3)=-1
       jmg=mg(iny2,iny3)
       IF (jmg.EQ.0) mg(iny2,iny3)=-1
    ENDDO
    DO j=1,qr3s
       DO i=1,qr2s
          IF (mg(i,j).LT.0) THEN
             img=img+1
             mg(i,j)=img
          ENDIF
       ENDDO
    ENDDO
    mfrays=img
    ! 
    jhgl=0
    jgwl=0
    CALL zeroing(spm(:,:,ipool))!,9*maxcpu+9)
    spm(1,parai%mepos,ipool)=jhg
    spm(2,parai%mepos,ipool)=jhgl
    spm(3,parai%mepos,ipool)=jgw
    spm(4,parai%mepos,ipool)=jgwl
    spm(5,parai%mepos,ipool)=lr1
    spm(6,parai%mepos,ipool)=lr2
    spm(7,parai%mepos,ipool)=lr3
    spm(8,parai%mepos,ipool)=mfrays
    spm(9,parai%mepos,ipool)=msrays
    ! 
    nzff => nzffp(:,ipool)
    inzf => inzfp(:,ipool)
    nzfs => nzfsp(:,ipool)
    inzs => inzsp(:,ipool)
    DO i=0,parai%nproc-1
       ipro=i
       CALL mp_bcast(spm(:,i,ipool),9,ipro,parai%allgrp)
    ENDDO
    ! MAXIMUM OF LR1, MSRAYS AND MFRAYS FOR MP_INDEX
    lr1m = 0
    lfrm = 0
    lsrm = 0
    DO i=0,parai%nproc-1
       lr1m = MAX(lr1m,spm(5,i,ipool))
       lfrm = MAX(lfrm,spm(8,i,ipool))
       lsrm = MAX(lsrm,spm(9,i,ipool))
    ENDDO
    qr1m=MAX(lr1m+MOD(lr1m+1,2),qr1)
    lmsq=MAX(lfrm,lsrm)
    ! SCATTER ARRAY FOR FFT ALONG X
    ALLOCATE(mq(lmsq*2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(mq)!,2*lmsq)
    DO i=1,qr2s
       DO j=1,qr3s
          ij=mg(i,j)
          IF (ij.GT.0) THEN
             mq(ij)=i
             mq(lmsq+ij)=j
          ENDIF
       ENDDO
    ENDDO
    ! CONCATENATE GATHER/SCATTER ARRAYS
    msglen = lfrm * 8/2
    CALL my_concat(mq(1),msqs,msglen,parai%allgrp)
    CALL my_concat(mq(lmsq+1),msqf,msglen,parai%allgrp)
    ! TRANSLATE I,J TO A SINGLE G/S INDEX
    DO ip=0,parai%nproc-1
       mxrp=parap%sparm(8,ip)
       DO ixf=1,mxrp
          jj=ixf+ip*lmsq
          i=msqs(jj)
          j=msqf(jj)
          msqfpool(ixf,ip+1,ipool)=i+(j-1)*qr2s
          IF (ixf.LE.parap%sparm(9,ip))&
               msqspool(ixf,ip+1,ipool)=i+(j-qr3min)*qr2s
       ENDDO
    ENDDO
    ! REDEFINE NZH AND INDZ FOR COMPRESSED STORAGE
    !$omp parallel do private(IG,NY1,NY2,NY3,INY1,INY2,INY3)
    DO ig=1,jhg
       ny1=inzh(1,ig)
       ny2=inzh(2,ig)
       ny3=inzh(3,ig)
       iny1=-ny1+2*nh1
       iny2=-ny2+2*nh2
       iny3=-ny3+2*nh3
       nzff(ig)=ny1 + (mg(ny2,ny3)-1)*qr1s
       inzf(ig)=iny1 + (mg(iny2,iny3)-1)*qr1s
    ENDDO
    !$omp parallel do private(IG)
    DO ig=1,jgw
       nzfs(ig)=nzff(ig)
       inzs(ig)=inzf(ig)
    ENDDO
    DEALLOCATE(mg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(my,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mz,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mq,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE setrays
  ! ==================================================================
  SUBROUTINE Prep_fft_com( comm_send, comm_recv, sendsize, sendsize_rem, nodes_numb, mype, my_node, my_node_rank, node_task_size, buffer_size, comm_sendrecv, do_comm, WAVE )
    IMPLICIT NONE
  
    INTEGER, INTENT(IN)                                 :: sendsize, sendsize_rem, nodes_numb, mype, my_node, my_node_rank, node_task_size, buffer_size
    COMPLEX(DP), INTENT(IN)                             :: comm_send( : , : )
    COMPLEX(DP), INTENT(INOUT)                          :: comm_recv( : , : )
    INTEGER, INTENT(OUT)                                :: comm_sendrecv( : )
    LOGICAL, INTENT(OUT)                                :: do_comm
    INTEGER, INTENT(IN)                                 :: WAVE
  
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
  
!    IF( WAVE ) THEN
!       IF( first ) THEN
!          first = .false.
!          buffer_size_save = buffer_size
!       ELSE
!          DO i = 1, buffer_size_save
!             DO j = 1, comm_sendrecv(1)
!                CALL MPI_REQUEST_FREE( parai%send_handle( j , i , 1 ) )
!                CALL MPI_REQUEST_FREE( parai%send_handle( j , i , 2 ) )
!             ENDDO
!             DO j = 1, comm_sendrecv(2)
!                CALL MPI_REQUEST_FREE( parai%recv_handle( j , i , 1 ) )
!                CALL MPI_REQUEST_FREE( parai%recv_handle( j , i , 2 ) )
!             ENDDO
!          ENDDO
!          buffer_size_save = buffer_size
!       END IF
!    END IF
    
    IF( ALLOCATED( parai%send_handle ) .and. WAVE .eq. 1 )       DEALLOCATE( parai%send_handle )
    IF( ALLOCATED( parai%recv_handle ) .and. WAVE .eq. 1 )       DEALLOCATE( parai%recv_handle )
    
    IF( .not. ALLOCATED(parai%send_handle) ) ALLOCATE( parai%send_handle( comm_sendrecv(1), buffer_size, 2, 2 ) )
    IF( .not. ALLOCATED(parai%recv_handle) ) ALLOCATE( parai%recv_handle( comm_sendrecv(2), buffer_size, 2, 2 ) )
    
    DO i = 1, buffer_size !INITIALIZE SENDING AND RECEIVING
  
       DO j = 1, comm_sendrecv( 1 )
  
          target_node = (comm_info_send(j,mype+1)-1) / node_task_size
  
          CALL mp_send_init_complex( comm_send(:,i), target_node*sendsize, sendsize, comm_info_send( j , mype+1 ) - 1, mype, parai%allgrp, parai%send_handle( j , i, 1, WAVE ) )
  
       ENDDO        
  
       DO j = 1, comm_sendrecv( 2 )
  
          origin_node = (comm_info_recv(j,mype+1)-1) / node_task_size
  
          CALL mp_recv_init_complex( comm_recv(:,i), origin_node*sendsize, sendsize, comm_info_recv( j , mype+1 ) - 1, parai%allgrp, parai%recv_handle( j , i, 1, WAVE ) )
  
       ENDDO        
  
    ENDDO
  
    IF( sendsize_rem .ne. 0 ) THEN
  
       DO i = 1, buffer_size !INITIALIZE SENDING AND RECEIVING
    
          DO j = 1, comm_sendrecv( 1 )
  
             target_node = (comm_info_send(j,mype+1)-1) / node_task_size
    
             CALL mp_send_init_complex( comm_send(:,i), target_node*sendsize_rem, sendsize_rem, comm_info_send( j , mype+1 ) - 1, mype, parai%allgrp, parai%send_handle( j , i, 2, WAVE ) )
   
          ENDDO        
    
          DO j = 1, comm_sendrecv( 2 )
    
             origin_node = (comm_info_recv(j,mype+1)-1) / node_task_size
    
             CALL mp_recv_init_complex( comm_recv(:,i), origin_node*sendsize_rem, sendsize_rem, comm_info_recv( j , mype+1 ) - 1, parai%allgrp, parai%recv_handle( j , i, 2, WAVE ) )
   
          ENDDO        
    
       ENDDO
  
    END IF
  
  END SUBROUTINE Prep_fft_com
  ! ==================================================================
  SUBROUTINE Make_Manual_Maps( tfft, batch_size, rem_size, nss, nr1s, ngs, which )
    IMPLICIT NONE
  
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    INTEGER, INTENT(IN) :: batch_size, rem_size
    INTEGER, INTENT(IN) :: nss(:), nr1s, ngs

    CHARACTER(*), PARAMETER                     :: procedureN = 'Prep_fft_com'
  
    INTEGER :: i, j, overlap_cor, eff_nthreads, ierr, which

    ierr = 0
  
    IF( cntl%overlapp_comm_comp .and. parai%ncpus .gt. 1 .and. tfft%do_comm(which) .and. which .eq. 1 ) THEN
       eff_nthreads = parai%ncpus - 1
    ELSE
       eff_nthreads = parai%ncpus
    END IF
    overlap_cor = parai%ncpus - eff_nthreads
  
  ! z things
    
    DO j = 1, parai%node_nproc * parai%nnode
  
       DO i = 1+overlap_cor, parai%ncpus
          tfft%thread_z_sticks( i, 1, j, which ) = ( nss( j ) * batch_size ) / eff_nthreads
       ENDDO
       DO i = 1+overlap_cor, mod( nss( j ) * batch_size, eff_nthreads ) + overlap_cor
          tfft%thread_z_sticks( i, 1, j, which ) = tfft%thread_z_sticks( i, 1, j, which ) + 1
       ENDDO
     
       tfft%thread_z_start( 1+overlap_cor, 1, j, which ) = 1
       DO i = 2+overlap_cor, parai%ncpus
          tfft%thread_z_start( i, 1, j, which ) = tfft%thread_z_start( i-1, 1, j, which ) + tfft%thread_z_sticks( i-1, 1, j, which )
       ENDDO
       DO i = 1+overlap_cor, parai%ncpus
          tfft%thread_z_end( i, 1, j, which ) = tfft%thread_z_start( i, 1, j, which ) + tfft%thread_z_sticks( i, 1, j, which ) - 1
       ENDDO
     
       IF( rem_size .ne. 0 ) THEN
          DO i = 1+overlap_cor, parai%ncpus
             tfft%thread_z_sticks( i, 2, j, which ) = ( nss( j ) * rem_size ) / eff_nthreads
          ENDDO
          DO i = 1+overlap_cor, mod( nss( j ) * rem_size, eff_nthreads ) + overlap_cor
             tfft%thread_z_sticks( i, 2, j, which ) = tfft%thread_z_sticks( i, 2, j, which ) + 1
          ENDDO
        
          tfft%thread_z_start( 1+overlap_cor, 2, j, which ) = 1
          DO i = 2+overlap_cor, parai%ncpus
             tfft%thread_z_start( i, 2, j, which ) = tfft%thread_z_start( i-1, 2, j, which ) + tfft%thread_z_sticks( i-1, 2, j, which )
          ENDDO
          DO i = 1+overlap_cor, parai%ncpus
             tfft%thread_z_end( i, 2, j, which ) = tfft%thread_z_start( i, 2, j, which ) + tfft%thread_z_sticks( i, 2, j, which ) - 1
          ENDDO
       END IF
  
    ENDDO
  
  ! y things
  
    DO i = 1+overlap_cor, parai%ncpus
       tfft%thread_y_sticks( i, 1, which ) = ( nr1s * tfft%my_nr3p * batch_size ) / eff_nthreads
    ENDDO
    DO i = 1+overlap_cor, mod( nr1s * tfft%my_nr3p * batch_size, eff_nthreads ) + overlap_cor
       tfft%thread_y_sticks( i, 1, which ) = tfft%thread_y_sticks( i, 1, which ) + 1
    ENDDO
  
    tfft%thread_y_start( 1+overlap_cor, 1, which ) = 1
    DO i = 2+overlap_cor, parai%ncpus
       tfft%thread_y_start( i, 1, which ) = tfft%thread_y_start( i-1, 1, which ) + tfft%thread_y_sticks( i-1, 1, which )
    ENDDO
    DO i = 1+overlap_cor, parai%ncpus
       tfft%thread_y_end( i, 1, which ) = tfft%thread_y_start( i, 1, which ) + tfft%thread_y_sticks( i, 1, which ) - 1
    ENDDO
  
    IF( rem_size .ne. 0 ) THEN
       DO i = 1+overlap_cor, parai%ncpus
          tfft%thread_y_sticks( i, 2, which ) = ( nr1s * tfft%my_nr3p * rem_size ) / eff_nthreads
       ENDDO
       DO i = 1+overlap_cor, mod( nr1s * tfft%my_nr3p * rem_size, eff_nthreads ) + overlap_cor
          tfft%thread_y_sticks( i, 2, which ) = tfft%thread_y_sticks( i, 2, which ) + 1
       ENDDO
     
       tfft%thread_y_start( 1+overlap_cor, 2, which ) = 1
       DO i = 2+overlap_cor, parai%ncpus
          tfft%thread_y_start( i, 2, which ) = tfft%thread_y_start( i-1, 2, which ) + tfft%thread_y_sticks( i-1, 2, which )
       ENDDO
       DO i = 1+overlap_cor, parai%ncpus
          tfft%thread_y_end( i, 2, which ) = tfft%thread_y_start( i, 2, which ) + tfft%thread_y_sticks( i, 2, which ) - 1
       ENDDO
    END IF
  
  ! x things
    
    DO i = 1+overlap_cor, parai%ncpus
       tfft%thread_x_sticks( i, 1, which ) = ( tfft%my_nr3p * tfft%nr2 * batch_size ) / eff_nthreads
    ENDDO
    DO i = 1+overlap_cor, mod( tfft%my_nr3p * tfft%nr2 * batch_size, eff_nthreads ) + overlap_cor
       tfft%thread_x_sticks( i, 1, which ) = tfft%thread_x_sticks( i, 1, which ) + 1
    ENDDO
  
    tfft%thread_x_start( 1+overlap_cor, 1, which ) = 1
    DO i = 2+overlap_cor, parai%ncpus
       tfft%thread_x_start( i, 1, which ) = tfft%thread_x_start( i-1, 1, which ) + tfft%thread_x_sticks( i-1, 1, which )
    ENDDO
    DO i = 1+overlap_cor, parai%ncpus
       tfft%thread_x_end( i, 1, which ) = tfft%thread_x_start( i, 1, which ) + tfft%thread_x_sticks( i, 1, which ) - 1
    ENDDO
  
    IF( rem_size .ne. 0 ) THEN
       DO i = 1+overlap_cor, parai%ncpus
          tfft%thread_x_sticks( i, 2, which ) = ( tfft%my_nr3p * tfft%nr2 * rem_size ) / eff_nthreads
       ENDDO
       DO i = 1+overlap_cor, mod( tfft%my_nr3p * tfft%nr2 * rem_size, eff_nthreads ) + overlap_cor
          tfft%thread_x_sticks( i, 2, which ) = tfft%thread_x_sticks( i, 2, which ) + 1
       ENDDO
     
       tfft%thread_x_start( 1+overlap_cor, 2, which ) = 1
       DO i = 2+overlap_cor, parai%ncpus
          tfft%thread_x_start( i, 2, which ) = tfft%thread_x_start( i-1, 2, which ) + tfft%thread_x_sticks( i-1, 2, which )
       ENDDO
       DO i = 1+overlap_cor, parai%ncpus
          tfft%thread_x_end( i, 2, which ) = tfft%thread_x_start( i, 2, which ) + tfft%thread_x_sticks( i, 2, which ) - 1
       ENDDO
    END IF
  
  ! gspace things
  
    DO i = 1+overlap_cor, parai%ncpus
       tfft%thread_ngms( i ) = ( ngs ) / eff_nthreads
    ENDDO
    DO i = 1+overlap_cor, mod( ngs, eff_nthreads ) + overlap_cor
       tfft%thread_ngms( i ) = tfft%thread_ngms( i ) + 1
    ENDDO
  
    tfft%thread_ngms_start( 1+overlap_cor ) = 1
    DO i = 2+overlap_cor, parai%ncpus
       tfft%thread_ngms_start( i ) = tfft%thread_ngms_start( i-1 ) + tfft%thread_ngms( i-1 )
    ENDDO
    DO i = 1+overlap_cor, parai%ncpus
       tfft%thread_ngms_end( i ) = tfft%thread_ngms_start( i ) + tfft%thread_ngms( i ) - 1
    ENDDO
  
  ! rspace things
    
    DO i = 1+overlap_cor, parai%ncpus
       tfft%thread_rspace( i ) = ( tfft%my_nr3p * tfft%nr2 * tfft%nr1 ) / eff_nthreads
    ENDDO
    DO i = 1+overlap_cor, mod( tfft%my_nr3p * tfft%nr2 * tfft%nr1, eff_nthreads ) + overlap_cor
       tfft%thread_rspace( i ) = tfft%thread_rspace( i ) + 1
    ENDDO
  
    tfft%thread_rspace_start( 1+overlap_cor ) = 1
    DO i = 2+overlap_cor, parai%ncpus
       tfft%thread_rspace_start( i ) = tfft%thread_rspace_start( i-1 ) + tfft%thread_rspace( i-1 )
    ENDDO
    DO i = 1+overlap_cor, parai%ncpus
       tfft%thread_rspace_end( i ) = tfft%thread_rspace_start( i ) + tfft%thread_rspace( i ) - 1
    ENDDO
    
  END SUBROUTINE Make_Manual_Maps
  ! ==================================================================
  SUBROUTINE Make_inv_yzCOM_Maps( tfft, map_acinv, batch_size, ir1s, nss, my_nr1s, small_chunks, big_chunks, zero_start, zero_end )
    IMPLICIT NONE
  
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    INTEGER, INTENT(IN)  :: batch_size, my_nr1s, small_chunks, big_chunks
    INTEGER, INTENT(OUT) :: map_acinv(:)
    INTEGER, INTENT(IN)  :: ir1s(:), nss(:)
    INTEGER, OPTIONAL, INTENT(OUT) :: zero_start(:), zero_end(:)
  
    LOGICAL :: l_map( tfft%my_nr3p * my_nr1s * tfft%nr2 )
    LOGICAL :: first
    INTEGER :: ibatch, j, l, i, k
    INTEGER :: iproc, offset, it, mc, m1, m2, i1
    INTEGER :: ierr
  
    l_map = .false.
    map_acinv = 0
    
    !$omp parallel private( ibatch, j, l, iproc, offset, i, it, mc, m1, m2, i1, k)
    DO ibatch = 1, batch_size
       DO j = 1, parai%nnode
          DO l = 1, parai%node_nproc
             iproc = (j-1)*parai%node_nproc + l
             offset = ( parai%node_me*parai%node_nproc + (l-1) ) * small_chunks + ( (j-1)*batch_size + (ibatch-1) ) * big_chunks
             !$omp do
             DO i = 1, nss( iproc )
                it = offset + tfft%nr3px * (i-1) 
                mc = tfft%ismap( i + tfft%iss(iproc) ) ! this is  m1+(m2-1)*nr1x  of the  current pencil
                m1 = mod ( mc-1, tfft%nr1 ) + 1
                m2 = (mc-1)/tfft%nr1 + 1
                i1 = m2 + ( ir1s(m1) - 1 ) * tfft%nr2 + (ibatch-1)*tfft%my_nr3p*my_nr1s*tfft%nr2
                DO k = 1, tfft%my_nr3p
                   IF( ibatch .eq. 1 ) l_map( i1 ) = .true.
                   map_acinv( i1 ) = k + it
                   i1 = i1 + tfft%nr2*my_nr1s
                ENDDO
             ENDDO
             !$omp end do
          ENDDO
       ENDDO
    ENDDO
    !$omp end parallel
  
    IF( present( zero_start ) ) THEN
     
       zero_start = 0
       zero_end = tfft%nr2
       first = .true.
       
!       IF( parai%me .eq. 0 ) THEN 
          DO j = 1, my_nr1s
             first = .true.
             DO l = 1, tfft%nr2
       
                IF( l_map( (j-1)*tfft%nr2 + l ) .eqv. .true. ) THEN
                   IF( first .eqv. .false. ) THEN
                      zero_end( j ) = l-1
                      first = .true.
                   END IF
                ELSE
                   IF( first .eqv. .true. ) THEN
                      zero_start( j ) = l
                      first = .false.
                   END IF
                END IF
       
             ENDDO
          ENDDO
!       END IF
!       IF( parai%node_me .eq. 0 .and. parai%nnode .ne. 1 ) THEN
!          CALL MP_BCAST( zero_start, my_nr1s, 0, tfft%inter_node_comm )
!          CALL MP_BCAST( zero_end  , my_nr1s, 0, tfft%inter_node_comm )
!       END IF
!       CALL MP_BCAST( zero_start, my_nr1s, 0, tfft%node_comm )
!       CALL MP_BCAST( zero_end  , my_nr1s, 0, tfft%node_comm )
  
    END IF
    
  END SUBROUTINE Make_inv_yzCOM_Maps
  ! ==================================================================

END MODULE fftnew_utils
