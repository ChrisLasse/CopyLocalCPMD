MODULE vtaupsi_utils
  USE cppt,                            ONLY: gk,&
                                             indzs,&
                                             nzhs
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE parac,                           ONLY: parai
  USE spin,                            ONLY: lspin2,&
                                             spin_mod
  USE system,                          ONLY: fpar,&
                                             ncpw,&
                                             nkpt,&
                                             parm
  USE tauf,                            ONLY: itaur,&
                                             vtau
  USE tauofr_utils,                    ONLY: dpsisc
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vtaupsi
  PUBLIC :: taupot
  PUBLIC :: ftauadd

CONTAINS

  ! ==================================================================
  SUBROUTINE vtaupsi(c0,c2,f,psi,nstate,ispin)
    ! ==================================================================
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate)
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,nstate), &
                                                c0(nkpt%ngwk,nstate)
    INTEGER                                  :: ispin

    INTEGER                                  :: is, is1, is2, isub
    REAL(real_8), POINTER                    :: vpot(:,:)

    CALL tiset('   VTAUPSI',isub)
    IF (tkpts%tkpnt) CALL stopgm("VTAUPSI","KPNT NOT IMPLEMENTED",& 
         __LINE__,__FILE__)
    IF (lspin2%tlse) CALL stopgm("VTAUPSI","LSE NOT IMPLEMENTED",& 
         __LINE__,__FILE__)
    vpot => vtau(:, itaur:)
    ! ==--------------------------------------------------------------==
    ! Loop over the electronic states
    DO is=1,nstate,2
       is1=is
       is2=is+1
       ! x direction
       CALL zeroing(psi)!,maxfft)
       CALL dpsisc(c0,psi,1,is1,is2,nstate)
       CALL  invfftn(psi,.TRUE.,parai%allgrp)
       CALL taupot(vpot,psi,is1,is2,ispin)
       CALL  fwfftn(psi,.TRUE.,parai%allgrp)
       CALL ftauadd(c2,psi,f,1,is1,is2,nstate)
       ! y direction
       CALL zeroing(psi)!,maxfft)
       CALL dpsisc(c0,psi,2,is1,is2,nstate)
       CALL  invfftn(psi,.TRUE.,parai%allgrp)
       CALL taupot(vpot,psi,is1,is2,ispin)
       CALL  fwfftn(psi,.TRUE.,parai%allgrp)
       CALL ftauadd(c2,psi,f,2,is1,is2,nstate)
       ! z direction
       CALL zeroing(psi)!,maxfft)
       CALL dpsisc(c0,psi,3,is1,is2,nstate)
       CALL  invfftn(psi,.TRUE.,parai%allgrp)
       CALL taupot(vpot,psi,is1,is2,ispin)
       CALL  fwfftn(psi,.TRUE.,parai%allgrp)
       CALL ftauadd(c2,psi,f,3,is1,is2,nstate)
    ENDDO
    CALL tihalt('   VTAUPSI',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vtaupsi
  ! ==================================================================
  SUBROUTINE taupot(vpot,psi,is1,is2,ispin)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: is1, is2, ispin
    REAL(real_8)                             :: vpot(fpar%nnr1,ispin)

    INTEGER                                  :: ir
    REAL(real_8)                             :: r1, r2

    IF (ispin.EQ.1) THEN
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          psi(ir)=psi(ir)*vpot(ir,1)
       ENDDO
    ELSE
       IF (is2.LE.spin_mod%nsup) THEN
          !$omp parallel do private(IR)
          DO ir=1,fpar%nnr1
             psi(ir)=psi(ir)*vpot(ir,1)
          ENDDO
       ELSEIF (is1.GT.spin_mod%nsup) THEN
          !$omp parallel do private(IR)
          DO ir=1,fpar%nnr1
             psi(ir)=psi(ir)*vpot(ir,2)
          ENDDO
       ELSE
          !$omp parallel do private(IR,R1,R2)
          DO ir=1,fpar%nnr1
             r1=AIMAG(psi(ir))
             r2=REAL(psi(ir))
             psi(ir)=CMPLX(r2*vpot(ir,2),r1*vpot(ir,1),kind=real_8)
          ENDDO
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE taupot
  ! ==================================================================
  SUBROUTINE ftauadd(c2,psi,f,k,is1,is2,nstate)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: f(*)
    INTEGER                                  :: k, is1, is2, nstate
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,nstate)

    COMPLEX(real_8)                          :: fm, fp
    INTEGER                                  :: ig
    REAL(real_8)                             :: fi1, fi2

    IF (is2.GT.nstate) THEN
       fi1=0.25_real_8*f(is1)*parm%tpiba2
       !$omp parallel do private(IG,FP,FM)
       DO ig=1,ncpw%ngw
          fp=psi(nzhs(ig))+psi(indzs(ig))
          fm=psi(nzhs(ig))-psi(indzs(ig))
          c2(ig,is1)=c2(ig,is1)-&
               fi1*gk(k,ig)*CMPLX(REAL(fm),AIMAG(fp),kind=real_8)
       ENDDO
    ELSE
       fi1=0.25_real_8*f(is1)*parm%tpiba2
       fi2=0.25_real_8*f(is2)*parm%tpiba2
       !$omp parallel do private(IG,FP,FM)
       DO ig=1,ncpw%ngw
          fp=psi(nzhs(ig))+psi(indzs(ig))
          fm=psi(nzhs(ig))-psi(indzs(ig))
          c2(ig,is1)=c2(ig,is1)-&
               fi1*gk(k,ig)*CMPLX(REAL(fm),AIMAG(fp),kind=real_8)
          c2(ig,is2)=c2(ig,is2)-&
               fi2*gk(k,ig)*CMPLX(AIMAG(fm),-REAL(fp),kind=real_8)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ftauadd
  ! ==================================================================

END MODULE vtaupsi_utils
