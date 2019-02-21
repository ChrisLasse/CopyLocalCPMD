#include "cpmd_global.h"

MODULE noforce_utils
  USE csize_utils,                     ONLY: csize
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: blocking_fft
  USE fnonloc_utils,                   ONLY: fnonloc,&
                                             give_scr_fnonloc
  USE geq0mod,                         ONLY: geq0
  USE gsize_utils,                     ONLY: gsize
  USE hnlmat_utils,                    ONLY: hnlmat
  USE hubbardu,                        ONLY: c2u0
  USE hubbardu_utils,                  ONLY: add_hubbardu
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlforce_utils,                   ONLY: give_scr_nlforce,&
                                             nlforce
  USE nlps,                            ONLY: imagp
  USE nlsl_utils,                      ONLY: nlsl
  USE norm,                            ONLY: cnorm,&
                                             gemax,&
                                             gnmax,&
                                             gnorm
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE puttau_utils,                    ONLY: taucl
  USE reigs_utils,                     ONLY: reigs
  USE rgsvan_utils,                    ONLY: rgsvan
  USE rnlfl_utils,                     ONLY: rnlfl
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE ropt,                            ONLY: ropt_mod
  USE rotate_utils,                    ONLY: rotate,&
                                             rottr
  USE rscpot_utils,                    ONLY: give_scr_rscpot,&
                                             rscpot
  USE rswfmod,                         ONLY: rsactive
  USE spin,                            ONLY: clsd,&
                                             lspin2,&
                                             spin_mod
  USE summat_utils,                    ONLY: summat
  USE symtrz_utils,                    ONLY: symvec
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean
  USE vpsi_utils,                      ONLY: vpsi,vpsi_blocking
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: noforce
  PUBLIC :: give_scr_noforce
CONTAINS

  ! ==================================================================
  SUBROUTINE noforce(c0,c2,sc0,tau0,fion,eigv,rhoe,psi,nstate,tfor)
    ! ==--------------------------------------------------------------==
    ! ==          COMPUTES FOR  A SET OF NONORTHOGONAL ORBITALS       ==
    ! ==                     THE TOTAL ENERGY                         ==
    ! ==                  THE ELECTRONIC FORCES                       ==
    ! ==                  THE FORCES ON THE IONS                      ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:), sc0(:,:)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:), &
                                                eigv(*), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    INTEGER                                  :: nstate
    LOGICAL                                  :: tfor

    CHARACTER(*), PARAMETER                  :: procedureN = 'noforce'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8)                          :: pab
    INTEGER                                  :: i, iabsl, ierr,  &
                                                il_gam(1), il_smat(2), &
                                                isub, lnoforce
    INTEGER, EXTERNAL                        :: izamax
    LOGICAL                                  :: debug
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: gam(:), smat(:,:)
#else
    REAL(real_8), ALLOCATABLE                :: gam(:), smat(:,:)
#endif
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    debug=.FALSE.
    IF (lspin2%tlse) CALL stopgm('NOFORCE','NO LSE ALLOWED HERE',& 
         __LINE__,__FILE__)
    pab = 0._real_8
    IF (imagp.EQ.2) CALL stopgm('NOFORCE','K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (geq0) CALL zclean(c0,nstate,ncpw%ngw)
    gnmax=0.0_real_8
    gnorm=0.0_real_8
    ! ==--------------------------------------------------------------==
    il_gam=imagp*nstate*nstate    
    il_smat(1)=nstate
    il_smat(2)=nstate
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_gam,gam,procedureN//'_gam')
    CALL request_scratch(il_smat,smat,procedureN//'_smat')
#else
    ALLOCATE(gam(il_gam(1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(smat(il_smat(1),il_smat(2)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
#endif
    CALL dcopy(2*ncpw%ngw*nstate,c0,1,sc0,1)
    ! ==--------------------------------------------------------------==
    ! ==   CALCULATE INVERSE CHOLESKY DECOMPOSITION                   ==
    ! ==   SCHMIDT ORTHOGONALIZATION    SC0 = U**(-1)*C0              ==
    ! ==--------------------------------------------------------------==
    call rgsvan(sc0,nstate,smat,.TRUE.,need_fnl=.TRUE.,rot_fnl=.TRUE.)
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    ! ==   CALCULATE THE POTENTIAL AND THE FORCE ON THE IONS          ==
    ! ==--------------------------------------------------------------==
    IF (pslo_com%tivan) THEN
       CALL rnlsm(sc0,nstate,1,1,tfor,rdst=.FALSE.,only_dfnl=.TRUE.,packed_dfnl=.TRUE.)
    ELSE
       CALL rnlsm(sc0,nstate,1,1,tfor)
    END IF
    rsactive = cntl%krwfn
    CALL rscpot(sc0,tau0,fion,rhoe,psi,tfor,ropt_mod%calste,nstate,1)
    ! ==--------------------------------------------------------------==
    ! ==   CALCULATE THE ELECTRONIC FORCE                             ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(c2)!,ngw*nstate)
    ! ==--------------------------------------------------------------==
    ! == Compute the force on the electronic degrees of freedom due   ==
    ! == to the local potential (stored in RHOE)                      ==
    ! ==--------------------------------------------------------------==

    IF (blocking_fft) THEN
       CALL vpsi_blocking(sc0,c2,crge%f(:,1),rhoe,psi(:,1),nstate,1,clsd%nlsd,.TRUE.)
    ELSE
       CALL vpsi(sc0,c2,crge%f(:,1),rhoe,psi(:,1),nstate,1,clsd%nlsd,.TRUE.)
    END IF
    ! c2u0 is calculated in uprho or rscpot
    IF (cntl%thubb) CALL add_hubbardu(c2,c2u0,nstate)
    IF (pslo_com%tivan) THEN
       CALL ovlap(nstate,gam,c2,sc0,rdst=.FALSE.,full=.FALSE.)
       CALL hnlmat(gam,crge%f,nstate)
       CALL summat(gam,nstate,.TRUE.,.TRUE.,parai%cp_grp)
       CALL nlforce(c2,crge%f,gam,nstate)
       CALL rotate(-1.0_real_8,sc0,1.0_real_8,c2,gam,nstate,2*ncpw%ngw,cntl%tlsd,&
            spin_mod%nsup,spin_mod%nsdown)
       IF (tfor) CALL rnlfl(fion,gam,nstate,1)
       IF (ropt_mod%calste) CALL nlsl(gam,nstate)
    ELSE
       ! ==--------------------------------------------------------------==
       ! == Compute the force on the electronic degrees of freedom due   ==
       ! == to the non-local part of the potential, and add it to the    ==
       ! == other piece, coming from the local contribution.             ==
       ! ==--------------------------------------------------------------==
       CALL fnonloc(c2,crge%f,nstate,1,clsd%nlsd,.TRUE.)
       CALL ovlap(nstate,gam,c2,sc0,rdst=.FALSE.,full=.FALSE.)
       CALL summat(gam,nstate,.TRUE.,.TRUE.,parai%cp_grp)
       ! ==--------------------------------------------------------------==
       ! ==   C2(I) = C2(I) - SUM(J) <SC(I) | H | SC(J)> SC(J)           ==
       ! ==--------------------------------------------------------------==
       CALL rotate(-1.0_real_8,sc0,1.0_real_8,c2,gam,nstate,2*ncpw%ngw,cntl%tlsd,&
            spin_mod%nsup,spin_mod%nsdown)
    ENDIF
    IF (ropt_mod%prteig.AND.paral%parent) THEN
       CALL dscal(nstate*nstate,-1.0_real_8,gam,1)
       CALL reigs(nstate,gam,crge%f)
       CALL dscal(nstate*nstate,-1.0_real_8,gam,1)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==   ROTATE ELECTRONIC FORCE BACK INTO NONORTHOGONAL BASIS      ==
    ! ==--------------------------------------------------------------==
    CALL rottr(1._real_8,c2,smat,"T",nstate,ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    IF (geq0) CALL zclean(c2,nstate,ncpw%ngw)
    ! ==--------------------------------------------------------------==
    CALL csize(c2,nstate,gemax,cnorm)
    ! ==--------------------------------------------------------------==
    rsactive = .FALSE.
    IF (tfor) THEN
       CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%allgrp)
       IF (paral%parent) THEN
          CALL symvec(fion)
          CALL taucl(fion)
          CALL gsize(fion,gnmax,gnorm)
       ENDIF
    ENDIF
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_smat,smat,procedureN//'_smat')
    CALL free_scratch(il_gam,gam,procedureN//'_gam')
#else
    DEALLOCATE(gam,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(smat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#endif
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE noforce
  ! ==================================================================
  SUBROUTINE give_scr_noforce(lnoforce,il_gam,il_auxc,il_smat,&
       il_ddia,tag,nstate,tfor)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lnoforce, il_gam, il_auxc, &
                                                il_smat, il_ddia
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: tfor

    INTEGER                                  :: lrnlsm, lrscpot

    il_smat=0
    IF (pslo_com%tivan) THEN
       CALL give_scr_nlforce(il_gam,il_auxc,il_ddia,nstate)
    ELSE
       CALL give_scr_fnonloc(il_auxc,il_ddia,nstate)
       il_gam=imagp*nstate*nstate
    ENDIF
    IF (ropt_mod%prteig.AND.paral%parent) THEN
       il_auxc=MAX(il_auxc,nstate*(nstate+1)/2+3*nstate)! REIGS
    ENDIF
    il_auxc=MAX(il_auxc,nstate*nstate)
    il_smat=MAX(il_smat,nstate*nstate)
    il_ddia=MAX(il_ddia,2*maxsys%nax*nstate)
    CALL give_scr_rnlsm(lrnlsm,tag,nstate,tfor)
    CALL give_scr_rscpot(lrscpot,tag,ropt_mod%calste)
    lnoforce=il_gam+il_auxc+il_smat+il_ddia+100
    lnoforce=MAX(lnoforce,lrnlsm,lrscpot)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_noforce

END MODULE noforce_utils
