#include "cpmd_global.h"

MODULE newd_utils
  USE cppt,                            ONLY: gk,&
                                             inyh
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms
  USE distribution_utils,              ONLY: dist_size
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: nlps_com
  !$ USE omp_lib,                       ONLY: omp_get_thread_num
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE qvan2_utils,                     ONLY: qvan2
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE system,                          ONLY: cntl,&
                                             iatpe,&
                                             ipept,&
                                             maxsys,&
                                             ncpw,&
                                             parm,&
                                             paraw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: newd
  !public :: cftemp
  PUBLIC :: give_scr_newd

CONTAINS

    ! ==--------------------------------------------------------------==
    ! ==         THIS ROUTINE CALCULATES THE ARRAY DEEQ               ==
    ! ==                                                              ==
    ! ==  DEEQ_I,IJ = OMEGA( V(G=0) Q_I,IJ(G=0) +                     ==
    ! ==       2 SUM_G> RE[ V*(G) Q_I,IJ(G) E^-IG.R_I ] )             ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
  SUBROUTINE newd(fnla,deeq,f,vpot,fion,nstate,tfor)
    REAL(real_8)                             :: fnla(:,:,:), &
                                                deeq(ions1%nat,maxsys%nhxs,*)
    COMPLEX(real_8)                          :: vpot(ncpw%nhg)
    REAL(real_8)                             :: fion(:,:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate)
    LOGICAL                                  :: tfor
    INTEGER                                  :: i, ii, j, is, ia, isa, isa0, ijv, iv, jv, &
                                                ig, k, nhh, num_orb, methread, &
                                                na_grp(2,ions1%nsp,parai%cp_nogrp), &
                                                il_fnlat(3), il_gktemp(2), il_vtmp(3), &
                                                il_ylm1(2), il_ylm2(3), il_qg(2), ierr, isub
    REAL(real_8)                             :: ftmp, otr, rhovan, tom, omtpiba, fiont(3),tmp
    REAL(real_8), EXTERNAL                   :: dotp
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: fnlat(:,:,:),gktemp(:,:),ylm1(:,:),&
                                                ylm2(:,:,:)
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: qg(:,:), vtmp(:,:,:)    
#else
    REAL(real_8), ALLOCATABLE                :: fnlat(:,:,:),gktemp(:,:),ylm1(:,:),&
                                                ylm2(:,:,:)
    COMPLEX(real_8), ALLOCATABLE             :: qg(:,:), vtmp(:,:,:)    
#endif
#ifdef _VERBOSE_FORCE_DBG
    REAL(real_8),ALLOCATABLE                 :: dbg_forces(:,:,:)
#endif
    CHARACTER(*), PARAMETER                  :: procedureN='newd'

    CALL tiset(procedureN,isub)
    !generate 'new' nst12 mapping -> nwa12
    CALL dist_size(nstate,parai%nproc,paraw%nwa12,nblock=1,nbmax=num_orb,fw=1)
    num_orb=paraw%nwa12(2,parai%mepos)-paraw%nwa12(1,parai%mepos)+1
    !no need to shift, because fnla is a pointer to the right states of fnl
    
    tom=2.0_real_8*parm%omega
    IF (cntl%bigmem) THEN
       IF (tfor) THEN
          IF (parai%cp_nogrp.GT.1 .AND. parai%cp_inter_me .GT. 0) THEN
             !$omp parallel do private(is,ia)
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   fion(1:3,ia,is)=0._real_8
                END DO
             END DO
          END IF
          OMTPIBA=parm%omega*parm%tpiba

          il_gktemp(1)=ncpw%nhg_cp
          il_gktemp(2)=3
#ifdef _USE_SCRATCHLIBRARY
          CALL request_scratch(il_gktemp,gktemp,procedureN//'_gktemp')
#else
          ALLOCATE(gktemp(il_gktemp(1),il_gktemp(2)), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate gktemp)',& 
               __LINE__,__FILE__)
#endif
          !$OMP do private(ig,k) collapse(2)
          DO ig=1,ncpw%nhg_cp
             DO k=1,3
                gktemp(ig,k)=gk(k,ncpw%nhg_start-1+ig)
             ENDDO
          ENDDO
          CALL cp_grp_split_atoms(na_grp)
       END IF

       isa0=0
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is)) THEN
             nhh=(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2
             il_qg(1)=ncpw%nhg_cp
             il_qg(2)=nhh
             il_ylm1(1)=ions0%na(is)
             il_ylm1(2)=ncpw%nhg_cp
#ifdef _USE_SCRATCHLIBRARY
             CALL request_scratch(il_qg,qg,procedureN//'_qg')
             CALL request_scratch(il_ylm1,ylm1,procedureN//'_ylm1')            
#else
             ALLOCATE(qg(il_qg(1),il_qg(2)), stat=ierr)
             IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate qg)',& 
                  __LINE__,__FILE__)
             ALLOCATE(ylm1(il_ylm1(1),il_ylm1(2)), stat=ierr)
             IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ylm1)',& 
                  __LINE__,__FILE__)
#endif
             IF (tfor) THEN
                il_ylm2(1)=3
                il_ylm2(2)=nhh
                il_ylm2(3)=ions0%na(is)
                il_vtmp(1)=ncpw%nhg_cp
                il_vtmp(2)=3
                il_vtmp(3)=nhh
                il_fnlat(1)=num_orb
                il_fnlat(2)=nlps_com%ngh(is)
                il_fnlat(3)=parai%ncpus
#ifdef _USE_SCRATCHLIBRARY
                CALL request_scratch(il_ylm2,ylm2,procedureN//'_ylm2')
                CALL request_scratch(il_vtmp,vtmp,procedureN//'_vtmp')
                CALL request_scratch(il_fnlat,fnlat,procedureN//'_fnlat')
#else
             ALLOCATE(ylm2(il_ylm2(1),il_ylm2(2),il_ylm2(3)), stat=ierr)
             IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ylm2)',& 
                  __LINE__,__FILE__)
             ALLOCATE(vtmp(il_vtmp(1),il_vtmp(2),il_vtmp(3)), stat=ierr)
             IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate vtmp)',& 
                  __LINE__,__FILE__)
             ALLOCATE(fnlat(il_fnlat(1),il_fnlat(2),il_fnlat(3)), stat=ierr)
             IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlat)',& 
                  __LINE__,__FILE__)
#endif
             END IF
             ijv=0
             DO iv=1,nlps_com%ngh(is)
                DO jv=iv,nlps_com%ngh(is)
                   ijv=ijv+1
                   CALL qvan2(iv,jv,is,qg(1:ncpw%nhg_cp,ijv))
#if defined(__VECTOR)
                   !$omp parallel do private(IG)
#else
                   !$omp parallel do private(IG) schedule(static)
#endif
                   DO ig=1,ncpw%nhg_cp
                      qg(ig,ijv)=CONJG(qg(ig,ijv))*vpot(ncpw%nhg_start-1+ig)
                   ENDDO
                ENDDO
             ENDDO
             CALL dgemm('T','N',ions0%na(is),nhh,2*ncpw%nhg_cp,tom,&
                  eigrb(ncpw%nhg_start,isa0+1),2*ncpw%nhg,&
                  qg(1,1),2*ncpw%nhg_cp,&
                  0.0_real_8,ylm1(1,1),ions0%na(is))
             IF (geq0 .AND. ncpw%nhg_start .EQ. 1 ) then
                !$omp parallel do private(iv,jv,ijv,i,j,ia,isa) schedule(static,1)
                DO iv=1,nlps_com%ngh(is)
                   DO jv=iv,nlps_com%ngh(is)
                      ijv=0
                      DO i=1,iv-1
                         DO j=i,nlps_com%ngh(is)
                            ijv=ijv+1
                         END DO
                      END DO
                      DO i=iv,jv
                         ijv=ijv+1
                      END DO
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         ylm1(ia,ijv)=ylm1(ia,ijv)-parm%omega*REAL(qg(1,ijv)*eigrb(1,isa))
                      ENDDO
                   ENDDO
                ENDDO
             END IF
             CALL mp_sum(ylm1,ions0%na(is)*nhh,parai%cp_grp)
             !$omp parallel do private(iv,jv,ijv,i,j,ia,isa) 
             DO iv=1,nlps_com%ngh(is)
                DO jv=iv,nlps_com%ngh(is)
                   ijv=0
                   DO i=1,iv-1
                      DO j=i,nlps_com%ngh(is)
                         ijv=ijv+1
                      END DO
                   END DO
                   DO i=iv,jv
                      ijv=ijv+1
                   END DO
                   DO ia=1,ions0%na(is)
                      isa=isa0+ia
                      deeq(isa,iv,jv)=ylm1(ia,ijv)
                   ENDDO
                   DO ia=1,ions0%na(is)
                      isa=isa0+ia
                      deeq(isa,jv,iv)=ylm1(ia,ijv)
                   ENDDO
                ENDDO
             ENDDO
             IF (tfor) then  
                ijv=0
                IF (cntl%tfdist) THEN
                   CALL zeroing(ylm1)!,maxsys%nax)
                   DO iv=1,nlps_com%ngh(is)
                      DO jv=iv,nlps_com%ngh(is)
                         ijv=ijv+1
                         !$omp parallel do private(IA,ISA,II,I,tmp,ftmp) 
                         DO ia=1,ions0%na(is)
                            isa=isa0+ia
                            IF (iatpe(isa).EQ.parai%mepos) THEN
                               ii=isa-ipept(1,parai%mepos)+1
                               tmp=0._real_8
                               DO i=1,nstate
                                  tmp=tmp+f(i)*fnla(ii,iv,i)*&
                                       fnla(ii,jv,i)
                               END DO
                               FTMP=OMTPIBA
                               IF(IV.NE.JV) FTMP=FTMP*2._real_8
                               ylm1(ia,ijv)=tmp*ftmp
                            ENDIF
                         ENDDO
                      END DO
                   END DO
                   CALL mp_sum(ylm1,ions0%na(is)*nhh,parai%allgrp)
                ELSE
                   methread=1
                   !$omp parallel private(ia,isa,iv,i,ii,ijv,jv,tmp,ftmp,methread) 
                   !$ methread=omp_get_thread_num()+1
                   DO ijv=1,nhh
                      !$omp DO 
                      DO ia=1,ions0%na(is)
                         ylm1(ia,ijv)=0.0_real_8
                      END DO
                   END DO
                   !$omp do 
                   DO ia=na_grp(1,is,parai%cp_inter_me +1),na_grp(2,is,parai%cp_inter_me +1)
                      isa=isa0+ia
                      !$omp simd
                      DO iv=1,nlps_com%ngh(is)
                         DO i=paraw%nwa12(1,parai%mepos),paraw%nwa12(2,parai%mepos)
                            ii=i-paraw%nwa12(1,parai%mepos)+1
                            fnlat(ii,iv,methread)=fnla(isa,iv,i)
                         END DO
                      END DO
                      IJV=0
                      DO iv=1,nlps_com%ngh(is)
                         DO jv=iv,nlps_com%ngh(is)
                            ijv=ijv+1
                            tmp=0._real_8
                            !$omp simd reduction(+:tmp)
                            DO i=paraw%nwa12(1,parai%mepos),paraw%nwa12(2,parai%mepos)
                               ii=i-paraw%nwa12(1,parai%mepos)+1
                               tmp=tmp+f(i)*fnlat(ii,iv,methread)*fnlat(ii,jv,methread)
                            ENDDO
                            FTMP=OMTPIBA
                            IF(IV.NE.JV) FTMP=FTMP*2._real_8
                            YLM1(ia,ijv)=tmp*FTMP
                         END DO
                      END DO
                   END DO
                   !$omp end parallel
                   CALL mp_sum(ylm1,nhh*ions0%na(is),parai%cp_grp)
                ENDIF
                isa=isa0+1
                !$omp parallel private(k,ig,ijv) 
                DO ijv=1,nhh
                   DO k=1,3
                      !$omp do
                      DO ig=1,ncpw%nhg_cp
                         vtmp(ig,k,ijv)=&
                              CMPLX(gktemp(ig,k)*aimag(qg(ig,ijv)),-gktemp(ig,k)*real(qg(ig,ijv),kind=real_8))
                      END DO
                      !$omp end do nowait
                   END DO
                END DO
                !$omp end parallel
                CALL dgemm('T','N',3*nhh,ions0%na(is),2*ncpw%nhg_cp,2._real_8,&
                     vtmp(1,1,1),2*ncpw%nhg_cp, &
                     eigrb(ncpw%nhg_start,isa),2*ncpw%nhg, &
                     0.0_real_8,ylm2(1,1,1),3*nhh)
                IF (geq0 .AND. ncpw%nhg_start .EQ. 1) &
                     CALL dger(3*nhh,ions0%na(is),-1.0_real_8,&
                     vtmp(1,1,1),2*ncpw%nhg_cp,&
                     eigrb(1,isa),2*ncpw%nhg &
                     ,ylm2(1,1,1),3*nhh)
                !$omp parallel do private(ia,fiont,iv,jv,ijv,i,j,otr,k)
                DO ia=1,ions0%na(is)
                   fiont=0._real_8
                   ijv=0
                   DO iv=1,nlps_com%ngh(is)
                      DO jv=iv,nlps_com%ngh(is)
                         ijv=ijv+1
                         OTR=YLM1(ia,ijv)
                         DO k=1,3
                            fiont(k)=fiont(k)+OTR*YLM2(k,ijv,ia)
                         ENDDO
                      END DO
                   END DO
                   DO k=1,3
                      fion(k,ia,is)=fion(k,ia,is)+fiont(k)
                   END DO
                END DO
             ENDIF
          END IF
          IF (tfor) THEN
#ifdef _USE_SCRATCHLIBRARY
             CALL free_scratch(il_fnlat,fnlat,procedureN//'_fnlat')
             CALL free_scratch(il_vtmp,vtmp,procedureN//'_vtmp')
             CALL free_scratch(il_ylm2,ylm2,procedureN//'_ylm2')
#else
             DEALLOCATE(ylm2, stat=ierr)
             IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate ylm2)',& 
                  __LINE__,__FILE__)
             DEALLOCATE(vtmp, stat=ierr)
             IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate vtmp)',& 
                  __LINE__,__FILE__)
             DEALLOCATE(fnlat, stat=ierr)
             IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fnlat)',& 
                  __LINE__,__FILE__)
#endif
          END IF
#ifdef _USE_SCRATCHLIBRARY
          CALL free_scratch(il_ylm1,ylm1,procedureN//'_ylm1')
          CALL free_scratch(il_qg,qg,procedureN//'_qg')
#else
          DEALLOCATE(qg, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate qg)',& 
               __LINE__,__FILE__)
          DEALLOCATE(ylm1, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate ylm1)',& 
               __LINE__,__FILE__)
#endif
          isa0=isa0+ions0%na(is)
       END DO
       IF (tfor) THEN
#ifdef _USE_SCRATCHLIBRARY
          CALL free_scratch(il_gktemp,gktemp,procedureN//'_gktemp')
#else
          DEALLOCATE(gktemp, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate ylm1)',& 
               __LINE__,__FILE__)
#endif
       END IF
       IF (parai%cp_nogrp.gt.1 ) then
          CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%cp_inter_grp)
       END IF
    ELSE
       ! ==--------------------------------------------------------------==
       il_qg(1)=ncpw%nhg
       il_qg(2)=1
       il_vtmp(1)=ncpw%nhg
       il_vtmp(2)=1
       il_vtmp(3)=1
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_qg,qg,procedureN//'_qg')
       CALL request_scratch(il_vtmp,vtmp,procedureN//'_vtmp')
#else
       ALLOCATE(qg(il_qg(1),il_qg(2)), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate qg)',& 
            __LINE__,__FILE__)
       ALLOCATE(vtmp(il_vtmp(1),il_vtmp(2),il_vtmp(3)), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ylm1)',& 
            __LINE__,__FILE__)
#endif
       isa0=0
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is)) THEN
             ijv=0
             DO iv=1,nlps_com%ngh(is)
                DO jv=iv,nlps_com%ngh(is)
                   ijv=ijv+1
                   CALL qvan2(iv,jv,is,qg(1:ncpw%nhg,1))
                   !$omp parallel do private(IG) shared(QG) schedule(static)
                   DO ig=1,ncpw%nhg
                      qg(ig,1)=CONJG(qg(ig,1))*vpot(ig)
                   ENDDO
                   DO ia=1,ions0%na(is)
                      isa=isa0+ia
                      !$omp parallel do private(IG) shared(VTMP) schedule(static)
                      DO ig=1,ncpw%nhg
                         vtmp(ig,1,1)=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                              ei3(isa,inyh(3,ig))
                      ENDDO
                      deeq(isa,iv,jv)=parm%omega*dotp(ncpw%nhg,qg(1,1),vtmp(1,1,1))
                      deeq(isa,jv,iv)=deeq(isa,iv,jv)
                      IF (tfor) THEN
                         rhovan=0.0_real_8
                         IF (cntl%tfdist) THEN
                            IF (iatpe(isa).EQ.parai%mepos) THEN
                               ii=isa-ipept(1,parai%mepos)+1
                               !$omp parallel do private(I) reduction(+:RHOVAN)
                               DO i=1,nstate
                                  rhovan=rhovan+f(i)*fnla(ii,iv,i)*fnla(ii,jv,i)
                               ENDDO
                            ENDIF
                            CALL mp_sum(rhovan,parai%allgrp)
                         ELSE
                            !$omp parallel do private(I) reduction(+:RHOVAN)
                            DO i=1,nstate
                               rhovan=rhovan+f(i)*fnla(isa,iv,i)*fnla(isa,jv,i)
                            ENDDO
                         ENDIF
                         otr=parm%omega*parm%tpiba*rhovan
                         DO k=1,3
                            CALL cftemp(ncpw%nhg,qg(1,1),vtmp(1,1,1),gk,k,ftmp)
                            IF (iv.NE.jv) ftmp=2.0_real_8*ftmp
                            fion(k,ia,is)=fion(k,ia,is)+ftmp*otr
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          isa0=isa0+ions0%na(is)
       ENDDO
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_vtmp,vtmp,procedureN//'_vtmp')
       CALL free_scratch(il_qg,qg,procedureN//'_qg')
#else
       DEALLOCATE(qg, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate qg)',& 
            __LINE__,__FILE__)
       DEALLOCATE(vtmp, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ylm1)',& 
            __LINE__,__FILE__)
#endif
       CALL mp_sum(deeq,ions1%nat*maxsys%nhxs*maxsys%nhxs,parai%allgrp)
    ENDIF

#ifdef _VERBOSE_FORCE_DBG
    ALLOCATE(dbg_forces(3,maxsys%nax,maxsys%nsx), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate dbg_forces',& 
         __LINE__,__FILE__)
    dbg_forces=fion
    CALL mp_sum(dbg_forces,3*maxsys%nax*maxsys%nsx,parai%allgrp)
    IF (paral%io_parent) THEN
       WRITE(6,*) "===================================="
       WRITE(6,*) "DEBUG FORCES", procedureN
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             WRITE(6,*) dbg_forces(1:3,ia,is),ia,is
          END DO
       END DO
    END IF
    DEALLOCATE(dbg_forces,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#endif
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE newd
  ! ==================================================================
  SUBROUTINE give_scr_newd(lnewd,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lnewd
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: nhh

    IF (cntl%bigmem) THEN
       nhh = (maxsys%nhxs*(maxsys%nhxs+1))/2
       lnewd = 2*ncpw%nhg*nhh + 6*ncpw%nhg + maxsys%nax*(4+nhh) + 100
    ELSE
       lnewd = 2*ncpw%nhg + 6*ncpw%nhg + maxsys%nax + 100
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_newd
  ! ==================================================================


END MODULE newd_utils

! ==================================================================
SUBROUTINE cftemp(nhg,qg,vtmp,gk,k,ftmp)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: nhg
  REAL(real_8)                               :: qg(2,nhg), vtmp(2,nhg), &
                                                gk(3,nhg)
  INTEGER                                    :: k
  REAL(real_8)                               :: ftmp

  INTEGER                                    :: ig

! Variables
! ==--------------------------------------------------------------==
! G=0 TERM IS ZERO !

  ftmp=0.0_real_8
#if defined(__VECTOR)
  !$omp parallel do private(IG) reduction(+:FTMP)
#else
  !$omp parallel do private(IG) reduction(+:FTMP) schedule(static)
#endif
  DO ig=1,nhg
     ftmp=ftmp+gk(k,ig)*(vtmp(1,ig)*qg(2,ig)-vtmp(2,ig)*qg(1,ig))
  ENDDO
  ftmp=2.0_real_8*ftmp
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE cftemp
