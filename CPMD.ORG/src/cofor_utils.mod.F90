#include "cpmd_global.h"

MODULE cofor_utils
  USE cppt,                            ONLY: gk,&
                                             inyh
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions1,&
                                             ions0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlcc,                            ONLY: corel,&
                                             rhoc,&
                                             vnlcc,&
                                             vnlt
  USE parac,                           ONLY: parai, &
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: ei1t,&
                                             ei2t,&
                                             ei3t,&
                                             eigrb
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             maxsys,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cofor

CONTAINS

  ! ==================================================================
  SUBROUTINE cofor(fion,vpot)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == FORCES ON IONS DUE TO CORE CHARGES                           ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)
    COMPLEX(real_8)                          :: vpot(ncpw%nhg,*)

    COMPLEX(real_8)                          :: ei123, vxc
    INTEGER                                  :: ia, ig, is, isa, isub, k, isa0, ierr
    REAL(real_8)                             :: omtp, vcgs, fiont(3)
#ifdef _VERBOSE_FORCE_DBG
    REAL(real_8),ALLOCATABLE                 :: dbg_forces(:,:,:)
#endif
    CHARACTER(*), PARAMETER                  :: procedureN = 'cofor'

    CALL tiset(procedureN,isub)

    if (parai%cp_nogrp.gt.1 .and. parai%cp_inter_me .gt. 0) then
       !$omp parallel do private(is,ia)
       do is=1,ions1%nsp
          do ia=1,ions0%na(is)
             fion(1:3,ia,is)=0._real_8
          end do
       end do
    end if

    omtp=2._real_8*parm%omega*parm%tpiba
    IF (cntl%bigmem) THEN
       if(cntl%tlsd) then
          isa0=0
          DO is=1,ions1%nsp
             IF (corel%tnlcc(is)) THEN
                !$omp parallel do private(ia,isa,ig,vxc,vcgs,k,fiont)
                do ia=1,ions0%na(is)
                   fiont=0._real_8
                   isa=isa0+ia
                   DO ig=ncpw%nhg_start,ncpw%nhg_last
                      vxc=0.5_real_8*(vpot(ig,1)+vpot(ig,2)-2*vnlt(ig)-vnlcc(ig,1)-vnlcc(ig,2))
                      vcgs=-AIMAG(CONJG(vxc)*eigrb(ig,isa)*rhoc(ig,is))                  
                      do k=1,3
                         fiont(k)=fiont(k)+gk(k,ig)*vcgs*omtp
                      end do
                   ENDDO
                   do k=1,3
                      fion(k,ia,is)=fion(k,ia,is)+fiont(k)
                   end do
                END do
             end if
             isa0=isa0+ions0%na(is)
          ENDDO
       ELSE
          isa0=0
          DO is=1,ions1%nsp
             IF (corel%tnlcc(is)) THEN
                !$omp parallel do private(ia,isa,ig,vxc,vcgs,k,fiont)
                do ia=1,ions0%na(is)
                   fiont=0._real_8
                   isa=isa0+ia
                   DO ig=ncpw%nhg_start,ncpw%nhg_last
                      vxc=vpot(ig,1)-vnlt(ig)-vnlcc(ig,1)
                      vcgs=-AIMAG(CONJG(vxc)*eigrb(ig,isa)*rhoc(ig,is))                  
                      do k=1,3
                         fiont(k)=fiont(k)+gk(k,ig)*vcgs*omtp
                      end do
                   ENDDO
                   do k=1,3
                      fion(k,ia,is)=fion(k,ia,is)+fiont(k)
                   end do
                END do
             end if
             isa0=isa0+ions0%na(is)
          ENDDO
       END IF
    ELSE
       !$omp parallel do private(ISA,IA,IS,IG,EI123,VXC,VCGS) shared(OMTP)
       DO isa=1,ions1%nat
          ia=iatpt(1,isa)
          is=iatpt(2,isa)
          IF (corel%tnlcc(is)) THEN
             DO ig=ncpw%nhg_start,ncpw%nhg_last
                ei123=ei1t(inyh(1,ig),isa)*ei2t(inyh(2,ig),isa)*&
                     ei3t(inyh(3,ig),isa)
                vxc=vpot(ig,1)-vnlt(ig)-vnlcc(ig,1)
                IF (cntl%tlsd) vxc=0.5_real_8*(vxc+vpot(ig,2)-vnlt(ig)-vnlcc(ig,2))
                vcgs=-AIMAG(CONJG(vxc)*ei123*rhoc(ig,is))
                fion(1,ia,is)=fion(1,ia,is)+gk(1,ig)*vcgs*omtp
                fion(2,ia,is)=fion(2,ia,is)+gk(2,ig)*vcgs*omtp
                fion(3,ia,is)=fion(3,ia,is)+gk(3,ig)*vcgs*omtp
             ENDDO
          ENDIF
       ENDDO
    END IF

    if (parai%cp_nogrp.gt.1 ) then
       CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%cp_inter_grp)
    end if

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
  END SUBROUTINE cofor
  ! ==================================================================

END MODULE cofor_utils
