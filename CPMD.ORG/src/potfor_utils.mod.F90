#include "cpmd_global.h"

#if defined(__SR11000)
!option OPT(O(ss))
#endif

MODULE potfor_utils
  USE cppt,                            ONLY: gk,&
                                             inyh,&
                                             nzh,&
                                             rhops,&
                                             scg,&
                                             vps
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             iatpt,&
                                             maxsys,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: potfor
  PUBLIC :: potabfor

CONTAINS

  ! ==================================================================
  SUBROUTINE potfor(fion,v,eirop)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == POTENTIAL ENERGY CONTRIBUTIONS TO THE FORCES ON THE IONS     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)
    COMPLEX(real_8)                          :: v(*), eirop(ncpw%nhg)

    COMPLEX(real_8)                          :: ei123, gx, gy, gz, rhet, &
                                                rhets, rhog, rhogs, rp, txx, &
                                                tyy, tzz, vcgs, fiont(3)
    INTEGER                                  :: ia, ig, ig1, is, isa, isub, ierr
    REAL(real_8)                             :: omtp
#ifdef _VERBOSE_FORCE_DBG
    REAL(real_8),ALLOCATABLE                 :: dbg_forces(:,:,:)
#endif
    CHARACTER(*), PARAMETER                  :: procedureN = 'potfor'

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
    ig1=ncpw%nhg_start
    IF (geq0 .and. ncpw%nhg_start .eq. 1) ig1=2
!#if defined (__VECTOR)
    IF (cntl%bigmem) THEN
       !$omp parallel do private(ISA,IA,IS,IG,EI123,RP,RHET,RHOG,RHETS,RHOGS, &
       !$omp  GX,GY,GZ,VCGS,fiont)
!#ifdef _vpp_
!       !OCL NOALIAS
!#endif
       DO isa=1,ions1%nat
          ia=iatpt(1,isa)
          is=iatpt(2,isa)
          fiont=0._real_8
          DO ig=ig1,ncpw%nhg_last
             ei123=eigrb(ig,isa)
             rp=eirop(ig)
             rhet=v(nzh(ig))
             rhog=rhet+rp
             rhets=CONJG(rhet)
             rhogs=CONJG(rhog)
             gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
             gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
             gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
             vcgs=scg(ig)*rhogs
             ei123=ei123*(rhops(is,ig)*vcgs+vps(is,ig)*rhets)
             fiont(1)=fiont(1)+REAL(ei123*gx)
             fiont(2)=fiont(2)+REAL(ei123*gy)
             fiont(3)=fiont(3)+REAL(ei123*gz)

!             fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*gx)*omtp
!             fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*gy)*omtp
!             fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*gz)*omtp
          ENDDO          
             fion(1,ia,is)=fion(1,ia,is)+fiont(1)*omtp
             fion(2,ia,is)=fion(2,ia,is)+fiont(2)*omtp
             fion(3,ia,is)=fion(3,ia,is)+fiont(3)*omtp

       ENDDO
    ELSE
       !$omp parallel do private(ISA,IA,IS,IG,EI123,RP,RHET,RHOG,RHETS,RHOGS, &
       !$omp  GX,GY,GZ,VCGS,fiont)
!#ifdef _vpp_
!       !OCL NOALIAS
!#endif
       DO isa=1,ions1%nat
          ia=iatpt(1,isa)
          is=iatpt(2,isa)
          fiont=0._real_8
          DO ig=ig1,ncpw%nhg_last
             ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                  ei3(isa,inyh(3,ig))
             rp=eirop(ig)
             rhet=v(nzh(ig))
             rhog=rhet+rp
             rhets=CONJG(rhet)
             rhogs=CONJG(rhog)
             gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
             gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
             gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
             vcgs=scg(ig)*rhogs
             ei123=ei123*(rhops(is,ig)*vcgs+vps(is,ig)*rhets)
!             fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*gx)*omtp
!             fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*gy)*omtp
!             fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*gz)*omtp
             fiont(1)=fiont(1)+REAL(ei123*gx)
             fiont(2)=fiont(2)+REAL(ei123*gy)
             fiont(3)=fiont(3)+REAL(ei123*gz)
          ENDDO
          fion(1,ia,is)=fion(1,ia,is)+fiont(1)*omtp
          fion(2,ia,is)=fion(2,ia,is)+fiont(2)*omtp
          fion(3,ia,is)=fion(3,ia,is)+fiont(3)*omtp
       ENDDO
    ENDIF
!#else                
!    ! mb   For scalar version & pseudo-vector machines
!    IF (cntl%bigmem) THEN
!#ifdef __SR8000
!       !poption parallel
!       !poption tlocal(RP,RHET,RHOG,RHETS,RHOGS,GX,GY,GZ,VCGS)
!       !poption tlocal(ISA,TXX,TYY,TZZ,EI123)
!       !poption psum(FION)
!#else
!#ifdef __HPC
!       ! no-OMP - segmentation fault if OMP are used here
!#else
!       !$omp parallel do private(IG,IS,IA,ISA) &
!       !$omp  private(RP,RHET,RHOG,RHETS,RHOGS,GX,GY,GZ,VCGS) &
!       !$omp  private(TXX,TYY,TZZ,EI123) &
!       !$omp  reduction(+:FION)
!#endif
!#endif
!       DO ig=ig1,ncpw%nhg
!          rp=eirop(ig)
!          rhet=v(nzh(ig))
!          rhog=rhet+rp
!          rhets=CONJG(rhet)
!          rhogs=CONJG(rhog)
!          gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
!          gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
!          gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
!          vcgs=scg(ig)*rhogs
!          isa=0
!          DO is=1,ions1%nsp
!             txx=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)*gx
!             tyy=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)*gy
!             tzz=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)*gz
!             DO ia=1,ions0%na(is)
!                isa=isa+1
!                ei123=eigrb(ig,isa)
!                fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
!                fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
!                fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
!             ENDDO
!          ENDDO
!       ENDDO
!    ELSE
!#ifdef __SR8000
!       !poption parallel
!       !poption tlocal(RP,RHET,RHOG,RHETS,RHOGS,GX,GY,GZ,VCGS)
!       !poption tlocal(ISA,TXX,TYY,TZZ,EI123)
!       !poption psum(FION)
!#else
!#ifdef __HPC
!       ! no-OMP - segmentation fault if OMP are used here
!#else
!       !$omp parallel do private(IG,IS,IA,ISA) &
!       !$omp  private(RP,RHET,RHOG,RHETS,RHOGS,GX,GY,GZ,VCGS) &
!       !$omp  private(TXX,TYY,TZZ,EI123) &
!       !$omp  reduction(+:FION)
!#endif
!#endif
!       DO ig=ig1,ncpw%nhg
!          rp=eirop(ig)
!          rhet=v(nzh(ig))
!          rhog=rhet+rp
!          rhets=CONJG(rhet)
!          rhogs=CONJG(rhog)
!          gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
!          gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
!          gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
!          vcgs=scg(ig)*rhogs
!          isa=0
!          DO is=1,ions1%nsp
!             txx=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)*gx
!             tyy=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)*gy
!             tzz=(rhops(is,ig)*vcgs+vps(is,ig)*rhets)*gz
!             DO ia=1,ions0%na(is)
!                isa=isa+1
!                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
!                     ei3(isa,inyh(3,ig))
!                fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
!                fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
!                fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
!             ENDDO
!          ENDDO
!       ENDDO
!    ENDIF
!#endif

    if (parai%cp_nogrp.gt.1) then
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
  END SUBROUTINE potfor
  ! ==================================================================
  SUBROUTINE potabfor(fion,rhog)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == POTENTIAL ENERGY CONTRIBUTIONS TO THE FORCES ON THE IONS     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)
    COMPLEX(real_8)                          :: rhog(fpar%nnr1)

    CHARACTER(*), PARAMETER                  :: procedureN = 'potabfor'

    COMPLEX(real_8)                          :: ei123, gx, gy, gz, rg, txx, &
                                                tyy, tzz
    INTEGER                                  :: ia, ig, ig1, is, isa, isub
    REAL(real_8)                             :: omtp

    CALL tiset(procedureN,isub)
    omtp=2._real_8*parm%omega*parm%tpiba
    ig1=1
    IF (geq0) ig1=2
#if defined (__VECTOR) && (! (__NEC))
    isa=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          isa=isa+1
          IF (cntl%bigmem) THEN
             !$omp parallel do private(IG,EI123,RG,GX,GY,GZ,TXX,TYY,TZZ) &
             !$omp  shared(ISA,RHOG)
             DO ig=ig1,ncpw%nhg
                ei123=eigrb(ig,isa)
                rg=CONJG(rhog(ig))
                ! mb              GX=CMPLX(0._real_8,GK(1,IG))
                ! mb              GY=CMPLX(0._real_8,GK(2,IG))
                ! mb              GZ=CMPLX(0._real_8,GK(3,IG))
                txx=vps(is,ig)*rg*CMPLX(0._real_8,gk(1,ig),kind=real_8)
                tyy=vps(is,ig)*rg*CMPLX(0._real_8,gk(2,ig),kind=real_8)
                tzz=vps(is,ig)*rg*CMPLX(0._real_8,gk(3,ig),kind=real_8)
                fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
                fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
                fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
             ENDDO
          ELSE
             !$omp parallel do private(IG,EI123,RG,GX,GY,GZ,TXX,TYY,TZZ) &
             !$omp  shared(ISA,RHOG)
             DO ig=ig1,ncpw%nhg
                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
                rg=CONJG(rhog(ig))
                ! mb              GX=CMPLX(0._real_8,GK(1,IG))
                ! mb              GY=CMPLX(0._real_8,GK(2,IG))
                ! mb              GZ=CMPLX(0._real_8,GK(3,IG))
                txx=vps(is,ig)*rg*CMPLX(0._real_8,gk(1,ig),kind=real_8)
                tyy=vps(is,ig)*rg*CMPLX(0._real_8,gk(2,ig),kind=real_8)
                tzz=vps(is,ig)*rg*CMPLX(0._real_8,gk(3,ig),kind=real_8)
                fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
                fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
                fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
             ENDDO
          ENDIF
       ENDDO
    ENDDO
#else
    DO ig=ig1,ncpw%nhg
       rg=CONJG(rhog(ig))
       gx=CMPLX(0._real_8,gk(1,ig),kind=real_8)
       gy=CMPLX(0._real_8,gk(2,ig),kind=real_8)
       gz=CMPLX(0._real_8,gk(3,ig),kind=real_8)
       isa=0
       DO is=1,ions1%nsp
          txx=vps(is,ig)*rg*gx
          tyy=vps(is,ig)*rg*gy
          tzz=vps(is,ig)*rg*gz
          DO ia=1,ions0%na(is)
             isa=isa+1
             IF (cntl%bigmem) THEN
                ei123=eigrb(ig,isa)
             ELSE
                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
             ENDIF
             fion(1,ia,is)=fion(1,ia,is)+REAL(ei123*txx)*omtp
             fion(2,ia,is)=fion(2,ia,is)+REAL(ei123*tyy)*omtp
             fion(3,ia,is)=fion(3,ia,is)+REAL(ei123*tzz)*omtp
          ENDDO
       ENDDO
    ENDDO
#endif
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE potabfor
  ! ==================================================================

END MODULE potfor_utils
