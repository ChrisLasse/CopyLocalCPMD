MODULE rnlrh_utils
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms
  USE cvan,                            ONLY: dvan
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_c,&
                                             ener_com
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: wk
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlm,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: fnl,&
                                             fnl2,&
                                             fnlgp,&
                                             fnla
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE spin,                            ONLY: clsd,&
                                             lspin2
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             ipept,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlrh
  !public :: rnlcas
  PUBLIC :: rnlrhg

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlrh(enl,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==     THE NON-LOCAL POTENTIAL CONTRIBUTION TO THE TOTAL        ==
    ! ==     ENERGY, I.E. ENL                                         ==
    ! ==     K-POINTS IMPLEMENTED (FNL IS COMPLEX -> IMAGP=2)         ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: enl
    INTEGER                                  :: nstate, nkpoint

    INTEGER                                  :: i, ia, ii, ik, is, isa, isa0, &
                                                isub, iv, jv, ki, kj, l, l2, &
                                                li, lj, isa_start, isa_end, &
                                                na_grp(2,ions1%nsp,parai%cp_nogrp)
    REAL(real_8)                             :: sum, weight, enl_priv

    ! split atoms between cp groups
    call cp_grp_split_atoms(na_grp)

    enl=0.0_real_8

    ! If no non-local components -> return.
    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('     RNLRH',isub)
    ! ==--------------------------------------------------------------==
    ! == Compute the non-local contribution to the total energy (ENL) ==
    ! ==--------------------------------------------------------------==
    IF (pslo_com%tivan) THEN
       ! VANDERBILT PSEUDOPOTENTIAL
       IF (imagp.EQ.2)&
            CALL stopgm('RNLRH','K-POINTS NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       !$omp parallel private (ik,isa0,is,i,ii,iv,jv,sum,ia,isa,isa_end,isa_start) reduction(+:enl) 
       DO ik=1,nkpoint
          isa0=0
          DO is=1,ions1%nsp
             IF (pslo_com%tvan(is)) THEN
                isa_start=isa0+na_grp(1,is,parai%cp_inter_me +1)
                isa_end=isa_start+na_grp(2,is,parai%cp_inter_me +1)-na_grp(1,is,parai%cp_inter_me +1)
                IF(isa_end-isa_start+1.GT.0) THEN
                   !$omp do 
                   DO i=parap%nst12(1,parai%mepos),parap%nst12(2,parai%mepos)
                      ii=i-parap%nst12(1,parai%mepos)+1
                      DO iv=1,nlps_com%ngh(is)
                         DO jv=1,nlps_com%ngh(is)
                            sum=0.0_real_8
                            IF (cntl%tfdist) THEN
                               DO ia=1,ions0%na(is)
                                  isa=isa0+ia
                                  sum=sum+crge%f(i,ik)*fnl2(1,isa,iv,ii,ik)*&
                                       fnl2(1,isa,jv,ii,ik)
                               ENDDO
                            ELSE
                               !$omp simd reduction(+:sum)
                               DO isa=isa_start,isa_end
                                  sum=sum+fnla(isa,iv,i)*fnla(isa,jv,i)
                               ENDDO
                            ENDIF
                            enl=enl+dvan(jv,iv,is)*sum*crge%f(i,ik)
                         ENDDO
                      ENDDO
                   ENDDO
                   !$omp end do nowait
                END IF
             END IF
             isa0=isa0+ions0%na(is)
          END DO
       END DO
       !$omp end parallel
    ELSE
       DO ik=1,nkpoint
          isa0=0
          DO is=1,ions1%nsp
             IF (sgpp1%tsgp(is)) THEN
                ! Stefan Goedecker pp
                DO iv=1,nlps_com%ngh(is)
                   l=nghtol(iv,is)+1
                   ki=sgpp2%lfval(iv,is)
                   li=sgpp2%lpval(iv,is)
                   DO jv=iv,nlps_com%ngh(is)
                      l2=nghtol(jv,is)+1
                      lj=sgpp2%lpval(jv,is)
                      IF (l2.EQ.l.AND.li.EQ.lj) THEN
                         kj=sgpp2%lfval(jv,is)
                         sum=0.0_real_8
                         DO i=parap%nst12(1,parai%mepos),parap%nst12(2,parai%mepos)
                            weight=wk(ik)*crge%f(i,ik)
                            IF (weight.EQ.0._real_8) GOTO 2000
                            ii=i-parap%nst12(1,parai%mepos)+1
                            IF (imagp.EQ.1) THEN
                               IF (cntl%tfdist) THEN
                                  DO ia=1,ions0%na(is)
                                     isa=isa0+ia
                                     sum=sum+weight*fnl2(1,isa,iv,ii,ik)*&
                                          fnl2(1,isa,jv,ii,ik)
                                  ENDDO
                               ELSE
                                  DO ia=na_grp(1,is,parai%cp_inter_me +1),na_grp(2,is,parai%cp_inter_me +1)
                                     isa=isa0+ia
                                     sum=sum+weight*fnl2(1,isa,iv,i,ik)*&
                                          fnl2(1,isa,jv,i,ik)
                                  ENDDO
                               ENDIF
                            ELSE
                               IF (cntl%tfdist) THEN
                                  DO ia=1,ions0%na(is)
                                     isa=isa0+ia
                                     sum=sum+weight*&
                                          (fnl2(1,isa,iv,ii,ik)*fnl2(1,isa,jv,ii,ik)&
                                          +fnl2(2,isa,iv,ii,ik)*fnl2(2,isa,jv,ii,ik))
                                  ENDDO
                               ELSE
                                  DO ia=na_grp(1,is,parai%cp_inter_me +1),na_grp(2,is,parai%cp_inter_me +1)
                                     isa=isa0+ia
                                     sum=sum+weight*&
                                          (fnl2(1,isa,iv,i,ik)*fnl2(1,isa,jv,i,ik)&
                                          +fnl2(2,isa,iv,i,ik)*fnl2(2,isa,jv,i,ik))
                                  ENDDO
                               ENDIF
                            ENDIF
2000                        CONTINUE
                         ENDDO
                         IF (iv.NE.jv) sum=2._real_8*sum
                         enl=enl+sum*sgpp2%hlsg(ki,kj,l,is)
                      ENDIF
                   ENDDO
                ENDDO
             ELSE
                ! BHS AND RELATED 
                DO iv=1,nlps_com%ngh(is)
                   sum=0.0_real_8
                   DO i=parap%nst12(1,parai%mepos),parap%nst12(2,parai%mepos)
                      ii=i-parap%nst12(1,parai%mepos)+1
                      IF (imagp.EQ.1) THEN
                         IF (cntl%tfdist) THEN
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               sum=sum+crge%f(i,ik)*wk(ik)*&
                                    fnl2(1,isa,iv,ii,ik)*fnl2(1,isa,iv,ii,ik)
                            ENDDO
                         ELSE
                            DO ia=na_grp(1,is,parai%cp_inter_me +1),na_grp(2,is,parai%cp_inter_me +1)
                               isa=isa0+ia
                               sum=sum+crge%f(i,ik)*wk(ik)*&
                                    fnl2(1,isa,iv,i,ik)*fnl2(1,isa,iv,i,ik)
                            ENDDO
                         ENDIF
                      ELSE
                         IF (cntl%tfdist) THEN
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               sum=sum+crge%f(i,ik)*wk(ik)*&
                                    (fnl2(1,isa,iv,ii,ik)*fnl2(1,isa,iv,ii,ik)&
                                    +fnl2(2,isa,iv,ii,ik)*fnl2(2,isa,iv,ii,ik))
                            ENDDO
                         ELSE
                            DO ia=na_grp(1,is,parai%cp_inter_me +1),na_grp(2,is,parai%cp_inter_me +1)
                               isa=isa0+ia
                               sum=sum+crge%f(i,ik)*wk(ik)*&
                                    (fnl2(1,isa,iv,i,ik)*fnl2(1,isa,iv,i,ik)&
                                    +fnl2(2,isa,iv,i,ik)*fnl2(2,isa,iv,i,ik))
                            ENDDO
                         ENDIF
                      ENDIF
                   ENDDO
                   enl=enl+wsg(is,iv)*sum
                ENDDO
             ENDIF
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO
    END IF
    if (parai%cp_nogrp .gt. 1) call mp_sum(enl,parai%cp_inter_grp)
    IF (lspin2%tlse .AND. lspin2%tcas22) CALL rnlcas
    CALL tihalt('     RNLRH',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnlrh
  ! ==================================================================
  SUBROUTINE rnlcas
    ! ==--------------------------------------------------------------==
    ! == Calculates NL PP contribution to CAS22 Energies              ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ia, iat, iatl, is, iv, jv, &
                                                ki, kj, l, l2, li, lj
    REAL(real_8)                             :: sm

! ==--------------------------------------------------------------==

    ener_c%enl_ab = 0._real_8
    ener_c%enl_a  = ener_com%enl
    DO iat=ipept(1,parai%mepos),ipept(2,parai%mepos)
       ia=iatpt(1,iat)
       is=iatpt(2,iat)
       IF (pslo_com%tvan(is)) THEN
          ! VANDERBILT PSEUDOPOTENTIAL
       ELSEIF (sgpp1%tsgp(is)) THEN
          ! Stefan Goedecker pp
          DO iv=1,nlps_com%ngh(is)
             l=nghtol(iv,is)+1
             li=sgpp2%lpval(iv,is)
             ki=sgpp2%lfval(iv,is)
             DO jv=iv,nlps_com%ngh(is)
                IF (iv.EQ.jv) THEN
                   sm=1._real_8
                ELSE
                   sm=2._real_8
                ENDIF
                l2=nghtol(jv,is)+1
                lj=sgpp2%lpval(jv,is)
                IF (l2.EQ.l.AND.li.EQ.lj) THEN
                   kj=sgpp2%lfval(jv,is)
                   IF (cntl%tfdist) THEN
                      iatl=iat-ipept(1,parai%mepos)+1
                   ELSE
                      iatl=iat
                   ENDIF
                   ener_c%enl_ab=ener_c%enl_ab+sgpp2%hlsg(ki,kj,l,is)*0.5_real_8*sm*&
                        (fnl(1,iatl,iv,clsd%ialpha,1)*fnl(1,iatl,jv,clsd%ibeta,1)+&
                        fnl(1,iatl,jv,clsd%ialpha,1)*fnl(1,iatl,iv,clsd%ibeta,1))
                   ener_c%enl_a=ener_c%enl_a+sgpp2%hlsg(ki,kj,l,is)*sm*&
                        (fnl(1,iatl,iv,clsd%ialpha,1)*fnl(1,iatl,jv,clsd%ialpha,1)-&
                        fnl(1,iatl,iv,clsd%ibeta,1)*fnl(1,iatl,jv,clsd%ibeta,1))
                ENDIF
             ENDDO
          ENDDO
       ELSE
          ! BHS AND RELATED
          DO iv=1,nlps_com%ngh(is)
             IF (cntl%tfdist) THEN
                iatl=iat-ipept(1,parai%mepos)+1
             ELSE
                iatl=iat
             ENDIF
             ener_c%enl_ab=ener_c%enl_ab+wsg(is,iv)*fnl(1,iatl,iv,clsd%ialpha,1)*&
                  fnl(1,iatl,iv,clsd%ibeta,1)
             ener_c%enl_a=ener_c%enl_a+wsg(is,iv)*(fnl(1,iatl,iv,clsd%ialpha,1)**2&
                  -fnl(1,iatl,iv,clsd%ibeta,1)**2)
          ENDDO
       ENDIF
    ENDDO
    ener_c%enl_2  = ener_com%enl - (ener_c%enl_a - ener_com%enl)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnlcas
  ! EHR[
  ! ==================================================================
  SUBROUTINE rnlrhg(enl,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==     THE NON-LOCAL POTENTIAL CONTRIBUTION TO THE TOTAL        ==
    ! ==     ENERGY, I.E. ENL                                         ==
    ! ==     K-POINTS IMPLEMENTED (FNL IS COMPLEX -> IMAGP=2)         ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: enl
    INTEGER                                  :: nstate, nkpoint

    INTEGER                                  :: i, ia, ii, ik, is, isa, isa0, &
                                                isub, iv
    REAL(real_8)                             :: sum

    enl=0._real_8
    ! If no non-local components -> return.
    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('     RNLRH',isub)
    ! ==--------------------------------------------------------------==
    ! == Compute the non-local contribution to the total energy (ENL) ==
    ! ==--------------------------------------------------------------==
    DO ik=1,nkpoint
       isa0=0
       DO is=1,ions1%nsp
          ! BHS AND RELATED
          DO iv=1,nlps_com%ngh(is)
             sum=0.0_real_8
             DO i=parap%nst12(1,parai%mepos),parap%nst12(2,parai%mepos)
                ii=i-parap%nst12(1,parai%mepos)+1
                IF (imagp.EQ.1) THEN
                   IF (cntl%tfdist) THEN
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         sum=sum+crge%f(i,ik)*wk(ik)*&
                              fnl2(1,isa,iv,ii,ik)*fnl2(1,isa,iv,ii,ik)
                      ENDDO
                   ELSE
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         sum=sum+crge%f(i,ik)*wk(ik)*&
                              fnl2(1,isa,iv,i,ik)*fnl2(1,isa,iv,i,ik)
                      ENDDO
                   ENDIF
                ELSE
                   IF (cntl%tfdist) THEN
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         sum=sum+crge%f(i,ik)*wk(ik)*&
                              (fnl2(1,isa,iv,ii,ik)*fnl2(1,isa,iv,ii,ik)&
                              +fnl2(2,isa,iv,ii,ik)*fnl2(2,isa,iv,ii,ik))
                      ENDDO
                   ELSE
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         sum=sum+crge%f(i,ik)*wk(ik)* 2._real_8 *&
                              ( (fnl(1,isa,iv,i,ik)*fnlgp(1,isa,iv,i,ik)&
                              +fnl(2,isa,iv,i,ik)*fnlgp(2,isa,iv,i,ik))&
                              )
                      ENDDO
                   ENDIF
                ENDIF
             ENDDO
             enl=enl+wsg(is,iv)*sum
          ENDDO
          isa0=isa0+ions0%na(is)
       ENDDO
    ENDDO
    IF (lspin2%tlse .AND. lspin2%tcas22) CALL rnlcas
    CALL tihalt('     RNLRH',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnlrhg
  ! ==================================================================
  ! EHR]

END MODULE rnlrh_utils
