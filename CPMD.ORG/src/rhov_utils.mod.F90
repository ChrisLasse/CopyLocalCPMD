#include "cpmd_global.h"

MODULE rhov_utils
  USE cppt,                            ONLY: indz,&
                                             inyh,&
                                             nzh
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms, &
                                             cp_grp_redist_array
  USE distribution_utils,              ONLY: dist_size
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE fft_maxfft,                      ONLY: maxfftn
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_sum  
  USE nlps,                            ONLY: imagp,&
                                             nlps_com
  !$ USE omp_lib,                       ONLY: omp_get_thread_num
  USE parac,                           ONLY: parai,paral
  USE part_1d,                         ONLY: part_1d_get_blk_bounds
  USE pslo,                            ONLY: pslo_com
  USE qvan2_utils,                     ONLY: qvan2
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb,&
                                             fnl,&
                                             fnla
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw, &
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

  PUBLIC :: rhov
  PUBLIC :: give_scr_rhov

CONTAINS

  ! ==================================================================
  SUBROUTINE rhov(nstate,is1,is2,rsumv,psi)
    ! ==--------------------------------------------------------------==
    ! ==         THIS ROUTINE CALCULATES THE VANDERBILT DENSITY       ==
    ! ==                                                              ==
    ! == N_V(G) = SUM_I,IJ RHO_I,IJ Q_I,JI(G) E^-IG.R_I               ==
    ! == RHO_I,IJ = SUM_N < BETA_I,I | PSI_N >< PSI_N | BETA_I,J >    ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate, is1, is2
    REAL(real_8)                             :: rsumv
    COMPLEX(real_8),OPTIONAL                 :: psi(:)

    INTEGER                                  :: i, ia, ierr, ig,  ijv, &
                                                is, isa, isa0, isub, iv, jv, ii, &
                                                isub1, il_ctmp(2), il_deltar(1), &
                                                il_qg(1), il_ddia(2), il_fnlt(3), &
                                                methread, nhh0, nhh, ia_sum, &
                                                na_grp(2,ions1%nsp,parai%cp_nogrp),num_orb
    REAL(real_8)                             :: fac, sum, dtmp
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: ddia(:,:), fnlt(:,:,:)
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: ctmp(:,:), deltar(:), qg(:)
#else
    REAL(real_8), ALLOCATABLE                :: ddia(:,:), fnlt(:,:,:)
    COMPLEX(real_8), ALLOCATABLE             :: ctmp(:,:), deltar(:), qg(:)
#endif
    LOGICAL                                  :: zero
    CHARACTER(*), PARAMETER                  :: procedureN='rhov'
    CALL tiset(procedureN,isub)

    ! split atoms between cp groups
    call cp_grp_split_atoms(na_grp)
    ! generate 'new' nst12 mapping -> nwa12
    CALL dist_size(is2-is1+1,parai%nproc,paraw%nwa12,nblock=1,nbmax=num_orb,fw=1)
    !shift everything to starting state
    DO i=0,parai%nproc-1
       paraw%nwa12(1,i)=paraw%nwa12(1,i)+is1-1
       paraw%nwa12(2,i)=paraw%nwa12(2,i)+is1-1
    END DO
    num_orb=paraw%nwa12(2,parai%mepos)-paraw%nwa12(1,parai%mepos)+1
    CALL setfftn(0)
    IF (cntl%tfdist) CALL stopgm('RHOV','TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (imagp.EQ.2)CALL stopgm('RHOV','K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
   
!    ti=0.0d0
    nhh0=0
    do is=1,ions1%nsp
       if ( pslo_com%tvan(is) ) then
          nhh0=nhh0+nlps_com%ngh(is)*(nlps_com%ngh(is)+1)/2
       end if
    end do

    il_deltar(1)=ncpw%nhg
    il_qg(1)=ncpw%nhg
    IF(cntl%bigmem)THEN
       il_qg(1)=ncpw%nhg_cp
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_deltar,deltar,procedureN//'_deltar')
    CALL request_scratch(il_qg,qg,procedureN//'_qg')
#else
    ALLOCATE(deltar(il_deltar(1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(qg(il_qg(1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
#endif
   ! ==--------------------------------------------------------------==
    CALL zeroing(deltar)!,nhg)
    isa0=0
    IF (cntl%bigmem) THEN

       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is)) THEN

             nhh=(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2
             il_ctmp(1)=ncpw%nhg_cp
             il_ctmp(2)=nhh
             il_ddia(1)=ions0%na(is)
             il_ddia(2)=nhh
             il_fnlt(1)=num_orb
             il_fnlt(2)=nlps_com%ngh(is)
             il_fnlt(3)=parai%ncpus                      
#ifdef _USE_SCRATCHLIBRARY
             CALL request_scratch(il_ctmp,ctmp,procedureN//'_ctmp')
             CALL request_scratch(il_ddia,ddia,procedureN//'_ddia')
             CALL request_scratch(il_fnlt,fnlt,procedureN//'_fnlt')
#else
             ALLOCATE(ctmp(il_ctmp(1),il_ctmp(2)),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)
             ALLOCATE(ddia(il_ddia(1),il_ddia(2)),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)
             ALLOCATE(fnlt(il_fnlt(1),il_fnlt(2),il_fnlt(3)),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)
#endif
!             ti(1)=m_walltime()
             !$omp parallel private(methread,ijv,ia,isa,iv,i,ii,jv,dtmp,fac)
             methread=1
             !$ methread=omp_get_thread_num()+1
             !$omp do 
             DO ijv=1,nhh
                DO ia=1,ions0%na(is)
                   ddia(ia,ijv)=0.0_real_8
                END DO
             END DO
             ia_sum=na_grp(2,is,parai%cp_inter_me +1)-na_grp(2,is,parai%cp_inter_me +1)+1
             IF(ia_sum.GT.0) THEN
                !$omp do 
                DO ia=na_grp(1,is,parai%cp_inter_me +1),na_grp(2,is,parai%cp_inter_me +1)
                   isa=isa0+ia
                   !$omp simd
                   DO iv=1,nlps_com%ngh(is)
                      DO i=paraw%nwa12(1,parai%mepos),paraw%nwa12(2,parai%mepos)
                         ii=i-paraw%nwa12(1,parai%mepos)+1
                         fnlt(ii,iv,methread)=fnla(isa,iv,i)
                      END DO
                   END DO
                   ijv=0
                   DO iv=1,nlps_com%ngh(is)
                      DO jv=iv,nlps_com%ngh(is)
                         ijv=ijv+1
                         dtmp=0._real_8
                         !$omp simd reduction(+:dtmp)
                         DO i=paraw%nwa12(1,parai%mepos),paraw%nwa12(2,parai%mepos)
                            ii=i-paraw%nwa12(1,parai%mepos)+1
                            dtmp=dtmp+crge%f(i,1)*fnlt(ii,iv,methread)*fnlt(ii,jv,methread)
                         END DO
                         fac=1.0_real_8
                         IF (iv.NE.jv) fac=2.0_real_8
                         ddia(ia,ijv)=dtmp*fac
                      ENDDO
                   ENDDO
                END DO
             END IF
             !$omp end parallel
             CALL mp_sum(ddia,ions0%na(is)*nhh,parai%cp_grp)
             CALL dgemm('n','n',2*ncpw%nhg_cp,nhh,ions0%na(is) &
                  ,1.0d0,EIGRB(ncpw%nhg_start,isa0+1),2*ncpw%nhg &
                  ,DDIA(1,1),ions0%na(is),0.0_real_8,CTMP(1,1),2*ncpw%nhg_cp)
             ijv=0
             DO iv=1,nlps_com%ngh(is)
                DO jv=iv,nlps_com%ngh(is)
                   ijv=ijv+1
                   CALL qvan2(iv,jv,is,qg(:))
                   !$omp parallel do private(IG)
                   !CDIR NODEP
                   DO ig=1,ncpw%nhg_cp
                      deltar(ncpw%nhg_start-1+ig)=deltar(ncpw%nhg_start-1+ig)+qg(ig)*ctmp(ig,ijv)
                   ENDDO
                END DO
             END DO
#ifdef _USE_SCRATCHLIBRARY
             CALL free_scratch(il_fnlt,fnlt,procedureN//'_fnlt')
             CALL free_scratch(il_ddia,ddia,procedureN//'_ddia')
             CALL free_scratch(il_ctmp,ctmp,procedureN//'_ctmp')
#else
             DEALLOCATE(ctmp,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)
             DEALLOCATE(ddia,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)
             DEALLOCATE(fnlt,STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)
#endif
          END IF
          isa0=isa0+ions0%na(is)
       end do
       IF (parai%CP_NOGRP .GT. 1 .AND. cntl%nonort .and. pslo_com%tivan) then
          CALL TISET(procedureN//'_grpsb',isub1)
          CALL cp_grp_redist_array(deltar,ncpw%nhg)
          CALL TIHALT(procedureN//'_grpsb',isub1)
       END IF

    ELSE
       
       il_ctmp(1)=ncpw%nhg
       il_ctmp(2)=1
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_ctmp,ctmp,procedureN//'_ctmp')
#else
       ALLOCATE(ctmp(il_ctmp(1),il_ctmp(2)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
#endif
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is)) THEN
             DO iv=1,nlps_com%ngh(is)
                DO jv=iv,nlps_com%ngh(is)
                   CALL qvan2(iv,jv,is,qg(:))
                   DO ia=1,ions0%na(is)
                      isa=isa0+ia
                      !$omp parallel do private(IG)
#ifdef __SR8000
                      !poption parallel
#endif
                      DO ig=1,ncpw%nhg
                         ctmp(ig,1)=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                              ei3(isa,inyh(3,ig))
                      ENDDO
                      sum=0.0_real_8
                      !$omp parallel do private(I) reduction(+:SUM)
                      DO i=is1,is2
                         sum = sum +&
                              crge%f(i,1)*fnl(1,isa,iv,i,1)*fnl(1,isa,jv,i,1)
                      ENDDO
                      IF (iv.NE.jv) sum=2._real_8*sum
                      !$omp parallel do private(IG)
                      !CDIR NODEP
                      DO ig=1,ncpw%nhg
                         deltar(ig)=deltar(ig)+qg(ig)*sum*ctmp(ig,1)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          isa0=isa0+ions0%na(is)
       ENDDO
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_ctmp,ctmp,procedureN//'_ctmp')
#else
       DEALLOCATE(ctmp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
#endif
    ENDIF
    IF (geq0) THEN
       rsumv=REAL(deltar(1))
    ELSE
       rsumv=0.0_real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (PRESENT(psi)) THEN
       CALL zeroing(psi)!,maxfft)
       !$omp parallel do private(IG) shared(PSI)
       !CDIR NODEP
#ifdef __SR8000
       !poption parallel
#endif
       DO ig=1,ncpw%nhg
          psi(nzh(ig))=deltar(ig)
          psi(indz(ig))=CONJG(deltar(ig))
       ENDDO
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_qg,qg,procedureN//'_qg')
    CALL free_scratch(il_deltar,deltar,procedureN//'_deltar')
#else
    DEALLOCATE(qg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(deltar,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#endif
    IF (PRESENT(psi)) THEN
       CALL invfftn(psi, .FALSE.,parai%allgrp)
    END IF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rhov
  ! ==================================================================
  SUBROUTINE give_scr_rhov(lrhov,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrhov
    CHARACTER(len=30)                        :: tag

    lrhov=6*ncpw%nhg+maxsys%nax
    tag='6*NHG+maxsys%nax'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rhov
  ! ==================================================================

END MODULE rhov_utils
