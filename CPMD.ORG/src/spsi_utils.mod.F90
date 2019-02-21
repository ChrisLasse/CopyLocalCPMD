#include "cpmd_global.h"

MODULE spsi_utils
  USE cp_grp_utils,                    ONLY: cp_grp_redist,&
                                             cp_grp_split_atoms  
  USE cppt,                            ONLY: twnl
  USE cvan,                            ONLY: qq
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlps_com
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: eigr,&
                                             fnla
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             maxsys,&
                                             ncpw,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: spsi
  PUBLIC :: give_scr_spsi

CONTAINS

  ! ==================================================================
  SUBROUTINE spsi(nstate,sc0)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8), TARGET                  :: sc0(ncpw%ngw,nstate)
    CHARACTER(*), PARAMETER                  :: procedureN = 'spsi'

    COMPLEX(real_8)                          :: ci
    COMPLEX(real_8), ALLOCATABLE             :: work(:,:)
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: temp(:,:)
    REAL(real_8), ALLOCATABLE                :: aa(:,:)
    INTEGER                                  :: ia, ierr, ig, is, isa, isa0, &
                                                isub, iv, i, jv, idx, niv, ld_aa, &
                                                iac, iac0, il_work(2),il_temp(2),il_aa(2),&
                                                na_grp(2,ions1%nsp,parai%cp_nogrp),&
                                                revcnt(parai%nproc),displs(parai%nproc)
    REAL(real_8)                             :: t1, t2, t3

    CALL tiset(procedureN,isub)
    IF (cntl%tfdist) CALL stopgm('SPSI','TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (imagp.EQ.2) CALL stopgm('SPSI','K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==

    ! split atoms between cp_groups
    CALL cp_grp_split_atoms(na_grp)   
    
    ld_aa=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is) ) THEN
          ld_aa=ld_aa+ &
               (na_grp(2,is,parai%cp_inter_me +1)-na_grp(1,is,parai%cp_inter_me +1) +1)*nlps_com%ngh(is)
       END IF
    END DO

    il_temp(1)=ncpw%ngw
    il_temp(2)=nstate
    il_aa(1)=ld_aa
    il_aa(2)=nstate
    il_work(1)=ncpw%ngw
    il_work(2)=ld_aa

    IF(parai%cp_nogrp.gt.1)THEN
       ALLOCATE(temp(il_temp(1),il_temp(2)), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate temp',& 
            __LINE__,__FILE__)
       !$omp parallel do private(i,ig)
       DO i=1,nstate
          DO ig=1,ncpw%ngw
             temp(ig,i)=(0.0_real_8,0.0_real_8)
          END DO
       END DO
    ELSE
       temp=>sc0   
    END IF
    ALLOCATE(aa(il_aa(1),il_aa(2)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate aa',& 
         __LINE__,__FILE__)
    ALLOCATE(work(il_work(1),il_work(2)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate work',& 
         __LINE__,__FILE__)

    !this loop may be parallelized using nst12, however, with reducing work by skipping qq=0 
    !this seems not to be benificial at all

    !$omp parallel private(i,isa0,iac0,is,iv,iac,isa,ia,ci,ig,t1,t2,t3)     
    !$omp do 
    DO i=1,nstate
       DO iac=1,ld_aa
          aa(iac,i)=0.0_real_8
       END DO
       isa0=0
       iac0=0
       DO is=1,ions1%nsp
          DO iv=1,nlps_com%ngh(is)
             DO jv=1,nlps_com%ngh(is)
                IF (ABS(qq(jv,iv,is)).GT.1.e-5_real_8) THEN
                   iac=iac0
                   isa=isa0+na_grp(1,is,parai%cp_inter_me +1)-1
                   DO ia=1,na_grp(2,is,parai%cp_inter_me +1)-na_grp(1,is,parai%cp_inter_me +1) +1
                      iac=iac+1
                      isa=isa+1
                      aa(iac,i)=aa(iac,i)+qq(iv,jv,is)*fnla(isa,jv,i)
                   END DO
                END IF
             END DO
             iac0=iac0+na_grp(2,is,parai%cp_inter_me +1)-na_grp(1,is,parai%cp_inter_me +1) +1
          END DO
          isa0=isa0+ions0%na(is)            
       END DO
    END DO
    !$omp end do nowait

    iac0=0
    isa0=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          DO iv=1,nlps_com%ngh(is)
             ci=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
             !       MAKE USE OF THE SPECIAL STRUCTURE OF CI
             IF (ABS(REAL(ci)).GT.0.5_real_8) THEN
                !           CI IS REAL               
                !$omp do 
                DO ia=na_grp(1,is,parai%cp_inter_me +1),na_grp(2,is,parai%cp_inter_me +1)
                   isa=isa0+ia
                   iac=iac0+ia-na_grp(1,is,parai%cp_inter_me +1)+1
                   DO ig=1,ncpw%ngw
                      t1=REAL(eigr(ig,isa,1))
                      t2=AIMAG(eigr(ig,isa,1))
                      t3=twnl(ig,iv,is,1)*REAL(ci)
                      work(ig,iac)=CMPLX(t1*t3,t2*t3,kind=real_8)
                   END DO
                END DO
                !$omp end do nowait
             ELSE
                !          CI IS IMAGINARY
                !$omp do
                DO ia=na_grp(1,is,parai%cp_inter_me +1),na_grp(2,is,parai%cp_inter_me +1)
                   isa=isa0+ia
                   iac=iac0+ia-na_grp(1,is,parai%cp_inter_me +1)+1
                   DO ig=1,ncpw%ngw
                      t1=REAL(eigr(ig,isa,1))
                      t2=AIMAG(eigr(ig,isa,1))
                      t3=twnl(ig,iv,is,1)*AIMAG(ci)
                      work(ig,iac)=CMPLX(-t2*t3,t1*t3,kind=real_8)
                   END DO
                END DO
                !$omp end do nowait
             END IF
             iac0=iac0+na_grp(2,is,parai%cp_inter_me +1)-na_grp(1,is,parai%cp_inter_me +1) +1
          END DO
       END IF
       isa0=isa0+ions0%na(is)
    END DO
    !$omp end parallel 

    IF(ncpw%ngw /= 0) THEN
        CALL DGEMM('n','n',2*ncpw%ngw,nstate, &
            ld_aa, 1.0_real_8,work(1,1),2*ncpw%ngw,&
            aa(1,1),ld_aa,&
            1.0_real_8,temp(1,1),2*ncpw%ngw)
    ENDIF
    IF (parai%cp_nogrp .gt. 1) THEN
       CALL cp_grp_redist(temp,ncpw%ngw,nstate)
       !$omp parallel do private(i,ig)
       DO i=1,nstate
          DO ig=1,ncpw%ngw
             sc0(ig,i)=sc0(ig,i)+temp(ig,i)
          END DO
       END DO
       DEALLOCATE(temp, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate temp',& 
            __LINE__,__FILE__)
    END IF
    DEALLOCATE(aa, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate aa',& 
         __LINE__,__FILE__)
    DEALLOCATE(work, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate work',& 
         __LINE__,__FILE__)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE spsi
  ! ==================================================================
  SUBROUTINE give_scr_spsi(lspsi,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lspsi
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    lspsi=ncpw%ngw*2
    tag ='NGW*2'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_spsi
  ! ==================================================================

END MODULE spsi_utils
