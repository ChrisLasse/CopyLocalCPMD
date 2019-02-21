#include "cpmd_global.h"

MODULE csmat_utils
  USE cvan,                            ONLY: qq
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE fnl_utils,                       ONLY: split_atoms_btw_buffers
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: nlps_com
  USE nort,                            ONLY: nort_com,nort_ovlap
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE spin,                            ONLY: spin_mod
  USE summat_utils,                    ONLY: summat
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxsys,&
                                             ncpw,&
                                             paraw,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE
  INTEGER,ALLOCATABLE,SAVE                 :: na_buff(:,:,:,:)

  PRIVATE

  PUBLIC :: csmat

CONTAINS

  ! ==================================================================
  ! == FOR TKPNT=.TRUE. FNL IS COMPLEX -- CSMAT_C WILL BE WRITTEN   ==
  ! ==================================================================
  SUBROUTINE csmat(a,c0,fnl,nstate,ikind,store_nonort,parent,rdst)
    ! ==--------------------------------------------------------------==
    ! ==         COMPUTES THE OVERLAP MATRIX A = < C0 |S| C0 >        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8),INTENT(IN)                  :: fnl(ions1%nat,maxsys%nhxs,nstate,*)
    COMPLEX(real_8),INTENT(IN)               :: c0(ncpw%ngw,nstate)
    REAL(real_8),INTENT(OUT)                 :: a(nstate,*)
    INTEGER,INTENT(IN)                       :: ikind, nstate
    INTEGER                                  :: i, ia, is, isa, isa0, isub, iv,&
                                                j, jv, ierr, count, ia0, ia_start,&
                                                ia_end, ia_sum, nmin, nchunk, nmin_2,&
                                                nchunk_2, il_fnlat(3), il_fnlatj(4) 
    LOGICAL                                  :: second,only_parent,redist
    REAL(real_8)                             :: fractions(parai%nproc),selem, temp
    CHARACTER(*), PARAMETER                  :: procedureN = 'csmat' 
    REAL(real_8),ALLOCATABLE                 :: fnlat(:,:,:),fnlatj(:,:,:)    
    LOGICAL,INTENT(IN)                       :: store_nonort
    LOGICAL,INTENT(IN),OPTIONAL              :: parent,rdst
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)

    IF (cntl%tfdist) CALL stopgm('CSMAT','TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)

    IF(PRESENT(parent))THEN
       only_parent=parent
    ELSE
       only_parent=.FALSE.
    END IF
    
    IF(PRESENT(rdst))THEN
       redist=rdst
    ELSE
       redist=.TRUE.
    END IF

  
    !if needed allocate nort_ovlap for later use in crotwf
    IF (store_nonort) THEN
       IF (.NOT. ALLOCATED(nort_ovlap)) THEN
          ALLOCATE(nort_ovlap(nstate,nstate), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate nort_ovlap',& 
               __LINE__,__FILE__)
          !$omp parallel do private (i,j)
          DO i=1,nstate
             DO j=1,nstate
                nort_ovlap(j,i)=0.0_real_8
             END DO
          END DO
       END IF      
    END IF

    !we only need the upper part here
    CALL ovlap(nstate,a,c0,c0,rdst=redist,full=.FALSE.)
    IF (store_nonort) THEN
       !if needed save nort_ovlap for later use in crotwf
       IF(cntl%tlsd)THEN
          !$omp parallel private(i,j) 
          !$omp do 
          DO i=1,spin_mod%nsup
             DO j=1,i
                nort_ovlap(j,i)=a(j,i)
             END DO
          END DO
          !$omp end do nowait
          !$omp  do
          DO i=spin_mod%nsup+1,nstate
             DO j=spin_mod%nsup+1,i
                nort_ovlap(j,i)=a(j,i)
             END DO
          END DO
          !$omp end parallel
       ELSE
          !$omp parallel do private(i,j) schedule(static,1)
          DO i=1,nstate
             DO j=1,i
                nort_ovlap(j,i)=a(j,i)
             END DO
          END DO
       END IF
#ifdef _HAS_LIBELPA
       !elpa always needs the full matrix
       !if ovlap already redistributed, sum over allgrp
       IF (redist) THEN
          CALL summat(nort_ovlap,nstate,symmetrization=.TRUE.,lsd=.TRUE.,gid=parai%allgrp,&
               parent=.FALSE.)
       ELSE
          CALL summat(nort_ovlap,nstate,symmetrization=.TRUE.,lsd=.TRUE.,gid=parai%cp_grp,&
               parent=.FALSE.)
       END IF
#else
       !lapack only needs upper or lower part
       !if ovlap already redistributed, sum over allgrp
       IF (redist) THEN
          CALL summat(nort_ovlap,nstate,symmetrization=.FALSE.,lsd=.TRUE.,gid=parai%allgrp,&
               parent=.FALSE.)
       ELSE
          CALL summat(nort_ovlap,nstate,symmetrization=.FALSE.,lsd=.TRUE.,gid=parai%cp_grp,&
               parent=.FALSE.)
       END IF
#endif
       !always check threshold.
       !needs rethinking, in spin case there should be two seperate thresholds...
       temp=0.0_real_8
       IF(cntl%tlsd)THEN
          !$omp parallel private(i,j,selem)reduction(max:temp)
          !$omp do
          DO i=1,spin_mod%nsup
             DO j=1,i-1
                selem=ABS(nort_ovlap(j,i))
                IF (i.EQ.j) selem=ABS(nort_ovlap(j,i)-1.0_real_8)
                IF (temp.LT.selem) temp=selem
             ENDDO
          ENDDO
          !$omp end do nowait
          !$omp do
          DO i=spin_mod%nsup+1,nstate
             DO j=spin_mod%nsup+1,i-1
                selem=ABS(nort_ovlap(j,i))
                IF (i.EQ.j) selem=ABS(nort_ovlap(j,i)-1.0_real_8)
                IF (temp.LT.selem) temp=selem
             ENDDO
          ENDDO
          !$omp end parallel         
       ELSE
          !$omp parallel do private(i,j,selem)reduction(max:temp) schedule(static,8)
          DO i=1,nstate
             DO j=1,i-1
                selem=ABS(nort_ovlap(j,i))
                IF (i.EQ.j) selem=ABS(nort_ovlap(j,i)-1.0_real_8)
                IF (temp.LT.selem) temp=selem
             ENDDO
          ENDDO
       END IF
       nort_com%scond=temp
    END IF
        
    IF (pslo_com%tivan) then
       !preparation
       IF (.NOT. ALLOCATED(na_buff))THEN
          !distribute atoms between procs, take into account cp_groups, filter out any non 
          !uspp atom
          ALLOCATE(na_buff(2,ions1%nsp,parai%nproc,parai%cp_nogrp), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate na_buff',& 
               __LINE__,__FILE__)
          fractions=1.0_real_8/REAL(parai%nproc,KIND=real_8)
          CALL split_atoms_btw_buffers(parai%nproc,fractions,na_buff,only_uspp=.TRUE.)
       END IF

       !spin settings
       nmin=1
       nchunk=nstate
       second=.FALSE.
       nmin_2=0
       nchunk_2=0
       IF (cntl%tlsd) THEN
          nmin=1
          nchunk=spin_mod%nsup
          second=.TRUE.
          nmin_2=spin_mod%nsup+1
          nchunk_2=spin_mod%nsdown
       ENDIF
       !end preparation

       isa0=0
       DO count=1,maxsys%nhxs
          ia_sum=0
          !count all atoms with ngh(is) .EQ. count
          DO is=1,ions1%nsp
             IF (nlps_com%ngh(is).NE.count) CYCLE
             IF (.NOT.pslo_com%tvan(is)) CYCLE
             ia_sum=ia_sum+na_buff(2,is,parai%mepos+1,parai%cp_inter_me +1)&
                  -na_buff(1,is,parai%mepos+1,parai%cp_inter_me +1)+1
          END DO
          IF (ia_sum.EQ.0) CYCLE
          !get work space for this group of atoms

          il_fnlat(1)=ia_sum
          il_fnlat(2)=count
          il_fnlat(3)=nstate

          il_fnlatj(1)=ia_sum
          il_fnlatj(2)=count
          il_fnlatj(3)=nstate

          ALLOCATE(fnlat(il_fnlat(1),il_fnlat(2),il_fnlat(3)), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlati',& 
               __LINE__,__FILE__)
          ALLOCATE(fnlatj(il_fnlatj(1),il_fnlatj(2),il_fnlatj(3)), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlatj',& 
               __LINE__,__FILE__)

          !$omp parallel private(i,j,is,jv,iv,isa,ia_start,ia_end,ia,ia0,isa0)

          !prepare fnl chunk
          ia0=0
          DO is=1,ions1%nsp
             IF (nlps_com%ngh(is).NE.count) CYCLE
             IF (.NOT.pslo_com%tvan(is)) CYCLE
             isa0=0
             DO i=1,is-1
                isa0=isa0+ions0%na(i)
             END DO
             ia_start=ia0+1
             ia_end=ia0+na_buff(2,is,parai%mepos+1,parai%cp_inter_me +1)&
                        -na_buff(1,is,parai%mepos+1,parai%cp_inter_me +1) +1
             IF(ia_end-ia_start+1.GT.0) THEN
                !$omp do schedule(static)
                DO i=1,nstate
                   DO iv=1,nlps_com%ngh(is)
                      isa=isa0+na_buff(1,is,parai%mepos+1,parai%cp_inter_me +1) -1
                      DO ia=ia_start,ia_end
                         isa=isa+1
                         fnlat(ia,iv,i)=fnl(isa,iv,i,ikind)
                      END DO
                   END DO
                END DO
                !$omp end do nowait
                ia0=ia0+na_buff(2,is,parai%mepos+1,parai%cp_inter_me +1)&
                     -na_buff(1,is,parai%mepos+1,parai%cp_inter_me +1) +1
             END IF
          END DO
          !prepare fnl chunk, 'j' states
          !$omp do schedule(static)
          DO j=1,nstate
             DO iv=1,count
                DO ia=1,ia_sum
                   fnlatj(ia,iv,j)=0.0_real_8
                END DO
             END DO
          END DO
          !$omp end do nowait
          ia0=0
          DO is=1,ions1%nsp
             IF (nlps_com%ngh(is).NE.count) CYCLE
             IF (.NOT.pslo_com%tvan(is)) CYCLE
             ia_start=ia0+1
             ia_end=ia0+na_buff(2,is,parai%mepos+1,parai%cp_inter_me +1)&
                              -na_buff(1,is,parai%mepos+1,parai%cp_inter_me +1) +1
             IF(ia_end-ia_start+1.GT.0) THEN
                !$omp do schedule(static)
                DO j=1,nstate
                   DO iv=1,count
                      DO jv=1,count
                         IF (ABS(qq(jv,iv,is)).GT.1.e-5_real_8) THEN
                            DO ia=ia_start,ia_end
                               fnlatj(ia,iv,j)=fnlatj(ia,iv,j)+fnlat(ia,jv,j)*qq(jv,iv,is)
                            END DO
                         END IF
                      END DO
                   END DO
                END DO
                !$omp end do nowait
                ia0=ia0+na_buff(2,is,parai%mepos+1,parai%cp_inter_me +1)&
                     -na_buff(1,is,parai%mepos+1,parai%cp_inter_me +1) +1
             END IF
          END DO

          !$omp end parallel

          CALL dgemm('T','N',nchunk,nchunk,ia_sum*count,1.0_real_8,&
               fnlat(1,1,nmin),ia_sum*count,fnlatj(1,1,nmin),ia_sum*count,1.0_real_8,&
               a(nmin,nmin),nstate)

          IF (second) THEN
             CALL dgemm('T','N',nchunk_2,nchunk_2,ia_sum*count,1.0_real_8,&
                  fnlat(1,1,nmin_2),ia_sum*count,fnlatj(1,1,nmin_2),ia_sum*count,1.0_real_8,&
                  a(nmin_2,nmin_2),nstate)
          END IF
          DEALLOCATE(fnlat,fnlatj, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fnlat,fnlatj',& 
               __LINE__,__FILE__)
       END DO

       IF(redist)THEN
          CALL summat(a,nstate,symmetrization=.FALSE.,lsd=.TRUE.,gid=parai%allgrp,parent=only_parent)
       ELSE
          CALL summat(a,nstate,symmetrization=.FALSE.,lsd=.TRUE.,gid=parai%cp_grp,parent=only_parent)
       END IF
    END IF

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE csmat
  ! ==================================================================

END MODULE csmat_utils
