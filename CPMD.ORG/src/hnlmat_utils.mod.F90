#include "cpmd_global.h"

MODULE hnlmat_utils
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms
  USE cvan,                            ONLY: deeq,&
                                             dvan
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE fnl_utils,                       ONLY: split_atoms_btw_buffers
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai,paral
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: fnla
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             parap,&
                                             maxsys
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE
  INTEGER,ALLOCATABLE,SAVE                 :: na_buff(:,:,:,:)

  PRIVATE

  PUBLIC :: hnlmat

CONTAINS
  ! ==================================================================
  SUBROUTINE hnlmat(hmat,f,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate
    REAL(real_8),INTENT(IN)                  :: f(nstate)
    REAL(real_8),INTENT(INOUT),TARGET        :: hmat(nstate,nstate)
    INTEGER                                  :: i, ia, ind, ind1, is, isa, isa0, ispin, isub,&
                                                iv, j, jmax, jv, ierr, count, spin, ia_start,&
                                                ia_end, ia0, ia_sum, nmin, nchunk, nmin_2,&
                                                nchunk_2, il_fnlat(3), il_fnlatj(4),&
                                                il_deeqt(4), il_hmat_loc(2)
    LOGICAL                                  :: need_hmat_loc, non_uspp, second
    REAL(real_8)                             :: ffi, fac, sum , fractions(parai%nproc)
    REAL(real_8),ALLOCATABLE                 :: fnlat(:,:,:), deeqt(:,:,:,:), fnlatj(:,:,:,:)
    REAL(real_8),POINTER __CONTIGUOUS        :: hmat_loc(:,:)
    CHARACTER(*), PARAMETER                  :: procedureN = 'hnlmat'
    ! ==--------------------------------------------------------------==
    ! == Compute the non-local contribution to the Hamilton matrix    ==
    ! == optimized version for uspp only, calls old hnlmat in case of ==
    ! == other pseudo potentials                                      ==
    ! ==--------------------------------------------------------------==
    
    CALL tiset(procedureN,isub)
    IF (cntl%tfdist) CALL stopgm('HNLMAT','TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (imagp.EQ.2)CALL stopgm('HNLMAT','K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    
    IF (pslo_com%tivan) THEN
       
       ! VANDERBILT PP
       spin=1
       nmin=1
       nchunk=nstate
       fac=-2.0_real_8
       second=.FALSE.
       nmin_2=0
       nchunk_2=0
       IF (cntl%tlsd) THEN
          spin=2   
          nmin=1
          nchunk=spin_mod%nsup
          fac=-1.0_real_8
          second=.TRUE.
          nmin_2=spin_mod%nsup+1
          nchunk_2=spin_mod%nsdown
       ENDIF

       IF (.NOT. ALLOCATED(na_buff))THEN
          !distribute atoms between procs, take into account cp_group distributed fnl, filter out any non uspp atom
          ALLOCATE(na_buff(2,ions1%nsp,parai%nproc,parai%cp_nogrp), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate na_buff',& 
               __LINE__,__FILE__)
          fractions=1.0_real_8/REAL(parai%nproc,KIND=real_8)
          CALL split_atoms_btw_buffers(parai%nproc,fractions,na_buff,only_uspp=.TRUE.)
       END IF

       !do we have equal occupation numbers?
       need_hmat_loc=.FALSE.
       ffi=f(1)
       DO i=1,nstate
          IF (ffi.NE.f(i)) THEN
             need_hmat_loc=.TRUE.
             EXIT
          ENDIF
       ENDDO
       IF(need_hmat_loc)THEN
          il_hmat_loc(1)=nstate
          il_hmat_loc(2)=nstate

          ALLOCATE(hmat_loc(il_hmat_loc(1),il_hmat_loc(2)), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate hmat_loc',& 
               __LINE__,__FILE__)
          !$omp parallel do private(i,j)
          DO i=1,nstate
             DO j=1,nstate
                hmat_loc(j,i)=0.0_real_8
             END DO
          END DO
          fac=-1.0_real_8
       ELSE
          hmat_loc=>hmat
       END IF

       isa0=0
       DO count=1,maxsys%nhxs
          ia_sum=0
          !count all atoms with ngh(is) .EQ. count
          DO is=1,ions1%nsp
             IF (nlps_com%ngh(is).NE.count) CYCLE
             IF (.NOT. pslo_com%tvan(is)) CYCLE
             ia_sum=ia_sum+na_buff(2,is,parai%mepos+1,parai%cp_inter_me +1)&
                  -na_buff(1,is,parai%mepos+1,parai%cp_inter_me +1)+1
          END DO
          IF (ia_sum .EQ. 0) CYCLE
          !get work space for this group of atoms

          il_fnlat(1)=ia_sum
          il_fnlat(2)=count
          il_fnlat(3)=nstate

          il_fnlatj(1)=ia_sum
          il_fnlatj(2)=count
          il_fnlatj(3)=nstate
          il_fnlatj(4)=spin

          il_deeqt(1)=ia_sum
          il_deeqt(2)=count
          il_deeqt(3)=count
          il_deeqt(4)=spin

          ALLOCATE(fnlat(il_fnlat(1),il_fnlat(2),il_fnlat(3)), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlat',& 
               __LINE__,__FILE__)
          ALLOCATE(fnlatj(il_fnlatj(1),il_fnlatj(2),il_fnlatj(3),il_fnlatj(4)), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlatj',& 
               __LINE__,__FILE__)
          ALLOCATE(deeqt(il_deeqt(1),il_deeqt(2),il_deeqt(3),il_deeqt(4)), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate deeqt',& 
               __LINE__,__FILE__)
          !$omp parallel private(i,j,ispin,is,jv,iv,isa,ia,ia_end,ia_start,ia0,isa0)

          !fill work space
          ia0=0
          DO is=1,ions1%nsp
             IF (nlps_com%ngh(is).NE.count) CYCLE
             IF (.NOT. pslo_com%tvan(is)) CYCLE
             isa0=0         
             DO i=1,is-1
                isa0=isa0+ions0%na(i)
             END DO
             ia_start=ia0+1
             ia_end=ia0+na_buff(2,is,parai%mepos+1,parai%cp_inter_me +1)&
                  -na_buff(1,is,parai%mepos+1,parai%cp_inter_me +1) +1
             IF(ia_end-ia_start+1.GT.0) THEN
                DO ispin=1,spin
                   !$omp do
                   DO iv=1,count
                      DO jv=1,count
                         isa=isa0+na_buff(1,is,parai%mepos+1,parai%cp_inter_me +1) -1
                         DO ia=ia_start,ia_end
                            isa=isa+1
                            deeqt(ia,jv,iv,ispin)=deeq(isa,jv,iv,ispin)+dvan(jv,iv,is)
                         END DO
                      END DO
                   END DO
                   !$omp end do nowait
                END DO
                ia0=ia0+na_buff(2,is,parai%mepos+1,parai%cp_inter_me +1)&
                     -na_buff(1,is,parai%mepos+1,parai%cp_inter_me +1) +1
             END IF
          END DO
          !deeqt needs to be complete
          !$omp barrier
          !prepare fnl chunk
          ia0=0
          DO is=1,ions1%nsp
             IF (nlps_com%ngh(is).NE.count) CYCLE
             IF (.NOT. pslo_com%tvan(is)) CYCLE
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
                         fnlat(ia,iv,i)=fnla(isa,iv,i)
                      END DO
                   END DO
                END DO
                !$omp end do nowait
                ia0=ia0+na_buff(2,is,parai%mepos+1,parai%cp_inter_me +1)&
                     -na_buff(1,is,parai%mepos+1,parai%cp_inter_me +1) +1
             END IF
          END DO
          !prepare fnl chunk, 'j' states,
          DO ispin=1,spin
             !$omp do schedule(static)
             DO j=1,nstate
                DO iv=1,count
                   DO ia=1,ia_sum
                      fnlatj(ia,iv,j,ispin)=0.0_real_8
                   END DO
                END DO
                DO iv=1,count
                   DO jv=1,count
                      DO ia=1,ia_sum
                         fnlatj(ia,iv,j,ispin)=fnlatj(ia,iv,j,ispin)+deeqt(ia,jv,iv,ispin)*fnlat(ia,jv,j)
                      ENDDO
                   END DO
                END DO
             END DO
             !$omp end do nowait
          END DO
          !$omp end parallel

          ispin=1
          CALL dgemm('T','N',nchunk,nchunk,ia_sum*count,fac,&
               fnlat(1,1,nmin),ia_sum*count,fnlatj(1,1,nmin,ispin),ia_sum*count,1.0_real_8,&
               hmat_loc(nmin,nmin),nstate)
          IF (second) THEN
             ispin=2
             CALL dgemm('T','N',nchunk_2,nchunk_2,ia_sum*count,fac,&
                  fnlat(1,1,nmin_2),ia_sum*count,fnlatj(1,1,nmin_2,ispin),ia_sum*count,1.0_real_8,&
                  hmat_loc(nmin_2,nmin_2),nstate)
          END IF
          DEALLOCATE(fnlat, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fnlat',& 
               __LINE__,__FILE__)
          DEALLOCATE(fnlatj, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fnlatj',& 
               __LINE__,__FILE__)
          DEALLOCATE(deeqt, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate deeqt',& 
               __LINE__,__FILE__)
       END DO

       IF(need_hmat_loc)THEN
          !$omp parallel do private(i,j)
          DO i=1,nstate
             DO j=1,nstate
                IF (f(i) .GE. 1.e-5_real_8) THEN
                   hmat(j,i)=hmat(j,i)+hmat_loc(j,i)*f(i)
                ELSE 
                   hmat(j,i)=hmat(j,i)+hmat_loc(j,i)
                END IF
             END DO
          END DO
          DEALLOCATE(hmat_loc, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate deeqt',& 
               __LINE__,__FILE__)
       END IF
       
    END IF

    !!END VANDERBILT OPTIMIZED!!
    !lets check if there is something else...

    non_uspp=.FALSE.
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is))CYCLE
       non_uspp=.TRUE.
    END DO

    IF (non_uspp) CALL hnlmat_old(hmat,f,nstate)
    CALL tihalt(procedureN,isub)

  END SUBROUTINE hnlmat
  ! ==================================================================
  SUBROUTINE hnlmat_old(hmat,f,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate), hmat(nstate,nstate)

    INTEGER                                  :: chunk, i, ia, ind, ind1, is, &
                                                isa, isa0, ispin, isub, iv, &
                                                j, jmax, jv, ki, kj, l, l2, &
                                                li, lj, mxmypos, mxrank, mypos, &
                                                idx, npairs, &
                                                na_grp(2,ions1%nsp,parai%cp_nogrp)
    REAL(real_8)                             :: dd, fdd, ffi, sum
    INTEGER,ALLOCATABLE                      :: map(:,:)
    CHARACTER(*), PARAMETER                  :: procedureN = 'hnlmat_old'    
! Variables
! ==--------------------------------------------------------------==
! == Compute the non-local contribution to the Hamilton matrix    ==
! == fallback for non Vanderbilt species                          ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF (cntl%tfdist) CALL stopgm('HNLMAT','TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (imagp.EQ.2)CALL stopgm('HNLMAT','K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)

    ! split atoms between cp_groups
    CALL cp_grp_split_atoms(na_grp)   

    chunk = 0
    mxrank = 0
    IF (cntl%tlsd) THEN
       npairs=(spin_mod%nsup*(spin_mod%nsup+1)+(nstate-spin_mod%nsup)*((nstate-spin_mod%nsup)+1))/2 
    ELSE
       npairs=nstate*(nstate+1)/2
    ENDIF

    IF (parap%nst12(1,parai%mepos).LE.parap%nst12(2,parai%mepos)) THEN
       ind1 = 0
       DO i=0, parai%nproc-1
          IF (parap%nst12(1,i).LE.parap%nst12(2,i)) THEN
             ind1 = ind1 + 1
             IF (i.GT.mxrank) mxrank=i
             IF (parai%mepos.EQ.i) mxmypos = ind1
          ENDIF
       ENDDO

       chunk = npairs / ind1
       IF(parai%mepos.EQ.mxrank) chunk=npairs-(ind1-1)*chunk
       ind = 1
       mypos = (mxmypos-1)*chunk+1
       IF (parai%mepos.EQ.mxrank) THEN
          mypos=(ind1-1)*((npairs)/ind1)+1
       ENDIF
       allocate(map(2,chunk))
       IDX=0
       DO i=1, nstate
          jmax=nstate
          IF(cntl%tlsd.AND.i.LE.spin_mod%nsup) jmax=spin_mod%nsup
          DO j=i, jmax
             IF (ind.GE.mypos) THEN
                idx=idx+1
                map(1,idx)=i
                map(2,idx)=j
             ENDIF
             ind = ind +1
             IF(idx.GE.chunk) GOTO 999
          ENDDO
       ENDDO
999    CONTINUE

    ENDIF
   
    !$omp parallel do private(ind,i,j,ispin,sum,isa0,is,jv,iv,isa,ia,fdd,l,ki,li,l2,lj,kj,dd)
    DO ind = 1, chunk!chunk_start, chunk_end
       i=map(1,ind)
       j=map(2,ind)
       ispin=1
       IF(cntl%tlsd.AND.i.GT.spin_mod%nsup) ispin=2
       sum=0.0_real_8
       isa0=0
       DO is=1,ions1%nsp
          IF(pslo_com%tvan(is)) THEN
             ! VANDERBILT PP moved to new hnlmat above
          ELSEIF (sgpp1%tsgp(is)) THEN
             ! Stefan Goedecker PP
             DO iv=1,nlps_com%ngh(is)
                l=nghtol(iv,is)+1
                ki=sgpp2%lfval(iv,is)
                li=sgpp2%lpval(iv,is)
                DO jv=1,nlps_com%ngh(is)
                   l2=nghtol(jv,is)+1
                   lj=sgpp2%lpval(jv,is)
                   IF (l2.EQ.l.AND.li.EQ.lj) THEN
                      kj=sgpp2%lfval(jv,is)
                      DO ia=na_grp(1,is,parai%cp_inter_me +1),na_grp(2,is,parai%cp_inter_me +1)
                         isa=isa0+ia
                         fdd=sgpp2%hlsg(ki,kj,l,is)*fnla(isa,iv,i)*&
                              fnla(isa,jv,j)
                         sum=sum-fdd
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ELSE
             ! BACHELET HAMANN SCHLUTER
             DO iv=1,nlps_com%ngh(is)
                dd=wsg(is,iv)
                DO ia=na_grp(1,is,parai%cp_inter_me +1),na_grp(2,is,parai%cp_inter_me +1)
                   isa=isa0+ia
                   fdd=dd*fnla(isa,iv,i)*fnla(isa,iv,j)
                   sum=sum-fdd
                ENDDO
             END DO
          ENDIF
          isa0=isa0+ions0%na(is)
       ENDDO
       IF (f(i) .GE. 1.e-5_real_8) THEN
          hmat(j,i)=hmat(j,i)+sum*f(i)
          IF(i.NE.j) hmat(i,j)=hmat(i,j)+sum*f(i)
       ELSE
          hmat(j,i)=hmat(j,i)+sum
          IF(i.NE.j) hmat(i,j)=hmat(i,j)+sum
       END IF
    ENDDO
    IF (parap%nst12(1,parai%mepos).LE.parap%nst12(2,parai%mepos)) DEALLOCATE(map)

    CALL tihalt(procedureN,isub)

  END SUBROUTINE hnlmat_old
  ! ==================================================================

END MODULE hnlmat_utils
