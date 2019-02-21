#include "cpmd_global.h"

MODULE nlforce_utils
  USE cppt,                            ONLY: twnl
  USE cp_grp_utils,                    ONLY: cp_grp_redist,&
                                             cp_grp_split_atoms  
  USE cvan,                            ONLY: deeq,&
                                             dvan,&
                                             qq
  USE distribution_utils,              ONLY: dist_size
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE fnl_utils,                       ONLY: split_atoms_btw_buffers
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlps_com,&
                                             wsg
  !$ USE omp_lib,                      ONLY: omp_set_nested,&
  !$                                         omp_get_thread_num,&
  !$                                         omp_set_num_threads
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE reshaper,                        ONLY: reshape_inplace
  USE sfac,                            ONLY: eigr,&
                                             fnla
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             paraw,&
                                             norbpe
  USE tbxc,                            ONLY: toldcode
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nlforce
  PUBLIC :: give_scr_nlforce

CONTAINS

  ! ==================================================================
  SUBROUTINE nlforce(c2,f,gam,nstate,redist)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),TARGET                   :: c2(ncpw%ngw,nstate)
    REAL(real_8)                             :: f(*),gam(nstate,*)
    INTEGER                                  :: NSTATE
    LOGICAL,INTENT(IN),OPTIONAL              :: redist

    character(*), PARAMETER                  :: proceduren = 'nlforce'
    COMPLEX(real_8)                          :: ct, ci,ctm,cir,cii
    INTEGER                                  :: i, ia, ig, is, isa, isa0, iac, iac0, iv, jv, ld_dai,&
                                                nmax, nmin, nmax_2, nmin_2, nchunk, nchunk_2, &
                                                ispin, ki, kj, l, l2, li, lj, &
                                                buff, last_buff, methread, nthreads, nested_threads,&
                                                revcnt(parai%nproc), displs(parai%nproc),&
                                                isub, isub1, ierr, start_dai, start_eiscr, buffcount,&
                                                na_grp(2,ions1%nsp,parai%cp_nogrp),iac_start,iac_end,isa_start
    INTEGER,ALLOCATABLE,SAVE                 :: map_isa(:,:,:),na_buff(:,:,:,:)

    REAL(real_8)                             :: dd, ffi, t1, t2, t3, fractions(4)
    !$ LOGICAL                               :: locks(4,2)

    REAL(real_8),POINTER __CONTIGUOUS        :: dai_ptr(:,:)
    COMPLEX(real_8),POINTER __CONTIGUOUS     :: auxc(:,:)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: dai(:,:)
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: eiscr(:,:)
#else
    REAL(real_8),ALLOCATABLE                 :: dai(:,:)
    COMPLEX(real_8),ALLOCATABLE              :: eiscr(:,:)
#endif
    INTEGER                                  :: il_dai(2),il_eiscr(2),il_auxc(2)
    LOGICAL                                  :: need_scratch, redst, ready, second

    ! split atoms between cp groups
    CALL cp_grp_split_atoms(na_grp)

    IF (PRESENT(redist)) THEN
       redst=redist
    ELSE
       redst=.TRUE.
    END IF

    ! ==--------------------------------------------------------------==
    ! == compute the non-local contributions to the electron gradient ==
    ! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF (imagp.EQ.2) call stopgm(procedureN,'k-point not implemented',& 
         __LINE__,__FILE__)
    IF (cntl%tfdist) call stopgm(procedureN,'fnl dist. not implemented',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! fixme: ak 2008/02/14: this was broken by changes from 2006/05/12.
    ! not sure whether those were the real cause.
    ! undoing those changes does not make it work
    ! again. so we stop here until somebody fixes it
    ! update ak 2008/05/24: it looks as IF it can be worked around using 
    ! oldcode. it seems to be faster, too.
    IF (pslo_com%tivan.and.cntl%tlsd) THEN
       IF (.NOT.toldcode)&
            CALL stopgm(procedureN,'vanderbilt with lsd requires USE of'&
            // ' "oldcode" flag in &dft section.',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL dist_size(nstate,parai%nproc,paraw%nwa12,nblock=1,fw=1)
    IF (pslo_com%tivan) THEN
       ! ==--------------------------------------------------------==
       ! ==  vanderbilt pp                                         ==
       ! ==--------------------------------------------------------==
       ld_dai=0
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is) ) THEN
             ld_dai=ld_dai+ &
                  (na_grp(2,is,parai%cp_inter_me +1)-na_grp(1,is,parai%cp_inter_me +1) +1)*&
                  nlps_com%ngh(is)
          END IF
       END DO     

       il_dai(1)=ld_dai*nstate
       il_dai(2)=1
       il_eiscr(1)=nkpt%ngwk
       il_eiscr(2)=ld_dai
       il_auxc(1)=nkpt%ngwk
       il_auxc(2)=nstate

#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_dai,dai,procedureN//'_dai')
       CALL request_scratch(il_eiscr,eiscr,procedureN//'_eiscr')
#else
       ALLOCATE(dai(il_dai(1),il_dai(2)), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate dai',& 
            __LINE__,__FILE__)
       ALLOCATE(eiscr(il_eiscr(1),il_eiscr(2)), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate eiscr',& 
            __LINE__,__FILE__)
#endif
       IF(parai%cp_nogrp.GT.1)THEN
#ifdef _USE_SCRATCHLIBRARY
          CALL request_scratch(il_auxc,auxc,procedureN//'_auxc')
#else
          ALLOCATE(auxc(il_auxc(1),il_auxc(2)), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate auxc',& 
               __LINE__,__FILE__)
#endif
          !$omp parallel do private(i,ig)
          DO i=1,nstate
             DO ig=1,ncpw%ngw
                auxc(ig,i)=(0.0_real_8,0.0_real_8)
             END DO
          END DO
       ELSE
          auxc=>c2
       END IF


       !use overlapping com/comp
       buffcount=4
       
       IF (.NOT. ALLOCATED(map_isa) ) THEN

          ALLOCATE(map_isa(ions1%nat+1,buffcount,parai%cp_nogrp), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate map_isa',& 
               __LINE__,__FILE__)

          ALLOCATE(na_buff(2,ions1%nsp,buffcount,parai%cp_nogrp), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate na_buff',& 
               __LINE__,__FILE__)

          fractions=1.0_real_8/REAL(buffcount,KIND=real_8)

          CALL split_atoms_btw_buffers(buffcount,fractions,na_buff,map_isa,.true.)

       END IF

       !set up spin settings for dai dgemms
       IF (cntl%tlsd) THEN
          IF (paraw%nwa12(1,parai%mepos).LE.spin_mod%nsup.AND.&
               paraw%nwa12(2,parai%mepos).LE.spin_mod%nsup) THEN
             !only nsup states
             nmin=1
             nmax=spin_mod%nsup
             second=.FALSE.
             nchunk=norbpe
             nmin_2=0
             nmax_2=-1
             nchunk_2=0
          ELSE IF (paraw%nwa12(1,parai%mepos).GT.spin_mod%nsup.AND.&
               paraw%nwa12(2,parai%mepos).GT.spin_mod%nsup) THEN
             !only nsdown states
             nmin=spin_mod%nsup+1
             nmax=nstate
             nchunk=norbpe
             second=.FALSE.
             nmin_2=0
             nmax_2=-1
             nchunk_2=0
          ELSE IF (paraw%nwa12(1,parai%mepos).LE.spin_mod%nsup.AND.&
               paraw%nwa12(2,parai%mepos).GT.spin_mod%nsup) THEN
             !unfortunately nst12 coveres nsup and nsdown
             nmin=1
             nmax=spin_mod%nsup
             nchunk=spin_mod%nsup-paraw%nwa12(1,parai%mepos)+1
             nmin_2=spin_mod%nsup+1
             nmax_2=nstate
             nchunk_2=paraw%nwa12(2,parai%mepos)-(spin_mod%nsup+1)+1
             second=.TRUE.
          ENDIF
       ELSE
          nmin=1
          nmax=nstate
          nchunk=paraw%nwa12(2,parai%mepos)-paraw%nwa12(1,parai%mepos)+1
          second=.FALSE.
          nmin_2=0
          nmax_2=-1
          nchunk_2=0
       ENDIF

       !last buffer should always be done with all threads
       last_buff=buffcount
       !If com finished before last_buff is reached, exit parallel region
       ready=.FALSE.

       methread=0
#ifdef _USE_OVERLAPPING_COM_COMP
       nested_threads=MAX(parai%ncpus-1,1)
       nthreads=MIN(2,parai%ncpus)
#else
       nested_threads=1
       nthreads=parai%ncpus
#endif
       !$ locks=.TRUE.
       !$OMP parallel num_threads(nthreads) private(methread,buff,start_dai,start_eiscr,ld_dai,i) shared(locks,nthreads) proc_bind(close)
       !$ methread = omp_get_thread_num()
       !$ IF (nested_threads .NE. parai%ncpus) THEN
       !$    IF (methread.EQ.1.OR.nthreads.EQ.1) THEN
       !$       CALL omp_set_num_threads(nested_threads)
       !$       CALL omp_set_nested(1)
       !$    END IF
       !$ END IF
       !$OMP barrier
             
       IF (methread.EQ.1.OR.nthreads.EQ.1) THEN
          start_dai=0
          DO buff=1,buffcount
             ld_dai=map_isa(ions1%nat+1,buff,parai%cp_inter_me +1)
             IF (ld_dai.GT.0.AND.nchunk.GT.0) THEN
                call reshape_inplace(dai(start_dai+1:start_dai+ld_dai*nstate,1), (/ld_dai,nstate/),dai_ptr)
                !$omp parallel private (i,ffi,ispin,iac,iac0,isa0,is,iv,jv,isa,ia,iac_start,iac_end)
                !$omp do 
                DO i=paraw%nwa12(1,parai%mepos),paraw%nwa12(2,parai%mepos)
                   !setup spin settings
                   ffi=f(i)
                   IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
                   ispin=1
                   IF (cntl%tlsd.AND.i.GT.spin_mod%nsup) ispin=2

                   !clear dai
                   DO iac=1,ld_dai
                      dai_ptr(iac,i)=0.0_real_8
                   END DO

                   iac0=0
                   isa0=0
                   !fill local part of dai
                   DO is=1,ions1%nsp
                      IF(pslo_com%tvan(is)) THEN
                         ia=na_buff(2,is,buff,parai%cp_inter_me+1)-&
                              na_buff(1,is,buff,parai%cp_inter_me+1)+1
                         IF(ia.GT.0)THEN
                            DO iv=1,nlps_com%ngh(is)
                               DO jv=1,nlps_com%ngh(is)
                                  iac_start=iac0+1
                                  iac_end=iac0+ia
                                  isa=isa0+na_buff(1,is,buff,parai%cp_inter_me +1)-1
                                  DO iac=iac_start,iac_end
                                     isa=isa+1
                                     dai_ptr(iac,i)=dai_ptr(iac,i)&
                                          -ffi*fnla(isa,jv,i)*&
                                          (deeq(isa,jv,iv,ispin)+dvan(jv,iv,is))
                                  ENDDO
                               END DO
                               iac0=iac0+ia
                            END DO
                         END IF
                      END IF
                      isa0=isa0+ions0%na(is)
                   END DO
                END DO
                !$omp end do

                iac0=1
                isa0=0
                DO is=1,ions1%nsp
                   IF(pslo_com%tvan(is)) THEN
                      ia=(na_buff(2,is,buff,parai%cp_inter_me +1)-&
                           na_buff(1,is,buff,parai%cp_inter_me +1)+1)
                      IF(ia.GT.0)THEN
                         !$omp do schedule(dynamic)
                         DO iv=1,nlps_com%ngh(is)
                            iac=iac0+(iv-1)*ia
                            DO jv=1,nlps_com%ngh(is)
                               IF (ABS(qq(jv,iv,is)).GT.1.e-5_real_8) THEN
                                  CALL dgemm('N','N',ia,nchunk,nmax-nmin+1,-qq(jv,iv,is),&
                                       fnla(isa0+na_buff(1,is,buff,parai%cp_inter_me +1),jv,nmin),&
                                       ions1%nat*maxsys%nhxs, &
                                       gam(nmin,paraw%nwa12(1,parai%mepos)),nstate,1.0_real_8,&
                                       dai_ptr(iac,paraw%nwa12(1,parai%mepos)),ld_dai)
                                  IF (second) THEN
                                     CALL dgemm('N','N',ia,nchunk_2,nmax_2-nmin_2+1,-qq(jv,iv,is),&
                                          fnla(isa0+na_buff(1,is,buff,parai%cp_inter_me +1),jv,nmin_2),&
                                          ions1%nat*maxsys%nhxs, &
                                          gam(nmin_2,nmin_2),nstate,1.0_real_8,&
                                          dai_ptr(iac,nmin_2),ld_dai)
                                  END IF
                               END IF
                            END DO
                         END DO
                         !$omp end do nowait
                         iac0=iac0+ia*nlps_com%ngh(is)
                      END IF
                   END IF
                   isa0=isa0+ions0%na(is)
                END DO
                !$omp end parallel
             END IF
             !$ locks(buff,1) = .FALSE.
             !$omp flush(locks)
             start_dai=start_dai+ld_dai*nstate
          END DO
          !$omp parallel private (iac,iac0,isa0,is,iv,jv,isa,ia,ig,t1,t2,t3,ci,start_eiscr,&
          !$omp ld_dai,buff,i,iac_start,iac_end,isa_start)
          start_eiscr=0
          DO buff=1,buffcount
             ld_dai=map_isa(ions1%nat+1,buff,parai%cp_inter_me +1)
             IF (ld_dai.GT.0) THEN
                iac0=start_eiscr
                isa0=0
                DO is=1,ions1%nsp
                   IF (pslo_com%tvan(is)) THEN
                      isa_start=isa0+na_buff(1,is,buff,parai%cp_inter_me +1)-1
                      ia=na_buff(2,is,buff,parai%cp_inter_me +1)-&
                           na_buff(1,is,buff,parai%cp_inter_me +1)+1
                      IF(ia.GT.0)THEN
                         DO iv=1,nlps_com%ngh(is)
                            ci=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
                            iac_start=iac0+1
                            iac_end=iac0+ia
                            !MAKE USE OF THE SPECIAL STRUCTURE OF CI
                            IF (ABS(REAL(ci)).GT.0.5_real_8) THEN
                               !CI IS REAL
                               !$omp do
                               DO iac=iac_start,iac_end
                                  isa=isa_start+iac-iac_start+1
                                  DO ig=1,ncpw%ngw
                                     t1=REAL(eigr(ig,isa,1))
                                     t2=AIMAG(eigr(ig,isa,1))
                                     t3=twnl(ig,iv,is,1)*REAL(ci)
                                     eiscr(ig,iac)=CMPLX(t1*t3,t2*t3,kind=real_8)
                                  END DO
                               END DO
                               !$omp end do nowait                     
                            ELSE
                               !CI IS IMAGINARY
                               !$omp do
                               DO iac=iac_start,iac_end
                                  isa=isa_start+iac-iac_start+1
                                  DO ig=1,ncpw%ngw
                                     t1=REAL(eigr(ig,isa,1))
                                     t2=AIMAG(eigr(ig,isa,1))
                                     t3=twnl(ig,iv,is,1)*AIMAG(ci)
                                     eiscr(ig,iac)=CMPLX(-t2*t3,t1*t3,kind=real_8)
                                  END DO
                               END DO
                               !$omp end do nowait
                            END IF
                            iac0=iac0+ia
                         END DO
                      END IF
                   END IF
                   isa0=isa0+ions0%na(is)
                END DO
                start_eiscr=start_eiscr+ld_dai
             END IF
          END DO
          !$omp end parallel
       END IF

       IF (methread.EQ.0.OR.nthreads.EQ.1) THEN
          start_dai=0
          DO buff=1,buffcount
             ld_dai=map_isa(ions1%nat+1,buff,parai%cp_inter_me +1)
             IF (ld_dai.GT.0) THEN
                DO I=0,parai%NPROC-1
                   displs(i+1)=ld_dai*(paraw%nwa12(1,i)-1)
                   revcnt(i+1)=ld_dai*(paraw%nwa12(2,i)-paraw%nwa12(1,i)+1)
                ENDDO
             END IF
             !$omp flush(locks)
             !$ DO WHILE ( locks(buff,1) )
             !$omp flush(locks)
             !$ END DO
             !$ locks(buff,1) = .TRUE.
             !$omp flush(locks)
             IF (ld_dai.GT.0) THEN
                CALL MY_CONCATV_INPLACE(dai(start_dai+1,1),revcnt,displs,parai%allgrp)
             END IF
             !$ locks(buff,1)=.FALSE.
             !$omp flush(locks)
             !$ locks(buff,2)=.FALSE.
             !$omp flush(locks)           
             start_dai=start_dai+ld_dai*nstate
          END DO
          ready=.TRUE.
          !$omp flush(ready)
       END IF

       IF (methread.EQ.1.OR.nthreads.EQ.1) THEN
#ifdef _INTEL_MKL
          !$ CALL mkl_set_dynamic(0)
#endif
          start_dai=0
          start_eiscr=0
          DO buff=1,buffcount
             !$omp flush(locks)
             !$ DO WHILE ( locks(buff,2) )
             !$omp flush(locks)
             !$ END DO               
             !$omp flush(locks)
             !$ DO WHILE ( locks(buff,1) )
             !$omp flush(locks)
             !$ END DO
             IF (buff.LT.last_buff) THEN
                !$omp flush (ready)
                IF (ready) THEN
                   last_buff=buff         
                   CYCLE
                END IF
                ld_dai=map_isa(ions1%nat+1,buff,parai%cp_inter_me +1)
                IF (ld_dai.GT.0) THEN
                   CALL DGEMM('N','N',2*ncpw%ngw,nstate,ld_dai &
                        ,1._real_8,eiscr(1,start_eiscr+1),2*ncpw%ngw&
                        ,dai(start_dai+1,1),ld_dai,1._real_8,auxc(1,1),2*ncpw%ngw)
                   start_dai=start_dai+ld_dai*nstate
                   start_eiscr=start_eiscr+ld_dai
                END IF
             END IF
          END DO
#ifdef _INTEL_MKL
          !$ CALL mkl_set_dynamic(1)
#endif
       END IF
       !$OMP barrier

       !$ IF (nested_threads.NE.parai%ncpus) THEN
       !$    IF (methread.EQ.1.OR.nthreads.EQ.1) THEN
       !$       CALL omp_set_num_threads(parai%ncpus)
       !$       CALL omp_set_nested(0)
       !$    END IF
       !$ END IF

       !$omp end parallel

       start_dai=0
       start_eiscr=0
       DO buff=1,buffcount        
          ld_dai=map_isa(ions1%nat+1,buff,parai%cp_inter_me +1)
          IF (ld_dai.GT.0) THEN
             IF (buff.GE.last_buff) THEN
                CALL DGEMM('N','N',2*ncpw%ngw,nstate,ld_dai &
                     ,1._real_8,eiscr(1,start_eiscr+1),2*ncpw%ngw&
                     ,dai(start_dai+1,1),ld_dai,1._real_8,auxc(1,1),2*ncpw%ngw)
             END IF
             start_dai=start_dai+ld_dai*nstate
             start_eiscr=start_eiscr+ld_dai
          END IF
       END DO

       IF (parai%cp_nogrp .gt. 1) THEN
          CALL tiset(procedureN//'_grps_b',isub1)
          IF (redst) THEN
             CALL cp_grp_redist(auxc,nkpt%ngwk,nstate)
          END IF
          !$omp parallel do private(i,ig)
          DO i=1,nstate
             DO ig=1,ncpw%ngw
                c2(ig,i)=c2(ig,i)+auxc(ig,i)
             END DO
          END DO
#ifdef _USE_SCRATCHLIBRARY
          CALL free_scratch(il_auxc,auxc,procedureN//'_auxc')          
#else
          DEALLOCATE(auxc, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate auxc',& 
               __LINE__,__FILE__)
#endif
          CALL tihalt(procedureN//'_grps_b',isub1)
       END IF
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_eiscr,eiscr,procedureN//'_eiscr')
       CALL free_scratch(il_dai,dai,procedureN//'_dai')
#else
       DEALLOCATE(dai, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate dai',& 
            __LINE__,__FILE__)
       DEALLOCATE(eiscr, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate eiscr',& 
            __LINE__,__FILE__)
#endif
    ENDIF
    need_scratch=.FALSE.
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is) )CYCLE
       need_scratch=.TRUE.
    END DO
    
    IF (need_scratch) THEN

       il_dai(1)=maxsys%nax
       il_dai(2)=nstate
       il_auxc(1)=nkpt%ngwk
       il_auxc(2)=nstate
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_dai,dai,procedureN//'_dai')
       CALL request_scratch(il_auxc,auxc,procedureN//'_auxc')
#else
       ALLOCATE(dai(il_dai(1),il_dai(2)), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate dai',& 
            __LINE__,__FILE__)
       ALLOCATE(auxc(il_auxc(1),il_auxc(2)), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate auxc',& 
            __LINE__,__FILE__)
#endif
    END IF

    isa0=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN 
          ! see newcode above
       ELSEIF (sgpp1%tsgp(is)) THEN
          ! ==--------------------------------------------------------==
          ! ==  Stefan Goedecker PP                                   ==
          ! ==--------------------------------------------------------==
          DO iv=1,nlps_com%ngh(is)
             CALL zeroing(auxc)!,ngw*nstate)
             CALL zeroing(dai)!,maxsys%nax*nstate)
             DO i=paraw%nwa12(1,parai%mepos),paraw%nwa12(2,parai%mepos)
                ffi=f(i)
                IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
                DO ia=1,ions0%na(is)
                   isa=isa0+ia
                   l=nghtol(iv,is)+1
                   ki=sgpp2%lfval(iv,is)
                   li=sgpp2%lpval(iv,is)
                   dd=0.0_real_8
                   DO jv=1,nlps_com%ngh(is)
                      l2=nghtol(jv,is)+1
                      lj=sgpp2%lpval(jv,is)
                      IF (l.EQ.l2.AND.li.EQ.lj) THEN
                         kj=sgpp2%lfval(jv,is)
                         dai(ia,i)=dai(ia,i)+&
                              fnla(isa,jv,i)*sgpp2%hlsg(ki,kj,l,is)
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             CALL mp_sum(dai,nstate*maxsys%nax,parai%allgrp)
             IF (ncpw%ngw.GT.0)&
                  CALL dgemm("N","N",2*ncpw%ngw,nstate,ions0%na(is),1._real_8,&
                  eigr(1,isa0+1,1),2*ncpw%ngw,dai(1,1),maxsys%nax,0._real_8,auxc,2*ncpw%ngw)
             ct=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
             !$omp parallel DO private(IG,I,FFI,CTM)
             DO i=1,nstate
                ffi=f(i)
                IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
                ctm=-ffi*ct
                DO ig=1,ncpw%ngw
                   auxc(ig,i)=ctm*auxc(ig,i)
                   c2(ig,i)=c2(ig,i)+twnl(ig,iv,is,1)*auxc(ig,i)
                ENDDO
             ENDDO
          ENDDO
       ELSE
          ! ==--------------------------------------------------------==
          ! ==  BACHELET HAMANN SCHLUTER                              ==
          ! ==--------------------------------------------------------==
          DO iv=1,nlps_com%ngh(is)
             CALL zeroing(dai)!,maxsys%nax*nstate)
             DO i=paraw%nwa12(1,parai%mepos),paraw%nwa12(2,parai%mepos)
                CALL dcopy(ions0%na(is),fnla(isa0+1,iv,i),1,dai(1,i),1)
             ENDDO
             CALL mp_sum(dai,nstate*maxsys%nax,parai%allgrp)
             IF (ncpw%ngw.GT.0)&
                  CALL dgemm('N','N',2*ncpw%ngw,nstate,ions0%na(is),1._real_8,&
                  eigr(1,isa0+1,1),2*ncpw%ngw,dai,maxsys%nax,0.0_real_8,auxc,2*ncpw%ngw)
             ct=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)*wsg(is,iv)
             !$omp parallel DO private(IG,I,FFI,CTM)
             DO i=1,nstate
                ffi=f(i)
                IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
                ctm=-ffi*ct
                DO ig=1,ncpw%ngw
                   auxc(ig,i)=ctm*auxc(ig,i)
                   c2(ig,i)=c2(ig,i)+twnl(ig,iv,is,1)*auxc(ig,i)
                ENDDO
             ENDDO
          ENDDO
          ! ==--------------------------------------------------------------==
       ENDIF
       isa0=isa0+ions0%na(is)
    ENDDO
    IF (need_scratch) THEN
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_auxc,auxc,procedureN//'_auxc')
       CALL free_scratch(il_dai,dai,procedureN//'_dai')
#else
       DEALLOCATE(dai, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate dai',& 
            __LINE__,__FILE__)
       DEALLOCATE(auxc, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate auxc',& 
            __LINE__,__FILE__)
#endif
    END IF

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nlforce
  ! ==================================================================
  SUBROUTINE give_scr_nlforce(il_gam,il_auxc,il_ddia,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: il_gam, il_auxc, il_ddia, &
                                                nstate

! ==--------------------------------------------------------------==

    il_gam =imagp*nstate*nstate
    il_auxc=2*nkpt%ngwk*nstate+2
    il_ddia=2*maxsys%nax*nstate+2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_nlforce
  ! ==================================================================


END MODULE nlforce_utils
