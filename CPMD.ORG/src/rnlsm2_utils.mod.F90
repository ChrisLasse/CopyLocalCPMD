#include "cpmd_global.h"

MODULE rnlsm2_utils
  USE cppt,                            ONLY: gk,&
                                             twnl
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms  
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE fnl_utils,                       ONLY: unpack_dfnl
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: eigkr,&
                                             rk
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_walltime
  USE mm_dimmod,                       ONLY: mmdim
  USE mp_interface,                    ONLY: mp_sum,mp_bcast
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlm,&
                                             nlps_com
  !$ USE omp_lib,                      ONLY: omp_set_nested,&
  !$                                         omp_get_thread_num,&
  !$                                         omp_set_num_threads
  USE parac,                           ONLY: parai,paral
  USE sfac,                            ONLY: dfnl,&
                                             dfnla,&
                                             dfnl_packed,&
                                             eigr,&
                                             il_dfnl_packed
  USE system,                          ONLY: ncpw,&
                                             nkpt,&
                                             cnti,&
                                             cntr,&
                                             parap,&
                                             parm,iatpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch,&
                                             save_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlsm2
  PUBLIC :: give_scr_rnlsm2

CONTAINS

  ! ==================================================================
  SUBROUTINE RNLSM2(c0,nstate,ikpt,ikind,store_dfnl_packed)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ARRAY DFNL WHICH IS USED IN THE SUBROUTINE RNLFOR       ==
    ! ==  K-POINT VERSION IS IMPLEMENTED                              ==
    ! ==          NOT IMPLEMENTED FOR TSHEL(IS)=TRUE                  ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate)
    INTEGER                                  :: ikpt, ikind
    LOGICAL,OPTIONAL,INTENT(IN)              :: store_dfnl_packed

    LOGICAL                                  :: store,unpack
    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlsm2'
    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8) ,&
                                                zzero = (0._real_8,0._real_8)
    COMPLEX(real_8)                          :: ci, zfac
    REAL(real_8)                             :: cii, cir, ei, er, tfac, temp
    INTEGER                                  :: i, ia, ierr, ig, ii, is, iv, k, isa, isa0, &
                                                isub, isub2, offset, offset0,  start_state,&
                                                methread, nthreads, nested_threads, &
                                                nst_buffer(2,15), nst, tot_work,&
                                                im, na_grp(2,ions1%nsp,parai%cp_nogrp), &
                                                il_gktemp(2),il_eiscr(2),&
                                                il_t(1), ia_sum,start_ia,end_ia,&
                                                 buffcount, buff
#ifdef _USE_SCRATCHLIBRARY    
    REAL(real_8), POINTER __CONTIGUOUS       :: t(:),gktemp(:,:)
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: eiscr(:,:)
#else
    REAL(real_8),ALLOCATABLE                 :: t(:),gktemp(:,:)
    COMPLEX(real_8),ALLOCATABLE              :: eiscr(:,:)
#endif
    REAL(real_8), SAVE                       :: timings(2)=0.0_real_8
    INTEGER, SAVE                            :: autotune_it=0
    !$ LOGICAL                               :: locks(cnti%rnlsm2_bc,2)
    ! ==--------------------------------------------------------------==
    CALL TISET(procedureN,ISUB)
    ! ==--------------------------------------------------------------==
    ! IF no non-local components -> return.
    IF(NLM.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    ! If we store packed dfnl, we are able to skip glosum
    IF(PRESENT(store_dfnl_packed))THEN
       store=store_dfnl_packed
    ELSE
       store=.FALSE.
    END IF

    ! split atoms between cp_groups
    CALL cp_grp_split_atoms(na_grp)

    IF(tkpts%tkpnt) THEN
       zfac=parm%tpiba*zone
    ELSE
       tfac=2._real_8*parm%tpiba
    ENDIF
  
    !divede work between buffers along states   
    nst_buffer(1,:)=0
    nst_buffer(2,:)=-1
    nst_buffer(1,1)=1
    IF(store)THEN
       !disable overlapping when no glosum is requiered
       buffcount=1
       nst_buffer(2,1)=nstate
    ELSE
       IF (autotune_it.EQ.0.AND.cnti%rnlsm_autotune_maxit.GT.0) THEN
          timings=0.0_real_8
          cnti%rnlsm2_bc=1
          cntr%rnlsm2_b1=1.0_real_8
          cntr%rnlsm2_b2=0.0_real_8
       END IF
       buffcount=cnti%rnlsm2_bc       
       nst_buffer(2,1)=CEILING(cntr%rnlsm2_b1*REAL(nstate))
       DO i=2,buffcount
          nst_buffer(1,i)=nst_buffer(2,i-1)+1
          nst_buffer(2,i)=nst_buffer(1,i)+CEILING(cntr%rnlsm2_b2*REAL(nstate))
          IF(i.EQ.buffcount) nst_buffer(2,i)=nstate
       END DO
    END IF
    
    temp=0
    !sanity check
    DO i=1,buffcount
       temp=temp+nst_buffer(2,i)-nst_buffer(1,i)+1
    END DO
    
    IF(temp.GT.nstate) THEN
       nst_buffer(2,buffcount)=nst_buffer(2,buffcount)-(temp-nstate)
    ELSE IF(temp.LT.nstate) THEN
       nst_buffer(2,buffcount)=nst_buffer(2,buffcount)+(nstate-temp)
    END IF

    tot_work=0
    DO is=1,ions1%nsp
       tot_work=tot_work+&
            (na_grp(2,is,parai%cp_inter_me +1)-na_grp(1,is,parai%cp_inter_me +1)+1)*&
            nlps_com%ngh(is)
    END DO
    tot_work=tot_work*3

    !use scratch
    il_dfnl_packed(1)=tot_work*imagp
    il_dfnl_packed(2)=nstate
    il_eiscr(1)=nkpt%ngwk
    il_eiscr(2)=tot_work
    il_gktemp(1)=nkpt%ngwk
    il_gktemp(2)=3
    il_t(1)=nkpt%ngwk    
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_dfnl_packed,dfnl_packed,'DFNL_packed')
    CALL request_scratch(il_gktemp,gktemp,procedureN//'_gktemp')
    CALL request_scratch(il_t,t,procedureN//'_t')
    CALL request_scratch(il_eiscr,eiscr,procedureN//'_eiscr')
#else
    ALLOCATE(dfnl_packed(il_dfnl_packed(1),il_dfnl_packed(2)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate dfnl_packed',& 
         __LINE__,__FILE__)
    ALLOCATE(eiscr(il_eiscr(1),il_eiscr(2)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate eiscr',& 
         __LINE__,__FILE__)
    ALLOCATE(gktemp(il_gktemp(1),il_gktemp(2)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate gktemp',& 
         __LINE__,__FILE__)
    ALLOCATE(t(il_t(1)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate t',& 
         __LINE__,__FILE__)
#endif
    !$omp parallel private(ig,k,offset0,isa0,is,iv,ci,cir,cii,ia,offset,isa,er,ei,ia_sum,&
    !$omp start_ia,end_ia)
    IF (tkpts%tkpnt) THEN
       !$OMP do collapse(2)
       DO ig=1,nkpt%ngwk
          DO k=1,3
             gktemp(ig,k)=gk(k,ig)+RK(K,IKIND)
          ENDDO
       ENDDO
       !$omp end do nowait
    ELSE
       !$OMP do collapse(2)
       DO ig=1,nkpt%ngwk
          DO k=1,3
             gktemp(ig,k)=gk(k,ig)
          ENDDO
       ENDDO
       !$omp end do nowait
    ENDIF

    !offset is the counter
    offset0=0
    isa0=0
    DO is=1,ions1%nsp
       DO iv=1,nlps_com%ngh(is)
          ci=(0.0_real_8,-1.0_real_8)**(nghtol(iv,is)+1)
          cir=REAL(ci,kind=real_8)
          cii=AIMAG(ci)
          !$omp barrier
          !Make use of the special structure of CI
          IF (ABS(cir).GT.0.5_real_8) THEN
             !CI is real
             !$omp do 
             DO ig=1,nkpt%ngwk
                t(ig)=twnl(ig,iv,is,ikind)*cir
             END DO
             start_ia=na_grp(1,is,parai%cp_inter_me +1)
             end_ia=na_grp(2,is,parai%cp_inter_me +1)
             ia_sum=end_ia-start_ia+1            
             DO k=1,3
                !$OMP do 
                DO ia=start_ia,end_ia
                   isa=isa0+ia
                   offset=offset0+(k-1)*ia_sum+ia-start_ia+1
                   DO ig=1,nkpt%ngwk
                      er=REAL(eigkr(ig,isa,ikind),kind=real_8)
                      ei=AIMAG(eigkr(ig,isa,ikind))
                      eiscr(ig,offset) = CMPLX(t(ig)*er*gktemp(ig,k),t(ig)*ei*&
                           gktemp(ig,k),kind=real_8)
                   ENDDO
                   IF(geq0) THEN 
                      IF(tkpts%tkpnt) THEN
                         eiscr(ncpw%ngw+1,offset)=&
                              cmplx(0._real_8,0._real_8,kind=real_8)
                      ELSE
                         eiscr(1,offset)=&
                              0.5_real_8*eiscr(1,offset)
                      ENDIF
                   ENDIF
                ENDDO
                !$omp end do nowait
             END DO
          ELSE
             !CI is imaginary
             !$omp do 
             DO ig=1,nkpt%ngwk
                t(ig)=twnl(ig,iv,is,ikind)*cii
             END DO
             start_ia=na_grp(1,is,parai%cp_inter_me +1)
             end_ia=na_grp(2,is,parai%cp_inter_me +1)
             ia_sum=end_ia-start_ia+1            
             DO k=1,3
                !$OMP do 
                DO ia=start_ia,end_ia
                   isa=isa0+ia
                   offset=offset0+(k-1)*ia_sum+ia-start_ia+1
                   DO ig=1,nkpt%ngwk
                      er=real(eigkr(ig,isa,ikind),kind=real_8)
                      ei=aimag(eigkr(ig,isa,ikind))
                      eiscr(ig,offset) =CMPLX(-t(ig)*ei*gktemp(ig,k),t(ig)*er*&
                           gktemp(ig,k),kind=real_8)
                   ENDDO
                   IF(geq0) THEN 
                      IF(tkpts%tkpnt) THEN
                         eiscr(ncpw%ngw+1,offset)=cmplx(0._real_8,0._real_8)
                      ELSE
                         eiscr(1,offset)=0.5_real_8*eiscr(1,offset)
                      ENDIF
                   ENDIF
                ENDDO
                !$omp end do nowait
             END DO
          ENDIF
          offset0=offset0+ia_sum*3
       ENDDO
       isa0=isa0+ions0%na(is)
    ENDDO
    !$omp end parallel

    buff=1
    IF (.NOT.store)THEN
       IF(autotune_it.LT.cnti%rnlsm_autotune_maxit) THEN
          autotune_it=autotune_it+1
          temp=m_walltime()
       END IF
    END IF
    nst=nst_buffer(2,buff)-nst_buffer(1,buff)+1
    IF(TKPTS%TKPNT)THEN
       CALL ZGEMM('C','N',tot_work,nst,nkpt%ngwk,zfac &
            ,eiscr(1,1),nkpt%ngwk &
            ,c0(1,1),nkpt%ngwk,zzero &
            ,dfnl_packed(1,1),tot_work)
    ELSE
       CALL DGEMM('T','N',tot_work,nst,2*nkpt%ngwk,tfac &
            ,eiscr(1,1),2*nkpt%ngwk &
            ,c0(1,1),2*nkpt%ngwk,0._real_8 &
            ,dfnl_packed(1,1),tot_work)
    END IF
    !finished here for "store" option
    IF(.NOT.store)THEN
       IF(autotune_it.GT.0)timings(1)=timings(1)+m_walltime()-temp

#ifdef _USE_OVERLAPPING_COM_COMP
       nthreads=MIN(2,parai%ncpus)
       nested_threads=(MAX(parai%ncpus-1,1))
#else
       nthreads=1
       nested_threads=parai%ncpus
#endif
       IF(buffcount.eq.1)THEN
          nthreads=1
          nested_threads=parai%ncpus
       END IF
       methread=0
       !$ locks=.TRUE.
       !$ locks(1,1)=.FALSE.
       !$omp parallel if(nthreads.eq.2) num_threads(nthreads) &
       !$omp private(methread,buff,nst,start_state)&
       !$omp shared(locks,nthreads) proc_bind(close)
       !$ methread = omp_get_thread_num()
       !$ IF (nested_threads .NE. parai%ncpus) THEN
       !$    IF (methread .EQ. 1 .OR. nthreads .EQ. 1) THEN
       !$       CALL omp_set_num_threads(nested_threads)
#ifdef _INTEL_MKL
       !$       CALL mkl_set_dynamic(0)
#endif
       !$       CALL omp_set_nested(1)
       !$    END IF
       !$ END IF
       !$OMP barrier

       IF(methread.EQ.1.OR.nthreads.EQ.1) THEN
          DO buff=2,buffcount
             start_state=nst_buffer(1,buff)
             nst=nst_buffer(2,buff)-nst_buffer(1,buff)+1
             IF(tkpts%tkpnt) THEN
                CALL ZGEMM('C','N',tot_work,nst,nkpt%ngwk,zfac &
                     ,eiscr(1,1),2*nkpt%ngwk &
                     ,c0(1,start_state),nkpt%ngwk,zzero &
                     ,dfnl_packed(1,start_state),tot_work)
             ELSE
                CALL DGEMM('T','N',tot_work,nst,2*nkpt%ngwk,tfac &
                     ,eiscr(1,1),2*nkpt%ngwk &
                     ,c0(1,start_state),2*nkpt%ngwk,0._real_8 &
                     ,dfnl_packed(1,start_state),tot_work)
             ENDIF
             !$ locks(buff,1) = .FALSE.
             !$omp flush(locks)
          ENDDO
       ENDIF

       IF(methread.EQ.0.or.nthreads.EQ.1) THEN
          DO buff=1,buffcount
             nst=nst_buffer(2,buff)-nst_buffer(1,buff)+1          
             start_state=nst_buffer(1,buff)
             CALL TISET(procedureN//'_barrier',ISUB2)
             !$omp flush(locks)
             !$ DO WHILE ( locks(buff,1) )
             !$omp flush(locks)
             !$ END DO
             !$ locks(buff,1) = .TRUE.
             !$omp flush(locks)
             CALL TIHALT(procedureN//'_barrier',ISUB2)
             IF(autotune_it.GT.0) temp=m_walltime()
             CALL mp_sum(dfnl_packed(:,start_state),nst*tot_work,parai%allgrp)
             IF(autotune_it.GT.0) timings(2)=timings(2)+m_walltime()-temp
             !$ locks(buff,1)=.FALSE.
             !$omp flush(locks)
             !$ locks(buff,2)=.FALSE.
             !$omp flush(locks)
          ENDDO
       ENDIF

       !$omp barrier
       !$ IF (nested_threads .NE. parai%ncpus) THEN
       !$    IF (methread .EQ. 1 .OR. nthreads .EQ. 1) THEN
       !$       CALL omp_set_num_threads(parai%ncpus)
#ifdef _INTEL_MKL
       !$       CALL mkl_set_dynamic(1)
#endif
       !$       CALL omp_set_nested(0)
       !$    END IF
       !$ END IF

       !$OMP end parallel

       CALL unpack_dfnl(dfnl_packed,dfnl,dfnla,na_grp,tkpts%tkpnt,ikind)

       !set buffcount/buffsizes
       IF(cnti%rnlsm_autotune_maxit.GT.0)THEN
          IF(autotune_it.EQ.cnti%rnlsm_autotune_maxit)THEN
             autotune_it=autotune_it+1
             IF(paral%parent)THEN
                !first buffersize:
                IF(timings(1).GT.timings(2))THEN
                   !total work = nstate
                   tot_work=nstate
                   !work for overlapping part
                   tot_work=CEILING(REAL(tot_work,kind=real_8)*timings(2)/timings(1))
                   !take at least 100 in each buffer
                   cnti%rnlsm2_bc=tot_work/100
                   IF(cnti%rnlsm2_bc.EQ.0) cnti%rnlsm2_bc=2
                   IF(cnti%rnlsm2_bc.GT.15) cnti%rnlsm2_bc=15
                   cntr%rnlsm2_b1=1.0_real_8-timings(2)/timings(1)
                   IF(cnti%rnlsm2_bc.GT.1) THEN
                      cntr%rnlsm2_b2=timings(2)/timings(1)/real(cnti%rnlsm2_bc-1,kind=real_8)
                   ELSE
                      cntr%rnlsm2_b2=0.0_real_8
                   END IF
                ELSE
                   tot_work=nstate
                   cnti%rnlsm2_bc=tot_work/100
                   IF(cnti%rnlsm2_bc.GT.15) cnti%rnlsm2_bc=15
                   cntr%rnlsm2_b1=1.0_real_8/real(cnti%rnlsm2_bc)
                   cntr%rnlsm2_b2=1.0_real_8/real(cnti%rnlsm2_bc)
                END IF
                WRITE(6,*) "####rnlsm2_autotuning results####"
                WRITE(6,*) cnti%rnlsm2_bc
                WRITE(6,*) cntr%rnlsm2_b1
                WRITE(6,*) cntr%rnlsm2_b2
                WRITE(6,*) timings
             END IF
             CALL mp_bcast(cnti%rnlsm2_bc,parai%source,parai%allgrp)
             CALL mp_bcast(cntr%rnlsm2_b1,parai%source,parai%allgrp)
             CALL mp_bcast(cntr%rnlsm2_b2,parai%source,parai%allgrp)
          END IF
       END IF
    END IF

#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_eiscr,eiscr,procedureN//'_eiscr')
    CALL free_scratch(il_t,t,procedureN//'_t')
    CALL free_scratch(il_gktemp,gktemp,procedureN//'_gktemp')
    IF(store)THEN
       CALL save_scratch(il_dfnl_packed,dfnl_packed,'DFNL_packed')
    ELSE
       CALL free_scratch(il_dfnl_packed,dfnl_packed,'DFNL_packed')
       il_dfnl_packed(1)=0
    END IF
#else
    DEALLOCATE(eiscr, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate eiscr',& 
         __LINE__,__FILE__)
    DEALLOCATE(gktemp, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate gktemp',& 
         __LINE__,__FILE__)
    DEALLOCATE(t, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate t',& 
         __LINE__,__FILE__)
    
    IF(.NOT.store)THEN
       DEALLOCATE(dfnl_packed, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate dfnl_packed',& 
            __LINE__,__FILE__)
       il_dfnl_packed(1)=0
    END IF
#endif       
    CALL TIHALT(procedureN,ISUB)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE RNLSM2
  ! ==================================================================
  SUBROUTINE give_scr_rnlsm2(lrnlsm2,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrnlsm2
    CHARACTER(LEN=30)                        :: tag
    INTEGER                                  :: nstate,il_dfnl_packed,il_eiscr(2),&
                                                il_gktemp(2),tot_work,is, il_t,&
                                                na_grp(2,ions1%nsp,parai%cp_nogrp)
    
    lrnlsm2=0
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rnlsm2
  ! ==================================================================

END MODULE rnlsm2_utils
