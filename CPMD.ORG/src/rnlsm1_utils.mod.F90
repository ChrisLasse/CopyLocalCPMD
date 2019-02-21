#include "cpmd_global.h"

MODULE rnlsm1_utils
  USE cppt,                            ONLY: twnl
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms
  USE error_handling,                  ONLY: stopgm
  USE fnl_utils,                       ONLY: unpack_fnl
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: eigkr
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_walltime
  USE mm_dimmod,                       ONLY: mmdim
  USE mp_interface,                    ONLY: mp_sum,mp_bcast
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlm,&
                                             nlps_com
  USE nvtx_utils
  !$ USE omp_lib,                      ONLY: omp_set_nested,&
  !$                                         omp_get_thread_num,&
  !$                                         omp_set_num_threads
  USE parac,                           ONLY: parai,paral
  USE sfac,                            ONLY: fnl,&
                                             fnla,&
                                             fnl_packed,&
                                             il_fnl_packed
  USE system,                          ONLY: cntl,&
                                             cnti,&
                                             cntr,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,iatpt,parap,cnti
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

  PUBLIC :: rnlsm1
  PUBLIC :: give_scr_rnlsm1

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlsm1(c0,nstate,ikind,store_fnl_packed)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ARRAY FNL (ALSO CALLED BEC IN SOME VANDERBILT ROUTINES) ==
    ! ==  K-POINT VERSION IS IMPLEMENTED                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate)
    INTEGER                                  :: ikind
    LOGICAL,OPTIONAL,INTENT(IN)              :: store_fnl_packed
    LOGICAL                                  :: store
    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlsm1'
    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8) ,&
                                                zzero = (0._real_8,0._real_8)
    COMPLEX(real_8)                          :: ci, zfac
    REAL(real_8)                             :: cii, cir, ei, er, tfac, temp
    INTEGER                                  :: i, ia, iv, k, ig, is, isa, isa0, im,&
                                                isub, isub2, ierr, buffcount, buff,&
                                                offset, offset0, start_state, methread,&
                                                nthreads, nested_threads, nst_buffer(2,15),&
                                                nst, tot_work,&
                                                na_grp(2,ions1%nsp,parai%cp_nogrp),&
                                                il_eiscr(2),il_t(1)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: t(:)
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: eiscr(:,:)
#else
    REAL(real_8),ALLOCATABLE                 :: t(:)
    COMPLEX(real_8),ALLOCATABLE              :: eiscr(:,:)
#endif
    REAL(real_8), SAVE                       :: timings(2)=0.0_real_8
    INTEGER, SAVE                            :: autotune_it=0
    !$ LOGICAL                               :: locks(cnti%rnlsm1_bc,2)
   

    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )

    ! ==--------------------------------------------------------------==
    ! IF no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    IF (present(store_fnl_packed) ) THEN
       store=store_fnl_packed
    ELSE
       store=.FALSE.
    END IF

    ! split states between cp groups
    CALL cp_grp_split_atoms(na_grp)

    IF (autotune_it.EQ.0.AND.cnti%rnlsm_autotune_maxit.GT.0) THEN
       timings=0.0_real_8
       cnti%rnlsm1_bc=1
       cntr%rnlsm1_b1=1.0_real_8
       cntr%rnlsm1_b2=0.0_real_8
    END IF  

    buffcount=cnti%rnlsm1_bc       

    !divede work between buffers along states
    nst_buffer(1,:)=0
    nst_buffer(2,:)=-1
    nst_buffer(1,1)=1
    nst_buffer(2,1)=CEILING(cntr%rnlsm1_b1*REAL(nstate))
    DO i=2,buffcount
       nst_buffer(1,i)=nst_buffer(2,i-1)+1
       nst_buffer(2,i)=nst_buffer(1,i)+CEILING(cntr%rnlsm1_b2*REAL(nstate))
       IF(i.EQ.buffcount) nst_buffer(2,i)=nstate
    END DO
    
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
       tot_work=tot_work+(na_grp(2,is,parai%cp_inter_me +1)-na_grp(1,is,parai%cp_inter_me +1)+1)*nlps_com%ngh(is)
    END DO

    !50 iterations for autotuning to collect averaged timings
    IF (autotune_it.LT.cnti%rnlsm_autotune_maxit) autotune_it=autotune_it+1
    IF(nstate.lt.100)THEN
       buffcount=1
       nst_buffer(1,:)=0
       nst_buffer(2,:)=-1
       nst_buffer(1,1)=1
       nst_buffer(2,1)=nstate
    END IF
    !use scratch

    il_fnl_packed(1)=tot_work*imagp
    il_fnl_packed(2)=nstate
    il_eiscr(1)=nkpt%ngwk
    il_eiscr(2)=tot_work*imagp
    il_t(1)=nkpt%ngwk   
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_fnl_packed,fnl_packed,'fnl_packed')
    CALL request_scratch(il_t,t,procedureN//'_t')
    CALL request_scratch(il_eiscr,eiscr,procedureN//'_eiscr')
#else
    ALLOCATE(fnl_packed(il_fnl_packed(1),il_fnl_packed(2)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnl_packed',& 
         __LINE__,__FILE__)
    ALLOCATE(eiscr(il_eiscr(1),il_eiscr(2)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate eiscr',& 
         __LINE__,__FILE__)
    ALLOCATE(t(il_t(1)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate t',& 
         __LINE__,__FILE__)
#endif
    !$omp parallel private(offset0,isa0,is,iv,ci,cir,cii,ig,ia,offset,isa,er,ei)
    offset0=0
    isa0=0
    DO is=1,ions1%nsp
       DO iv=1,nlps_com%ngh(is)
          ci=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
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
             !$OMP do 
             DO ia=na_grp(1,is,parai%cp_inter_me +1),na_grp(2,is,parai%cp_inter_me +1)
                offset=offset0+ia-na_grp(1,is,parai%cp_inter_me +1)+1
                isa=isa0+ia
                DO ig=1,nkpt%ngwk
                   er=REAL(eigkr(ig,isa,ikind),kind=real_8)
                   ei=AIMAG(eigkr(ig,isa,ikind))
                   eiscr(ig,offset) = CMPLX(T(ig)*ER,T(ig)*EI,kind=real_8)
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
          ELSE
             !CI is imaginary
             !$omp do 
             DO ig=1,nkpt%ngwk
                t(ig)=twnl(ig,iv,is,ikind)*cii
             END DO
             !$OMP do
             DO ia=na_grp(1,is,parai%cp_inter_me +1),na_grp(2,is,parai%cp_inter_me +1)
                offset=offset0+ia-na_grp(1,is,parai%cp_inter_me +1)+1
                isa=isa0+ia
                DO ig=1,nkpt%ngwk
                   er=real(eigkr(ig,isa,ikind),kind=real_8)
                   ei=aimag(eigkr(ig,isa,ikind))
                   eiscr(ig,offset) =CMPLX(-t(ig)*ei,t(ig)*er,kind=real_8)
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
          ENDIF
          offset0=offset0+&
               na_grp(2,is,parai%cp_inter_me +1)-na_grp(1,is,parai%cp_inter_me +1)+1
       ENDDO
       isa0=isa0+ions0%na(is)
    ENDDO
    !$omp end parallel

    buff=1   
    IF (autotune_it.GT.0) temp=m_walltime()
    nst=nst_buffer(2,buff)-nst_buffer(1,buff)+1
    IF(TKPTS%TKPNT) THEN       
       CALL ZGEMM('C','N',tot_work,nst,nkpt%ngwk,zone &
            ,eiscr(1,1),nkpt%ngwk &
            ,c0(1,1),nkpt%ngwk,zzero &
            ,fnl_packed(1,1),tot_work)
    ELSE
       CALL DGEMM('T','N',tot_work,nst,2*nkpt%ngwk,2._real_8 &
            ,eiscr(1,1),2*nkpt%ngwk &
            ,c0(1,1),2*nkpt%ngwk,0._real_8 &
            ,fnl_packed(1,1),tot_work)
    ENDIF
    IF(autotune_it.GT.0)timings(1)=timings(1)+m_walltime()-temp

    !now we split up the threads, thread=0 is used to communicate,
    !thread=1 is used to do the dgemms and copy the data to its destionation.
    !dim 1 of buff contains locks for the dgemm and dim 2 contains locks for
    !the communication,since we already computed buff 1 and 3 we skip setting
    !the lock here


#ifdef _USE_OVERLAPPING_COM_COMP
    nthreads=MIN(2,parai%ncpus)
    nested_threads=(MAX(parai%ncpus-1,1))
#else
    nthreads=1
    nested_threads=parai%ncpus
#endif
    methread=0
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
             CALL ZGEMM('C','N',tot_work,nst,nkpt%ngwk,zone &
                  ,eiscr(1,1),2*nkpt%ngwk &
                  ,c0(1,start_state),nkpt%ngwk,zzero &
                  ,fnl_packed(1,start_state),tot_work)
          ELSE
             CALL DGEMM('T','N',tot_work,nst,2*nkpt%ngwk,2._real_8 &
                  ,eiscr(1,1),2*nkpt%ngwk &
                  ,c0(1,start_state),2*nkpt%ngwk,0._real_8 &
                  ,fnl_packed(1,start_state),tot_work)
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
          CALL mp_sum(fnl_packed(:,start_state),nst*tot_work,parai%allgrp)
          IF(autotune_it.GT.0) timings(2)=timings(2)+m_walltime()-temp
          !$ locks(buff,1)=.FALSE.
          !$omp flush(locks)
          !$ locks(buff,2)=.FALSE.
          !$omp flush(locks)
       ENDDO
    ENDIF
    
    IF (methread.EQ.1.or.nthreads.EQ.1) THEN
       DO buff=1,BUFFCOUNT
          !$omp flush(locks)
          !$ DO WHILE ( locks(buff,2) )
          !$omp flush(locks)
          !$ END DO               
          !$omp flush(locks)
          !$ DO WHILE ( locks(buff,1) )
          !$omp flush(locks)
          !$ END DO
          CALL unpack_fnl(fnl_packed,fnl,fnla,na_grp,nst_buffer(1,buff),nst_buffer(2,buff),&
               tkpts%tkpnt,ikind)
       END DO
    END IF
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

#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_eiscr,eiscr,procedureN//'_eiscr')
    CALL free_scratch(il_t,t,procedureN//'_t')
    IF(STORE)THEN
       CALL save_scratch(il_fnl_packed,fnl_packed,'fnl_packed')
    ELSE
       CALL free_scratch(il_fnl_packed,fnl_packed,'fnl_packed')
    END IF
#else
    IF(.NOT.store)THEN
       DEALLOCATE(fnl_packed, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fnl_packed',& 
            __LINE__,__FILE__)
    END IF
    DEALLOCATE(eiscr, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate eiscr',& 
         __LINE__,__FILE__)
    DEALLOCATE(t, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate t',& 
         __LINE__,__FILE__)
#endif
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
                cnti%rnlsm1_bc=tot_work/100
                IF(cnti%rnlsm1_bc.EQ.0) cnti%rnlsm1_bc=2
                IF(cnti%rnlsm1_bc.GT.15) cnti%rnlsm1_bc=15
                cntr%rnlsm1_b1=1.0_real_8-timings(2)/timings(1)
                IF(cnti%rnlsm1_bc.GT.1) THEN
                   cntr%rnlsm1_b2=timings(2)/timings(1)/real(cnti%rnlsm1_bc-1,kind=real_8)
                ELSE
                   cntr%rnlsm1_b2=0.0_real_8
                END IF
             ELSE
                tot_work=nstate
                cnti%rnlsm1_bc=tot_work/100
                IF(cnti%rnlsm1_bc.GT.15) cnti%rnlsm1_bc=15
                cntr%rnlsm1_b1=1.0_real_8/real(cnti%rnlsm1_bc)
                cntr%rnlsm1_b2=1.0_real_8/real(cnti%rnlsm1_bc)
             END IF
             WRITE(6,*) "####rnlsm1_autotuning results####"
             WRITE(6,*) cnti%rnlsm1_bc
             WRITE(6,*) cntr%rnlsm1_b1
             WRITE(6,*) cntr%rnlsm1_b2
             WRITE(6,*) timings
          END IF
          CALL mp_bcast(cnti%rnlsm1_bc,parai%source,parai%allgrp)
          CALL mp_bcast(cntr%rnlsm1_b1,parai%source,parai%allgrp)
          CALL mp_bcast(cntr%rnlsm1_b2,parai%source,parai%allgrp)
       END IF
    END IF

    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==

  END SUBROUTINE rnlsm1
  ! ==================================================================
  SUBROUTINE give_scr_rnlsm1(lrnlsm1,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrnlsm1
    CHARACTER(LEN=30)                        :: tag
    INTEGER                                  :: nstate,il_fnl_packed,il_eiscr(2),&
                                                tot_work,is, il_t,&
                                                na_grp(2,ions1%nsp,parai%cp_nogrp)
    
    lrnlsm1=0
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rnlsm1
  ! ==================================================================

END MODULE rnlsm1_utils
