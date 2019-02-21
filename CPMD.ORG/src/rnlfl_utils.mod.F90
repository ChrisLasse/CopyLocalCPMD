#include "cpmd_global.h"

MODULE rnlfl_utils
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms  
  USE cvan,                            ONLY: qq
  USE distribution_utils,              ONLY: dist_size
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum,&
                                             mp_win_alloc_shared_mem,&
                                             mp_win_sync,&
                                             mp_sync
  USE nlps,                            ONLY: imagp,&
                                             nlps_com
  !$ USE omp_lib,                       ONLY: omp_get_thread_num
  USE part_1d,                         ONLY: part_1d_get_blk_bounds
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: fnla,&
                                             dfnl_packed,&
                                             il_dfnl_packed
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             parap,&
                                             maxsys,&
                                             norbpe
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE, INTRINSIC :: ISO_C_BINDING,     ONLY: C_PTR,&
                                             C_F_POINTER
  
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch,&
                                             request_saved_scratch
#endif
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlfl

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlfl(fion,gamma,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE FORCE ON THE IONIC DEGREES OF FREEDOM DUE TO THE        ==
    ! ==  ORTHOGONALITY CONSTRAINED IN VANDERBILT PSEUDO-POTENTIALS   ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: GAMMA(nstate,*)
    INTEGER                                  :: nkpoint
    INTEGER                                  :: i, ia, ii, ik, is, isa, isa0, &
                                                iv, j, jv, k , ia_sum, iac, isub, &
                                                il_fnlt(3), il_fiont(4),&
                                                na_grp(2,ions1%nsp,parai%cp_nogrp),&
                                                ierr, tot_work, proc, &
                                                work_proc(2,0:parai%node_nproc-1),&
                                                start_ia, end_ia, offset1, offset2,&
                                                offset3, offset4, offset_base, offset, &
                                                loc_work, start_work, end_work, methread, &
                                                l_chunks(2,0:parai%node_nproc),&
                                                nmin, nmin_2, nchunk, nchunk_2, arrayshape(2)
    
    LOGICAL                                  :: second
    TYPE(C_PTR)                              :: baseptr(0:parai%node_nproc-1)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: fnlt(:,:),fiont(:,:,:,:)
#else
    REAL(real_8), ALLOCATABLE                :: fnlt(:,:),fiont(:,:,:,:)
#endif
    REAL(real_8)                             :: temp,temp1,temp2,ft(maxsys%nax,3)
    TYPE CONTAINER
       REAL(real_8),  POINTER __CONTIGUOUS :: tk(:,:)
    END TYPE CONTAINER
    TYPE(CONTAINER),ALLOCATABLE :: iproc(:)
#ifdef _VERBOSE_FORCE_DBG
    REAL(real_8),ALLOCATABLE                 :: dbg_forces(:,:,:)
#endif

    real(real_8),external :: ddot
    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlfl'
    ! split atoms between cp groups
    CALL cp_grp_split_atoms(na_grp)
    
    CALL tiset(procedureN,isub)

    IF (cntl%tfdist) CALL stopgm('FNL_RSPACE','TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (imagp.EQ.2)CALL stopgm('RNLFL','K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)

    IF(pslo_com%tivan) THEN
       
       if (parai%cp_nogrp.gt.1 .and. parai%cp_inter_me .gt. 0) then
          !$omp parallel do private(is,ia)
          do is=1,ions1%nsp
             do ia=1,ions0%na(is)
                fion(1:3,ia,is)=0._real_8
             end do
          end do
       end if
       
       tot_work=0
       DO is=1,ions1%nsp
          tot_work=tot_work+(na_grp(2,is,parai%cp_inter_me +1)-&
               na_grp(1,is,parai%cp_inter_me +1)+1)*nlps_com%ngh(is)
       END DO

       CALL dist_size(tot_work,parai%node_nproc,work_proc)     
       loc_work=work_proc(2,parai%node_me)-work_proc(1,parai%node_me)+1
       start_work=work_proc(1,parai%node_me)
       end_work=work_proc(2,parai%node_me)

       il_fnlt(1)=loc_work
       il_fnlt(2)=nstate

       il_fiont(1)=3
       il_fiont(2)=maxsys%nax
       il_fiont(3)=ions1%nsp+10 !padding
       il_fiont(4)=parai%ncpus

       ALLOCATE(iproc(0:parai%node_nproc-1), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate iproc',& 
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL request_saved_scratch(il_dfnl_packed,dfnl_packed,'DFNL_packed')
       CALL request_scratch(il_fiont,fiont,procedureN//'_fiont')
       CALL request_scratch(il_fnlt,fnlt,procedureN//'_fnlt')
#else
       ALLOCATE(fiont(il_fiont(1),il_fiont(2),il_fiont(3),il_fiont(4)), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fiont',& 
            __LINE__,__FILE__)
       ALLOCATE(fnlt(il_fnlt(1),il_fnlt(2)), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlt',& 
            __LINE__,__FILE__)
#endif
       !parallelize dgemms on node level
       CALL mp_win_alloc_shared_mem('r',tot_work,nstate,baseptr,parai%node_nproc,&
            parai%node_me,parai%node_grp)
       arrayshape(1)=tot_work
       arrayshape(2)=nstate
       DO proc=0,parai%node_nproc-1
          CALL C_F_POINTER(baseptr(proc), iproc(proc)%tk,arrayshape)
       END DO

       !set up spin settings
       IF (cntl%tlsd) THEN
          nmin=1
          nchunk=spin_mod%nsup
          nmin_2=spin_mod%nsup+1
          nchunk_2=spin_mod%nsdown
          second=.TRUE.
       ELSE
          nmin=1
          nchunk=nstate
          nmin_2=0
          nchunk_2=0
          second=.FALSE.
       ENDIF

       !pack fnl for efficient dgemm
       !$omp parallel do private(i,isa0,offset,is,iv,ia,isa,start_ia,end_ia,ia_sum) &
       !$omp proc_bind(close)
       states: DO i=1,nstate
          isa0=0
          offset=0
          DO is=1,ions1%nsp            
             start_ia=na_grp(1,is,parai%cp_inter_me +1)
             end_ia=na_grp(2,is,parai%cp_inter_me +1)
             ia_sum=end_ia-start_ia+1            
             IF (offset+nlps_com%ngh(is)*ia_sum.LT.start_work) THEN
                offset=offset+nlps_com%ngh(is)*ia_sum
             ELSE
                DO iv=1,nlps_com%ngh(is)
                   IF (offset+ia_sum.LT.start_work) THEN
                      offset=offset+ia_sum
                   ELSE
                      isa=isa0+start_ia-1
                      DO ia=1,ia_sum
                         isa=isa+1
                         offset=offset+1
                         IF (offset.GT.end_work) CYCLE states
                         IF (offset.GE.start_work) fnlt(offset-start_work+1,i)=fnla(isa,IV,I)
                      ENDDO
                   END IF
                ENDDO
             END IF
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO states

       !apply gamma to fnl 
       !something similar is done in nlforce... think about combining!      
       CALL dgemm('N','N',loc_work,nchunk,nchunk,1.0_real_8,&
            fnlt(1,nmin),loc_work,&
            gamma(nmin,nmin),nstate,&
            0.0_real_8,iproc(parai%node_me)%tk(start_work,nmin),tot_work)
       IF(second)THEN
          CALL dgemm('N','N',loc_work,nchunk_2,nchunk_2,1.0_real_8,&
               fnlt(1,nmin_2),loc_work,&
               gamma(nmin_2,nmin_2),nstate,&
               0.0_real_8,iproc(parai%node_me)%tk(start_work,nmin_2),tot_work)
       END IF
       
       !sync shared memory window
       CALL mp_win_sync(parai%node_grp)
       CALL mp_sync(parai%node_grp)
       !get data from other procs
       DO proc=0,parai%node_nproc-1
          IF (proc.NE.parai%node_me) THEN
             !$omp parallel do private(i,j)
             DO i=1,nstate
                !$omp simd
                DO j=work_proc(1,proc),work_proc(2,proc)
                   iproc(parai%node_me)%tk(j,i)=iproc(proc)%tk(j,i)
                END DO
             END DO
          END IF
       END DO

       methread=1
       !$omp parallel private (methread,is,ia,k,isa0,offset_base,start_ia,end_ia,ia_sum,ft,&
       !$omp i,iv,jv,temp,offset1,offset2,offset3,offset4,iac)
       !$ methread=omp_get_thread_num()+1
       fiont(:,:,:,methread)=0.0_real_8

       isa0=0
       offset_base=0
       DO is=1,ions1%nsp
          IF (.NOT. pslo_com%tvan(is)) THEN 
             isa0=isa0+ions0%na(is)
             CYCLE
          END IF
          
          start_ia=na_grp(1,is,parai%cp_inter_me +1)
          end_ia=na_grp(2,is,parai%cp_inter_me +1)
          ia_sum=end_ia-start_ia+1
          IF (ia_sum.GT.0) THEN            
             ft=0.0_real_8
             !$omp do 
             do i=1,nstate
                DO jv=1,nlps_com%ngh(is)
                   temp=qq(jv,jv,is)*2.0_real_8
                   do k=1,3
                      offset2=offset_base+(jv-1)*ia_sum
                      offset3=offset_base*3+(jv-1)*ia_sum*3+(k-1)*ia_sum
                      !$omp simd
                      do ia=1,ia_sum
                         ft(ia,k)=ft(ia,k)+iproc(parai%node_me)%tk(offset2+ia,i)*&
                              dfnl_packed(offset3+ia,i)*temp
                      end do
                   end DO
                end DO
                DO iv=1,nlps_com%ngh(is)
                   DO jv=iv+1,nlps_com%ngh(is)
                      temp=qq(jv,iv,is)*2.0_real_8
                      IF (ABS(temp).GT.2.e-5_real_8) THEN
                         do k=1,3
                            offset1=offset_base+(iv-1)*ia_sum
                            offset2=offset_base+(jv-1)*ia_sum
                            offset3=offset_base*3+(jv-1)*ia_sum*3+(k-1)*ia_sum
                            offset4=offset_base*3+(iv-1)*ia_sum*3+(k-1)*ia_sum
                            !$omp simd
                            do ia=1,ia_sum
                               ft(ia,k)=ft(ia,k)+iproc(parai%node_me)%tk(offset1+ia,i)*&
                                    dfnl_packed(offset3+ia,i)*temp
                               ft(ia,k)=ft(ia,k)+iproc(parai%node_me)%tk(offset2+ia,i)*&
                                    dfnl_packed(offset4+ia,i)*temp
                            end do
                         end do
                      end IF
                   end DO
                end DO
             ENDDO
             !$omp end do nowait
             DO ia=1,ia_sum
                do k=1,3
                   fiont(k,ia,is,methread)=ft(ia,k)
                end do
             end do
          END IF
          offset_base=offset_base+nlps_com%ngh(is)*ia_sum
          isa0=isa0+ions0%na(is)
       END DO
       !$omp end parallel
       do methread=1,parai%ncpus
          DO is=1,ions1%nsp
             do ia=na_grp(1,is,parai%cp_inter_me+1),na_grp(2,is,parai%cp_inter_me+1)
                iac=ia-na_grp(1,is,parai%cp_inter_me+1)+1
                do k=1,3
                   fion(k,ia,is)=fion(k,ia,is)-fiont(k,iac,is,methread)
                end do
             end do
          end DO
       end do

#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_fnlt,fnlt,procedureN//'_fnlt')
       CALL free_scratch(il_fiont,fiont,procedureN//'_fiont')
       CALL free_scratch(il_dfnl_packed,dfnl_packed,'DFNL_packed')
#else
       DEALLOCATE(fiont, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fiont',& 
            __LINE__,__FILE__)
       DEALLOCATE(fnlt, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fnlt',& 
            __LINE__,__FILE__)
       DEALLOCATE(dfnl_packed, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate dfnl_packed',& 
            __LINE__,__FILE__)
#endif
       DEALLOCATE(iproc, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate iproc',& 
            __LINE__,__FILE__)

       IF (parai%cp_nogrp.GT.1 ) THEN
          CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%cp_inter_grp)
       END IF

    END IF

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
  END SUBROUTINE rnlfl
  ! ==================================================================

END MODULE rnlfl_utils
