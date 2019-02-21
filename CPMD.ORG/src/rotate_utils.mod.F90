#include "cpmd_global.h"

MODULE rotate_utils
  USE error_handling,                  ONLY: stopgm
  USE cp_grp_utils,                    ONLY: cp_grp_redist_array,&
                                             cp_grp_split_states  
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE nvtx_utils
  USE parac,                           ONLY: parai
  USE part_1d,                         ONLY: part_1d_get_blk_bounds
  USE system,                          ONLY: parap
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rotate
  !!public :: rotate_c
  !!public :: rotate_da
  PUBLIC :: rottr
  PUBLIC :: rottr_fnl

CONTAINS
  ! ==================================================================
  SUBROUTINE rottr_fnl(a,fnl,gam,transa,nstate,n,tlsd,na,nb,gid,np,me)
    ! ==--------------------------------------------------------------==
    ! ==         FNL <= A*FNL*GAM                                     ==
    ! ==   SPECIAL CASE FOR GAM UPPER TRIAGONAL                       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8),INTENT(INOUT)               :: fnl(n,*)
    REAL(real_8),INTENT(IN)                  :: a, gam(:,:)
    CHARACTER(len=*),INTENT(IN)              :: transa
    LOGICAL,INTENT(IN)                       :: tlsd
    INTEGER,INTENT(IN)                       :: nstate, n, na, nb, gid, np, me

    INTEGER                                  :: isub, naa, i, n_loc,&
                                                l_chunks(2,0:np-1), n_start, n_end, ierr,&
                                                t_chunks(2,0:np-1)
    CHARACTER(*), PARAMETER                  :: procedureN = 'rottr_fnl'
! ==--------------------------------------------------------------==
    IF (n.EQ.0 .OR. nstate.EQ.0) RETURN
    CALL tiset(procedureN,isub)
    IF(np.GT.1)THEN
       DO i = 0,np-1
          CALL part_1d_get_blk_bounds(n,i,np,l_chunks(1,i),l_chunks(2,i))
       END DO
       n_end=l_chunks(2,me)
       n_start=l_chunks(1,me)
    ELSE
       n_end=n
       n_start=1
    END IF
    n_loc=n_end-n_start+1
    IF (tlsd) THEN
       naa=na+1
       CALL dtrmm('R','U',transa,'N',n_loc,na,a,gam,nstate,fnl(n_start,1),n)
       CALL dtrmm('R','U',transa,'N',n_loc,nb,a,gam(naa,naa),nstate,&
            fnl(n_start,naa),n)
    ELSE
       CALL dtrmm('R','U',transa,'N',n_loc,nstate,a,gam,nstate,fnl(n_start,1),n)
    END IF
    IF (np .GT. 1) THEN
       t_chunks(1,:)=1
       t_chunks(2,:)=nstate
       call mp_allgatherv_intranode_inplace_real(fnl,n,nstate,l_chunks,t_chunks,np,me,gid)
    END IF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rottr_fnl
  ! ==================================================================
  SUBROUTINE rottr(a,c1,gam,transa,nstate,n,tlsd,na,nb)
    ! ==--------------------------------------------------------------==
    ! ==         C1 <= A*C1*GAM                                       ==
    ! ==   SPECIAL CASE FOR GAM UPPER TRIAGONAL                       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8),intent(in)                  :: a
    COMPLEX(real_8),TARGET                   :: c1(:,:)
    REAL(real_8)                             :: gam(:,:)
    CHARACTER(len=*)                         :: transa
    INTEGER                                  :: nstate, n
    LOGICAL                                  :: tlsd
    INTEGER                                  :: na, nb

    INTEGER                                  :: isub, naa, il_c1_loc(3), ig, ig1, i, n_loc,&
                                                n_proc(3,0:parai%cp_nogrp-1), grp, n_start, ierr
    CHARACTER(*), PARAMETER                  :: procedureN = 'rottr'
    COMPLEX(real_8),ALLOCATABLE              :: c1_loc(:,:,:)
!(nstate,nstate)
! Variables
! ==--------------------------------------------------------------==
    IF (n.EQ.0 .OR. nstate.EQ.0) RETURN
    CALL tiset('     ROTTR',isub)
    IF (tlsd) THEN
       naa=na+1
       CALL dtrmm('R','U',transa,'N',2*n,na,a,gam,nstate,c1,2*n)
       CALL dtrmm('R','U',transa,'N',2*n,nb,a,gam(naa,naa),nstate,&
            c1(1,naa),2*n)
    ELSE
       IF (parai%cp_nogrp .gt. 1) THEN
          DO i = 0,parai%cp_nogrp-1
             CALL part_1d_get_blk_bounds(n,i,parai%cp_nogrp,n_proc(1,i),n_proc(2,i))
             n_proc(3,i)=n_proc(2,i)- n_proc(1,i)+1
          END DO

          il_c1_loc(1)=MAXVAL(n_proc(3,:))
          il_c1_loc(2)=nstate
          il_c1_loc(3)=parai%cp_nogrp
          ALLOCATE(c1_loc(il_c1_loc(1),il_c1_loc(2),il_c1_loc(3)),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          n_loc=n_proc(3,parai%cp_inter_me)
          n_start=n_proc(1,parai%cp_inter_me)
       ELSE
          n_loc=n
          n_start=1
       END if

       
       CALL dtrmm('R','U',transa,'N',2*n_loc,&
            nstate,a,gam,nstate,c1(n_start,1),2*n)

       IF (parai%cp_nogrp .GT. 1) THEN
          !$omp parallel do private (i,ig,ig1) 
          DO i=1,nstate
             DO ig=n_proc(1,parai%cp_inter_me),n_proc(2,parai%cp_inter_me)
                ig1=ig-n_proc(1,parai%cp_inter_me)+1
                c1_loc(ig1,i,parai%cp_inter_me+1)=c1(ig,i)
             END DO
          END DO

          CALL my_concat_inplace(c1_loc,il_c1_loc(1)*nstate*2,parai%cp_inter_grp)

          !$omp parallel private(i,ig,ig1,grp)
          DO grp=1,parai%cp_nogrp
             IF (grp .NE. parai%cp_inter_me+1) THEN
                !$omp do 
                DO i=1,nstate
                   DO ig=n_proc(1,grp-1),n_proc(2,grp-1)
                      ig1=ig-n_proc(1,grp-1)+1
                      c1(ig,i)=c1_loc(ig1,i,grp)
                   END DO
                END DO
                !$omp end do nowait
             END IF
          END DO
          !$omp end parallel
          DEALLOCATE(c1_loc,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
       END IF
    END IF
    CALL tihalt('     ROTTR',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rottr
  ! ==================================================================
  SUBROUTINE rotate(a,c1,b,c2,gam,nstate,n,tlsd,na,nb,redist)
    ! ==--------------------------------------------------------------==
    ! ==         C2 <= B*C2 + A*C1*GAM                                ==
    ! ==   ALSO FOR LSD                                               ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: a
    COMPLEX(real_8)                          :: c1(:,:)
    REAL(real_8)                             :: b
    COMPLEX(real_8)                          :: c2(:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: gam(nstate,*)
    INTEGER                                  :: n
    LOGICAL                                  :: tlsd
    INTEGER                                  :: na, nb
    LOGICAL,INTENT(IN),OPTIONAL              :: redist

    CHARACTER(*), PARAMETER                  :: procedureN = 'rotate'

    INTEGER                                  :: isub, isub1, naa
    INTEGER                                  :: first_state,nstate_grp
    LOGICAL                                  :: redst
    ! split states between cp groups


! Variables
! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )

    IF (PRESENT(redist)) THEN
       redst=redist
    ELSE
       redst=.true.
    END IF

    IF (n>0) THEN
       IF (tlsd) THEN
          naa=na+1
          CALL dgemm('N','N',n,na,na,a,c1(1,1),n,gam(1,1),nstate,b,c2(1,1),n)
          CALL dgemm('N','N',n,nb,nb,a,c1(1,naa),n,gam(naa,naa),nstate,b,c2(1,naa),n)
       ELSE
          call cp_grp_split_states(nstate,nstate_grp=nstate_grp,first_state=first_state)
          CALL dgemm('N','N',n,nstate_grp,nstate,a,c1,n,gam(1,first_state),nstate,b,c2(1,first_state),n)
       ENDIF
    ENDIF
    
    IF (parai%cp_nogrp.GT.1.AND..NOT.tlsd) THEN
       !fixme: tlsd requieres parap%first_state_up and parap%first_state_down
       IF (redst) THEN
          CALL tiset(procedureN//'_grpsb',isub1)
          CALL cp_grp_redist_array(c2,n/2,nstate)
          CALL tihalt(procedureN//'_grpsb',isub1)
       END IF
    END IF
    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rotate



END MODULE rotate_utils

! ==================================================================
SUBROUTINE rotate_da(a,c1,b,c2,gam,ld,n,nstate,ndd1,ndd2,nddx,mepos,pgroup,nproc,grp,tlsd,na,nb)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_bcast
  USE cuda_types,                      ONLY: cuda_memory_t
  USE cuda_utils,                      ONLY: cuda_alloc_bytes, cuda_dealloc, &
       cuda_memcpy_host_to_device, cuda_memcpy_device_to_host, cuda_d_points_to, cuda_mem_zero_bytes
  USE cublas_types,                    ONLY: cublas_handle_t
  USE cublas_utils,                    ONLY: cublas_create, cublas_destroy, cublas_dgemm, cublas_dger
  USE sizeof_kinds,                    ONLY: sizeof_real_8
  USE cp_cuwfn_types, ONLY: cp_cuwfn_device_get_ptrs, cp_cuwfn, cp_cuwfn_get
  USE nvtx_utils
  IMPLICIT NONE
  REAL(real_8)                               :: a, b
  INTEGER                                    :: ld, n, nstate
  REAL(real_8)                               :: c1(ld,nstate), c2(ld,nstate), &
                                                gam(nstate,nstate)
  INTEGER                                    :: ndd1(0:*), ndd2(0:*), nddx, &
                                                mepos, pgroup(0:*), nproc, grp
  LOGICAL                                    :: tlsd
  INTEGER                                    :: na, nb

  CHARACTER(*), PARAMETER                    :: procedureN = 'rotate_da'
  LOGICAL, PARAMETER                         :: use_gpu = .FALSE.

  INTEGER                                    :: i, i1, ierr, ii, ip, isub, j, &
                                                n1, nmax, nmin
  INTEGER(int_8)                             :: n_bytes
  REAL(real_8), ALLOCATABLE                  :: aux(:)
  REAL(real_8), ALLOCATABLE, DIMENSION(:, :) :: c2_tmp
  TYPE(cublas_handle_t)                      :: blas_handle
  TYPE(cuda_memory_t)                        :: c1_d, c2_d, gam_d

! ==--------------------------------------------------------------==

  CALL tiset(procedureN,isub)
  __NVTX_TIMER_START ( procedureN )

  IF (tlsd) THEN
     IF (ndd1(mepos).LT.na) THEN
        nmax=MIN(na,ndd2(mepos))
        DO i=ndd1(mepos),nmax
           ii=i-ndd1(mepos)+1
           DO j=na+1,nstate
              gam(j,ii)=0._real_8
           ENDDO
        ENDDO
     ENDIF
     IF (ndd2(mepos).GT.na) THEN
        nmin=MAX(na+1,ndd1(mepos))
        DO i=nmin,ndd2(mepos)
           ii=i-ndd1(mepos)+1
           DO j=1,na
              gam(j,ii)=0._real_8
           ENDDO
        ENDDO
     ENDIF
  ENDIF
  CALL dscal(nstate*n,b,c2,1)
  ALLOCATE(aux(nstate*nddx),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  CALL dcopy(nstate*nddx,gam,1,aux,1)

  IF( use_gpu ) THEN
     WRITE(*,*) 'rotate_utils.mod.F90: use_gpu',use_gpu
     CALL cublas_create ( blas_handle, 0 )

     n_bytes = SIZE( c1, KIND=int_8 ) * INT( sizeof_real_8, int_8 )
     CALL cuda_alloc_bytes ( c1_d, n_bytes, 0 )
     CALL cuda_memcpy_host_to_device ( c1, c1_d )

     n_bytes = SIZE( c1, KIND=int_8 ) * INT( sizeof_real_8, int_8 )
     CALL cuda_alloc_bytes ( c2_d, n_bytes, 0 )
     CALL cuda_mem_zero_bytes ( c2_d, n_bytes )

     n_bytes = SIZE( gam, KIND=int_8 ) * INT( sizeof_real_8, int_8 )
     CALL cuda_alloc_bytes ( gam_d, n_bytes, 0 )
  ENDIF

  DO ip=0,nproc-1
     i1=ndd1(ip)
     n1=ndd2(ip)-ndd1(ip)+1
     IF (n1.GT.0) THEN
        CALL dcopy(nstate*nddx,aux,1,gam,1)
        CALL mp_bcast(gam,nstate*nddx,pgroup(ip+1),grp)
        IF (ld>0) THEN
           IF( use_gpu ) THEN
              CALL cuda_memcpy_host_to_device ( gam, gam_d )
              CALL cublas_dgemm ( blas_handle, 'N', 'N', n, n1, nstate, &
                   & a, c1_d, ld, &
                   & gam_d, nstate, &
                   & 1.0_real_8, cuda_d_points_to(c2_d,(i1-1)*ld+1), ld )
           ELSE
              CALL dgemm('N','N',n,n1,nstate,a,c1(1,1),ld,&
                   gam(1,1),nstate,1._real_8,c2(1,i1),ld)
           ENDIF
        ENDIF
     ENDIF
  ENDDO

  IF( use_gpu ) THEN
     ALLOCATE(c2_tmp(SIZE(c2,1),SIZE(c2,2)),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
          __LINE__,__FILE__)
     CALL cuda_memcpy_device_to_host ( c2_d, c2_tmp )
     CALL daxpy ( SIZE(c2), 1.0_real_8, c2_tmp, 1, c2, 1 )
     DEALLOCATE(c2_tmp,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
          __LINE__,__FILE__)
     CALL cuda_dealloc ( gam_d )
     CALL cuda_dealloc ( c2_d )
     CALL cuda_dealloc ( c1_d )
     CALL cublas_destroy ( blas_handle )
  ENDIF

  CALL dcopy(nstate*nddx,aux,1,gam,1)
  DEALLOCATE(aux,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)

  __NVTX_TIMER_STOP
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE rotate_da
! ==================================================================
SUBROUTINE rotate_c(a,c1,b,c2,gam,nstate)
  ! ==--------------------------------------------------------------==
  ! ==         C2 <= B*C2 + A*C1*GAM                                ==
  ! ==   ALSO FOR LSD                                               ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:cntl,nkpt
  USE parac, ONLY : paral,parai
  USE spin , ONLY:spin_mod
  IMPLICIT NONE
  COMPLEX(real_8)                            :: a, c1(nkpt%ngwk,*), b, &
                                                c2(nkpt%ngwk,*)
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: gam(nstate,*)

  INTEGER                                    :: isub, naa

  IF (nkpt%ngwk.EQ.0 .OR. nstate.EQ.0) RETURN
  CALL tiset('  ROTATE_C',isub)
  IF (cntl%tlsd) THEN
     naa=spin_mod%nsup+1
     IF (nkpt%ngwk>0) THEN
        CALL zgemm('N','N',nkpt%ngwk,spin_mod%nsup,spin_mod%nsup,a,c1(1,1),nkpt%ngwk,gam(1,1),&
             nstate,b,c2(1,1),nkpt%ngwk)
        CALL zgemm('N','N',nkpt%ngwk,spin_mod%nsdown,spin_mod%nsdown,a,c1(1,naa),nkpt%ngwk,&
             gam(naa,naa),nstate,b,c2(1,naa),nkpt%ngwk)
     ENDIF
  ELSE
     IF (nkpt%ngwk>0) THEN
        CALL zgemm('N','N',nkpt%ngwk,nstate,nstate,a,c1(1,1),nkpt%ngwk,&
             gam(1,1),nstate,b,c2(1,1),nkpt%ngwk)
     ENDIF
  ENDIF
  CALL tihalt('  ROTATE_C',isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE rotate_c
