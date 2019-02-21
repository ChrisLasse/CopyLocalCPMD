#include "cpmd_global.h"

MODULE crotwf_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE part_1d,                         ONLY: part_1d_get_blk_bounds
  USE pslo,                            ONLY: pslo_com
  USE nort,                            ONLY: nort_ovlap,&
                                             nort_com
  USE rotate_utils,                    ONLY: rotate
  USE spin,                            ONLY: spin_mod
  USE summat_utils,                    ONLY: summat
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE utils,                           ONLY: dsyevx_driver,&
                                             dsyevd_driver
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE
#ifdef _HAS_LIBELPA
  PUBLIC :: crotwf,crotwf_para
#else
  PUBLIC :: crotwf
#endif
  PUBLIC :: give_scr_crotwf

CONTAINS
#ifdef _HAS_LIBELPA
  SUBROUTINE crotwf_para(c0,cm,c2,sc0,nstate,gam)
    ! ==--------------------------------------------------------------==
    USE elpa
    USE elpa_utils,                           ONLY: elpa_ob_create
    IMPLICIT NONE
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: sc0(ncpw%ngw,nstate), c2(ncpw%ngw,nstate), &
      cm(ncpw%ngw,nstate), c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: gam(nstate,nstate)

    INTEGER                                  :: i, ic1b, ierr, iopt, j, k,success
    REAL(real_8), ALLOCATABLE                :: aux(:,:), c1(:,:), w(:)
    
    class(elpa_t),SAVE,pointer :: el1,el2
    INTEGER, SAVE, allocatable               :: idx1(:,:,:),idx2(:,:,:)
    INTEGER,SAVE                             :: na_row1,na_row2,na_col1,na_col2
    integer                                  :: isub
    CHARACTER(*), PARAMETER                  :: procedureN = 'crotwf_para'    
    
    CALL tiset(procedureN,isub)

!already done in csmat -> nort_ovlap    
!    IF (pslo_com%tivan) THEN
!        CALL ovlap(nstate,gam,c0,c0,.false.)
!        CALL summat(gam,nstate,.true.,parai%cp_grp)
!    ELSE
!        CALL ovlap(nstate,gam,c0,c0,.true.)
!        CALL summat(gam,nstate,.true.,parai%allgrp)
!    END IF

    call zeroing(gam)
    
    IF (.NOT.cntl%tlsd) THEN
      
      if(.not.ALLOCATED(idx1)) call elpa_ob_create(nstate,el1,idx1,na_row1,na_col1)    

      if(na_row1.GT.0) then
        ALLOCATE(c1(na_row1,na_col1),w(nstate),aux(na_row1,na_col1),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
              __LINE__,__FILE__)
        
        !&OMP parallel do private(i,j,k,l)
        do j=1,na_col1
            DO i=1,na_row1
              c1(i,j)=nort_ovlap(idx1(i,j,1), idx1(i,j,2))
            end do
        end DO 
        call el1%eigenvectors(c1, w, aux, success)

        !$OMP parallel do private(i,j)
        do j=1,na_col1
          do i=1,na_row1
            gam(idx1(i,j,1),idx1(i,j,2))=aux(i,j)
          enddo
        enddo
    
      END IF
    ELSE
       
      if(.not.ALLOCATED(idx1))  then
        call elpa_ob_create(spin_mod%nsup  ,el1,idx1,na_row1,na_col1) 
        call elpa_ob_create(spin_mod%nsdown,el2,idx2,na_row2,na_col2)
        ! shift idx2 to the beginning of the down spins
        idx2(:,:,:)=idx2(:,:,:)+spin_mod%nsup
      endif
      if(na_row1.GT.0) then
        ALLOCATE(c1(na_row1,na_col1),w(spin_mod%nsup),aux(na_row1,na_col1),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
              __LINE__,__FILE__)
        !&OMP parallel do private(i,j,k,l)
        do j=1,na_col1
            DO i=1,na_row1
              c1(i,j)=nort_ovlap(idx1(i,j,1), idx1(i,j,2))
            end do
        end DO         
        call el1%eigenvectors(c1, w, aux, success)
        !$OMP parallel do private(i,j)
        do j=1,na_col1
          do i=1,na_row1
            gam(idx1(i,j,1),idx1(i,j,2))=aux(i,j)
          enddo
        enddo
        DEALLOCATE(c1,w,aux,STAT=ierr)
        ALLOCATE(c1(na_row2,na_col2),w(spin_mod%nsdown),aux(na_row2,na_col2),STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
              __LINE__,__FILE__)
          
        !&OMP parallel do private(i,j,k,l)
        do j=1,na_col2
            do i=1,na_row2
              c1(i,j)=nort_ovlap(idx2(i,j,1),idx2(i,j,2))
            end do
        end do          
        call el2%eigenvectors(c1, w, aux, success)
        !$OMP parallel do private(i,j)
        do j=1,na_col2
          do i=1,na_row2
            gam(idx2(i,j,1),idx2(i,j,2))=aux(i,j)
          enddo
        enddo
      endif

    ENDIF
    call mp_sum(gam,nstate*nstate,parai%cp_grp)
    ! to avoid problems we bcast the result (we dont need the eigvals)
!    CALL mp_bcast(gam,nstate**2,parai%io_source,parai%cp_grp)

    CALL rotate(1.0_real_8,c0,0.0_real_8,sc0,gam,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    CALL dcopy(2*ncpw%ngw*nstate,sc0(1,1),1,c0(1,1),1)
    CALL rotate(1.0_real_8,cm,0.0_real_8,sc0,gam,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    CALL dcopy(2*ncpw%ngw*nstate,sc0(1,1),1,cm(1,1),1)
    CALL rotate(1.0_real_8,c2,0.0_real_8,sc0,gam,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    CALL dcopy(2*ncpw%ngw*nstate,sc0(1,1),1,c2(1,1),1)

!    DEALLOCATE(nort_ovlap,STAT=ierr)
!    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
!         __LINE__,__FILE__)
    if(na_row1.GT.0) then
    ! ==--------------------------------------------------------------==
      DEALLOCATE(c1,w,aux,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
          __LINE__,__FILE__)
    endif
    
    CALL tihalt(procedureN,isub)

    RETURN  
  END SUBROUTINE
          
#endif
  ! ==================================================================
  SUBROUTINE crotwf(c0,cm,c2,sc0,nstate,gam)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: sc0(ncpw%ngw,nstate), c2(ncpw%ngw,nstate), &
      cm(ncpw%ngw,nstate), c0(ncpw%ngw,nstate)
    REAL(real_8)                             :: gam(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'crotwf'

    INTEGER                                  :: i, ierr, iopt, j, isub, &
                                                il_eigval(1), il_temp(2), &
                                                first,last,&
                                                chunks(2,0:parai%cp_nproc-1),&
                                                recvcnt(0:parai%cp_nproc-1),&
                                                displ(0:parai%cp_nproc-1)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: eigval(:),temp(:,:),temp1(:,:)
#else
    REAL(real_8), ALLOCATABLE                :: eigval(:),temp(:,:),temp1(:,:)
#endif
    CALL tiset(procedureN,isub)

    IF(cntl%tlsd)THEN
       il_eigval(1)=max(spin_mod%nsup,spin_mod%nsdown)
    ELSE
       il_eigval(1)=nstate
    END IF

#ifdef _USE_SCRATCHLIBRARY    
    CALL request_scratch(il_eigval,eigval,procedureN//'_eigval')
#else
    ALLOCATE(eigval(il_eigval(1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
#endif
!ovlap already done in csmat => nort_ovlap
!    CALL ovlap(nstate,gam,c0,c0)
!    IF (.NOT.cntl%tlsd) THEN
!       k=1
!       DO j=1,nstate
!          DO i=j,nstate
!             c1(k)=gam(i,j)
!             k=k+1
!          ENDDO
!       ENDDO
!       CALL mp_sum(c1,nstate*nstate,parai%allgrp)
    iopt=21
    IF (.NOT.cntl%tlsd) THEN
       IF(nort_com%scond.LT.1.e-9_real_8.OR.parai%cp_nproc.LT.17) THEN
          !parallelization using multiple dsyevx/r does not work => fall back to dsyevd on root
          IF(paral%io_parent) CALL dsyevd_driver(iopt,nort_ovlap,eigval,nstate)
          CALL mp_bcast(nort_ovlap,nstate**2,parai%io_source,parai%cp_grp)
       ELSE
          recvcnt=-1
          displ=0
          DO i = 0,parai%cp_nproc-1
             CALL part_1d_get_blk_bounds(nstate,i,parai%cp_nproc,chunks(1,i),chunks(2,i))       
             recvcnt(i)=(chunks(2,i)-chunks(1,i)+1)*nstate
             IF (i.GT.0) displ(i)=displ(i-1)+recvcnt(i-1)
          END DO
          
          first=chunks(1,parai%cp_me)
          last=chunks(2,parai%cp_me)
          call dsyevx_driver(iopt,nort_ovlap,gam,eigval,nstate,first,last,-1.0_real_8)
          CALL my_concatv(gam,nort_ovlap,(last-first+1)*nstate,recvcnt,displ,parai%cp_grp)
       END if
    ELSE
       !for the moment we copy both spins out in a temporary buffer...
       !spin up
       il_temp(1)=spin_mod%nsup
       il_temp(2)=spin_mod%nsup
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_temp,temp,procedureN//'_temp')
       CALL request_scratch(il_temp,temp1,procedureN//'_temp1')
#else
       ALLOCATE(temp(il_temp(1),il_temp(2)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(temp1(il_temp(1),il_temp(2)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)      
#endif
       !$omp parallel do private(i,j)
       DO i=1,spin_mod%nsup
          DO j=1,i
             temp(j,i)=nort_ovlap(j,i)
          END DO
       END DO
       IF(nort_com%scond.LT.1.e-9_real_8.OR.parai%cp_nproc.LT.17) THEN
          !parallelization using multiple dsyevx/r does not work => fall back to dsyevd on root
          IF(paral%io_parent) CALL dsyevd_driver(iopt,temp,eigval,spin_mod%nsup)
          CALL mp_bcast(temp,spin_mod%nsup**2,parai%io_source,parai%cp_grp)
       ELSE
          recvcnt=-1
          displ=0
          DO i = 0,parai%cp_nproc-1
             CALL part_1d_get_blk_bounds(spin_mod%nsup,i,parai%cp_nproc,chunks(1,i),chunks(2,i))       
             recvcnt(i)=(chunks(2,i)-chunks(1,i)+1)*spin_mod%nsup
             IF (i.GT.0) displ(i)=displ(i-1)+recvcnt(i-1)
          END DO
          
          first=chunks(1,parai%cp_me)
          last=chunks(2,parai%cp_me)
          CALL dsyevx_driver(iopt,temp,temp1,eigval,spin_mod%nsup,first,last,-1.0_real_8)
          CALL my_concatv(temp1,temp,(last-first+1)*spin_mod%nsup,recvcnt,displ,parai%cp_grp)
       END if
       !$omp parallel do private(i,j)
       DO i=1,spin_mod%nsup
          DO j=1,spin_mod%nsup
             nort_ovlap(j,i)=temp(j,i)
          END DO
       END DO
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_temp,temp1,procedureN//'_temp1')
       CALL free_scratch(il_temp,temp,procedureN//'_temp')
#else
       DEALLOCATE(temp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(temp1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
#endif
       !spin down
       il_temp(1)=spin_mod%nsdown
       il_temp(2)=spin_mod%nsdown
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_temp,temp,procedureN//'_temp')
       CALL request_scratch(il_temp,temp1,procedureN//'_temp1')
#else
       ALLOCATE(temp(il_temp(1),il_temp(2)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(temp1(il_temp(1),il_temp(2)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)      
#endif
       !$omp parallel do private(i,j)
       DO i=spin_mod%nsup+1,nstate
          DO j=spin_mod%nsup+1,i
             temp(j-spin_mod%nsup,i-spin_mod%nsup)=nort_ovlap(j,i)
          END DO
       END DO
       IF(nort_com%scond.LT.1.e-9_real_8.OR.parai%cp_nproc.LT.17) THEN
          !parallelization using multiple dsyevx/r does not work => fall back to dsyevd on root
          IF(paral%io_parent) CALL dsyevd_driver(iopt,temp,eigval,spin_mod%nsdown)
          CALL mp_bcast(temp,spin_mod%nsdown**2,parai%io_source,parai%cp_grp)
       ELSE
          recvcnt=-1
          displ=0
          DO i = 0,parai%cp_nproc-1
             CALL part_1d_get_blk_bounds(spin_mod%nsdown,i,parai%cp_nproc,chunks(1,i),chunks(2,i))       
             recvcnt(i)=(chunks(2,i)-chunks(1,i)+1)*spin_mod%nsdown
             IF (i.GT.0) displ(i)=displ(i-1)+recvcnt(i-1)
          END DO
          
          first=chunks(1,parai%cp_me)
          last=chunks(2,parai%cp_me)
          CALL dsyevx_driver(iopt,temp,temp1,eigval,spin_mod%nsdown,first,last,-1.0_real_8)
          CALL my_concatv(temp1,temp,(last-first+1)*spin_mod%nsdown,recvcnt,displ,parai%cp_grp)
       END if
       !$omp parallel do private(i,j)
       DO i=spin_mod%nsup+1,nstate
          DO j=spin_mod%nsup+1,nstate
             nort_ovlap(j,i)=temp(j-spin_mod%nsup,i-spin_mod%nsup)
          END DO
       END DO
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_temp,temp1,procedureN//'_temp1')
       CALL free_scratch(il_temp,temp,procedureN//'_temp')
#else
       DEALLOCATE(temp,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(temp1,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
#endif
    ENDIF

    CALL rotate(1.0_real_8,c0,0.0_real_8,sc0,nort_ovlap,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    CALL dcopy(2*ncpw%ngw*nstate,sc0(1,1),1,c0(1,1),1)
    CALL rotate(1.0_real_8,cm,0.0_real_8,sc0,nort_ovlap,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    CALL dcopy(2*ncpw%ngw*nstate,sc0(1,1),1,cm(1,1),1)
    CALL rotate(1.0_real_8,c2,0.0_real_8,sc0,nort_ovlap,nstate,2*ncpw%ngw,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
    CALL dcopy(2*ncpw%ngw*nstate,sc0(1,1),1,c2(1,1),1)
    ! ==--------------------------------------------------------------==
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_eigval,eigval,procedureN//'_eigval')
#else
    DEALLOCATE(eigval,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#endif
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE crotwf
  ! ==================================================================
  SUBROUTINE give_scr_crotwf(lcrotwf,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lcrotwf
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

! ==--------------------------------------------------------------==

    lcrotwf=nstate*nstate+4*nstate
    tag   ='NSTATE*NSTATE+4*NSTATE'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_crotwf
  ! ==================================================================

END MODULE crotwf_utils
