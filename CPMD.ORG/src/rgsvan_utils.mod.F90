#include "cpmd_global.h"

MODULE rgsvan_utils
  USE csmat_utils,                     ONLY: csmat
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms  
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE fnl_utils,                       ONLY: unpack_fnl
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE rgs_utils,                       ONLY: uinv
  USE rotate_utils,                    ONLY: rottr,&
                                             rottr_fnl
  USE rnlsm_utils,                     ONLY: rnlsm
  USE sfac,                            ONLY: ldf1,&
                                             fnla,&
                                             fnl,&
                                             fnl_packed,&
                                             il_fnl_packed
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean,&
                                             symmat_pack,&
                                             symmat_unpack

#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch,&
                                             request_saved_scratch
#endif
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rgsvan

CONTAINS

  ! ==================================================================
  SUBROUTINE rgsvan(c0,nstate,smat,store_nonort,need_fnl,rot_fnl)
    ! ==--------------------------------------------------------------==
    ! ==  Gram-Schmidt orthogonalization for Vanderbilt pp            ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: smat(nstate,nstate)
    INTEGER                                  :: il_smatpacked(1), ierr, isub,&
                                                na_grp(2,ions1%nsp,parai%cp_nogrp)
    CHARACTER(*), PARAMETER                  :: procedureN = 'rgsvan'
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: smatpacked(:)
#else
    REAL(real_8), ALLOCATABLE                :: smatpacked(:)
#endif
    LOGICAL, INTENT(IN)                      :: store_nonort
    LOGICAL,INTENT(IN),OPTIONAL              :: rot_fnl,need_fnl
    LOGICAL                                  :: fnl_rot
! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    CALL cp_grp_split_atoms(na_grp)
    IF(PRESENT(rot_fnl))THEN
       fnl_rot=rot_fnl
    ELSE
       fnl_rot=.TRUE.
    END IF
    IF(PRESENT(need_fnl))THEN
       IF(need_fnl) CALL rnlsm(c0,nstate,1,1,.FALSE.,rdst=.FALSE.,packed_fnl=rot_fnl)
    ELSE
       CALL rnlsm(c0,nstate,1,1,.FALSE.,rdst=.FALSE.,packed_fnl=fnl_rot)
    END IF
    CALL csmat(smat,c0,fnla,nstate,1,store_nonort,rdst=.FALSE.,parent=.TRUE.)
    IF(cntl%tlsd) THEN
       il_smatpacked=spin_mod%nsup*(spin_mod%nsup+1)/2+&
            spin_mod%nsdown*(spin_mod%nsdown+1)/2
    ELSE
       il_smatpacked=nstate*(nstate+1)/2
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_smatpacked,smatpacked,procedureN//'_smatpacked')
#else
    ALLOCATE(smatpacked(il_smatpacked(1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
#endif
    IF(paral%io_parent)THEN
       IF (cntl%tlsd) THEN
          CALL uinv('U',smat(1,1),nstate,spin_mod%nsup)
          CALL uinv('U',smat(spin_mod%nsup+1,spin_mod%nsup+1),nstate,spin_mod%nsdown)
          CALL symmat_pack(smat,smatpacked,nstate,spin_mod%nsup,spin_mod%nsdown)
       ELSE
          CALL uinv('U',smat,nstate,nstate)
          CALL symmat_pack(smat,smatpacked,nstate,nstate,0)
       ENDIF
    END IF
    CALL mp_bcast(smatpacked,il_smatpacked(1),parai%io_source,parai%cp_grp)
    IF(cntl%tlsd)THEN
       CALL symmat_unpack(smat,smatpacked,nstate,spin_mod%nsup,spin_mod%nsdown,.FALSE.)
    ELSE
       CALL symmat_unpack(smat,smatpacked,nstate,nstate,0,.FALSE.)
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_smatpacked,smatpacked,procedureN//'_smatpacked')
#else
    DEALLOCATE(smatpacked,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#endif
    IF (ncpw%ngw.GT.0)&
         CALL rottr(1._real_8,c0,smat,"N",nstate,ncpw%ngw,cntl%tlsd,spin_mod%nsup,&
         spin_mod%nsdown)

    IF(fnl_rot)THEN
#ifdef _USE_SCRATCHLIBRARY
       CALL request_saved_scratch(il_fnl_packed,fnl_packed,'fnl_packed')
#endif
       CALL rottr_fnl(1._real_8,fnl_packed,smat,"N",nstate,size(fnl_packed,1),cntl%tlsd,&
            spin_mod%nsup,spin_mod%nsdown,parai%node_grp,parai%node_nproc,parai%node_me)
       CALL unpack_fnl(fnl_packed,fnl,fnla,na_grp,1,nstate,.FALSE.,1)
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_fnl_packed,fnl_packed,'fnl_packed')
#else
       DEALLOCATE(fnl_packed,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
#endif
    END IF
    IF (geq0) CALL zclean(c0,nstate,ncpw%ngw)

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rgsvan
  ! ==================================================================

END MODULE rgsvan_utils
