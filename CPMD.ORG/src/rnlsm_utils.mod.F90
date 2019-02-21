MODULE rnlsm_utils
  USE cp_grp_utils,                    ONLY: cp_grp_redist_dfnl_fnl
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nlm
  USE rnlsm1_utils,                    ONLY: give_scr_rnlsm1,&
                                             rnlsm1
  USE rnlsm2_utils,                    ONLY: give_scr_rnlsm2,&
                                             rnlsm2
  USE rnlsmd_utils,                    ONLY: rnlsmd
  USE system,                          ONLY: cntl
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlsm
  PUBLIC :: give_scr_rnlsm

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlsm(c0,nstate,ikpt,ikind,tfor,rdst,only_dfnl,packed_fnl,packed_dfnl)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(:,:)
    INTEGER                                  :: nstate, ikpt, ikind
    LOGICAL,INTENT(IN)                       :: tfor
    LOGICAL, INTENT(IN), OPTIONAL            :: rdst,only_dfnl,packed_fnl,packed_dfnl
    LOGICAL                                  :: redist,dfnl,p_dfnl,p_fnl
    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlsm'

    INTEGER                                  :: isub

! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF(nlm.EQ.0) RETURN
    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    IF(PRESENT(rdst))THEN
       redist=rdst
    ELSE
       redist=.TRUE.
    END IF
    IF(PRESENT(packed_fnl))THEN
       p_fnl=packed_fnl
    ELSE
       p_fnl=.FALSE.
    END IF
    IF(PRESENT(packed_dfnl))THEN
       p_dfnl=packed_dfnl
    ELSE
       p_dfnl=.FALSE.
    END IF
    IF(PRESENT(only_dfnl))THEN
       dfnl=only_dfnl
    ELSE
       dfnl=.FALSE.
    END IF
    IF(cntl%tfdist) THEN
       CALL rnlsmd(c0,nstate,ikind)
    ELSE      
       IF(.NOT.dfnl) CALL rnlsm1(c0,nstate,ikind,store_fnl_packed=p_fnl)
    ENDIF
    IF(tfor) CALL rnlsm2(c0,nstate,ikpt,ikind,store_dfnl_packed=p_dfnl)
    ! cp_group redistribution if necessary (defaults to yes)
    IF(redist)THEN
       IF(tfor.AND..NOT.dfnl) CALL cp_grp_redist_dfnl_fnl(.TRUE.,.TRUE.,nstate,ikind)
       IF(.NOT.tfor) CALL cp_grp_redist_dfnl_fnl(.TRUE.,.FALSE.,nstate,ikind)
       IF(tfor.AND.dfnl) CALL cp_grp_redist_dfnl_fnl(.FALSE.,.TRUE.,nstate,ikind)
    END IF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)

  END SUBROUTINE rnlsm
  ! ==================================================================
  SUBROUTINE give_scr_rnlsm(lrnlsm,tag,nstate,tfor)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrnlsm
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: tfor

    INTEGER                                  :: lrnlsm1, lrnlsm2

    CALL give_scr_rnlsm1(lrnlsm1,tag,nstate)
    IF (tfor) THEN
       CALL give_scr_rnlsm2(lrnlsm2,tag,nstate)
    ELSE
       lrnlsm2=0
    ENDIF
    lrnlsm=MAX(lrnlsm1,lrnlsm2)
    tag   ='MAX(LRNLSM1,LRNLSM2)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rnlsm
  ! ==================================================================


END MODULE rnlsm_utils
