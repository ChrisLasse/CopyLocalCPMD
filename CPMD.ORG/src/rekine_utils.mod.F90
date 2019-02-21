MODULE rekine_utils
  USE cp_grp_utils,                    ONLY: cp_grp_split_states  
  USE dotp_utils,                      ONLY: dotp
  USE geq0mod,                         ONLY: geq0
  USE harm,                            ONLY: xmu
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             ncpw,parap
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rekine

CONTAINS

  ! ==================================================================
  SUBROUTINE rekine(cm,nstate,ekinc)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: cm(ncpw%ngw,nstate)
    REAL(real_8)                             :: ekinc

    INTEGER                                  :: i, ig , isub
    REAL(real_8)                             :: ax, bx, pf
    CHARACTER(*), PARAMETER                  :: procedureN = 'rekine'        
    INTEGER                                  :: first_state,last_state

! ==--------------------------------------------------------------==
! ==  COMPUTE FICTITIOUS KINETIC ENERGY OF THE ELECTRONS          ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    
    ! split states between cp groups
    call cp_grp_split_states(nstate,first_state=first_state,last_state=last_state)
    
    ekinc=0._real_8
    IF (cntl%tmass) THEN
       !$omp parallel do private(i,ig,pf,ax,bx) reduction(+:ekinc)
       DO i=first_state,last_state
          DO ig=1,ncpw%ngw
             pf=2.0_real_8*xmu(ig)
             ax=REAL(cm(ig,i))
             bx=AIMAG(cm(ig,i))
             ekinc=ekinc+pf*(ax*ax+bx*bx)
          ENDDO
          IF (geq0) ekinc=ekinc-xmu(1)*REAL(cm(1,i))*REAL(cm(1,i))
       ENDDO
    ELSE
       !$omp parallel do private(i) reduction(+:ekinc)
       DO i=first_state,last_state
          ekinc=ekinc+dotp(ncpw%ngw,cm(:,i),cm(:,i))
       ENDDO
       ekinc=ekinc*cntr%emass
    ENDIF
    CALL mp_sum(ekinc,parai%cp_grp)
    CALL tihalt(procedureN,isub)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rekine
  ! ==================================================================

END MODULE rekine_utils
