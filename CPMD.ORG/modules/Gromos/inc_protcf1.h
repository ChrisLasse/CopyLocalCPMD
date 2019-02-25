C hello emacs, this is -*- fortran -*-
#ifdef TCF_READREF
                    IF (I.EQ.0) THEN
                      IF (LRDCNF) THEN
C read reference configuration
                        PRINT *
                        PRINT *,'READING REFERENCE COORDINATE FILE'
                        PRINT *,'UNIT 99'
                        PRINT *
                        CALL RDCONF (99,ICOORD,LFORM2,LRED2,.FALSE.,
     $                       COOTIT,NCOOTI,
     $                       NRP,XP,.FALSE.,DUM,LEOF2)
                        CALL PRTIT (' ',NCOOTI,COOTIT )
C translate reference configuration
                        CALL CENMAS (MAXNAT,1,NRP,0,0,0,3,3,XP,XP,
     $                       TMASP,
     $                       0,WMAS,WMASS,DUM,XCM,XCM,XCM,DUM,XCM,-1)
                        J3=0
                        DO 9025 J=1,NRP
                          DO 9024 M=1,3
                            J3=J3+1
                            XP(J3)=XP(J3)-XCM(M)
 9024                     CONTINUE
 9025                   CONTINUE
                      ELSE
C reference configuration is the first configuration
!$OMP parallel do private(J3)
                        DO J3=1,NRP3
                          XP(J3)=F(J3)
                        ENDDO
                      ENDIF
                      IROT=1
                    ENDIF
C fit current configuration to reference configuration
C if first step, do that only if the configuration was read from file
                    IF (I.NE.0 .OR. LRDCNF) THEN
                      CALL LSQSTR (NRP,W,XP,F,DUM,IROT,3)
                      IF (IROT.EQ.0) STOP
                      RMS = 0.0d0
                      SUMW = 0.0d0
                      DO II = 1,NRP
                        II3 = 3*II-2
                        RMS = RMS + W(II)*DSQRT(
     .                       + (XP(II3  )-F(II3  ))**2
     .                       + (XP(II3+1)-F(II3+1))**2
     .                       + (XP(II3+2)-F(II3+2))**2
     .                       )
                        SUMW = SUMW + W(II)
                      ENDDO
                      IF (RMS.GT.1.0d-10) THEN
                        RMS = RMS/SUMW
                      ELSE
                        RMS = 0.0d0
                      ENDIF
                    ELSE
C the first configuration is the reference configuration
                      RMS = 0.0d0
                    ENDIF
#endif
