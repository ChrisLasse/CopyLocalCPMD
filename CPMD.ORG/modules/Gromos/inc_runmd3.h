C hello emacs, this is -*- fortran -*-
#ifdef EWALD
         IF (NGHTUP.LT.0) THEN
           NGHTUP = -NGHTUP
           LGHOPT = .TRUE.
           LGHWRT = .TRUE.
C does not really matter
           LGHUPD = .FALSE.
         ELSE
           LGHOPT = .FALSE.
           LGHWRT = .FALSE.
           LGHUPD = (MOD(NSTEP+1,NGHTUP) .EQ. 0)
         ENDIF
#endif
#ifdef READ_COORD
C while loop 
      IF (NRDCOR.NE.0) THEN
        IF (NSTEP.EQ.0) THEN
C always read the first configuration
          NRECRD = 1
        ELSE
          NRECRD = NRDCOR 
        ENDIF
        CALL RDBHDR(IURDCO,LFORMR,IDBLKR,NRECR)
        LGTCRD = .FALSE.
        LGTTIM = .FALSE.
 8000   IF (IDBLKR .NE. IDERR) THEN        
          IF (IDBLKR .EQ. IPORID) THEN
            CALL RDXRED(IURDCO,LFORMR,.FALSE.,NRECR,NATTOT,NDIM,X)
            LGTCRD = .TRUE.
          ELSEIF (IDBLKR .EQ. ITIMID) THEN
            CALL RDTIME(IURDCO,LFORMR,NRECR,NSTEPR,TIMER)
            LGTTIM = .TRUE.
          ELSE
            CALL SKPBLK(IURDCO,LFORMR,NRECR)
          ENDIF
          IF ( LGTCRD .AND. LGTTIM) THEN
            LGTCRD = .FALSE.
            LGTTIM = .FALSE.
            NRECRD = NRECRD - 1
          ENDIF
          IF ( NRECRD .NE. 0 ) THEN
            CALL RDBHDR(IURDCO,LFORMR,IDBLKR,NRECR)
            GOTO 8000
          ENDIF
        ELSE
          PRINT '(A,I3)',
     $         'ATTEMPTING TO READ CONFIGURATIONS PAST EOF ON UNIT ',
     $         IURDCO
          STOP
        ENDIF
        IF (NSTEP+NSTINC.GT.NSTLIM-1) THEN
          CALL CLSFIL(IURDCO)
          PRINT *,'ALL POSSIBLE CONFIGURATIONS READ ON UNIT ',IURDCO
        ENDIF        
        IF ( NSTEPR .NE. NSTEP .OR. 
     $       DABS(TIMER-TIME).GT.1.0D-4 ) THEN
          PRINT *,'TIME AND STEP READ FROM TRAJECTORY ARE INCONSISTENT'
          PRINT *,'WITH LOCAL VALUE'
          PRINT '(A,I6,A,F10.3)',
     $         'READ STEP ',NSTEPR,', READ TIME ',TIMER
          PRINT '(A,I6,A,F10.3)',
     $         'CURRENT STEP ',NSTEP,', CURRENT TIME ',TIME
CCC          STOP
        ENDIF
      ENDIF
#endif
