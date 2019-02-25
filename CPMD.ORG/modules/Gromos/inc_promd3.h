C hello emacs, this is -*- fortran -*-
#ifdef READ_COORD
      IF (NRDCOR .EQ. 0) THEN
        CALL GTCOOR(PRGSTR,NATTOT,MXEWRT,DLSUM,DMSUM)
      ELSE        
        CALL DETOPN('RDCO',IURDCO,LFORMR)
        IF (IURDCO .LT.0) THEN
          PRINT *,'Error opening position trajectory file'
          STOP
        ENDIF
        CALL RDBHDR(IURDCO,LFORMR,IDBLKR,NRECR)
        IF (IDBLKR .NE. ITITID ) THEN
          PRINT *,'TITLE block expected unit ',IURDCO
          STOP
        ENDIF
        CALL RDTIT(IURDCO,LFORMR,NRECR,MAXTIT,RDCTIT,NRDCTI)
        CALL PRTIT('INPUT COORDINATE TRAJECTORY',NRDCTI,RDCTIT)
        PRINT *
        PRINT *,'VELOCITIES WILL BE SET TO ZERO'
        PRINT *
      ENDIF   
#endif
