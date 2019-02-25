C hello emacs, this is -*- fortran -*-
#ifdef WRITE_FORCES
C forces
      IF (NTWFOR .NE. 0) THEN
        CALL OPNFIL('TRJF',FRMSTR,STATST,IUTRJF)
        IF (IUTRJF .LT. 0)THEN
          PRINT *,PRGSTR,
     $         ' :failed to open force trajectory file'
          STOP
        ENDIF
        CALL WRTIT(IUTRJF,LFORM,MDTLNS,MDTITL)
        CALL WRFMT(IUTRJF,LFORM)
      ENDIF   
#endif
