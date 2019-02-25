C hello emacs, this is -*- fortran -*-
#ifdef WRITEFORCES
      IF (NTWFOR .LT. 0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTWFOR'
         PRINT FMNII,VSTR,NTWFOR
         PRINT FMGEI,VSTR
         STOP
      ENDIF
#endif
