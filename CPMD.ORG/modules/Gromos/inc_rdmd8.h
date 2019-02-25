C hello emacs, this is -*- fortran -*-
#ifdef EWALD
C incompatible with RF
#endif
#ifdef READ_COORD
      IF (NRDCOR.NE.0) THEN
        IF (NSNB.NE.1) THEN
          PRINT *,'NRDCOR.NE.0 BUT NSNB.NE.1'
          STOP
        ENDIF
        IF (NSTCOR.LE.0) THEN
          PRINT  *,'NSTCOR.LE.0'
          STOP
        ENDIF
      ENDIF
#endif
