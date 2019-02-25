C hello emacs, this is -*- fortran -*-
#ifdef TCF_READREF
      IF (NTTR.EQ.3) THEN
        LRDCNF = .TRUE.
        NTTR = 1
      ELSEIF (NTTR.EQ.-3) THEN
        LRDCNF = .TRUE.
        NTTR = -1
      ELSE
        LRDCNF = .FALSE.
      ENDIF
#endif
