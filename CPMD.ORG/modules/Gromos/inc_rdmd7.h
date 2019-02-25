C hello emacs, this is -*- fortran -*-
#ifdef READ_COORD
      IF (NRDCOR.NE.0) THEN
        IF (INIT.NE.INITCO) THEN
          PRINT '(A,I2)','NRDCOR.NE.0 REQUIRES INIT = ',INITCO
          STOP
        ENDIF
        NSTINC = NRDCOR * NSTCOR
        PRINT *
        PRINT *,'SUCCESSIVE CONFIGURATIONS WILL BE READ FROM FILE'
        PRINT *,'VELOCITIES WILL BE SET TO ZERO'
        PRINT '(A,I5,A,I5,A)','FILE WRITTEN WITH NTWX = ',
     $       NSTCOR,', ONE RECORD OVER ',NRDCOR,' WILL BE USED'
        PRINT *
      ELSE
        NSTINC = 1
      ENDIF
#endif
