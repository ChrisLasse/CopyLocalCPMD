C hello emacs, this is -*- fortran -*-
#ifdef WRITE_FORCES
C write forces if required
         IF (NTWFOR .NE. 0) THEN
           IF (MOD(NSTEP,NTWFOR) .EQ. 0) THEN
             CALL WRTIME(IUTRJF,LFORM,NSTEP,TIME)
             IF (NTWFOR .GT. 0) THEN
               CALL WRFRED(IUTRJF,LFORM,LWR4,NATTOT,NDIM,F)
             ELSE
C only write solute coords
               CALL WRFRED(IUTRJF,LFORM,LWR4,NRPT,NDIM,F)
             ENDIF
C             IF (NTB .NE. NTBVAC) THEN
C               CALL WRBOX(IUTRJF,LFORM,BOX)
C             ENDIF
           ENDIF
         ENDIF
#endif
