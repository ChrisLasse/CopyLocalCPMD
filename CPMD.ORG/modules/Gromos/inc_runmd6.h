C hello emacs, this is -*- fortran -*-
#ifdef READ_COORD
Cmb - Revised on 4 Oct. 2005
         ELSE
C IF (NRDCOR.NE.0)
C
C set velocities to zero
           CALL mm_AZZERO(V,NATTOT) 
C set kinetic energies to zero
           ENER(IKSLU3) = 0.0d0
           ENER(IKSLU4) = 0.0d0
           ENER(IKSLV3) = 0.0d0
           ENER(IKSLV4) = 0.0d0
           ENER(IKTOT4) = 0.0d0
           ENER(IKTOT3) = 0.0d0
           ENER(IK4THD) = 0.0d0
         ENDIF
#endif
