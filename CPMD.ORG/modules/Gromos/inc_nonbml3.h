C hello emacs, this is -*- fortran -*-
#ifdef EWALD
#ifdef EWATCUT
                  IF (RIJ2 .LE. RCUTP2) THEN
#endif

#define REAL_SPACE
#include "inc_latsum9.h"
#undef REAL_SPACE

#ifdef HAT_SHAPE
#ifdef IONDERIV_READY
#ifdef IONDERIV
                    DERIR = DERIR + 
     .                DSQRT(138.93569d0)*CGAJ(IPQ)*AMULT(IPQ)
     .                            *(RIJINV+FACG3*RIJ2/RIJINV
     .                            + FACG2*RIJ2 + FACG0 )
#endif  
#endif
#endif

#ifdef EWATCUT
                  ELSE
#ifdef GAU_SHAPE
C since here rcutp2 = alpha2, dont do atomic truncation of the LJ 
C for the hat and parabolic functions !
                    VLJA4D = 0.0d0
#endif
                    VELA4D = 0.0d0
                    DFA4D = 0.0d0
                  ENDIF
#endif
#endif


#ifdef EWATCUT
#ifndef EWALD
                  IF (RIJ2 .GT. RCUTP2) THEN
#ifdef GAU_SHAPE
C since here rcutp2 = alpha2, dont do atomic truncation of the LJ 
C for the hat and parabolic functions !
                    VLJA4D = 0.0d0
#endif
                    VELA4D = 0.0d0
                    DFA4D = 0.0d0
                  ENDIF
#endif
#endif
