C hello emacs, this is -*- fortran -*-
#ifdef EWATCUT
C arg
      real*8 RCUTP,RCUTP2
#endif
#ifdef EWALD
C args
      real*8 ALPHA
      real*8 DF3LJ,DF3EL
#ifdef HUM_VIR
      real*8 FELEC(NDIM*NATTOT),FELECI(MAXLAT)
#endif
C
#include "inc_latsum3.h"
#endif


#ifdef EWATCUT
#ifndef GAU_SHAPE
      RCUTP2 = ALPHA**2

#ifdef GAS_SHAPE
#ifdef GAU_TRUNC
      RCUTP2 = GAU_TRUNC*GAU_TRUNC
#else
      PRINT *,'GAS SHAPE: DEFINE GAU TRUNC'
      STOP
#endif
#endif

#else
      RCUTP2 = RCUTP**2
#endif
#endif
