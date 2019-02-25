C hello emacs, this is -*- fortran -*-
#define RIJ2        VRIJ2(IPQ)
#define RIJINV      (1.0/VRIJ(IPQ))
C because VRINV(IPQ) is zero if amult=0
#define RIJIN2      VRINV2(IPQ)
#define RIJIN6      VRINV6(IPQ)
#define VELA4D      VEELJ(IPQ)
#define VLJA4D      VLJJ(IPQ)
#define DFA4D       DF3(IPQ)
#define CRA1        2.0E0*VRINV6(IPQ)*VC12A(IPQ)-VC6A(IPQ)
#define CRA12       0.0
#include "inc_nonbml3.h"
#undef RIJ2
#undef RIJINV
#undef RIJIN2
#undef RIJIN6
#undef VELA4D
#undef VLJA4D
#undef DFA4D
#undef CRA1
#undef CRA12
