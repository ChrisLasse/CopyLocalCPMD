C hello emacs, this is -*- fortran -*-
Cmb - inc_latsum9.h: Revised on 26 April 2005
#ifdef RECI_SPACE_VAR_DECL
      real*8 ALPHAC,ALPHAK,FKCUT2
C and for TSC
      real*8 SALPHAK2,SALPHAK3
#endif


C********** SPHERICAL HAT 
#ifdef HAT_SHAPE

#ifdef RECI_SPACE_VAR_SET
      ALPHAC = 1.2d+1/ALPHA**4
#endif
#ifdef RECI_SPACE
      ALPHAK = ALPHA*DSQRT(K2)
      GAMHAT = ALPHAC/(K2*K2)*(2.0d0-2.0d0*DCOS(ALPHAK)
     .                         -ALPHAK*DSIN(ALPHAK))
#endif
#ifdef REAL_SPACE_VAR
      real*8 FACG0,FACG2,FACG3,FACGF1,FACGF2,ALPHA2
C
      FACG0  = -2.0d0/ALPHA
      FACG2  =  2.0d0/ALPHA**3
      FACG3  = -1.0d0/ALPHA**4
      FACGF1 = -2.0d0*FACG2
      FACGF2 = -3.0d0*FACG3
      ALPHA2 = ALPHA**2
#endif
#ifdef REAL_SPACE
                   VELA4D=QIQJA*(RIJINV 
     .                          + FACG3*RIJ2/RIJINV 
     .                          + FACG2*RIJ2+FACG0) 
                   DFA4D =QIQJA*RIJINV
     .                   *(RIJIN2+FACGF1/RIJINV
     .                   + FACGF2*RIJ2)+RIJIN2
     .            *(6.0d0*(CRA1+CRA12)* RIJIN6)
#endif
#ifdef EXCL_CORR
            IF (RIJ2.LE.ALPHA2) THEN
C exclusion distance shorter than alpha, add eta(r)-1/r
              VELA4D=QIQJA
     .             *(FACG3*RIJ2/RIJINV+FACG2*RIJ2+FACG0 )
              DFA4D =QIQJA*RIJINV
     .             *(FACGF1/RIJINV+FACGF2*RIJ2)
            ELSE
C exclusion distance larger than alpha, add -1/r
              VELA4D = -QIQJA*RIJINV
              DFA4D  = VELA4D*RIJIN2
cmb              DFA4D  = -QIQJA*RIJINV*RIJIN2
            ENDIF
#endif
#ifdef SELF_TERM
        A1=-4.0d0*(4.0d0*DATAN(1.0d0))*ALPHA**2/(1.5d+1*VOL)
        A3=-2.0d0/ALPHA
#endif

#endif


C********** SPHERICAL PARABOLA (CENTERED)
#ifdef PAC_SHAPE

#ifdef RECI_SPACE_VAR_SET
      ALPHAC = 1.5d+1/ALPHA**5
#endif
#ifdef RECI_SPACE
      ALPHAK = ALPHA*DSQRT(K2)
      GAMHAT = ALPHAC/(K2*K2*DSQRT(K2))
     .       *(-3.0d0*ALPHAK*DCOS(ALPHAK)
     .       +(3.0d0-ALPHAK**2)*DSIN(ALPHAK) )
#endif
#ifdef REAL_SPACE_VAR
      real*8 FACG0,FACG3,FACG4,FACGF2,FACGF3,ALPHA2
C
      FACG0 = -1.5d+1/(8.0d0*ALPHA)
      FACG2 =  5.0d0/(4.0d0*ALPHA**3)
      FACG4 = -3.0d0/(8.0d0*ALPHA**5)
      FACGF1 = -2.0d0*FACG2
      FACGF3 = -4.0d0*FACG4      
cmb      FACGF1 = -5.0d0/(2.0d0*ALPHA**3)
cmb      FACGF3 =  3.0d0/(2.0d0*ALPHA**5)
      ALPHA2 = ALPHA**2
#endif
#ifdef REAL_SPACE
                    VELA4D = QIQJA*(RIJINV 
     .                     + FACG4*RIJ2*RIJ2
     .                     + FACG2*RIJ2+FACG0)
                    DFA4D = QIQJA*RIJINV
     .                    * (RIJIN2+FACGF1/RIJINV
     .                    + FACGF3*RIJ2/RIJINV)
     .                    + RIJIN2 
     .                    *(6.0d0*(CRA1+CRA12)*RIJIN6)
#endif
#ifdef EXCL_CORR
            IF (RIJ2.LE.ALPHA2) THEN
C exclusion distance shorter than alpha, add eta(r)-1/r
              VELA4D=QIQJA*((FACG4*RIJ2+FACG2)*RIJ2+FACG0)
              DFA4D =QIQJA*(FACGF1+FACGF3*RIJ2)
cmb              DFA4D =QIQJA*RIJINV*(FACGF1/RIJINV+FACGF3*RIJ2/RIJINV)
            ELSE
C exclusion distance larger than alpha, add -1/r
              VELA4D = -QIQJA*RIJINV
              DFA4D  = VELA4D*RIJIN2
cmb              DFA4D  = - QIQJA*RIJINV*RIJIN2
            ENDIF
#endif
#ifdef SELF_TERM
        A1=-2.0d0*(4.0d0*DATAN(1.0d0))*ALPHA**2/(7.0d0*VOL)
        A3=-1.875d0/ALPHA
cmb        A3=-1.5d+1/(8.0d0*ALPHA)
#endif

#endif


C********** SPHERICAL PARABOLA (BERENDSEN)
#ifdef PARAB_SHAPE

#ifdef RECI_SPACE_VAR_SET
      ALPHAC = 2.0d+1/ALPHA**5
#endif
#ifdef RECI_SPACE
      ALPHAK = ALPHA*DSQRT(K2)
      GAMHAT = ALPHAC/(K2*K2*DSQRT(K2))
     .       *(-2.0d0*ALPHAK*(1.0d0+2.0d0*DCOS(ALPHAK))
     .         +(6.0d0-ALPHAK**2)*DSIN(ALPHAK))
#endif
#ifdef REAL_SPACE_VAR
      real*8 FACG0,FACG3,FACG4,FACGF2,FACGF3,ALPHA2
C
      FACG0 = -5.0d0/(3.0d0*ALPHA)
      FACG3 =  5.0d0/(3.0d0*ALPHA**4)
      FACG4 = -1.0d0/ALPHA**5
      FACGF2 = -3.0d0*FACG3
      FACGF3 = -4.0d0*FACG4
      ALPHA2 = ALPHA**2
#endif
#ifdef REAL_SPACE
                    VELA4D = QIQJA*(RIJINV 
     .                    + RIJ2*(FACG4*RIJ2
     .                    + FACG3/RIJINV) 
     .                    + FACG0)
                    DFA4D  = QIQJA*RIJINV
     .                    * (RIJIN2+FACGF2*RIJ2
     .                    + FACGF3*RIJ2/RIJINV)
     .                    + RIJIN2 
     .              *(6.0d0*(CRA1+CRA12)*RIJIN6)
#endif
#ifdef EXCL_CORR
            IF (RIJ2.LE.ALPHA2) THEN
C exclusion distance shorter than alpha, add eta(r)-1/r
              VELA4D=QIQJA*((FACG4*RIJ2+FACG3/RIJINV)*RIJ2+FACG0)
              DFA4D =QIQJA*RIJINV*(FACGF2*RIJ2+FACGF3*RIJ2/RIJINV)
            ELSE
C exclusion distance larger than alpha, add -1/r
              VELA4D = -QIQJA*RIJINV
              DFA4D  = VELA4D*RIJIN2
cmb              DFA4D  = - QIQJA*RIJINV*RIJIN2
            ENDIF
#endif
#ifdef SELF_TERM
        A1=-2.0d+1*(4.0d0*DATAN(1.0d0))*ALPHA**2/(6.3d+1*VOL)
        A3=-5.0d0/(3.0d0*ALPHA)
#endif

#endif


C********** SPHERICAL SHELL
#ifdef SSH_SHAPE

#ifdef RECI_SPACE
                      ALPHAK = ALPHA*DSQRT(K2)
                      GAMHAT = DSIN(ALPHAK)/ALPHAK
#endif
#ifdef REAL_SPACE_VAR
      real*8 FACG0,ALPHA2
C
      FACG0 = -1.0d0/ALPHA
      ALPHA2 = ALPHA**2
#endif
#ifdef REAL_SPACE
                   VELA4D = QIQJA*(RIJINV+FACG0) 
                   DFA4D  = QIQJA*RIJINV*RIJIN2
     .                    + RIJIN2 
     .              *(6.0d0*(CRA1+CRA12)*RIJIN6)
#endif
#ifdef EXCL_CORR
            IF (RIJ2.LE.ALPHA2) THEN
C exclusion distance shorter than alpha, add eta(r)-1/r
              VELA4D = QIQJA*FACG0
              DFA4D  = 0.0d0
            ELSE
C exclusion distance larger than alpha, add -1/r
              VELA4D = -QIQJA*RIJINV
              DFA4D  = VELA4D*RIJIN2
cmb              DFA4D  = - QIQJA*RIJINV*RIJIN2
            ENDIF
#endif
#ifdef SELF_TERM
        A1=-2.0d0*(4.0d0*DATAN(1.0d0))*ALPHA**2/(3.0d0*VOL)
        A3=-1.0d0/ALPHA
#endif

#endif


C********** SPHERICAL CONSTANT
#ifdef SCO_SHAPE

#ifdef RECI_SPACE_VAR_SET
      ALPPAC = 3.0d0/ALPHA**3
#endif
#ifdef RECI_SPACE
      ALPHAK = ALPHA*DSQRT(K2)
      GAMHAT = ALPPAC/(K2*DSQRT(K2))
     .       *(DSIN(ALPHAK)-ALPHAK*DCOS(ALPHAK))
#endif
#ifdef REAL_SPACE_VAR
      real*8 FACG0, FACG2, FACGF1,ALPHA2
C
      FACG0 = -1.5d0/ALPHA
      FACG2 =  0.5d0/(ALPHA**3)
      FACGF1 = -2.0d0*FACG2
      ALPHA2 = ALPHA**2
#endif
#ifdef REAL_SPACE
                    VELA4D = QIQJA*(RIJINV 
     .                     + FACG2*RIJ2+FACG0)
                    DFA4D  = QIQJA*RIJINV
     .                     *(RIJIN2+FACGF1/RIJINV)
     .                     + RIJIN2 
     .                *(6.0d0*(CRA1+CRA12)*RIJIN6)
#endif
#ifdef EXCL_CORR
            IF (RIJ2.LE.ALPHA2) THEN
C exclusion distance shorter than alpha, add eta(r)-1/r
              VELA4D = QIQJA*(FACG2*RIJ2+FACG0)
              DFA4D  = QIQJA*RIJINV*(FACGF1/RIJINV)
            ELSE
C exclusion distance larger than alpha, add -1/r
              VELA4D = -QIQJA*RIJINV
              DFA4D  = VELA4D*RIJIN2
cmb              DFA4D  = - QIQJA*RIJINV*RIJIN2
            ENDIF
#endif
#ifdef SELF_TERM
        A1=-2.0d0*(4.0d0*DATAN(1.0d0))*ALPHA**2/(5.0d0*VOL)
        A3=-1.5d0/ALPHA
#endif

#endif


C********** SPHERICAL INVERSE WITH OFFSET
#ifdef SIO_SHAPE

#ifdef RECI_SPACE_VAR_SET
      ALPHAC = 6.0d0/ALPHA**3
#endif
#ifdef RECI_SPACE
      ALPHAK = ALPHA*DSQRT(K2)
      GAMHAT = ALPHAC/(K2*DSQRT(K2))*(ALPHAK-DSIN(ALPHAK))
#endif
#ifdef REAL_SPACE_VAR
      real*8 FACG0, FACG2, FACGF1,ALPHA2
C
      FACG0 = -3.0d0/(ALPHA)
      FACG1 =  3.0d0/(ALPHA**2)
      FACG2 = -1.0d0/(ALPHA**3)
      FACGF0 = -FACG1
      FACGF1 = -2.0d0*FACG2
cmb      FACGF0 = -3.0d0/(ALPHA**2)
cmb      FACGF1 =  2.0d0/(ALPHA**3)
      ALPHA2 = ALPHA**2
#endif
#ifdef REAL_SPACE
                    VELA4D = QIQJA*(RIJINV
     .                     + FACG2*RIJ2
     .                     + FACG1/RIJINV
     .                     + FACG0)
                    DFA4D  = QIQJA*RIJINV*(RIJIN2
     .                     + FACGF1/RIJINV
     .                     + FACGF0 )
     .                     + RIJIN2 
     .                *(6.0d0*(CRA1+CRA12)*RIJIN6)
#endif
#ifdef EXCL_CORR
            IF (RIJ2.LE.ALPHA2) THEN
C exclusion distance shorter than alpha, add eta(r)-1/r
              VELA4D = QIQJA
     .             * (FACG2*RIJ2+FACG1/RIJINV+FACG0)
              DFA4D  = QIQJA * RIJINV
     .             * (FACGF1/RIJINV+FACGF0)
            ELSE
C exclusion distance larger than alpha, add -1/r
              VELA4D = -QIQJA*RIJINV
              DFA4D  = VELA4D*RIJIN2
cmb              DFA4D  = - QIQJA*RIJINV*RIJIN2
            ENDIF
#endif
#ifdef SELF_TERM
        A1=-(4.0d0*ATAN(1.0d0))*ALPHA**2/(5.0d0*VOL)
        A3=-3.0d0/ALPHA
#endif

#endif


C********** SPHERICAL INVERSE WITH *NO* OFFSET
#ifdef SIN_SHAPE

#ifdef RECI_SPACE_VAR_SET
      ALPHAC = 2.0d0/ALPHA**2
#endif
#ifdef RECI_SPACE
      ALPHAK = ALPHA*DSQRT(K2)
      GAMHAT = ALPHAC/K2*(1.0d0-DCOS(ALPHAK))
#endif
#ifdef REAL_SPACE_VAR
      real*8 FACG0, FACG2, FACGF1,ALPHA2
C
      FACG0 = -2.0d0/(ALPHA)
      FACG1 =  1.0d0/(ALPHA**2)
      FACGF0 = -FACG1
cmb      FACGF0 = -1.0d0/(ALPHA**2)
      ALPHA2 = ALPHA**2
#endif
#ifdef REAL_SPACE
                    VELA4D = QIQJA*(RIJINV
     .                     + FACG1/RIJINV+FACG0)
                    DFA4D  = QIQJA*RIJINV*(RIJIN2
     .                     +FACGF0)+RIJIN2
     .               *(6.0d0*(CRA1+CRA12)*RIJIN6)
#endif
#ifdef EXCL_CORR
            IF (RIJ2.LE.ALPHA2) THEN
C exclusion distance shorter than alpha, add eta(r)-1/r
              VELA4D = QIQJA*(FACG1/RIJINV+FACG0)
              DFA4D  = QIQJA*RIJINV*FACGF0
            ELSE
C exclusion distance larger than alpha, add -1/r
              VELA4D = -QIQJA*RIJINV
              DFA4D  = VELA4D*RIJIN2
cmb              DFA4D  = - QIQJA*RIJINV*RIJIN2
            ENDIF
#endif
#ifdef SELF_TERM
        A1=-(4.0d0*DATAN(1.0d0))*ALPHA**2/(3.0d0*VOL)
        A3=-2.0d0/ALPHA
#endif

#endif



C********** SPHERICAL GAUSSIAN
#ifdef GAU_SHAPE

#ifdef RECI_SPACE_VAR_SET
      ALPHAC = -0.25d0/ALPHA**2
C
C     PREVENT UNDERFLOW IN THE EXP FUNCTION (FOR RHAT) BY TRUNCATING 
C     THE SQUARE OF THE RECIPROCAL VECTORS K AT A CUT-OFF FKCUT2. 
C     HERE FOR AN ACCURACY OF 10E-20
C
      FKCUT2 = -2.0d+1*DLOG(1.0d+1)/ALPHAC
#endif
#ifdef RECI_SPACE
      GAMHAT =  EXP(ALPHAC*K2)
#endif
#ifdef REAL_SPACE_VAR
      real*8 PEW,AEW1,AEW2,AEW3,AEW4,AEW5
      PARAMETER ( AEW1 = 0.254829592d0, AEW2 = -0.284496736d0)
      PARAMETER ( AEW3 = 1.421413741d0, AEW4 = -1.453152027d0)
      PARAMETER ( AEW5 = 1.061405429d0, PEW  =  0.3275911d0  )
      real*8 ALPHAR,TMPEW,EXPAR2,ERFC,FACERW
C
      FACERW = 2.0d0*ALPHA/DSQRT(4.0d0*DATAN(1.0d0))
#endif
#ifdef REAL_SPACE
                   ALPHAR = ALPHA/RIJINV
                   TMPEW  = 1.0d0/(1.0d0+PEW*ALPHAR)
                   EXPAR2 = DEXP(-ALPHAR**2)
                   ERFC   = ((((AEW5*TMPEW+AEW4)*TMPEW+AEW3)
     .                    *TMPEW+AEW2)*TMPEW+AEW1)*TMPEW*EXPAR2
                   VELA4D = QIQJA*RIJINV * ERFC
                   DFA4D  = QIQJA*RIJINV * RIJIN2
     .                    * (FACERW/RIJINV*EXPAR2+ERFC)
     .                    + RIJIN2*(6.0d0*(CRA1+CRA12)*RIJIN6)
#endif
#ifdef EXCL_CORR
C add eta(r)-1/r = [erfc(alpha*r)-1]/r = -erf(alpha*r)/r
            ALPHAR = ALPHA / RIJINV
            TMPEW  = 1.0d0/(1.0d0+PEW*ALPHAR)
            EXPAR2 = DEXP(-ALPHAR**2)
            ERFC   = ((((AEW5*TMPEW+AEW4)*TMPEW+AEW3)
     .             *TMPEW+AEW2)*TMPEW+AEW1)*TMPEW*EXPAR2
            VELA4D =  QIQJA*RIJINV*(ERFC-1.0d0)
            DFA4D  =  QIQJA*RIJINV*RIJIN2 
     .             * ((FACERW/RIJINV*EXPAR2+ERFC)-1.0d0)
#endif
#ifdef SELF_TERM
CCC        GAMSLF = - ALPHA/SQRT(3.0D0*ATAN(1.0D0))
C that was wrong !
        A1=-(4.0d0*ATAN(1.0d0))/(VOL*ALPHA**2)
        A3=-2.0d0*ALPHA/DSQRT(4.0d0*DATAN(1.0d0))
#endif

#endif


C********** SPHERICAL GAUSSIAN - we use this to enforce
C strict truncation of electrostatics at given distance 
C -DGAU_TRUNC xxx
C otherwise identical to GAU_SHAPE
#ifdef GAS_SHAPE

#ifdef RECI_SPACE_VAR_SET
      ALPHAC = -0.25d0/ALPHA**2
C
C     PREVENT UNDERFLOW IN THE EXP FUNCTION (FOR RHAT) BY TRUNCATING 
C     THE SQUARE OF THE RECIPROCAL VECTORS K AT A CUT-OFF FKCUT2. 
C     HERE FOR AN ACCURACY OF 10E-20
C
      FKCUT2 = - 2.0d+1*DLOG(1.0d+1)/ALPHAC
#endif
#ifdef RECI_SPACE
      GAMHAT = DEXP(ALPHAC*K2)
#endif
#ifdef REAL_SPACE_VAR
      real*8 PEW,AEW1,AEW2,AEW3,AEW4,AEW5
      PARAMETER ( AEW1 = 0.254829592d0, AEW2 = -0.284496736d0)
      PARAMETER ( AEW3 = 1.421413741d0, AEW4 = -1.453152027d0)
      PARAMETER ( AEW5 = 1.061405429d0, PEW  =  0.3275911d0  )
      real*8 ALPHAR,TMPEW,EXPAR2,ERFC,FACERW
C
      FACERW = 2.0d0*ALPHA/DSQRT(4.0d0*DATAN(1.0d0))
#endif
#ifdef REAL_SPACE
                    ALPHAR = ALPHA/RIJINV
                    TMPEW  = 1.0d0/(1.0d0+PEW*ALPHAR)
                    EXPAR2 = DEXP(-ALPHAR**2)
                    ERFC   = ((((AEW5*TMPEW+AEW4)*TMPEW+AEW3)
     .                     *TMPEW+AEW2)*TMPEW+AEW1)*TMPEW*EXPAR2
                    VELA4D =  QIQJA*RIJINV*ERFC
                    DFA4D  =  QIQJA*RIJINV*RIJIN2
     .                     * (FACERW/RIJINV*EXPAR2+ERFC)
     .                     + RIJIN2*(6.0d0*(CRA1+CRA12)*RIJIN6)
#endif
#ifdef EXCL_CORR
C add eta(r)-1/r = [erfc(alpha*r)-1]/r = -erf(alpha*r)/r
            ALPHAR = ALPHA / RIJINV
            TMPEW  = 1.0d0/(1.0d0+PEW*ALPHAR)
            EXPAR2 = DEXP(-ALPHAR**2)
            ERFC   = ((((AEW5*TMPEW+AEW4)*TMPEW+AEW3)
     .             *TMPEW+AEW2)*TMPEW+AEW1)*TMPEW*EXPAR2
            VELA4D =  QIQJA*RIJINV*(ERFC-1.0d0)
            DFA4D  =  QIQJA*RIJINV*RIJIN2 
     .            * ((FACERW/RIJINV*EXPAR2+ERFC)-1.0d0)
#endif
#ifdef SELF_TERM
CCC        GAMSLF = - ALPHA/SQRT(3.0D0*ATAN(1.0D0))
C that was wrong !
        A1=-(4.0d0*DATAN(1.0d0))/(VOL*ALPHA**2)
        A3=-2.0d0*ALPHA/DSQRT(4.0d0*ATAN(1.0d0))
#endif

#endif


C********** TSC SHAPE
#ifdef TSC_SHAPE

#ifdef RECI_SPACE_VAR_SET
      ALPHAC = 243.0d0/ALPHA**5
#endif
#ifdef RECI_SPACE
                      ALPHAK = ALPHA*DSQRT(K2)
                      SALPHAK3 = SIN(ALPHAK/3.0d0)                      
                      GAMHAT = ALPHAC/(K2*K2*DSQRT(K2))
     .                       *SALPHAK3**2*(3.0d0*SALPHAK3
     .                       -ALPHAK*DCOS(ALPHAK/3.0d0))
#endif
#ifdef REAL_SPACE_VAR
      real*8 FACGIP, FACG0M, FACG0P, 
     $     FACG2M, FACG2P, FACG3P, FACG4M, FACG4P  
      real*8 FACGFIP, FACGF1M, FACGF1P, 
     $     FACGF2P, FACGF3M, FACGF3P
      real*8 ALPHA32
C
      FACG0M = -39.0d0/(16.0d0*ALPHA)
      FACG2M =  27.0d0/(8.0d0*ALPHA**3)
      FACG4M = -243.0d0/(80.0d0*ALPHA**5)
C
      FACGIP =  81.0d0/80.0d0
      FACG0P = -81.0d0/(32.0d0*ALPHA)
      FACG2P =  81.0d0/(16.0d0*ALPHA**3)
      FACG3P = -81.0d0/(16.0d0*ALPHA**4)
      FACG4P =  243.0d0/(160.0d0*ALPHA**5)
C
      FACGF1M = -2.0d0*FACG2M
      FACGF3M = -4.0d0*FACG4M
C
      FACGFIP = FACGIP
      FACGF1P = -2.0d0*FACG2P
      FACGF2P = -3.0d0*FACG3P
      FACGF3P = -4.0d0*FACG4P
C
      ALPHA2 = ALPHA**2
      ALPHA32 = (ALPHA/3.0d0)**2
#endif
#ifdef REAL_SPACE
      IF (RIJ2.LE.ALPHA32) THEN
                    VELA4D = QIQJA*(RIJINV 
     .                     + RIJ2*(FACG4M*RIJ2
     .                     + FACG2M) +FACG0M )
                    DFA4D  = QIQJA*RIJINV*(RIJIN2
     .                     + FACGF3M*RIJ2/RIJINV
     .                     + FACGF1M/RIJINV)
     .                     + RIJIN2 
     .               *(6.0d0*(CRA1+CRA12)*RIJIN6)
      ELSE
                    VELA4D = QIQJA*(FACGIP*RIJINV
     .                     + RIJ2*(FACG4P*RIJ2
     .                     + FACG3P/RIJINV) 
     .                     + FACG2P*RIJ2+FACG0P)
                    DFA4D  = QIQJA*RIJINV
     .                     *(FACGFIP*RIJIN2
     .                     + RIJ2*(FACGF3P/RIJINV
     .                     + FACGF2P) 
     .                     + FACGF1P/RIJINV )
     .                     + RIJIN2 
     .               *(6.0d0*(CRA1+CRA12)*RIJIN6)
      ENDIF
#endif
#ifdef EXCL_CORR
            IF (RIJ2.LE.ALPHA2) THEN
C exclusion distance shorter than alpha, add eta(r)-1/r
              IF (RIJ2.LE.ALPHA32) THEN
                VELA4D=QIQJA*((FACG4M*RIJ2+FACG2M)*RIJ2+FACG0M)
                DFA4D =QIQJA*(FACGF3M*RIJ2+FACGF1M)
cmb                DFA4D =QIQJA*RIJINV*(FACGF3M*RIJ2/RIJINV
cmb     .                +FACGF1M/RIJINV)
              ELSE
                VELA4D=QIQJA*((FACGIP-1.0d0)*RIJINV+FACG4P*RIJ2*RIJ2
     .                +FACG3P*RIJ2/RIJINV+FACG2P*RIJ2+FACG0P)
                DFA4D =QIQJA*RIJINV*((FACGFIP-1.0d0)*RIJIN2
     .                +FACGF3P*RIJ2/RIJINV+FACGF2P*RIJ2 
     .                +FACGF1P/RIJINV)
              ENDIF
            ELSE
C exclusion distance larger than alpha, add -1/r
              VELA4D = -QIQJA*RIJINV
              DFA4D  = VELA4D*RIJIN2
cmb              DFA4D  = - QIQJA*RIJINV*RIJIN2
            ENDIF
#endif
#ifdef SELF_TERM
        A1=-26.0d0*(4.0d0*DATAN(1.0d0))*ALPHA**2/(135.0d0*VOL)
        A3=-39.0d0/(16.0d0*ALPHA)
#endif

#endif
