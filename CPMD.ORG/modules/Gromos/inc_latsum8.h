C hello emacs, this is -*- fortran -*-
C periodic copy of the point on the grid
            IX = MOD(IPX + NGRDX, NGRDX) + 1
            IY = MOD(IPY + NGRDY, NGRDY) + 1
            IZ = MOD(IPZ + NGRDZ, NGRDZ) + 1
C square distance and distance from the ion
            DIST2 = ((IPX-NGRDX/2)*HX)**2
     $           +  ((IPY-NGRDY/2)*HY)**2
     $           +  ((IPZ-NGRDZ/2)*HZ)**2
            DIST = SQRT(DIST2)
#ifdef X_FACE
            W3 = W1*W2 * ABS(IPX-NGRDX/2)*HX / (DIST*DIST2)
#endif
#ifdef Y_FACE
            W3 = W1*W2 * ABS(IPY-NGRDY/2)*HY / (DIST*DIST2)
#endif
#ifdef Z_FACE
            W3 = W1*W2 * ABS(IPZ-NGRDZ/2)*HZ / (DIST*DIST2)
#endif


#ifdef DISTRIB_CHARGE
            SUMCHG = SUMCHG + W3
#ifdef COMPLEX_GRID
            GRID(IX,IY,IZ) = GRID(IX,IY,IZ) 
     $           + CMPLX(W3,0.0)
#else 
            ITMP = IX+(IY-1)*NGRDX+(IZ-1)*NGRDXY
            GRIDR(ITMP) = GRIDR(ITMP) + W3
#endif
#endif


#ifdef GET_DERIV
#ifdef COMPLEX_GRID
            DERIK2 = DERIK2 + CMPLX(W3,0.0)*
     $           GRID(IX,IY,IZ)
#else 
            ITMP = IX+(IY-1)*NGRDX
     $           +(IZ-1)*NGRDXY
            DERIK2 = DERIK2 + W3*GRIDR(ITMP)
#endif
#endif
 


