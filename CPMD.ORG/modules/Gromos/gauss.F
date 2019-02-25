C     ==================================================================
      SUBROUTINE GAUSS(AM,SD,NUM,V,IG)
C     ==--------------------------------------------------------------==
C     W.F. VAN GUNSTEREN, GRONINGEN, JULY 1987
C     This Polar Box-Muller code contributed by Ilario Tironi
C     replaces the GROMOS87 version of GAUSS.
C
C     SUBROUTINE GAUSS (AM,SD,NUM,V,IG)
C
C     GAUSS will supply NUM normally distributed random numbers with
C     a given mean and standard deviation. It uses a polar Box-Muller
C     method to compute normal random numbers.
C     (Black, S.C., Computers in Physics, 1989)
C
C     AM    = the desired mean of the normal distribution
C     SD    = the desired standard deviation of the normal distribution
C     NUM   = the desired number of elements in V to produce
C     V(NUM)= delivered with NUM normal random variables
C     IG    = random number generator seed. This variable is modified
C             by subroutine RANDOM.
C
C     GAUSS uses a random number generator L<RANDOM>.
C
C     This can be replaced by subroutine RANDU(IG,IY,Y) on IBM machines,
C     function RANDU(IG) on VAX machines, or intrinsic function RAN(IG)
C     on CDC machines.
Cmb - Revised on 16 Sept. 2005
C     ==--------------------------------------------------------------==
C args
      INTEGER NUM,IG
      REAL*8 AM,SD,V(NUM)

C local vars
      INTEGER I
      LOGICAL LNEW
      REAL*8 R,STORE,SQLN,W1,W2
      SAVE LNEW,STORE
      DATA LNEW/.FALSE./

C begin
      DO I=1,NUM
         IF (LNEW) THEN
            V(I) = AM + SD*STORE
            LNEW = .FALSE.
         ELSE
C repeat until (r<1.0 and r<>0.0)
 10         CONTINUE
               CALL GR_RANDOM(W1,IG)
               CALL GR_RANDOM(W2,IG)
               W1 = 2.0d0*W1-1.0d0
               W2 = 2.0d0*W2-1.0d0
               R = W1*W1 + W2*W2
            IF ((R.GE.1.0d0).OR.(R.EQ .0.0d0)) GO TO 10
C end of repeat until loop
            SQLN = DSQRT(-2.0D0*DLOG(R)/R)
            V(I) = AM+SD*(W1*SQLN)
            STORE = W2*SQLN
            LNEW = .TRUE.
         ENDIF
      ENDDO
C end gauss
C     ==--------------------------------------------------------------==
      RETURN
      END
