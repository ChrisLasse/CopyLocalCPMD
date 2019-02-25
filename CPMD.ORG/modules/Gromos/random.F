C     ==================================================================
      SUBROUTINE GR_RANDOM(RAND,IG)
C     ==--------------------------------------------------------------==
CCCCCC R. GEURTSEN, GRONINGEN, WFVG, JULY 1987 CCCCCCCCCCCCCCCCCCCCCCCCC
Cmb Revised 15 Sept. 2005
C
C     SUBROUTINE RANDOM (RAND,IG)
C
C     RANDOM GENERATES A RANDOM NUMBER RAND, USING A LINEAR
C     CONGRUENTIAL METHOD. THE RECURSION FORMULA
C
C         IRAND = MOD(IRAND * B + 1, A)
C
C     IS USED WITH  B = 31415821  AND  A = 100000000. THE LAST DIGIT
C     FROM THE RANDOM INTEGER IRAND IS CHOPPED OF, AND THE NUMBER
C     IS SCALED TO A real*8 VALUE RAND BETWEEN 0 AND 1, INCLUDING 0 BUT
C     EXCLUDING 1.
C
C     RAND = DELIVERED WITH RANDOM NUMBER BETWEEN 0 AND 1
C     IG = RANDOM NUMBER GENERATOR SEED, IS DELIVERED WITH A RANDOM
C          INTEGER.
C
C     THIS VERSION WHICH STORES A NUMBER OF RANDOM NUMBERS LOCALLY
C     BY W. SCOTT, ZUERICH, 1995
C     ==--------------------------------------------------------------==
C args
      real*8 RAND
      INTEGER IG
C local params
      INTEGER M,M1,MULT
      PARAMETER (M=100000000,M1=10000,MULT=31415821)
C     how many random numbers to store locally
      INTEGER MAXLOC
      PARAMETER (MAXLOC = 20)
C local vars
      LOGICAL LNEW
      INTEGER IRAND
      SAVE LNEW,IRAND
C
      INTEGER IRANDH,IRANDL,MULTH,MULTL,IUSED,I
      INTEGER ISAVE(MAXLOC)
      real*8 R,RSAVE(MAXLOC)
      SAVE IUSED,ISAVE,RSAVE
C data
      DATA IRAND/0/
      DATA LNEW/.TRUE./
C begin
      IF (LNEW) THEN
         LNEW = .FALSE.
         IRAND = MOD(IABS(IG),M)
         IUSED = MAXLOC
      ENDIF

      IF (IUSED .EQ. MAXLOC) THEN
C     have to generate new numbers which we store
C     into ISAVE and RSAVE
C
         DO I=1,MAXLOC
C*****MULTIPLY IRAND BY MULT, BUT TAKE INTO ACCOUNT THAT OVERFLOW MUST
C*****BE DISCARDED, AND DO NOT GENERATE AN ERROR.
C
            IRANDH = IRAND/M1
            IRANDL = MOD(IRAND,M1)
            MULTH = MULT/M1
            MULTL = MOD(MULT,M1)
C
            IRAND = MOD(IRANDH*MULTL+IRANDL*MULTH,M1)*M1 + IRANDL*MULTL
            IRAND = MOD(IRAND+1,M)
C
C*****CONVERT IRAND TO A real*8 RANDOM NUMBER BETWEEN 0 AND 1.
C
            R = DBLE(IRAND/10)*10/DBLE(M)
            IF ((R.LE.0.0d0).OR.(R.GT.1.0d0)) R = 0.0d0
            RSAVE(I) = R
            ISAVE(I) = IRAND
         ENDDO
         IUSED = 0
      ENDIF

C     now give out a number
      IUSED = IUSED + 1
      RAND = RSAVE(IUSED)
      IG   = ISAVE(IUSED)
C end random
C     ==--------------------------------------------------------------==
      RETURN
      END
