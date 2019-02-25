C     ==================================================================
      SUBROUTINE RESTX(NRC,NATTOT,JRC,CXC,XC,X,F,EC,
     $     NFRSTX,XC0,EC0,LEVERY)
C     ==================================================================
COMMSUBR RESTX W.F. VAN GUNSTEREN, GRONINGEN, APR. 1985 
C
C     RESTX will supply energies and forces for a harmonic
C     potential
C     V(R) = sum over all specified atom of
C            0.5 * CXC * (X-XC)**2.
C     Periodic boundary conditions can be taken into
C     account.
C     NRC   number of atoms for which corrections are to be made
C           NRC .LE. NATTOT
C     JRC(NRC)
C           sequence numbers for these atoms
C     XC(NDIM*NATTOT)
C           reference atom cartesian coordinates
C     CXC(NATTOT)
C           force constants
COMMEND
C     ==--------------------------------------------------------------==
cmb - Revised on 13 Oct. 2005
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'restx.h'
C arguments
      INTEGER NRC,NATTOT,NFRSTX
      INTEGER JRC(NRC)
      LOGICAL LEVERY
      REAL*8 XC(NDIM*NATTOT),X(NDIM*NATTOT),F(NDIM*NATTOT)
      REAL*8 EC,CXC(NATTOT),XC0(NRC),EC0(NRC)
C local vars
      LOGICAL LFIRST
      SAVE LFIRST

      LOGICAL LMONO,LOCTO,LVAC,LDOVIR,LDOTRA
      SAVE LMONO,LOCTO,LVAC,LDOVIR,LDOTRA

      REAL*8 BETAR,COSB,COSB2,BOXOH,BOXOQ,BOXHM2
      SAVE COSB,COSB2,BOXOH,BOXOQ,BOXHM2

      INTEGER I,I3,M,PASS
      LOGICAL LDOFRC,LDOCMP
      REAL*8 XIJDST(MAXDIM),DSTTMP,RIJ2,TMP,BOXHMI
      REAL*8 TRVEC(MAXDIM),AVDF
      REAL*8 DF,ELOC
C data
      DATA LFIRST /.TRUE./
C begin
      IF (LEVERY .OR. LFIRST) THEN
         LFIRST = .FALSE.

         LMONO = (NTB .GT. 0)
         LOCTO = (NTB. LT. 0)
         LVAC  = (NTB .EQ. 0)
         LDOVIR= (ABS(NTB) .EQ. NTBVIR)

         LDOTRA = .FALSE.
         IF (LMONO)THEN
            BETAR = BETA*DATAN(1.D0)/45.0D+00
            COSB = DCOS(BETAR)
            LDOTRA = (DABS(COSB).GE.1.D-4)
            COSB2 = 2.0D0*COSB
         ENDIF

         BOXHMI = BOXH(1)
         DO 2 M=2,NDIM
            IF (BOXH(M) .LT. BOXHMI) THEN
               BOXHMI = BOXH(M)
            ENDIF   
 2       CONTINUE
         BOXHM2 = BOXHMI**2

         IF (LOCTO) THEN
            BOXOH = BOXH(1)
            BOXOQ = BOX(1)*0.75D0
            BOXHM2 = BOXHM2*0.75D0
         ELSE
            BOXOH = 0.0D0
            BOXOQ = 0.0D0
         ENDIF
      ENDIF

      EC = 0.0D0
      IF (NFRSTX .EQ. NFRDNR .OR. NFRSTX .EQ. NFRENR) THEN
C NRC is the constraint to calculate
         I3 = NDIM*(JRC(NRC)-1)
C calc distances
         RIJ2 = 0.0D0
         DO M=1,NDIM
            TMP = X(I3+M) - XC(I3+M)
            RIJ2 = RIJ2 + TMP**2
            XIJDST(M) = TMP
         ENDDO
         IF (.NOT. LVAC .AND. RIJ2 .GE. BOXHM2) THEN
            IF (LDOTRA) THEN
               CALL TRACO(1,0,XIJDST,BETA,1,LEVERY)
            ENDIF
C Phil : Error #13
            RIJ2 = 0.0D0
C
            DO M=1,NDIM
               DSTTMP = XIJDST(M)
               IF (DSTTMP .GE. BOXH(M)) THEN
                  DSTTMP = DSTTMP - BOX(M)
               ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
                  DSTTMP = DSTTMP + BOX(M)
               ENDIF
               XIJDST(M) = DSTTMP
               RIJ2 = RIJ2 + DSTTMP**2
            ENDDO
            IF (LDOTRA) THEN
C Phil : Error #21
               RIJ2 = RIJ2 + 
     $              COSB2*XIJDST(1)*XIJDST(3)
               CALL TRACO(1,0,XIJDST,BETA,-1,LEVERY)
            ENDIF
C
            IF (LOCTO) THEN
               DSTTMP = BOXOQ - DABS(XIJDST(1))-
     $              DABS(XIJDST(2))-DABS(XIJDST(3))
               IF (DSTTMP .LT. 0.0D0) THEN
                  RIJ2 = RIJ2 + DSTTMP*BOX(1)
                  DO M=1,NDIM
                     XIJDST(M)=XIJDST(M)-SIGN(BOXOH,XIJDST(M))
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
C--- end dist calc

         IF (NFRSTX .EQ. NFRDNR) THEN
            EC = DSQRT(RIJ2)
         ELSE
            DF = CXC(JRC(NRC))
            EC = 0.5D0*DF*RIJ2
         ENDIF
         RETURN
      ENDIF


      LDOFRC = (NFRSTX .EQ. NFRALL .OR. NFRSTX .EQ. NFREF)
      LDOCMP = (NFRSTX .EQ. NFRALL .OR. NFRSTX .EQ. NFREDT)

      DO 500 PASS=1,2
         IF (PASS.EQ.1) THEN
            AVDF=0.D0
cmb            DO M=1,NDIM
cmb               TRVEC(M)=0.0D0
cmb            ENDDO
            CALL mm_AZZERO(TRVEC,NDIM)
         ELSE
            DO M=1,NDIM
               TRVEC(M)=TRVEC(M)/AVDF
            ENDDO
            DO I=1,NRC
               I3 = NDIM*(JRC(I)-1)
               DO M=1,NDIM
                  XC(I3+M)=XC(I3+M)+TRVEC(M) 
               ENDDO
            ENDDO
         ENDIF
      DO 101 I=1,NRC
         I3 = NDIM*(JRC(I)-1)
C calc distances
         RIJ2 = 0.0D0
         DO M=1,NDIM
            TMP = X(I3+M) - XC(I3+M)
            RIJ2 = RIJ2 + TMP**2
            XIJDST(M) = TMP
         ENDDO
         IF (.NOT. LVAC .AND. RIJ2 .GE. BOXHM2) THEN
            IF (LDOTRA) THEN
               CALL TRACO(1,0,XIJDST,BETA,1,LEVERY)
            ENDIF
C Phil : Error #13
            RIJ2 = 0.0D0
C
            DO M=1,NDIM
               DSTTMP = XIJDST(M)
               IF (DSTTMP .GE. BOXH(M)) THEN
                  DSTTMP = DSTTMP - BOX(M)
               ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
                  DSTTMP = DSTTMP + BOX(M)
               ENDIF
               XIJDST(M) = DSTTMP
               RIJ2 = RIJ2 + DSTTMP**2
            ENDDO
            IF (LDOTRA) THEN
C Phil : Error #21
               RIJ2 = RIJ2 + 
     $              COSB2*XIJDST(1)*XIJDST(3)
               CALL TRACO(1,0,XIJDST,BETA,-1,LEVERY)
            ENDIF
C
            IF (LOCTO) THEN
               DSTTMP = BOXOQ - DABS(XIJDST(1))-
     $              DABS(XIJDST(2))-DABS(XIJDST(3))
               IF (DSTTMP .LT. 0.0D0) THEN
                  RIJ2 = RIJ2 + DSTTMP*BOX(1)
                  DO M=1,NDIM
                     XIJDST(M)=XIJDST(M)-SIGN(BOXOH,XIJDST(M))
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
C--- end dist calc

         DF = CXC(JRC(I))
         IF (PASS.EQ.1) THEN
            AVDF=AVDF+DF
            DO M=1,NDIM
               TRVEC(M)=TRVEC(M)+DF*XIJDST(M)
            ENDDO
         ELSE
            ELOC = DF*RIJ2
            EC = EC + ELOC
            IF (LDOCMP) THEN
               EC0(I) = ELOC*0.5D0
               XC0(I) = DSQRT(RIJ2)
            ENDIF
            IF (LDOFRC) THEN
               DO M=1,NDIM
                  F(I3+M) = F(I3+M) - DF*XIJDST(M)
               ENDDO
            ENDIF
         ENDIF
 101  CONTINUE
 500  CONTINUE 

      EC = EC*0.5D0
C end restx
C     ==--------------------------------------------------------------==
      RETURN
      END
