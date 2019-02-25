#ifdef __SR11000
*option MP(P(0)), LANGLVL(SAVE(0))
#endif
C     ==================================================================
      SUBROUTINE C06FCF(X,Y,PTS,WORK,IFAIL)
C     ==--------------------------------------------------------------==
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
Cmb - Revised on 21 June 2006
C     COMPLEX FOURIER TRANSFORM
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(PTS), X(PTS), Y(PTS)
C     .. Local Scalars ..
      DOUBLE PRECISION  SQPTS
      INTEGER           IERROR, IPTS, PMAX, TWOGRP
C     .. Local Arrays ..
      INTEGER           RFACT(21), TFACT(21)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
c      EXTERNAL          P01ABF
C     .. External Subroutines ..
c      EXTERNAL          C06ECW, C06FAY, C06FAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE,DSQRT
C     .. Data statements ..
      DATA              PMAX/19/
      DATA              TWOGRP/8/
C     .. Executable Statements ..
      IF (PTS.LE.1) GO TO 40
      CALL C06FAZ(PTS,PMAX,TWOGRP,TFACT,RFACT,IERROR)
      IF (IERROR.NE.0) GO TO 60
      CALL C06ECW(X,Y,PTS,TFACT)
      CALL C06FAY(X,PTS,RFACT,WORK)
      CALL C06FAY(Y,PTS,RFACT,WORK)
      SQPTS = DSQRT(DBLE(PTS))
!$OMP parallel do private(IPTS) shared(SQPTS)
#ifdef __SR11000
*poption parallel, tlocal(IPTS)
#endif
      DO IPTS = 1, PTS
         X(IPTS) = X(IPTS)/SQPTS
         Y(IPTS) = Y(IPTS)/SQPTS
      ENDDO
      IFAIL = 0
      GO TO 80
C
   40 IERROR = 3
   60 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
   80 RETURN
      END
C     ==================================================================
      SUBROUTINE C06FAZ(PTS,PMAX,TWOGRP,TFACT,RFACT,IERROR)
C     ==--------------------------------------------------------------==
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COPY REORDERING FACTORING PROGRAMME
C     .. Scalar Arguments ..
      INTEGER           IERROR, PMAX, PTS, TWOGRP
C     .. Array Arguments ..
      INTEGER           RFACT(21), TFACT(21)
C     .. Local Scalars ..
      INTEGER           F, J, JJ, N, NEST, P, PTWO, Q
C     .. Local Arrays ..
      INTEGER           PP(10), QQ(20)
C     .. Data statements ..
      DATA              NEST/20/
C     .. Executable Statements ..
      N = PTS
      F = 2
      P = 0
      Q = 0
   20 CONTINUE
      IF (N.LE.1) GO TO 100
      DO 40 J = F, PMAX 
         IF (N.EQ.(N/J)*J) GO TO 60
   40 CONTINUE
      GO TO 280
   60 CONTINUE
      IF (2*P+Q.GE.NEST) GO TO 300
      F = J
      N = N/F
      IF (N.EQ.(N/F)*F) GO TO 80
      Q = Q + 1
      QQ(Q) = F
      GO TO 20
   80 CONTINUE
      N = N/F
      P = P + 1
      PP(P) = F
      GO TO 20
C     
  100 CONTINUE
      IF (P.LT.1) GO TO 140
#ifdef _vpp_
!OCL SCALAR
#endif 
      DO J = 1, P
         JJ = P + 1 - J
         TFACT(J) = PP(JJ)
         JJ = P + Q + J
         TFACT(JJ) = PP(J)
      ENDDO
  140 CONTINUE
      IF (Q.LT.1) GO TO 180
!$OMP parallel do private(J,JJ)
#ifdef _vpp_
!OCL SCALAR
#endif
      DO J = 1, Q
         JJ = P + J
         TFACT(JJ) = QQ(J)
      ENDDO
  180 CONTINUE
      JJ = 2*P + Q
      TFACT(JJ+1) = 0
      RFACT(JJ+1) = 0
!$OMP parallel do private(J)
#ifdef __SR11000
*poption parallel, tlocal(J)
#endif
      DO J = 1, JJ
         RFACT(J) = TFACT(J)
      ENDDO
      IF (JJ.EQ.1) RFACT(1) = 0
      PTWO = 1
      J = 0
  220 CONTINUE
      J = J + 1
      IF (TFACT(J).EQ.0) GO TO 260
      IF (TFACT(J).NE.2) GO TO 220
      PTWO = PTWO*2
      TFACT(J) = 1
      IF (PTWO.GE.TWOGRP) GO TO 240
      IF (TFACT(J+1).EQ.2) GO TO 220
  240 CONTINUE
      TFACT(J) = PTWO
      PTWO = 1
      GO TO 220
  260 CONTINUE
      IERROR = 0
      GO TO 320
  280 IERROR = 1
      GO TO 320
C
  300 IERROR = 2
  320 RETURN
      END
C     ==================================================================
      SUBROUTINE C06FAY(X,PTS,FACTOR,WORK)
C     ==--------------------------------------------------------------==
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11 REVISED. IER-444 (FEB 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     SINGLE COPY REORDERING PROGRAMME
C     EQUIVALENCE (I1,I(1)), (K1,K(1)), (L1,L(1))
C     EQUIVALENCE (I2,I(2)), (K2,K(2)), (L2,L(2))
C     EQUIVALENCE (I3,I(3)), (K3,K(3)), (L3,L(3))
C     EQUIVALENCE (I4,I(4)), (K4,K(4)), (L4,L(4))
C     EQUIVALENCE (I5,I(5)), (K5,K(5)), (L5,L(5))
C     EQUIVALENCE (I6,I(6)), (K6,K(6)), (L6,L(6))
C     EQUIVALENCE (I7,I(7)), (K7,K(7)), (L7,L(7))
C     EQUIVALENCE (I8,I(8)), (K8,K(8)), (L8,L(8))
C     EQUIVALENCE (I9,I(9)), (K9,K(9)), (L9,L(9))
C     EQUIVALENCE (I10,I(10)), (K10,K(10)), (L10,L(10))
C     EQUIVALENCE (K11,K(11))
C     .. Scalar Arguments ..
      INTEGER           PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(PTS), X(PTS)
      INTEGER           FACTOR(21)
C     .. Local Scalars ..
      INTEGER           I1, I10, I2, I3, I4, I5, I6, I7, I8, I9, J, JJ,
     *                  K1, K10, K11, K2, K3, K4, K5, K6, K7, K8, K9,
     *                  KK, L1, L10, L2, L3, L4, L5, L6, L7, L8, L9,
     *                  LEVEL, LOOP, NEST
C     .. Local Arrays ..
      INTEGER           I(20), K(20), L(20)
C     .. Data statements ..
      DATA              NEST/20/
      DATA              LOOP/10/
C     .. Executable Statements ..
      IF (FACTOR(1).EQ.0) GO TO 380
!$OMP parallel do private(J)
#ifdef __SR11000
*poption parallel, tlocal(J)
#endif
      DO J = 1, NEST
         L(J) = 1
         I(J) = 1
      ENDDO
      KK = PTS
      DO 40 J = 1, NEST
         IF (FACTOR(J).EQ.0) GO TO 60
         L(J) = KK
         KK = KK/FACTOR(J)
         I(J) = KK
   40 CONTINUE
   60 CONTINUE
C
      L1 = L(1)
      L2 = L(2)
      L3 = L(3)
      L4 = L(4)
      L5 = L(5)
      L6 = L(6)
      L7 = L(7)
      L8 = L(8)
      L9 = L(9)
      L10 = L(10)
      I1 = I(1)
      I2 = I(2)
      I3 = I(3)
      I4 = I(4)
      I5 = I(5)
      I6 = I(6)
      I7 = I(7)
      I8 = I(8)
      I9 = I(9)
      I10 = I(10)
C
      KK = 0
      LEVEL = NEST
      K(LEVEL) = 1
      GO TO 100
   80 CONTINUE
      IF (LEVEL.GE.NEST) GO TO 340
      LEVEL = LEVEL + 1
      K(LEVEL) = K(LEVEL) + I(LEVEL)
      IF (K(LEVEL).GT.L(LEVEL)) GO TO 80
  100 CONTINUE
      LEVEL = LEVEL - 1
      DO J = LOOP, LEVEL
         JJ = LEVEL + LOOP - J
         K(JJ) = K(JJ+1)
      ENDDO
      K11 = K(11)
#ifdef _vpp_
!OCL SCALAR
#endif 
      DO 320 K10 = K11, L10, I10
         DO 300 K9 = K10, L9, I9
            DO 280 K8 = K9, L8, I8
               DO 260 K7 = K8, L7, I7
                  DO 240 K6 = K7, L6, I6
                     DO 220 K5 = K6, L5, I5
                        DO 200 K4 = K5, L4, I4
                           DO 180 K3 = K4, L3, I3
                              DO 160 K2 = K3, L2, I2
                                 DO 140 K1 = K2, L1, I1
                                    KK = KK + 1
                                    WORK(KK) = X(K1)
  140                            CONTINUE
  160                         CONTINUE
  180                      CONTINUE
  200                   CONTINUE
  220                CONTINUE
  240             CONTINUE
  260          CONTINUE
  280       CONTINUE
  300    CONTINUE
  320 CONTINUE
      LEVEL = LOOP
      GO TO 80
  340 CONTINUE
!$OMP parallel do private(J)
#ifdef __SR11000
*poption parallel, tlocal(J)
#endif
      DO J = 1, PTS
         X(J) = WORK(J)
      ENDDO
C     
  380 CONTINUE
      RETURN
      END 
C     ==================================================================
      SUBROUTINE C06ECW(X,Y,PTS,FACTOR)
C     ==--------------------------------------------------------------==
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COMPLEX FOURIER TRANSFORM KERNEL DRIVER
C     .. Scalar Arguments ..
      INTEGER           PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X(PTS), Y(PTS)
      INTEGER           FACTOR(21)
C     .. Local Scalars ..
      INTEGER           F, M, M1, M2, M3, M4, M5, M6, M7, P
C     .. External Subroutines ..
c      EXTERNAL          C06ECQ, C06ECR, C06ECS, C06ECT, C06ECU, C06ECV
C     .. Executable Statements ..
      F = 0
      M = PTS
   20 CONTINUE
      F = F + 1
      P = FACTOR(F)
      IF (P.EQ.0) RETURN
      IF (P.EQ.1) GO TO 20
      M = M/P
      M1 = PTS - M
      M2 = M1 - M
      M3 = M2 - M
      M4 = M3 - M
      M5 = M4 - M
      M6 = M5 - M
      M7 = M6 - M
      IF (P.EQ.2) GO TO 40
      IF (P.EQ.3) GO TO 60
      IF (P.EQ.4) GO TO 80
      IF (P.EQ.5) GO TO 100
      IF (P.EQ.8) GO TO 120
      GO TO 140
C
   40 CONTINUE
      CALL C06ECV(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,M)
      GO TO 20
C
   60 CONTINUE
      CALL C06ECU(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1),M2,M)
      GO TO 20
C
   80 CONTINUE
      CALL C06ECT(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1)
     *            ,M2,X(3*M+1),Y(3*M+1),M3,M)
      GO TO 20
C
  100 CONTINUE
      CALL C06ECS(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1)
     *            ,M2,X(3*M+1),Y(3*M+1),M3,X(4*M+1),Y(4*M+1),M4,M)
      GO TO 20
C
  120 CONTINUE
      CALL C06ECR(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1)
     *            ,M2,X(3*M+1),Y(3*M+1),M3,X(4*M+1),Y(4*M+1),M4,X(5*M+1)
     *            ,Y(5*M+1),M5,X(6*M+1),Y(6*M+1),M6,X(7*M+1),Y(7*M+1)
     *            ,M7,M)
      GO TO 20
C
  140 CONTINUE
      CALL C06ECQ(X,Y,PTS,M,P)
      GO TO 20
C
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE C06ECV(X0,Y0,PTS,X1,Y1,M1,M)
C     ==--------------------------------------------------------------==
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX TWO COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), Y0(PTS), Y1(M1)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANGLE, C, IS, IU, RS, RU, S, TWOPI
      INTEGER           J, K, K0, M2, MOVER2
      LOGICAL           FOLD, ZERO
C     .. External Functions ..
cmb       DOUBLE PRECISION  X01AAF
c      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         DCOS, DBLE, DSIN
C     .. Executable Statements ..
      M2 = M*2          
      MOVER2 = M/2 + 1  
cmb      TWOPI = 2.0D0*X01AAF(0.0D0)
      TWOPI = 8.0D0*DATAN(1.0D0)
C
cmbike... J=1
!$OMP parallel do private(K,RS,IS,RU,IU)
#ifdef __SR11000
*poption parallel, tlocal(K,RS,IS,RU,IU)
#endif
      DO K = 1, PTS, M2
        RS = X0(K) + X1(K)
        IS = Y0(K) + Y1(K)
        RU = X0(K) - X1(K)
        IU = Y0(K) - Y1(K)
        X0(K) = RS
        Y0(K) = IS
        X1(K) = RU
        Y1(K) = IU
      ENDDO
cmbike... J=2,...,MOVER2
      DO J = 2, MOVER2
         FOLD = 2*J .LT. M + 2
         K0 = J         
         ANGLE = TWOPI*DBLE(J-1)/DBLE(M2)
         C = DCOS(ANGLE)
         S = DSIN(ANGLE)
         GO TO 40
   20    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
         C = -C
   40    CONTINUE
cmb
!$OMP parallel do private(K,RS,IS,RU,IU)
#ifdef __SR11000
*poption parallel, tlocal(K,RS,IS,RU,IU)
#endif
         DO K = K0, PTS, M2
           RS = X0(K) + X1(K)
           IS = Y0(K) + Y1(K)
           RU = X0(K) - X1(K)
           IU = Y0(K) - Y1(K)
           X0(K) = RS
           Y0(K) = IS
           X1(K) = RU*C + IU*S
           Y1(K) = IU*C - RU*S
         ENDDO
         IF (FOLD) GO TO 20
      ENDDO
C
      RETURN
      END
C     ==================================================================
      SUBROUTINE C06ECU(X0,Y0,PTS,X1,Y1,M1,X2,Y2,M2,M)
C     ==--------------------------------------------------------------==
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX THREE COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, M2, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), X2(M2), Y0(PTS), Y1(M1), Y2(M2)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, ANGLE, B, C1, C2, I0, I1, I2, IA, IB, IS, R0,
     *                  R1, R2, RA, RB, RS, S1, S2, T, TWOPI
      INTEGER           J, K, K0, M3, MOVER2
      LOGICAL           FOLD, ZERO
C     .. External Functions ..
cmb      DOUBLE PRECISION  X01AAF
c      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         DCOS, DBLE, DSIN, DSQRT
C     .. Executable Statements ..
      M3 = M*3
      MOVER2 = M/2 + 1
cmb      TWOPI = 2.0D0*X01AAF(0.0D0)
      TWOPI = 8.0D0*DATAN(1.0D0)
C     A = COS(TWOPI/3.0) 
C     B = SIN(TWOPI/3.0) 
      A = -0.5D0
      B = DSQRT(0.75D0)
C
      DO 120 J = 1, MOVER2
         FOLD = J .GT. 1 .AND. 2*J .LT. M + 2
         K0 = J
         ANGLE = TWOPI*DBLE(J-1)/DBLE(M3)
         ZERO = ANGLE .EQ. 0.0D0
         C1 = DCOS(ANGLE)
         S1 = DSIN(ANGLE)
         C2 = C1*C1 - S1*S1
         S2 = S1*C1 + C1*S1
         GO TO 40
   20    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
         T = C1*A + S1*B
         S1 = C1*B - S1*A
         C1 = T
         T = C2*A - S2*B
         S2 = -C2*B - S2*A
         C2 = T
   40    CONTINUE
cmb
         IF (ZERO) THEN
!$OMP parallel do private(K,R0,I0,RS,IS,RA,IA,RB,IB)
#ifdef __SR11000
*poption parallel, tlocal(K,R0,I0,RS,IS,RA,IA,RB,IB)
#endif
         DO K = K0, PTS, M3
            R0 = X0(K)
            I0 = Y0(K)
            RS = X1(K) + X2(K)
            IS = Y1(K) + Y2(K)
            X0(K) = R0 + RS
            Y0(K) = I0 + IS
            RA = R0 + RS*A
            IA = I0 + IS*A
            RB = (X1(K)-X2(K))*B
            IB = (Y1(K)-Y2(K))*B
            X1(K) = RA + IB
            Y1(K) = IA - RB
            X2(K) = RA - IB
            Y2(K) = IA + RB
         ENDDO
         ELSE
!$OMP parallel do private(K,R0,I0,RS,IS,RA,IA,RB,IB,R1,I1,R2,I2)
#ifdef __SR11000
*poption parallel, tlocal(K,R0,I0,RS,IS,RA,IA,RB,IB,R1,I1,R2,I2)
#endif
         DO K = K0, PTS, M3
            R0 = X0(K)
            I0 = Y0(K)
            RS = X1(K) + X2(K)
            IS = Y1(K) + Y2(K)
            X0(K) = R0 + RS
            Y0(K) = I0 + IS
            RA = R0 + RS*A
            IA = I0 + IS*A
            RB = (X1(K)-X2(K))*B
            IB = (Y1(K)-Y2(K))*B
            R1 = RA + IB
            I1 = IA - RB
            R2 = RA - IB
            I2 = IA + RB
            X1(K) = R1*C1 + I1*S1
            Y1(K) = I1*C1 - R1*S1
            X2(K) = R2*C2 + I2*S2
            Y2(K) = I2*C2 - R2*S2
         ENDDO
cmb
         ENDIF
         IF (FOLD) GO TO 20
  120 CONTINUE
C
      RETURN
      END
C     ==================================================================
      SUBROUTINE C06ECT(X0,Y0,PTS,X1,Y1,M1,X2,Y2,M2,X3,Y3,M3,M)
C     ----------------------------------------------------------------==
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX FOUR COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, M2, M3, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), X2(M2), X3(M3), Y0(PTS),
     *                  Y1(M1), Y2(M2), Y3(M3)
C     .. Local Scalars .. 
      DOUBLE PRECISION  ANGLE, C1, C2, C3, I1, I2, I3, IS0, IS1, IU0,
     *                  IU1, R1, R2, R3, RS0, RS1, RU0, RU1, S1, S2, S3,
     *                  T, TWOPI
      INTEGER           J, K, K0, M4, MOVER2
      LOGICAL           FOLD, ZERO
C     .. External Functions ..
cmb      DOUBLE PRECISION  X01AAF
c      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         DCOS, DBLE, DSIN
C     .. Executable Statements ..
      M4 = M*4
      MOVER2 = M/2 + 1
cmb      TWOPI = 2.0D0*X01AAF(0.0D0)
      TWOPI = 8.0D0*DATAN(1.0D0)
cmbike... J=1
!$OMP parallel do private(K,RS0,IS0,RU0,IU0,RS1,IS1,RU1,IU1)
#ifdef __SR11000
*poption parallel, tlocal(K,RS0,IS0,RU0,IU0,RS1,IS1,RU1,IU1)
#endif
      DO K = 1, PTS, M4
        RS0 = X0(K) + X2(K)
        IS0 = Y0(K) + Y2(K)
        RU0 = X0(K) - X2(K)
        IU0 = Y0(K) - Y2(K)
        RS1 = X1(K) + X3(K)
        IS1 = Y1(K) + Y3(K)
        RU1 = X1(K) - X3(K)
        IU1 = Y1(K) - Y3(K)
        X0(K) = RS0 + RS1
        Y0(K) = IS0 + IS1
        X2(K) = RU0 + IU1
        Y2(K) = IU0 - RU1
        X1(K) = RS0 - RS1
        Y1(K) = IS0 - IS1
        X3(K) = RU0 - IU1
        Y3(K) = IU0 + RU1
      ENDDO
cmbike... J=2,...,MOVER2
      DO J = 2, MOVER2
         FOLD = 2*J .LT. M + 2
         K0 = J
         ANGLE = TWOPI*DBLE(J-1)/DBLE(M4)
         C1 = DCOS(ANGLE)
         S1 = DSIN(ANGLE)
         C2 = C1*C1 - S1*S1
         S2 = S1*C1 + C1*S1
         C3 = C2*C1 - S2*S1
         S3 = S2*C1 + C2*S1
         GO TO 40
   20    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
         T = C1
         C1 = S1
         S1 = T
         C2 = -C2
         T = C3
         C3 = -S3
         S3 = -T
   40    CONTINUE
C
!$OMP parallel do private(K,RS0,IS0,RU0,IU0,RS1,IS1,RU1,IU1)
!$OMP+ private(R1,I1,R2,I2,R3,I3)
#ifdef __SR11000
*poption parallel
*poption tlocal(K,RS0,IS0,RU0,IU0,RS1,IS1,RU1,IU1)
*poption tlocal(R1,I1,R2,I2,R3,I3)
#endif
         DO K = K0, PTS, M4
           RS0 = X0(K) + X2(K)
           IS0 = Y0(K) + Y2(K)
           RU0 = X0(K) - X2(K)
           IU0 = Y0(K) - Y2(K)
           RS1 = X1(K) + X3(K)
           IS1 = Y1(K) + Y3(K)
           RU1 = X1(K) - X3(K)
           IU1 = Y1(K) - Y3(K)
           X0(K) = RS0 + RS1
           Y0(K) = IS0 + IS1
           R1 = RU0 + IU1
           I1 = IU0 - RU1
           R2 = RS0 - RS1
           I2 = IS0 - IS1
           R3 = RU0 - IU1
           I3 = IU0 + RU1
           X2(K) = R1*C1 + I1*S1
           Y2(K) = I1*C1 - R1*S1
           X1(K) = R2*C2 + I2*S2
           Y1(K) = I2*C2 - R2*S2
           X3(K) = R3*C3 + I3*S3
           Y3(K) = I3*C3 - R3*S3
         ENDDO
         IF (FOLD) GO TO 20
      ENDDO
C  
      RETURN
      END
C     ==================================================================
      SUBROUTINE C06ECS(X0,Y0,PTS,X1,Y1,M1,X2,Y2,M2,X3,Y3,M3,X4,Y4,M4,M)
C     ----------------------------------------------------------------==
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX FIVE COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, M2, M3, M4, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), X2(M2), X3(M3), X4(M4),
     *                  Y0(PTS), Y1(M1), Y2(M2), Y3(M3), Y4(M4)
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2, ANGLE, AS, AU, B1, B2, C1, C2, C3, C4,
     *                  I0, I1, I2, I3, I4, IA1, IA2, IAS, IAU, IB1,
     *                  IB2, IS1, IS2, ISS, IU1, IU2, R0, R1, R2, R3,
     *                  R4, RA1, RA2, RAS, RAU, RB1, RB2, RS1, RS2, RSS,
     *                  RU1, RU2, S1, S2, S3, S4, T, TWOPI
      INTEGER           J, K, K0, M5, MOVER2
      LOGICAL           FOLD, ZERO
C     .. External Functions ..
cmb      DOUBLE PRECISION  X01AAF
c      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         DCOS, DBLE, DSIN, DSQRT
C     .. Executable Statements ..
      M5 = M*5
      MOVER2 = M/2 + 1
cmb      TWOPI = 2.0D0*X01AAF(0.0D0)
      TWOPI = 8.0D0*DATAN(1.0D0)
      A1 = DCOS(0.2D0*TWOPI)
      B1 = DSIN(0.2D0*TWOPI)
      A2 = DCOS(0.4D0*TWOPI)
      B2 = DSIN(0.4D0*TWOPI)
      AS = -0.25D0
      AU = DSQRT(5.0D0)/4.0D0
C
      DO 120 J = 1, MOVER2
         FOLD = J .GT. 1 .AND. 2*J .LT. M + 2
         K0 = J
         ANGLE = TWOPI*DBLE(J-1)/DBLE(M5)
         ZERO = ANGLE .EQ. 0.0D0
         C1 = DCOS(ANGLE)
         S1 = DSIN(ANGLE)
         C2 = C1*C1 - S1*S1
         S2 = S1*C1 + C1*S1
         C3 = C2*C1 - S2*S1
         S3 = S2*C1 + C2*S1
         C4 = C2*C2 - S2*S2
         S4 = S2*C2 + C2*S2
         GO TO 40
   20    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
         T = C1*A1 + S1*B1
         S1 = C1*B1 - S1*A1
         C1 = T
         T = C2*A2 + S2*B2
         S2 = C2*B2 - S2*A2
         C2 = T
         T = C3*A2 - S3*B2
         S3 = -C3*B2 - S3*A2
         C3 = T
         T = C4*A1 - S4*B1
         S4 = -C4*B1 - S4*A1
         C4 = T
   40    CONTINUE
cmb
         IF (ZERO) THEN
!$OMP parallel do private(K,R0,I0,RS1,IS1,RU1,IU1,RS2,IS2)
!$OMP+    private(RU2,IU2,RSS,ISS,RAS,IAS,RAU,IAU) 
!$OMP+    private(RA1,IA1,RA2,IA2,RB1,IB1,RB2,IB2)
#ifdef __SR11000
*poption parallel
*poption tlocal(K,R0,I0,RS1,IS1,RU1,IU1,RS2,IS2)
*poption tlocal(RU2,IU2,RSS,ISS,RAS,IAS,RAU,IAU)
*poption tlocal(RA1,IA1,RA2,IA2,RB1,IB1,RB2,IB2)
#endif
          DO K = K0, PTS, M5
            R0 = X0(K)
            I0 = Y0(K)
            RS1 = X1(K) + X4(K)
            IS1 = Y1(K) + Y4(K)
            RU1 = X1(K) - X4(K)
            IU1 = Y1(K) - Y4(K)
            RS2 = X2(K) + X3(K)
            IS2 = Y2(K) + Y3(K)
            RU2 = X2(K) - X3(K)
            IU2 = Y2(K) - Y3(K)
            RSS = RS1 + RS2
            ISS = IS1 + IS2
            RAS = R0 + RSS*AS
            IAS = I0 + ISS*AS
            RAU = (RS1-RS2)*AU
            IAU = (IS1-IS2)*AU
            RA1 = RAS + RAU
            IA1 = IAS + IAU
            RA2 = RAS - RAU
            IA2 = IAS - IAU
            RB1 = RU1*B1 + RU2*B2
            IB1 = IU1*B1 + IU2*B2
            RB2 = RU1*B2 - RU2*B1
            IB2 = IU1*B2 - IU2*B1
            X0(K) = R0 + RSS
            Y0(K) = I0 + ISS
            X1(K) = RA1 + IB1
            Y1(K) = IA1 - RB1
            X2(K) = RA2 + IB2
            Y2(K) = IA2 - RB2
            X3(K) = RA2 - IB2
            Y3(K) = IA2 + RB2
            X4(K) = RA1 - IB1
            Y4(K) = IA1 + RB1
          ENDDO
         ELSE 
!$OMP parallel do private(K,R0,I0,RS1,IS1,RU1,IU1,RS2,IS2)
!$OMP+ private(RU2,IU2,RSS,ISS,RAS,IAS,RAU,IAU) 
!$OMP+ private(RA1,IA1,RA2,IA2,RB1,IB1,RB2,IB2)
!$OMP+ private(R1,I1,R2,I2,R3,I3,R4,I4)
#ifdef __SR11000
*poption parallel
*poption tlocal(K,R0,I0,RS1,IS1,RU1,IU1,RS2,IS2)
*poption tlocal(RU2,IU2,RSS,ISS,RAS,IAS,RAU,IAU)
*poption tlocal(RA1,IA1,RA2,IA2,RB1,IB1,RB2,IB2)
*poption tlocal(R1,I1,R2,I2,R3,I3,R4,I4)
#endif
          DO K = K0, PTS, M5
            R0 = X0(K)
            I0 = Y0(K)
            RS1 = X1(K) + X4(K)
            IS1 = Y1(K) + Y4(K)
            RU1 = X1(K) - X4(K)
            IU1 = Y1(K) - Y4(K)
            RS2 = X2(K) + X3(K)
            IS2 = Y2(K) + Y3(K)
            RU2 = X2(K) - X3(K)
            IU2 = Y2(K) - Y3(K)
            RSS = RS1 + RS2
            ISS = IS1 + IS2
            RAS = R0 + RSS*AS
            IAS = I0 + ISS*AS
            RAU = (RS1-RS2)*AU
            IAU = (IS1-IS2)*AU
            RA1 = RAS + RAU
            IA1 = IAS + IAU
            RA2 = RAS - RAU
            IA2 = IAS - IAU
            RB1 = RU1*B1 + RU2*B2
            IB1 = IU1*B1 + IU2*B2
            RB2 = RU1*B2 - RU2*B1
            IB2 = IU1*B2 - IU2*B1
            X0(K) = R0 + RSS
            Y0(K) = I0 + ISS
            R1 = RA1 + IB1
            I1 = IA1 - RB1
            R2 = RA2 + IB2
            I2 = IA2 - RB2
            R3 = RA2 - IB2
            I3 = IA2 + RB2
            R4 = RA1 - IB1 
            I4 = IA1 + RB1 
            X1(K) = R1*C1 + I1*S1
            Y1(K) = I1*C1 - R1*S1
            X2(K) = R2*C2 + I2*S2
            Y2(K) = I2*C2 - R2*S2
            X3(K) = R3*C3 + I3*S3
            Y3(K) = I3*C3 - R3*S3
            X4(K) = R4*C4 + I4*S4
            Y4(K) = I4*C4 - R4*S4
          ENDDO
         ENDIF
cmb
         IF (FOLD) GO TO 20
  120 CONTINUE
C
      RETURN
      END
C     ==================================================================
      SUBROUTINE C06ECR(X0,Y0,PTS,X1,Y1,M1,X2,Y2,M2,X3,Y3,M3,X4,Y4,M4,
     *                  X5,Y5,M5,X6,Y6,M6,X7,Y7,M7,M)
C     ----------------------------------------------------------------==
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX EIGHT COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, M2, M3, M4, M5, M6, M7, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), X2(M2), X3(M3), X4(M4), X5(M5),
     *                  X6(M6), X7(M7), Y0(PTS), Y1(M1), Y2(M2), Y3(M3),
     *                  Y4(M4), Y5(M5), Y6(M6), Y7(M7)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANGLE, C1, C2, C3, C4, C5, C6, C7, E, I1, I2,
     *                  I3, I4, I5, I6, I7, IS0, IS1, IS2, IS3, ISS0,
     *                  ISS1, ISU0, ISU1, IU0, IU1, IU2, IU3, IUS0,
     *                  IUS1, IUU0, IUU1, R1, R2, R3, R4, R5, R6, R7,
     *                  RS0, RS1, RS2, RS3, RSS0, RSS1, RSU0, RSU1, RU0,
     *                  RU1, RU2, RU3, RUS0, RUS1, RUU0, RUU1, S1, S2,
     *                  S3, S4, S5, S6, S7, T, TWOPI
      INTEGER           J, K, K0, M8, MOVER2
      LOGICAL           FOLD, ZERO
C     .. External Functions ..
cmb      DOUBLE PRECISION  X01AAF
c      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         DCOS, DBLE, DSIN
C     .. Executable Statements ..
      M8 = M*8
      MOVER2 = M/2 + 1
cmb      TWOPI = 2.0D0*X01AAF(0.0D0)
      TWOPI = 8.0D0*DATAN(1.0D0)
      E = DCOS(TWOPI/8.0D0)
C
cmbike... J=1
      C1 = 1.0D0
      S1 = 0.0D0
      C2 = C1*C1 - S1*S1
      S2 = S1*C1 + C1*S1
      C3 = C2*C1 - S2*S1
      S3 = S2*C1 + C2*S1
      C4 = C2*C2 - S2*S2
      S4 = S2*C2 + C2*S2
      C5 = C4*C1 - S4*S1
      S5 = S4*C1 + C4*S1
      C6 = C4*C2 - S4*S2
      S6 = S4*C2 + C4*S2
      C7 = C4*C3 - S4*S3
      S7 = S4*C3 + C4*S3
!$OMP parallel do private(K,RS0,IS0,RU0,IU0,RS1,IS1)
!$OMP+ private(RU1,IU1,RS2,IS2,RU2,IU2,RS3,IS3)
!$OMP+ private(RU3,IU3,RSS0,ISS0,RSU0,ISU0)
!$OMP+ private(RSS1,ISS1,RSU1,ISU1,RUS0,IUS0)
!$OMP+ private(RUU0,IUU0,RUS1,IUS1,RUU1,IUU1)
!$OMP+ private(T) shared(E)
#ifdef __SR11000
*poption parallel
*poption tlocal(K,RS0,IS0,RU0,IU0,RS1,IS1)
*poption tlocal(RU1,IU1,RS2,IS2,RU2,IU2,RS3,IS3)
*poption tlocal(RU3,IU3,RSS0,ISS0,RSU0,ISU0)
*poption tlocal(RSS1,ISS1,RSU1,ISU1,RUS0,IUS0)
*poption tlocal(RUU0,IUU0,RUS1,IUS1,RUU1,IUU1,T)
#endif
      DO K = 1, PTS, M8
        RS0 = X0(K) + X4(K)
        IS0 = Y0(K) + Y4(K)
        RU0 = X0(K) - X4(K)
        IU0 = Y0(K) - Y4(K)
        RS1 = X1(K) + X5(K)
        IS1 = Y1(K) + Y5(K)
        RU1 = X1(K) - X5(K)
        IU1 = Y1(K) - Y5(K)
        RS2 = X2(K) + X6(K)
        IS2 = Y2(K) + Y6(K)
        RU2 = X2(K) - X6(K)
        IU2 = Y2(K) - Y6(K)
        RS3 = X3(K) + X7(K)
        IS3 = Y3(K) + Y7(K)
        RU3 = X3(K) - X7(K)
        IU3 = Y3(K) - Y7(K)
        RSS0 = RS0 + RS2
        ISS0 = IS0 + IS2
        RSU0 = RS0 - RS2
        ISU0 = IS0 - IS2
        RSS1 = RS1 + RS3
        ISS1 = IS1 + IS3
        RSU1 = RS1 - RS3
        ISU1 = IS1 - IS3
        RUS0 = RU0 - IU2
        IUS0 = IU0 + RU2
        RUU0 = RU0 + IU2
        IUU0 = IU0 - RU2
        RUS1 = RU1 - IU3
        IUS1 = IU1 + RU3
        RUU1 = RU1 + IU3
        IUU1 = IU1 - RU3
        T = (RUS1+IUS1)*E
        IUS1 = (IUS1-RUS1)*E
        RUS1 = T
        T = (RUU1+IUU1)*E
        IUU1 = (IUU1-RUU1)*E
        RUU1 = T
        X0(K) = RSS0 + RSS1
        Y0(K) = ISS0 + ISS1
        X4(K) = RUU0 + RUU1
        Y4(K) = IUU0 + IUU1
        X2(K) = RSU0 + ISU1
        Y2(K) = ISU0 - RSU1
        X6(K) = RUS0 + IUS1
        Y6(K) = IUS0 - RUS1
        X1(K) = RSS0 - RSS1
        Y1(K) = ISS0 - ISS1
        X5(K) = RUU0 - RUU1
        Y5(K) = IUU0 - IUU1
        X3(K) = RSU0 - ISU1
        Y3(K) = ISU0 + RSU1
        X7(K) = RUS0 - IUS1
        Y7(K) = IUS0 + RUS1
      ENDDO
cmbike... J=2,...,MOVER2
      DO J = 2, MOVER2
        FOLD = 2*J .LT. M + 2
        K0 = J
        ANGLE = TWOPI*DBLE(J-1)/DBLE(M8)
        C1 = DCOS(ANGLE)
        S1 = DSIN(ANGLE)
        C2 = C1*C1 - S1*S1
        S2 = S1*C1 + C1*S1
        C3 = C2*C1 - S2*S1
        S3 = S2*C1 + C2*S1
        C4 = C2*C2 - S2*S2
        S4 = S2*C2 + C2*S2
        C5 = C4*C1 - S4*S1
        S5 = S4*C1 + C4*S1
        C6 = C4*C2 - S4*S2
        S6 = S4*C2 + C4*S2
        C7 = C4*C3 - S4*S3
        S7 = S4*C3 + C4*S3
        GO TO 40
   20   CONTINUE
        FOLD = .FALSE.
        K0 = M + 2 - J
        T = (C1+S1)*E
        S1 = (C1-S1)*E
        C1 = T
        T = S2
        S2 = C2
        C2 = T 
        T = (-C3+S3)*E
        S3 = (C3+S3)*E
        C3 = T
        C4 = -C4
        T = -(C5+S5)*E
        S5 = (-C5+S5)*E
        C5 = T
        T = -S6
        S6 = -C6
        C6 = T 
        T = (C7-S7)*E
        S7 = -(C7+S7)*E
        C7 = T
   40   CONTINUE
!$OMP parallel do private(K,RS0,IS0,RU0,IU0,RS1,IS1,RU1,IU1)
!$OMP+ private(RS2,IS2,RU2,IU2,RS3,IS3,RU3,IU3,RSS0,ISS0)
!$OMP+ private(RSU0,ISU0,RSS1,ISS1,RSU1,ISU1,RUS0,IUS0)
!$OMP+ private(RUU0,IUU0,RUS1,IUS1,RUU1,IUU1,T)
!$OMP+ private(R1,I1,R2,I2,R3,I3,R4,I4,R5,I5,R6,I6,R7,I7)
!$OMP+ shared(E)
#ifdef __SR11000
*poption parallel
*poption tlocal(K,RS0,IS0,RU0,IU0,RS1,IS1,RU1,IU1)
*poption tlocal(RS2,IS2,RU2,IU2,RS3,IS3,RU3,IU3,RSS0,ISS0)
*poption tlocal(RSU0,ISU0,RSS1,ISS1,RSU1,ISU1,RUS0,IUS0)
*poption tlocal(RUU0,IUU0,RUS1,IUS1,RUU1,IUU1,T)
*poption tlocal(R1,I1,R2,I2,R3,I3,R4,I4,R5,I5,R6,I6,R7,I7)
#endif
        DO K = K0, PTS, M8
          RS0 = X0(K) + X4(K)
          IS0 = Y0(K) + Y4(K)
          RU0 = X0(K) - X4(K)
          IU0 = Y0(K) - Y4(K)
          RS1 = X1(K) + X5(K)
          IS1 = Y1(K) + Y5(K)
          RU1 = X1(K) - X5(K)
          IU1 = Y1(K) - Y5(K)
          RS2 = X2(K) + X6(K)
          IS2 = Y2(K) + Y6(K)
          RU2 = X2(K) - X6(K)
          IU2 = Y2(K) - Y6(K)
          RS3 = X3(K) + X7(K)
          IS3 = Y3(K) + Y7(K)
          RU3 = X3(K) - X7(K)
          IU3 = Y3(K) - Y7(K)
          RSS0 = RS0 + RS2
          ISS0 = IS0 + IS2
          RSU0 = RS0 - RS2
          ISU0 = IS0 - IS2
          RSS1 = RS1 + RS3
          ISS1 = IS1 + IS3
          RSU1 = RS1 - RS3
          ISU1 = IS1 - IS3
          RUS0 = RU0 - IU2
          IUS0 = IU0 + RU2
          RUU0 = RU0 + IU2
          IUU0 = IU0 - RU2
          RUS1 = RU1 - IU3
          IUS1 = IU1 + RU3
          RUU1 = RU1 + IU3
          IUU1 = IU1 - RU3
          T = (RUS1+IUS1)*E
          IUS1 = (IUS1-RUS1)*E
          RUS1 = T
          T = (RUU1+IUU1)*E
          IUU1 = (IUU1-RUU1)*E
          RUU1 = T
          X0(K) = RSS0 + RSS1
          Y0(K) = ISS0 + ISS1
          R1 = RUU0 + RUU1
          I1 = IUU0 + IUU1
          R2 = RSU0 + ISU1
          I2 = ISU0 - RSU1
          R3 = RUS0 + IUS1
          I3 = IUS0 - RUS1
          R4 = RSS0 - RSS1
          I4 = ISS0 - ISS1
          R5 = RUU0 - RUU1
          I5 = IUU0 - IUU1
          R6 = RSU0 - ISU1
          I6 = ISU0 + RSU1
          R7 = RUS0 - IUS1
          I7 = IUS0 + RUS1
          X4(K) = R1*C1 + I1*S1
          Y4(K) = I1*C1 - R1*S1
          X2(K) = R2*C2 + I2*S2
          Y2(K) = I2*C2 - R2*S2
          X6(K) = R3*C3 + I3*S3
          Y6(K) = I3*C3 - R3*S3
          X1(K) = R4*C4 + I4*S4
          Y1(K) = I4*C4 - R4*S4
          X5(K) = R5*C5 + I5*S5
          Y5(K) = I5*C5 - R5*S5
          X3(K) = R6*C6 + I6*S6
          Y3(K) = I6*C6 - R6*S6
          X7(K) = R7*C7 + I7*S7
          Y7(K) = I7*C7 - R7*S7
        ENDDO
        IF (FOLD) GO TO 20
      ENDDO
C           
      RETURN
      END   
C     ==================================================================
      SUBROUTINE C06ECQ(X,Y,PTS,M,P)
C     ----------------------------------------------------------------==
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX PRIME COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, P, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X(PTS), Y(PTS)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANGLE, IS, IU, RS, RU, T, TWOPI, XT, YT
      INTEGER           J, JJ, K, K0, KS1, KS2, MOVER2, MP, PM, PP, U, V
     .                , JJ2
      LOGICAL           FOLD, ZERO
C     .. Local Arrays ..
      DOUBLE PRECISION  A(18), AA(9,9), B(18), BB(9,9), C(18), IA(9),
     *                  IB(9), RA(9), RB(9), S(18)
C     .. External Functions ..
cmb      DOUBLE PRECISION  X01AAF
c      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         DCOS, DBLE, DSIN
C     .. Executable Statements ..
cmb      TWOPI = 2.0D0*X01AAF(0.0D0)
      TWOPI = 8.0D0*DATAN(1.0D0)
      MOVER2 = M/2 + 1
      MP = M*P 
      PP = P/2 
      PM = P - 1
!$OMP parallel do private(U,V,JJ,JJ2,ANGLE)
      DO U = 1, PP
         JJ = P - U
         ANGLE = TWOPI*DBLE(U)/DBLE(P)
         A(U) = DCOS(ANGLE)
         B(U) = DSIN(ANGLE)
         A(JJ) = A(U)
         B(JJ) = -B(U)
         DO V = 1, PP
            JJ2 = U*V - ((U*V)/P)*P
            AA(V,U) = A(JJ2)
            BB(V,U) = B(JJ2)
         ENDDO
      ENDDO
cmb      DO U = 1, PP
cmb         DO V = 1, PP
cmb            JJ = U*V - ((U*V)/P)*P
cmb            AA(V,U) = A(JJ)
cmb            BB(V,U) = B(JJ)
cmb         ENDDO
cmb      ENDDO
C
      DO 300 J = 1, MOVER2
         FOLD = J .GT. 1 .AND. 2*J .LT. M + 2
         K0 = J
         ANGLE = TWOPI*DBLE(J-1)/DBLE(MP)
         ZERO = ANGLE .EQ. 0.0D0
         C(1) = DCOS(ANGLE)
         S(1) = DSIN(ANGLE)
!$OMP parallel do private(U)
         DO U = 2, PM
            C(U) = C(U-1)*C(1) - S(U-1)*S(1)
            S(U) = S(U-1)*C(1) + C(U-1)*S(1)
         ENDDO
         GO TO 140
  100    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
!$OMP parallel do private(U,T)
         DO U = 1, PM
            T = C(U)*A(U) + S(U)*B(U)
            S(U) = -S(U)*A(U) + C(U)*B(U)
            C(U) = T
         ENDDO
  140    CONTINUE
C
         DO 280 K = K0, PTS, MP
            XT = X(K)
            YT = Y(K)
            KS1 = M + K
            KS2 = (P-1)*M + K
            RS = X(KS1) + X(KS2)
            IS = Y(KS1) + Y(KS2)
            RU = X(KS1) - X(KS2)
            IU = Y(KS1) - Y(KS2)
!$OMP parallel do private(U) shared(XT,YT,RS,IS,RU,IU)
#ifdef __SR11000
*poption parallel, tlocal(U)
#endif
            DO U = 1, PP
               RA(U) = XT + RS*AA(U,1)
               IA(U) = YT + IS*AA(U,1)
               RB(U) = RU*BB(U,1)
               IB(U) = IU*BB(U,1)
            ENDDO
            XT = XT + RS
            YT = YT + IS
            DO U = 2, PP
               JJ = P - U
               KS1 = U*M + K
               KS2 = JJ*M + K
               RS = X(KS1) + X(KS2)
               IS = Y(KS1) + Y(KS2)
               RU = X(KS1) - X(KS2)
               IU = Y(KS1) - Y(KS2)
               XT = XT + RS
               YT = YT + IS
!$OMP parallel do private(V) shared(RS,IS,RU,IU)
#ifdef __SR11000
*poption parallel, tlocal(V)
#endif
               DO V = 1, PP
                  RA(V) = RA(V) + RS*AA(V,U)
                  IA(V) = IA(V) + IS*AA(V,U)
                  RB(V) = RB(V) + RU*BB(V,U)
                  IB(V) = IB(V) + IU*BB(V,U)
               ENDDO
            ENDDO
            X(K) = XT 
            Y(K) = YT
            IF (ZERO) THEN  ! cmb
             DO U = 1, PP
               JJ = P - U
               KS1 = U*M + K
               X(KS1) = RA(U) + IB(U)
               Y(KS1) = IA(U) - RB(U)
               KS1 = JJ*M + K
               X(KS1) = RA(U) - IB(U)
               Y(KS1) = IA(U) + RB(U)
             ENDDO 
            ELSE  ! cmb
             DO U = 1, PP
               JJ = P - U
               XT = RA(U) + IB(U)
               YT = IA(U) - RB(U)
               KS1 = U*M + K
               X(KS1) = XT*C(U) + YT*S(U)
               Y(KS1) = YT*C(U) - XT*S(U)
               XT = RA(U) - IB(U)
               YT = IA(U) + RB(U)
               KS1 = JJ*M + K
               X(KS1) = XT*C(JJ) + YT*S(JJ)
               Y(KS1) = YT*C(JJ) - XT*S(JJ)
             ENDDO  ! cmb
            ENDIF
  280    CONTINUE
         IF (FOLD) GO TO 100
  300 CONTINUE
C           
      RETURN
      END   
C     ==================================================================
      INTEGER FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C     ==--------------------------------------------------------------==
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-621 (APR 1988).
C     MARK 13B REVISED. IER-668 (AUG 1988).
C
C     P01ABF is the error-handling routine for the NAG Library.
C
C     P01ABF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ABF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ABF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR
      CHARACTER*72            MESS
C     .. External Subroutines ..
c      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   20       CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     *  ' =',I6)
      END
C     ==================================================================
      SUBROUTINE P01ABZ
C     ==--------------------------------------------------------------==
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     Terminates execution when a hard failure occurs.
C
C     ******************** IMPLEMENTATION NOTE ********************
C     The following STOP statement may be replaced by a call to an
C     implementation-dependent routine to display a message and/or
C     to abort the program.
C     *************************************************************
C     .. Executable Statements ..
      WRITE(6,'(''!!! Execution terminated due to hard failure !!!'')')
      WRITE(6,'(''Program stopped at P01ABZ in c06xxx.F'')')
      CALL STOPGM('P01ABZ','FATAL ERROR')
      RETURN
      END
C     ==================================================================
      SUBROUTINE X04AAF(I,NERR)
C     ==--------------------------------------------------------------==
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-829 (DEC 1989).
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     .. Scalar Arguments ..
      INTEGER           I, NERR
C     .. Local Scalars ..
      INTEGER           NERR1
C     .. Save statement ..
      SAVE              NERR1
C     .. Data statements ..
      DATA              NERR1/0/
C     .. Executable Statements ..
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
C     ==================================================================
      SUBROUTINE X04BAF(NOUT,REC)
C     ==--------------------------------------------------------------==
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAF writes the contents of REC to the unit defined by NOUT.
C
C     Trailing blanks are not output, except that if REC is entirely
C     blank, a single blank character is output.
C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
C     then no output occurs.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*(*)     REC
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (NOUT.GE.0) THEN
C        Remove trailing blanks 
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
C        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
C
99999 FORMAT (A)
      END
C     ==================================================================
cmb      DOUBLE PRECISION FUNCTION X01AAF(X)
C     ==--------------------------------------------------------------==
C     MARK 8 RE-ISSUE. NAG COPYRIGHT 1980.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     RETURNS THE VALUE OF THE MATHEMATICAL CONSTANT PI.
C
C     X IS A DUMMY ARGUMENT
C
C     IT MAY BE NECESSARY TO ROUND THE real*8 CONSTANT IN THE
C     ASSIGNMENT STATEMENT TO A SMALLER NUMBER OF SIGNIFICANT
C     DIGITS IN ORDER TO AVOID COMPILATION PROBLEMS.  IF SO, THEN
C     THE NUMBER OF DIGITS RETAINED SHOULD NOT BE LESS THAN
C     .     2 + INT(FLOAT(IT)*ALOG10(IB))
C     WHERE  IB  IS THE BASE FOR THE REPRESENTATION OF FLOATING-
C     .             POINT NUMBERS
C     . AND  IT  IS THE NUMBER OF IB-ARY DIGITS IN THE MANTISSA OF
C     .             A FLOATING-POINT NUMBER.
C
C     .. Scalar Arguments ..
cmb      DOUBLE PRECISION                 X
C     .. Executable Statements ..
cmb      X01AAF = 3.14159265358979323846264338328D0
cmb      RETURN
cmb      END
C     ==================================================================
