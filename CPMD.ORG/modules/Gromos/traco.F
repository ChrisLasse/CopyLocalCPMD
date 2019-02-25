C     ==================================================================
      SUBROUTINE TRACO(NR,NSKIP,X,BBETA,IT,LNEW)
C     ==--------------------------------------------------------------==
C     W.F. VAN GUNSTEREN, GRONINGEN, DEC. 1985
C
C     SUBROUTINE TRACO (NR,NPA,X,BETA,IT,LNEW)
C
C     TRACO PERFORMS THE TRANSFORMATION FROM CARTESIAN COORDINATES
C     TO OBLIQUE (MONOCLINIC) CONTRAVARIANT COORDINATES OR VICE VERSA.
C
C     NR = NUMBER OF ATOMS FOR WHICH THE TRANSFORMATION IS TO BE
C          PERFORMED (NR>0)
C     NPA = NUMBER OF ATOMS (IN X) PRECEEDING THE FIRST ATOM FOR WHICH
C           THE TRANSFORMATION IS TO BE PERFORMED
C     X(3*NPA+1..3*(NPA+NR)) = COORDINATES TO BE TRANSFORMED,
C                              DELIVERED WITH TRANSFORMED COORDINATES
C     BETA = ANGLE BETWEEN THE X- AND Z-AXES OF THE OBLIQUE COORDINATE
C            SYSTEM (DEGREES, 0<BETA<180)
C            THE Y-AXIS IS ASSUMED TO BE ORTHOGONAL TO THE OTHER ONES
C     IT > 0 : TRANSFORMATION FROM CARTESIAN TO OBLIQUE COORDINATES
C        = 0 : NO TRANSFORMATION IS PERFORMED
C        < 0 : TRANSFORMATION FROM OBLIQUE TO CARTESIAN COORDINATES
C
C     LNEW = .FALSE. :
C                QUANTITIES DEPENDING ON BETA ARE ONLY CALCULATED AT
C                THE FIRST SUBR. CALL
C          = .TRUE.  : THEY ARE CALCULATED AT EVERY SUBR. CALL
C     ==--------------------------------------------------------------==
Cmb - Revised on 15 Sept. 2005
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
C args
      INTEGER NR,NSKIP,IT
      real*8 BBETA,X(NDIM*(NR+NSKIP))
      LOGICAL LNEW
C local vars
      LOGICAL LFIRST
      SAVE LFIRST
      INTEGER I3,J
      real*8 BETAR,COSB,SINB,SINBIN,XH
      SAVE BETAR,COSB,SINB,SINBIN
C data
      DATA LFIRST/.TRUE./
C begin
      IF (LNEW .OR. LFIRST) THEN
         LFIRST = .FALSE.
         BETAR = BBETA*DATAN(1.0D0)/45.d0
         COSB = DCOS(BETAR)
         SINB = DSIN(BETAR)
         SINBIN = 1.0D0/SINB
      ENDIF

      IF (IT .GT. 0) THEN
         I3 = NDIM*NSKIP + 1
         DO J = 1,NR
            XH = X(I3+2)*SINBIN
            X(I3+2) = XH
            X(I3) = X(I3) - COSB*XH
            I3 = I3 + NDIM
         ENDDO
      ELSEIF(IT .LT. 0) THEN
         I3 = NDIM*NSKIP + 1
         DO J = 1,NR
            X(I3) = X(I3) + COSB*X(I3+2)
            X(I3+2) = X(I3+2)*SINB
            I3 = I3 + NDIM
         ENDDO
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
