C     ==================================================================
      SUBROUTINE CONAT(NATTOT,NRP,NSKIP,XREF,X,LEVERY)
C     ==--------------------------------------------------------------==
C     W.F. VAN GUNSTEREN, GRONINGEN, JAN. 1986
Cmb - Revised on 21 June 2006
C     SUBROUTINE CONAT(NATTOT,NRP,NSKIP,XREF,X,LEVERY)
C
C     CONAT gathers the atoms of a  molecule by applying periodic
C     boundary conditions, such that each atom lies within BOX/2 or
C     BOX*SQRT(3)/4 of its predecessor. The first atom of the molecule
C     is compared to the reference position XREF.
C
C     NOTE: ONLY THE FIRST THREE DIMENSIONS ARE TREATED,
C     REGARDLESS OF NDIM.
C
C     The periodic box can be:
C     a truncated octahedron,
C     rectangular (L<BETA> =90 degrees) or
C     monoclinic (L<BETA> .NE. 90 degrees)
C
C     NATTOT the total number of atoms in the X array
C     NRP   number of molecular atoms
C     NSKIP number of atoms in X preceeding the first atom of the
C           molecule
C     XREF(NDIM)
C           cartesian coordinates of the reference position to
C           which the first molecular atom is compared.
C     X(NDIM*NATTOT)
C           atom cartesian coordinates
C     LEVERY
C      .FALSE.
C           quantities depending on BETA are only calculated at
C           the first SUBR. call.
C      .TRUE.
C           they are calculated every SUBR. call.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
C args
      INTEGER NRP,NSKIP,NATTOT
      LOGICAL LEVERY
      real*8 XREF(NDIM),X(NDIM*NATTOT)
C local vars
      LOGICAL LFIRST,LDOTRA,LMONO,LOCTO,LVAC,LDOVIR
      SAVE LFIRST,LDOTRA,LMONO,LOCTO,LVAC,LDOVIR

      real*8 BOXOH,BOXOQ
      SAVE BOXOH,BOXOQ
C
      INTEGER J,M,I3
      real*8 XAUX,XAUY,XAUZ
      real*8 BETAR,COSB
      real*8 XI(MAXDIM),XIJ(MAXDIM)
C data
      DATA LFIRST/.TRUE./
C begin
      IF (LFIRST .OR. LEVERY) THEN
         LFIRST = .FALSE.

         LMONO = (NTB .GT. 0)
         LOCTO = (NTB. LT. 0)
         LVAC  = (NTB .EQ. 0)
         LDOVIR= (ABS(NTB) .EQ. NTBVIR)

         LDOTRA = .FALSE.
         IF (LMONO)THEN
            BETAR=BETA*DATAN(1.0d+00)/45.0d+00
            COSB=DCOS(BETAR)
            LDOTRA=(DABS(COSB).GE.1.d-04)
         ENDIF

         BOXOH = BOXH(1)
         BOXOQ = BOX(1)*0.75d+00
      ENDIF

      IF (LVAC) THEN
         PRINT *,'CONAT: called with vacuum!'
         RETURN
      ENDIF

      IF (LDOTRA) THEN
         CALL TRACO(1,0,XREF,BETA,1,LEVERY)
         CALL TRACO(NRP,NSKIP,X,BETA,1,LEVERY)
      ENDIF

C gather the atoms
!$OMP parallel do private(M)
#ifdef __SR11000
*poption parallel, tlocal(M)
#endif
      DO M=1,NDIM
         XI(M) = XREF(M)
      ENDDO

      IF (LOCTO) THEN
        I3 = NDIM*NSKIP
        DO J=1,NRP
!$OMP parallel do private(M,XAUX)
          DO M=1,NDRMAX
            XAUX = (X(I3+M)-XI(M))*BOXINV(M)
            XIJ(M) = (XAUX-NINT(XAUX))*BOX(M)
          ENDDO
          XAUZ =DABS(XIJ(1))+DABS(XIJ(2))+DABS(XIJ(3))
          IF (XAUZ .GE. BOXOQ) THEN
!$OMP parallel do private(M)
            DO M=1,NDRMAX
                  XIJ(M) = XIJ(M) - SIGN(BOXOH,XIJ(M))
            ENDDO
          ENDIF
!$OMP parallel do private(M,XAUX)         
          DO M=1,NDRMAX
            XAUX = XI(M) + XIJ(M)
            XI(M) = XAUX
            X(I3+M) = XAUX
          ENDDO
          I3 = I3 + NDIM
        ENDDO
      ELSE
        I3 = NDIM*NSKIP
        DO J=1,NRP
!$OMP parallel do private(M,XAUX,XAUY) shared(I3)
          DO M=1,NDRMAX
            XAUX = (X(I3+M)-XI(M))*BOXINV(M)
            XIJ(M) = (XAUX-NINT(XAUX))*BOX(M)
            XAUY = XI(M) + XIJ(M)
            XI(M) = XAUY
            X(I3+M) = XAUY
          ENDDO
          I3 = I3 + NDIM
        ENDDO
      ENDIF

      IF (LDOTRA) THEN
         CALL TRACO(1,0,XREF,BETA,-1,LEVERY)
         CALL TRACO(NRP,NSKIP,X,BETA,-1,LEVERY)
      ENDIF
C end conat
C     ==--------------------------------------------------------------==
      RETURN
      END
