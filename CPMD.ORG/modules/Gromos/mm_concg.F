C     ==================================================================
      SUBROUTINE mm_SHIAG(NATTOT,NPM,NSM,XMIN,X,LEVERY)
C     ==--------------------------------------------------------------==
C     W.F. VAN GUNSTEREN, GRONINGEN, APR. 1985
C     Id: concg.f,v 1.10 1996/06/26 11:24:41 wscott Exp
C     
C     SUBROUTINE SHIAG(NATTOT,NPM,NSM,XMIN,X,LEVERY)
C
C     SHIAG translates solute charge group atoms and solvent molecules,
C     applying periodic boundary conditions, such that the first atom
C     of a solute charge group or of a solvent molecule lies within
C     the specified periodic box. This box lies in the positive
C     quadrant with respect to an origin XMIN(1..NDIM).
C     It can be a truncated octahedron, rectangular or monoclinic.
C
C     The atoms of a charge group or solvent molecule must lie within
C     BOX/2 or BOX*SQRT(3)/4 of each other.
C
C     NOTE: ONLY THE FIRST THREE DIMENSIONS ARE TREATED, I.E.
C     ATOMS ARE ONLY PLACED BACK INTO THE BOX IN THE FIRST
C     THREE DIMENSIONS REGARDLESS OF L<NDIM>.
C
C     NATTOT the total number of atoms
C     NPM    the number of solute molecules
C     NSM    the number os solvent molecules
C     XMIN(MAXDIM)
C            the origin of the box specified.
C     X(NDIM*NATTOT)
C            the coordinates to be translated.
C            delivered with the translated coordinates
C
C     LEVERY
C      .FALSE.
C           quantities depending on BETA are only calculated at
C           the first SUBR. call.
C      .TRUE.
C           they are calculated every SUBR. call.
C
C     Topology information is taken from common blocks defined in
C     toposz.h and topoar.h.
C     The box information and type of boundary condition is taken
C     from common blocks defined in box.h.
C
C     Rewritten by Walter Scott, Zuerich 1995/6. Special thanks to
C     Ton Linssen (RU Groningen) for pointing out a small error in the
C     truncated octahedron case.
C     ==--------------------------------------------------------------==
cmb - Revised on 13 Oct. 2005
C includes
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
C args
      INTEGER NATTOT,NPM,NSM
      real*8 XMIN(*)
      real*8 X(*)
      LOGICAL LEVERY
C local vars
      LOGICAL LFIRST,LDOTRA,LMONO,LOCTO,LVAC
      SAVE LFIRST,LDOTRA,LMONO,LOCTO,LVAC

      INTEGER I3,NN,NFIRST,NLAST,NATCG,NCG,M,J
      LOGICAL LSHIFT

      real*8 BOXOH,BOXOHH(MAXDIM),BOXOQ
      SAVE BOXOH,BOXOHH,BOXOQ

      real*8 BETAR,COSB
      real*8 XI(MAXDIM),XAUX
      logical change
      integer change_count
C data
      DATA LFIRST/.TRUE./
C begin
      IF (LFIRST .OR. LEVERY) THEN
         LFIRST = .FALSE.

         LMONO = (NTB .GT. 0)
         LOCTO = (NTB. LT. 0)
         LVAC  = (NTB .EQ. 0)

         LDOTRA = .FALSE.
         IF (LMONO)THEN
C we have a monoclinic box
            BETAR = BETA*DATAN(1.0d0)/45.0D0
            COSB = DCOS(BETAR)
            LDOTRA = (DABS(COSB).GE.1.D-4)
         ENDIF

         IF (LOCTO) THEN
            BOXOH = BOXH(1)
            BOXOQ = BOX(1)*0.75d0
!$OMP parallel do private(M)
            DO M=NDRMIN,NDRMAX
               BOXOHH(M)= BOXH(M) + XMIN(M)
            ENDDO
         ELSE
            BOXOH = 0.0d0
            BOXOQ = 0.0d0
!$OMP parallel do private(M)
            DO M=NDRMIN,NDRMAX
               BOXOHH(M) = 0.0d0
            ENDDO
         ENDIF
      ENDIF

      IF (LVAC) THEN
         WRITE(6,*)'SHIAG CALLED IN VACUUM!!'
         RETURN
      ENDIF

      IF (NATTOT .NE. NPM*NRP + NSM*NRAM) THEN
         WRITE(6,*)'SHIAG: inconsistent number of atoms!'
         WRITE(6,*)'NATTOT =',NATTOT,'  but NPM*NRP+NSM*NRAM =',
     .              NPM*NRP+NSM*NRAM
         CALL STOPGM('MM_SHIAG','FATAL ERROR')
      ENDIF

      IF (LDOTRA) THEN
         CALL TRACO(NATTOT,0,X,BETA,1,LEVERY)
      ENDIF

      change_count=0
117   continue
      change=.false.
      change_count=change_count + 1
      if(MOD(change_count,100).eq.0)
     &    write(6,*) 'mm_SHIAG| WARNING!',change_count,' translations'
      if(change_count.gt.500)
     &     CALL STOPGM('mm_SHIAG','too much translation? NaN?')

C translate the solute atoms
C NATCG: the number of atoms in the current CG
      I3 = 0
      DO 10 NN=1,NPM
         NFIRST = 1
         DO 20 NCG=1, NCAG
            NLAST = INC(NCG)
            NATCG = NLAST - NFIRST +1

C see if we need to shift first atom
            LSHIFT = .FALSE.
            DO M=NDRMIN,NDRMAX
               XAUX = X(I3+M)
               XI(M)= XAUX
               IF (XAUX .GE. XMIN(M)+BOX(M)) THEN
                  X(I3+M) = XAUX - BOX(M)
                  change=.true.
                  LSHIFT = .TRUE.
               ELSEIF (XAUX .LT. XMIN(M)) THEN
                  X(I3+M) = XAUX + BOX(M)
                  change=.true.
                  LSHIFT = .TRUE.
               ENDIF
            ENDDO

            IF (LOCTO) THEN
               XAUX = DABS(X(I3+1)-BOXOHH(1))
     .              + DABS(X(I3+2)-BOXOHH(2))
     .              + DABS(X(I3+3)-BOXOHH(3))
               IF (XAUX .GE. BOXOQ) THEN
                  LSHIFT = .TRUE.
!$OMP parallel do private(M)
                  DO M=NDRMIN,NDRMAX
                     X(I3+M) = X(I3+M) - SIGN(BOXOH,X(I3+M)-BOXOHH(M))
                     change=.true.
                  ENDDO
               ENDIF
            ENDIF

C shift other atoms of this CG if we have to and
C if there are any 
            IF (NATCG .GT. 1 .AND. LSHIFT) THEN
C XI contains the original (unshifted) coords of first atom
C X(I3) contains the new coords of first atom
C now calculate 'shift vector' into XI
!$OMP parallel do private(M)
               DO M=NDRMIN,NDRMAX
                  XI(M) = X(I3+M) - XI(M)
               ENDDO
C add the shift vector to all other atoms of this cg
               I3 = I3+NDIM
               DO J=2,NATCG
!$OMP parallel do private(M)
                  DO M=NDRMIN,NDRMAX
                     X(I3+M) = X(I3+M) + XI(M)
                     change=.true.
                  ENDDO
                  I3 = I3 + NDIM
               ENDDO
            ELSE
               I3 = I3 + NDIM*NATCG
            ENDIF
            NFIRST = NLAST + 1
 20      CONTINUE
 10   CONTINUE

C do solvent molecules
      DO 100 NN=1,NSM
C see if we need to shift first atom
         LSHIFT = .FALSE.
         DO M=NDRMIN,NDRMAX
            XAUX = X(I3+M)
            XI(M)= XAUX
            IF (XAUX .GE. XMIN(M)+BOX(M)) THEN
               X(I3+M) = XAUX - BOX(M)
               change=.true.
               LSHIFT = .TRUE.
            ELSEIF (XAUX .LT. XMIN(M)) THEN
               X(I3+M) = XAUX + BOX(M)
               change=.true.
               LSHIFT = .TRUE.
            ENDIF
         ENDDO

         IF (LOCTO) THEN
            XAUX = DABS(X(I3+1)-BOXOHH(1))
     .           + DABS(X(I3+2)-BOXOHH(2))
     .           + DABS(X(I3+3)-BOXOHH(3))
            IF (XAUX .GE. BOXOQ) THEN
               LSHIFT = .TRUE.
!$OMP parallel do private(M)
               DO M=NDRMIN,NDRMAX
                  X(I3+M) = X(I3+M) - SIGN(BOXOH,X(I3+M)-BOXOHH(M))
                  change=.true.
               ENDDO
            ENDIF
         ENDIF

C shift other atoms of this solvent molecule if we have to and
C if there are any 
         IF (NRAM .GT. 1 .AND. LSHIFT) THEN
C XI contains the original (unshifted) coords of first atom
C X(I3) contains the new coords of first atom
C now calculate 'shift vector' into XI
!$OMP parallel do private(M)
            DO M=NDRMIN,NDRMAX
               XI(M) = X(I3+M) - XI(M)
            ENDDO
C add the shift vector to all other atoms of this cg
            I3 = I3 + NDIM
            DO J=2,NRAM
!$OMP parallel do private(M)
               DO M=NDRMIN,NDRMAX
                  X(I3+M) = X(I3+M) + XI(M)
                  change=.true.
               ENDDO
               I3 = I3 + NDIM
            ENDDO
         ELSE
            I3 = I3 + NDIM*NRAM
         ENDIF
 100  CONTINUE

      if(change) goto 117

      IF (LDOTRA) THEN
         CALL TRACO(NATTOT,0,X,BETA,-1,LEVERY)
      ENDIF
C end shiag
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE CONCG(NRPT,NPM,X,LEVERY)
C     ==--------------------------------------------------------------==
COMMSUBR CONCG
C     SUBROUTINE CONCG(NRPT,NPM,X,LEVERY)
C     
C     CONCG gathers the atoms forming a charge group by applying
C     periodic boundary conditions, such that the atoms of one
C     charge group lie within BOX/2 or BOX*SQRT(3)/4 of each other.
C     CONCG does this for a number of solute molecules.
C     The charge group information is taken from topology.
C
C     The periodic box can be a truncated octahedron, rectangular
C     (L<BETA> = 90.0 degrees) or monoclinic (L<BETA> .NE. 90.0 degrees)
C
C     NOTE: CONCG only performs its coordinate transformations
C     in three dimensions, even if L<NDIM> = 4.
C
C     NRPT: total number of solute atoms
C     NPM: number of solute molecules
C     X(NDIM*NRPT): atom coordinates
C     LEVERY
C      .FALSE.
C           quantities depending on BETA are only calculated at
C           the first SUBR. call.
C      .TRUE.
C           they are calculated avery SUBR. call.
COMMEND
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
C args
      INTEGER NRPT,NPM
      LOGICAL LEVERY
      real*8 X(NRPT*NDIM)
C local vars
      LOGICAL LFIRST,LDOTRA,LMONO,LOCTO,LVAC
      SAVE LFIRST,LDOTRA,LMONO,LOCTO,LVAC

      INTEGER I3,NN,NFIRST,NLAST,NATCG,NCG,M,J
      LOGICAL LSHIFT

      real*8 BOXOH,BOXOQ
      SAVE BOXOH,BOXOQ

      real*8 BETAR,COSB
      real*8 XI(MAXDIM),XIJ(MAXDIM),XAUX
C data
      DATA LFIRST/.TRUE./
C begin
      IF (LFIRST .OR. LEVERY) THEN
         LFIRST = .FALSE.

         LMONO = (NTB .GT. 0)
         LOCTO = (NTB. LT. 0)
         LVAC  = (NTB .EQ. 0)

         LDOTRA = .FALSE.
         IF (LMONO)THEN
C we have a monoclinic box
            BETAR = BETA*DATAN(1.0d0)/45.0D0
            COSB = DCOS(BETAR)
            LDOTRA = (DABS(COSB).GE.1.D-4)
         ENDIF

         IF (LOCTO) THEN
            BOXOH = BOXH(1)
            BOXOQ = BOX(1)*0.75D0
         ELSE
            BOXOH = 0.0D0
            BOXOQ = 0.0D0
         ENDIF
      ENDIF

      IF (LVAC) THEN
         WRITE(6,*)'CONCG CALLED IN VACUUM!!'
         RETURN
      ENDIF

      IF (LDOTRA) THEN
         CALL TRACO(NRPT,0,X,BETA,1,LEVERY)
      ENDIF

      I3 = 0
      DO NN=1,NPM
         NFIRST = 1
         DO NCG = 1,NCAG
            NLAST = INC(NCG)
            NATCG = NLAST - NFIRST +1
            IF (NATCG .GT. 1) THEN
!$OMP parallel do private(M)
               DO M = NDRMIN,NDRMAX
                  XI(M) = X(I3+M)
               ENDDO
               I3 = I3 + NDIM
C see if any neighbours are too far away
               DO J = 2,NATCG
                  LSHIFT = .FALSE.
                  DO M = NDRMIN,NDRMAX
                     XAUX = X(I3+M) - XI(M)
                     IF (XAUX .GT. BOXH(M)) THEN
                        XAUX = XAUX - BOX(M)
                        LSHIFT = .TRUE.
                     ELSEIF (XAUX .LT. -BOXH(M)) THEN
                        XAUX = XAUX + BOX(M)
                        LSHIFT = .TRUE.
                     ENDIF
                     XIJ(M) = XAUX
                  ENDDO
                  IF (LOCTO) THEN
                     IF (DABS(XIJ(1))+DABS(XIJ(2))+DABS(XIJ(3))
     .                  .GE.BOXOQ) THEN
!$OMP parallel do private(M)
                        DO M = NDRMIN,NDRMAX
                           XIJ(M) = XIJ(M) - SIGN(BOXOH,XIJ(M))
                        ENDDO
                        LSHIFT = .TRUE.
                     ENDIF
                  ENDIF
C if distance XIJ is too great, shift atom J
                  IF (LSHIFT) THEN
!$OMP parallel do private(M)
                     DO M = NDRMIN,NDRMAX
                        X(I3+M) = XI(M) + XIJ(M)
                     ENDDO
                  ENDIF
                  I3 = I3 + NDIM
               ENDDO
            ELSE
               I3 = I3 + NDIM
            ENDIF
            NFIRST = NLAST + 1
         ENDDO
      ENDDO

      IF (LDOTRA) THEN
         CALL TRACO(NRPT,0,X,BETA,-1,LEVERY)
      ENDIF
C end concg
      END
C     ==================================================================
      SUBROUTINE SHIAG(NATTOT,NPM,NSM,XMIN,X,LEVERY)
C     ==--------------------------------------------------------------==
COMMSUBR SHIAG
C     W.F. VAN GUNSTEREN, GRONINGEN, APR. 1985
C     
C     SUBROUTINE SHIAG(NATTOT,NPM,NSM,XMIN,X,LEVERY)
C
C     SHIAG translates solute charge group atoms and solvent molecules,
C     applying periodic boundary conditions, such that the first atom
C     of a solute charge group or of a solvent molecule lies within
C     the specified periodic box. This box lies in the positive
C     quadrant with respect to an origin XMIN(1..NDIM).
C     It can be a truncated octahedron, rectangular or monoclinic.
C
C     The atoms of a charge group or solvent molecule must lie within
C     BOX/2 or BOX*SQRT(3)/4 of each other.
C
C     NOTE: ONLY THE FIRST THREE DIMENSIONS ARE TREATED, I.E.
C     ATOMS ARE ONLY PLACED BACK INTO THE BOX IN THE FIRST
C     THREE DIMENSIONS REGARDLESS OF L<NDIM>.
C
C     NATTOT the total number of atoms
C     NPM    the number of solute molecules
C     NSM    the number os solvent molecules
C     XMIN(MAXDIM)
C            the origin of the box specified.
C     X(NDIM*NATTOT)
C            the coordinates to be translated.
C            delivered with the translated coordinates
C
C     LEVERY
C      .FALSE.
C           quantities depending on BETA are only calculated at
C           the first SUBR. call.
C      .TRUE.
C           they are calculated every SUBR. call.
C
C     Topology information is taken from common blocks defined in
C     toposz.h and topoar.h.
C     The box information and type of boundary condition is taken
C     from common blocks defined in box.h.
C
C     Rewritten by Walter Scott, Zuerich 1995/6. Special thanks to
C     Ton Linssen (RU Groningen) for pointing out a small error in the
C     truncated octahedron case.
COMMEND
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
C args
      INTEGER NATTOT,NPM,NSM
      real*8 XMIN(MAXDIM)
      real*8 X(NATTOT*NDIM)
      LOGICAL LEVERY
C local vars
      LOGICAL LFIRST,LDOTRA,LMONO,LOCTO,LVAC
      SAVE LFIRST,LDOTRA,LMONO,LOCTO,LVAC

      INTEGER I3,NN,NFIRST,NLAST,NATCG,NCG,M,J
      LOGICAL LSHIFT

      real*8 BOXOH,BOXOHH(MAXDIM),BOXOQ
      SAVE BOXOH,BOXOHH,BOXOQ

      real*8 BETAR,COSB
      real*8 XI(MAXDIM),XAUX
C data
      DATA LFIRST/.TRUE./
C begin
      IF (LFIRST .OR. LEVERY) THEN
         LFIRST = .FALSE.

         LMONO = (NTB .GT. 0)
         LOCTO = (NTB. LT. 0)
         LVAC  = (NTB .EQ. 0)

         LDOTRA = .FALSE.
         IF (LMONO)THEN
C we have a monoclinic box
            BETAR = BETA*DATAN(1.0d0)/45.0d0
            COSB = DCOS(BETAR)
            LDOTRA = (DABS(COSB).GE.1.D-4)
         ENDIF

         IF (LOCTO) THEN
            BOXOH = BOXH(1)
            BOXOQ = BOX(1)*0.75D0
!$OMP parallel do private(M)
            DO M=NDRMIN,NDRMAX
               BOXOHH(M)= BOXH(M) + XMIN(M)
            ENDDO
         ELSE
            BOXOH = 0.0D0
            BOXOQ = 0.0D0
!$OMP parallel do private(M)
            DO M=NDRMIN,NDRMAX
               BOXOHH(M) = 0.0D0
            ENDDO
         ENDIF
      ENDIF

      IF (LVAC) THEN
         WRITE(6,*)'SHIAG CALLED IN VACUUM!!'
         RETURN
      ENDIF

      IF (NATTOT .NE. NPM*NRP + NSM*NRAM) THEN
         WRITE(6,*)'SHIAG: inconsistent number of atoms!'
         WRITE(6,*)'NATTOT =',NATTOT,' but NPM*NRP+NSM*NRAM =',
     .              NPM*NRP+NSM*NRAM
         CALL STOPGM('MM_SHIAG','FATAL ERROR')
      ENDIF

      IF (LDOTRA) THEN
         CALL TRACO(NATTOT,0,X,BETA,1,LEVERY)
      ENDIF

C translate the solute atoms
C NATCG: the number of atoms in the current CG
      I3 = 0
      DO NN=1,NPM
         NFIRST = 1
         DO NCG=1, NCAG
            NLAST = INC(NCG)
            NATCG = NLAST - NFIRST +1
C see if we need to shift first atom
            LSHIFT = .FALSE.
            DO M=NDRMIN,NDRMAX
               XAUX = X(I3+M)
               XI(M)= XAUX
               IF (XAUX .GE. XMIN(M)+BOX(M)) THEN
                  X(I3+M) = XAUX - BOX(M)
                  LSHIFT = .TRUE.
               ELSEIF (XAUX .LT. XMIN(M)) THEN
                  X(I3+M) = XAUX + BOX(M)
                  LSHIFT = .TRUE.
               ENDIF
            ENDDO
            IF (LOCTO) THEN
               XAUX = DABS(X(I3+1)-BOXOHH(1))
     .              + DABS(X(I3+2)-BOXOHH(2))
     .              + DABS(X(I3+3)-BOXOHH(3))
               IF (XAUX .GE. BOXOQ) THEN
                  LSHIFT = .TRUE.
!$OMP parallel do private(M)
                  DO M=NDRMIN,NDRMAX
                     X(I3+M) = X(I3+M) - SIGN(BOXOH,X(I3+M)-BOXOHH(M))
                  ENDDO
               ENDIF
            ENDIF
C shift other atoms of this CG if we have to and
C if there are any 
            IF (NATCG .GT. 1 .AND. LSHIFT) THEN
C XI contains the original (unshifted) coords of first atom
C X(I3) contains the new coords of first atom
C now calculate 'shift vector' into XI
!$OMP parallel do private(M)
               DO M=NDRMIN,NDRMAX
                  XI(M) = X(I3+M) - XI(M)
               ENDDO
C add the shift vector to all other atoms of this cg
               I3 = I3+NDIM
               DO J=2,NATCG
!$OMP parallel do private(M)
                  DO M=NDRMIN,NDRMAX
                     X(I3+M) = X(I3+M) + XI(M)
                  ENDDO
                  I3 = I3 + NDIM
               ENDDO
            ELSE
               I3 = I3 + NDIM*NATCG
            ENDIF
            NFIRST = NLAST + 1
         ENDDO
      ENDDO

C do solvent molecules
      DO NN=1,NSM
C see if we need to shift first atom
         LSHIFT = .FALSE.
         DO M=NDRMIN,NDRMAX
            XAUX = X(I3+M)
            XI(M)= XAUX
            IF (XAUX .GE. XMIN(M)+BOX(M)) THEN
               X(I3+M) = XAUX - BOX(M)
               LSHIFT = .TRUE.
            ELSEIF (XAUX .LT. XMIN(M)) THEN
               X(I3+M) = XAUX + BOX(M)
               LSHIFT = .TRUE.
            ENDIF
         ENDDO
         IF (LOCTO) THEN
            XAUX = DABS(X(I3+1)-BOXOHH(1))
     .           + DABS(X(I3+2)-BOXOHH(2))
     .           + DABS(X(I3+3)-BOXOHH(3))
            IF (XAUX .GE. BOXOQ) THEN
               LSHIFT = .TRUE.
!$OMP parallel do private(M)
               DO M=NDRMIN,NDRMAX
                  X(I3+M) = X(I3+M) - SIGN(BOXOH,X(I3+M)-BOXOHH(M))
               ENDDO
            ENDIF
         ENDIF
C shift other atoms of this solvent molecule if we have to and
C if there are any 
         IF (NRAM .GT. 1 .AND. LSHIFT) THEN
C XI contains the original (unshifted) coords of first atom
C X(I3) contains the new coords of first atom
C now calculate 'shift vector' into XI
!$OMP parallel do private(M)
            DO M=NDRMIN,NDRMAX
               XI(M) = X(I3+M) - XI(M)
            ENDDO
C add the shift vector to all other atoms of this cg
            I3 = I3 + NDIM
            DO J=2,NRAM
!$OMP parallel do private(M)
               DO M=NDRMIN,NDRMAX
                  X(I3+M) = X(I3+M) + XI(M)
               ENDDO
               I3 = I3 + NDIM
            ENDDO
         ELSE
            I3 = I3 + NDIM*NRAM
         ENDIF
      ENDDO

      IF (LDOTRA) THEN
         CALL TRACO(NATTOT,0,X,BETA,-1,LEVERY)
      ENDIF
C end shiag
      END
