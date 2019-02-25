C hello emacs, this is -*- fortran -*-
#ifdef EWALD
COMMSUBR PRPVIR
C     SUBROUTINE PRPVIR2(NATTOT,NPM,NSM,X,V,XR,XCEN,TMASS,
C    $                   NSPM,NSP,SUBMAS,EKCM,EKCMTO,LEVERY)
C
C     PRPVIR does some preparatory calculations for calculating
C     the virial.
C     Essentially, what we do in this routine is calculate the
C     centre of masses for later use.
C Phil
C     XR(1..NATTOT) delivered with atom coordinates with respect to the center
C                   of mass
C     XCEN(1..3*(NPM*NSPM+NSM)) delivered with center of mass of solute and solvent
C                   (sub)molecules in the central unit cell
Cmb - Revised on 21 June 2006
COMMEND

      SUBROUTINE PRPVIR2(NATTOT,NPM,NSM,X,V,XR,XCEN,TMASS,
     .                   NSPM,NSP,SUBMAS,EKCM,EKCMTO,LEVERY)
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'box.h'
      INCLUDE 'cenmas.h'
C args
      INTEGER NATTOT,NPM,NSM,NSPM
      INTEGER NSP(NSPM)
      LOGICAL LEVERY
      real*8 TMASS,SUBMAS(NSPM),X(NATTOT*NDIM),V(NATTOT*NDIM)
      real*8 XR(NATTOT*NDIM)
      real*8 EKCM(NDIM),EKCMTO
C local vars
      LOGICAL LSCLCM
      INTEGER NATTO3,N,M,NSUB,NFIRST,NLAST,NUMAT,NNP
      INTEGER NSKIPA,NSKPA3,NSKIPT, NSUBTO
      real*8 XREF(MAXDIM),TM
      real*8 XCM(MAXDIM),VCM(MAXDIM),ACMDUM(MAXDIM),
     c    OCMDUM(MAXDIM),EKRDUM
      real*8 ESCLCM,DSTTMP
C Phil
      real*8 XCEN(3*(NPM*NSPM+NSM))
C begin
      LSCLCM = .FALSE.
      EKCMTO = 0.0D0
      ESCLCM = 0.0D0

      NATTO3 = NDIM*NATTOT
      CALL mm_AZZERO(EKCM,NDIM)

!$OMP parallel do private(N)
#ifdef __SR11000
*poption parallel, tlocal(N)
#endif
      DO N=1,NATTO3
         XR(N) = X(N)
      ENDDO

C NSKIPA : number of atoms to skip in the coord array
C NSKIPT : number of atoms to skip in the topology(WINV array)
      NSKIPA = 0
      NSKPA3 = 0
C Phil
      NSUBTO = 1
C done Phil
      DO 30 NNP=1,NPM
         NFIRST = 1
         NSKIPT = 0
         DO 40 NSUB=1,NSPM
            NLAST = NSP(NSUB)
            NUMAT = NLAST - NFIRST + 1
C XREF is always the first atom in the current submolecule
            DO M=1,NDIM
               XREF(M) = XR(NSKPA3+M)
            ENDDO
            CALL CONAT(NATTOT,NUMAT,NSKIPA,XREF,XR,LEVERY)
            TM = SUBMAS(NSUB)
!           CALL CENMAS(NATTOT,1,NUMAT,0,0,NSKIPA,NDIM,NDRMAX,XR,V,TM,
!    $           NSKIPT,WMAS,WMASS,
!    $           EKCMTO,XCM,VCM,ACMDUM,EKRDUM,OCMDUM, -ICMVEL)
C
!$OMP parallel do private(M)
            DO M=1,NDIM
               EKCM(M) = EKCM(M) + TM*VCM(M)**2
            ENDDO
            DO N=1,NUMAT
               DO M=1, NDIM
                  XR(NSKPA3+M) = XR(NSKPA3+M) - XCM(M)
               ENDDO
               NSKPA3 = NSKPA3 + NDIM
            ENDDO
C Phil
            DO M=1,3
              DSTTMP = XCM(M)
 1000         CONTINUE
              IF (DSTTMP .GE. BOX(M)) THEN
                DSTTMP = DSTTMP - BOX(M)
                GOTO 1000
C              ELSEIF (DSTTMP .LT. -BOX(M)) THEN
C                DSTTMP = DSTTMP + BOX(M)
C Phil 12/19/97
              ELSEIF (DSTTMP .LT. 0.0) THEN
                DSTTMP = DSTTMP + BOX(M)
                GOTO 1000
              ENDIF
C
              if (dsttmp.gt.box(m).or.dsttmp.lt.0.0) then
                print *,'prpvir error solute'
                stop
              endif
C
              XCM(M) = DSTTMP
C
            ENDDO
C
            XCEN(3*NSUBTO-2) = XCM(1)
            XCEN(3*NSUBTO-1) = XCM(2)
            XCEN(3*NSUBTO  ) = XCM(3)
            NSUBTO = NSUBTO + 1
C done Phil
            NSKIPT = NSKIPT + NUMAT
            NSKIPA = NSKIPA + NUMAT
            NFIRST = NLAST + 1
 40      CONTINUE
 30   CONTINUE

C center of mass heat bath: prepare ESCLCM
      IF (LSCLCM) THEN
!$OMP parallel do private(M) reduction(+:ESCLCM)
         DO M = 1,NDIM
            ESCLCM = ESCLCM + EKCM(M)
         ENDDO
      ENDIF

C do solvent
C NSKIPT stays 0
      NSKIPT = 0
      NSKIPA = NPM*NRP
      NSKPA3 = NDIM*NSKIPA
      DO 100 NNP=1,NSM
!        CALL CENMAS(NATTOT,0,0,1,NRAM,NSKIPA,NDIM,NDRMAX,XR,V,
!    $        TMASS,NSKIPT,WMAS,WMASS,
!    $        EKCMTO,XCM,VCM,ACMDUM,EKRDUM,OCMDUM, -ICMVEL)
cmb - the above routine (CENMAS) is no longer used anywhere !
!$OMP parallel do private(M)
         DO M=1,NDIM
            EKCM(M) = EKCM(M) + TMASS*VCM(M)**2
         ENDDO
         DO N=1,NRAM
            DO M=1,NDIM
               XR(NSKPA3+M) = XR(NSKPA3+M) - XCM(M)
            ENDDO
            NSKPA3 = NSKPA3 + NDIM
         ENDDO
C Phil
         DO M=1,3
           DSTTMP = XCM(M)
 1100      CONTINUE
           IF (DSTTMP .GE. BOX(M)) THEN
             DSTTMP = DSTTMP - BOX(M)
             GOTO 1100
           ELSEIF (DSTTMP .LT. 0.0) THEN
             DSTTMP = DSTTMP + BOX(M)
             GOTO 1100
           ENDIF
           XCM(M) = DSTTMP
C
           IF (dsttmp.gt.box(m).or.dsttmp.lt.0.0d0) THEN
             print *,'inc_runmd2: PRPVIR ERROR SOLVENT'
             stop
           ENDIF
         ENDDO
C
         XCEN(3*NSUBTO-2) = XCM(1)
         XCEN(3*NSUBTO-1) = XCM(2)
         XCEN(3*NSUBTO  ) = XCM(3)
         NSUBTO = NSUBTO + 1
C done Phil
         NSKIPA = NSKIPA + NRAM
 100  CONTINUE
C end PRPVIR2
      END
#endif
