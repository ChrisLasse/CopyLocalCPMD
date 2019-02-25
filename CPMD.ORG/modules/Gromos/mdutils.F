C Id: mdutils.f,v 1.7 1996/11/06 19:57:19 wscott Exp 
C a collection of utility and checking routines used in PROMD and PROFEE.
Cmb - Revised on 28 January 2013
#if defined(__WITH_OMP3)
#define __COLLAPSE2 collapse(2)
#else
#define __COLLAPSE2
#endif
C     ==--------------------------------------------------------------==
      SUBROUTINE GTCOOR(PRGSTR,NR,NDSZE,DLSUM,DMSUM)
C     ==--------------------------------------------------------------==
C     GTCOOR sets up the initial coordinates and velocities
C     either by reading them from file using L<RDCFG>
C     or generating them, depending on L<NTX> and L<TEMPI>.
C
C     GTCOOR reads in or generates initial coordinates L<X>, initial
C     velocities L<V> (depending on L<NTX> and L<TEMPI>), 
C     the simulation box L<BOX>, position re(con)straint
C     reference positions L<XC>, and L<RIIAVE> time averaged
C     distances for distance restraining.
C
C     Depending on the input switches defined, L<GTCOOR> calls
C     L<RDCFG> to read some of the quantities from file, and generates
C     others. Subroutine L<GAUSS> is called in order to produce
C     the initial velocities.
C
C     If the simulation is performed in 4D, L<GT4NDX> is called
C     in order read in the L<C4D> array.
C
C     NR     the total number of atoms in the system
C     Other arrays effected by this routine:
C     X,V,F,SX,XC,XR
C     BOX, C4D, RIIAVE, CXC, JRC
C
C     Subroutine L<RDXRST> can be called which reads in L<CXC>
C     (the atomic B factors) and/or the position re(con)straint
C     reference positions L<XC>.
C     The logic for reading in position restraining 
C     reference  coords and restraining pot. coefficients
C     is the following:
COMMVERB
C     IF (NTR  .NE. NTROFF) THEN
C     position restraining is turned on
C        IF (NTX < 8) THEN
C           read pos. rest. coords using L<RDCFG>
C           IF (NTR .EQ. NTRCHO) THEN
C              set constants according to CH0
C           ELSEIF (NTR .EQ. NTRBFA) THEN
C              read atomic B factors from fort.24     CASE 1
C           ENDIF
C        ELSE
C           IF (NTR .EQ. NTRCHO) THEN
C              read just coords using L<GTXRST>       CASE 2
C           ELSEIF (NTR .EQ. NTRBFA) THEN
C              read coords and atomic B
C              factors from fort.24                   CASE 3
C           ENDIF
C        ENDIF
C     ENDIF
C
C     fort.24 can thus contain in
COMMLIST
CITM  CASE 1:
C     just abf
CITM  CASE 2:
C     just coords
CITM  CASE 3:
C     coords and abf
C
C     and must contain the correct block types accordingly.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
      INCLUDE 'toposz.h'
      INCLUDE 'units.h'
      INCLUDE 'topoar.h'
      INCLUDE 'md.h'
      INCLUDE 'box.h'
C args
      INTEGER NR,NDSZE
      CHARACTER PRGSTR*(*)
      real*8 DLSUM(NDSZE),DMSUM(NDSZE)
C local vars
      INTEGER I,II,J,M,NN,IUNIT,NR3,NRPT3,NTIT
      LOGICAL LFORM,LRDVEL,LRDXC,LRDDIS,LRDSTC,LRDPER,LDOEM
      LOGICAL LRD4DX,LRD4DV,LRDBOX
      LOGICAL LRDJV,LRDLE,LSTOP
      real*8 BOLFAC,AHEAT,SD
      CHARACTER TITLE*(MXCOTI)
      DIMENSION TITLE(MAXLNS)
C begin
      LDOEM = (IDOPRO .EQ. IDOEM)

C set the whole arrays to 0.0 for better debugging
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
      DO I=1,MAXXCO
         X(I) = 0.0d+00
         V(I) = 0.0d+00
         F(I) = 0.0d+00
         XC(I)= 0.0d+00
         XR(I)= 0.0d+00
         SX(I)= 0.0d+00
      ENDDO
      CALL mm_AZZERO(RIIAVE,MAXTDR)

      NR3 = NDIM*NR
      NRPT3= NDIM*(NRP+NPM)

!  QMMM: already opened
!     CALL DETOPN('XVI',IUNIT,LFORM)
!     IF (IUNIT .LT. 0) THEN
!        PRINT *,PRGSTR,': failed to open file for ',
!    $      'reading initial coordinates'
!        STOP
!     ENDIF
      IUNIT=IOXVI
      lform=.true.

C for setting these logicals, see 'md.h'
      IF (LDOEM) THEN
         LRDVEL = .FALSE.
      ELSE
         LRDVEL = (TEMPI .LE. EPS .AND.
     $        (NTX .EQ. NTXV .OR. NTX .EQ. NTSX))
      ENDIF

      LRDXC  = (NTR .NE. NTROFF .AND. NRDRX .EQ. NRDRXY)

C do I read 4th dimensional components ?
      LRD4DX =(NT4DIM .NE.  NT4OFF .AND.
     $     (NT4XI .EQ. NT4XRM .OR. NT4XI .EQ. NT4XRR))
      IF (LDOEM) THEN
         LRD4DV = .FALSE.
      ELSE
         LRD4DV =(LRDVEL .AND. NT4DIM .NE.  NT4OFF .AND.
     $        (NT4XI .EQ. NT4XRR))
      ENDIF

C read stochastic integrals ?
      LRDSTC = (IDOPRO .EQ. IDOSD .AND. NTX .EQ. NTSX)

C read box
      LRDBOX  = (NTB .NE. NTBVAC .AND. NRDBOX .EQ. NRDBXY)

C read perturbation data ?
      LRDPER = (NTG .NE. NTGOFF .AND. NRDGL .EQ. NRDGLY)

C read distance restraint averages ?
      LRDDIS = (NTDR .LT. NTDROF .AND. LRDDR)

C read j-val restraint data ?
      LRDJV = (NTJR .LT. NTJROF .AND. NRDJR .EQ. NRDJRY)

C local elevation memory ?
      LRDLE = (NTLE .NE. NTLEOF .AND. NRDLE .EQ. NRDLEY)

      CALL RDCFG(IUNIT,LFORM,
     $     MAXLNS,NTIT,TITLE,
     $     LRD4DX,LRD4DV,
     $     NR,NDIM,NDR*NPM,NDJV*NPM,
     $     .TRUE.,X,
     $     LRDVEL,V,
     $     LRDBOX,BOX,
     $     LRDXC,XC,
     $     LRDSTC,SX,IG,
     $     LRDPER,NDSZE,RLAM,DLSUM,RMU,DMSUM,
     $     LRDDIS,RIIAVE,
     $     LRDJV,COSQAV,COSIAV,
     $     LRDLE,NDLE,NLECFG,NLEMEM,ILEMEM,NLEVST)

C read any other blocks from initial configuration here
      CALL CLSFIL(IUNIT)

C We have to check the box (CHKBOX) and cutoff distance (CHKRCT) here 
C and set the other BOX related vars if we have read the box from file.
C Otherwise, this is done in RDMD.
      IF (LRDBOX) THEN
         CALL CHKBOX(PRGSTR,.FALSE.)
         CALL CHKRCT(PRGSTR,.FALSE.)
      ENDIF

C calculate "derived energies" of DLSUM and DMSUM if needed
C as they are not saved to file but can be recalculated
C in CLCDRE in file force.f
C
      IF (LRDPER) THEN
         CALL CLCDRE(DLSUM)
         CALL CLCDRE(DMSUM)
      ENDIF

C     print something to screen
      CALL PRTIT('INITIAL CONFIGURATION TITLE',NTIT,TITLE)

C take velocities from a Maxwellian if required
      IF (.NOT. LDOEM) THEN
         IF (TEMPI .GT. EPS) THEN
            PRINT 560,3
            BOLFAC = BOLTZ*TEMPI
            I = 1
            DO NN=1,NPM
               DO J=1,NRP
                  SD=DSQRT(BOLFAC*WINV(J))
                  CALL GAUSS(0.d0,SD,NDRMAX,V(I),IG)
                  I = I+NDIM
               ENDDO
            ENDDO
            DO NN=1,NSM
               DO J=1,NRAM
                  SD=DSQRT(BOLFAC*WINVS(J))
                  CALL GAUSS(0.d0,SD,NDRMAX,V(I),IG)
                  I = I+NDIM
               ENDDO
            ENDDO
         ELSE
            IF (DABS(HEAT).GT.EPS) THEN
               IF (HEAT .LE. 0.d0) THEN
                  WRITE(6,500)
                  AHEAT=DABS(HEAT)
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
                  DO I=1,NRPT3
                     V(I)=V(I)*AHEAT
                  ENDDO
               ELSE
                  WRITE(6,501)
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
                  DO I=1,NR3
                     V(I)=V(I)*HEAT
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
      ENDIF

C initialize 4D position and velocities if necessary
      IF (NT4DIM .NE. NT4OFF) THEN
C IF we perform
C a) the nonbonded interactions in 4D AND
C b) a Mu perturbation THEN
C --> all 4D atoms must be in the perturbation list
C     
         IF (((NTG .EQ. NTGMU) .OR. (NTG .EQ. NTGBOT)) .AND.
     $        NDO4D(N4DNBD) .EQ. MAXDIM) THEN
            I = 0
            LSTOP = .FALSE.
            DO NN=1,NPM
               DO J=1,NRP
                  I = I + 1
                  IF (C4D(I).GE.0.d0 .AND. IPERT(J).EQ.NOPERT) THEN
                     LSTOP = .TRUE.
                     WRITE(6,570)PRGSTR,NTG
                     WRITE(6,580)I,J,NN
                  ENDIF
               ENDDO
            ENDDO
C now check solvent: a solvent atom can t be in the perturbation
C     list, so we can t have any solvent atom in 4D in this case.
            DO NN=1,NSM
               DO J=1,NRAM
                  I = I + 1
                  IF (C4D(I).GE.0.d0) THEN
                     LSTOP = .TRUE.
                     WRITE(6,571)PRGSTR,NTG
                     WRITE(6,580)I,J,NN
                  ENDIF
               ENDDO
            ENDDO

            IF (LSTOP) THEN
               WRITE(6,590)
               CALL STOPGM('MDUTILS','FATAL ERROR')
            ENDIF
         ENDIF

C first check for velocities
         IF (.NOT. LDOEM .AND.
     .        (NT4XI .EQ. NT4X0V .OR. NT4XI .EQ. NT4XBM .OR.
     .        NT4XI .EQ. NT4XRM)) THEN
            PRINT 560,4
            BOLFAC = BOLTZ*TEMP4I
            I = 0
            II= 1
            DO NN=1,NPM
               DO J=1,NRP
                  IF (C4D(II).GE.0.d0) THEN
                     SD=DSQRT(BOLFAC*WINV(J))
                     CALL GAUSS(0.d0,SD,1,V(I+NDHMIN),IG)
                  ELSE
!$OMP parallel do private(M)
                     DO M=NDHMIN,NDHMAX
                        V(I+M) = 0.d0
                     ENDDO
                  ENDIF
                  II = II+1
                  I = I+NDIM
               ENDDO
            ENDDO
C
            DO NN=1,NSM
               DO J=1,NRAM
                  IF (C4D(II).GE.0.d0) THEN
                     SD=DSQRT(BOLFAC*WINVS(J))
                     CALL GAUSS(0.d0,SD,1,V(I+NDHMIN),IG)
                  ELSE
!$OMP parallel do private(M)
                     DO M=NDHMIN,NDHMAX
                        V(I+M) = 0.d0
                     ENDDO
                  ENDIF
                  II = II+1
                  I = I+NDIM
               ENDDO
            ENDDO
         ENDIF
      
C     now positions
         IF (NT4XI .EQ. NT4X0V) THEN
C     set to zero
            I = 0
            DO J=1,NR
!$OMP parallel do private(M)
               DO M=NDHMIN,NDHMAX
                  X(I+M) = 0.d0
               ENDDO
               I = I+NDIM
            ENDDO
         ELSEIF (NT4XI .EQ. NT4XBM) THEN
            WRITE(6,*)'BOLTZMANN DIST OF 4D POS NOT IMPLEMENTED!'
C     boltzmann dist of positions
            I = 0
            DO J=1,NR
!$OMP parallel do private(M)
               DO M=NDHMIN,NDHMAX
                  X(I+M) = 0.d0
               ENDDO
               I = I+NDIM
            ENDDO
         ENDIF
      ENDIF

C******READ OR GENERATE FRICTION COEFFS IF NEEDED,
C     AND INITIALIZE SWINK and SWINKS arrays.
C     DO THIS --AFTER-- RDPERT, AS RDPERT MAY HAVE
C     MODIFIED THE MASS ARRAYS
      IF (IDOPRO .EQ. IDOSD) THEN
         CALL GTFRIC(NR,NDIM,GAM,X,
     .        CC1,CC2,CC3,CC4,CC5,CC6,CC7,CC8,CC9)
         CALL CLCSWI(NRP,NRAM,WINV,WINVS,BOLTZ*TEMPSD,
     .        SWINK,SWINKS)
      ENDIF

C position re(con)straining reference positions
      IF (NTR .NE. NTROFF) THEN
         WRITE(6,2110)
         IF (LRDXC) THEN
            WRITE(6,505)
         ELSE
            WRITE(6,506)
         ENDIF
C read pos and/or b factors from file iff necessary
         CALL RDXRST(NR,NDIM,.NOT. LRDXC,NTR .EQ. NTRBFA)

         IF (NTR .EQ. NTRCHO) THEN
C     set b factors using CHO
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
            DO I=1,NR
               CXC(I) = CHO
            ENDDO
         ELSEIF (NTR .EQ. NTRBFA) THEN
!$OMP parallel do private(I)
            DO I=1,NR
               IF (CXC(I).GT.0.d0) THEN
                  CXC(I) = CHO/CXC(I)
               ENDIF
            ENDDO
         ENDIF
      ENDIF

C check restriction: a 4D atom may not be position constrained
      IF (NTR .EQ. NTRCON .AND. NT4DIM .NE.  NT4OFF) THEN
         LSTOP = .FALSE.
         DO I=1,NRCON
            IF (C4D(JRC(I)).GE.0.d0) THEN
               WRITE(6,520)PRGSTR,JRC(I)
               LSTOP = .TRUE.
            ENDIF
         ENDDO
         IF (LSTOP) THEN
            CALL STOPGM('MDUTILS','FATAL ERROR')
         ENDIF
      ENDIF

C j-value constraints
      IF (NTJR .NE. NTJROF) THEN
         IF (ABS(NTJR) .EQ. NTJCJR) THEN
C set CPJV to CJR
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
            DO I=1,NDJV
               CPJV(I) = CJR
            ENDDO
         ELSEIF (ABS(NTJR) .EQ. NTJRED) THEN
C     set CPJV to CJR*CPJV
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
            DO I=1,NDJV
               CPJV(I) = CPJV(I)*CJR
            ENDDO
         ELSE
            WRITE(6,*)'ILLEGAL VALUE FOR NTJR = ',NTJR
            CALL STOPGM('MDUTILS','FATAL ERROR')
         ENDIF

         CALL SETJAV(NTJR .LT. NTJROF,NRDJR .EQ. NRDJRY,
     .        COSQAV,COSIAV,NPM,NRP,NRP*NPM,NDIM,NDJV,
     .        AJV,BJV,CJV,PJR0,PSJR,IPJV,JPJV,KPJV,LPJV,X)

      ENDIF

C format statements
 500  FORMAT (1X,'VELOCITIES OF SOLUTE ARE RESCALED')
 501  FORMAT (1X,'VELOCITIES OF ALL ATOMS ARE RESCALED')
 505  FORMAT (1X,'POSITION RESTRAINING POSITIONS ARE READ ',
     $        'FROM INITIAL CONFIG. FILE')
 506  FORMAT (1X,'POSITION RESTRAINING POSITIONS ARE READ ',
     $        'FROM A SEPARATE FILE')
 560  FORMAT (1X,I1,'D VELOCITIES ARE TAKEN FROM A MAXWELLIAN',
     $     ' DISTRIBUTION',/)

 520  FORMAT (1X,A6,': IMPLEMENTATION RESTRICTION: 4D ATOM MAY',
     $     'NOT BE POSITIONALLY CONSTRAINED:',I5)
 570  FORMAT (1X,A6,': INCONSISTENT INPUT: NON-PERTURBED 4D ATOM',
     $     ' WITH NTG =',I5)
 571  FORMAT (1X,A6,': INCONSISTENT INPUT: 4D SOLVENT ATOM WITH',
     $     ' NTG =',I5)
 580  FORMAT (1X,'ATOM NUMBER :',I5,'(ATOM ',I5,' IN MOLECULE ',I5,')')
 590  FORMAT (1X,'IF THE NONBONDED INTERACTION IS IN 4D, THEN',
     $     ' EVERY 4D ATOM MUST BE IN THE PERTURBATION LIST')

 2110 FORMAT (//,1X,'R E F E R E N C E   C O O R D I N A T E S',
     $     '   RE ( C O N ) S T R A I N I N G',/)
C end GTCOOR
C     ==--------------------------------------------------------------==
      RETURN
      END 
C     ==================================================================
      SUBROUTINE SETJAV(LTAVE,LRED,COSQAV,COSIAV,NPM,NRP,NRPT,NDIM,
     $     NDJV,AJV,BJV,CJV,PJR0,PSJR,IPJV,JPJV,KPJV,LPJV,
     $     X)
C     ==--------------------------------------------------------------==
C     SETJAV sets the initial values for L<COSQAV>, L<COSIAV> which
C     are used in subroutine L<RESTJ> when calculating time
C     averaged J-value restraints.
C     The initial values for the averages are set to the optimal
C     value corresponding to PJR0. See subroutine L<RESTJ> for
C     details.
C
C     <COSQAV>, L<COSIAV> are set to zero if time averaging is not used.
C
C     SETJAV is called by L<GTCOOR> if L<NTJR> .NE. L<NTJROF>.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'cobond.h'
C args
      LOGICAL LTAVE,LRED
      INTEGER NPM,NRP,NDJV,NRPT,NDIM
      INTEGER IPJV(NDJV),JPJV(NDJV),KPJV(NDJV),LPJV(NDJV)
      real*8 COSQAV(MAXJ),COSIAV(MAXJ)
      real*8 AJV(NDJV),BJV(NDJV),CJV(NDJV),PJR0(NDJV),PSJR(NDJV)
      real*8 X(NDIM*NRPT)
C local parameters
C     we have to be able to store coordinates and force of
C     maximum 4 atoms
      INTEGER NTMPAT
      PARAMETER (NTMPAT = 4)
      INTEGER MAXLCX
      PARAMETER (MAXLCX = NTMPAT*MAXDIM)
      INTEGER NCALCD
      PARAMETER (NCALCD = 3)
      LOGICAL LEVERY
      PARAMETER (LEVERY = .FALSE.)
cmb      real*8 ONE
cmb      PARAMETER (ONE = 1.0d+0)
C local vars
      INTEGER I
      real*8 COSVAL

#include "inc_mdutils1.h"

C      INTEGER IP,JP,KP,LP,IDUM,NN
C      real*8 RMOLIN,XDUM,PHI,TETA
C      real*8 PYE, CONV, PYETWO
C      real*8 DPHIDR(MAXLCX)
C begin
      IF (LTAVE) THEN
         IF (.NOT. LRED) THEN
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
            DO I=1,MAXJ
               COSQAV(I) = 0.0d+00
               COSIAV(I) = 0.0d+00
            ENDDO

C     if we have time averaging AND didn t read in the initial values
C     for COSQAV,COSIAV from file in RDCFG, then initialise them here:
C
C     assign to COSIAV the optimal value corresponding to PJR0.
C     In this case, it is wrong to divide by NPM: remember to
C     comment out the loop 200 when using this part of the code.
C     This code is as described in the errata to the GROMOS96 manual.
C
C     Solve quadratic equation: A*x*x + B*x + C - pjr0 = 0 for x, where
C     x = cos(teta)
!$OMP parallel do private(I,COSVAL)
            DO I=1,NDJV
               COSVAL=(-BJV(I)+DSQRT(BJV(I)**2-4.0d0*AJV(I)*
     .                 (CJV(I)-PJR0(I))))/(2.0d0*AJV(I))
               IF(DABS(COSVAL) .GT. 1.0d+0) THEN
                 COSVAL=(-BJV(I)-DSQRT(BJV(I)**2-4.0d0*AJV(I)*
     .                   (CJV(I)-PJR0(I))))/(2.0d0*AJV(I))
               ENDIF
               COSIAV(I) = COSIAV(I) + COSVAL
               COSQAV(I) = COSQAV(I) + COSVAL*COSVAL
            ENDDO
C-----------------------------------------------------------------
C     assign to COSIAV the instantaneous value of the coordinates.
C     In the case of multiple molecules, we get the average instantaneous
C     value over all molecules.
C     This is as described in the GROMOS96 manual, but corrected
C     in an errata. Consequently, this part of the code is NOT used.
C            PYE = 4.0 * ATAN(ONE)
C            CONV = 180.0 / PYE
C            PYETWO = 2.0 * PYE
C
C            DO 110 NN=1,NPM
C               DO 120 I=1,NDJV
C                  IP = IPJV(I)
C                  JP = JPJV(I)
C                  KP = KPJV(I)
C                  LP = LPJV(I)
C     calculate the dihedral angle PHI (in degrees)
C     and derivative DPHIDR using DIHANG
C                  CALL DIHANG(NCALCD,NN,NRP,NFDER,IDINCS,LEVERY,
C     $                 1,IP,JP,KP,LP,IDUM,
C     $                 0,XDUM,XDUM,IDUM,
C     $                 NRPT,X,DPHIDR,PHI,XDUM,XDUM,
C     $                 .FALSE.,XDUM,0,XDUM,XDUM,XDUM,XDUM,
C     $                 IDUM,IDUM,XDUM,
C     $                 .FALSE.,XDUM,XDUM,XDUM)
C                  PHI = PHI / CONV
C                  TETA = PHI + PSJR(I)
C                  IF (TETA .GT.  PYE) TETA = TETA - PYETWO
C                  IF (TETA .LT. -PYE) TETA = TETA + PYETWO
C                  COSVAL = COS(TETA)
C                  COSIAV(I) = COSIAV(I) + COSVAL
C                  COSQAV(I) = COSQAV(I) + COSVAL*COSVAL
C 120           CONTINUE
C 110        CONTINUE
C            RMOLIN = 1.0/REAL(NPM)
C            DO 200 I=1,NDJV
C               COSIAV(I) = COSIAV(I)*RMOLIN
C               COSQAV(I) = COSQAV(I)*RMOLIN
C 200        CONTINUE
C-----------------------------------------------------------------
         ENDIF
      ELSE
C no time averaging, set to zero
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
         DO I=1,MAXJ
            COSQAV(I) = 0.d0
            COSIAV(I) = 0.d0
         ENDDO
      ENDIF
C setjav
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GTFRIC(NATTOT,NDIM,GAM,X,
     $     CC1,CC2,CC3,CC4,CC5,CC6,CC7,CC8,CC9)
C     ==--------------------------------------------------------------==
C     Initializes GAM and the CC? auxiliary arrays for SD.
C     Other arguments not modified.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'md.h'
C  args
       INTEGER NATTOT,NDIM
       real*8 GAM(NATTOT),X(NDIM*NATTOT)
       real*8 CC1(NATTOT),CC2(NATTOT),CC3(NATTOT),CC4(NATTOT)
       real*8 CC5(NATTOT),CC6(NATTOT),CC7(NATTOT)
       real*8 CC8(NATTOT),CC9(NATTOT)
C local vars
       INTEGER I
       LOGICAL LCLCGA
C begin
       LCLCGA = .FALSE.
       IF (NTFR .EQ. NTFR0) THEN
          CALL mm_AZZERO(GAM,NATTOT)
       ELSEIF (NTFR .EQ. NTFRFR) THEN
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
          DO I=1,NATTOT
             GAM(I) = CFRIC
          ENDDO
       ELSEIF (NTFR .EQ. NTFRGA) THEN
C read from file, then multiply by CFRIC
          CALL RDGAM(NATTOT,GAM)
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
          DO I=1,NATTOT
             GAM(I) = GAM(I)*CFRIC
          ENDDO
       ELSEIF (NTFR .EQ. NTFRCA ) THEN
C calc using FRIC
          LCLCGA = .TRUE.
       ELSE
          WRITE(6,*)'GETFRIC: invalid value for NTFR =',NTFR
          CALL STOPGM('MDUTILS','FATAL ERROR')
       ENDIF

C     now initialize all the CC variables, and possibly
C     GAM if not set above
!      CALL FRIC(NATTOT,NDIM,LCLCGA,NBREF,RCUTF,CFRIC,DT,
!    $      X,GAM,
!    $      CC1,CC2,CC3,CC4,CC5,CC6,CC7,CC8,CC9)
C end gtfric
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE CLCSWI(NRP,NRAM,WINV,WINVS,BOLTZT,
     $     SWINK,SWINKS)
C     ==--------------------------------------------------------------==
C     CLCSWI calculates the SWINK arrays used in L<RUNMD>
C     when performing SD.
C     SWINK = SQRT(KT/MASS) for solute
C     SWINKS= SQRT(KT/MASS) for solvent
C     All other variables are not modified.
C     ==--------------------------------------------------------------==
C args
      INTEGER NRP,NRAM
      real*8 WINV(NRP),WINVS(NRAM),BOLTZT
      real*8 SWINK(NRP),SWINKS(NRAM)
C local vars
      INTEGER J
C begin
!$OMP parallel do private(J)
#ifdef __SR11000
*poption parallel, tlocal(J)
#endif
      DO J=1,NRP
         SWINK(J) = DSQRT(BOLTZT*WINV(J))
      ENDDO
!$OMP parallel do private(J)
#ifdef __SR11000
*poption parallel, tlocal(J)
#endif
      DO J=1,NRAM
         SWINKS(J) = DSQRT(BOLTZT*WINVS(J))
      ENDDO
C end clcswi
      RETURN
      END
C     ==================================================================
      SUBROUTINE RDIHRS(PRGSTR,NRPT)
C     ==--------------------------------------------------------------==
C     RDIHRS reads in and converts dihedral angle restraining
C     specifications.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
      INCLUDE 'md.h'
C args
      CHARACTER PRGSTR*(*)
      INTEGER NRPT
C local vars
      INTEGER K
      REAL*8 DEG2RD,DGRDSQ
C begin
      DEG2RD = DATAN(1.0d+00)/4.5d+01
      DGRDSQ = DEG2RD*DEG2RD

      WRITE(6,1510)
      CALL GTDHRS(NRPT)
      IF (NTDLR .EQ. NTDLRC) THEN
C     overwrite CPLR with CDLR
!$OMP parallel do private(K)
#ifdef __SR11000
*poption parallel, tlocal(K)
#endif
         DO K=1,NDLR
            CPLR(K) = CDLR
         ENDDO
      ELSEIF (NTDLR .EQ. NTDLRF) THEN
C     scale CPLR with CDLR
!$OMP parallel do private(K)
#ifdef __SR11000
*poption parallel, tlocal(K)
#endif
         DO K=1,NDLR
            CPLR(K) = CPLR(K)*CDLR
         ENDDO
      ELSE
         WRITE(6,*)PRGSTR,'invalid value for NTDLR = ',NTDLR
         CALL STOPGM('MDUTILS','FATAL ERROR')
      ENDIF
C     Now convert the force constant from energy per degrees **2 to energy
C     per radians**2 and of the minimum energy angle from degress to
C     radians.

!$OMP parallel do private(K)
#ifdef __SR11000
*poption parallel, tlocal(K)
#endif
      DO K=1,NDLR
         CPLR(K) = CPLR(K)/DGRDSQ
         PDLR(K) = PDLR(K)*DEG2RD
      ENDDO

C formats
 1510 FORMAT (//,1X,'R E S T R A I N E D   D I H E D R A L S')
C end rdihrs
      RETURN
      END
C     ==================================================================
      SUBROUTINE CHKPRT(PRGSTR)
C     ==--------------------------------------------------------------==
C     Subroutine CHKPRT:
C     Make sure every atom pair IEB, JEB in the perturbed atom pair list
C     is excluded in the topology in order to avoid "double counting"
C     of the interaction in the non bonded routines
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'pertsz.h'
      INCLUDE 'pertar.h'
C args
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER NN,IMIN,JMAX,ITMP,IFRST,ILAST
      LOGICAL LERR,LFOUND
C begin
      LERR = .FALSE.
      DO 10 NN=1,NEB
         IMIN = IEB(NN)
         JMAX = JEB(NN)
         IF (IMIN .GT. JMAX) THEN
            ITMP = IMIN
            IMIN = JMAX
            JMAX = ITMP
         ENDIF
C look for JMAX in exclusion list of IMIN
         IFRST = KNE(IMIN)+1
         ILAST = KNE(IMIN)+INE(IMIN)
         LFOUND = .FALSE.
 20      IF (IFRST .LE. ILAST .AND.  .NOT. LFOUND) THEN
            LFOUND = (JSNE(IFRST) .EQ. JMAX)
            IFRST = IFRST + 1
            GOTO 20
         ENDIF
         IF (.NOT. LFOUND) THEN
            WRITE(6,500)PRGSTR
            WRITE(6,510)IEB(NN)
            WRITE(6,520)JEB(NN)
            LERR = .TRUE.
         ENDIF
 10   CONTINUE

      IF (LERR) THEN
         CALL STOPGM('MDUTILS','FATAL ERROR')
      ENDIF

 500  FORMAT (1X,A,': ERROR: Perturbed atom pair must',
     .     ' be excluded atoms in the topology !')
 510  FORMAT (1X,'Atom IEB:',I5)
 520  FORMAT (1X,'Atom JEB:',I5,/)
C end chkprt
      RETURN
      END
C     ==================================================================
      SUBROUTINE CHKCG(PRGSTR)
C     ==--------------------------------------------------------------==
C     Subroutine CHKCG:
C     It is checked whether no charge group has more than
C     MXATCG atoms.
C     If the topology read in has more than MXATCG atoms, then
C     an error is printed and STOP is invoked.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
C args
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER IFTOPA,ILTOPA,INUMAT,ITOPCG
      LOGICAL LERR
C begin
      LERR = .FALSE.
      IFTOPA = 1
      DO ITOPCG=1,NCAG
         ILTOPA = INC(ITOPCG)
         INUMAT = ILTOPA - IFTOPA + 1
         IF (INUMAT .GT. MXATCG) THEN
            WRITE(6,500)PRGSTR,ITOPCG
            WRITE(6,510)INUMAT
            WRITE(6,520)MXATCG
            LERR = .TRUE.
         ENDIF
         IFTOPA = ILTOPA + 1
      ENDDO

      IF (LERR) THEN
         CALL STOPGM('MDUTILS','FATAL ERROR')
      ENDIF
 500  FORMAT (1X,A,': too many atoms in charge group number',I4)
 510  FORMAT (1X,'number of atoms       :',I5)
 520  FORMAT (1X,'maximum allowed MXATCG=',I5)
C chkcg
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE CHKPCN(PRGSTR,NRCON,JRC)
C     ==--------------------------------------------------------------==
C     Subroutine CHKPCN:
C     Check that we don t perturb the mass of a
C     positionally constrained atom.
C     ==--------------------------------------------------------------==
C args
      CHARACTER PRGSTR*(*)
      INTEGER NRCON,JRC(NRCON)
C includes
      INCLUDE 'pertsz.h'
      INCLUDE 'pertar.h'
C local params
      real*8 EPS
      PARAMETER (EPS=1.0D-05)
C local vars
      INTEGER N,II
      LOGICAL LMAPT,LERR
C begin
      LERR = .FALSE.
      DO 10 N=1,NJLA
         LMAPT = (DABS(WMA(N) - WMB(N)).GT.EPS)
         IF (LMAPT) THEN
C look for N in constrained list
            II = 1
 20         IF (II .LE. NRCON) THEN
               IF (JRC(II) .EQ. N) THEN
                  WRITE(6,500)PRGSTR
                  WRITE(6,510)N
                  LERR = .TRUE.
               ENDIF
               II = II+1
               GOTO 20
            ENDIF
         ENDIF
 10   CONTINUE

      IF (LERR) THEN
         CALL STOPGM('MDUTILS','FATAL ERROR')
      ENDIF
 500  FORMAT (1X,A,': The mass of a positionally constrained atom may',
     .     'NOT be perturbed!')
 510  FORMAT (1X,'Conflict on atom number:',I5)
C end chkpcon
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE CHKNRE(PRGSTR,NATTOT)
C     ==--------------------------------------------------------------==
C     Subroutine CHKNRE:
C     CHKNRE checks the value for L<NEGR> and L<NRE> input and fills in
C     the L<IAGRP> with the L<NRELKP> ,L<NUSNRE> and L<NUNRE2> variables.
C     It also creates formatting strings L<EPTFMT> used in SUBR. 
C     L<PRNRG> for writing out the energy matrices.
C     ==--------------------------------------------------------------==
C args
      CHARACTER PRGSTR*(*)
      INTEGER NATTOT
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
      INCLUDE 'md.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'forcear.h'
C local vars
      INTEGER I,J,IND,IVAL,IFIRST,ILAST
      INTEGER IRES,ICLNUM,ISPACE
      CHARACTER CHARO
C begin
      IF (NEGR .EQ. 0) THEN
C the user has specified no energy groups;
C init IAGRP and NRELKP to use one group which
C is not printed out.
         CHARO = CHAR(1)
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
         DO I=1,MAXNAT
            IAGRP(I) = CHARO
         ENDDO
         NUSNRE = 0
         NUNRE2 = 1

!$OMP parallel do private(I,J) __COLLAPSE2
#ifdef __SR11000
*poption parallel, tlocal(I,J)
#endif
         DO I=1,MAXNRE
            DO J=I,MAXNRE
               NRELKP(I,J) = 1
               NRELKP(J,I) = 1
            ENDDO
         ENDDO

      ELSE
         CHARO = CHAR(0)
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
         DO I=1,MAXNAT
            IAGRP(I) = CHARO
         ENDDO
         NUSNRE = 0
         NUNRE2 = 0
C check input
         DO 5 I=1,NEGR-1
            IF (NRE(I) .LT. 0 .OR. NRE(I) .GT. NRE(I+1)) THEN
               WRITE(6,500)PRGSTR
               CALL STOPGM('MDUTILS','FATAL ERROR')
            ENDIF
 5       CONTINUE

         IF (NRE(NEGR) .NE. NATTOT) THEN
            WRITE(6,500)PRGSTR
            WRITE(6,'('' NEGR ='',I6,4x,''NRE(NEGR) ='',I10)') 
     .            NEGR,NRE(NEGR) ! cmb
            PRINT *,'The last entry must be = NATTOT =',NATTOT
            CALL STOPGM('MDUTILS','FATAL ERROR')
         ENDIF

         IFIRST = 1
         DO I=1,NEGR
            ILAST = NRE(I)
            IF (ILAST .GE. IFIRST) THEN
               NUSNRE = NUSNRE+1
!$OMP parallel do private(IVAL)
               DO IVAL=IFIRST,ILAST
                  IAGRP(IVAL) = CHAR(NUSNRE)
               ENDDO
               
C generate strings for heading of matrix
               CALL SCOIN2(IFIRST,ILAST,6,EPTSTR(NUSNRE))
               IFIRST = ILAST + 1
            ENDIF
         ENDDO

         NUNRE2 = NUSNRE*(NUSNRE+1)/2

C     generate lookup table
         IND = 1
         DO I=1,NUSNRE
            DO J=I,NUSNRE
               NRELKP(I,J) = IND
               NRELKP(J,I) = IND
               IND = IND + 1
            ENDDO
         ENDDO
         
C generate format string table for writing out
C the matrix with headings
         ICLNUM = NUSNRE
         ISPACE = 0
         DO I = 1, NUSNRE
C zero times a space gives me one space !!
C have to do it this way
            IF (ISPACE .EQ. 0)THEN
               WRITE (EPTFMT(I),IOSTAT = IRES,FMT = 601)ICLNUM
            ELSE
               WRITE (EPTFMT(I),IOSTAT = IRES,FMT = 600)
     .              ISPACE,ICLNUM
            ENDIF
            IF (IRES .NE. 0)THEN
               WRITE(6,*)'FMTSTR FORMAT ERROR !!!',ISPACE
               CALL STOPGM('MDUTILS','FATAL ERROR')
            ENDIF

            ICLNUM = ICLNUM-1
            ISPACE = ISPACE+1
         ENDDO
      ENDIF

 500  FORMAT (1X,A,': invalid input for NRE !')
 600  FORMAT ('(1X,A12,',I2,'(13X),',I2,'(1X,E12.5))')
 601  FORMAT ('(1X,A12,',I2,'(1X,E12.5))')
C end chknre
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE CHKNSP(PRGSTR,NDIM,NPM,NRP,NCAG,INC)
C     ==--------------------------------------------------------------==
C     Subroutine CHKNSP:
C
C     CHKNSP checks consistency between the topology file
C     (L<NRP> = number of solute atoms) and the L<PROMD> input:
C
C     L<NSPM> and L<NSP>(K) = last atom of the K-th submolecule
C     which are read in a SUBMOLECULE input block.
C
C     Atoms in a given charge group must be in the same submolecule.
C     The total number of solute submolecule coordinates must
C     be < L<MXSBCO>.
C
C     NOTE: his routine must not be called with NSPM=0, as is the case
C     when there was no SUBMOLECULE block on the PROMD input file.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
C args
      CHARACTER PRGSTR*(*)
      INTEGER NDIM,NPM,NRP,NCAG
      INTEGER INC(NCAG)
C local params
      INTEGER MAXFLD
      PARAMETER (MAXFLD = 5)
      INTEGER MAXSTR
      PARAMETER (MAXSTR = 2*MAXFLD+1)
C local vars
      INTEGER K,KCG
      INTEGER KFRST,KLST
C begin
C check for centre of mass velocity coordinates
      IF (NDIM*NSPM*NPM .GT. MXSBCO) THEN
         WRITE(6,600)PRGSTR
         WRITE(6,610)NSPM
         WRITE(6,620)NPM
         WRITE(6,630)NDIM
         WRITE(6,640)NDIM*NSPM*NPM
         WRITE(6,650)MXSBCO
         CALL STOPGM('MDUTILS','FATAL ERROR')
      ENDIF

      KFRST = 1
      DO 5 K=1,NSPM
         KLST = NSP(K)

         IF (KLST .GT. NRP) THEN
            WRITE(6,490)PRGSTR
            WRITE(6,500)NRP
            WRITE(6,510)K,KLST
            WRITE(6,*)'NSP(K) must be <= NRP'
            CALL STOPGM('MDUTILS','FATAL ERROR')
         ENDIF

C KLST must be in the INC array
         KCG=1
C while loop
 7       IF (KCG .LT. NCAG) THEN
            IF (KLST .GT. INC(KCG)) THEN
               KCG = KCG + 1
               GOTO 7
            ENDIF
C end of while loop
         ENDIF

         IF (KLST .NE. INC(KCG)) THEN  
            WRITE(6,490)PRGSTR
            WRITE(6,*)'Charge groups may not span submolecules'
            WRITE(6,480)KLST
            WRITE(6,485)K
            CALL STOPGM('MDUTILS','FATAL ERROR')
         ENDIF

         KFRST = KLST + 1
 5    CONTINUE

      IF (NSP(NSPM) .NE. NRP) THEN
         WRITE(6,*)PRGSTR,': Omitted atoms in submolecule definition!'
         WRITE(6,500)NRP
         WRITE(6,510)NSPM,NSP(NSPM)
         WRITE(6,*)'NSP(NSPM) must be equal to NRP'
         CALL STOPGM('MDUTILS','FATAL ERROR')
      ENDIF

 480  FORMAT (1X,'Atom number',I5,' defines the boundary of')
 485  FORMAT (1X,'submolecule',I5,' but not a charge group boundary')
 490  FORMAT (1X,A,': Inconsistent definition of submolecules!')
 500  FORMAT (1X,'Number of atoms per solute molecule (NRP)=',I4)
 510  FORMAT (1X,'But NSP(',I2,') = ',I4)
C
 600  FORMAT (1X,A,': Number of submolecule coordinates too large!')
 610  FORMAT (1X,'Number of submolecules per solute (NSPM):',I5)
 620  FORMAT (1X,'Number of solute molecules         (NPM):',I5)
 630  FORMAT (1X,'Dimension of simulation           (NDIM):',I5)
 640  FORMAT (1X,'Space needed               NSPM*NPM*NDIM:',I5)
 650  FORMAT (1X,'but only have                     MXSBCO=',I5)
C end chknsp
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SETCON
C     ==--------------------------------------------------------------==
C     Subroutine SETCON:
C
C     SETCON calculates the bond constraints for L<SHAKE>
C     from bond length data in the topology.
C     The data is read from L<IB>, L<JB> and L<ICB> for bonds not
C     involving hydrogen atoms and from L<IBH>, L<JBH> and L<ICBH> for
C     bonds involving hydrogen atoms.
C
C     Note: this routine checks that that the perturbation file
C     has been read in before this routine is called.
C     Every bond is checked whether it has been put into the
C     list already by L<RDPERT> before being added here in order to
C     avoid double counting of perturbed bonds.
C
C     The data is written to L<ICOG>, L<JCOG> and L<CONP>.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'pertsz.h'
      INCLUDE 'pertar.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
      INCLUDE 'md.h'
C externals
!     LOGICAL LPTOK
!     EXTERNAL LPTOK
C local vars
      INTEGER I,IAUX
      LOGICAL LADD,LERR
      CHARACTER PRGSTR*(6)
C data
      DATA PRGSTR /'SETCON'/
C begin
      LERR = .FALSE.
!     IF ((NTG .EQ. NTGLAM .OR. NTG .EQ. NTGBOT)
!    $     .AND. .NOT. LPTOK()) THEN
!        PRINT *,PRGSTR,': ERROR! perturbation file must be read before'
!        PRINT *,'constraints are set.'
!        STOP
!     ENDIF

C constraint handling of bonds
      IF (NTC .NE. NTCNON) THEN
         DO I=1,NBONH
C find I in NCBHG
            LADD = .TRUE.
            IAUX=0
C while loop
 80         IF (IAUX .LT. NBONHG .AND. LADD) THEN
               IAUX = IAUX+1
               LADD = (NCBHG(IAUX) .NE. I)
               GOTO 80
C end of while loop
            ENDIF

            IF (LADD) THEN
C not found --> add I to constraint list
               NCONG = NCONG+1
               IF (NCONG .GT. MAXCON) THEN
                  LERR = .TRUE.
               ELSE
                  CONP(NCONG) = B0(ICBH(I))**2
                  ICOG(NCONG) = IBH(I)
                  JCOG(NCONG) = JBH(I)
               ENDIF
            ENDIF
         ENDDO
      ENDIF

      IF (NTC .EQ. NTCDOB) THEN
         DO I=1,NBON
            LADD = .TRUE.
            IAUX = 0
C while loop
 100        IF (IAUX .LT. NBONG .AND. LADD) THEN
               IAUX = IAUX+1
               LADD = (NCBG(IAUX) .NE. I)
               GOTO 100
            ENDIF
C end of while loop

            IF (LADD) THEN
               NCONG = NCONG+1
               IF (NCONG .GT. MAXCON) THEN
                  LERR = .TRUE.
               ELSE
                  CONP(NCONG) = B0(ICB(I))**2
                  ICOG(NCONG) = IB(I)
                  JCOG(NCONG) = JB(I)
               ENDIF
            ENDIF
         ENDDO
      ENDIF

      IF (LERR) THEN
         WRITE(6,500)PRGSTR
         WRITE(6,501)NCONG
         WRITE(6,502)MAXCON
         CALL STOPGM('MDUTILS','FATAL ERROR')
      ENDIF

 500  FORMAT (1X,A,': too many distance constraints!')
 501  FORMAT (1X,'Need to use     : NCONG =',I5)
 502  FORMAT (1X,'maximum allowed : MAXCON=',I5)
C end SETCON
      END
C     ==================================================================
      SUBROUTINE mm_IAZZERO(IA,N)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER N
      INTEGER IA(N)
C     Variables
      INTEGER I
C     ==--------------------------------------------------------------==
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
#ifdef _vpp_
!OCL NOALIAS
#endif
      DO I=1,N
        IA(I)=0
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE mm_AZZERO(A,N)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
*     Arguments
      INTEGER N
      REAL*8 A(N)
#if defined(__VECTOR)
*     Variables
      INTEGER I
C     ==--------------------------------------------------------------==
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
#ifdef _vpp_
!OCL NOALIAS
#endif
      DO I=1,N
        A(I)=0.0D0
      ENDDO
#else
*     Variables
      INTEGER I,II
*     ==--------------------------------------------------------------==
      II=MOD(N,4)
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
      DO I=1,II
        A(I)=0.0D0
      ENDDO
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
      DO I=1+II,N,4
        A(I)  =0.0D0
        A(I+1)=0.0D0
        A(I+2)=0.0D0
        A(I+3)=0.0D0
      ENDDO
#endif
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE mm_ZAZZERO(A,N)
C     ==--------------------------------------------------------------==
#if defined(__VECTOR)
      IMPLICIT NONE
C     Arguments
      INTEGER N
      COMPLEX*16 A(N)
C     Variables
      INTEGER I
      COMPLEX*16 ZZERO
      PARAMETER  (ZZERO=(0.D0,0.D0))
C     ==--------------------------------------------------------------==
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
#ifdef _vpp_
!OCL NOVREC
#endif
      DO I=1,N
        A(I)=ZZERO
      ENDDO
#else
      IMPLICIT NONE
C     Arguments
      INTEGER N
      COMPLEX*16 A(N)
C     Variables
      INTEGER I,II
      COMPLEX*16 ZZERO
      PARAMETER  (ZZERO=(0.D0,0.D0))
C     ==--------------------------------------------------------------==
      II=MOD(N,4)
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
      DO I=1,II
        A(I)=ZZERO
      ENDDO
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
      DO I=1+II,N,4
        A(I)  =ZZERO
        A(I+1)=ZZERO
        A(I+2)=ZZERO
        A(I+3)=ZZERO
      ENDDO
#endif
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
