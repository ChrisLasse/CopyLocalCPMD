C     ==================================================================
COMMSUBR COBOND
C W.F. VAN GUNSTEREN, GRONINGEN, JUNE 1985
C The GROMOS96 functional form versions supplied by Thomas Beutler.
C Revised/documented/GROMOS87 compatibility P.Hunenberger, sept. 95
Cmb - Revised on 17 October 2006
C
C      SUBROUTINE COBOND (NCALCD,NMOL,NAM,NFCTR,IDIBND,LEVERY,
C     $     NUMBND,IB,JB,ICB,
C     $     NBTY,CB,B0,
C     $     NATTOT,X,F,EB,XB0,EB0,
C     $     LDOPER,RLAM,NPBND,CBA,CBB,BA0,BB0,EBGL,
C     $     LDO34,RMUE,C4D,EBGM)
C
C     COBOND will supply the covalent bond energies and forces
C     according to a potential energy function depending on IDIBND
C     for bonds in NCALCD dimensions (NCALCD.LE.NDIM).
C     Peridoc boundary conditions can be taken into account
C     when NTB.NE.NTBVAC (with some restrictions mentionned below
C     for the 4D case)
C     More than one identical molecule can be considered.
C
C     The cobond.h header file defines reference values for the
C     NFCTR and IDIBND switches.
C     COBOND assumes that NDIM, NTB, BOX(3), BOXH(3) and BETA are 
C     defined and transmitted through the common blocks of the box.h 
C     header file.
C     The file pertsz.h (perturbation) is optional (LDOPER.EQ.TRUE.)
C    
C     When LDOPER. is .FALSE. the standard potential,
C     V_stand(B;CB,B0) = sum over all specified bonds of
C     1.       0.25 * CB * (B**2 - B0**2)**2  if IDIBND .EQ. L<IBND4>
C              (intended for covalent bonds, GROMOS96)
C     2.       0.5 * CB * (B)**2              if IDIBND .EQ. L<IBND2>
C              (intended for path integral calculations ) 
C     3.       0.5 * (CB*2*B0**2) * (B-B0)**2 if IDIBND .EQ. L<IBNDH>
C              (limited implementation, intended for compatibility
C               with GROMOS87 program)
C     is evaluated, where B0 denotes the bond length at minimum energy 
C     and CB the force constant (according to GROMOS96 functional form).
C
C     When LDOPER is .TRUE., a LAMDA (RLAM) dependent perturbation 
C     potential, 
C     V_pert(B;RLAM;CBA,CBB,BA0,BB0) = sum over all specified bonds of
C     1.      0.25 * [ (1-RLAM)*CBA+RLAM*CBB) ] * 
C                    [ B**2 - ((1-RLAM)*BA0+RLAM*BB0)**2 ]**2
C            -0.25 * CBA *(B**2-BA0**2)**2
C            if IDIBND = L<IBND4>
C     2.      0.5  * [ (1-RLAM)*CHA+RLAM*CHB) ] *
C                    [ B - ((1-RLAM)*BA0+RLAM*BB0) ]**2
C            -0.5 *  CHA*(B-BA0)**2
C            with CHA=CBA*2*BA0**2 and CHB=CBB*2*BB0**2
C            if IDIBND = L<IBNDH>
C
C     is evaluated, where B0A/B0B denote the bond lengths at minimum 
C     energy in states A/B and CBA/CBB the force constants (according 
C     to GROMOS96 functional form) in states A/B.
C     The second term in 1. and 2. is intented to cancel the
C     V_pert(RLAM=0) contribution already taken into account by a previous
C     call to cobond (all bonds, LDOPER=.FALSE.) i.e. to avoid double 
C     counting. dV_pert/d(RLAM) (EBGL) is also evaluated.
C     Currently no perturbations are possible if IDIBND = IBND2
C     It is not implemented for NFCTR.EQ.NFNDXL.
C
C     If LDO34 is .TRUE., a MU (RMUE) dependent 3D-4D perturbation 
C     potential,
C     V_stand34(B;RMUE;CB,B0) = sum over all specified bonds of
C             + RMUE*V(B_4D;CB,B0) - RMUE*V(B_3D;CB,B0) 
C             if LDOPER is .FALSE.
C     or
C     V_pert34(B;RMUE;RLAM;CBA,CBB,BA0,BB0) = sum over all specified 
C     bonds of
C             + RMUE * V_pert(B_4D;RLAM;CBA,CBB,BA0,BB0)
C             - RMUE * V_pert(B_3D;RLAM;CBA,CBB,BA0,BB0)
C             if LDOPER is .TRUE.
C     is evaluated.
C
C     Note that the contribution of the 3D component is -RMUE and 
C     not (1-RMUE) to avoid double counting of the V_34(RMUE=0) component
C     already taken into account by previous call(s) to cobond.
C     dV_pert/d(RLAM) (EBGL) and dV_pert/d(RMUE) (EBGM) are also evaluated.
C     LDO34 requires NDIM.EQ.4 .AND. NCALCD.EQ.4.
C     It is not implemented for NFCTR.EQ.NFNDXL.OR.NFCTR.EQ.NFNDXE.
C     It is only implemented for IDIBND.EQ.IBND4.
C     Only a rectangular or truncated-octahedron box (BETA.EQ.90) can be 
C     used. In the latter case, only rectangular periodic gathering is
C     applied in the NCALCD-3 last dimensions.
C     Bond length values are returned in a mixed 3D-4D RMUE dependant form
C     too, if required:
C     RIJ = SQRT( RMUE*RIJ(4D)**2 + (1.0-RMUE)*RIJ(3D)**2 ).
C
C     Arguments:
C
C     NCALCD      the number of dimensions in which to consider
C                 all vectors, coordinate/forces (NCALCD.LE.NDIM).
C     NMOL        The meaning of this argument depends on NFCTR: it
C                 is either the number of identical molecules with 
C                 NUMBND bonds or the sequence number of the molecule 
C                 for which bond number NUMBND is specified. 
C                 See the table for NFCTR below.
C     NAM         number of atoms per molecule (NAM.GT.1), only used
C                 if NMOL.GT.1
C     NFCTR       determines the actions performed.
C                 NFCTR also influences the meanings of NUMBND and NMOL.
C                 (NFCTR.GE.NFCMIN.AND.NFCTR.LE.NFCMAX)
C
C          NUMBND means   NMOL means       vars delivered
C   ------------------------------------------------------------------------
C   L<NFALL>  num of bonds   num of mol.      EB,F,EB0,XB0
C   L<NFEF>   num of bonds   num of mol.      EB,F
C   L<NFEL>   num of bonds   num of mol.      EB,EB0,XB0
C   L<NFE>    num of bonds   num of mol.      EB
C   L<NFNDXL> seq num        seq num of mol.  EB(with length of NUMBND-th bond)
C   L<NFNDXE> seq num        seq num of mol.  EB(with energy of NUMBND-th bond)
C   ------------------------------------------------------------------------
C
C     IDIBND      Functional form selected for the potential, see above.
C                 (IDIBND.GT.IBNMIN.AND.IDIBND.LT.IBNMAX)
C     LEVERY
C         .FALSE.    quantities depending on BETA are only calculated
C                    at the first subroutine call.
C         .TRUE.     they are calculated at every subroutine call.
C     NUMBND      the meaning of NUMBND depends on NFCTR:
C                 it is either the number of bonds
C                 per molecule or the bond sequence number.
C                 See the table for NFCTR above.
C                 Note that if LDOPER.EQ.TRUE, NUMBND should be equal
C                 to NPBND
C     IB,JB(NUMBND) atom topology indices of the atoms forming the bonds
C                 ( accessed whatever LDOPER )
C     ICB(NUMBND)   bond-type code or pair code, denoting the bond index
C                 in the arrays CB and B0 for bond I-J. Only used if
C                 LDOPER.EQ..FALSE. 
C                 Otherwise the bond parameters are taken sequentially 
C                 from the arrays CBA,CBB,BA0,BB0(NPBND).
C     NBTY        number of bond types
C     CB(NBTY)    force constants - GROMOS96 functional form,
C                 (only accessed if LDOPER.EQ..FALSE.)
C     B0(NBTY)    bond lengths at minimum energy
C                 (only accessed if LDOPER.EQ..FALSE.,
C                  ignored if  IDIBND.EQ.IBND2)
C     NATTOT      the total number of atoms to be considered in the 
C                 X and F arrays.
C     X(NDIM*NATTOT)   atom cartesian coordinates.
C     F(NDIM*NATTOT)   delivered with the bond forces added if
C                 NFCTR.EQ.NFALL.OR.NFCTR.EQ.NFEF
C     EB          what is delivered in this variable depends on
C                 NFCTR according to the table below.
C
C   NFCTR                    EB delivered with
C   ------------------------------------------------------------------------
C   NFALL,NFEF,NFEL,NFE      total bond energy
C   NFNDXL                   length of the NUMBND-th
C                            bond of the NMOL-th molecule
C   NFNDXE                   energy of the NUMBND-th
C                            bond of the NMOL-th molecule
C   ------------------------------------------------------------------------
C
C     XB0(NUMBND*NMOL)   delivered with the bond lengths if
C                 NFCTR.EQ.NFEL.OR.NFCTR.EQ.NFALL.
C                 If LDO34.EQ..TRUE., see 3D-4D mixing rule above.
C     EB0(NUMBND*NMOL)   delivered with the bond energies.
C                 NFCTR.EQ.NFEL.OR.NFCTR.EQ.NFALL
C     LDOPER      determines whether a lamda dependent potential
C                 is calculated or not.
C            LDOPER .EQ. .FALSE.
C                 the potential V_stand is evaluated using ICB sequentially
C                 and the corresponding CB,B0.
C            LDOPER .EQ. .TRUE.
C                 LAMDA dependent potential V_pert is evaluated
C                 using CBA, CBB, BA0 and BB0 sequentially.
C     RLAM        StateA - StateB Hamiltonian coupling parameter lambda.
C                 (only accessed if LDOPER.EQ.TRUE)
C     NPBND       number of perturbed bonds.
C                 (only accessed if LDOPER.EQ..TRUE.)
C     CBA,CBB(NPBND) force constants - GROMOS96 functional form, 
C                 state A and B
C                 (only accessed if LDOPER.EQ..TRUE.)
C     BA0,BB0(NPBND) bond lengths at minimum energy - GROMOS96 functional 
C                 form, state A and B
C                 (only accessed if LDOPER.EQ..TRUE.)
C     EBGL        delivered with the value of the derivative
C                 V_pert (or V_pert34) with respect to RLAM, 
C                 if LDOPER.EQ..TRUE.
C     LDO34       determines whether a mu dependent potential
C                 is calculated or not.
C            LDO34 .EQ. .FALSE.
C                 the potential V_stand or V_pert are evaluated.
C            LDO34 .EQ. .TRUE.
C                 the potential V_stand34 or V_pert34 are evaluated.
C     RMUE        3D - 4D Hamiltonian coupling parameter mu.
C                 (only accessed if LDO34.EQ..TRUE.)
C     C4D(NATTOT) force constant for the 4th dimentional harmonic
C                 spring. Only used here as C4D(I).GT.EPSSM (EPSSM=-1.0E-20)
C                 implies that the atom is in 4D. A bond is treated in 4D
C                 if any of the atoms forming the bond is in 4D.
C                 (only accessed if LDO34.EQ..TRUE.).
C     EBGM        delivered with the value of the derivative
C                 V_stand34 or V_pert34 with respect to RMUE, if
C                 LDO34.EQ..TRUE.
C
C     ==--------------------------------------------------------------==
      SUBROUTINE COBOND (NCALCD,NMOL,NAM,NFCTR,IDIBND,LEVERY,
     $     NUMBND,IB,JB,ICB,
     $     NBTY,CB,B0,
     $     NATTOT,X,F,EB,XB0,EB0,
     $     LDOPER,RLAM,NPBND,CBA,CBB,BA0,BB0,EBGL,
     $     LDO34,RMUE,C4D,EBGM)
      IMPLICIT NONE
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'pertsz.h'
      INCLUDE 'cobond.h'
C args
      LOGICAL LDOPER,LEVERY
      INTEGER NUMBND,NCALCD,NMOL,NAM,NFCTR,IDIBND,NATTOT
      INTEGER IB(NUMBND),JB(NUMBND),ICB(NUMBND)
C bond types
      INTEGER NBTY
      real*8 CB(NBTY),B0(NBTY)
      real*8 X(NDIM*NATTOT),F(NDIM*NATTOT)
      real*8 XB0(NUMBND*NMOL),EB0(NUMBND*NMOL),EB
C perturbation vars
      INTEGER NPBND
      real*8 CBA(NPBND),BA0(NPBND),CBB(NPBND),BB0(NPBND),EBGL
      real*8 RLAM
C 3D-4D coupling stuff
      LOGICAL LDO34
      real*8 C4D(NATTOT),RMUE,EBGM

C local params
      real*8 ONE
      PARAMETER (ONE = 1.0d0)
C local vars
      LOGICAL LFIRST
      SAVE LFIRST

      LOGICAL LMONO,LOCTO,LVAC,LDOVIR,LDOTRA
      SAVE LMONO,LOCTO,LVAC,LDOVIR,LDOTRA

      real*8 CONV,BETAR,COSB,COSB2,BOXOH,BOXOQ
      SAVE COSB,COSB2,BOXOH,BOXOQ
C
      LOGICAL LFORCE,LBNDS
      INTEGER I3,J3,IATOFF,IBOFF,MOL,NBN,M,IPC
      INTEGER I,J
      real*8 XIJDST(MAXDIM),DSTTMP,RIJ2,RIJ
      real*8 DB,BA,BB,DC,DB0,DRB,DRB2,DBL,DBL2,DCL,DRBL,DRBL2
      real*8 EBH,EBGH,DF,XH,XH2,XH3
      real*8 EBH3,EBH4,EBGL3,EBGL4,DF3D,DF4D,RIJ24D
      real*8 BA2,B02,CBLA
      real*8 EPSSM
      CHARACTER PRGSTR*(6)
C data
      DATA EPSSM /-1.0D-20/
      DATA LFIRST/.TRUE./
      DATA PRGSTR /'COBOND'/
C begin
C******** CHECK ARGUMENTS FOR CONSISTENCY
C
      IF (IDIBND.LT.IBNMIN .OR. IDIBND.GT.IBNMAX ) THEN
         PRINT *,PRGSTR,': invalid choice of IDIBND' 
         CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF (NFCTR.LT.NFCMIN .OR. NFCTR.GT.NFCMAX) THEN
        PRINT *,PRGSTR,': invalid choice of NFCTR'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF (NCALCD.GT.NDIM) THEN
        PRINT *,PRGSTR,': NCALCD.GT.NDIM'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( NAM.LE.1 .AND. NMOL.GT.1 ) THEN
        PRINT *,PRGSTR,': NMOL.GT.1 requires NAM.GT.1'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF (LDOPER .AND. NUMBND.NE.NPBND) THEN
        PRINT *,PRGSTR,': LDOPER requires NUMBND.EQ.NPBND'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( LDO34 .AND. ( NDIM .NE. 4 .OR. NCALCD .NE. 4 ) ) THEN
         PRINT 500,PRGSTR
         PRINT *,'3D - 4D only possible if NDIM .EQ. 4 .AND. ',
     $           'NCALCD .EQ. 4'
         CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( LDO34 .AND. IDIBND .NE. IBND4 ) THEN
         PRINT 500,PRGSTR
         PRINT *,'3D - 4D only possible if IDIBND .EQ. ', IBND4
         CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C         
      IF ( LDO34 .AND. NTB.NE.NTBVAC .AND. BETA .NE. 90.0D0 ) THEN
         PRINT 500,PRGSTR
         PRINT *,'3D - 4D only possible if  BETA .EQ. 90.0 '
         CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( LDOPER.AND.IDIBND.EQ.IBND2 ) THEN
        PRINT 500,PRGSTR
        PRINT *,'LDOPER not implemented for IDIBND .EQ. IBND2'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( LDO34.AND.NFCTR.EQ.NFNDXE ) THEN
        PRINT 500,PRGSTR
        PRINT *,'3D - 4D incompatible with NFCTR .EQ. NFNDXE'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( (LDOPER.OR.LDO34).AND.NFCTR.EQ.NFNDXL ) THEN
        PRINT 500,PRGSTR
        PRINT *,'3D - 4D or LDOPER incompatible with NFCTR',
     $          ' .EQ. NFNDXL'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
C******** DONE CHECK ARGUMENTS FOR CONSISTENCY
C******** FIRST CALL SETTINGS
C      
      IF (LEVERY .OR. LFIRST) THEN
         LFIRST = .FALSE.
C
         CONV = 4.5D1/DATAN(ONE)
C
         LMONO = (NTB .GT. 0)
         LOCTO = (NTB. LT. 0)
         LVAC  = (NTB .EQ. 0)
         LDOVIR= (ABS(NTB) .EQ. NTBVIR)
         LDOTRA = .FALSE.
         IF (LMONO)THEN
            BETAR = BETA/CONV
            COSB = DCOS(BETAR)
            LDOTRA = (DABS(COSB).GE.1.D-4)
            COSB2 = 2.0D0*COSB
         ENDIF
C
         IF (LOCTO) THEN
            BOXOH = BOXH(1)
            BOXOQ = BOX(1)*0.75D0
         ELSE
            BOXOH = 0.0D0
            BOXOQ = 0.0D0
         ENDIF
      ENDIF
C
C******** DONE FIRST CALL SETTINGS
C
      EB = 0.0D0
      EBGL = 0.0D0
      EBGM = 0.0D0
C
C******** ONE BOND CALCULATION
C
      IF (NFCTR .EQ. NFNDXL .OR. NFCTR .EQ. NFNDXE) THEN
C
C* NFCTR.EQ.NFNDXL -> EB returned with length of NUMBND-th bond 
C  of the NMOL-th molecule
C  NFCTR.EQ.NFNDXE -> EB returned with energy of NUMBND-th bond
C  of the NMOL-th molecule
C
         IATOFF = (NMOL-1)*NAM -1
         I3 = NDIM*(IB(NUMBND)+IATOFF)
         J3 = NDIM*(JB(NUMBND)+IATOFF)
C
C** CALCULATE DISTANCE SQUARED, CORRECT FOR BOUNDARY CONDITION IF REQUIRED
C
!$OMP parallel do private(M)
         DO M=1,NCALCD
            XIJDST(M) = X(I3+M) - X(J3+M)
         ENDDO
         RIJ2 = 0.0D0
         IF (LVAC) THEN
            DO M=1,NCALCD
               RIJ2 = RIJ2 + XIJDST(M)**2
            ENDDO
         ELSE
            IF (LDOTRA) THEN
               CALL TRACO(1,0,XIJDST,BETA,1,LEVERY)
            ENDIF
            DO 550 M=1,NCALCD
               DSTTMP = XIJDST(M)
               IF (DSTTMP .GE. BOXH(M)) THEN
                  DSTTMP = DSTTMP - BOX(M)
               ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
                  DSTTMP = DSTTMP + BOX(M)
               ENDIF
               XIJDST(M) = DSTTMP
               RIJ2 = RIJ2 + DSTTMP**2
 550        CONTINUE
C
            IF (LOCTO) THEN
               DSTTMP = BOXOQ - DABS(XIJDST(1))-
     .              DABS(XIJDST(2))-DABS(XIJDST(3))
               IF (DSTTMP .LT. 0.0D0) THEN
                  RIJ2 = RIJ2 + DSTTMP*BOX(1)
                  DO M=1,NCALCD
                    XIJDST(M)=XIJDST(M)-SIGN(BOXOH,XIJDST(M)) 
                  ENDDO
               ENDIF
            ELSEIF (LDOTRA) THEN
               RIJ2=RIJ2+COSB2*XIJDST(1)*XIJDST(3)
            ENDIF
         ENDIF
C
C** DONE DISTANCE SQUARED CALCULATION
C
         IF (NFCTR .EQ. NFNDXL) THEN
C
C* if NFCTR.EQ.NFNDXL -> EB returned with length of NUMBND-th bond 
C  of the NMOL-th molecule
C
            EB = DSQRT(RIJ2)
         ELSE
C
C* if NFCTR.EQ.NFNDXE -> EB returned with energy of NUMBND-th bond
C  of the NMOL-th molecule
C
            IF (IDIBND .EQ. IBND2) THEN
C
C** PATH INTEGRAL FUNCTIONAL FORM
C
               IPC = ICB(NUMBND)
               EB  = 0.5D0*CB(IPC)*RIJ2
            ELSEIF (IDIBND .EQ. IBNDH) THEN
C
C** OLD HARMONIC FUNCTIONAL FORM
C
               RIJ = DSQRT(RIJ2)
               IF (LDOPER)THEN
                  BA    = BA0(NUMBND)
                  BB    = BB0(NUMBND)
                  CBLA  = CBA(NUMBND)*2.0D0*BA**2
                  DC    = CBB(NUMBND)*2.0D0*BB**2 - CBLA
                  DB0   = BB - BA
                  DRB   = RIJ - BA
                  DRB2  = DRB**2
                  DBL   = RLAM*DB0 + BA
                  DCL   = RLAM*DC + CBLA
                  DRBL  = RIJ - DBL
                  DRBL2 = DRBL**2
C
                  EB    = 0.5D0*(DCL*DRBL2-CBLA*DRB2)
                  EBGL  = 0.5D0*(DC *DRBL2-2.0D0*DCL*DRBL*DB0)
C
               ELSE
                  IPC   = ICB(NUMBND)
                  EB    = 0.5D0*(CB(IPC)*2.0D0*B0(IPC)**2) 
     .                  * (RIJ-B0(IPC))**2
               ENDIF
            ELSE
C
C** FOURTH POWER FUNCTIONAL FORM
C
               IF (LDOPER)THEN
                  BA    = BA0(NUMBND)
                  BB    = BB0(NUMBND)
                  CBLA  = CBA(NUMBND)
                  DC    = CBB(NUMBND) - CBLA
                  DB0   = BB - BA
                  DRB   = RIJ2 - BA**2
                  DRB2  = DRB**2
                  DBL   = RLAM*DB0 + BA
                  DBL2  = DBL**2
                  DCL   = RLAM*DC + CBLA
                  DRBL  = RIJ2 - DBL2
                  DRBL2 = DRBL**2
C
                  EB    = 0.25D0*(DCL*DRBL2-CBLA*DRB2 )
                  EBGL  = 0.25D0*(DC *DRBL2-4.0D0*DCL*DRBL*DBL*DB0 ) 
C
               ELSE
                  IPC   = ICB(NUMBND)
                  EB    = 0.25D0*CB(IPC)*(RIJ2-B0(IPC)**2)**2
               ENDIF           
            ENDIF
         ENDIF
C
         RETURN
      ENDIF
C
C******** DONE ONE BOND CALCULATION
C******** MULTIPLE BOND CALCULATION
C
      IF (NFCTR .EQ. NFALL) THEN
         LBNDS = .TRUE.
         LFORCE= .TRUE.
      ELSEIF (NFCTR .EQ. NFEF) THEN
         LBNDS = .FALSE.
         LFORCE = .TRUE.
      ELSEIF (NFCTR .EQ. NFEL) THEN
         LBNDS = .TRUE.
         LFORCE = .FALSE.
      ELSEIF (NFCTR .EQ. NFE) THEN
         LBNDS = .FALSE.
         LFORCE = .FALSE.
      ENDIF
C
      IATOFF = -1
      IBOFF = 0
C
C**** PATH INTEGRAL FUNCTIONAL FORM
C
      IF (IDIBND. EQ. IBND2) THEN

         DO 100 MOL=1,NMOL
            DO 200 NBN=1,NUMBND
               I3 = NDIM*(IB(NBN)+IATOFF)
               J3 = NDIM*(JB(NBN)+IATOFF)
C
C** CALCULATE DISTANCE SQUARED, CONNECTING VECTOR,
C   CORRECT FOR BOUNDARY CONDITION IF REQUIRED
C
!$OMP parallel do private(M)
               DO M=1,NCALCD
                  XIJDST(M) = X(I3+M) - X(J3+M)
               ENDDO
               RIJ2 = 0.0D0
               IF (LVAC) THEN
                  DO M=1,NCALCD
                     RIJ2 = RIJ2 + XIJDST(M)**2
                  ENDDO
               ELSE
                  IF (LDOTRA) THEN
                     CALL TRACO(1,0,XIJDST,BETA,1,LEVERY)
                  ENDIF
                  DO M=1,NCALCD
                     DSTTMP = XIJDST(M)
                     IF (DSTTMP .GE. BOXH(M)) THEN
                        DSTTMP = DSTTMP - BOX(M)
                     ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
                        DSTTMP = DSTTMP + BOX(M)
                     ENDIF
                     XIJDST(M) = DSTTMP
                     RIJ2 = RIJ2 + DSTTMP**2
                  ENDDO
C     
                  IF (LOCTO) THEN
                     DSTTMP = BOXOQ - DABS(XIJDST(1))-
     .                    DABS(XIJDST(2))-DABS(XIJDST(3))
                     IF (DSTTMP .LT. 0.0D0) THEN
                        RIJ2 = RIJ2 + DSTTMP*BOX(1)
                        DO M=1,NCALCD
                          XIJDST(M)=XIJDST(M)-SIGN(BOXOH,XIJDST(M)) 
                        ENDDO
                     ENDIF
                  ELSEIF (LDOTRA) THEN
                     RIJ2=RIJ2+COSB2*XIJDST(1)*XIJDST(3)
                  ENDIF
               ENDIF
C
C** DONE DISTANCE SQUARED CALCULATION
C   CALCULATE ENERGY/FORCE
C
               IPC = ICB(NBN)
               DF  = CB(IPC)
               EBH = 0.5D0*DF*RIJ2
               EB  = EB + EBH
C               
               IF (LBNDS) THEN
                  EB0(NBN+IBOFF) = EBH
                  XB0(NBN+IBOFF) = DSQRT(RIJ2)
               ENDIF
C
               IF (LFORCE) THEN
                  IF (LDOTRA) THEN
                     CALL TRACO(1,0,XIJDST,BETA,-1,LEVERY)
                  ENDIF
!$OMP parallel do private(M,XH)
                  DO M = 1,NCALCD
                     XH = XIJDST(M)*DF
                     F(I3+M) = F(I3+M) - XH
                     F(J3+M) = F(J3+M) + XH
                  ENDDO
               ENDIF
 200        CONTINUE
            IATOFF = IATOFF + NAM
            IBOFF  = IBOFF  + NUMBND
 100     CONTINUE
C
C**** DONE PATH INTEGRAL FUNCTIONAL FORM
C
      ELSEIF ( IDIBND.EQ.IBNDH ) THEN
C
C**** OLD HARMONIC FUNCTIONAL FORM
C
         IF (.NOT.LDO34 ) THEN
C
C** CALCULATION IN 3D (OLD HARMONIC FUNCTIONAL FORM)
C
            DO 1101 MOL=1,NMOL
               DO 1201 NBN=1,NUMBND
                  I3 = NDIM*(IB(NBN)+IATOFF)
                  J3 = NDIM*(JB(NBN)+IATOFF)
C
C** CALCULATE DISTANCE SQUARED, CONNECTING VECTOR,
C   CORRECT FOR BOUNDARY CONDITION IF REQUIRED
C
!$OMP parallel do private(M)
                  DO M=1,NCALCD
                     XIJDST(M) = X(I3+M) - X(J3+M)
                  ENDDO
                  RIJ2 = 0.0D0
                  IF (LVAC) THEN
                     DO M=1,NCALCD
                        RIJ2 = RIJ2 + XIJDST(M)**2
                     ENDDO
                  ELSE
                     DO M=1,NCALCD
                        DSTTMP = XIJDST(M)
                        IF (DSTTMP .GE. BOXH(M)) THEN
                           DSTTMP = DSTTMP - BOX(M)
                        ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
                           DSTTMP = DSTTMP + BOX(M)
                        ENDIF
                        XIJDST(M) = DSTTMP
                        RIJ2 = RIJ2 + DSTTMP**2
                     ENDDO
C     
                     IF (LOCTO) THEN
                        DSTTMP = BOXOQ - DABS(XIJDST(1))-
     .                       DABS(XIJDST(2))-DABS(XIJDST(3))
                        IF (DSTTMP .LT. 0.0D0) THEN
                           RIJ2 = RIJ2 + DSTTMP*BOX(1)
                           DO M=1,NCALCD
                             XIJDST(M)=XIJDST(M)-SIGN(BOXOH,XIJDST(M)) 
                           ENDDO
                        ENDIF
                     ELSEIF (LDOTRA) THEN
                        RIJ2=RIJ2+COSB2*XIJDST(1)*XIJDST(3) 
                     ENDIF
                  ENDIF
C
C** DONE DISTANCE CALCULATION
C   CALCULATE ENERGY/FORCES/LAMBDA DERIVATIVE
C
                  RIJ = DSQRT(RIJ2)
C we need it here     
                  IF (LDOPER)THEN
                     BA    = BA0(NBN)
                     BB    = BB0(NBN)
                     CBLA  = CBA(NBN)*2.0D0*BA**2
                     DC    = CBB(NBN)*2.0D0*BB**2 - CBLA
                     DB0   = BB - BA
                     DRB   = RIJ - BA
                     DRB2  = DRB**2
                     DBL   = RLAM*DB0 + BA
                     DCL   = RLAM*DC + CBLA
                     DRBL  = RIJ - DBL
                     DRBL2 = DRBL**2
C     
                     EBH   = 0.5D0*(DCL*DRBL2 - CBLA*DRB2)
                     EB    = EB + EBH
                     EBGH  = 0.5D0*(DC*DRBL2 - 2.0D0*DCL*DRBL*DB0) 
                     EBGL  = EBGL + EBGH
                     DF    = 1.0D0/RIJ * (DCL*DRBL - CBLA*DRB)
C     
                  ELSE
C
                     IPC   = ICB(NBN)
                     BA    = B0(IPC)
                     CBLA  = (CB(IPC)*2.0D0*BA**2)
                     DB    = RIJ - BA
                     DF    = CBLA * DB 
                     EBH   = 0.5D0 * DF * DB
                     EB    = EB + EBH
                     DF    = DF/RIJ
                  ENDIF
C                 
                  IF (LBNDS) THEN
                     EB0(NBN+IBOFF) = EBH
                     XB0(NBN+IBOFF) = RIJ
                  ENDIF
C
                  IF (LFORCE) THEN
                     IF (LDOTRA) THEN
                        CALL TRACO(1,0,XIJDST,BETA,-1,LEVERY)
                     ENDIF
!$OMP parallel do private(M,XH)
                     DO M = 1,NCALCD
                        XH = XIJDST(M)*DF
                        F(I3+M) = F(I3+M) - XH
                        F(J3+M) = F(J3+M) + XH
                     ENDDO
                  ENDIF
 1201          CONTINUE
C             
               IATOFF = IATOFF + NAM
               IBOFF  = IBOFF  + NUMBND
 1101       CONTINUE
C
C** DONE CALCULATION IN 3D (OLD HARMONIC FUNCTIONAL FORM)
C
         ELSE
            PRINT *,PRGSTR,'3D-4D not implemented IDIBND.EQ.IBNDH'
            CALL STOPGM('COBOND','FATAL ERROR')
         ENDIF
C
C**** DONE OLD HARMONIC FUNCTIONAL FORM
C**** FOURTH POWER FUNCTIONAL FORM
C
      ELSEIF ( .NOT.LDO34 ) THEN
C
C** CALCULATION IN 3D (FOURTH POWER FUNCTIONAL FORM)
C
         DO 101 MOL=1,NMOL
            DO 201 NBN=1,NUMBND
               I3 = NDIM*(IB(NBN)+IATOFF)
               J3 = NDIM*(JB(NBN)+IATOFF)
C
C** CALCULATE DISTANCE SQUARED, CONNECTING VECTOR,
C   CORRECT FOR BOUNDARY CONDITION IF REQUIRED
C
!$OMP parallel do private(M)
               DO M=1,NCALCD
                  XIJDST(M) = X(I3+M) - X(J3+M)
               ENDDO
               RIJ2 = 0.0D0
               IF (LVAC) THEN
                  DO M=1,NCALCD
                     RIJ2 = RIJ2 + XIJDST(M)**2
                  ENDDO
               ELSE
                  IF (LDOTRA) THEN
                     CALL TRACO(1,0,XIJDST,BETA,1,LEVERY)
                  ENDIF
                  DO M=1,NCALCD
                     DSTTMP = XIJDST(M)
                     IF (DSTTMP .GE. BOXH(M)) THEN
                        DSTTMP = DSTTMP - BOX(M)
                     ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
                        DSTTMP = DSTTMP + BOX(M)
                     ENDIF
                     XIJDST(M) = DSTTMP
                     RIJ2 = RIJ2 + DSTTMP**2
                  ENDDO
C     
                  IF (LOCTO) THEN
                     DSTTMP = BOXOQ - DABS(XIJDST(1))-
     .                    DABS(XIJDST(2))-DABS(XIJDST(3))
                     IF (DSTTMP .LT. 0.0D0) THEN
                        RIJ2 = RIJ2 + DSTTMP*BOX(1)
                        DO M=1,NCALCD
                          XIJDST(M)=XIJDST(M)-SIGN(BOXOH,XIJDST(M)) 
                        ENDDO
                     ENDIF
                  ELSEIF (LDOTRA) THEN
                     RIJ2=RIJ2+COSB2*XIJDST(1)*XIJDST(3)
                  ENDIF
               ENDIF
C
C** DONE DISTANCE CALCULATION
C   CALCULATE ENERGY/FORCES/LAMBDA DERIVATIVE
C
               IF (LDOPER) THEN
                  BA    = BA0(NBN)
                  BB    = BB0(NBN)
                  CBLA  = CBA(NBN)
                  DC    = CBB(NBN) - CBLA
                  DB0   = BB - BA
                  DRB   = RIJ2 - BA**2
                  DRB2  = DRB**2
                  DBL   = RLAM*DB0 + BA
                  DBL2  = DBL**2
                  DCL   = RLAM*DC + CBLA
                  DRBL  = RIJ2 - DBL2
                  DRBL2 = DRBL**2
C     
                  EBH   = 0.25D0*(DCL*DRBL2 - CBLA*DRB2)
                  EB    = EB + EBH
                  EBGH  = 0.25D0*(DC*DRBL2 - 4.0D0*DCL*DRBL*DBL*DB0) 
                  EBGL  = EBGL + EBGH
                  DF    = DCL*DRBL - CBLA*DRB
               ELSE
                  IPC   = ICB(NBN)
                  DB    = RIJ2 - B0(IPC)**2
                  DF    = CB(IPC)*DB
                  EBH   = 0.25D0*DF*DB
                  EB    = EB + EBH
               ENDIF
C
               IF (LBNDS) THEN
                  EB0(NBN+IBOFF) = EBH
                  XB0(NBN+IBOFF) = DSQRT(RIJ2)
               ENDIF
C
               IF (LFORCE) THEN
                  IF (LDOTRA) THEN
                     CALL TRACO(0,1,XIJDST,BETA,-1,LEVERY)
                  ENDIF
!$OMP parallel do private(M,XH)
                  DO M = 1,NCALCD
                     XH = XIJDST(M)*DF
                     F(I3+M) = F(I3+M) - XH
                     F(J3+M) = F(J3+M) + XH
                  ENDDO
               ENDIF
 201        CONTINUE
C
            IATOFF = IATOFF + NAM
            IBOFF  = IBOFF  + NUMBND
 101     CONTINUE
C
C** DONE CALCULATION IN 3D (FOURTH POWER FUNCTIONAL FORM)
C
      ELSE
C
C** CALCULATION WITH 3D-4D COUPLING (FOURTH POWER FUNCTIONAL FORM)
C
         DO 102 MOL=1,NMOL
            DO 202 NBN=1,NUMBND
               I = IB(NBN)+IATOFF
               J = JB(NBN)+IATOFF
               IF ( C4D(I+1) .GE. EPSSM .OR. C4D(J+1) .GE. EPSSM ) THEN
                 I3 = NDIM*I
                 J3 = NDIM*J
C
C** CALCULATE DISTANCE SQUARED, CONNECTING VECTOR,
C   CORRECT FOR BOUNDARY CONDITION IF REQUIRED, IN 3D AND 4D
C
                 XIJDST(1) = X(I3+1) - X(J3+1)
                 XIJDST(2) = X(I3+2) - X(J3+2)
                 XIJDST(3) = X(I3+3) - X(J3+3)
                 XIJDST(4) = X(I3+4) - X(J3+4)
                 RIJ2 = 0.0D0
                 IF (LVAC) THEN
                   RIJ2 = XIJDST(1)*XIJDST(1)
     .                  + XIJDST(2)*XIJDST(2)
     .                  + XIJDST(3)*XIJDST(3)
                   RIJ24D = RIJ2 + XIJDST(4)**2
                 ELSE
C no monoclinic in 4D
cmb                   DO 652 M=1,4
cmb                     DSTTMP = XIJDST(M)
cmb                     IF (DSTTMP .GE. BOXH(M)) THEN
cmb                       DSTTMP = DSTTMP - BOX(M)
cmb                     ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
cmb                       DSTTMP = DSTTMP + BOX(M)
cmb                     ENDIF
cmb                     XIJDST(M) = DSTTMP
cmb                     IF ( M .LE. 3 ) THEN
cmb                       RIJ2 = RIJ2 + DSTTMP**2
cmb                     ELSE
cmb                       RIJ24D = RIJ2 + XIJDST(4)**2
cmb                     ENDIF
cmb 652               CONTINUE
cmb                   
                   DO M=1,3
                     DSTTMP = XIJDST(M)
                     IF (DSTTMP .GE. BOXH(M)) THEN
                       DSTTMP = DSTTMP - BOX(M)
                     ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
                       DSTTMP = DSTTMP + BOX(M)
                     ENDIF
                     XIJDST(M) = DSTTMP
                     RIJ2 = RIJ2 + DSTTMP*DSTTMP
                   ENDDO
cmb                M=4
                   DSTTMP = XIJDST(4)
                   IF (DSTTMP .GE. BOXH(4)) THEN
                     DSTTMP = DSTTMP - BOX(4)
                   ELSEIF (DSTTMP .LT. -BOXH(4)) THEN
                     DSTTMP = DSTTMP + BOX(4)
                   ENDIF
                   XIJDST(4) = DSTTMP
                   RIJ24D = RIJ2 + XIJDST(4)*XIJDST(4)
cmb
                   IF (LOCTO) THEN
                     DSTTMP = BOXOQ - DABS(XIJDST(1))-
     .                    DABS(XIJDST(2))-DABS(XIJDST(3))
                     IF (DSTTMP .LT. 0.0D0) THEN
                       RIJ2 = RIJ2 + DSTTMP*BOX(1)
                       DO M=1,NCALCD
                         XIJDST(M)=XIJDST(M)-SIGN(BOXOH,XIJDST(M)) 
                       ENDDO
                     ENDIF
                   ENDIF
C no monoclinic in 4D                     
                 ENDIF
C     
C** DONE DISTANCE CALCULATION
C   CALCULATE ENERGY/FORCES/LAMBDA/MU DERIVATIVE
C
                 IF (LDOPER)THEN
                   BA    = BA0(NBN)
                   BB    = BB0(NBN)
                   BA2   = BA**2
                   CBLA  = CBA(NBN)
                   DC    = CBB(NBN) - CBLA
                   DB0   = BB - BA
                   DBL   = RLAM*DB0 + BA
                   DBL2  = DBL**2
                   DCL   = RLAM*DC + CBLA
C* 3D
                   DRB   = RIJ2 - BA2
                   DRB2  = DRB**2
                   DRBL  = RIJ2 - DBL2
                   DRBL2 = DRBL**2
C
                   EBH3  =  0.25D0*(DCL*DRBL2-CBLA*DRB2)
                   EBGL3 =  0.25D0*(DC *DRBL2-4.0D0*DCL*DRBL*DBL*DB0) 
                   DF3D  = -RMUE * (DCL*DRBL -CBLA*DRB)
C* 4D
                   DRB   = RIJ24D - BA2
                   DRB2  = DRB**2
                   DRBL  = RIJ24D - DBL2
                   DRBL2 = DRBL**2
C     
                   EBH4  =  0.25D0*(DCL*DRBL2-CBLA*DRB2) 
                   EBGL4 =  0.25D0*(DC *DRBL2-4.0D0*DCL*DRBL*DBL*DB0) 
                   DF4D  =  RMUE * (DCL*DRBL  - CBLA*DRB)
C* sum
                   EB    = EB   + RMUE * ( EBH4  - EBH3 )
                   EBGL  = EBGL + RMUE * ( EBGL4 - EBGL3 )
                   EBGM  = EBGM + EBH4 - EBH3
C
                 ELSE
C
                   IPC   = ICB(NBN)
                   B02   = B0(IPC)**2
C* 3D
                   DB    = RIJ2 - B02
                   DF3D  = CB(IPC)*DB
C
                   EBH3  = 0.25D0*DF3D*DB
                   DF3D  = -DF3D*RMUE
C* 4D
                   DB    = RIJ24D - B02
                   DF4D  = CB(IPC)*DB
C
                   EBH4  = 0.25D0*DF4D*DB
                   DF4D  = DF4D*RMUE
C* sum
                   EB    = EB   + RMUE * ( EBH4  - EBH3 )
                   EBGM  = EBGM + EBH4 - EBH3
                 ENDIF
C                  
                 IF (LBNDS) THEN
                   EB0(NBN+IBOFF) = RMUE * ( EBH4 - EBH3 )
                   XB0(NBN+IBOFF) = 
     $                 DSQRT(RMUE*RIJ24D + (1.0D0-RMUE)*RIJ2)
                 ENDIF
C
                 IF (LFORCE) THEN
C non monoclinic in 4D
                   XH  = XIJDST(1)*DF3D
                   XH2 = XIJDST(2)*DF3D
                   XH3 = XIJDST(3)*DF3D
                   F(I3+1) = F(I3+1) - XH
                   F(I3+2) = F(I3+2) - XH2
                   F(I3+3) = F(I3+3) - XH3
                   F(J3+1) = F(J3+1) + XH 
                   F(J3+2) = F(J3+2) + XH2
                   F(J3+3) = F(J3+3) + XH3
!$OMP parallel do private(M,XH)
                   DO M = 1,NCALCD
                     XH = XIJDST(M)*DF4D
                     F(I3+M) = F(I3+M) - XH
                     F(J3+M) = F(J3+M) + XH
                   ENDDO
                 ENDIF
               ENDIF
 202         CONTINUE
C
             IATOFF = IATOFF + NAM
             IBOFF  = IBOFF  + NUMBND
C     
 102       CONTINUE
C
C** DONE CALCULATION WITH 3D-4D COUPLING (FOURTH POWER FUNCTIONAL FORM)
C**** DONE FOURTH POWER FUNCTIONAL FORM
C
      ENDIF
C
C******** DONE MULTIPLE BOND CALCULATION
C
C     end cobond
 500  FORMAT (1X,A6,': IMPLEMENTATION RESTRICTION')
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
COMMSUBR ANGLE
C W.F. VAN GUNSTEREN, GRONINGEN, JUNE 1985
C The GROMOS96 functional form versions supplied by Thomas Beutler.
C Revised/documented/GROMOS87 compatibility P.Hunenberger, sept. 95
C
C      SUBROUTINE ANGLE(NCALCD,NMOL,NAM,NFCTR,IDIANG,LEVERY,
C     $     NUMANG,IT,JT,KT,ICT,
C     $     NANGTY,CT,T0,
C     $     NATTOT,X,F,EBANG,XT0,ET0,
C     $     LDOPER,RLAM,NPANG,CTA,CTB,TA0,TB0,EBANGL,
C     $     LDO34,RMUE,C4D,EBANGM)
C
C
C     ANGLE will supply the covalent bond angle (T) energies and forces
C     according to a potential energy function depending on IDIANG
C     for angles in NCALCD dimensions (NCALCD .LE. NDIM).
C     Peridoc boundary conditions can be taken into account
C     when NTB .NE. L<NTBVAC> (with some restrictions mentioned below
C     for the 4D case)
C     More than one identical molecule can be considered.
C
C     The cobond.h header file defines reference values for the
C     NFCTR and IDIANG switches.
C     ANGLE assumes that NDIM, NTB, BOX(3), BOXH(3) and BETA are 
C     defined and transmitted through the common blocks of the box.h 
C     header file.
C     The file pertsz.h (perturbation) is optional (LDOPER.EQ.TRUE.)
C
C     When LDOPER. is .FALSE. the standard potential,
C     V_stand(T;T0,CT) = sum over all specified bond angles of
C     1.       0.5 * CT * (COS(T)-T0)**2.         if IDIANG = L<IANGC>
C              (intended for covalent bond angles, GROMOS96)
C     2.       0.5 * F(CT) * (T-ACOS(T0))**2      if IDIANG = L<IANGH>
C              (limited implementation, intended for compatibility
C               with GROMOS87 program)
C     is evaluated, where T0 denotes the cosine of the bond angle at 
C     minimum energy and CT the force constant (according to GROMOS96 
C     functional form). F(CT) is the inverse of the g function 
C     (see manual table 2.5.3.1)
C
C     When LDOPER is .TRUE., a LAMDA (RLAM) dependent perturbation
C     potential, 
C     V_pert(T;RLAM;CTA,CTB,TA0,TB0) = sum over all specified angles of
C     1.       0.5 * [ (1-RLAM)*CTA+RLAM*CTB) ] *
C                    [ COS(T) - ((1-RLAM)*TA0+RLAM*TB0) ]**2
C             -0.5 * CTA * (COS(T)-TA0)**2
C             if IDIANG = L<IANGC>
C     2.       0.5 * [ (1-RLAM)*F(CTA)+RLAM*F(CTB)) ] *
C                    [ T - ((1-RLAM)*ACOS(TA0)+RLAM*ACOS(TB0)) ]**2
C             -0.5 * F(CTA) * (T-ACOS(TA0))**2
C             if IDIANG = L<IANGH>
C     is evaluated, where TA0/TB0 denote the bond angle cosines at
C     minimum energy in states A/B and CTA/CTB the force constants
C     (according to GROMOS96 functional form) in states A/B.
C     The second term in 1. and 2. is intented to cancel the
C     V_pert(RLAM=0) contribution already taken into account by a previous
C     call to angle (all angles, LDOPER=.FALSE.) i.e. to avoid double 
C     counting. dV_pert/d(RLAM) (EBGL) is also evaluated.
C     It is not implemented for NFCTR.EQ.NFNDXL.
C
C     If LDO34 is .TRUE., a MU (RMUE) dependent 3D-4D perturbation 
C     potential,
C     V_stand34(T;RMUE;CT,T0) = sum over all specified angles of
C             + RMUE * V(T_4D;CT,T0) - RMUE * V(T_3D;CT,T0) 
C             if LDOPER.EQ..FALSE.
C     or
C     V_pert34(T;RMUE;RLAM;CTA,CTB,TA0,TB0) = sum over all specified angles of
C             + RMUE * V_pert(T_4D;RLAM;CTA,CTB,TA0,TB0)
C             - RMUE * V_pert(T_3D;RLAM;CTA,CTB,TA0,TB0)
C     is evaluated.
C
C     Note that the contribution of the 3D component is -RMUE and 
C     not (1-RMUE) to avoid double counting of the V_34(RMUE=0) component
C     already taken into account by previous call(s) to angle.
C     dV_pert/d(RLAM) (EBGL) and dV_pert/d(RMUE) (EBGM) are also evaluated.
C     LDO34 requires NDIM.EQ.4 .AND. NCALCD.EQ.4.
C     It is not implemented for NFCTR.EQ.NFNDXL.OR.NFCTR.EQ.NFNDXE.
C     It is only implemented for IDIBND.EQ.IBNDC.
C     Only a rectangular or truncated-octahedron box (BETA.EQ.90) can be 
C     used. In the latter case, only rectangular periodic gathering is
C     applied in the NCALCD-3 last dimensions.
C     Bond angle values are returned in a mixed 3D-4D RMUE dependant form
C     too, if required:
C     ANGLE = ACOS( RMUE*COSINE(4D) + (1.0-RMUE)*COSINE(3D) )
C
C     Arguments:
C
C     NCALCD      the number of dimensions in which to consider
C                 all vectors, coordinate/forces (NCALCD.LE.NDIM).
C     NMOL        The meaning of this argument depends on NFCTR: it
C                 is either the number of identical molecules with 
C                 NUMANG angles or the sequence number of the molecule 
C                 for which angle number NUMANG is specified.
C                 See the table for NFCTR below.
C     NAM         number of atoms per molecule (NAM.GT.2), only used
C                 if NMOL.GT.1
C     NFCTR       determines the actions performed.
C                 NFCTR also influences the meanings of NUMBND and NMOL.
C                 (NFCTR.GE.NFCMIN.AND.NFCTR.LE.NFCMAX)
C
C          NUMANG means   NMOL means       vars delivered
C   ------------------------------------------------------------------------
C   L<NFALL>  num of angles  num of mol.      EBANG,F,ET0,XT0
C   L<NFEF>   num of angles  num of mol.      EBANG,F
C   L<NFEL>   num of angles  num of mol.      EBANG,ET0,XT0
C   L<NFE>    num of angles  num of mol.      EBANG
C   L<NFNDXL> seq num        seq num of mol.  EBANG(with value NUMANG-th angle)
C   L<NFNDXE> seq num        seq num of mol.  EBANG(with energy NUMANG-th angle)
C   ------------------------------------------------------------------------
C
C     IDIANG      Functional form selected for the potential, see above.
C                 (IDIANG.GT.IANMIN.AND.IDIANG.LT.IANMAX)
C     LEVERY
C         .FALSE.    quantities depending on BETA are only calculated
C                    at the first subroutine call.
C         .TRUE.     they are calculated at every subroutine call.
C     NUMANG      the meaning of NUMANG depends on NFCTR:
C                 it is either the number of angles
C                 per molecule or the angle sequence number.
C                 See the table for NFCTR above. 
C                 Note that if LDOPER.EQ.TRUE, NUMANG should be equal
C                 to NPANG
C     IT,JT,KT(NUMANG) atom topology indices of the atoms forming the angles.
C                 ( accessed whatever LDOPER )
C     ITB(NUMANG)   angle-type code denoting the angle index
C                 in the arrays CT and T0 for angle I-J-K. Only used if
C                 LDOPER.EQ..FALSE. 
C                 Otherwise the angle parameters are taken sequentially 
C                 from the arrays CTA,CTB,TA0,TB0(NPANG).
C     NANGTY      number of angle types
C     CT(NANGTY)  force constants - GROMOS96 functional form,
C                 (only accessed if LDOPER.EQ..FALSE.)
C     T0(NANGTY)  angle cosines at minimum energy
C                 (only accessed if LDOPER.EQ..FALSE.,
C                  ignored if  IDIBND.EQ.IBND2)
C     NATTOT      the total number of atoms to be considered in the 
C                 X and F arrays.
C     X(NDIM*NATTOT)   atom cartesian coordinates.
C     F(NDIM*NATTOT)   delivered with the angle forces added if
C                 NFCTR.EQ.NFALL.OR.NFCTR.EQ.NFEF
C     EBANG       what is delivered in this variable depends on
C                 NFCTR according to the table below.
C
C   NFCTR                    EBANG delivered with
C   ------------------------------------------------------------------------
C   NFALL,NFEF,NFEL,NFE      total angle energy
C   NFNDXL                   value of the NUMANG-th
C                            angle of the NMOL-th molecule
C   NFNDXE                   energy of the NUMANG-th
C                            angle of the NMOL-th molecule
C   ------------------------------------------------------------------------
C
C
C     XT0(NUMANG*NMOL)   delivered with the angle values if
C                 NFCTR.EQ.NFEL.OR.NFCTR.EQ.NFALL.
C                 If LDO34.EQ..TRUE., see 3D-4D mixing rule above.
C     ET0(NUMANG*NMOL)   delivered with the angle energies.
C                 NFCTR.EQ.NFEL.OR.NFCTR.EQ.NFALL
C     LDOPER      determines whether a lamda dependent potential
C                 is calculated or not.
C            LDOPER .EQ. .FALSE.
C                 the potential V_stand is evaluated using ITB sequentially
C                 and the corresponding CT,T0.
C            LDOPER .EQ. .TRUE.
C                 LAMDA dependent potential V_pert is evaluated
C                 using CTA, CTB, TA0 and TB0 sequentially.
C     RLAM        StateA - StateB Hamiltonian coupling parameter lambda.
C                 (only accessed if LDOPER.EQ.TRUE)
C     NPANG       number of perturbed angles.
C                 (only accessed if LDOPER.EQ..TRUE.)
C     CTA,CTB(NPANG) force constants - GROMOS96 functional form, 
C                 state A and B
C                 (only accessed if LDOPER.EQ..TRUE.)
C     TA0,TB0(NPBND) bond angle cosines at minimum energy - GROMOS96 
C                 functional form, state A and B
C                 (only accessed if LDOPER.EQ..TRUE.)
C     EBGL        delivered with the value of the derivative
C                 V_pert (or V_pert34) with respect to RLAM, 
C                 if LDOPER.EQ..TRUE.
C     LDO34       determines whether a mu dependent potential
C                 is calculated or not.
C            LDO34 .EQ. .FALSE.
C                 the potential V_stand or V_pert are evaluated.
C            LDO34 .EQ. .TRUE.
C                 the potential V_stand34 or V_pert34 are evaluated.
C     RMUE        3D - 4D Hamiltonian coupling parameter mu.
C                 (only accessed if LDO34.EQ..TRUE.)
C     C4D(NATTOT) force constant for the 4th dimentional harmonic
C                 spring. Only used here as C4D(I).GT.EPSSM (EPSSM=-1E-20)
C                 implies that the atom is in 4D. An angle is treated in 4D
C                 if any of the atoms forming the angle is in 4D.
C                 (only accessed if LDO34.EQ..TRUE.).
C     EBGM        delivered with the value of the derivative
C                 V_stand34 or V_pert34 with respect to RMUE, if
C                 LDO34.EQ..TRUE.
COMMEND
C     ==--------------------------------------------------------------==
      SUBROUTINE ANGLE(NCALCD,NMOL,NAM,NFCTR,IDIANG,LEVERY,
     $                 NUMANG,IT,JT,KT,ICT,NANGTY,CT,T0,
     $                 NATTOT,X,F,EBANG,XT0,ET0,LDOPER,
     $                 RLAM,NPANG,CTA,CTB,TA0,TB0,EBANGL,
     $                 LDO34,RMUE,C4D,EBANGM)
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'pertsz.h'
      INCLUDE 'cobond.h'
C args
      LOGICAL LDOPER,LEVERY
      INTEGER NUMANG,NCALCD,NFCTR,NMOL,NAM,NANGTY,IDIANG,NATTOT
      INTEGER IT(NUMANG),JT(NUMANG),KT(NUMANG),ICT(NUMANG)
      real*8 CT(NANGTY),T0(NANGTY),EBANG,EBANGL
      real*8 X(NDIM*NATTOT),F(NDIM*NATTOT)
      real*8 XT0(NUMANG*NMOL),ET0(NUMANG*NMOL)
C perturbation stuff
      INTEGER NPANG
      real*8 CTA(NPANG),CTB(NPANG),TA0(NPANG),TB0(NPANG)
      real*8 RLAM
C 3D-4D coupling stuff
      LOGICAL LDO34
      real*8 C4D(NATTOT),RMUE,EBANGM

C local params
      real*8 ONE
      PARAMETER (ONE = 1.0D0)
C local vars
      LOGICAL LFIRST
      SAVE LFIRST

      LOGICAL LMONO,LOCTO,LVAC,LDOVIR,LDOTRA
      SAVE LMONO,LOCTO,LVAC,LDOVIR,LDOTRA

      real*8 BETAR,COSB,COSB2,BOXOH,BOXOQ,CONV
      SAVE COSB,COSB2,BOXOH,BOXOQ,CONV
C
      LOGICAL LBNDS,LFORCE
      CHARACTER PRGSTR*(5)
      INTEGER I3,J3,K3,IATOFF,IBOFF,M,IC,NNA,MOL
      INTEGER I,J,K
      real*8 BOXM,BOXHM,RIJ2,RKJ2,XTMP,DSTTMP
      real*8 XXIJ(MAXDIM),XXKJ(MAXDIM)
      real*8 CST,AT,RKJ,RIJ,RIJRKJ
      real*8 DCB,DB0,DBA,DF,DB2,DB3,EBAH,DB4,EBGH,DA
      real*8 DB1,DBL,DB5
      real*8 DFI,DFJ,DFK
      real*8 RIJ4D,RKJ4D,RIJ24D,RKJ24D,RIJRK4,CST4D
      real*8 EBAH3,EBAH4,DF3D,DF4D
      real*8 EBGLH3,EBGLH4
      real*8 TREFA,TREFB,CTTA,CTTB
      real*8 EPSSM,EPSSIN
C
      real*8 TREF,CSQ,COST2,CTHA,SDET,SNT
      real*8 EKT
      PARAMETER (EKT = 2.5d0)
C     Ekt for bond angle force constant conversion
C     2.5 [kJ/mol] at 300K
C data
C
      DATA EPSSM /-1.0d-20/
      DATA EPSSIN /1.0d-12/
      DATA LFIRST/.TRUE./
      DATA PRGSTR /'ANGLE'/
C begin
C******** CHECK ARGUMENTS FOR CONSISTENCY
C

      IF (IDIANG.LT.IANMIN .OR. IDIANG.GT.IANMAX ) THEN
         PRINT *,PRGSTR,': invalid choice of IDIANG'
         CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF (NFCTR.LT.NFCMIN .OR. NFCTR.GT.NFCMAX) THEN
        PRINT *,PRGSTR,': invalid choice of NFCTR'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF (NCALCD.GT.NDIM) THEN
        PRINT *,PRGSTR,': NCALCD.GT.NDIM'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( NAM.LE.2 .AND. NMOL.GT.1 ) THEN
        PRINT *,PRGSTR,': NMOL.GT.1 requires NAM.GT.2'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF (LDOPER .AND. NUMANG.NE.NPANG) THEN
        PRINT *,PRGSTR,': LDOPER requires NUMANG.EQ.NPANG'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( LDO34 .AND. ( NDIM .NE. 4 .OR. NCALCD .NE. 4 ) ) THEN
         PRINT 500,PRGSTR
         PRINT *,'3D - 4D only possible if NDIM .EQ. 4 .AND. ',
     $           'NCALCD .EQ. 4'
         CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( LDO34 .AND. IDIANG .NE. IANGC ) THEN
         PRINT 500,PRGSTR
         PRINT *,'3D - 4D only possible if IDIANG .EQ. IANGC'
         CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( LDO34 .AND. NTB.NE.NTBVAC .AND. BETA .NE. 90.0D0 ) THEN
         PRINT 500,PRGSTR
         PRINT *,'3D - 4D only possible if  BETA .EQ. 90.0 '
         CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( LDO34.AND.NFCTR.EQ.NFNDXE ) THEN
        PRINT 500,PRGSTR
        PRINT *,'3D - 4D incompatible with NFCTR .EQ. NFNDXE'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( (LDOPER.OR.LDO34).AND.NFCTR.EQ.NFNDXL ) THEN
        PRINT 500,PRGSTR
        PRINT *,'3D - 4D and LDOPER incompatible with NFCTR',
     $          ' .EQ. NFNDXL'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
C******** DONE CHECK ARGUMENTS FOR CONSISTENCY
C******** FIRST CALL SETTINGS
C 
      IF (LEVERY .OR. LFIRST) THEN
         LFIRST = .FALSE.
C
         CONV = 4.5D1/DATAN(ONE)
         LMONO = (NTB .GT. 0)
         LOCTO = (NTB. LT. 0)
         LVAC  = (NTB .EQ. 0)
         LDOVIR= (ABS(NTB) .EQ. NTBVIR)
         LDOTRA = .FALSE.
         IF (LMONO)THEN
            BETAR = BETA/CONV
            COSB = DCOS(BETAR)
            LDOTRA = (DABS(COSB).GE.1.d-4)
            COSB2 = 2.0d0*COSB
         ENDIF
C
         IF (LOCTO) THEN
            BOXOH = BOXH(1)
            BOXOQ = BOX(1)*0.75d0
         ELSE
            BOXOH = 0.0d0
            BOXOQ = 0.0d0
         ENDIF
      ENDIF
C
C******** DONE FIRST CALL SETTINGS
C
      EBANG = 0.0d0
      EBANGL = 0.0d0
      EBANGM = 0.0d0
C
C******** ONE ANGLE CALCULATION
C
      IF (NFCTR .EQ. NFNDXL .OR. NFCTR .EQ. NFNDXE) THEN
C
C* NFCTR.EQ.NFNDXL -> EBANG returned with value of NUMANG-th angle
C  of the NMOL-th molecule
C  NFCTR.EQ.NFNDXE -> EBANG returned with energy of NUMANG-th angle
C  of the NMOL-th molecule
C
        IATOFF = (NMOL-1)*NAM -1
        I3 = NDIM*(IT(NUMANG)+IATOFF)
        J3 = NDIM*(JT(NUMANG)+IATOFF)
        K3 = NDIM*(KT(NUMANG)+IATOFF)
C
C** CALCULATE ANGLE COSINE, CORRECT FOR BOUNDARY CONDITION IF REQUIRED
C
!$OMP parallel do private(M,XTMP)
         DO M = 1,NCALCD
            XTMP = X(J3+M)
            XXIJ(M) = X(I3+M) - XTMP
            XXKJ(M) = X(K3+M) - XTMP
         ENDDO
         RIJ2 = 0.0d0
         RKJ2 = 0.0d0
         IF (LVAC) THEN
            DO M=1,NCALCD
               RIJ2 = RIJ2 + XXIJ(M)**2
               RKJ2 = RKJ2 + XXKJ(M)**2
            ENDDO
         ELSE
            IF (LDOTRA) THEN
               CALL TRACO(1,0,XXIJ,BETA,1,LEVERY)
               CALL TRACO(1,0,XXKJ,BETA,1,LEVERY)
            ENDIF
            DO M=1,NCALCD
               BOXM = BOX(M)
               BOXHM= BOXH(M)
C
               DSTTMP = XXIJ(M)
               IF (DSTTMP .GE. BOXHM) THEN
                  DSTTMP = DSTTMP - BOXM
               ELSEIF (DSTTMP .LT. -BOXHM) THEN
                  DSTTMP = DSTTMP + BOXM
               ENDIF
               XXIJ(M) = DSTTMP
               RIJ2 = RIJ2 + DSTTMP*DSTTMP
C
               DSTTMP = XXKJ(M)
               IF (DSTTMP .GE. BOXHM) THEN
                  DSTTMP = DSTTMP - BOXM
               ELSEIF (DSTTMP .LT. -BOXHM) THEN
                  DSTTMP = DSTTMP + BOXM
               ENDIF
               XXKJ(M) = DSTTMP
               RKJ2 = RKJ2 + DSTTMP*DSTTMP
            ENDDO
C
            IF (LOCTO) THEN
               DSTTMP = BOXOQ - DABS(XXIJ(1))-
     .            DABS(XXIJ(2))-DABS(XXIJ(3))
               IF (DSTTMP .LT. 0.0d0) THEN
                  RIJ2 = RIJ2 + DSTTMP*BOX(1)
                  DO M=1,NCALCD
                    XXIJ(M)=XXIJ(M)-SIGN(BOXOH,XXIJ(M)) 
                  ENDDO
               ENDIF
C
               DSTTMP = BOXOQ - DABS(XXKJ(1))-
     .            DABS(XXKJ(2))-DABS(XXKJ(3))
               IF (DSTTMP .LT. 0.0d0) THEN
                  RKJ2 = RKJ2 + DSTTMP*BOX(1)
                  DO M=1,NCALCD
                    XXKJ(M)=XXKJ(M)-SIGN(BOXOH,XXKJ(M)) 
                  ENDDO
               ENDIF
            ELSEIF (LDOTRA) THEN
               RIJ2=RIJ2+COSB2*XXIJ(1)*XXIJ(3)
               RKJ2=RKJ2+COSB2*XXKJ(1)*XXKJ(3)
            ENDIF
         ENDIF
C
         RKJ = DSQRT(RKJ2)
         RIJ = DSQRT(RIJ2)
         RIJRKJ = 0.0d0
         DO M=1, NCALCD
            RIJRKJ = RIJRKJ + XXIJ(M)*XXKJ(M)
         ENDDO         
         CST = RIJRKJ/RKJ/RIJ
         IF (CST.GT. ONE) CST = 1.d0
         IF (CST.LT.-ONE) CST = -1.d0
C
C** DONE ANGLE COSINE CALCULATION
C
        IF (NFCTR .EQ. NFNDXL) THEN
C
C* if NFCTR.EQ.NFNDXL -> EBANG returned with value of NUMBND-th bond 
C  of the NMOL-th molecule
C
          EBANG = DACOS(CST)*CONV
        ELSE
C
C* if NFCTR.EQ.NFNDXE -> EBANG returned with energy of NUMBND-th bond
C  of the NMOL-th molecule
C
          IF (IDIANG.EQ.IANGH) THEN
C
C** OLD ANGLE/HARMONIC FUNCTIONAL FORM
C
            AT = DACOS(CST)
C we need it here               
            IF (LDOPER) THEN
              TREFA  = DACOS(TA0(NUMANG))
              TREFB  = DACOS(TB0(NUMANG))
C convert back new functional form force constant
              COST2  = TA0(NUMANG)*TA0(NUMANG)
              SDET   = COST2*COST2-(2.0d0*COST2-1.0d0)
     .             * (1.0d0-EKT/CTA(NUMANG))
              CSQ    = (COST2-DSQRT(SDET))/(2.0d0*COST2-1.0d0)
              CTTA   =  EKT/DACOS(CSQ)**2
              COST2  = TB0(NUMANG)*TB0(NUMANG)
              SDET   = COST2*COST2-(2.0d0*COST2-1.0d0)
     .             * (1.0d0-EKT/CTB(NUMANG))
              CSQ    = (COST2-DSQRT(SDET))/(2.0d0*COST2-1.0d0)
              CTTB   =  EKT/DACOS(CSQ)**2
C     
              DCB    = CTTB - CTTA
              DB0    = TREFB - TREFA
              DBA    = AT - TREFA
              DBL    = RLAM*DB0
              DB1    = DBA - DBL
              DB2    = DCB*DB1
              DB3    = DB0*CTTA
              DF     = (DB2-DB3)*RLAM
              DB5    = DB0*DB3*RLAM
              DB4    = DB1*DB2 + DB5 - 2.d0*DBA*DB3
              EBANG  = 0.5d0 * RLAM*DB4
              EBANGL = 0.5d0 * (DB4 + DB5 - 2.d0*DBL*DB2)
C     
            ELSE
              IC     = ICT(NUMANG)
              TREF   = DACOS(T0(IC))
              DA     = AT - TREF
              IF(DABS(CT(IC)).lt.EKT) THEN
                CTHA=0.d0
              ELSE
C     convert back new functional form force constant
                COST2 = T0(IC)*T0(IC)
                SDET  = COST2*COST2-(2.0d0*COST2-1.0d0)*
     .                  (1.0d0-EKT/CT(IC))
                CSQ   = (COST2-DSQRT(SDET))/(2.0d0*COST2-1.0d0)
                CTHA  =  EKT/DACOS(CSQ)**2
              ENDIF
C     
              EBANG  = 0.5d0*CTHA*DA*DA
            ENDIF
C     
          ELSE
C
C** NEW ANGLE COSINE/HARMONIC FUNCTIONAL FORM
C
            IF (LDOPER) THEN
              DCB    = CTB(NUMANG) - CTA(NUMANG)
              DB0    = TB0(NUMANG) - TA0(NUMANG)
              DBA    = CST - TA0(NUMANG)
              DBL    = RLAM*DB0
              DB1    = DBA - DBL
              DB2    = DCB*DB1
              DB3    = DB0*CTA(NUMANG)
              DF     = (DB2-DB3)*RLAM
              DB5    = DB0*DB3*RLAM
              DB4    = DB1*DB2+DB5-2.d0*DBA*DB3
              EBANG  = 0.5d0*RLAM*DB4
              EBANGL = 0.5d0*(DB4+DB5-2.D0*DBL*DB2)
            ELSE
              IC     = ICT(NUMANG)
              DA     = CST - T0(IC)
              DF     = CT(IC)*DA
              EBANG  = 0.5d0*DF*DA
            ENDIF
          ENDIF
        ENDIF
C
        RETURN
      ENDIF
C     
C******** DONE ONE ANGLE CALCULATION
C******** MULTIPLE ANGLE CALCULATION
C
      IF (NFCTR .EQ. NFALL) THEN
         LBNDS = .TRUE.
         LFORCE= .TRUE.
      ELSEIF (NFCTR .EQ. NFEF) THEN
         LBNDS = .FALSE.
         LFORCE = .TRUE.
      ELSEIF (NFCTR .EQ. NFEL) THEN
         LBNDS = .TRUE.
         LFORCE = .FALSE.
      ELSEIF (NFCTR .EQ. NFE) THEN
         LBNDS = .FALSE.
         LFORCE = .FALSE.
      ENDIF
C
      IATOFF = -1
      IBOFF = 0
C
      IF (IDIANG.EQ.IANGH) THEN
C
C**** OLD ANGLE/HARMONIC FUNCTIONAL FORM
C
        IF (.NOT.LDO34 ) THEN
C
C** CALCULATION IN 3D (OLD ANGLE/HARMONIC FUNCTIONAL FORM)
C        
          DO 100 MOL=1,NMOL
            DO 200 NNA=1,NUMANG

               I3 = NDIM*(IT(NNA)+IATOFF)
               J3 = NDIM*(JT(NNA)+IATOFF)
               K3 = NDIM*(KT(NNA)+IATOFF)
C
C** CALCULATE ANGLE COSINE, AND DERIVATIVE, CORRECT FOR BOUNDARY 
C   CONDITION IF REQUIRED
C
!$OMP parallel do private(M,XTMP)
               DO M = 1,NCALCD
                  XTMP = X(J3+M)
                  XXIJ(M) = X(I3+M) - XTMP
                  XXKJ(M) = X(K3+M) - XTMP
               ENDDO
               RIJ2 = 0.0d0
               RKJ2 = 0.0d0
               IF (LVAC) THEN
                  DO M=1,NCALCD
                     RIJ2 = RIJ2 + XXIJ(M)**2
                     RKJ2 = RKJ2 + XXKJ(M)**2
                  ENDDO
               ELSE
                  IF (LDOTRA) THEN
                     CALL TRACO(1,0,XXIJ,BETA,1,LEVERY)
                     CALL TRACO(1,0,XXKJ,BETA,1,LEVERY)
                  ENDIF
                  DO M=1,NCALCD
                     BOXM = BOX(M)
                     BOXHM= BOXH(M)
C     
                     DSTTMP = XXIJ(M)
                     IF (DSTTMP .GE. BOXHM) THEN
                        DSTTMP = DSTTMP - BOXM
                     ELSEIF (DSTTMP .LT. -BOXHM) THEN
                        DSTTMP = DSTTMP + BOXM
                     ENDIF
                     XXIJ(M) = DSTTMP
                     RIJ2 = RIJ2 + DSTTMP*DSTTMP
C     
                     DSTTMP = XXKJ(M)
                     IF (DSTTMP .GE. BOXHM) THEN
                        DSTTMP = DSTTMP - BOXM
                     ELSEIF (DSTTMP .LT. -BOXHM) THEN
                        DSTTMP = DSTTMP + BOXM
                     ENDIF
                     XXKJ(M) = DSTTMP
                     RKJ2 = RKJ2 + DSTTMP*DSTTMP
                  ENDDO
C     
                  IF (LOCTO) THEN
                     DSTTMP = BOXOQ - DABS(XXIJ(1))-
     $                  DABS(XXIJ(2))-DABS(XXIJ(3))
                     IF (DSTTMP .LT. 0.0d0) THEN
                        RIJ2 = RIJ2 + DSTTMP*BOX(1)
                        DO M=1,NCALCD
                          XXIJ(M)=XXIJ(M)-SIGN(BOXOH,XXIJ(M)) 
                        ENDDO
                     ENDIF
C
                     DSTTMP = BOXOQ - DABS(XXKJ(1))-
     $                  DABS(XXKJ(2))-DABS(XXKJ(3))
                     IF (DSTTMP .LT. 0.0d0) THEN
                        RKJ2 = RKJ2 + DSTTMP*BOX(1)
                        DO M=1,NCALCD
                          XXKJ(M)=XXKJ(M)-SIGN(BOXOH,XXKJ(M)) 
                        ENDDO
                     ENDIF
                  ELSEIF (LDOTRA) THEN
                     RIJ2=RIJ2+COSB2*XXIJ(1)*XXIJ(3)
                     RKJ2=RKJ2+COSB2*XXKJ(1)*XXKJ(3)
                  ENDIF
               ENDIF
C
               RKJ = DSQRT(RKJ2)
               RIJ = DSQRT(RIJ2)
               RIJRKJ = 0.0d0
               DO M=1, NCALCD
                  RIJRKJ = RIJRKJ + XXIJ(M)*XXKJ(M)
               ENDDO
C               
               CST = RIJRKJ/RKJ/RIJ
               IF (CST.GT. ONE) CST = 1.d0
               IF (CST.LT. -ONE) CST = -1.d0
C
C** DONE ANGLE COSINE AND DERIVATIVE CALCULATION
C   CALCULATE ENERGY/FORCES/LAMBDA DERIVATIVE
C
               AT = DACOS(CST)
               SNT= DSIN(AT)
               IF (DABS(SNT).LT.EPSSIN) SNT=1.D-12
C we need it here               
               IF (LDOPER) THEN
                 TREFA = DACOS(TA0(NNA))
                 TREFB = DACOS(TB0(NNA))
C convert back new functional form force constant
                 COST2 = TA0(NNA)*TA0(NNA)
                 SDET  = COST2*COST2-(2.0d0*COST2-1.0d0)
     .                 *(1.0d0-EKT/CTA(NNA))
                 CSQ   = (COST2-DSQRT(SDET))/(2.0d0*COST2-1.0d0)
                 CTTA  =  EKT/DACOS(CSQ)**2
                 COST2 = TB0(NNA)*TB0(NNA)
                 SDET  = COST2*COST2-(2.0d0*COST2-1.0d0)
     .                 *(1.0d0-EKT/CTB(NNA))
                 CSQ   = (COST2-DSQRT(SDET))/(2.0d0*COST2-1.0d0)
                 CTTB  =  EKT/DACOS(CSQ)**2
C
                 DCB    = CTTB - CTTA
                 DB0    = TREFB - TREFA
                 DBA    = AT - TREFA
                 DBL    = RLAM*DB0
                 DB1    = DBA - DBL
                 DB2    = DCB*DB1
                 DB3    = DB0*CTTA
                 DF     = (DB2-DB3)*RLAM
                 DB5    = DB0*DB3*RLAM
                 DB4    = DB1*DB2 + DB5 - 2.D0*DBA*DB3
                 EBAH   = 0.5d0*RLAM*DB4
                 EBANG  = EBANG + EBAH
                 EBGH   = 0.5d0*(DB4 + DB5 - 2.D0*DBL*DB2)
                 EBANGL = EBANGL + EBGH
                 DF = -DF/SNT
C                
               ELSE
C
                  IC = ICT(NNA)
                  TREF = DACOS(T0(IC))
                  DA = AT - TREF
                  IF(DABS(CT(IC)).lt.EKT) THEN
                    CTHA=0.d0 
                  ELSE
C convert back new functional form force constant
                    COST2 = T0(IC)*T0(IC)
                    SDET  = COST2*COST2-(2.0d0*COST2-1.0d0)
     .                    *(1.0d0-EKT/CT(IC))
                    CSQ   = (COST2-DSQRT(SDET))
     .                    /(2.0d0*COST2-1.0d0)
                    CTHA =  EKT/DACOS(CSQ)**2
                  ENDIF
                  DF = CTHA * DA
                  EBAH = 0.5D0*DF*DA
                  EBANG = EBANG + EBAH
                  DF = -DF/SNT
C
               ENDIF
C
               IF (LBNDS) THEN
                  ET0(IBOFF+NNA) = EBAH
                  XT0(IBOFF+NNA) = AT*CONV
               ENDIF
C
               IF (LFORCE) THEN
                  IF (LDOTRA) THEN
                     CALL TRACO(1,0,XXIJ,BETA,-1,LEVERY)
                     CALL TRACO(1,0,XXKJ,BETA,-1,LEVERY)
                  ENDIF
!$OMP parallel do private(M,DFI,DFK,DFJ)
                  DO M=1, NCALCD
                     DFI = DF*(XXKJ(M)/RKJ - XXIJ(M)/RIJ*CST)/RIJ
                     DFK = DF*(XXIJ(M)/RIJ - XXKJ(M)/RKJ*CST)/RKJ
                     DFJ = -DFI - DFK
                     F(I3+M) = F(I3+M) - DFI
                     F(J3+M) = F(J3+M) - DFJ
                     F(K3+M) = F(K3+M) - DFK
                  ENDDO
               ENDIF
 200        CONTINUE
C
            IATOFF = IATOFF + NAM
            IBOFF  = IBOFF  + NUMANG
C
 100     CONTINUE
C
C** DONE CALCULATION IN 3D (OLD ANGLE/HARMONIC FUNCTIONAL FORM)
C
         ELSE
           PRINT *,PRGSTR,'3D-4D not implemented IDIANG.EQ.IANGH'
           CALL STOPGM('COBOND','FATAL ERROR')
         ENDIF
C
C**** DONE OLD ANGLE/HARMONIC FUNCTIONAL FORM
C**** NEW ANGLE COSINE/HARMONIC FUNCTIONAL FORM
C         
      ELSEIF ( .NOT.LDO34 ) THEN
C
C** CALCULATION IN 3D (NEW ANGLE COSINE/HARMONIC FUNCTIONAL FORM)
C
         DO 1100 MOL=1,NMOL
            DO 1200 NNA=1,NUMANG
C
               I3 = NDIM*(IT(NNA)+IATOFF)
               J3 = NDIM*(JT(NNA)+IATOFF)
               K3 = NDIM*(KT(NNA)+IATOFF)
C
C** CALCULATE ANGLE COSINE, AND DERIVATIVE, CORRECT FOR BOUNDARY 
C   CONDITION IF REQUIRED
C
!$OMP parallel do private(M,XTMP)
               DO M = 1,NCALCD
                  XTMP = X(J3+M)
                  XXIJ(M) = X(I3+M) - XTMP
                  XXKJ(M) = X(K3+M) - XTMP
               ENDDO
               RIJ2 = 0.0D0
               RKJ2 = 0.0D0
               IF (LVAC) THEN
                  DO M=1,NCALCD
                     RIJ2 = RIJ2 + XXIJ(M)**2
                     RKJ2 = RKJ2 + XXKJ(M)**2
                  ENDDO
               ELSE
                  IF (LDOTRA) THEN
                     CALL TRACO(1,0,XXIJ,BETA,1,LEVERY)
                     CALL TRACO(1,0,XXKJ,BETA,1,LEVERY)
                  ENDIF
                  DO M=1,NCALCD
                     BOXM = BOX(M)
                     BOXHM= BOXH(M)
C     
                     DSTTMP = XXIJ(M)
                     IF (DSTTMP .GE. BOXHM) THEN
                        DSTTMP = DSTTMP - BOXM
                     ELSEIF (DSTTMP .LT. -BOXHM) THEN
                        DSTTMP = DSTTMP + BOXM
                     ENDIF
                     XXIJ(M) = DSTTMP
                     RIJ2 = RIJ2 + DSTTMP*DSTTMP
C     
                     DSTTMP = XXKJ(M)
                     IF (DSTTMP .GE. BOXHM) THEN
                        DSTTMP = DSTTMP - BOXM
                     ELSEIF (DSTTMP .LT. -BOXHM) THEN
                        DSTTMP = DSTTMP + BOXM
                     ENDIF
                     XXKJ(M) = DSTTMP
                     RKJ2 = RKJ2 + DSTTMP*DSTTMP
                  ENDDO
C     
                  IF (LOCTO) THEN
                     DSTTMP = BOXOQ - DABS(XXIJ(1))-
     .                  DABS(XXIJ(2))-DABS(XXIJ(3))
                     IF (DSTTMP .LT. 0.0d0) THEN
                        RIJ2 = RIJ2 + DSTTMP*BOX(1)
                        DO M=1,NCALCD
                          XXIJ(M)=XXIJ(M)-SIGN(BOXOH,XXIJ(M)) 
                        ENDDO
                     ENDIF
C
                     DSTTMP = BOXOQ - DABS(XXKJ(1))-
     .                  DABS(XXKJ(2))-DABS(XXKJ(3))
                     IF (DSTTMP .LT. 0.0d0) THEN
                        RKJ2 = RKJ2 + DSTTMP*BOX(1)
                        DO M=1,NCALCD
                          XXKJ(M)=XXKJ(M)-SIGN(BOXOH,XXKJ(M)) 
                        ENDDO
                     ENDIF
                  ELSEIF (LDOTRA) THEN
                     RIJ2=RIJ2+COSB2*XXIJ(1)*XXIJ(3)
                     RKJ2=RKJ2+COSB2*XXKJ(1)*XXKJ(3)
                  ENDIF
               ENDIF
C
               RKJ = DSQRT(RKJ2)
               RIJ = DSQRT(RIJ2)
               RIJRKJ = 0.0d0
               DO M=1, NCALCD
                 RIJRKJ = RIJRKJ + XXIJ(M)*XXKJ(M)
               ENDDO
C               
               CST = RIJRKJ/RKJ/RIJ
               IF (CST.GT. ONE) CST = 1.D0
               IF (CST.LT. -ONE) CST = -1.D0
C
C** DONE ANGLE COSINE AND DERIVATIVE CALCULATION
C   CALCULATE ENERGY/FORCES/LAMBDA DERIVATIVE
C
               IF (LDOPER) THEN
                  DCB    = CTB(NNA) - CTA(NNA)
                  DB0    = TB0(NNA) - TA0(NNA)
                  DBA    = CST - TA0(NNA)
                  DBL    = RLAM*DB0
                  DB1    = DBA - DBL
                  DB2    = DCB*DB1
                  DB3    = DB0*CTA(NNA)
                  DF     = (DB2-DB3)*RLAM
                  DB5    = DB0*DB3*RLAM
                  DB4    = DB1*DB2 + DB5 - 2.D0*DBA*DB3
                  EBAH   = 0.5d0*RLAM*DB4
                  EBANG  = EBANG + EBAH
                  EBGH   = 0.5d0*(DB4 + DB5 - 2.D0*DBL*DB2)
                  EBANGL = EBANGL + EBGH
               ELSE
                  IC     = ICT(NNA)
                  DA     = CST - T0(IC)
                  DF     = CT(IC)*DA
                  EBAH   = 0.5d0*DF*DA
                  EBANG  = EBANG + EBAH
               ENDIF
C
               IF (LBNDS) THEN
                 ET0(IBOFF+NNA) = EBAH
                 XT0(IBOFF+NNA) = DACOS(CST)*CONV
               ENDIF
C
               IF (LFORCE) THEN
                  IF (LDOTRA) THEN
                     CALL TRACO(1,0,XXIJ,BETA,-1,LEVERY)
                     CALL TRACO(1,0,XXKJ,BETA,-1,LEVERY)
                  ENDIF
!$OMP parallel do private(M,DFI,DFK,DFJ) 
                  DO M=1, NCALCD
                     DFI = DF*(XXKJ(M)/RKJ - XXIJ(M)/RIJ*CST)/RIJ
                     DFK = DF*(XXIJ(M)/RIJ - XXKJ(M)/RKJ*CST)/RKJ
                     DFJ = -DFI - DFK
                     F(I3+M) = F(I3+M) - DFI
                     F(J3+M) = F(J3+M) - DFJ
                     F(K3+M) = F(K3+M) - DFK
                  ENDDO
                ENDIF
 1200        CONTINUE
C
             IATOFF = IATOFF + NAM
             IBOFF  = IBOFF  + NUMANG
 1100     CONTINUE
C         
C** DONE CALCULATION IN 3D (NEW ANGLE COSINE/HARMONIC FUNCTIONAL FORM)
C
      ELSE
C
C** CALCULATION WITH 3D-4D COUPLING (NEW ANGLE COSINE/HARMONIC FUNCTIONAL FORM)
C
         DO 101 MOL=1,NMOL
            DO 201 NNA=1,NUMANG
C
               I = IT(NNA)+IATOFF
               J = JT(NNA)+IATOFF
               K = KT(NNA)+IATOFF
C
               IF ( C4D(I+1) .GE. EPSSM .OR. C4D(J+1) .GE. EPSSM 
     $              .OR. C4D(K+1) .GE. EPSSM ) THEN                  
                  I3 = NDIM*I
                  J3 = NDIM*J
                  K3 = NDIM*K
C
C** CALCULATE ANGLE COSINE, AND DERIVATIVE, CORRECT FOR BOUNDARY 
C   CONDITION IF REQUIRED, IN 3D AND 4D
C
                  DO M = 1,4
                     XTMP = X(J3+M)
                     XXIJ(M) = X(I3+M) - XTMP
                     XXKJ(M) = X(K3+M) - XTMP
                  ENDDO
C
                  IF (LVAC) THEN
                     RIJ2 = XXIJ(1)**2+XXIJ(2)**2+XXIJ(3)**2
                     RKJ2 = XXKJ(1)**2+XXKJ(2)**2+XXKJ(3)**2
                     RIJ24D = RIJ2 + XXIJ(4)**2
                     RKJ24D = RKJ2 + XXKJ(4)**2
                  ELSE
C no monoclinic case
cmb                     DO 351 M=1,4
cmb                        BOXM = BOX(M)
cmb                        BOXHM= BOXH(M)
cmbC     
cmb                        DSTTMP = XXIJ(M)
cmb                        IF (DSTTMP .GE. BOXHM) THEN
cmb                           DSTTMP = DSTTMP - BOXM
cmb                        ELSEIF (DSTTMP .LT. -BOXHM) THEN
cmb                           DSTTMP = DSTTMP + BOXM
cmb                        ENDIF
cmb                        XXIJ(M) = DSTTMP
cmb                        IF ( M.LE.3 ) THEN
cmb                           RIJ2 = RIJ2 + DSTTMP**2
cmb                        ELSE
cmb                           RIJ24D = RIJ2 + DSTTMP**2
cmb                        ENDIF
cmbC     
cmb                        DSTTMP = XXKJ(M)
cmb                        IF (DSTTMP .GE. BOXHM) THEN
cmb                           DSTTMP = DSTTMP - BOXM
cmb                        ELSEIF (DSTTMP .LT. -BOXHM) THEN
cmb                           DSTTMP = DSTTMP + BOXM
cmb                        ENDIF
cmb                        XXKJ(M) = DSTTMP
cmb                        IF ( M.LE.3 ) THEN
cmb                           RKJ2 = RKJ2 + DSTTMP**2
cmb                        ELSE
cmb                           RKJ24D = RKJ2 + DSTTMP**2
cmb                        ENDIF
cmb 351                 CONTINUE
cmb
                     DO M=1,3
                        BOXM = BOX(M)
                        BOXHM= BOXH(M)
                        DSTTMP = XXIJ(M)
                        IF (DSTTMP .GE. BOXHM) THEN
                           DSTTMP = DSTTMP - BOXM
                        ELSEIF (DSTTMP .LT. -BOXHM) THEN
                           DSTTMP = DSTTMP + BOXM
                        ENDIF
                        XXIJ(M) = DSTTMP
                        RIJ2 = RIJ2 + DSTTMP*DSTTMP
                        DSTTMP = XXKJ(M)
                        IF (DSTTMP .GE. BOXHM) THEN
                           DSTTMP = DSTTMP - BOXM
                        ELSEIF (DSTTMP .LT. -BOXHM) THEN
                           DSTTMP = DSTTMP + BOXM
                        ENDIF
                        XXKJ(M) = DSTTMP
                        RKJ2 = RKJ2 + DSTTMP*DSTTMP
                     ENDDO
cmb                  M=4
                     BOXM = BOX(4)
                     BOXHM= BOXH(4)
                     DSTTMP = XXIJ(4)
                     IF (DSTTMP .GE. BOXHM) THEN
                        DSTTMP = DSTTMP - BOXM
                     ELSEIF (DSTTMP .LT. -BOXHM) THEN
                        DSTTMP = DSTTMP + BOXM
                     ENDIF
                     XXIJ(4) = DSTTMP
                     RIJ24D = RIJ2 + DSTTMP*DSTTMP
                     DSTTMP = XXKJ(4)
                     IF (DSTTMP .GE. BOXHM) THEN
                        DSTTMP = DSTTMP - BOXM
                     ELSEIF (DSTTMP .LT. -BOXHM) THEN
                        DSTTMP = DSTTMP + BOXM
                     ENDIF
                     XXKJ(4) = DSTTMP
                     RKJ24D = RKJ2 + DSTTMP*DSTTMP
cmb
                     IF (LOCTO) THEN
                        DSTTMP = BOXOQ -DABS(XXIJ(1))-
     .                    DABS(XXIJ(2))-DABS(XXIJ(3))
                        IF (DSTTMP .LT. 0.0d0) THEN
                           RIJ2 = RIJ2 + DSTTMP*BOX(1)
                           DO M=1,NCALCD
                             XXIJ(M)=XXIJ(M)-SIGN(BOXOH,XXIJ(M)) 
                           ENDDO
                        ENDIF
C
                        DSTTMP = BOXOQ - DABS(XXKJ(1))-
     .                     DABS(XXKJ(2))-DABS(XXKJ(3))
                        IF (DSTTMP .LT. 0.0d0) THEN
                           RKJ2 = RKJ2 + DSTTMP*BOX(1)
                           DO M=1,NCALCD
                             XXKJ(M)=XXKJ(M)-SIGN(BOXOH,XXKJ(M)) 
                           ENDDO
                        ENDIF
C no monoclinic case
                     ENDIF
                  ENDIF
C
                  RKJ = DSQRT(RKJ2)
                  RIJ = DSQRT(RIJ2)
                  RKJ4D = DSQRT(RKJ24D)
                  RIJ4D = DSQRT(RIJ24D)
C
                  RIJRKJ = XXIJ(1)*XXKJ(1)
     .                   + XXIJ(2)*XXKJ(2)
     .                   + XXIJ(3)*XXKJ(3)
                  RIJRK4 = RIJRKJ + XXIJ(4)*XXKJ(4)
C                  
                  CST = RIJRKJ/RKJ/RIJ
                  IF (CST.GT. ONE) CST = 1.d0
                  IF (CST.LT. -ONE) CST = -1.d0
C
                  CST4D = RIJRK4/RKJ4D/RIJ4D
                  IF (CST4D.GT. ONE) CST4D = 1.d0
                  IF (CST4D.LT. -ONE) CST4D = -1.d0
C
C** DONE ANGLE COSINE AND DERIVATIVE CALCULATION
C   CALCULATE ENERGY/FORCES/LAMBDA/MU DERIVATIVE
C
                  IF (LDOPER) THEN
                     DCB    = CTB(NNA) - CTA(NNA)
                     DB0    = TB0(NNA) - TA0(NNA)
                     DB3    = DB0*CTA(NNA)
                     DBL    = RLAM*DB0
C* 3D
                     DBA    = CST - TA0(NNA)
                     DB1    = DBA - DBL
                     DB2    = DCB*DB1
                     DB5    = DB0*DB3*RLAM
                     DB4    = DB1*DB2 + DB5 - 2.d0*DBA*DB3
C
                     EBAH3  = 0.5d0*RLAM*DB4
                     EBGLH3 = 0.5d0*(DB4 + DB5 - 2.d0*DBL*DB2)
                     DF3D   = -(DB2-DB3)*RLAM*RMUE
C* 4D
                     DBA    = CST4D - TA0(NNA)
                     DB1    = DBA - DBL
                     DB2    = DCB*DB1
                     DB5    = DB0*DB3*RLAM
                     DB4    = DB1*DB2 + DB5 - 2.d0*DBA*DB3
C
                     EBAH4  = 0.5d0*RLAM*DB4
                     EBGLH4 = 0.5d0*(DB4 + DB5 - 2.d0*DBL*DB2)
                     DF4D   = (DB2-DB3)*RLAM*RMUE
C* sum
                     EBANG  = EBANG  + RMUE * (EBAH4  - EBAH3)
                     EBANGL = EBANGL + RMUE * (EBGLH4 - EBGLH3)
                     EBANGM = EBANGM + EBAH4 - EBAH3
C
                  ELSE
C
                     IC     = ICT(NNA)
C 3D
                     DA     = CST - T0(IC)
                     DF3D   = CT(IC)*DA
                     EBAH3  = 0.5d0*DF3D*DA
                     DF3D   = -DF3D*RMUE
C 4D
                     DA     = CST4D - T0(IC)
                     DF4D   = CT(IC)*DA
                     EBAH4  = 0.5d0*DF4D*DA
                     DF4D   = DF4D*RMUE
C sum
                     EBANG  = EBANG  + RMUE * (EBAH4 - EBAH3)
                     EBANGM = EBANGM + EBAH4 - EBAH3
                  ENDIF
C
                  IF (LBNDS) THEN
                     AT = DACOS(RMUE*CST4D + (1.0d0-RMUE)*CST4D)
                     ET0(IBOFF+NNA) = RMUE * (EBAH4 - EBAH3)
                     XT0(IBOFF+NNA) = AT*CONV
                  ENDIF
C
                  IF (LFORCE) THEN
C no monoclinic case
                     DO M=1, 3
                        DFI = DF3D*(XXKJ(M)/RKJ - 
     $                       XXIJ(M)/RIJ*CST)/RIJ
                        DFK = DF3D*(XXIJ(M)/RIJ - 
     $                       XXKJ(M)/RKJ*CST)/RKJ
                        DFJ = -DFI - DFK
                        F(I3+M) = F(I3+M) - DFI
                        F(J3+M) = F(J3+M) - DFJ
                        F(K3+M) = F(K3+M) - DFK
                     ENDDO
                     DO M=1, 4
                        DFI = DF4D*(XXKJ(M)/RKJ4D - 
     $                       XXIJ(M)/RIJ4D*CST4D)/RIJ4D
                        DFK = DF4D*(XXIJ(M)/RIJ4D - 
     $                       XXKJ(M)/RKJ4D*CST4D)/RKJ4D
                        DFJ = -DFI - DFK
                        F(I3+M) = F(I3+M) - DFI
                        F(J3+M) = F(J3+M) - DFJ
                        F(K3+M) = F(K3+M) - DFK
                     ENDDO
                  ENDIF
               ENDIF
 201        CONTINUE
C
            IATOFF = IATOFF + NAM
            IBOFF  = IBOFF  + NUMANG
C
 101     CONTINUE
C         
C
C** DONE CALCULATION WITH 3D-4D COUPLING (NEW ANGLE COSINE/HARMONIC 
C   FUNCTIONAL FORM)
C**** DONE NEW ANGLE COSINE/HARMONIC FUNCTIONAL FORM
C
      ENDIF
C
C******** DONE MULTIPLE ANGLE CALCULATION
C
C end angle
 500  FORMAT (1X,A6,': IMPLEMENTATION RESTRICTION')
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
COMMSUBR DIHANG
C W.F. VAN GUNSTEREN, GRONINGEN, JAN. 1987
C The GROMOS96 functional form versions supplied by Thomas Beutler.
C Revised/documented/GROMOS87 compatibility P.Hunenberger, sept. 95
C
C      SUBROUTINE DIHANG(NCALCD,NMOL,NAM,NFCTR,IDIDA,LEVERY,
C     $     NUMDA,IP,JP,KP,LP,ICP,
C     $     NDATY,CP,PD,NP,
C     $     NATTOT,X,F,EDA,XP0,EP0,
C     $     LDOPER,RLAM,NPDA,CPA,CPB,PDA,PDB,NPA,NPB,EDAGL,
C     $     LDO34,RMUE,C4D,EDAGM)
C
C     DIHANG will supply the dihedral angle (P), energies and forces
C     according to a potential energy function depending on IDIDA
C     for dihedral in NCALCD dimensions (NCALCD.LE.NDIM).
C     when NTB.NE.NTBVAC (with some restrictions mentionned below
C     for the 4D case)
C     More than one identical molecule can be considered.
C
C     The cobond.h header file defines reference values for the
C     NFCTR and IDIDA switches.
C     DIHANG assumes that NDIM, NTB, BOX(3), BOXH(3) and BETA are 
C     defined and transmitted through the common blocks of the box.h 
C     header file.
C     The file pertsz.h (perturbation) is optional (LDOPER.EQ.TRUE.)
C
C     When LDOPER. is .FALSE. the standard potential,
C     V_stand(P;PD,NP,CP) = sum over all specified dihedral angles of
C     1.         CP * (1+PD*COS(NP*P))  if IDIDA=IDICOS (NCALCD = 3 or 4)
C     2.        0.5 * CP * (P-PD)**2    if IDIDA=IDINCS (NCALCD = 3 only )
C     is evaluated, where PD is the cosine of the phase shift (IDIDA.EQ.
C     IDICOS) or the dihedral at minimal energy (IDIDA.EQ.IDINCS), 
C     CP the force constant, NP the multiplicity (IDIDA.EQ.IDICOS only) and
C     the difference P-PD in 2. is taken MODULO 2*PI (360 degrees).
C
C     In GROMOS, the 1. functional form above is typically used for
C     dihedrals and the 2. form is used for improper dihedrals.
C
C     The dihedral angle (P) is defined according to the IUPAC-IUB
C     convention (BIOCHEMISTRY 9 (1970) 3471) when IDIDA = IDICOS or
C     IDIDA = IDINCS.
C
C     When LDOPER is .TRUE., a LAMBDA (RLAM) dependent perturbation 
C     potential,
C     V_pert(P;RLAM;PD,CP,NP) = sum over all specified dihedral angles of
C     1.       RLAM * (CPB*(1+PDB*COS(NPB*P)) - CPA*(1+PDA*COS(NPA*P)))
C              if IDIDA=IDICOS
C     2.       0.5 * [ (1-RLAM)*CPA+RLAM*CPB ] *
C                    [ P - ((1-RLAM)*PDA+RLAM*PDB) ]**2
C             -0.5 * CPA * (P-PDA)**2
C              if IDIDA=IDINCS
C     is evaluated, where PDA/PDB denote the cosines of the phase shifts 
C     (IDIDA.EQ.IDICOS) or the dihedral at minimal energy (IDIDA.EQ.IDINCS),
C     CPA/CPB the force constants and NPA/NPB the multiplicities (IDIDA.EQ.
C     only), in states A/B. The difference P-PD in 2. is taken MODULO 2*PI 
C     (360 degrees).
C     The second term in 2. (already substracted in 1.)is intented to 
C     cancel the V_pert(RLAM=0) contribution already taken into account 
C     by a previous call to dihangle (all dihedrals, LDOPER=.FALSE.) 
C     i.e. to avoid double counting. 
C     dV_pert/d(RLAM) (EBGL) is also evaluated.
C     It is not implemented for NFCTR.EQ.NFNDXL.
C
C     If LDO34 is .TRUE., a MU (RMUE) dependent 3D-4D perturbation 
C     potential,
C     V_stand34(P;RMUE;CP,NP) = sum over all specified angles of
C             + RMUE * V(P_4D;CP,PD,NP) - RMUE * V(P_3D;CP,PD,NP) 
C             if LDOPER.EQ..FALSE.
C     or
C     V_pert34(T;RMUE;RLAM;CPA,CPB,PDA,PDB,NPA,NPB) = sum over all 
C     specified dihedral angles of
C             + RMUE * V_pert(P_4D;RLAM;CPA,CPB,PDA,PDB,NPA,NPB)
C             - RMUE * V_pert(P_3D;RLAM;CPA,CPB,PDA,PDB,NPA,NPB)
C     is evaluated.
C
C     This is only allowed when IDIDA.EQ.IDICOS
C     Note that the contribution of the 3D component is -RMUE and 
C     not (1-RMUE) to avoid double counting of the V_34(RMUE=0) component
C     already taken into account by previous call(s) to dihang.
C     dV_pert/d(RLAM) (EBGL) and dV_pert/d(RMUE) (EBGM) are also evaluated.
C     LDO34 requires NDIM.EQ.4 .AND. NCALCD.EQ.4.
C     It is not implemented for NFCTR.EQ.NFNDXL.OR.NFCTR.EQ.NFNDXE.
C     Only a rectangular or truncated-octahedron box (BETA.EQ.90) can be 
C     used. In the latter case, only rectangular periodic gathering is
C     applied in the NCALCD-3 last dimensions.
C     Dihedral angle values are not returned in this case.
C
C     Arguments:
C
C     NCALCD      the number of dimensions in which to consider
C                 all vectors, coordinate/forces (NCALCD.LE.NDIM).
C     NMOL        The meaning of this argument depends on NFCTR: it
C                 is either the number of identical molecules with 
C                 NUMDA dihedrals or the sequence number of the molecule 
C                 for which dihedral number NUMDA is specified.
C                 See the table for NFCTR below.
C     NAM         number of atoms per molecule (NAM.GT.3), only used
C                 if NMOL.GT.1
C     NFCTR       determines the actions performed.
C                 NFCTR also influences the meanings of NUMDA and NMOL.
C                 (NFCTR.GE.NFCDMN.AND.NFCTR.LE.NFCMAX)
C
C          NUMDA  means   NMOL means       vars delivered
C   ------------------------------------------------------------------------
C   L<NFALL>  num of dihedrals  num of mol.      EDA,F,EP0,XP0
C   L<NFEF>   num of dihedrals  num of mol.      EDA,F
C   L<NFEL>   num of dihedrals  num of mol.      EDA,EP0,XP0
C   L<NFE>    num of dihedrals  num of mol.      EDA
C   L<NFNDXL> seq num           seq num of mol.  EDA(with angle of the
C                                             NUMDA-th dihedral)
C   L<NFNDXE> seq num           seq num of mol.  EDA(with energy of the
C                                             NUMDA-th dihedral)
C   L<NFDER>  seq num           seq num of mol.  EDA(with angle of the
C                                             NUMDA-th dihedral), 
C                                             F with derivative
C   ------------------------------------------------------------------------
C
C     IDIDA       Functional form selected for the potential, see above.
C                 (IDIDA .GT. L<IDIMIN> .AND. IDIDA .LT. L<IDIMAX>)
C     LEVERY
C         .FALSE.    quantities depending on BETA are only calculated
C                    at the first subroutine call.
C         .TRUE.     they are calculated at every subroutine call.
C     NUMDA       the meaning of NUMDA depends on NFCTR:
C                 it is either the number of dihedral angles
C                 per molecule or the dihedral angle sequence number.
C                 See the table for NFCTR above. 
C                 Note that if LDOPER.EQ.TRUE, NUMDA should be equal
C                 to NPDA
C     IP,JP,KP,LP(NUMDA) atom topology indices of the atoms forming the 
C                 dihedrals.
C                 ( accessed whatever LDOPER )
C     IPB(NUMDA)  dihedral-type code denoting the index
C                 in the arrays CP,PD and NP for dihedral I-J-K-L. Only used 
C                 if LDOPER.EQ..FALSE. 
C                 Otherwise the angle parameters are taken sequentially 
C                 from the arrays CPA,CPB,PDA,PDB,NPA,NPB(NPDA).
C     NDATY       number of dihedral angle types
C     CP(NDATY)   force constants (only accessed if LDOPER.EQ..FALSE.)
C     PD(NDATY)   cosines of the phase shifts (IDIDA.EQ.IDICOS) or the 
C                 dihedral at minimal energy (IDIDA.EQ.IDINCS).
C                 (only accessed if LDOPER.EQ..FALSE.)
C     NP(NDATY)   multiplicities, restricted to 
C                 NP(I)=0,1,2,3,4,5 OR 6 
C                 (IDIDA.EQ.IDICOS, ignored otherwise)
C     NATTOT      the total number of atoms to be considered in the 
C                 X and F arrays.
C     X(NDIM*NATTOT)   atom cartesian coordinates.
C     F(NDIM*NATTOT)   delivered with the dihedral angle forces added if
C                 NFCTR.EQ.NFALL.OR.NFCTR.EQ.NFEFC, or 
C                 the dihedral angle derivative if NFCTR.EQ.NFDER.
C                 In the latter case NCALCD = 3 and IDIDA = IDINCS 
C                 are required.
C                 The twelve values are written into, not added into,
C                 the _beginning_ of the force array, regardless of the
C                 IP,JP,KP and LP and of NMOL and NAM.
C                 Thus, if NDIM = 3, the first 12 elements
C                 of F are filled with values. If NDIM = 4, the first
C                 16 elements are used, with the 4,8,12 and 16th elements
C                 set to zero.
C     EDA         what is delivered in this variable depends on
C                 NFCTR according to the table below.
C
C   NFCTR                    EDA delivered with
C   ------------------------------------------------------------------------
C   NFALL,NFEF,NFEL,NFE      total dihedral energy
C   NFNDXL                   value of the NUMDA-th
C                            dihedral of the NMOL-th molecule
C   NFNDXE                   energy of the NUMDA-th
C                            dihedral of the NMOL-th molecule
C   ------------------------------------------------------------------------
C
C
C     XP0(NUMDA*NMOL)  delivered with the dihedral angle values if
C                 NFCTR.EQ.NFEL.OR.NFCTR.EQ.NFALL.
C                 ( LDO34.EQ..FALSE. )
C     EP0(NUMDA*NMOL)  delivered with the dihedral angle energies.
C                 NFCTR.EQ.NFEL.OR.NFCTR.EQ.NFALL
C     LDOPER      determines whether a lamda dependent potential
C                 is calculated or not.
C            LDOPER .EQ. .FALSE.
C                 the potential V_stand is evaluated using IPB sequentially
C                 and the corresponding CP,PD,NP.
C            LDOPER .EQ. .TRUE.
C                 LAMDA dependent potential V_pert is evaluated
C                 using CPA,CPB,PDA,PDB,NPA,NPD sequentially.
C     RLAM        StateA - StateB Hamiltonian coupling parameter lambda.
C                 (only accessed if LDOPER.EQ.TRUE)
C     NPDA        number of perturbed dihedral angles.
C                 (only accessed if LDOPER.EQ..TRUE.)
C     CPA,CPB(NPDA) force constants state A and B
C                 (only accessed if LDOPER.EQ..TRUE.)
C     PDA,PDB(NPDA) cosines of the phase shifts (IDIDA.EQ.IDICOS) or the 
C                 dihedrals at minimal energy (IDIDA.EQ.IDINCS),state A and B
C                 (only accessed if LDOPER.EQ..FALSE.)
C     NPA,NPB(NPDA) multiplicities, restricted to 
C                 NP(I)=0,1,2,3,4,5 OR 6 
C                 (IDIDA.EQ.IDICOS, ignored otherwise), state A and B
C     EBGL        delivered with the value of the derivative
C                 V_pert (or V_pert34) with respect to RLAM, 
C                 if LDOPER.EQ..TRUE.
C     LDO34       determines whether a mu dependent potential
C                 is calculated or not.
C            LDO34 .EQ. .FALSE.
C                 the potential V_stand or V_pert are evaluated.
C            LDO34 .EQ. .TRUE.
C                 the potential V_stand34 or V_pert34 are evaluated.
C     RMUE        3D - 4D Hamiltonian coupling parameter mu.
C                 (only accessed if LDO34.EQ..TRUE.)
C     C4D(NATTOT) force constant for the 4th dimentional harmonic
C                 spring. Only used here as C4D(I).GT.EPSSM (EPSSM=-1E-20)
C                 implies that the atom is in 4D. A dihedral is treated in 4D
C                 if any of the atoms forming the dihedral is in 4D.
C                 (only accessed if LDO34.EQ..TRUE.).
C     EBGM        delivered with the value of the derivative
C                 V_stand34 or V_pert34 with respect to RMUE, if
C                 LDO34.EQ..TRUE.
COMMEND
C     ==--------------------------------------------------------------==
      SUBROUTINE DIHANG(NCALCD,NMOL,NAM,NFCTR,IDIDA,LEVERY,
     $     NUMDA,IP,JP,KP,LP,ICP,
     $     NDATY,CP,PD,NP,
     $     NATTOT,X,F,EDA,XP0,EP0,
     $     LDOPER,RLAM,NPDA,CPA,CPB,PDA,PDB,NPA,NPB,EDAGL,
     $     LDO34,RMUE,C4D,EDAGM)
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'pertsz.h'
      INCLUDE 'cobond.h'
C args
      LOGICAL LDOPER,LEVERY
Csco
C args
      INTEGER NUMDA,NCALCD,NMOL,NAM,NDATY,NFCTR,IDIDA,NATTOT
      INTEGER IP(NUMDA),JP(NUMDA),KP(NUMDA),LP(NUMDA),ICP(NUMDA)
      INTEGER NP(NDATY)
      real*8 CP(NDATY),PD(NDATY),EDA,EDAGL
      real*8 X(NDIM*NATTOT),F(NDIM*NATTOT)
      real*8 XP0(NUMDA*NMOL),EP0(NUMDA*NMOL)
C perturbation
      INTEGER NPDA
      INTEGER NPA(NPDA),NPB(NPDA)
      real*8 CPA(NPDA),PDA(NPDA),CPB(NPDA),PDB(NPDA)
      real*8 RLAM
C 3D-4D coupling staff
      LOGICAL LDO34
      real*8 C4D(NATTOT),RMUE,EDAGM

C local params
      real*8 ONE
      PARAMETER (ONE = 1.0d0)
      real*8 EPS
      PARAMETER (EPS = 1.D-6)
C local vars
      LOGICAL LFIRST
      SAVE LFIRST

      LOGICAL LMONO,LOCTO,LVAC,LDOVIR,LDOTRA
      SAVE LMONO,LOCTO,LVAC,LDOVIR,LDOTRA

      real*8 PYE,CONV,BETAR,COSB,COSB2,BOXOH,BOXOQ,PYETWO,PYETWI
      SAVE COSB,COSB2,BOXOH,BOXOQ,CONV,PYETWO,PYETWI
C
      INTEGER I3,J3,K3,L3,IATOFF,IDAOFF,ICOFF
      INTEGER II3,JJ3,KK3,LL3
      INTEGER M,IC,NDBA,NARG,MOL,NNA
      INTEGER I,J,K,L
      real*8 XJ,XK,XXIJ(MAXDIM),XXKJ(MAXDIM),XXKL(MAXDIM)
      real*8 XXMJ(MAXDIM),XXNK(MAXDIM),XXIM(MAXDIM),XXLN(MAXDIM)
      real*8 DPDRI(MAXDIM),DPDRJ(MAXDIM),DPDRK(MAXDIM),DPDRL(MAXDIM)
      real*8 DCPDRI(MAXDIM),DCPDRL(MAXDIM)
      real*8 DCPDRJ(MAXDIM),DCPDRK(MAXDIM)
      real*8 RIJ2,RKJ2,RKL2,RMJ2,RNK2,RIM2,RLN2
      real*8 RMJ,RNK,RIM,RLN,RKJ
      real*8 RRMJNK,RRIMLN
      real*8 RIJRNK,RMJRNK,RIJRKJ,RKLRKJ,RIMRLN
      real*8 BOXM,BOXHM,DSTTMP
      real*8 AP,CSP,EDAH,EDAD
      real*8 CSMP,CSMPA,CSMPB,DCMPDC,DCMPCA,DCMPCB
      real*8 SP,SPA,SPB,DCB,DB,DB0,DBA,DB1,DB2,DB3,DB4,DB5,DBL
      real*8 DUMA, DUMB
      real*8 ARG,EDGH
      real*8 CSP4D,SP3D,SP4D
      real*8 EDAH3,EDAH4
      real*8 EDGLH3,EDGLH4
      real*8 RIJ24D,RKL24D
      real*8 RIM4D,RLN4D,RIM24D,RLN24D
      real*8 RIJRK4,RKLRK4,RKJ24D,RIMRL4,RRIML4
      real*8 XXLN4(MAXDIM),XXIM4(MAXDIM)
      real*8 EPSSM
      LOGICAL LBNDS,LFORCE
      CHARACTER PRGSTR*(6)
C data
      DATA EPSSM /-1.0D-20/
      DATA LFIRST/.TRUE./
      DATA PRGSTR /'DIHANG'/
C begin
C******** CHECK ARGUMENTS FOR CONSISTENCY
C
      IF (IDIDA.LT.IDIMIN .OR. IDIDA.GT.IDIMAX ) THEN
         PRINT *,PRGSTR,': invalid choice of IDIDA'
         CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF (NFCTR.LT.NFCDMN .OR. NFCTR.GT.NFCMAX) THEN
        PRINT *,PRGSTR,': invalid choice of NFCTR'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( NAM.LE.3 .AND. NMOL.GT.1 ) THEN
        PRINT *,PRGSTR,': NMOL.GT.1 requires NAM.GT.3'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF (LDOPER .AND. NUMDA.NE.NPDA) THEN
        PRINT *,PRGSTR,': LDOPER requires NUMDA.EQ.NPDA'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF (NCALCD.LT.3.OR.NCALCD.GT.4) THEN
         PRINT 500,PRGSTR
         PRINT *,'DIHANG is defined only for NCALCD .EQ.3 or 4' 
         CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF (NCALCD .NE. 3 .AND. IDIDA .EQ. IDINCS) THEN
         PRINT *,PRGSTR,': IDIDA.EQ.IDINCS requires NCALCD.EQ.3'        
         CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF (LDO34 .AND. ( NDIM .NE. 4 .OR. NCALCD .NE. 4 ) ) THEN
         PRINT 500,PRGSTR
         PRINT *,'3D - 4D only possible if NDIM .EQ. 4 .AND. ',
     $           'NCALCD .EQ. 4'
         CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF (LDO34 .AND. NTB.NE.NTBVAC .AND. BETA .NE. 90.0D0 ) THEN
         PRINT *,PRGSTR,'3D - 4D only possible if  BETA .EQ. 90.0 '
         CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( (NFCTR .EQ. NFNDXL .OR. NFCTR .EQ. NFNDXE) .AND.
     $     NCALCD.NE.3 ) THEN
        PRINT 500,PRGSTR
        PRINT *,'NFCTR .EQ. NFNDXL or NFCTR .EQ. NFNDXE ',
     $       'require NCALCD .EQ. 3'        
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( (NFCTR.EQ.NFEL.OR.NFCTR.EQ.NFALL.OR.NFCTR.EQ.NFNDXL)
     $     .AND. NCALCD.NE.3 ) THEN
        PRINT 500,PRGSTR
        PRINT *,'NFCTR .EQ. NFEL or NFCTR .EQ. NFALL ',
     $       'NFCTR .EQ. NFNDXL require NCALCD .EQ. 3'
        PRINT *,'a dihedral angle cannot be calculated in 4D'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( LDO34.AND.NFCTR.EQ.NFNDXE ) THEN
        PRINT 500,PRGSTR
        PRINT *,'3D - 4D incompatible with NFCTR .EQ. NFNDXE'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF ( (LDOPER.OR.LDO34).AND.NFCTR.EQ.NFNDXL ) THEN
        PRINT 500,PRGSTR
        PRINT *,'3D - 4D and LDOPER incompatible with NFCTR',
     $          ' .EQ. NFNDXL'
        CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
      IF (NFCTR .EQ. NFDER .AND. IDIDA .NE. IDINCS) THEN
         PRINT *,PRGSTR,'invalid value combination NFCTR.EQ.NFDER ',
     $       'requires IDIDA = IDINCS'
         CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF
C
C******** DONE CHECK ARGUMENTS FOR CONSISTENCY
C******** FIRST CALL SETTINGS
C          
      IF (LEVERY .OR. LFIRST) THEN
C
         PYE = 4.D0*DATAN(ONE)
         PYETWO = 2.D0*PYE
         PYETWI = ONE/PYETWO
         CONV = 1.8D2/PYE
C
         LFIRST = .FALSE.
C
         LMONO = (NTB .GT. 0)
         LOCTO = (NTB. LT. 0)
         LVAC  = (NTB .EQ. 0)
         LDOVIR= (ABS(NTB) .EQ. NTBVIR)
C
         LDOTRA = .FALSE.
         IF (LMONO)THEN
            BETAR = BETA/CONV
            COSB = DCOS(BETAR)
            LDOTRA = (DABS(COSB).GE.1.D-4)
            COSB2 = 2.0d0*COSB
         ENDIF
C
         IF (LOCTO) THEN
            BOXOH = BOXH(1)
            BOXOQ = BOX(1)*0.75D0
         ELSE
            BOXOH = 0.0D0
            BOXOQ = 0.0D0
         ENDIF
      ENDIF
C
C******** DONE FIRST CALL SETTINGS
C
      EDA = 0.0D0
      EDAGL = 0.0D0
      EDAGM = 0.0D0
C
C******** ONE DIHEDRAL CALCULATION
C
      IF (NFCTR .EQ. NFNDXL .OR. NFCTR .EQ. NFNDXE 
     $     .OR. NFCTR .EQ. NFDER) THEN
C
C* NFCTR.EQ.NFNDXL -> EDA returned with value of NUMDA-th dihedral
C  of the NMOL-th molecule
C  NFCTR.EQ.NFDER  -> idem, F() delivered with the angle derivative
C  (requires IDIDA.EQ.IDINCS)
C  NFCTR.EQ.NFNDXE -> EDA returned with energy of NUMDA-th dihedral
C  of the NMOL-th molecule
C
         IATOFF = (NMOL-1)*NAM
         ICOFF  = NDIM*(IATOFF-1)
         I3 = NDIM*IP(NUMDA) + ICOFF
         J3 = NDIM*JP(NUMDA) + ICOFF
         K3 = NDIM*KP(NUMDA) + ICOFF
         L3 = NDIM*LP(NUMDA) + ICOFF
C
C** CALCULATE DIHEDRAL, AND DERIVATIVE, CORRECT FOR BOUNDARY CONDITION 
C   IF REQUIRED
C
         DO M=1,NCALCD
            XJ = X(J3+M)
            XK = X(K3+M)
            XXIJ(M) = X(I3+M) - XJ
            XXKJ(M) = XK - XJ
            XXKL(M) = XK - X(L3+M)
         ENDDO
C
         RIJ2 = 0.0D0
         RKJ2 = 0.0D0
         RKL2 = 0.0D0
         IF (LVAC) THEN
            DO M=1,NCALCD
               RIJ2 = RIJ2 + XXIJ(M)**2
               RKJ2 = RKJ2 + XXKJ(M)**2
               RKL2 = RKL2 + XXKL(M)**2
            ENDDO
         ELSE
            IF (LDOTRA) THEN
               CALL TRACO(1,0,XXIJ,BETA,1,LEVERY)
               CALL TRACO(1,0,XXKJ,BETA,1,LEVERY)
               CALL TRACO(1,0,XXKL,BETA,1,LEVERY)
            ENDIF
            DO 150 M=1,NCALCD
               BOXM = BOX(M)
               BOXHM= BOXH(M)
C
               DSTTMP = XXIJ(M)
               IF (DSTTMP .GE. BOXHM) THEN
                  DSTTMP = DSTTMP - BOXM
               ELSEIF (DSTTMP .LT. -BOXHM) THEN
                  DSTTMP = DSTTMP + BOXM
               ENDIF
               XXIJ(M) = DSTTMP
               RIJ2 = RIJ2 + DSTTMP*DSTTMP
C
               DSTTMP = XXKJ(M)
               IF (DSTTMP .GE. BOXHM) THEN
                  DSTTMP = DSTTMP - BOXM
               ELSEIF (DSTTMP .LT. -BOXHM) THEN
                  DSTTMP = DSTTMP + BOXM
               ENDIF
               XXKJ(M) = DSTTMP
               RKJ2 = RKJ2 + DSTTMP*DSTTMP
C
               DSTTMP = XXKL(M)
               IF (DSTTMP .GE. BOXHM) THEN
                  DSTTMP = DSTTMP - BOXM
               ELSEIF (DSTTMP .LT. -BOXHM) THEN
                  DSTTMP = DSTTMP + BOXM
               ENDIF
               XXKL(M) = DSTTMP
               RKL2 = RKL2 + DSTTMP*DSTTMP
 150        CONTINUE
C
            IF (LOCTO) THEN
               DSTTMP = BOXOQ - DABS(XXIJ(1))-
     .            DABS(XXIJ(2))-DABS(XXIJ(3))
               IF (DSTTMP .LT. 0.0d0) THEN
                  RIJ2 = RIJ2 + DSTTMP*BOX(1)
                  DO M=1,NCALCD
                    XXIJ(M)=XXIJ(M)-SIGN(BOXOH,XXIJ(M)) 
                  ENDDO
               ENDIF
C
               DSTTMP = BOXOQ - DABS(XXKJ(1))-
     .            DABS(XXKJ(2))-DABS(XXKJ(3))
               IF (DSTTMP .LT. 0.0d0) THEN
                  RKJ2 = RKJ2 + DSTTMP*BOX(1)
                  DO M=1,NCALCD
                    XXKJ(M)=XXKJ(M)-SIGN(BOXOH,XXKJ(M)) 
                  ENDDO
               ENDIF
C
               DSTTMP = BOXOQ - DABS(XXKL(1))-
     .            DABS(XXKL(2))-DABS(XXKL(3))
               IF (DSTTMP .LT. 0.0d0) THEN
                  RKL2 = RKL2 + DSTTMP*BOX(1)
                  DO M=1,NCALCD
                    XXKL(M)=XXKL(M)-SIGN(BOXOH,XXKL(M)) 
                  ENDDO
               ENDIF
            ELSEIF (LDOTRA) THEN
               RIJ2 = RIJ2 + COSB2*XXIJ(1)*XXIJ(3)
               RKJ2 = RKJ2 + COSB2*XXKJ(1)*XXKJ(3)
               RKL2 = RKL2 + COSB2*XXKL(1)*XXKL(3)
            ENDIF
         ENDIF
C
C** DONE DIHEDRAL CALCULATION
C
         IF (IDIDA .EQ. IDINCS .OR. NFCTR .EQ. NFNDXL ) THEN
C
C* if NFCTR.EQ.NFNDXL -> EDA returned with value of NUMBND-th bond 
C  of the NMOL-th molecule (calculated here whatever IDIDA)
C  if NFCTR.EQ.NFDER  -> idem, F() delivered with the angle derivative
C  (requires IDIDA.EQ.IDINCS)
C  NFCTR.EQ.NFNDXE -> EDA returned with energy of NUMDA-th dihedral
C  of the NMOL-th molecule (calculated here for IDIDA.EQ.IDINCS)
C
            XXMJ(1) = XXIJ(2)*XXKJ(3) - XXIJ(3)*XXKJ(2)
            XXMJ(2) = XXIJ(3)*XXKJ(1) - XXIJ(1)*XXKJ(3)
            XXMJ(3) = XXIJ(1)*XXKJ(2) - XXIJ(2)*XXKJ(1)
            XXNK(1) = XXKJ(2)*XXKL(3) - XXKJ(3)*XXKL(2)
            XXNK(2) = XXKJ(3)*XXKL(1) - XXKJ(1)*XXKL(3)
            XXNK(3) = XXKJ(1)*XXKL(2) - XXKJ(2)*XXKL(1)
C
            RMJ2 = 0.0D0
            RNK2 = 0.0D0
            DO M=1, NCALCD
               RMJ2 = RMJ2 + XXMJ(M)**2
               RNK2 = RNK2 + XXNK(M)**2
            ENDDO
            RRMJNK = DSQRT(RMJ2*RNK2)
            IF (RRMJNK .LT. EPS) THEN
C return zero if angle is not defined
               EDA = 0.0D0
               EDAGL = 0.0D0
            ELSE
C angle is defined              
               RIJRNK = 0.0D0
               RMJRNK = 0.0D0
               DO M=1, NCALCD
                  RIJRNK = RIJRNK + XXIJ(M)*XXNK(M)
                  RMJRNK = RMJRNK + XXMJ(M)*XXNK(M)
               ENDDO
C
               CSP = RMJRNK/RRMJNK
               IF (CSP.GT.1.D0) THEN
                  CSP = 1.D0
               ELSEIF (CSP.LT.-1.D0) THEN
                  CSP = -1.D0
               ENDIF
C
               IF (NFCTR .EQ. NFNDXL .OR. NFCTR .EQ. NFDER) THEN
C
C* if NFCTR.EQ.NFNDXL -> EDA returned with value of NUMBND-th bond 
C  of the NMOL-th molecule (calculated here whatever IDIDA)
C  if NFCTR.EQ.NFDER  -> idem, F() delivered with the angle derivative
C  (requires IDIDA.EQ.IDINCS)
C
                  DB1 = XXIJ(1)*(XXKJ(2)*XXKL(3)-XXKJ(3)*XXKL(2))
                  DB2 = XXIJ(2)*(XXKJ(3)*XXKL(1)-XXKJ(1)*XXKL(3))
                  DB3 = XXIJ(3)*(XXKJ(1)*XXKL(2)-XXKJ(2)*XXKL(1))
                  DB = DB1 + DB2 + DB3               
                  AP = SIGN(DACOS(CSP), DB)
                  EDA = AP*CONV
                  IF (NFCTR .EQ. NFDER) THEN
C
C* return dphi/dr in the force array if required (NFCTR.EQ.NFDER,
C  IDIDA.EQ.IDINCS)
C
                     RKJ = DSQRT(RKJ2)
                     RIJRKJ = 0.0D0
                     RKLRKJ = 0.0D0
                     DO M=1, NCALCD
                        RIJRKJ = RIJRKJ + XXIJ(M)*XXKJ(M)
                        RKLRKJ = RKLRKJ + XXKL(M)*XXKJ(M)
                     ENDDO
                     II3 = 0
                     JJ3 = NDIM
                     KK3 = 2*NDIM
                     LL3 = 3*NDIM
                     DO M=1, NCALCD
                        DPDRI(M) = RKJ/RMJ2 * XXMJ(M)
                        DPDRL(M) = -RKJ/RNK2 * XXNK(M)
                        DPDRJ(M) = (RIJRKJ/RKJ2-1.0D0)*DPDRI(M)
     .                       - RKLRKJ/RKJ2 * DPDRL(M)
                        DPDRK(M) = -DPDRI(M)-DPDRL(M)-DPDRJ(M)

                        F(II3+M) = DPDRI(M)
                        F(JJ3+M) = DPDRJ(M)
                        F(KK3+M) = DPDRK(M)
                        F(LL3+M) = DPDRL(M)
                     ENDDO
!$OMP parallel do private(M)
                     DO M=NDHMIN,NDHMAX
                        F(II3+M) = 0.0D0
                        F(JJ3+M) = 0.0D0
                        F(KK3+M) = 0.0D0
                        F(LL3+M) = 0.0D0
                     ENDDO
                  ENDIF
               ELSE
C
C* NFCTR.EQ.NFNDXE -> EDA returned with energy of NUMDA-th dihedral
C  of the NMOL-th molecule (calculated here for IDIDA.EQ.IDINCS)
C
                  AP = SIGN(DACOS(CSP), RIJRNK)
                  IF (LDOPER) THEN
                     DCB = CPB(NUMDA) - CPA(NUMDA)
                     DB0 = PDB(NUMDA) - PDA(NUMDA)
                     DBA = AP - PDA(NUMDA)
                     NDBA = NINT(DBA*PYETWI)
                     DBA = DBA - NDBA*PYETWO
                     DBL = RLAM*DB0
                     DB1 = DBA - DBL
                     DB2 = DCB*DB1
                     DB3 = DB0*CPA(NUMDA)
                     SP = (DB2-DB3)*RLAM
                     DB5 = DB0*DB3*RLAM
                     DB4 = DB1*DB2 + DB5 - 2.D0*DBA*DB3
                     EDA = 0.5d0*RLAM*DB4
                     EDGH = DB4 + DB5 - 2.D0*DBL*DB2
                     EDAGL = 0.5d0*EDGH
                  ELSE
                     IC = ICP(NUMDA)
                     ARG = AP - PD(IC)
                     NARG = NINT(ARG*PYETWI)
                     ARG = ARG - NARG*PYETWO
                     SP = CP(IC)*ARG
                     EDA = 0.5d0*SP*ARG
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (IDIDA .EQ. IDICOS) THEN
C
C  NFCTR.EQ.NFNDXE -> EDA returned with energy of NUMDA-th dihedral
C  of the NMOL-th molecule (calculated here for IDIDA.EQ.IDICOS)
C
            RIJRKJ = 0.0D0
            RKLRKJ = 0.0D0
            DO M=1, NCALCD
               RIJRKJ = RIJRKJ + XXIJ(M)*XXKJ(M)
               RKLRKJ = RKLRKJ + XXKL(M)*XXKJ(M)
            ENDDO
C
!$OMP parallel do private(M)
            DO M=1, NCALCD
               XXIM(M) = XXIJ(M) - RIJRKJ/RKJ2 * XXKJ(M)
               XXLN(M) =-XXKL(M) + RKLRKJ/RKJ2 * XXKJ(M)
            ENDDO
C
            RIM2 = 0.0D0
            RLN2 = 0.0D0
            DO M=1, NCALCD
               RIM2 = RIM2 + XXIM(M)**2
               RLN2 = RLN2 + XXLN(M)**2
            ENDDO
            RRIMLN = DSQRT(RIM2*RLN2)
C
            IF (RRIMLN.LT.EPS) THEN
C     return zero if angle is not defined               
               EDA = 0.0D0
               EDAGL = 0.0D0
            ELSE
C     angle is defined
               RIMRLN = 0.0D0
               DO M=1, NCALCD
                  RIMRLN = RIMRLN + XXIM(M)*XXLN(M)
               ENDDO
C
               CSP = RIMRLN/RRIMLN
               IF (CSP.GT.1.D0) THEN
                  CSP = 1.D0
               ELSEIF (CSP.LT.-1.D0) THEN
                  CSP = -1.D0
               ENDIF
C
               IF (LDOPER) THEN
                  CALL GETDMP(NPA(NUMDA),CSP,CSMPA,DCMPCA)
                  CALL GETDMP(NPB(NUMDA),CSP,CSMPB,DCMPCB)
                  DUMA = CPA(NUMDA)* (1.D0+PDA(NUMDA)*CSMPA)
                  DUMB = CPB(NUMDA)* (1.D0+PDB(NUMDA)*CSMPB)
                  EDAD = DUMB - DUMA
                  EDA  = RLAM*EDAD
                  EDAGL = EDAD
               ELSE
                  IC = ICP(NUMDA)
                  CALL GETDMP(NP(IC),CSP,CSMP,DCMPDC)                  
                  EDA = CP(IC)* (1.D0+PD(IC)*CSMP)
               ENDIF
            ENDIF
         ENDIF
C
         RETURN
      ENDIF
C
C******** DONE ONE DIHEDRAL CALCULATION
C******** MULTIPLE DIHEDRAL CALCULATION
C
      IF (NFCTR .EQ. NFALL) THEN
         LBNDS = .TRUE.
         LFORCE= .TRUE.
      ELSEIF (NFCTR .EQ. NFEF) THEN
         LBNDS = .FALSE.
         LFORCE = .TRUE.
      ELSEIF (NFCTR .EQ. NFEL) THEN
         LBNDS = .TRUE.
         LFORCE = .FALSE.
      ELSEIF (NFCTR .EQ. NFE) THEN
         LBNDS = .FALSE.
         LFORCE = .FALSE.
       ENDIF
C
      IDAOFF = 0
C
      IF ( .NOT.LDO34 ) THEN
         DO 100 MOL=1,NMOL
C
C** CALCULATION IN 3D (BOTH POTENTIALS)
C 
            IATOFF = (MOL-1)*NAM
            ICOFF  = NDIM*(IATOFF-1)
C
            DO 200 NNA=1,NUMDA
C
               I3 = NDIM*IP(NNA) + ICOFF
               J3 = NDIM*JP(NNA) + ICOFF
               K3 = NDIM*KP(NNA) + ICOFF
               L3 = NDIM*LP(NNA) + ICOFF
C
C** CALCULATE DIHEDRAL, AND DERIVATIVE, CORRECT FOR BOUNDARY 
C   CONDITION IF REQUIRED
C
               DO M=1,NCALCD
                  XJ = X(J3+M)
                  XK = X(K3+M)
                  XXIJ(M) = X(I3+M) - XJ
                  XXKJ(M) = XK - XJ
                  XXKL(M) = XK - X(L3+M)
               ENDDO
C
               RIJ2 = 0.0D0
               RKJ2 = 0.0D0
               RKL2 = 0.0D0
               IF (LVAC) THEN
                  DO M=1,NCALCD
                     RIJ2 = RIJ2 + XXIJ(M)**2
                     RKJ2 = RKJ2 + XXKJ(M)**2
                     RKL2 = RKL2 + XXKL(M)**2
                  ENDDO
               ELSE
                  IF (LDOTRA) THEN
                     CALL TRACO(1,0,XXIJ,BETA,1,LEVERY)
                     CALL TRACO(1,0,XXKJ,BETA,1,LEVERY)
                     CALL TRACO(1,0,XXKL,BETA,1,LEVERY)
                  ENDIF
                  DO M=1,NCALCD
                     BOXM = BOX(M)
                     BOXHM= BOXH(M)
C     
                     DSTTMP = XXIJ(M)
                     IF (DSTTMP .GE. BOXHM) THEN
                        DSTTMP = DSTTMP - BOXM
                     ELSEIF (DSTTMP .LT. -BOXHM) THEN
                        DSTTMP = DSTTMP + BOXM
                     ENDIF
                     XXIJ(M) = DSTTMP
                     RIJ2 = RIJ2 + DSTTMP*DSTTMP
C     
                     DSTTMP = XXKJ(M)
                     IF (DSTTMP .GE. BOXHM) THEN
                        DSTTMP = DSTTMP - BOXM
                     ELSEIF (DSTTMP .LT. -BOXHM) THEN
                        DSTTMP = DSTTMP + BOXM
                     ENDIF
                     XXKJ(M) = DSTTMP
                     RKJ2 = RKJ2 + DSTTMP*DSTTMP
C     
                     DSTTMP = XXKL(M)
                     IF (DSTTMP .GE. BOXHM) THEN
                        DSTTMP = DSTTMP - BOXM
                     ELSEIF (DSTTMP .LT. -BOXHM) THEN
                        DSTTMP = DSTTMP + BOXM
                     ENDIF
                     XXKL(M) = DSTTMP
                     RKL2 = RKL2 + DSTTMP*DSTTMP
                  ENDDO
C     
                  IF (LOCTO) THEN
                     DSTTMP = BOXOQ -DABS(XXIJ(1))-
     .                 DABS(XXIJ(2))-DABS(XXIJ(3))
                     IF (DSTTMP .LT. 0.0D0) THEN
                        RIJ2 = RIJ2 + DSTTMP*BOX(1)
                        DO M=1,NCALCD
                          XXIJ(M)=XXIJ(M)-SIGN(BOXOH,XXIJ(M)) 
                        ENDDO
                     ENDIF
C     
                     DSTTMP = BOXOQ -DABS(XXKJ(1))-
     .                 DABS(XXKJ(2))-DABS(XXKJ(3))
                     IF (DSTTMP .LT. 0.0D0) THEN
                        RKJ2 = RKJ2 + DSTTMP*BOX(1)
                        DO M=1,NCALCD
                          XXKJ(M)=XXKJ(M)-SIGN(BOXOH,XXKJ(M)) 
                        ENDDO
                     ENDIF
C     
                     DSTTMP = BOXOQ -DABS(XXKL(1))-
     .                 DABS(XXKL(2))-DABS(XXKL(3))
                     IF (DSTTMP .LT. 0.0D0) THEN
                        RKL2 = RKL2 + DSTTMP*BOX(1)
                        DO M=1,NCALCD
                          XXKL(M)=XXKL(M)-SIGN(BOXOH,XXKL(M)) 
                        ENDDO
                     ENDIF
C
                  ELSEIF (LDOTRA) THEN
                     RIJ2 = RIJ2 + COSB2*XXIJ(1)*XXIJ(3)
                     RKJ2 = RKJ2 + COSB2*XXKJ(1)*XXKJ(3)
                     RKL2 = RKL2 + COSB2*XXKL(1)*XXKL(3)
                  ENDIF
               ENDIF
C
C** DONE DIHEDRAL AND DERIVATIVE CALCULATION
C   CALCULATE ENERGY/FORCES/LAMBDA DERIVATIVE
C
               IF (IDIDA .EQ. IDINCS) THEN
C
C* HARMONIC IMPROPER DIHEDRAL
C
                  XXMJ(1) = XXIJ(2)*XXKJ(3) - XXIJ(3)*XXKJ(2)
                  XXMJ(2) = XXIJ(3)*XXKJ(1) - XXIJ(1)*XXKJ(3)
                  XXMJ(3) = XXIJ(1)*XXKJ(2) - XXIJ(2)*XXKJ(1)
                  XXNK(1) = XXKJ(2)*XXKL(3) - XXKJ(3)*XXKL(2)
                  XXNK(2) = XXKJ(3)*XXKL(1) - XXKJ(1)*XXKL(3)
                  XXNK(3) = XXKJ(1)*XXKL(2) - XXKJ(2)*XXKL(1)
C
                  RMJ2 = 0.0D0
                  RNK2 = 0.0D0
                  DO M=1, NCALCD
                     RMJ2 = RMJ2 + XXMJ(M)**2
                     RNK2 = RNK2 + XXNK(M)**2
                  ENDDO
C
                  RMJ = DSQRT(RMJ2)
                  RNK = DSQRT(RNK2)
                  RRMJNK = RMJ*RNK
C
                  IF (RRMJNK.GT.EPS) THEN
C     angle is defined, else dont calculate the interactions            
                     RIJRNK = 0.0D0
                     RMJRNK = 0.0D0
                     DO M=1, NCALCD
                        RIJRNK = RIJRNK + XXIJ(M)*XXNK(M)
                        RMJRNK = RMJRNK + XXMJ(M)*XXNK(M)
                     ENDDO
C
                     CSP = RMJRNK/RRMJNK
                     IF (CSP.GT.1.D0) THEN
                        CSP = 1.D0
                     ELSEIF (CSP.LT.-1.D0) THEN
                        CSP = -1.D0
                     ENDIF
C
                     AP = SIGN(DACOS(CSP), RIJRNK)
                     IF (LDOPER) THEN
                        DCB = CPB(NNA) - CPA(NNA)
                        DB0 = PDB(NNA) - PDA(NNA)
                        DBA = AP - PDA(NNA)
                        NDBA = NINT(DBA*PYETWI)
                        DBA = DBA - NDBA*PYETWO
                        DBL = RLAM*DB0
                        DB1 = DBA - DBL
                        DB2 = DCB*DB1
                        DB3 = DB0*CPA(NNA)
                        SP = (DB2-DB3)*RLAM
                        DB5 = DB0*DB3*RLAM
                        DB4 = DB1*DB2 + DB5 - 2.D0*DBA*DB3
                        EDAH = 0.5D0*RLAM*DB4
                        EDA = EDA + EDAH
                        EDGH = DB4 + DB5 - 2.D0*DBL*DB2
                        EDAGL = EDAGL + 0.5D0*EDGH
                     ELSE
                        IC = ICP(NNA)
                        ARG = AP - PD(IC)
                        NARG = NINT(ARG*PYETWI)
                        ARG = ARG - NARG*PYETWO
                        SP = CP(IC)*ARG
                        EDAH = 0.5D0*SP*ARG
                        EDA = EDA + EDAH
                     ENDIF
C
                     IF (LFORCE) THEN
                        IF (LDOTRA) THEN
                           CALL TRACO(1,0,XXIJ,BETA,-1,LEVERY)
                           CALL TRACO(1,0,XXKJ,BETA,-1,LEVERY)
                           CALL TRACO(1,0,XXKL,BETA,-1,LEVERY)
                        ENDIF
                        RKJ = DSQRT(RKJ2)
                        RIJRKJ = 0.0D0
                        RKLRKJ = 0.0D0
                        DO M=1, NCALCD
                           RIJRKJ = RIJRKJ + XXIJ(M)*XXKJ(M)
                           RKLRKJ = RKLRKJ + XXKL(M)*XXKJ(M)
                        ENDDO
                        DO M=1, NCALCD
                           DPDRI(M) = RKJ/RMJ2 * XXMJ(M)
                           DPDRL(M) = -RKJ/RNK2 * XXNK(M)
                           DPDRJ(M) = (RIJRKJ/RKJ2 - 1.D0)*DPDRI(M)
     .                          - RKLRKJ/RKJ2 * DPDRL(M)
                           DPDRK(M) = -DPDRI(M)-DPDRL(M)-DPDRJ(M)
C
                           F(I3+M) = F(I3+M) - SP*DPDRI(M)
                           F(J3+M) = F(J3+M) - SP*DPDRJ(M)
                           F(K3+M) = F(K3+M) - SP*DPDRK(M)
                           F(L3+M) = F(L3+M) - SP*DPDRL(M)
C
                        ENDDO
                     ENDIF
C
                     IF (LBNDS) THEN
                        EP0(IDAOFF+NNA) = EDAH
                        XP0(IDAOFF+NNA) = AP*CONV
                     ENDIF
                  ENDIF                  
               ELSEIF (IDIDA .EQ. IDICOS) THEN
C
C* DONE HARMONIC IMPROPER DIHEDRAL
C* PERIODIC TORTIONAL DIHEDRAL
C
                 XXNK(1) = XXKJ(2)*XXKL(3) - XXKJ(3)*XXKL(2)
                 XXNK(2) = XXKJ(3)*XXKL(1) - XXKJ(1)*XXKL(3)
                 XXNK(3) = XXKJ(1)*XXKL(2) - XXKJ(2)*XXKL(1)
C
                 RIJRKJ = 0.0D0
                 RKLRKJ = 0.0D0
                 DO M=1, NCALCD
                   RIJRKJ = RIJRKJ + XXIJ(M)*XXKJ(M)
                   RKLRKJ = RKLRKJ + XXKL(M)*XXKJ(M)
                 ENDDO
C
!$OMP parallel do private(M)
                 DO M=1, NCALCD
                   XXIM(M) = XXIJ(M) - RIJRKJ/RKJ2 * XXKJ(M)
                   XXLN(M) =-XXKL(M) + RKLRKJ/RKJ2 * XXKJ(M)
                 ENDDO
C
                 RIM2 = 0.0D0
                 RLN2 = 0.0D0
                 DO M=1, NCALCD
                   RIM2 = RIM2 + XXIM(M)**2
                   RLN2 = RLN2 + XXLN(M)**2
                 ENDDO
C
                 RIM = DSQRT(RIM2)
                 RLN = DSQRT(RLN2)
                 RRIMLN = RIM*RLN
C                 
                 IF (RRIMLN.GT.EPS) THEN
                   RIJRNK = 0.0D0
                   DO M=1, NCALCD
                     RIJRNK = RIJRNK + XXIJ(M)*XXNK(M)
                   ENDDO
C     angle is defined, else dont calculate the interactions             
                   RIMRLN = 0.0D0
                   DO M=1, NCALCD
                     RIMRLN = RIMRLN + XXIM(M)*XXLN(M)
                   ENDDO
C                   
                   CSP = RIMRLN/RRIMLN
                   IF (CSP.GT.1.D0) THEN
                      CSP = 1.D0
                   ELSEIF (CSP.LT.-1.D0) THEN
                      CSP = -1.D0
                   ENDIF
                   AP = SIGN(DACOS(CSP),RIJRNK )
C
                   IF (LDOPER) THEN
                     CALL GETDMP(NPA(NNA),CSP,CSMPA,DCMPCA)
                     CALL GETDMP(NPB(NNA),CSP,CSMPB,DCMPCB)
                     DUMA  = CPA(NNA)* (1.D0+PDA(NNA)*CSMPA)
                     DUMB  = CPB(NNA)* (1.D0+PDB(NNA)*CSMPB)
                     EDAD  = DUMB - DUMA
                     EDAH  = RLAM*EDAD
                     EDA   = EDA + EDAH
                     EDAGL = EDAGL + EDAD
                     SPA   = CPA(NNA)*PDA(NNA)*DCMPCA
                     SPB   = CPB(NNA)*PDB(NNA)*DCMPCB
                     SP    = RLAM*(SPB - SPA)
                   ELSE
                     IC    = ICP(NNA)
                     CALL GETDMP(NP(IC),CSP,CSMP,DCMPDC)
                     EDAH  = CP(IC)* (1.D0+PD(IC)*CSMP)
                     EDA   = EDA + EDAH
                     SP    = CP(IC)*PD(IC)*DCMPDC
                   ENDIF
C
                   IF (LFORCE) THEN
                     DO M=1,NCALCD
                       DCPDRI(M) = (XXLN(M)/RLN -
     .                      XXIM(M)/RIM*CSP)/RIM
                       DCPDRL(M) = (XXIM(M)/RIM - 
     .                      XXLN(M)/RLN*CSP)/RLN
                       DCPDRJ(M) = 
     .                      (RIJRKJ/RKJ2 - 1.D0)*DCPDRI(M)
     .                      -(RKLRKJ/RKJ2)*DCPDRL(M)
                       DCPDRK(M) = -DCPDRI(M)
     .                      - DCPDRL(M) - DCPDRJ(M)
C                       
                       F(I3+M) = F(I3+M) - SP*DCPDRI(M)
                       F(J3+M) = F(J3+M) - SP*DCPDRJ(M)
                       F(K3+M) = F(K3+M) - SP*DCPDRK(M)
                       F(L3+M) = F(L3+M) - SP*DCPDRL(M)
C                       
                     ENDDO
                   ENDIF
C
                   IF (LBNDS) THEN
                      EP0(IDAOFF+NNA) = EDAH
                      XP0(IDAOFF+NNA) = AP*CONV
                   ENDIF                  
                  ENDIF
               ENDIF
C
C* DONE PERIODIC TORTIONAL DIHEDRAL
C
 200        CONTINUE
C
            IATOFF = IATOFF + NAM
            IDAOFF  = IDAOFF  + NUMDA
C
 100     CONTINUE
C
C** DONE CALCULATION IN 3D (BOTH POTENTIALS)
C
      ELSE
C
C** CALCULATION WITH 3D-4D COUPLING (PERIODIC TORTIONAL DIHEDRAL ONLY)
C
         DO 101 MOL=1,NMOL
C
            IATOFF = (MOL-1)*NAM
            ICOFF  = IATOFF-1
C
            DO 201 NNA=1,NUMDA
C
               I = IP(NNA) + ICOFF
               J = JP(NNA) + ICOFF
               K = KP(NNA) + ICOFF
               L = LP(NNA) + ICOFF
               IF ( C4D(I+1) .GE. EPSSM .OR. C4D(J+1) .GE. EPSSM .OR.
     $              C4D(K+1) .GE. EPSSM .OR. C4D(L+1) .GE. EPSSM  ) THEN
                  I3 = NDIM*I
                  J3 = NDIM*J
                  K3 = NDIM*K
                  L3 = NDIM*L
C
C** CALCULATE DIHEDRAL, AND DERIVATIVE, CORRECT FOR BOUNDARY 
C   CONDITION IF REQUIRED, IN 3D AND 4D
C
                  DO M=1,4
                     XJ = X(J3+M)
                     XK = X(K3+M)
                     XXIJ(M) = X(I3+M) - XJ
                     XXKJ(M) = XK - XJ
                     XXKL(M) = XK - X(L3+M)
                  ENDDO
C
                  RIJ2 = 0.0D0
                  RKJ2 = 0.0D0
                  RKL2 = 0.0D0
                  IF (LVAC) THEN
                     RIJ2 = XXIJ(1)*XXIJ(1)
     .                    + XXIJ(2)*XXIJ(2)
     .                    + XXIJ(3)*XXIJ(3)
                     RKJ2 = XXKJ(1)*XXKJ(1)
     .                    + XXKJ(2)*XXKJ(2)
     .                    + XXKJ(3)*XXKJ(3)
                     RKL2 = XXKL(1)*XXKL(1)
     .                    + XXKL(2)*XXKL(2)
     .                    + XXKL(3)*XXKL(3)
                     RIJ24D = RIJ2 + XXIJ(4)*XXIJ(4)
                     RKJ24D = RKJ2 + XXKJ(4)*XXKJ(4)
                     RKL24D = RKL2 + XXKL(4)*XXKL(4)
                  ELSE
C no monoclinic case
cmb                     DO 152 M=1,4
cmb                        BOXM = BOX(M)
cmb                        BOXHM= BOXH(M)
cmbC     
cmb                        DSTTMP = XXIJ(M)
cmb                        IF (DSTTMP .GE. BOXHM) THEN
cmb                           DSTTMP = DSTTMP - BOXM
cmb                        ELSEIF (DSTTMP .LT. -BOXHM) THEN
cmb                           DSTTMP = DSTTMP + BOXM
cmb                        ENDIF
cmb                        XXIJ(M) = DSTTMP
cmb                        IF ( M .LE. 3 ) THEN
cmb                           RIJ2 = RIJ2 + DSTTMP**2
cmb                        ELSE
cmb                           RIJ24D = RIJ2 + DSTTMP**2
cmb                        ENDIF
cmbC     
cmb                        DSTTMP = XXKJ(M)
cmb                        IF (DSTTMP .GE. BOXHM) THEN
cmb                           DSTTMP = DSTTMP - BOXM
cmb                        ELSEIF (DSTTMP .LT. -BOXHM) THEN
cmb                           DSTTMP = DSTTMP + BOXM
cmb                        ENDIF
cmb                        XXKJ(M) = DSTTMP
cmbC                        
cmb                        IF ( M .LE. 3 ) THEN
cmb                           RKJ2 = RKJ2 + DSTTMP**2
cmb                        ELSE
cmb                           RKJ24D = RKJ2 + DSTTMP**2
cmb                        ENDIF
cmbC     
cmb                        DSTTMP = XXKL(M)
cmb                        IF (DSTTMP .GE. BOXHM) THEN
cmb                           DSTTMP = DSTTMP - BOXM
cmb                        ELSEIF (DSTTMP .LT. -BOXHM) THEN
cmb                           DSTTMP = DSTTMP + BOXM
cmb                        ENDIF
cmb                        XXKL(M) = DSTTMP
cmb                        
cmb                        IF ( M .LE. 3 ) THEN
cmb                           RKL2 = RKL2 + DSTTMP**2
cmb                        ELSE
cmb                           RKL24D = RKL2 + DSTTMP**2
cmb                        ENDIF
cmb 152                 CONTINUE
cmb
                     DO M=1,3
                        BOXM = BOX(M)
                        BOXHM= BOXH(M)
                        DSTTMP = XXIJ(M)
                        IF (DSTTMP .GE. BOXHM) THEN
                           DSTTMP = DSTTMP - BOXM
                        ELSEIF (DSTTMP .LT. -BOXHM) THEN
                           DSTTMP = DSTTMP + BOXM
                        ENDIF
                        XXIJ(M) = DSTTMP
                        RIJ2 = RIJ2 + DSTTMP*DSTTMP
                        DSTTMP = XXKJ(M)
                        IF (DSTTMP .GE. BOXHM) THEN
                           DSTTMP = DSTTMP - BOXM
                        ELSEIF (DSTTMP .LT. -BOXHM) THEN
                           DSTTMP = DSTTMP + BOXM
                        ENDIF
                        XXKJ(M) = DSTTMP
                        RKJ2 = RKJ2 + DSTTMP*DSTTMP
                        DSTTMP = XXKL(M)
                        IF (DSTTMP .GE. BOXHM) THEN
                           DSTTMP = DSTTMP - BOXM
                        ELSEIF (DSTTMP .LT. -BOXHM) THEN
                           DSTTMP = DSTTMP + BOXM
                        ENDIF
                        XXKL(M) = DSTTMP
                        RKL2 = RKL2 + DSTTMP*DSTTMP
                     ENDDO
cmb                  M=4
                     BOXM = BOX(4)
                     BOXHM= BOXH(4)
                     DSTTMP = XXIJ(4)
                     IF (DSTTMP .GE. BOXHM) THEN
                        DSTTMP = DSTTMP - BOXM
                     ELSEIF (DSTTMP .LT. -BOXHM) THEN
                        DSTTMP = DSTTMP + BOXM
                     ENDIF
                     XXIJ(4) = DSTTMP
                     RIJ24D = RIJ2 + DSTTMP*DSTTMP
                     DSTTMP = XXKJ(4)
                     IF (DSTTMP .GE. BOXHM) THEN
                        DSTTMP = DSTTMP - BOXM
                     ELSEIF (DSTTMP .LT. -BOXHM) THEN
                        DSTTMP = DSTTMP + BOXM
                     ENDIF
                     XXKJ(4) = DSTTMP
                     RKJ24D = RKJ2 + DSTTMP*DSTTMP
                     DSTTMP = XXKL(4)
                     IF (DSTTMP .GE. BOXHM) THEN
                        DSTTMP = DSTTMP - BOXM
                     ELSEIF (DSTTMP .LT. -BOXHM) THEN
                        DSTTMP = DSTTMP + BOXM
                     ENDIF
                     XXKL(4) = DSTTMP
                     RKL24D = RKL2 + DSTTMP*DSTTMP
cmb
C     
                     IF (LOCTO) THEN
                        DSTTMP = BOXOQ -DABS(XXIJ(1))-
     .                    DABS(XXIJ(2))-DABS(XXIJ(3))
                        IF (DSTTMP .LT. 0.0D0) THEN
                           RIJ2 = RIJ2 + DSTTMP*BOX(1)
                           DO M=1,NCALCD
                             XXIJ(M)=XXIJ(M)-SIGN(BOXOH,XXIJ(M)) 
                           ENDDO
                        ENDIF
     
                        DSTTMP = BOXOQ -DABS(XXKJ(1))-
     .                    DABS(XXKJ(2))-DABS(XXKJ(3))
                        IF (DSTTMP .LT. 0.0D0) THEN
                           RKJ2 = RKJ2 + DSTTMP*BOX(1)
                           DO M=1,NCALCD
                              XXKJ(M)=XXKJ(M)-SIGN(BOXOH,XXKJ(M)) 
                           ENDDO
                        ENDIF
     
                        DSTTMP = BOXOQ -DABS(XXKL(1))-
     .                    DABS(XXKL(2))-DABS(XXKL(3))
                        IF (DSTTMP .LT. 0.0D0) THEN
                           RKL2 = RKL2 + DSTTMP*BOX(1)
                           DO M=1,NCALCD
                              XXKL(M)=XXKL(M)-SIGN(BOXOH,XXKL(M)) 
                           ENDDO
                        ENDIF
                     ENDIF
C     no monoclinic case
                  ENDIF
C
                  RIJRKJ = XXIJ(1)*XXKJ(1)
     .                   + XXIJ(2)*XXKJ(2)
     .                   + XXIJ(3)*XXKJ(3)
                  RKLRKJ = XXKL(1)*XXKJ(1)
     .                   + XXKL(2)*XXKJ(2)
     .                   + XXKL(3)*XXKJ(3)
                  RIJRK4 = RIJRKJ + XXIJ(4)*XXKJ(4)
                  RKLRK4 = RKLRKJ + XXKL(4)*XXKJ(4)
C
                  XXIM(1) = XXIJ(1) - RIJRKJ/RKJ2*XXKJ(1)
                  XXIM(2) = XXIJ(2) - RIJRKJ/RKJ2*XXKJ(2)
                  XXIM(3) = XXIJ(3) - RIJRKJ/RKJ2*XXKJ(3)
                  XXLN(1) =-XXKL(1) + RKLRKJ/RKJ2*XXKJ(1)
                  XXLN(2) =-XXKL(2) + RKLRKJ/RKJ2*XXKJ(2)
                  XXLN(3) =-XXKL(3) + RKLRKJ/RKJ2*XXKJ(3)
cmb                  DO M=1, 4
cmb                     XXIM4(M) = XXIJ(M) - RIJRK4/RKJ24D*XXKJ(M)
cmb                     XXLN4(M) =-XXKL(M) + RKLRK4/RKJ24D*XXKJ(M)
cmb                  ENDDO
C
                  RIM2 = XXIM(1)* XXIM(1)
     .                 + XXIM(2)* XXIM(2)
     .                 + XXIM(3)* XXIM(3)
                  RLN2 = XXLN(1)*XXLN(1)
     .                 + XXLN(2)*XXLN(2)
     .                 + XXLN(3)*XXLN(3)
                  RIM24D = 0.0D0
                  RLN24D = 0.0D0
                  DO M=1, 4
                     XXIM4(M) = XXIJ(M) - RIJRK4/RKJ24D*XXKJ(M)
                     XXLN4(M) =-XXKL(M) + RKLRK4/RKJ24D*XXKJ(M)
                     RIM24D = RIM24D + XXIM4(M)**2
                     RLN24D = RLN24D + XXLN4(M)**2
                  ENDDO
C
                  RIM = DSQRT(RIM2)
                  RLN = DSQRT(RLN2)
                  RRIMLN = RIM*RLN
C
                  RIM4D = DSQRT(RIM24D)
                  RLN4D = DSQRT(RLN24D)
                  RRIML4 = RIM4D*RLN4D
C
                  IF (RRIMLN.GT.EPS.AND.RRIML4.GT.EPS) THEN
C     angle is defined, else dont calculate the interactions             
                     RIMRLN = XXIM(1)*XXLN(1)
     .                      + XXIM(2)*XXLN(2)
     .                      + XXIM(3)*XXLN(3)
                     RIMRL4 = 0.0D0
                     DO M=1, 4
                        RIMRL4 = RIMRL4 + XXIM4(M)*XXLN4(M)
                     ENDDO
C
                     CSP = RIMRLN/RRIMLN
                     IF (CSP.GT.1.D0) THEN
                        CSP = 1.D0
                     ELSEIF (CSP.LT.-1.D0) THEN
                        CSP = -1.D0
                     ENDIF
                     CSP4D = RIMRL4/RRIML4
                     IF (CSP4D.GT.1.D0) THEN
                        CSP4D = 1.D0
                     ELSEIF (CSP4D.LT.-1.D0) THEN
                        CSP4D = -1.D0
                     ENDIF
C
C** DONE DIHEDRAL AND DERIVATIVE CALCULATION
C   CALCULATE ENERGY/FORCES/LAMBDA DERIVATIVE IN 3D AND 4D
C
                     IF (LDOPER) THEN
C* 3D
                        CALL GETDMP(NPA(NNA),CSP,CSMPA,DCMPCA)
                        CALL GETDMP(NPB(NNA),CSP,CSMPB,DCMPCB)
                        DUMA   = CPA(NNA)* (1.D0+PDA(NNA)*CSMPA)
                        DUMB   = CPB(NNA)* (1.D0+PDB(NNA)*CSMPB)
C
                        EDGLH3 = DUMB - DUMA
                        EDAH3  = RLAM*EDGLH3
                        SPA    = CPA(NNA)*PDA(NNA)*DCMPCA
                        SPB    = CPB(NNA)*PDB(NNA)*DCMPCB
                        SP3D   = -RLAM*(SPB - SPA)*RMUE
C* 4D
                        CALL GETDMP(NPA(NNA),CSP4D,CSMPA,DCMPCA)
                        CALL GETDMP(NPB(NNA),CSP4D,CSMPB,DCMPCB)
                        DUMA   = CPA(NNA)* (1.D0+PDA(NNA)*CSMPA)
                        DUMB   = CPB(NNA)* (1.D0+PDB(NNA)*CSMPB)
C
                        EDGLH4 = DUMB - DUMA
                        EDAH4  = RLAM*EDGLH4
                        SPA    = CPA(NNA)*PDA(NNA)*DCMPCA
                        SPB    = CPB(NNA)*PDB(NNA)*DCMPCB
                        SP4D   = RLAM*(SPB - SPA)*RMUE
C* sum
                        EDA    = EDA   + RMUE * ( EDAH4  - EDAH3 )
                        EDAGL  = EDAGL + RMUE * ( EDGLH4 - EDGLH3 )
                        EDAGM  = EDAGM + EDAH4 - EDAH3
                     ELSE
                        IC     = ICP(NNA)
C* 3D
                        CALL GETDMP(NP(IC),CSP,CSMP,DCMPDC) 
                        EDAH3  = CP(IC)* (1.D0+PD(IC)*CSMP)
                        SP3D   = -CP(IC)*PD(IC)*DCMPDC*RMUE
C* 4D
                        CALL GETDMP(NP(IC),CSP4D,CSMP,DCMPDC) 
                        EDAH4  = CP(IC)* (1.D0+PD(IC)*CSMP)
                        SP4D   = CP(IC)*PD(IC)*DCMPDC*RMUE
C* sum
                        EDAH  = RMUE * ( EDAH4 - EDAH3 )
                        EDA   = EDA + EDAH
                        EDAGM  = EDAGM + EDAH4 - EDAH3
C                        
                     ENDIF
                     IF (LFORCE) THEN
                        DO M=1,3
                           DCPDRI(M) = (XXLN(M)/RLN - 
     $                          XXIM(M)/RIM*CSP)/RIM
                           DCPDRL(M) = (XXIM(M)/RIM - 
     $                          XXLN(M)/RLN*CSP)/RLN
                           DCPDRJ(M) = (RIJRKJ/RKJ2 - 
     $                          1.D0)*DCPDRI(M)
     $                          -(RKLRKJ/RKJ2)*DCPDRL(M)
                           DCPDRK(M) = -DCPDRI(M) 
     $                          - DCPDRL(M) - DCPDRJ(M)
C                           
                           F(I3+M) = F(I3+M) - SP3D*DCPDRI(M)
                           F(J3+M) = F(J3+M) - SP3D*DCPDRJ(M)
                           F(K3+M) = F(K3+M) - SP3D*DCPDRK(M)
                           F(L3+M) = F(L3+M) - SP3D*DCPDRL(M)
C
                        ENDDO
                        DO M=1,4
                           DCPDRI(M) = (XXLN4(M)/RLN4D 
     $                          - XXIM4(M)/RIM4D*CSP4D)/RIM4D
                           DCPDRL(M) = (XXIM4(M)/RIM4D 
     $                          - XXLN4(M)/RLN4D*CSP4D)/RLN4D
                           DCPDRJ(M) = (RIJRK4/RKJ24D 
     $                          - 1.D0)*DCPDRI(M)
     $                          -(RKLRK4/RKJ24D)*DCPDRL(M)
                           DCPDRK(M) = -DCPDRI(M) 
     $                          - DCPDRL(M) - DCPDRJ(M)
C                           
                           F(I3+M) = F(I3+M) - SP4D*DCPDRI(M)
                           F(J3+M) = F(J3+M) - SP4D*DCPDRJ(M)
                           F(K3+M) = F(K3+M) - SP4D*DCPDRK(M)
                           F(L3+M) = F(L3+M) - SP4D*DCPDRL(M)
C
                        ENDDO
                     ENDIF
                     IF (LBNDS) THEN
                        AP = 0.0D0
C one could return ACOS(RMUE*CSP4D+(1-RMUE)*CSP) but the sign of
C the dihedral is undefined in 4D and one could get the wrong 3D sign
C for RMUE=0.0. Better return a zero.
                        EP0(IDAOFF+NNA) = EDAH
                        XP0(IDAOFF+NNA) = AP*CONV
                     ENDIF                  
                  ENDIF
               ENDIF
C
 201        CONTINUE
            IATOFF = IATOFF + NAM
            IDAOFF  = IDAOFF  + NUMDA
C
 101     CONTINUE
C
C** DONE CALCULATION WITH 3D-4D COUPLING (PERIODIC TORTIONAL DIHEDRAL ONLY)
C
      ENDIF
C
C******** DONE MULTIPLE DIHEDRAL CALCULATION
C
C     end dihang
 500  FORMAT (1X,A6,': IMPLEMENTATION RESTRICTION')
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GETDMP(NP,CSP,CSMP,DCMPDC) 
C     ==--------------------------------------------------------------==
COMSUBR GETDMP
C     
C     SUBROUTINE GETDMP(NP,CSP,CSMP,DCMPDC)
C     
C     GETDMP will supply the trigonometric functions
C     CSMP = cos(m*phi)
C     DCMPDC = d cos(m*phi) / d cos(phi)   and
C     
C     for a given multiplicity m=0,1,2,3,4,5 or 6
C     if NP is not in {m} set execution is stopped
C     
C     NP = multiplicity m
C     CSP = cos(phi)
C     CSMP,DCMPDC see above
C
C     
C     should probably not be a subroutine...
C     ==--------------------------------------------------------------==
C     args  
      INTEGER NP
      real*8 CSP,CSMP,DCMPDC
C     local params
      INTEGER ZERO, ONE,TWO,THREE,FOUR,FIVE,SIX
      PARAMETER (ZERO=0,ONE=1,TWO=2,THREE=3,FOUR=4,FIVE=5,SIX=6)
C     local vars
      real*8 CSP2,CSP4
C     begin
      
C     get cos(m*phi)

      IF (NP.EQ.ZERO) THEN
         CSMP = 1.D0
      ELSEIF (NP.EQ.ONE) THEN
         CSMP = CSP
      ELSEIF (NP.EQ.TWO) THEN
         CSMP = 2.D0*CSP*CSP - 1.D0
      ELSEIF (NP.EQ.THREE) THEN
         CSP2 = CSP*CSP
         CSMP = (4.D0*CSP2 - 3.D0)*CSP
      ELSEIF (NP.EQ.FOUR) THEN
         CSP2 = CSP*CSP
         CSMP = (8.D0*CSP2 - 8.D0)*CSP2 + 1.D0
      ELSEIF (NP.EQ.FIVE) THEN
         CSP2 = CSP*CSP
         CSP4 = CSP2*CSP2
         CSMP = (16.D0*CSP4 - 20.D0*CSP2 + 5.D0)*CSP
      ELSEIF (NP.EQ.SIX) THEN
         CSP2 = CSP*CSP
         CSP4 = CSP2*CSP2
         CSMP = (32.D0*CSP4 - 48.D0*CSP2 + 18.D0)*CSP2 - 1.D0
      ELSE
         PRINT *, 'GETDMP: unvalid multiplicity'
         CALL STOPGM('COBOND','FATAL ERROR')
      ENDIF

C     get d cos(m*phi) / d cos(phi)
      IF (NP.EQ.ZERO) THEN
         DCMPDC = 0.D0
      ELSEIF (NP.EQ.ONE) THEN
         DCMPDC = 1.D0
      ELSEIF (NP.EQ.TWO) THEN
         DCMPDC = 4.D0*CSP
      ELSEIF (NP.EQ.THREE) THEN
         DCMPDC = 12.D0*CSP*CSP - 3.D0
      ELSEIF (NP.EQ.FOUR) THEN
         CSP2 = CSP*CSP
         DCMPDC = (32.D0*CSP2-16.D0)*CSP
      ELSEIF (NP.EQ.FIVE) THEN
         CSP2 = CSP*CSP
         DCMPDC = (80.D0*CSP2 - 60.D0)*CSP2 + 5.D0
      ELSEIF (NP.EQ.SIX) THEN
         CSP2 = CSP*CSP
         CSP4 = CSP2*CSP2
         DCMPDC = (192.D0*CSP4-192.D0*CSP2+36.D0)*CSP
      ELSE
         PRINT *, 'GETDMP: unvalid multiplicity'
      ENDIF

C end getdmp
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==--------------------------------------------------------------==
