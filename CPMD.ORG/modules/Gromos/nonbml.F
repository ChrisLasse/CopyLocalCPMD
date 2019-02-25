Cmb Revised 28 January 2013
#if defined(__WITH_OMP3)
#define __COLLAPSE2 collapse(2)
#else
#define __COLLAPSE2
#endif
#ifndef __VECTOR 
COMMSUBR NONBML
C W.F. VAN GUNSTEREN, GRONINGEN, APR. 1985
C 4D/soft-core/reaction-field by Ilario Tironi & Philippe Huenenberger
C Path-integral code by Salomon Billeter, Zuerich, Aug. 95
C      SUBROUTINE NONBML(NATTOT,NPM,NRAGT,NCALCD,XCOORD,F,XR,
C     $     NSPM,NSP,VIR,
C     $     INB,NSZJNB,JNB,
C     $     MAXNRE,NRELKP,IAGRP,
C     $     RCRF,EPSRF,APPAK,
C     $     EGLEL,EGLRF,EGLRC,EGLLJ,
C     $     EG34EL,EG34RF,EG34RC,EG34LJ,
C     $     NUNRE2,EEL,ERF,ERC,ELJ,
C     $     LDOPER,L3D4D,RLAM,RMU,ALPHLJ,ALPHC,NLAM,MMU,
C     $     LPION,LEVERY,LFREE)
C     NONBML will supply the non-bonded energies and forces for a 
C     12-6 (LJ) and a Coulomb potential, corrected by a reaction-field
C     contribution, for a number of solute (protein) molecules plus solvent 
C     molecules.
C     The potential reads:
C     V(r_ij) = - C6/(r_ij**6) + C12/(r_ij**12)
C               + q_i*q_j/r_ij - (0.5*C1/RCRF**3)*q_i*q_j*(r_ij**2)
C                              - ((1-0.5*C1)/RCRF)*q_i*q_j
C     where RCRF is the radius of the reaction field boundary (in principle
C     L<RCUTL>) and 
C     C1 = ( (2-2*EPSRF)*(1+APPAK*RCRF) - EPSRF*(APPAK*RCRF)**2 ) /
C          ( (1+2*EPSRF)*(1+APPAK*RCRF) + EPSRF*(APPAK*RCRF)**2 )
C
C     is the coefficient for the reaction-field correction beyond RCRF
C     (with L<EPSRF> .GE. 1.0 and L<APPAK> .GE. 0.0).
C     If EPSRF .LE. 0.0, C1 is set to -1, corresponding to an infinite
C     continuum permittivity (a warning will be printed at first call)
C
C     If LDOPER.OR.L3D4D (L<NTG> .NE. L<NTGOFF>), a soft-core type, 
C     L<LAMDA> and L<MU> dependent perturbation potential is applied 
C     ( see manual )
C     It is evaluated for atom pairs involving the atoms II for which
C     L<IPERT>(II) .NE. L<NOPERT> .
C     The alpha coefficients for the soft-core perturbation are set to zero
C     except for atom pairs involving the atoms II for which
C     ISCLJ(II) .NE. 0  in which case alpha(LJ) = L<ALPHLJ>
C     and/or
C     ISCCB(II) .NE. 0  in which case alpha(Coulomb/rx field) = L<ALPHC>
C
C     When LDOPER.AND..NOT.L3D4D (L<NTG> .EQ. L<NTGLAM>), solely the lambda 
C     derivative of the potential will be returned. 
C     When .NOT.LDOPER.AND.L3D4D (L<NTG> .EQ. L<NTGMU>), solely the mu 
C     derivative of the potential will be returned. 
C     When LDOPER.AND.L3D4D (L<NTG> .EQ. L<NTGBOT>), both will be returned. 
C     Derivatives not explicitely required will be returned as zero.
C   NCALCD  LDOPER/L3D4D       calculated
C      
C   3  .F./.F.   Normal interactions in 3D
C   3  .T./.F.   Normal interactions in 3D
C                   For pairs involving perturbed atoms:
C                     Soft-core lambda dependent potential in 3D,
C                     mu has to be zero (if not, set locally to zero)
C   3  .F./.T.   - forbidden -
C   3  .T./.T.   - forbidden -
C   4  .F./.F.   Normal interactions in 4D
C   4  .T./.F.   Normal interactions in 4D
C                   For pairs involving perturbed atoms:
C                     Soft-core lambda dependent potential,
C                     3D-4D mixed at given rmu, lam deriv. returned.
C   4  .F./.T.   Normal interactions in 4D
C                     Soft-core lambda dependent potential,
C                     3D-4D mixed at given rmu, mu deriv. returned.
C   4  .T./.T.   Normal interactions in 4D
C                     Soft-core lambda dependent potential,
C                     3D-4D mixed at given rmu, both deriv. returned.
C  >4  .F./.F.   In principle allowed
C
C     IMPORTANT NOTE:
C
C     With NCALCD.EQ.4, no check is done whether an 
C     atom is in 4D or not according to 4D specification list (C4D array). 
C     This means that all atoms (solute/solvent) which are to be handled 
C     in 3D MUST HAVE ZERO 4thD COORDINATES !
C
C     This, however, wont prevent them from having a 4thD component
C     in the force acting on them, due to interaction with 4D atoms.
C
C     If ALPHLJ.NE.0.0 and for an atom II with ISCLJ(II).NE.0, the algorithm
C     assumes that (either in state A or B), if C6.EQ.0.0, then C12.EQ.0.0.
C     If this is not the case, a non soft-core perturbation of C12 is done.
C     A warning will be printed at first call in this case.
C
C     The forces are calculated in NCALCD dimensions (L<NCALCD>.LE.L<NDIM>).
C     The reaction-field correction is, however, always calculated in 3D.
C     If L<NCALCD>.EQ.3, whatever L<MU>, a local value of zero will be used.
C     A warning will be printed at first call in this case. Additionally,
C     L3D4D (L<NTG> .EQ. L<NTGMU> .OR. L<NTG> .EQ. L<NTGBOT>) will not be 
C     allowed in this case.    
C
C     Periodic boundary conditions can be applied, depending on L<NTB>,
C     where the periodic box can be a truncated octahedron, rectangular
C     or monoclinic. 
C
C     If the virial is to be calculated (ABS(NTB) = L<NTBVIR>), it must
C     be specified whether the solute consists of separate submolecules.
C     This can be done by setting the L<NSPM> > 1 and setting the corresponding
C     values into the L<NSP> array. The virial is only calculated
C     in 3D.
C
C     The solute atom charge groups and solvent molecule pairs are to
C     be specified. The pair-list must be ordered (NI < NJ).
C
C     The algorithm used assumes that the exclusions of atom II, 
C     L<JSNE>(L<KNE>(II)+1)...JSNE(KNE(II)+INE(II)), are ordered monotonously
C     in sequence (no inversions).
C     The same is assumed for the third neighbour list.
C     An error will occur at first call if the ordering is incorrect.
C     A warning will be printed if all excluded/third neighbours of a 
C     given atom are not scanned within the pair-list (which shouldnt
C     happen if RCUTP is reasonably large). The exception is if LFREE
C     is .TRUE. - in which case no warning will be written.
C
C     NATTOT   total number of atoms
C     NPM      number of (identical) solute molecules
C     NRAGT    total number of charge groups
C     NCALCD   the number of dimensions in which to calc forces.
C     XCOORD(NDIM*NATTOT)
C              atom cartesian coordinates.
C     F(NDIM*NATTOT)
C              delivered with the forces calculated added.
C     XR(NDIM*NATTOT)
C              atom cartesian coordinates relative to the
C              submolecular centres of mass.
C              This array is only used if the virial is to be
C              calculated, i.e. if (ABS(NTB)) .EQ. L<NTBVIR>).
C              It is still accessed (read only) in the opposite case 
C              and thus the dummy array used in the call should have 
C              the correct dimension (e.g. XCOORD).
C     NSPM     If the virial is to calculated, then L<NSPM> and L<NSP> are
C              accessed in order to calculate the molecular virial.
C              NSPM is the number of separate (sub)molecules forming 
C              one solute molecule
C     NSP(NSPM+1)
C              atom sequence numbers of the last atoms of the
C              submolecules.
C              Only elements 1..NSPM are considered to be filled
C              with data. We need NSPM+1 elements in order to avoid
C              a double IF statement when checking for submolecules
C     VIR(NDIM)
C              delivered with the components of the virial
C              (ABS(NTB)) .EQ. L<NTBVIR>) - only in 3D !
C     INB(NRAGT)
C              pointerlist of the charge group pairlist.
C              INB(I) specifies the number of charge groups
C              or solvent molecules J with J>I forming a pair
C              with solute charge group or solvent molecule I.
C     NSZJNB   the number of elements in the pairlist
C     JNB(NSZJNB)
C              charge group pairlist. It specifies solute charge
C              groups or solvent molecules J forming a pair with
C              solute charge group or solvent molecule I in 
C              natural order; for each I the pointer in JNB is to
C              incremented in INB(I)
C     MAXNRE   number of energy groups
C     NRELKP(MAXNRE,MAXNRE)   
C              index of pairs of energy groups in energy arrays
C     IAGRP(NATTOT) a CHARACTER array which is used to specify
C              the energy groups of the atoms involved
C     RCRF     radius of the reaction field boundary (in principle
C              RCUTL) and 
C     EPSRF    dielectric permittivity of the reaction field 
C              continuum outside the boundary
C     APPAK    inverse Debye screening length
C     EGLEL    delivered with the value of dV/d(L<LAMBDA>) at fix L<MU>
C              for the Coulomb term if LDOPER = .TRUE.
C     EGLRF    delivered with the value of dV/d(L<LAMBDA>) at fix L<MU>
C              for the distance dependent Poisson-Boltzmann reaction
C              field term if LDOPER = .TRUE.
C     EGLRC    delivered with the value of dV/d(L<LAMBDA>) at fix L<MU>
C              for the distance independent Poisson-Boltzmann reaction
C              field term if LDOPER = .TRUE.
C     EGLLJ    delivered with the value of dV/d(L<LAMBDA>) at fix L<MU>
C              for the Lennard-Jones term if LDOPER = .TRUE.
C     EG34EL   delivered with the 3D to 4D free energy contribution
C              for the Coulomb term if LDOPER = .TRUE. .AND. NCALCD = 4
C     EG34RF   delivered with the 3D to 4D free energy contribution
C              for the distant dependent Poisson-Boltzmann reaction field
C              term if LDOPER = .TRUE. .AND. NCALCD = 4
C     EG34RC   delivered with the 3D to 4D free energy contribution
C              for the distant independent Poisson-Boltzmann reaction field
C              term if LDOPER = .TRUE. .AND. NCALCD = 4
C     EG34LJ   delivered with the 3D to 4D free energy contribution
C              for the Lennard-Jones term if LDOPER = .TRUE. .AND. NCALCD = 4
C     NUNRE2   number of energy groups
C     EEL(NUNRE2)
C              delivered with the Coulomb energy components as
C              defined by L<NRE>
C     ERF(NUNRE2)
C              delivered with the distant dependent
C              Poisson-Boltzmann energy components as defined by L<NRE>
C     ERC(NUNRE2)
C              delivered with the distant independent
C              Poisson-Boltzmann energy components as defined by L<NRE>
C     ELJ(NUNRE2)
C              delivered with the Lennard-Jones energy components
C              as defined by L<NRE>
C     LDOPER   perturbation is applied, lambda derivative computed
C     LDO34D   perturbation is applied, mu derivative computed
C     RLAM     actual value of the lambda parameter
C     RMU      actual value of the mu parameter
C     ALPHLJ   soft-core alpha parameter for the Lennard-Jones component
C     ALPHC    soft-core alpha parameter for the electrostatics component
C     NLAM     exponent of lambda in the Hamiltonian coupling scheme
C     MMU      exponent of mu in the Hamiltonian coupling scheme
C     LPION    a path-integral calculation is performed 
C     LEVERY
C     .FALSE.  a number of quantities depending on BOX, BETA etc.
C              are only calculated at the first SUBR. call.
C     .TRUE.   they are calculated every SUBR. call.
C     LFREE
C     .TRUE.   the pair list and interaction is calculated solely 
C              for perturbed atoms. This is used in free energy 
C              extrapolation.
C
C     NOT passed as arguments, but accessed in common blocks are
C     the following:
C
C     If L<LDOPER> .EQ. .TRUE., then L<RLAM>, L<RMU>, L<IPERT>, L<CGB>
C     and L<ISCC>
C     are accessed in order to calculate a perturbation dependent
C     potential.
C     IPERT(MAXNRP)
C              If IPERT(II) =0, the charge of atom II is not perturbed.
C              otherwise it denotes where the new charge can be found
C              in L<CGB>.
C     CGB(MAXPAT)
C              charges determining the perturbation potential for
C              state B.
COMMEND
C     ==--------------------------------------------------------------==
      SUBROUTINE NONBML(NATTOT,NPM,NRAGT,NCALCD,XCOORD,F,
C     ==--------------------------------------------------------------==
#ifdef HUM_VIR
     $     FELEC,
#endif
     $     XR,
#ifdef EWALD
     $     ALPHA,
#endif
#ifdef EWATCUT
     $     RCUTP,
#endif
     $     NSPM,NSP,VIR,
     $     INB,NSZJNB,JNB,
     $     MAXNRE,NRELKP,IAGRP,
     $     RCRF,EPSRF,APPAK,
     $     EGLEL,EGLRF,EGLRC,EGLLJ,
     $     EG34EL,EG34RF,EG34RC,EG34LJ,
     $     NUNRE2,EEL,ERF,ERC,ELJ,
     $     LDOPER,L3D4D,RLAM,RMU,ALPHLJ,ALPHC,NLAM,MMU,
     $     LPION,LEVERY,LFREE)

      implicit none
C     
C     INCLUDES
C     
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'pertsz.h'
      INCLUDE 'pertar.h'
      INCLUDE 'forcesz.h'
      INCLUDE '../mm_input.inc'
C     
C     ARGUMENTS
C     
      LOGICAL LEVERY,LDOPER,L3D4D,LPION,LFREE
      INTEGER NATTOT,NPM,NRAGT,NCALCD
      INTEGER NSZJNB,INB(NRAGT),JNB(NSZJNB)
      INTEGER NLAM,MMU
      INTEGER MAXNRE,NRELKP(MAXNRE,*)
      INTEGER NSPM,NSP(NSPM+1)
      INTEGER NUNRE2
      real*8 XCOORD(NATTOT*NDIM),F(NATTOT*NDIM),XR(NATTOT*NDIM)
      real*8 VIR(NDIM)
      real*8 RCRF,EPSRF,APPAK
      real*8 RLAM,RMU,ALPHLJ,ALPHC
      real*8 EGLEL,EGLRF,EGLRC,EGLLJ
      real*8 EG34EL,EG34RF,EG34RC,EG34LJ
      real*8 EEL(*),ERF(*),ERC(*),ELJ(*)
      CHARACTER IAGRP(NATTOT)
C     
C     LOCAL PARAMS
C     
      REAL*8 ONE
      PARAMETER (ONE=1.0d+0)
      REAL*8 ZEROC6,ZEROAL
      PARAMETER (ZEROC6=1.0D-20,ZEROAL=1.0D-20)
C     
C     SAVED VARIABLES      
C     
      LOGICAL LMONO,LOCTO,LVAC,LDOVIR,LDOTRA,LFIRST,L4D
      SAVE    LMONO,LOCTO,LVAC,LDOVIR,LDOTRA,LFIRST,L4D
      LOGICAL LPIDOP
      SAVE    LPIDOP
      real*8  PININV
      SAVE    PININV
      real*8 COSB,COSB2,BOXOH,BOXOQ
      SAVE COSB,COSB2,BOXOH,BOXOQ
      real*8 RFF, RFE, RFC, C1, RCRF2
      SAVE RFF, RFE, RFC, C1, RCRF2
C     
C     LOCAL VARIABLES
C     
      LOGICAL LDOIT,LDO14,LJSOLV
      LOGICAL LIZERO, LJZERO
      LOGICAL LPERTL
      LOGICAL LEMIS,L14MIS
      LOGICAL LFULL4,LPART4,LFULL3,LPART3
      INTEGER NRAGP,NI,NJ,NNJ,NRPT
      INTEGER M,NDXGRP
      INTEGER NXTEX,NXTJEX,NENDEX
      INTEGER NXT14,NXTJ14,NEND14
      INTEGER ITMPA,ITMPB,JTYPE
      INTEGER ILPICI,MXNRE2
      REAL*8 BETAR
      real*8 XIJDST(MAXDIM)
      real*8 XH,DF3,DF4TH
      real*8 PIFACT
      CHARACTER PRGSTR*(6)
C     
C     LOCAL VARIABLES FOR IMPROVED PERFORMANCE
C     THE PARAMETERS OF AN ATOM ARE:
C     MASS, COORDINATES CHARGEA, CHARGEB, ITYPEA,ITYPEB,ENERGY GROUP,
C     LOCAL COPY OF IPERT, WHICH DETERMINES WHETHER AN ATOM
C     IS PERTURBED OR NOT.
C     
      LOGICAL LLIPER(MXATCG),LLJPER(MXATCG)
      INTEGER ITYPEA(MXATCG),ITYPEB(MXATCG),IGRP(MXATCG)
      INTEGER ISCLJL(MXATCG),ISCELL(MXATCG)
      INTEGER JSCLJL(MXATCG),JSCELL(MXATCG) 
      INTEGER JTYPEA(MXATCG),JTYPEB(MXATCG),JGRP(MXATCG)
      INTEGER ILNE(MXATCG),KLNE(MXATCG)
      INTEGER ILNE14(MXATCG),KLNE14(MXATCG)
      INTEGER ILPIC(MXATCG), ILPJC(MXATCG)
      real*8 CGIA(MXATCG),CGIB(MXATCG)
      real*8 CGJA(MXATCG),CGJB(MXATCG)
      real*8 XI(MAXLAT),XRI(MAXLAT),FI(MAXLAT)
C     
C     ATOM INTERACTION VARIABLES
C     
      LOGICAL LIPERT,LJPERT
      INTEGER IATOFF,IFTOPA,ILTOPA,IFCORA,ILCORA,IILOC,IILOC3
      INTEGER JATOFF,JFTOPA,JLTOPA,JFCORA,JLCORA,JJLOC,JJLOC3
      INTEGER INUMAT,INUMA3,II
      INTEGER JNUMAT,JNUMA3,JJ
      INTEGER ILATYP,ILBTYP,ILGRP
      INTEGER JLATYP,JLBTYP,JLGRP
      INTEGER IAUX,ITMP,IM,NUMIJ,NJPTR
      INTEGER ITOPCG,JTOPCG,NIPROT,NJPROT
      INTEGER ISCLJI, ISCELI
      INTEGER INTIJA,INTIJB
      INTEGER ICOPTR,JCOPTR
      INTEGER IFCG, ILCG, I, J, ILAST
      INTEGER, allocatable :: NSPT(:)
      SAVE    NSPT

      real*8 CGILA,CGILB
C     
C     PERTURBATION VARIABLES
C     
      real*8 RIJ3D2, VELA4D, VLJA4D, VRFA4D, DFA4D, DFRFA4
      real*8 DLLJA4, DLELA4, DLRFA4, VELB4D
      real*8 VLJB4D, VRFB4D, DFB4D, DFRFB4, DLLJB4, DLELB4
      real*8 DLRFB4
      real*8 VELA3D, VLJA3D, VRFA3D, DFA3D, DFRFA3, DLLJA3
      real*8 DLELA3, DLRFA3, VELB3D, VLJB3D, VRFB3D 
      real*8 DFB3D, DFRFB3, DLLJB3, DLELB3, DLRFB3
      real*8 VEL4D, VRF4D, VLJ4D, DF4D4D, DF4D, DLLJ4D, DLEL4D
      real*8 DLRF4D
      real*8 VEL3D, VRF3D, VLJ3D, DF3D, DLLJ3D, DLEL3D, DLRF3D
      real*8 DDM3EL, DDM3RF, DDM3LJ
      real*8 VRCA,VRCB,VRC,DLRC
      real*8 QRFR
      real*8 SOFTLJ, SOFTEL, CA126, CB126, SFEL, SFLJ
      real*8 RSOFT, RSOFT2, RIJ4, RIJ6, RRF, RRF2, RFSQRT
      real*8 RLMA, RL2A, RLMB, RL2B, RLLA, RLLB, RLB, RLA, RDLB, RDLA
      real*8 RLARMU, RLBRMU
      real*8 RM4D, RM3D, RDM3D, RDM4D
      real*8 RMCST, RMDER
      real*8 QIQJA, QIQJB 
      real*8 RIJ2,RIJINV,RIJIN2,RIJIN6
      real*8 DSTTMP
      real*8 CA12,CA6,CB12,CB6
      real*8 CRA1,CRA12
      real*8 RMULOC
C     
C     DATA
C     
      DATA LFIRST /.TRUE./
      DATA PRGSTR /'NONBML'/
#include "inc_nonbml2.h"
C     
C------------BEGIN NONBML  
C     
C------------DO FIRST CALL SETUP (SWITCHES,"MONOCLINIC" VARIABLES
C     REACTION FIELD VARIABLES, SIMPLE CONSISTENCY CHECKS)     
C

C
C---- IF L3D4D IS ON, WE HAVE TO GO THROUGH THE PERTURBATION CODE
C     ANYWAY
C
      LEMIS = .FALSE.
      L14MIS = .FALSE.
      MXNRE2 = MAXNRE*(MAXNRE+1)/2
C
      IF (LEVERY .OR. LFIRST) THEN

        LMONO = (NTB .GT. 0)
        LOCTO = (NTB. LT. 0)
        LVAC  = (NTB .EQ. 0)
        LDOVIR= (ABS(NTB) .EQ. NTBVIR)
C
        LDOTRA = .FALSE.
        IF (LMONO) THEN
          BETAR = BETA*DATAN(ONE)/45.0D0
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

      LPERTL = (LDOPER.OR.L3D4D)
C AK: 2005/07/22. we need this hack to make MAXINB dynamical.
      IF (LFIRST) THEN
         allocate(NSPT(MAXINB))
      ENDIF
      IF (LFIRST .OR. LEVERY) THEN
         LFIRST = .FALSE.
         L4D = (NCALCD .EQ. 4)
         LPIDOP = .TRUE.
         IF (NPID .NE. 0) THEN
            PININV = ONE / NPID
         ENDIF
C     
C---- PERFORM A FEW FIRST CALL CHECKS
C   
C--   IF LPERTL AND ALPHLJ.NE.0.0 CHECK THAT C6.EQ.0 => C12.EQ.0 FOR ATOMS II 
C     WITH ISCLJ(II).NE.0 AND ANY OTHER ATOM ( STATE A AND B )
C     
        IF (LPERTL .AND. DABS(ALPHLJ).GT.ZEROAL) THEN
          DO 26 II=1,NRP
            IF (IPERT(II) .NE. 0) THEN
              IF (ISCLJ(IPERT(II)) .NE. 0) THEN
                ITMPA = IAC(II)
                ITMPB = IACB(IPERT(II))
                DO 28 JTYPE=1,NRATT
                  INTIJA = MPAC(ITMPA,JTYPE)
                  INTIJB = MPAC(ITMPB,JTYPE)
                  IF (C6(INTIJA) .LT. ZEROC6.AND.
     $                 C12(INTIJA) .GT. ZEROC6) THEN
                     PRINT 900,PRGSTR,II
                     PRINT 910,'A',JTYPE
                  ENDIF
                  IF (C6(INTIJB) .LT. ZEROC6.AND.
     $                 C12(INTIJB) .GT. ZEROC6) THEN
                    PRINT 900,PRGSTR,II
                    PRINT 910,'B',JTYPE
                  ENDIF
                  IF (CS6(INTIJA) .LT. ZEROC6.AND.
     $                 CS12(INTIJA) .GT. ZEROC6) THEN
                    PRINT 900,PRGSTR,II
                    PRINT 915,'A',JTYPE
                  ENDIF
                  IF (CS6(INTIJB) .LT. ZEROC6.AND.
     $                 CS12(INTIJB) .GT. ZEROC6) THEN
                    PRINT 900,PRGSTR,II
                    PRINT 915,'B',JTYPE
                  ENDIF
 28             CONTINUE
              ENDIF
            ENDIF
 26       CONTINUE
C
C--   CALCULATING IN 4D WITH 3D COORDINATES IS WRONG
C     
          IF ( NCALCD.GT.NDIM ) THEN
             PRINT *,PRGSTR,' ERROR: NCALCD.GT.NDIM'
             CALL STOPGM('NONBML','FATAL ERROR')
          ENDIF
       ENDIF

C     
C---- DONE FIRST CALL CHECKS
C     NOW REACTION FIELD VARIABLES
C     
       IF (EPSRF .LE. 0.0D0) THEN
          C1 = -1.0D0
          PRINT *,PRGSTR, ' : WARNING, EPSRF.LE.0 , PERMITTIVITY ',
     $         'SET TO INFINITY'
       ELSE
          C1 = ((2.0D0 - 2.0D0*EPSRF) * (1.0D0 + APPAK*RCRF) 
     $         - EPSRF * (APPAK*RCRF)**2)
     $        /((1.0D0 + 2.0D0*EPSRF) * (1.0D0 + APPAK*RCRF) 
     $         + EPSRF * (APPAK * RCRF)**2)
       ENDIF
       RFF = C1 / RCRF**3
       RFE = -0.5D0*RFF
#ifdef EWALD
       RFF = 0.0D0
       RFE = 0.0D0
#endif       
       RFC = -(1.0D0-0.5D0*C1)/RCRF
       RCRF2 = RCRF*RCRF
       IF (LDOVIR) THEN
C     SET UP GLOBAL SUB-MOLECULE ARRAY
C     DENOTE SOLVENT MOLECULES BY -1
          IF (NSPM .GT. 0) THEN
             IFCG = 1
             ILCG = 1
             DO 1900 I = 1, NSPM
                ILAST = NSP(I)
 1880           IF (ILAST .GT. INC(ILCG)) THEN
                   ILCG = ILCG + 1
                   GOTO 1880
                ENDIF
                IF (ILAST .NE. INC(ILCG)) THEN
                   WRITE(6,*) 'NSP(I) ARRAY NOT CONSISTENT',
     $                  ' WITH CHARGE GROUP DEFINITION !'
                   WRITE(6,*) 'ERROR IN NONBML.F'
                   CALL STOPGM('NONBML','FATAL ERROR')
                ENDIF
                DO J = IFCG, ILCG
                   NSPT(J) = NSP(I)
                ENDDO
                IFCG = ILCG + 1
                ILCG = IFCG
 1900        CONTINUE
          ELSE
!$OMP parallel do private(I)
             DO I = 1, NCAG
                NSPT(I) = 1
             ENDDO
          ENDIF
          DO J = 2, NPM
             DO I = 1, NCAG
                NSPT((J-1)*NCAG+I)=NSPT(I)+ 
     .                  (NSPT(NCAG)*(J-1))
             ENDDO
          ENDDO
C     SOLVENT MOLECULES
          ILCG = NPM * NCAG
!$OMP parallel do private(I) shared(NRAGT,ILCG)
          DO I = 1, (NRAGT - ILCG)
             NSPT(ILCG + I) = -1
          ENDDO
          
       ELSE
!$OMP parallel do private(I)
          DO I = 1, NRAGT
             NSPT(I) = 1
          ENDDO

       ENDIF
       
      ENDIF

C     
C------------DONE FIRST CALL SETUP
C     CALCULATE RMU/RLAMBDA DEPENDANT VARIABLES AND A FEW MORE CHECKS
C
cmb-> Now fix some defaults and initialize variables
      LFULL3 = .FALSE.
      LPART3 = .FALSE.
      LFULL4 = .FALSE.
      LPART4 = .FALSE.
      RL2A = 0.0D0
      RLMA = 0.0D0
      RL2B = 0.0D0
      RLMB = 0.0D0
      RLB  = 0.0D0
      RDLB = 0.0D0
      RLA  = 0.0D0
      RDLA = 0.0D0
      RLLB = 0.0D0
      RLLA = 0.0D0
      RLARMU = 0.0D0
      RLBRMU = 0.0D0
      RM4D  = 0.0D0
      RDM4D = 0.0D0
      RM3D  = 0.0D0
      RDM3D = 0.0D0
      RMCST = 0.0D0
      RMDER = 0.0D0
cmb-> end defaults

      IF (LPERTL) THEN
C
C--   CALCULATING IN 3D DOESNT ALLOW THE CALCULATION OF MU DERIVATIVE
C
        IF (.NOT. L4D .AND. L3D4D) THEN
          PRINT *,PRGSTR,': NCALCD.EQ.3.AND.',
     $         '(NTG.EQ.NTGMU.OR.NTG.EQ.NTGBOT) NOT ALLOWED'
          CALL STOPGM('NONBML','FATAL ERROR')
        ENDIF
C
C--   CALCULATING IN 3D WITH RMU.NE.0 SHOULDNT BE DONE
C
        IF (.NOT. L4D .AND. RMU .NE. 0.0D0) THEN
          PRINT *,PRGSTR,' : WARNING, NCALCD.EQ.3.AND.RMU.NE.0.0'
          PRINT *,'         LOCAL RMU WILL BE SET TO ZERO'
          RMULOC = 0.0D0
        ELSE
          RMULOC = RMU
        ENDIF
C
C--   LOOKS LIKE SOME MACHINES DONT LIKE 0**0 = 1
C     ... WHICH WE UNDERSTAND ...
C
        RL2A = RLAM**2
        RLMA = RL2A + RMULOC**2
        RL2B = (1.0D0 - RLAM)**2
        RLMB = RL2B + RMULOC**2
C     
        IF (RLAM.EQ.0.0D0 .AND. NLAM.EQ.1) THEN
          RLB = 0.0D0
          RDLB = DBLE(NLAM)
        ELSE
          RLB = RLAM**NLAM
          RDLB = DBLE(NLAM)*RLB/RLAM
cmb          RDLB = RLAM**(NLAM - 1) * NLAM
        ENDIF
        IF (RLAM.EQ.1.0D0 .AND. NLAM.EQ.1) THEN
          RLA = 0.0D0
          RDLA = -DBLE(NLAM)
        ELSE
          RLA = (1.0D0 - RLAM)**NLAM          
          RDLA = -DBLE(NLAM)*RLA/(1.0D0-RLAM)
cmb          RDLA = -(1.0D0 - RLAM)**(NLAM-1) * NLAM
        ENDIF
C     
        RLLB = -(1.0D0 - RLAM) * RLB
        RLLA = RLAM * RLA
        RLARMU = RLA * RMULOC
        RLBRMU = RLB * RMULOC
C     
        IF (RMULOC.EQ.0.0D0 .AND. MMU.EQ.1) THEN
          RM4D = 0.0D0
          RDM4D = DBLE(MMU)
        ELSE
          RM4D = RMULOC**MMU
          RDM4D = DBLE(MMU)*RM4D/RMULOC
cmb          RDM4D = RMULOC**(MMU - 1) * MMU
        ENDIF
        IF (RMULOC.EQ.1.0D0 .AND. MMU.EQ.1) THEN
          RM3D = 0.0D0
          RDM3D = -DBLE(MMU)
        ELSE
          RM3D = (1.0D0 - RMULOC)**MMU
          RDM3D = -DBLE(MMU)*RM3D/(1.0D0-RMULOC)
cmb          RDM3D = - (1.0D0 - RMULOC)**(MMU - 1) * MMU
        ENDIF
C
        RMCST = RM3D + RM4D
        RMDER = RDM3D + RDM4D
C     
        IF (.NOT. L4D) THEN
          RDM3D = 0.0D0
          RMCST = 1.0D0
          RMDER = 0.0D0
        ENDIF     
C
C-- YET A FEW MORE LOGICALS FOR SPEED
C
        LFULL4 = (L4D .AND. (RMULOC .NE. 0.0D0))
        LPART4 = (RMULOC .EQ. 0.0D0 .AND. L3D4D .AND. MMU.EQ.1)
        LFULL3 = ((L4D .AND. (RMULOC .NE. 1.0D0)).OR. (.NOT. L4D))
        LPART3 = (RMULOC .EQ. 1.0D0 .AND. L3D4D .AND. MMU.EQ.1)
C
      ENDIF
C     
C------------DONE RMU/RLAMBDA DEPENDANT VARIABLES
C     TRANSFORM COORDINATES IF REQUIRED
C     
      IF (LDOTRA) THEN
        CALL TRACO(NATTOT,0,XCOORD,BETA,1,LEVERY)
        CALL TRACO(NATTOT,0,F,BETA,1,LEVERY)
      ENDIF
C     
C------------INITIALIZE ARRAYS AND GENERAL VARIABLES
C     
      EGLEL = 0.0D0
      EGLRF = 0.0D0
      EGLRC = 0.0D0
      EGLLJ = 0.0D0
C
      EG34EL = 0.0D0
      EG34RF = 0.0D0
      EG34RC = 0.0D0
      EG34LJ = 0.0D0
C     
!$OMP parallel do private(M)
#ifdef _vpp_
!OCL NOALIAS
#endif
      DO M=1,NUNRE2
        EEL(M) = 0.0D0
#ifndef EWALD
        ERF(M) = 0.0D0
        ERC(M) = 0.0D0
#endif
        ELJ(M) = 0.0D0
      ENDDO
C     
      CALL mm_AZZERO(VIR,NDIM)
C     
C------------INITIALIZE A FEW VARIABLES AND POINTERS
C     NRAGP: THE TOTAL NUMBER OF SOLUTE CHARGE GROUPS
C     NRPT : THE TOTAL NUMBER OF SOLUTE ATOMS
C     NJPTR IS USED TO ACCESS THE PAIRLIST
C     IATOFF IS THE ATOM SEQUENCE NUMBER OFFSET FOR 'I ATOMS'
C     ICOPTR IS THE COORDINATE ARRAY OFFSET GOT 'I ATOMS'
C     
      NRAGP = NPM*NCAG
      NRPT  = NPM*NRP
      NJPTR = 0
      IATOFF = 0
      ICOPTR = 0
      NI = 1
C     
C------------DONE INITIALIZE
C     START CALCULATION
C     
C------------LOOP OVER PROTEIN MOLECULES AND PRIMARY (PROTEIN ONLY)
C     CHARGE GROUPS      
C     
      DO 50 NIPROT=1,NPM        
        IFTOPA = 1
        IFCORA = IFTOPA + IATOFF
        DO 200 ITOPCG=1,NCAG
          ILTOPA = INC(ITOPCG)
          ILCORA = ILTOPA + IATOFF          
          INUMAT = ILTOPA - IFTOPA + 1
          INUMA3 = INUMAT*NDIM
C     
C---- PER PRIMARY CHARGE GROUP USE LOCAL ARRAYS
C     (TOPOLOGY INFORMATION, COORDINATES, FORCES )
C     DETERMINE SOLUTE SUBMOLECULES IF NECESSARY
C     
          DO II=1,INUMA3
            XI(II)  = XCOORD(ICOPTR+II)
            FI(II)= 0.0D0
#ifdef HUM_VIR
            FELECI(II) = 0.0D0
#endif
          ENDDO          
          IF (LDOVIR) THEN
            DO II=1,INUMA3
              XRI(II) = XR(ICOPTR+II)
            ENDDO
          ENDIF
          IM = 1
          ITMP = IFCORA-1
          DO II=IFTOPA,ILTOPA
            CGIA(IM)   = CG(II)
            ITYPEA(IM) = IAC(II)
            IGRP(IM)   = ICHAR(IAGRP(ITMP+IM))
C
            ILNE(IM) = INE(II)
            KLNE(IM) = KNE(II)
            ILNE14(IM) = INE14(II)
            KLNE14(IM) = KNE14(II)
            ILPIC(IM) = IPIC(II)
C
            IAUX    = IPERT(II)
            IF (LPERTL .AND. IAUX .NE. NOPERT) THEN
              LLIPER(IM) = .TRUE.
              CGIB(IM)   = CGB(IAUX)
              ITYPEB(IM) = IACB(IAUX)
              ISCLJL(IM) = 1 - ISCLJ(IAUX)
              ISCELL(IM) =  1 - ISCC(IAUX)
            ELSE
              LLIPER(IM) = .FALSE.
              CGIB(IM)   = CGIA(IM)
              ITYPEB(IM) = ITYPEA(IM)
              ISCLJL(IM) = 1
              ISCELL(IM) = 1
            ENDIF
            IM = IM+1
          ENDDO
C     
C---- DONE WITH SETTING UP THE LOCAL ARRAYS     
C     
C----------FIRST CALCULATE INTERACTIONS BETWEEN ATOMS WITHIN 
C     THE CHARGE GROUP NI
C     LOOP OVER ATOMS II OF CHARGE GROUP NI
C     
C     
          IILOC = 1
          IILOC3= 0
          DO 500 II=IFTOPA,ILTOPA-1
            CGILA = CGIA(IILOC)
            CGILB = CGIB(IILOC)
            ILATYP = ITYPEA(IILOC)
            ILBTYP = ITYPEB(IILOC)
            ILPICI = ILPIC(IILOC)
            LIPERT = LLIPER(IILOC)
            ILGRP  = IGRP(IILOC)
            ISCLJI = ISCLJL(IILOC)
            ISCELI = ISCELL(IILOC)
C     
C---- GET READY THE EXCLUSION AND 3RD NEIGHBOUR POINTERS OF ATOM II
C     

            NXTEX = KLNE(IILOC) + 1
            NENDEX = NXTEX + ILNE(IILOC)
            IF ( NXTEX.LT.NENDEX ) THEN
              NXTJEX = JSNE(NXTEX)
            ELSE
              NXTJEX = 0
            ENDIF
C
            NXT14 = KLNE14(IILOC) + 1
            NEND14 = NXT14 + ILNE14(IILOC)
            IF ( NXT14.LT.NEND14 ) THEN
              NXTJ14 = JSNE14(NXT14)
            ELSE
              NXTJ14 = 0
            ENDIF
C     
C---- DONE POINTERS FOR EXCLUSIONS AND 3RD NEIGHBOURS
C     
C--------LOOP OVER JJ ATOMS ( OF THE SAME CHARGE GROUP NI )
C     
            JJLOC = IILOC+1
            JJLOC3= IILOC3+NDIM
            DO 510 JJ=II+1,ILTOPA
C     
C--------CHECK FOR PATH INTEGRAL EXCLUSIONS
C     
              IF (LPION) THEN
                LIZERO = (ILPICI .EQ. 0)
                LJZERO = (ILPIC(JJLOC) .EQ. 0)
                IF (LIZERO .AND. LJZERO) THEN
                  PIFACT = ONE
                  LPIDOP = .TRUE.
                ELSEIF (LIZERO .OR. LJZERO) THEN
                  PIFACT = PININV
                  LPIDOP = .TRUE.
                ELSEIF (ILPICI .EQ. ILPIC(JJLOC)) THEN
                  PIFACT = PININV
                  LPIDOP = .TRUE.
                ELSE
                  PIFACT = 0.0D0
                  LPIDOP = .FALSE.
                ENDIF
              ENDIF
C
              IF ( JJ.EQ.NXTJEX ) THEN
C
C---- JJ EXCLUDED FROM II , SKIP INTERACTION
C     
                NXTEX = NXTEX + 1
                IF ( NXTEX.LT.NENDEX ) THEN
                  NXTJEX = JSNE(NXTEX)
                ELSE
                  NXTJEX = 0
                ENDIF
              ELSEIF (LPIDOP) THEN
C     
C---- JJ NOT EXCLUDED FROM II , CALCULATE INTERACTION...
C     CALCULATE DISTANCE AND CONNECTING VECTOR BETWEEN II AND JJ,
C     CORRECT FOR BOUNDARY CONDITIONS IF NECESSARY
C     
!$OMP parallel do private(M) shared(IILOC3,JJLOC3)
                DO M=1,NCALCD
                  XIJDST(M) = XI(IILOC3+M) - XI(JJLOC3+M)
                ENDDO
                RIJ2 = 0.0D0
                IF (LVAC) THEN
                  DO M=1,NCALCD
                    RIJ2 = RIJ2 + XIJDST(M)**2
                  ENDDO
                  IF (L4D) THEN
                    RIJ3D2 = RIJ2 - XIJDST(4)**2
                  ELSE
                    RIJ3D2 = RIJ2
                  ENDIF
                ELSE
                  DO M=1,NCALCD
                    DSTTMP = XIJDST(M)
                    IF (DSTTMP .GE. BOXH(M)) THEN
                      DSTTMP = DSTTMP - BOX(M)
                    ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
                      DSTTMP = DSTTMP + BOX(M)
                    ENDIF
                    XIJDST(M) = DSTTMP
                    RIJ2 = RIJ2 + DSTTMP*DSTTMP
                  ENDDO
C     
                  IF (L4D) THEN
                    RIJ3D2 = RIJ2 - XIJDST(4)**2
                  ELSE
                    RIJ3D2 = RIJ2
                  ENDIF
C                  
                  IF (LOCTO) THEN
                    DSTTMP = BOXOQ - DABS(XIJDST(1))-
     $                   DABS(XIJDST(2))-DABS(XIJDST(3))
                    IF (DSTTMP .LT. 0.0D0) THEN
                      RIJ2 = RIJ2 + DSTTMP*BOX(1) 
                      RIJ3D2 = RIJ3D2 + DSTTMP*BOX(1)
                      XIJDST(1) = XIJDST(1)-SIGN(BOXOH,XIJDST(1))
                      XIJDST(2) = XIJDST(2)-SIGN(BOXOH,XIJDST(2))
                      XIJDST(3) = XIJDST(3)-SIGN(BOXOH,XIJDST(3))
                    ENDIF
                  ELSEIF (LMONO .AND. LDOTRA) THEN
                    RIJ2 = RIJ2 + COSB2*XIJDST(1)*XIJDST(3)
                    RIJ3D2 = RIJ2
                  ENDIF
                ENDIF
C     
C---- DONE DISTANCE CALCULATION
C     NOW ACTUALLY CALCULATE THE INTERACTION BETWEEN II AND JJ
C     CHECK FOR SPECIAL 1,4 INTERACTIONS
C     
                LDO14 = .FALSE.
                IF (JJ .EQ. NXTJ14) THEN
                  LDO14 = .TRUE.
                  NXT14 = NXT14 + 1
                  IF (NXT14 .LT. NEND14) THEN
                    NXTJ14 = JSNE14(NXT14)
                  ELSE
                    NXTJ14 = 0
                  ENDIF
                ENDIF
C     
C---- DONE 1,4 INTERACTIONS CHECK
C     GET PARAMETERS FOR JJ FROM CG I LOCAL ARRAY (NOT CG J ARRAY!)
C     AND GET INTERACTION CODE
C     
                JLATYP = ITYPEA(JJLOC)
                LJPERT = LLIPER(JJLOC)
                JLGRP  = IGRP(JJLOC)                
                NDXGRP = NRELKP(ILGRP,JLGRP)

C AK FIXME: add check for energy group list.
                IF(NDXGRP.GT.MXNRE2) THEN
                  WRITE(6,*) 'NONBML| NDXGRP=',NDXGRP,
     &                          ' > MXNRE2=', MXNRE2
                  CALL STOPGM('NONBML','ARRAY SIZE MISMATCH')
                ENDIF
                INTIJA = MPAC(ILATYP,JLATYP)
                QIQJA  = CGILA * CGIA(JJLOC)
C     
                IF (LDO14) THEN
                  CA12 = CS12(INTIJA)
                  CA6  =  CS6(INTIJA)
                  QIQJA = QIQJA * scale_14_amber
                ELSE
                  CA12 = C12(INTIJA)
                  CA6  =  C6(INTIJA)
                ENDIF
C     
                IF (LIPERT .OR. LJPERT) THEN
C     
C---- IF EITHER II OR JJ IS PERTURBED, A SOFT-CORE PERTURBATION
C     CALCULATION IS DONE
C     GET STATE B PARAMETERS AND USE COMBINATION RULE FOR SOFT-CORE
C     ALPHA PARAMETERS
C     
                  JLBTYP = ITYPEB(JJLOC)
                  INTIJB = MPAC(ILBTYP,JLBTYP)
                  QIQJB  = CGILB * CGIB(JJLOC)
C     
                  IF (LDO14) THEN
                    CB12 = CS12(INTIJB)
                    CB6  =  CS6(INTIJB)
                  ELSE
                    CB12 = C12(INTIJB)
                    CB6  =  C6(INTIJB)
                  ENDIF
C     
                  SOFTLJ = ALPHLJ * (1 - ISCLJI * ISCLJL(JJLOC))
                  SOFTEL = ALPHC *  (1 - ISCELI * ISCELL(JJLOC))
C     
                  IF (DABS(CA6) .LT. ZEROC6) THEN
                    CA126 = 0.0D0
                  ELSE
                    CA126 = CA12 / CA6
                  ENDIF
                  IF (DABS(CB6) .LT. ZEROC6) THEN
                    CB126 = 0.0D0
                  ELSE
                    CB126 = CB12 / CB6
                  ENDIF
C     
                  IF (LFULL4) THEN
C     
C---- SOFT-CORE PERTURBATION, IF 0<RMULOC<=1.0 
C     CALCULATE FULL 4D CONTRIBUTION
C     
                    RIJ4 = RIJ2 * RIJ2
                    RIJ6 = RIJ4 * RIJ2
C     
C--   4D, STATE A
C     
C     COULOMB
                    SFEL   = SOFTEL * RL2A
                    RSOFT2 = 1.0D0 / (SFEL + RIJ2)
                    RSOFT  = DSQRT(RSOFT2)
                    VELA4D = QIQJA * RSOFT                    
C     REACTION FIELD
                    RRF2   = 1.0D0 / (SFEL + RCRF2)
                    RFSQRT = DSQRT(RRF2)
                    RRF    = RRF2 * RFSQRT
                    DFRFA4 = QIQJA * C1 * RRF
                    QRFR   = DFRFA4 * RIJ3D2
                    VRFA4D = -0.5D0 * QRFR
C     LENNARD-JONES
                    SFLJ   = SOFTLJ * CA126
                    CRA1   = 1.0D0 / (SFLJ * RL2A + RIJ6)
                    CRA12  = CA12 * CRA1
                    VLJA4D = CRA1 * (CRA12 - CA6)
                    DFA4D  = CRA1 * ( CRA1 * CRA12 + VLJA4D)
C     DV/DL
                    DLLJA4 = -2.0D0 * SFLJ * DFA4D
                    DLELA4 = -VELA4D * SOFTEL * RSOFT2
                    DLRFA4 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                    DFA4D = RSOFT2*VELA4D + 6.0D0*RIJ4*DFA4D
C     
C--   4D, STATE B
C     
C     COULOMB
                    SFEL   = SOFTEL * RL2B 
                    RSOFT2 = 1.0D0 / (SFEL + RIJ2)
                    RSOFT  = DSQRT(RSOFT2)
                    VELB4D = QIQJB * RSOFT
C     REACTION FIELD
                    RRF2   =  1.0D0 / (SFEL + RCRF2)
                    RFSQRT = DSQRT(RRF2)
                    RRF    = RRF2 * RFSQRT
                    DFRFB4 = QIQJB * C1 * RRF
                    QRFR   = DFRFB4 * RIJ3D2
                    VRFB4D = -0.5D0 * QRFR
C     LENNARD-JONES
                    SFLJ   = SOFTLJ * CB126
                    CRA1   = 1.0D0 / (SFLJ * RL2B + RIJ6)
                    CRA12  = CB12 * CRA1
                    VLJB4D = CRA1 * (CRA12 - CB6)
                    DFB4D  = CRA1 * ( CRA1 * CRA12 + VLJB4D)
C     DV/DL
                    DLLJB4 = -2.0D0 * SFLJ * DFB4D
                    DLELB4 = -VELB4D * SOFTEL * RSOFT2
                    DLRFB4 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                    DFB4D = RSOFT2*VELB4D + 6.0D0*RIJ4*DFB4D
C     
C--   COMBINE 4D STATE A AND B
C     
                    VEL4D  =    RLB  * VELB4D + RLA  * VELA4D
                    VRF4D  =    RLB  * VRFB4D + RLA  * VRFA4D
                    VLJ4D  =    RLB  * VLJB4D + RLA  * VLJA4D
                    DF4D4D =    RLB  * DFB4D  + RLA  * DFA4D
                    DF4D   =    DF4D4D 
     $                        + RLB  * DFRFB4 + RLA  * DFRFA4
                    DLLJ4D =    RDLB * VLJB4D + RLLB * DLLJB4
     $                        + RDLA * VLJA4D + RLLA * DLLJA4
                    DLEL4D =    RDLB * VELB4D + RLLB * DLELB4
     $                        + RDLA * VELA4D + RLLA * DLELA4
                    DLRF4D =    RDLB * VRFB4D + RLLB * DLRFB4
     $                        + RDLA * VRFA4D + RLLA * DLRFA4
C     
C---- DONE SOFT-CORE 4D (FULL CONTRIBUTION)
C     
                  ELSEIF (LPART4) THEN
C     
C---- SOFT-CORE PERTURBATION, IF 0=RMULOC 
C     CALCULATE JUST 4D ENERGIES FOR MU DERIVATIVE
C     
                    RIJ4 = RIJ2 * RIJ2
                    RIJ6 = RIJ4 * RIJ2
C     
C--   4D, STATE A
C     
C     COULOMB
                    SFEL   = SOFTEL * RL2A
                    VELA4D = QIQJA / DSQRT(SFEL + RIJ2)
C     REACTION FIELD
                    RRF2   = 1.0D0 / (SFEL + RCRF2)
                    RFSQRT = DSQRT(RRF2)
                    VRFA4D = -0.5D0*QIQJA*C1*RRF2*RFSQRT*RIJ3D2
C     LENNARD-JONES
                    CRA1   = 1.0D0 / (SOFTLJ*CA126*RL2A + RIJ6)
                    VLJA4D = CRA1 * (CA12 * CRA1 - CA6)
C     
C--   4D, STATE B
C     
C     COULOMB
                    SFEL   = SOFTEL * RL2B 
                    VELB4D = QIQJB / DSQRT(SFEL + RIJ2)
C     REACTION FIELD
                    RRF2   =  1.0D0 / (SFEL + RCRF2)
                    RFSQRT = DSQRT(RRF2)
                    VRFB4D = -0.5D0*QIQJB*C1*RRF2*RFSQRT*RIJ3D2
C     LENNARD-JONES
                    SFLJ   = SOFTLJ * CB126
                    CRA1   = 1.0D0 / (SFLJ * RL2B + RIJ6)
                    VLJB4D = CRA1 * (CB12 * CRA1 - CB6)
C     
C--   COMBINE 4D STATE A AND B
C     
                    VEL4D = RLB * VELB4D + RLA * VELA4D
                    VRF4D = RLB * VRFB4D + RLA * VRFA4D
                    VLJ4D = RLB * VLJB4D + RLA * VLJA4D
C     
C---- DONE SOFT-CORE 4D (PARTIAL CALCULATION)
C     
                  ENDIF
C     
                  IF (LFULL3) THEN
C     
C---- SOFT-CORE PERTURBATION, IF 0<=RMULOC<1.0
C     CALCULATE FULL 3D CONTRIBUTION
C     
                    RIJ4 = RIJ3D2 * RIJ3D2
                    RIJ6 = RIJ4 * RIJ3D2
C     
C--   3D, STATE A
C     
C     COULOMB
                    SFEL   = SOFTEL * RLMA
                    RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                    RSOFT  = DSQRT(RSOFT2)
                    VELA3D = QIQJA * RSOFT
C     REACTION FIELD
                    RRF2   = 1.0D0 / (SFEL + RCRF2)
                    RFSQRT = DSQRT(RRF2)
                    RRF    = RRF2 * RFSQRT
                    DFRFA3 = QIQJA * C1 * RRF
                    QRFR   = DFRFA3 * RIJ3D2
                    VRFA3D = -0.5D0 * QRFR
C     LENNARD-JONES
                    SFLJ   = SOFTLJ * CA126
                    CRA1   = 1.0D0 / (SFLJ * RLMA + RIJ6)
                    CRA12  = CA12 * CRA1
                    VLJA3D = CRA1 * (CRA12 - CA6)
                    DFA3D  = CRA1 * ( CRA1 * CRA12 + VLJA3D)
C     DV/DL
                    DLLJA3 = -2.0D0 * SFLJ * DFA3D
                    DLELA3 = -VELA3D * SOFTEL * RSOFT2
                    DLRFA3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                    DFA3D = RSOFT2*VELA3D + 6.0D0*RIJ4*DFA3D
C     
C--   3D, STATE B
C     
C     COULOMB
                    SFEL   = SOFTEL * RLMB 
                    RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                    RSOFT  = DSQRT(RSOFT2)
                    VELB3D = QIQJB * RSOFT
C     REACTION FIELD
                    RRF2   = 1.0D0 / (SFEL + RCRF2)
                    RFSQRT = DSQRT(RRF2)
                    RRF    = RRF2 * RFSQRT
                    DFRFB3 = QIQJB * C1 * RRF
                    QRFR   = DFRFB3 * RIJ3D2
                    VRFB3D = -0.5D0 * QRFR
C     LENNARD-JONES
                    SFLJ   = SOFTLJ * CB126
                    CRA1   = 1.0D0 / (SFLJ * RLMB + RIJ6)
                    CRA12  = CB12 * CRA1
                    VLJB3D = CRA1 * (CRA12 - CB6)
                    DFB3D  = CRA1 * ( CRA1 * CRA12 + VLJB3D)
C     DV/DL
                    DLLJB3 = -2.0D0 * SFLJ * DFB3D
                    DLELB3 = -VELB3D * SOFTEL * RSOFT2 
                    DLRFB3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                    DFB3D = RSOFT2*VELB3D + 6.0D0*RIJ4*DFB3D
C     
C--   COMBINE 3D, STATE A AND B
C     
                    VEL3D  =  RLB  * VELB3D + RLA  * VELA3D
                    VRF3D  =  RLB  * VRFB3D + RLA  * VRFA3D
                    VLJ3D  =  RLB  * VLJB3D + RLA  * VLJA3D
                    DF3D   =  RLB  * DFB3D  + RLA  * DFA3D
     $                      + RLB  * DFRFB3 + RLA  * DFRFA3
                    DLLJ3D =  RDLB * VLJB3D + RLLB * DLLJB3
     $                      + RDLA * VLJA3D + RLLA * DLLJA3
                    DLEL3D =  RDLB * VELB3D + RLLB * DLELB3 
     $                      + RDLA * VELA3D + RLLA * DLELA3
                    DLRF3D =  RDLB * VRFB3D + RLLB * DLRFB3
     $                      + RDLA * VRFA3D + RLLA * DLRFA3
                    DDM3EL =  RLBRMU * DLELB3 + RLARMU * DLELA3
                    DDM3RF =  RLBRMU * DLRFB3 + RLARMU * DLRFA3
                    DDM3LJ =  RLBRMU * DLLJB3 + RLARMU * DLLJA3
C     
C---- DONE SOFT-CORE 3D (FULL CONTRIBUTION)
C     
                  ELSEIF (LPART3) THEN
C     
C---- SOFT-CORE PERTURBATION, IF RMULOC=1.0
C     CALCULATE JUST 3D ENERGIES FOR MU DERIVATIVE
C     
                    RIJ4 = RIJ3D2 * RIJ3D2
                    RIJ6 = RIJ4 * RIJ3D2
C     
C--   3D, STATE A
C     
C     COULOMB
                    SFEL   = SOFTEL * RLMA
                    VELA3D = QIQJA / DSQRT(SFEL + RIJ3D2)
C     REACTION FIELD
                    RRF2   = 1.0D0 / (SFEL + RCRF2)
                    RFSQRT = DSQRT(RRF2)
                    VRFA3D = -0.5D0*QIQJA*C1*RRF2*RFSQRT*RIJ3D2
C     LENNARD-JONES
                    SFLJ   = SOFTLJ * CA126
                    CRA1   = 1.0D0 / (SFLJ * RLMA + RIJ6)
                    VLJA3D = CRA1 * (CA12 * CRA1 - CA6)
C     
C--   3D, STATE B
C     
C     COULOMB
                    SFEL   = SOFTEL * RLMB 
                    VELB3D = QIQJB / DSQRT(SFEL + RIJ3D2)
C     REACTION FIELD
                    RRF2   = 1.0D0 / (SFEL + RCRF2)
                    RFSQRT = DSQRT(RRF2)
                    VRFB3D = -0.5D0*QIQJB*C1*RRF2*RFSQRT*RIJ3D2
C     LENNARD-JONES
                    SFLJ   = SOFTLJ * CB126
                    CRA1   = 1.0D0 / (SFLJ * RLMB + RIJ6)
                    VLJB3D = CRA1 * (CB12 * CRA1 - CB6)
C     
C--   COMBINE 3D, STATE A AND B
C     
                    VEL3D  =  RLB  * VELB3D + RLA  * VELA3D
                    VRF3D  =  RLB  * VRFB3D + RLA  * VRFA3D
                    VLJ3D  =  RLB  * VLJB3D + RLA  * VLJA3D
C     
C---- DONE SOFT-CORE 3D (PARTIAL CALCULATION)
C     
                  ENDIF
C     
C---- THE DISTANCE INDEPENDANT REACTION FIELD COMPONENT IS
C     DIMENSIONALITY INDEPENDENT. CALCULATE IT HERE.
C     
                  VRCA = RFC*QIQJA
                  VRCB = RFC*QIQJB
                  VRC  = RLB*VRCB  + RLA*VRCA
                  DLRC = RDLB*VRCB + RDLA*VRCA
C     
C---- COMBINE 4D AND 3D COMPONENTS
C     NOTE THAT IF RMULOC.EQ.0, ONLY 4D ENERGIES ARE CALCULATED BUT
C     RM4D.EQ.0 AND RDM4D.EQ.0. 
C     SIMILARLY IF RMULOC.EQ.1, ONLY 3D ENERGIES ARE CALCULATED BUT
C     RM3D.EQ.0 AND RDM3D.EQ.0.
C AK FIXME: add check for energy group list.
              IF(NDXGRP.GT.MXNRE2) THEN
                WRITE(6,*) 'NONBML| NDXGRP=',NDXGRP,
     &                     ' > MXNRE2=', MXNRE2
                CALL STOPGM('NONBML','ARRAY SIZE MISMATCH')
              ENDIF
C     
                  EEL(NDXGRP) = EEL(NDXGRP) + 
     $                 RM3D * VEL3D + RM4D * VEL4D
#ifndef EWALD
                  ERF(NDXGRP) = ERF(NDXGRP) + 
     $                 RM3D * VRF3D + RM4D * VRF4D
                  ERC(NDXGRP) = ERC(NDXGRP) + 
     $                 RMCST * VRC
#endif
                  ELJ(NDXGRP) = ELJ(NDXGRP) + 
     $                 RM3D * VLJ3D + RM4D * VLJ4D
C
                  EGLEL = EGLEL + 
     $                 RM3D * DLEL3D + RM4D * DLEL4D
                  EGLRF = EGLRF + 
     $                 RM3D * DLRF3D + RM4D * DLRF4D
                  EGLRC = EGLRC +
     $                 RMCST * DLRC
                  EGLLJ = EGLLJ + 
     $                 RM3D * DLLJ3D + RM4D * DLLJ4D
C
                  EG34EL = EG34EL + 
     $                 RDM3D*VEL3D + RM3D*DDM3EL + RDM4D*VEL4D

                  EG34RF = EG34RF +
     $                 RDM3D*VRF3D + RM3D*DDM3RF + RDM4D*VRF4D
                  EG34RC = EG34RC +  RMDER * VRC
                  EG34LJ = EG34LJ +
     $                 RDM3D*VLJ3D + RM3D*DDM3LJ + RDM4D*VLJ4D
C     
                  DF3   = RM3D * DF3D + RM4D * DF4D
                  DF4TH = RM4D * DF4D4D
C     
                ELSE
C     
C---- DONE WITH THE PERTURBED CASE, NOW HANDLE THE CASE WHERE
C     NEITHER II NOR JJ IS PERTURBED
C
                  RIJIN2 = 1.0D0 / RIJ2
                  RIJINV = DSQRT(RIJIN2)
                  RIJIN6 = RIJIN2 * RIJIN2 * RIJIN2
C     
#ifndef EWALD
                  VELA4D = QIQJA * RIJINV
#endif
                  CRA1   = CA12 * RIJIN6
                  CRA12  = CRA1 - CA6
                  VLJA4D = CRA12 * RIJIN6
C     
#ifndef EWALD                 
                  DFA4D = RIJIN2*(VELA4D+6.0D0*(CRA1+CRA12)*RIJIN6) 
#endif
C
C Ewald REAL*8 space contribution - intra charge group
#include "inc_nonbml3.h"
C     
C---- PATH INTEGRAL: SCALE THE VALUES
C     
                  IF (LPION) THEN
                    QIQJA  = QIQJA  * PIFACT
                    VELA4D = VELA4D * PIFACT
                    VLJA4D = VLJA4D * PIFACT
                    DFA4D  = DFA4D  * PIFACT
                  ENDIF
C     
#ifndef EWALD
                  ERF(NDXGRP) = ERF(NDXGRP) +
     $                 QIQJA * RFE * RIJ3D2
                  ERC(NDXGRP) = ERC(NDXGRP) +
     $                 QIQJA * RFC
#endif
                  EEL(NDXGRP) = EEL(NDXGRP) +  VELA4D
                  ELJ(NDXGRP) = ELJ(NDXGRP) +  VLJA4D
C     
                  DF3   = DFA4D + RFF * QIQJA
                  DF4TH = DFA4D
C
                ENDIF
C     
C---- END TEST PERTURBATION
C     CALCULATE AND STORE THE FORCE INTO THE LOCAL FORCE ARRAY
C     NO VIRIAL COMPONENT HERE, WE HAVE INRAMOLECULAR FORCES
C     
                XH = XIJDST(1) * DF3
                FI(IILOC3 + 1) = FI(IILOC3 + 1) + XH
                FI(JJLOC3 + 1) = FI(JJLOC3 + 1) - XH
                XH = XIJDST(2) * DF3
                FI(IILOC3 + 2) = FI(IILOC3 + 2) + XH
                FI(JJLOC3 + 2) = FI(JJLOC3 + 2) - XH
                XH = XIJDST(3) * DF3
                FI(IILOC3 + 3) = FI(IILOC3 + 3) + XH
                FI(JJLOC3 + 3) = FI(JJLOC3 + 3) - XH
#ifdef HUM_VIR
                XH = XIJDST(1) * DF3EL
                FELECI(IILOC3 + 1) = FELECI(IILOC3 + 1) + XH
                FELECI(JJLOC3 + 1) = FELECI(JJLOC3 + 1) - XH
                XH = XIJDST(2) * DF3EL
                FELECI(IILOC3 + 2) = FELECI(IILOC3 + 2) + XH
                FELECI(JJLOC3 + 2) = FELECI(JJLOC3 + 2) - XH
                XH = XIJDST(3) * DF3EL
                FELECI(IILOC3 + 3) = FELECI(IILOC3 + 3) + XH
                FELECI(JJLOC3 + 3) = FELECI(JJLOC3 + 3) - XH
#endif
C
C-- THAT SHOULD WORK FOR 4D (PERT,NOPERT) AND >4D (NOPERT)
C
                DO M=4, NCALCD
                  XH = XIJDST(M) * DF4TH
                  FI(IILOC3 + M) = FI(IILOC3 + M) + XH
                  FI(JJLOC3 + M) = FI(JJLOC3 + M) - XH
                ENDDO
C
#ifdef ATVIR
C Phil, atomic virial - intra charge group contribution
                IF (LDOVIR) THEN                  
                  DO M=1,3
                    XH = XIJDST(M)
                    VIR(M) = VIR(M) - DF3* XH*XH
                  ENDDO
                ENDIF
#endif
C     
C---- END TEST IF EXCLUDED NEIGHBOUR
C     
              ENDIF
              JJLOC3 = JJLOC3 + NDIM
              JJLOC  = JJLOC+1
C     
C------END LOOP OVER JJ ATOMS
C     
 510        CONTINUE
C     
C---- SAVE POINTERS FOR EXCLUDED AND 3RD NEIGHBOURS TO KEEP TRACK
C     OF THE ONES ALREADY SCANNED BY EACH ATOM II OF CHARGE GROUP NI
C     
            KLNE(IILOC) = NXTEX - 1
            ILNE(IILOC) = NENDEX-NXTEX
            KLNE14(IILOC) = NXT14 - 1
            ILNE14(IILOC) = NEND14-NXT14
C     
C---- SAVE POINTERS DONE
C     
            IILOC3 = IILOC3 + NDIM
            IILOC  = IILOC+1
C     
C--------END LOOP OVER II ATOMS
C     
 500      CONTINUE
C     
C----------DONE WITH INTERACTIONS BETWEEN ATOMS WITHIN THE CHARGE GROUP NI
C     
C----------NOW CALCULATE INTERACTION OF CHARGE GROUP NI WITH THE OTHERS
C     CHARGE GROUPS NJ FROM THE PAIR-LIST
C     
          NUMIJ = INB(NI)
C     
C--------LOOP OVER SECONDARY CHARGE GROUPS (EITHER SOLUTE OR SOLVENT)
C     
          DO 600 NNJ= 1,NUMIJ
C     
C---- GET SECONDARY CHARGE GROUP FROM PAIR-LIST
C     PER SECONDARY CHARGE GROUP USE LOCAL ARRAYS
C     (TOPOLOGY INFORMATION), DISTINGUISH NJ SOLVENT OR SOLUTE
C     
            NJPTR = NJPTR + 1
            NJ = JNB(NJPTR)
            LJSOLV = (NJ .GT. NRAGP)
            LDOIT = (LDOVIR .AND.
     $           (NSPT(NI) .NE. NSPT(NJ)))
C            
            IF (LJSOLV) THEN
C     
C--   NJ IS A SOLVENT CHARGE GROUP
C     
              JFTOPA = 1
              JLTOPA = NRAM
              JFCORA = NRPT + NRAM*(NJ-NRAGP-1) + 1
              JLCORA = JFCORA + NRAM - 1
              JNUMAT = NRAM
              JNUMA3 = NDIM*JNUMAT
C     
              IM = 1
              ITMP = JFCORA-1
              DO JJ=JFTOPA,JLTOPA
                CGJA(IM)   = CGS(JJ)
                CGJB(IM)   = CGJA(IM)     
                JTYPEA(IM) = IACS(JJ)
                JTYPEB(IM) = JTYPEA(IM)
                JGRP(IM)   = ICHAR(IAGRP(ITMP+IM))
                LLJPER(IM) = .FALSE.
                JSCLJL(IM) = 1
                JSCELL(IM) = 1
                ILPJC(IM)  = 0
                IM = IM+1
              ENDDO
            ELSE
C     
C--   NJ IS A PROTEIN CHARGE GROUP
C     
              NJPROT = ((NJ-1) / NCAG) + 1
              JTOPCG = MOD(NJ-1,NCAG) + 1
              IF (JTOPCG .EQ. 1) THEN
                JFTOPA = 1
              ELSE
                JFTOPA = INC(JTOPCG-1) + 1
              ENDIF
              JLTOPA = INC(JTOPCG)
              JATOFF = (NJPROT-1)*NRP
              JFCORA = JFTOPA + JATOFF
              JLCORA = JLTOPA + JATOFF
              JNUMAT = JLTOPA - JFTOPA + 1
              JNUMA3 = NDIM* JNUMAT
C     
              IM = 1
              ITMP = JFCORA-1
              DO JJ=JFTOPA,JLTOPA
                CGJA(IM)   = CG(JJ)
                JTYPEA(IM) = IAC(JJ)
                JGRP(IM)   = ICHAR(IAGRP(ITMP+IM))                
                ILPJC(IM)  = IPIC(JJ)
                IAUX       = IPERT(JJ)
                IF (LPERTL .AND. IAUX .NE. NOPERT) THEN
                  LLJPER(IM) = .TRUE.
                  CGJB(IM)   = CGB(IAUX)
                  JTYPEB(IM) = IACB(IAUX)
                  JSCLJL(IM) = 1 - ISCLJ(IAUX)
                  JSCELL(IM) = 1 - ISCC(IAUX)
                ELSE
                  LLJPER(IM) = .FALSE.
                  CGJB(IM)   = CGJA(IM)
                  JTYPEB(IM) = JTYPEA(IM)
                  JSCLJL(IM) = 1
                  JSCELL(IM) = 1
                ENDIF
                IM = IM+1
              ENDDO
            ENDIF
C     
C---- DONE WITH SETTING UP THE LOCAL ARRAYS
C     
            JCOPTR = NDIM*(JFCORA-1)
C     
C--------LOOP OVER ATOMS II OF CHARGE GROUPS NI
C     
            IILOC = 1
            IILOC3= 0
            DO 640 II=IFTOPA,ILTOPA
              CGILA  = CGIA(IILOC)
              CGILB  = CGIB(IILOC)
              ILATYP = ITYPEA(IILOC)
              ILBTYP = ITYPEB(IILOC)
              ILPICI = ILPIC(IILOC)
              LIPERT = LLIPER(IILOC)
              ILGRP  = IGRP(IILOC)
              ISCLJI = ISCLJL(IILOC)
              ISCELI = ISCELL(IILOC)
C     
C---- GET THE EXCLUSION AND 3RD NEIGHBOUR POINTERS OF ATOM II READY
C     
              IF (.NOT. LJSOLV .AND. NIPROT .EQ. NJPROT) THEN
                NXTEX = KLNE(IILOC) + 1
                NENDEX = NXTEX + ILNE(IILOC)
 661            IF ( NXTEX.LT.NENDEX ) THEN
                  NXTJEX = JSNE(NXTEX)
                  IF ( NXTJEX.LT.JFTOPA ) THEN
C     
C--   FOR THE (WEIRD) CASE WHERE AN EXCLUDED ATOM IS NOT IN THE 
C     PAIR-LIST, SET A FLAG AND GO ON... THIS WILL HAPPEN IN
C     PROFEE...
C     
                    LEMIS = .TRUE.
                    NXTEX = NXTEX + 1
                    GOTO 661
                  ENDIF
                ELSE
                  NXTJEX = 0
                ENDIF
C     
                NXT14 = KLNE14(IILOC) + 1
                NEND14 = NXT14 + ILNE14(IILOC)
 662            IF ( NXT14.LT.NEND14 ) THEN
                  NXTJ14 = JSNE14(NXT14)
                  IF ( NXTJ14.LT.JFTOPA ) THEN
C     
C--   FOR THE (WEIRD) CASE WHERE 3RD NEIGHBOUR ATOM IS NOT IN THE 
C     PAIR-LIST, SET A FLAG AND GO ON... THIS WILL HAPPEN IN
C     PROFEE...
C     
                    L14MIS = .TRUE.
                    NXT14 = NXT14 + 1
                    GOTO 662
                  ENDIF
                ELSE
                  NXTJ14 = 0
                ENDIF
              ENDIF
C     
C---- DONE POINTERS FOR EXCLUSIONS AND 3RD NEIGHBOURS
C     
C------LOOP OVER ATOMS JJ OF CHARGE GROUP NJ
C     
              JJLOC  = 1
              JJLOC3 = JCOPTR
              DO 650 JJ=JFTOPA,JLTOPA                
C     
C--------CHECK FOR PATH INTEGRAL EXCLUSIONS
C     
                IF (LPION) THEN
                  LIZERO = (ILPICI .EQ. 0)
                  LJZERO = (ILPJC(JJLOC) .EQ. 0)
                  IF (LIZERO .AND. LJZERO) THEN
                    PIFACT = ONE
                    LPIDOP = .TRUE.
                  ELSEIF (LIZERO .OR. LJZERO) THEN
                    PIFACT = PININV
                    LPIDOP = .TRUE.
                  ELSEIF (ILPICI .EQ. ILPJC(JJLOC)) THEN
                    PIFACT = PININV
                    LPIDOP = .TRUE.
                  ELSE
                    PIFACT = 0.0D0
                    LPIDOP = .FALSE.
                  ENDIF
                ENDIF
C     
                IF (.NOT. LJSOLV .AND. NIPROT .EQ. NJPROT
     $               .AND. JJ.EQ.NXTJEX ) THEN
C     
C---- JJ EXCLUDED FROM II , SKIP INTERACTION
C     
                  NXTEX = NXTEX + 1
                  IF ( NXTEX.LT.NENDEX ) THEN
                    NXTJEX = JSNE(NXTEX)
                  ELSE
                    NXTJEX = 0
                  ENDIF
                ELSEIF (LPIDOP) THEN
C     
C---- JJ NOT EXCLUDED FROM II , CALCULATE INTERACTION...
C     CALCULATE DISTANCE AND CONNECTING VECTOR BETWEEN II AND JJ
C     CORRECT FOR BOUNDARY CONDITIONS IF NECESSARY
C     
!$OMP parallel do private(M) shared(IILOC3,JJLOC3)
                  DO M=1,NCALCD
                    XIJDST(M) = XI(IILOC3+M) - XCOORD(JJLOC3+M)
                  ENDDO
                  RIJ2 = 0.0D0
                  IF (LVAC) THEN
                    DO M=1,NCALCD
                      RIJ2 = RIJ2 + XIJDST(M)**2
                    ENDDO
                    IF (L4D) THEN
                      RIJ3D2 = RIJ2 - XIJDST(4)**2
                    ELSE
                      RIJ3D2 = RIJ2
                    ENDIF
                  ELSE
                    DO M=1,NCALCD
                      DSTTMP = XIJDST(M)
                      IF (DSTTMP .GE. BOXH(M)) THEN
                        DSTTMP = DSTTMP - BOX(M)
                      ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
                        DSTTMP = DSTTMP + BOX(M)
                      ENDIF
                      XIJDST(M) = DSTTMP
                      RIJ2 = RIJ2 + DSTTMP*DSTTMP
                    ENDDO
C     
                    IF (L4D) THEN
                      RIJ3D2 = RIJ2 - XIJDST(4)**2
                    ELSE
                      RIJ3D2 = RIJ2
                    ENDIF
C                    
                    IF (LOCTO) THEN
                      DSTTMP = BOXOQ - DABS(XIJDST(1))-
     $                     DABS(XIJDST(2))-DABS(XIJDST(3))
                      IF (DSTTMP .LT. 0.0D0) THEN
                        RIJ2 = RIJ2 + DSTTMP*BOX(1) 
                        RIJ3D2 = RIJ3D2 + DSTTMP*BOX(1)
                        XIJDST(1)=XIJDST(1)-SIGN(BOXOH,XIJDST(1)) 
                        XIJDST(2)=XIJDST(2)-SIGN(BOXOH,XIJDST(2)) 
                        XIJDST(3)=XIJDST(3)-SIGN(BOXOH,XIJDST(3)) 
                      ENDIF
                    ELSEIF (LMONO .AND. LDOTRA) THEN
                      RIJ2 = RIJ2 + COSB2*XIJDST(1)*XIJDST(3)
                      RIJ3D2 = RIJ2
                    ENDIF
                  ENDIF
C     
C---- DONE DISTANCE CALCULATION
C     NOW ACTUALLY CALCULATE THE INTERACTION BETWEEN II AND JJ
C     CHECK FOR SPECIAL 1,4 INTERACTIONS
C     
                  LDO14 = .FALSE.
                  IF (.NOT. LJSOLV .AND. NIPROT .EQ. NJPROT) THEN
                    IF ( JJ.EQ.NXTJ14 ) THEN
                      LDO14 = .TRUE.
                      NXT14 = NXT14 + 1
                      IF ( NXT14.LT.NEND14 ) THEN
                        NXTJ14 = JSNE14(NXT14)
                      ELSE
                        NXTJ14 = 0
                      ENDIF
                    ENDIF
                  ENDIF
C     
C---- DONE 1,4 INTERACTIONS CHECK
C     GET PARAMETERS FOR JJ FROM CG J LOCAL ARRAY
C     AND GET INTERACTION CODE
C     
                  JLATYP = JTYPEA(JJLOC)
                  LJPERT = LLJPER(JJLOC)
                  JLGRP  = JGRP(JJLOC)                  
                  NDXGRP = NRELKP(ILGRP,JLGRP)                  
                  INTIJA = MPAC(ILATYP,JLATYP)
                  QIQJA  = CGILA * CGJA(JJLOC) 
C     
                  IF (LDO14) THEN
                    CA12 = CS12(INTIJA)
                    CA6  =  CS6(INTIJA)
                    QIQJA = QIQJA * scale_14_amber
                  ELSE
                    CA12 = C12(INTIJA)
                    CA6  =  C6(INTIJA)
                  ENDIF
C     
                  IF (LIPERT .OR. LJPERT) THEN
C     
C---- IF EITHER II OR JJ IS PERTURBED, A SOFT-CORE PERTURBATION
C     CALCULATION IS DONE
C     GET STATE B PARAMETERS AND USE COMBINATION RULE FOR SOFT-CORE
C     ALPHA PARAMETERS
C     
                    JLBTYP = JTYPEB(JJLOC)
                    INTIJB = MPAC(ILBTYP,JLBTYP)
                    QIQJB  = CGILB * CGJB(JJLOC)
C     
                    IF (LDO14) THEN
                      CB12 = CS12(INTIJB)
                      CB6  =  CS6(INTIJB)
                    ELSE
                      CB12 = C12(INTIJB)
                      CB6  =  C6(INTIJB)
                    ENDIF
C     
                    SOFTLJ = ALPHLJ * (1 - ISCLJI * JSCLJL(JJLOC))
                    SOFTEL = ALPHC *  (1 - ISCELI * JSCELL(JJLOC))
C     
                    IF (DABS(CA6) .LT. ZEROC6) THEN
                      CA126 = 0.0D0
                    ELSE
                      CA126 = CA12 / CA6
                    ENDIF
                    IF (DABS(CB6) .LT. ZEROC6) THEN
                      CB126 = 0.0D0
                    ELSE
                      CB126 = CB12 / CB6
                    ENDIF
C     
                    IF (LFULL4) THEN
C     
C---- SOFT-CORE PERTURBATION, IF 0<RMULOC<=1.0 
C     CALCULATE FULL 4D CONTRIBUTION
C     
                      RIJ4 = RIJ2 * RIJ2
                      RIJ6 = RIJ4 * RIJ2
C     
C--   4D, STATE A
C     
C     COULOMB
                      SFEL   = SOFTEL * RL2A
                      RSOFT2 = 1.0D0 / (SFEL + RIJ2)
                      RSOFT  = DSQRT(RSOFT2)
                      VELA4D = QIQJA * RSOFT
C     REACTION FIELD
                      RRF2   = 1.0D0 / (SFEL + RCRF2)
                      RFSQRT = DSQRT(RRF2)
                      RRF    = RRF2 * RFSQRT
                      DFRFA4 = QIQJA * C1 * RRF
                      QRFR   = DFRFA4 * RIJ3D2
                      VRFA4D = -0.5D0 * QRFR
C     LENNARD-JONES
                      SFLJ   = SOFTLJ * CA126
                      CRA1   = 1.0D0 / (SFLJ * RL2A + RIJ6)
                      CRA12  = CA12 * CRA1
                      VLJA4D = CRA1 * (CRA12 - CA6)
                      DFA4D  = CRA1 * ( CRA1 * CRA12 + VLJA4D)
C     DV/DL
                      DLLJA4 = -2.0D0 * SFLJ * DFA4D
                      DLELA4 = -VELA4D * SOFTEL * RSOFT2
                      DLRFA4 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                      DFA4D = RSOFT2*VELA4D + 6.0D0*RIJ4*DFA4D
C     
C--   4D, STATE B
C     
C     COULOMB
                      SFEL   = SOFTEL * RL2B
                      RSOFT2 = 1.0D0 / (SFEL + RIJ2)
                      RSOFT  = DSQRT(RSOFT2)
                      VELB4D = QIQJB * RSOFT
C     REACTION FIELD
                      RRF2   = 1.0D0 / (SFEL + RCRF2)
                      RFSQRT = DSQRT(RRF2)
                      RRF    = RRF2 * RFSQRT
                      DFRFB4 = QIQJB * C1 * RRF
                      QRFR   = DFRFB4 * RIJ3D2
                      VRFB4D = -0.5D0 * QRFR
C     LENNARD-JONES
                      SFLJ   = SOFTLJ * CB126 
                      CRA1   = 1.0D0 / (SFLJ * RL2B + RIJ6)
                      CRA12  = CB12 * CRA1
                      VLJB4D = CRA1 * (CRA12 - CB6)
                      DFB4D  = CRA1 * ( CRA1 * CRA12 + VLJB4D)
C     DV/DL
                      DLLJB4 = -2.0D0 * SFLJ * DFB4D
                      DLELB4 = -VELB4D * SOFTEL * RSOFT2 
                      DLRFB4 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                      DFB4D = RSOFT2*VELB4D + 6.0D0*RIJ4*DFB4D
C     
C--   COMBINE 4D STATE A AND B
C     
                      VEL4D  = RLB * VELB4D + RLA * VELA4D
                      VRF4D  = RLB * VRFB4D + RLA * VRFA4D
                      VLJ4D  = RLB * VLJB4D + RLA * VLJA4D
                      DF4D4D = RLB * DFB4D  + RLA * DFA4D
                      DF4D   = DF4D4D 
     .                       + RLB * DFRFB4 + RLA * DFRFA4
                      DLLJ4D = RDLB * VLJB4D + RLLB * DLLJB4
     .                       + RDLA * VLJA4D + RLLA * DLLJA4
                      DLEL4D = RDLB * VELB4D + RLLB * DLELB4
     .                       + RDLA * VELA4D + RLLA * DLELA4
                      DLRF4D = RDLB * VRFB4D + RLLB * DLRFB4
     .                       + RDLA * VRFA4D + RLLA * DLRFA4
C     
C---- DONE SOFT-CORE 4D (FULL CONTRIBUTION)
C     
                    ELSEIF (LPART4) THEN
C     
C---- SOFT-CORE PERTURBATION, IF 0=RMULOC
C     CALCULATE JUST 4D ENERGIES FOR MU DERIVATIVE
C     
                      RIJ4 = RIJ2 * RIJ2
                      RIJ6 = RIJ4 * RIJ2
C     
C--   4D, STATE A
C     
C     COULOMB
                    SFEL   = SOFTEL * RL2A
                    VELA4D = QIQJA / DSQRT(SFEL + RIJ2)
C     REACTION FIELD
                    RRF2   = 1.0D0 / (SFEL + RCRF2)
                    RFSQRT = DSQRT(RRF2)
                    VRFA4D = -0.5D0*QIQJA*C1*RRF2*RFSQRT*RIJ3D2
C     LENNARD-JONES
                    CRA1   = 1.0D0 / (SOFTLJ*CA126*RL2A + RIJ6)
                    VLJA4D = CRA1 * (CA12*CRA1 - CA6)
C     
C--   4D, STATE B
C     
C     COULOMB
                    SFEL   = SOFTEL * RL2B 
                    VELB4D = QIQJB / DSQRT(SFEL + RIJ2)
C     REACTION FIELD
                    RRF2   =  1.0D0 / (SFEL + RCRF2)
                    RFSQRT = DSQRT(RRF2)
                    VRFB4D = -0.5D0*QIQJB*C1*RRF2*RFSQRT*RIJ3D2
C     LENNARD-JONES
                    SFLJ   = SOFTLJ * CB126
                    CRA1   = 1.0D0 / (SFLJ * RL2B + RIJ6)
                    VLJB4D = CRA1 * (CB12 * CRA1 - CB6)
C     
C--   COMBINE 4D STATE A AND B
C     
                    VEL4D = RLB  * VELB4D + RLA  * VELA4D
                    VRF4D = RLB  * VRFB4D + RLA  * VRFA4D
                    VLJ4D = RLB  * VLJB4D + RLA  * VLJA4D
C     
C---- DONE SOFT-CORE 4D (PARTIAL CALCULATION)
C     
                    ENDIF
C     
                    IF (LFULL3) THEN                 
C     
C---- SOFT-CORE PERTURBATION, IF 0<=RMULOC<1.0 
C     CALCULATE FULL 3D CONTRIBUTION
C     
                      RIJ4 = RIJ3D2 * RIJ3D2
                      RIJ6 = RIJ4 * RIJ3D2
C     
C--   3D, STATE A
C     
C     COULOMB
                      SFEL   = SOFTEL * RLMA
                      RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                      RSOFT  = DSQRT(RSOFT2)
                      VELA3D = QIQJA * RSOFT
C     REACTION FIELD
                      RRF2   = 1.0D0 / (SFEL + RCRF2)
                      RFSQRT = DSQRT(RRF2)
                      RRF    = RRF2 * RFSQRT
                      DFRFA3 = QIQJA * C1 * RRF
                      QRFR   = DFRFA3 * RIJ3D2
                      VRFA3D = -0.5D0 * QRFR
C     LENNARD-JONES
                      SFLJ   = SOFTLJ * CA126
                      CRA1   = 1.0D0 / (SFLJ * RLMA + RIJ6)
                      CRA12  = CA12 * CRA1
                      VLJA3D = CRA1 * (CRA12 - CA6)
                      DFA3D  = CRA1 * ( CRA1 * CRA12 + VLJA3D)
C     DV/DL
                      DLLJA3 = -2.0D0 * SFLJ * DFA3D
                      DLELA3 = -VELA3D * SOFTEL * RSOFT2 
                      DLRFA3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                      DFA3D = RSOFT2*VELA3D + 6.0D0*RIJ4*DFA3D
C     
C--   3D, STATE B
C     
C     COULOMB
                      SFEL   = SOFTEL * RLMB
                      RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                      RSOFT  = DSQRT(RSOFT2)
                      VELB3D = QIQJB * RSOFT
C     REACTION FIELD
                      RRF2   = 1.0D0 / (SFEL + RCRF2)
                      RFSQRT = DSQRT(RRF2)
                      RRF    = RRF2 * RFSQRT
                      DFRFB3 = QIQJB * C1 * RRF
                      QRFR   = DFRFB3 * RIJ3D2
                      VRFB3D = -0.5D0 * QRFR
C     LENNARD-JONES
                      SFLJ   = SOFTLJ * CB126 
                      CRA1   = 1.0D0 / (SFLJ * RLMB + RIJ6)
                      CRA12  = CB12 * CRA1
                      VLJB3D = CRA1 * (CRA12 - CB6)
                      DFB3D  = CRA1 * ( CRA1 * CRA12 + VLJB3D)
C     DV/DL
                      DLLJB3 = -2.0D0 * SFLJ * DFB3D
                      DLELB3 = -VELB3D * SOFTEL * RSOFT2 
                      DLRFB3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                      DFB3D = RSOFT2*VELB3D + 6.0D0*RIJ4*DFB3D

C     
C--   COMBINE 3D, STATE A AND B
C     
                      VEL3D =  RLB  * VELB3D + RLA  * VELA3D
                      VRF3D =  RLB  * VRFB3D + RLA  * VRFA3D
                      VLJ3D =  RLB  * VLJB3D + RLA  * VLJA3D
                      DF3D  =  RLB  * DFB3D  + RLA  * DFA3D
     $                       + RLB  * DFRFB3 + RLA  * DFRFA3
                      DLLJ3D = RDLB * VLJB3D + RLLB * DLLJB3
     $                       + RDLA * VLJA3D + RLLA * DLLJA3
                      DLEL3D = RDLB * VELB3D + RLLB * DLELB3 
     $                       + RDLA * VELA3D + RLLA * DLELA3
                      DLRF3D = RDLB * VRFB3D + RLLB * DLRFB3
     $                       + RDLA * VRFA3D + RLLA * DLRFA3
                      DDM3EL = RLBRMU *DLELB3 + RLARMU *DLELA3
                      DDM3RF = RLBRMU *DLRFB3 + RLARMU *DLRFA3
                      DDM3LJ = RLBRMU *DLLJB3 + RLARMU *DLLJA3

C     
C---- DONE SOFT-CORE 3D (FULL CONTRIBUTION)
C     
                    ELSEIF (LPART3) THEN
C     
C---- SOFT-CORE PERTURBATION, IF RMULOC=1.0 
C     CALCULATE JUST 3D ENERGIES FOR MU DERIVATIVE
C     
                      RIJ4 = RIJ3D2 * RIJ3D2
                      RIJ6 = RIJ4 * RIJ3D2
C     
C--   3D, STATE A
C     
C     COULOMB
                    SFEL   = SOFTEL * RLMA
                    VELA3D = QIQJA / DSQRT(SFEL + RIJ3D2)
C     REACTION FIELD
                    RRF2   = 1.0D0 / (SFEL + RCRF2)
                    RFSQRT = DSQRT(RRF2)
                    VRFA3D = -0.5D0*QIQJA*C1*RRF2*RFSQRT*RIJ3D2
C     LENNARD-JONES
                    SFLJ   = SOFTLJ * CA126
                    CRA1   = 1.0D0 / (SFLJ * RLMA + RIJ6)
                    VLJA3D = CRA1 * (CA12 * CRA1 - CA6)
C     
C--   3D, STATE B
C     
C     COULOMB
                    SFEL   = SOFTEL * RLMB 
                    VELB3D = QIQJB / DSQRT(SFEL + RIJ3D2)
C     REACTION FIELD
                    RRF2   = 1.0D0 / (SFEL + RCRF2)
                    RFSQRT = DSQRT(RRF2)
                    VRFB3D = -0.5D0*QIQJB*C1*RRF2*RFSQRT*RIJ3D2
C     LENNARD-JONES
                    SFLJ   = SOFTLJ * CB126
                    CRA1   = 1.0D0 / (SFLJ * RLMB + RIJ6)
                    VLJB3D = CRA1 * (CB12 * CRA1 - CB6)
C     
C--   COMBINE 3D, STATE A AND B
C     
                    VEL3D = RLB * VELB3D + RLA * VELA3D
                    VRF3D = RLB * VRFB3D + RLA * VRFA3D
                    VLJ3D = RLB * VLJB3D + RLA * VLJA3D
C     
C---- DONE SOFT-CORE 3D (PARTIAL CALCULATION)
C     
                    ENDIF
C     
C---- THE DISTANCE INDEPENDANT REACTION FIELD COMPONENT IS
C     DIMENSIONALITY INDEPENDANT. CALCULATE IT HERE.
C     
                    VRCA = RFC*QIQJA
                    VRCB = RFC*QIQJB
                    VRC  = RLB*VRCB  + RLA*VRCA
                    DLRC = RDLB*VRCB + RDLA*VRCA
C     
C---- COMBINE 4D AND 3D COMPONENTS
C     NOTE THAT IF RMULOC.EQ.0, ONLY 4D ENERGIES ARE CALCULATED BUT
C     RM4D.EQ.0 AND RDM4D.EQ.0. 
C     SIMILARLY IF RMULOC.EQ.1, ONLY 3D ENERGIES ARE CALCULATED BUT
C     RM3D.EQ.0 AND RDM3D.EQ.0.
C     
                    EEL(NDXGRP) = EEL(NDXGRP) +
     $                   RM3D * VEL3D + RM4D * VEL4D
#ifndef EWALD
                    ERF(NDXGRP) = ERF(NDXGRP) +
     $                   RM3D * VRF3D + RM4D * VRF4D
                    ERC(NDXGRP) = ERC(NDXGRP) +
     $                   RMCST * VRC
#endif
                    ELJ(NDXGRP) = ELJ(NDXGRP) +
     $                   RM3D * VLJ3D + RM4D * VLJ4D
C
                    EGLEL = EGLEL +
     $                   RM3D * DLEL3D + RM4D * DLEL4D
                    EGLRF = EGLRF +
     $                   RM3D * DLRF3D + RM4D * DLRF4D
                    EGLRC = EGLRC +
     $                   RMCST * DLRC
                    EGLLJ = EGLLJ +
     $                   RM3D * DLLJ3D + RM4D * DLLJ4D
C     
                    EG34EL = EG34EL +
     $                   RDM3D*VEL3D + RM3D*DDM3EL + RDM4D*VEL4D
                    EG34RF = EG34RF +
     $                   RDM3D*VRF3D + RM3D*DDM3RF + RDM4D*VRF4D
                    EG34RC = EG34RC +  RMDER * VRC
                    EG34LJ = EG34LJ +
     $                   RDM3D*VLJ3D + RM3D*DDM3LJ + RDM4D*VLJ4D
C     
                  DF3   = RM3D * DF3D + RM4D * DF4D
                  DF4TH = RM4D * DF4D4D
C 
                  ELSE
C     
C---- DONE WITH THE PERTURBED CASE, NOW HANDLE THE CASE WHERE
C     NEITHER II NOR JJ IS PERTURBED
C     
                    RIJIN2 = 1.0D0 / RIJ2
                    RIJINV = DSQRT(RIJIN2)
                    RIJIN6 = RIJIN2 * RIJIN2 * RIJIN2
C     
#ifndef EWALD
                    VELA4D = QIQJA * RIJINV                    
#endif
                    CRA1   = CA12 * RIJIN6
                    CRA12  = CRA1 - CA6
                    VLJA4D = CRA12 * RIJIN6
#ifndef EWALD
                    DFA4D = RIJIN2*(VELA4D+6.0D0*(CRA1+CRA12)*RIJIN6) 
#endif
C
C Ewald REAL*8 space contribution - solute-solute(\intra charge group),solute-solvent
#include "inc_nonbml3.h"
C     
C---- PATH INTEGRAL: SCALE THE VALUES
C     
                    IF (LPION) THEN
                      QIQJA  = QIQJA  * PIFACT
                      VELA4D = VELA4D * PIFACT
                      VLJA4D = VLJA4D * PIFACT
                      DFA4D  = DFA4D  * PIFACT
                    ENDIF
C     
#ifndef EWALD
                    ERF(NDXGRP) = ERF(NDXGRP) + 
     $                   QIQJA * RFE * RIJ3D2
                    ERC(NDXGRP) = ERC(NDXGRP) + 
     $                   QIQJA * RFC
#endif
                    EEL(NDXGRP) = EEL(NDXGRP) +  VELA4D
                    ELJ(NDXGRP) = ELJ(NDXGRP) +  VLJA4D
C     
                    DF3   = DFA4D + RFF * QIQJA
                    DF4TH = DFA4D
C
                  ENDIF
C
C---- END TEST PERTURBATION
C     CALCULATE AND STORE THE FORCE INTO THE LOCAL FORCE ARRAY
C     
                  XH = XIJDST(1) * DF3
                  FI(IILOC3 + 1) = FI(IILOC3 + 1) + XH
                  F(JJLOC3 + 1)  = F(JJLOC3 + 1) - XH
                  XH = XIJDST(2) * DF3
                  FI(IILOC3 + 2) = FI(IILOC3 + 2) + XH
                  F(JJLOC3 + 2)  = F(JJLOC3 + 2) - XH
                  XH = XIJDST(3) * DF3
                  FI(IILOC3 + 3) = FI(IILOC3 + 3) + XH
                  F(JJLOC3 + 3)  = F(JJLOC3 + 3) - XH
#ifdef HUM_VIR
                  XH = XIJDST(1) * DF3EL
                  FELECI(IILOC3 + 1) = FELECI(IILOC3 + 1) + XH
                  FELEC(JJLOC3 + 1)  = FELEC(JJLOC3 + 1) - XH
                  XH = XIJDST(2) * DF3EL
                  FELECI(IILOC3 + 2) = FELECI(IILOC3 + 2) + XH
                  FELEC(JJLOC3 + 2)  = FELEC(JJLOC3 + 2) - XH
                  XH = XIJDST(3) * DF3EL
                  FELECI(IILOC3 + 3) = FELECI(IILOC3 + 3) + XH
                  FELEC(JJLOC3 + 3)  = FELEC(JJLOC3 + 3) - XH
#endif
C
C-- THAT SHOULD WORK FOR 4D (PERT,NOPERT) AND >4D (NOPERT)
C
                  DO M=4, NCALCD
                    XH = XIJDST(M) * DF4TH
                    FI(IILOC3 + M) = FI(IILOC3 + M) + XH
                    F(JJLOC3 + M)  = F(JJLOC3 + M) - XH
                  ENDDO
C
C---- CALCULATE THE VIRIAL IF REQUIRED
C
C
#ifdef ATVIR
C Phil, atomic virial - solute-solute(\intra charge group),solute-solvent
C contribution, use LDOVIR, i.e. always a contribution
                  IF (LDOVIR) THEN                  
                    DO M=1,3
                      XH = XIJDST(M)
                      VIR(M) = VIR(M) - DF3 * XH*XH
                    ENDDO
                  ENDIF
#endif

#ifndef ATVIR
C use LDOIT, i.e. contribution only if in different submolecules
                  IF (LDOIT) THEN
                     DO M=1,3
                        XH = XIJDST(M)
#ifdef HUM_VIR
                        VIR(M) = VIR(M) - XH*DF3LJ*
     $                       (XH -XRI(IILOC3+M)+XR(JJLOC3+M))
#else
                        VIR(M) = VIR(M) - XH*DF3*
     $                       (XH -XRI(IILOC3+M)+XR(JJLOC3+M))
#endif
                     ENDDO
                  ENDIF
#endif
C     
C---- END TEST IF EXCLUDED NEIGHBOUR
C     
                ENDIF
                JJLOC3 = JJLOC3 + NDIM
                JJLOC = JJLOC+1
C     
C------END LOOP OVER JJ ATOMS
C     
 650          CONTINUE
C     
C---- SAVE POINTERS FOR EXCLUDED AND 3RD NEIGHBOURS TO KEEP TRACK
C     OF THE ONES ALREADY SCANNED BY EACH ATOM II OF CHARGE GROUP NI
C     
              KLNE(IILOC) = NXTEX - 1
              ILNE(IILOC) = NENDEX-NXTEX
              KLNE14(IILOC) = NXT14 - 1
              ILNE14(IILOC) = NEND14-NXT14
C     
C---- SAVE POINTERS DONE
C     
              IILOC3 = IILOC3 + NDIM
              IILOC = IILOC+1
C     
C--------END LOOP OVER II ATOMS
C     
 640        CONTINUE
C     
C----------END LOOP OVER NJ SECONDARY CHARGE GROUPS
C     
 600      CONTINUE
C     
C---- PERFORM AN FINAL CHECK IF ALL EXCLUDED AND THIRD NEIGHBOURS
C     HAVE BEEN SCANNED WITHIN THE PAIR-LIST - THIS REALLY COSTS
C     NOTHING
C     
          DO IM=1,INUMAT
            NXTEX = KLNE(IM) + 1
            NENDEX = NXTEX + ILNE(IM)
            NXT14 = KLNE14(IM) + 1
            NEND14 = NXT14 + ILNE14(IM)
            IF ( NXTEX.LT.NENDEX ) THEN
              LEMIS = .TRUE.
            ENDIF
            IF ( NXT14.LT.NEND14 ) THEN
              L14MIS = .TRUE.
            ENDIF
          ENDDO          
C     
C---- DONE NEIGHBOURS FINAL CHECK
C     ADD FORCE OF NI TO GLOBAL FORCE ARRAY
C     
!$OMP parallel do private(II) shared(ICOPTR)
          DO II=1,INUMA3
            F(ICOPTR+II) = F(ICOPTR+II) + FI(II)
#ifdef HUM_VIR
            FELEC(ICOPTR+II) = FELEC(ICOPTR+II) + FELECI(II)
#endif
          ENDDO
C     
C---- DONE ADD FORCE
C     
          IFTOPA = ILTOPA + 1
          IFCORA = IFTOPA + IATOFF
          ICOPTR = ICOPTR+INUMA3
          NI = NI + 1
 200    CONTINUE
        IATOFF = IATOFF + NRP
C     
C------------END LOOP OVER PROTEIN MOLECULES AND PRIMARY (PROTEIN ONLY)
C     CHARGE GROUPS      
C     
 50   CONTINUE
C     
C------------NOW LOOP OVER PRIMARY (SOLVENT ONLY) CHARGE GROUPS
C     SECONDARY CHARGE GROUP WILL ALSO BE SOLVENT I.E.
C     NO PERTURBATION, NO INTERNAL CHARGE GROUP INTERACTION,
C     NO EXCLUDED AND THIRD NEIGHBOURS
C     
      INUMA3 = NDIM*NRAM      
      IFCORA = NRPT+1
      ICOPTR = NRPT*NDIM
      DO 700 NI=NRAGP+1,NRAGT        
        IAUX = IFCORA-1
C     
        DO II=1,NRAM
          IGRP(II) = ICHAR(IAGRP(IAUX+II))
        ENDDO        
C     
        DO II=1,INUMA3
          XI(II)  = XCOORD(ICOPTR+II)
          FI(II)= 0.0D0
#ifdef HUM_VIR
          FELECI(II) = 0.0D0
#endif
        ENDDO
C     
        IF (LDOVIR) THEN
          DO II=1,INUMA3
            XRI(II) = XR(ICOPTR+II)
          ENDDO
        ENDIF
        NUMIJ = INB(NI)
C     
C--------LOOP OVER SECONDARY CHARGE GROUPS (SOLVENT)
C     
        DO 800 NNJ= 1,NUMIJ
C     
C---- GET SECONDARY CHARGE GROUP FROM PAIR-LIST
C     
          NJPTR = NJPTR + 1
          NJ = JNB(NJPTR)
          JFCORA = NRPT + NRAM*(NJ-NRAGP-1) + 1
          IAUX = JFCORA-1
          DO II=1,NRAM
            JGRP(II) = ICHAR(IAGRP(IAUX+II))
          ENDDO
          JCOPTR = NDIM*IAUX
C     
C--------LOOP OVER ATOMS II OF CHARGE GROUPS NI
C     
          IILOC3= 0
          DO 810 II=1,NRAM
            CGILA = CGS(II)
            ILATYP = IACS(II)
            ILGRP  = IGRP(II)
C     
C------LOOP OVER ATOMS JJ OF CHARGE GROUP NJ
C     
            JJLOC3 = JCOPTR
            DO 820 JJ=1,NRAM
C     
C---- CALC DISTANCE AND CONNECTING VECTOR BETWEEN II AND JJ
C     CORRECT FOR BOUNDARY CONDITIONS IF NECESSARY
C     
!$OMP parallel do private(M) shared(IILOC3,JJLOC3)
              DO M=1,NCALCD
                XIJDST(M) = XI(IILOC3+M) - XCOORD(JJLOC3+M)
              ENDDO
              RIJ2 = 0.0D0
              IF (LVAC) THEN
                DO M=1,NCALCD
                  RIJ2 = RIJ2 + XIJDST(M)**2
                ENDDO
                IF (L4D) THEN
                  RIJ3D2 = RIJ2 - XIJDST(4)**2
                ELSE
                  RIJ3D2 = RIJ2
                ENDIF
              ELSE
                DO M=1,NCALCD
                  DSTTMP = XIJDST(M)
                  IF (DSTTMP .GE. BOXH(M)) THEN
                    DSTTMP = DSTTMP - BOX(M)
                  ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
                    DSTTMP = DSTTMP + BOX(M)
                  ENDIF
                  XIJDST(M) = DSTTMP
                  RIJ2 = RIJ2 + DSTTMP*DSTTMP
                ENDDO
C     
                IF (L4D) THEN
                  RIJ3D2 = RIJ2 - XIJDST(4)**2
                ELSE
                  RIJ3D2 = RIJ2
                ENDIF
C     
                IF (LOCTO) THEN
                  DSTTMP = BOXOQ-DABS(XIJDST(1))-
     $                 DABS(XIJDST(2))-DABS(XIJDST(3))
                  IF (DSTTMP .LT. 0.0D0) THEN
                    RIJ2 = RIJ2 + DSTTMP*BOX(1) 
                    RIJ3D2 = RIJ3D2 + DSTTMP*BOX(1)
                    XIJDST(1) = XIJDST(1) - SIGN(BOXOH,XIJDST(1))
                    XIJDST(2) = XIJDST(2) - SIGN(BOXOH,XIJDST(2))
                    XIJDST(3) = XIJDST(3) - SIGN(BOXOH,XIJDST(3))
                  ENDIF
                ELSEIF (LMONO .AND. LDOTRA) THEN
                  RIJ2 = RIJ2 + COSB2*XIJDST(1)*XIJDST(3)
                  RIJ3D2 = RIJ2
                ENDIF
              ENDIF
C     
C---- DONE DISTANCE CALCULATION
C     NOW ACTUALLY CALCULATE THE INTERACTION BETWEEN II AND JJ
C     
              JLATYP = IACS(JJ)
              JLGRP  = JGRP(JJ)
C     
              INTIJA = MPAC(ILATYP,JLATYP)
              CA12   = C12(INTIJA)
              CA6    =  C6(INTIJA)
              NDXGRP = NRELKP(ILGRP,JLGRP)
C AK FIXME: add check for energy group list.
              IF(NDXGRP.GT.MXNRE2) THEN
                WRITE(6,*) 'NONBML| NDXGRP=',NDXGRP,
     &                     ' > MXNRE2=', MXNRE2
                CALL STOPGM('NONBML','ARRAY SIZE MISMATCH')
              ENDIF

              QIQJA  = CGILA * CGS(JJ)
C     
              RIJIN2 = 1.0D0 / RIJ2
              RIJINV = DSQRT(RIJIN2)
              RIJIN6 = RIJIN2 * RIJIN2 * RIJIN2
C     
#ifndef EWALD
              VELA4D = QIQJA * RIJINV
#endif
              CRA1   = CA12 * RIJIN6
              CRA12  = CRA1 - CA6
              VLJA4D = CRA12 * RIJIN6
C     
#ifndef EWALD
              DFA4D  = RIJIN2*(VELA4D+6.D0*(CRA1+CRA12)*RIJIN6) 
#endif
C     
C Ewald REAL*8 space contribution - solvent-solvent
#include "inc_nonbml3.h"
C
#ifndef EWALD
              ERF(NDXGRP) = ERF(NDXGRP) + QIQJA * RFE * RIJ3D2
              ERC(NDXGRP) = ERC(NDXGRP) + QIQJA * RFC
#endif
C AK FIXME: add check for energy group list.
              IF(NDXGRP.GT.MXNRE2) THEN
                WRITE(6,*) 'NONBML| NDXGRP=',NDXGRP,
     &                     ' > MXNRE2=', MXNRE2
                CALL STOPGM('NONBML','ARRAY SIZE MISMATCH')
              ENDIF
              EEL(NDXGRP) = EEL(NDXGRP) +  VELA4D
              ELJ(NDXGRP) = ELJ(NDXGRP) +  VLJA4D
C     
              DF3   = DFA4D + RFF * QIQJA
              DF4TH = DFA4D
C
              XH = XIJDST(1) * DF3
              FI(IILOC3 + 1) = FI(IILOC3 + 1) + XH
              F(JJLOC3 + 1)  = F(JJLOC3 + 1)  - XH
              XH = XIJDST(2) * DF3
              FI(IILOC3 + 2) = FI(IILOC3 + 2) + XH
              F(JJLOC3 + 2)  = F(JJLOC3 + 2)  - XH
              XH = XIJDST(3) * DF3
              FI(IILOC3 + 3) = FI(IILOC3 + 3) + XH
              F(JJLOC3 + 3)  = F(JJLOC3 + 3)  - XH
#ifdef HUM_VIR
              XH = XIJDST(1) * DF3EL
              FELECI(IILOC3 + 1) = FELECI(IILOC3 + 1) + XH
              FELEC(JJLOC3 + 1)  = FELEC(JJLOC3 + 1)  - XH
              XH = XIJDST(2) * DF3EL
              FELECI(IILOC3 + 2) = FELECI(IILOC3 + 2) + XH
              FELEC(JJLOC3 + 2)  = FELEC(JJLOC3 + 2)  - XH
              XH = XIJDST(3) * DF3EL
              FELECI(IILOC3 + 3) = FELECI(IILOC3 + 3) + XH
              FELEC(JJLOC3 + 3)  = FELEC(JJLOC3 + 3)  - XH            
#endif
C
C-- THAT SHOULD WORK FOR 4D (PERT,NOPERT) AND >4D (NOPERT)
C
              DO M = 4,NCALCD
                XH = XIJDST(M)*DF4TH
                FI(IILOC3+M) = FI(IILOC3+M) + XH
                F(JJLOC3+M) = F(JJLOC3+M) - XH
              ENDDO
C     
C---- CALCULATE THE VIRIAL IF REQUIRED
C
#ifdef ATVIR
C Phil, atomic virial - solvent - solvent
              IF (LDOVIR) THEN
                DO M=1,3
                  XH = XIJDST(M)                  
                  VIR(M) = VIR(M) - DF3 * XH*XH
                ENDDO
              ENDIF
#endif

#ifndef ATVIR
C molecular virial
              IF (LDOVIR) THEN
                DO M=1,3
                  XH = XIJDST(M)                  
#ifdef HUM_VIR
                  VIR(M) = VIR(M) - XH*DF3LJ*
     $                 (XH -XRI(IILOC3+M)+XR(JJLOC3+M))
#else
                  VIR(M) = VIR(M) - XH*DF3*
     $                 (XH -XRI(IILOC3+M)+XR(JJLOC3+M))
#endif
                ENDDO
              ENDIF
#endif
C     
C---- DONE WITH THE VIRIAL
C     
              JJLOC3 = JJLOC3+NDIM
C     
C------END LOOP OVER JJ ATOMS
C     
 820        CONTINUE
            IILOC3 = IILOC3+NDIM
C     
C--------END LOOP OVER II ATOMS
C     
 810      CONTINUE
C     
C----------END LOOP OVER NJ SECONDARY CHARGE GROUPS
C     
 800    CONTINUE
C     
C---- ADD FORCE OF NI TO GLOBAL FORCE ARRAY
C     
!$OMP parallel do private(II) shared(ICOPTR)
        DO II=1,INUMA3
          F(ICOPTR+II) = F(ICOPTR+II) + FI(II)
#ifdef HUM_VIR
          FELEC(ICOPTR+II) = FELEC(ICOPTR+II) + FELECI(II)
#endif
        ENDDO
C     
C---- DONE ADD FORCE
C     
        ICOPTR = ICOPTR + INUMA3
        IFCORA = IFCORA + NRAM
C     
C------------END LOOP OVER PRIMARY SOLVENT CHARGE GROUPS
C     
 700  CONTINUE
C     
C----------CORRECT VIRIAL IF NECESSARY
C     TRANSFORM COORDINATES BACK IF NECESSARY
C     
      IF (LDOVIR) THEN
!$OMP parallel do private(M)
#ifdef _vpp_
!OCL NOALIAS
#endif
        DO M=1,NDIM
          VIR(M) = VIR(M)*0.5D0
        ENDDO
      ENDIF
C     
      IF (LDOTRA) THEN
        CALL TRACO(NATTOT,0,XCOORD,BETA,-1,LEVERY)
        CALL TRACO(NATTOT,0,F,BETA,-1,LEVERY)
      ENDIF
C
C-- FLAG WEIRD CASES (OK FOR PROFEE)
C
      IF (LEMIS .AND. .NOT. LFREE) THEN
!  qmmm only
!       PRINT *,PRGSTR,': WARNING, NOT ALL EXCLUDED NEIGHBOUR ',
!    $       'WERE SCANNED IN THE PAIR-LIST'
      ENDIF
C
      IF (L14MIS .AND. .NOT. LFREE) THEN
!       PRINT *,PRGSTR,': WARNING, NOT ALL THIRD NEIGHBOUR ',
!    $       'WERE SCANNED IN THE PAIR-LIST'
      ENDIF

C     
C----------RETURN DERIVATIVE VARIABLES ACCORDING TO REQUEST
C     
      IF (.NOT.L3D4D) THEN
C     
C--   IN THIS CASE THESE WOULD BE *WRONG*
C     
        EG34EL = 0.0D0
        EG34RF = 0.0D0
        EG34RC = 0.0D0
        EG34LJ = 0.0D0
      ENDIF
      IF (.NOT.LDOPER) THEN
C     
C--   IN THIS CASE THESE WOULD BE OK, BUT SINCE BONDED ROUTINES
C     WONT RETURN THEM, RETURNING THEM MIGHT GENERATE CONFUSION... 
C     
        EGLEL = 0.0D0
        EGLRF = 0.0D0
        EGLRC = 0.0D0
        EGLLJ = 0.0D0
      ENDIF
C     
 900  FORMAT (1X,A6,': ATOM',I5,' SOFT-CORE PERTURBED')
 910  FORMAT (1X,'AND C6 .EQ. 0 BUT C12 .NE. 0 WITH STATE ',
     $     A1,' ATOM TYPE',I5)
 915  FORMAT (1X,'AND CS6 .EQ. 0 BUT CS12 .NE. 0 WITH STATE ',
     $     A1,' ATOM TYPE',I5)
C     
C------------END NONBML
C     
C     ==--------------------------------------------------------------==
      RETURN
      END
      
#endif
#ifdef __VECTOR 
C
COMMSUBR NONBML CCCCCC I.G. TIRONI, ZUERICH, JUNE 1996 CCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     SUBROUTINE NONBML(NATTOT,NPM,NRAGT,NCALCD,XCOORD,F,XR,
C     $     NSPM,NSP,VIR,
C     $     INB,NSZJNB,JNB,
C     $     MAXNRE,NRELKP,IAGRP,
C     $     RCRF,EPSRF,APPAK,
C     $     EGLEL,EGLRF,EGLRC,EGLLJ,
C     $     EG34EL,EG34RF,EG34RC,EG34LJ,
C     $     NUNRE2,EEL,ERF,ERC,ELJ,
C     $     LDOPER,L3D4D,RLAM,RMU,ALPHLJ,ALPHC,NLAM,MMU,
C     $     LPION,LEVERY,LFREE)
C
C     This is a special version suited for vectorprocessor machines
C     like a CRAY. It is derived from the old cray vector routine. 
C     Path-integral and free energy problems do not run as efficient
C     as a normal MD calculation.
C     Ilario Tironi June 1996.
C
C     
C     NONBML will supply the non-bonded energies and forces for a 
C     12-6 (LJ) and a Coulomb potential, corrected by a reaction-field
C     contribution, for a number of solute (protein) molecules plus solvent 
C     molecules.
C     
C     The potential reads:
C     V(r_ij) = - C6/(r_ij**6) + C12/(r_ij**12)
C     + q_i*q_j/r_ij - (0.5*C1/RCRF**3)*q_i*q_j*(r_ij**2)
C     - ((1-0.5*C1)/RCRF)*q_i*q_j
C     
C     where RCRF is the radius of the reaction field boundary (in principle
C     RCUTL) and 
C     
C     C1 = ( (2-2*EPSRF)*(1+APPAK*RCRF) - EPSRF*(APPAK*RCRF)**2 ) /
C     ( (1+2*EPSRF)*(1+APPAK*RCRF) + EPSRF*(APPAK*RCRF)**2 )
C     
C     is the coefficient for the reaction-field correction beyond RCRF
C     (with EPSRF.GE.1.0 and APPAK.GE.0.0).
C     If EPSRF.LE.0.0, C1 is set to -1, corresponding to an infinite
C     continuum permittivity ( a warning will be printed at first call )
C     
C     If LDOPER.OR.L3D4D (NTG.NE.L<NTGOFF>), a soft-core type, 
C     L<LAMDA> and L<MU> dependent perturbation potential is applied 
C     ( see manual )
C     It is evaluated for atom pairs involving the atoms II for which
C     IPERT(II) .NE. L<NOPERT> .
C     The alpha coefficients for the soft-core perturbation are set to zero
C     except for atom pairs involving the atoms II for which
C     ISCLJ(II) .NE. 0  in which case alpha(LJ) = L<ALPHLJ>
C     and/or
C     ISCCB(II) .NE. 0  in which case alpha(Coulomb/rx field) = L<ALPHC>
C     
C     When LDOPER.AND..NOT.L3D4D (NTG.EQ.L<NTGLAM>), solely the lambda 
C     derivative of the potential will be returned. 
C     When .NOT.LDOPER.AND.L3D4D (NTG.EQ.L<NTGMU>), solely the mu 
C     derivative of the potential will be returned. 
C     When LDOPER.AND.L3D4D (NTG.EQ.L<NTGBOT>), both will be returned. 
C     Derivatives not explicitely required will be returned as zero.
C     
C     NCALCD  LDOPER/L3D4D       calculated
C     
C     3  .F./.F.   Normal interactions in 3D
C     3  .T./.F.   Normal interactions in 3D
C     For pairs involving perturbed atoms:
C     Soft-core lambda dependent potential in 3D,
C     mu has to be zero (if not, set locally to zero)
C     3  .F./.T.   - forbidden -
C     3  .T./.T.   - forbidden -
C     4  .F./.F.   Normal interactions in 4D
C     4  .T./.F.   Normal interactions in 4D
C     For pairs involving perturbed atoms:
C     Soft-core lambda dependent potential,
C     3D-4D mixed at given rmu, lam deriv. returned.
C     4  .F./.T.   Normal interactions in 4D
C     Soft-core lambda dependent potential,
C     3D-4D mixed at given rmu, mu deriv. returned.
C     4  .T./.T.   Normal interactions in 4D
C     Soft-core lambda dependent potential,
C     3D-4D mixed at given rmu, both deriv. returned.
C     >4  .F./.F.   In principle allowed
C     
C     IMPORTANT NOTE:
C     
C     With NCALCD.EQ.4, no check is done whether an 
C     atom is in 4D or not according to 4D specification list (C4D array). 
C     This means that all atoms (solute/solvent) which are to be handled 
C     in 3D MUST HAVE ZERO 4thD COORDINATES !
C     
C     This, however, won't prevent them from having a 4thD component
C     in the force acting on them, due to interaction with 4D atoms.
C     
C     If ALPHLJ.NE.0.0 and for an atom II with ISCLJ(II).NE.0, the algorithm
C     assumes that (either in state A or B), if C6.EQ.0.0, then C12.EQ.0.0.
C     If this is not the case, a non soft-core perturbation of C12 is done.
C     A warning will be printed at first call in this case.
C     
C     The forces are calculated in NCALCD dimensions (L<NCALCD>.LE.L<NDIM>).
C     The reaction-field correction is, however, always calculated in 3D.
C     If L<NCALCD>.EQ.3, whatever L<MU>, a local value of zero will be used.
C     A warning will be printed at first call in this case. Additionally,
C     L3D4D (NTG.EQ.L<NTGMU>.OR.NTG.EQ.L<NTGBOT>) will not be 
C     allowed in this case.    
C     
C     Periodic boundary conditions can be applied, depending on L<NTB>,
C     where the periodic box can be a truncated octahedron, rectangular
C     or monoclinic. 
C     
C     If the virial is to be calculated (ABS(NTB) = L<NTBVIR>), it must
C     be specified whether the solute consists of separate submolecules.
C     This can be done by setting the NSPM > 1 and setting the corresponding
C     values into the NSP array. The virial is only calculated
C     in 3D.
C     
C     The solute atom charge groups and solvent molecule pairs are to
C     be specified. The pair-list must be ordered ( NI<NJ ).
C     
C     The algorithm used assumes that the exclusions of atom II, 
C     JSNE(KNE(II)+1)...JSNE(KNE(II)+INE(II)), are ordered monotonously
C     in sequence (no inversions).
C     The same is assumed for the third neighbour list.
C     An error will occur at first call if the ordering is incorrect.
C     A warning will be printed if all excluded/third neighbours of a 
C     given atom are not scanned within the pair-list (which shouldn't
C     happen if RCUTP is reasonably large). The exception is if LFREE
C     is .TRUE. - in which no warning will be written.
C     
C     NATTOT   total number of atoms
C     NPM      number of (identical) solute molecules
C     NRAGT    total number of charge groups
C     NCALCD   the number of dimensions in which to calc forces.
C     XCOORD(NDIM*NATTOT)
C     atom cartesian coordinates.
C     F(NDIM*NATTOT)
C     delivered with the forces calculated added.
C     XR(NDIM*NATTOT)
C     atom cartesian coordinates relative to the
C     submolecular centres of mass.
C     This array is only used if the virial is to be
C     calculated, i.e. if (ABS(NTB)) .EQ. L<NTBVIR>).
C     It is still accessed (read only) in the opposite case 
C     and thus the dummy array used in the call should have 
C     the correct dimension (e.g. XCOORD).
C     NSPM     If the virial is to calculated, then L<NSPM> and L<NSP> are
C     accessed in order to calculate the molecular virial.
C     NSPM is the number of separate (sub)molecules forming 
C     one solute molecule
C     NSP(NSPM+1)
C     atom sequence numbers of the last atoms of the
C     submolecules.
C     Only elements 1..NSPM are considered to be filled
C     with data. We need NSPM+1 elements in order to avoid
C     a double IF statement when checking for submolecules
C     VIR(NDIM)
C     delivered with the components of the virial
C     (ABS(NTB)) .EQ. L<NTBVIR>) - only in 3D !
C     INB(NRAGT)
C     pointerlist of the charge group pairlist.
C     INB(I) specifies the number of charge groups
C     or solvent molecules J with J>I forming a pair
C     with solute charge group or solvent molecule I.
C     NSZJNB   the number of elements in the pairlist
C     JNB(NSZJNB)
C     charge group pairlist. It specifies solute charge
C     groups or solvent molecules J forming a pair with
C     solute charge group or solvent molecule I in 
C     natural order; for each I the pointer in JNB is to
C     incremented in INB(I)
C     MAXNRE   number of energy groups
C     NRELKP(MAXNRE,MAXNRE)   
C     index of pairs of energy groups in energy arrays
C     IAGRP(NATTOT) a CHARACTER array which is used to specify
C     the energy groups of the atoms involved
C     RCRF     radius of the reaction field boundary (in principle
C     RCUTL) and 
C     EPSRF    dielectric permittivity of the reaction field 
C     continuum outside the boundary
C     APPAK    inverse Debye screening length
C     EGLEL    delivered with the value of dV/d(L<LAMBDA>) at fix L<MU>
C     for the Coulomb term if LDOPER = .TRUE.
C     EGLRF    delivered with the value of dV/d(L<LAMBDA>) at fix L<MU>
C     for the distance dependent Poisson-Boltzmann reaction
C     field term if LDOPER = .TRUE.
C     EGLRC    delivered with the value of dV/d(L<LAMBDA>) at fix L<MU>
C     for the distance independent Poisson-Boltzmann reaction
C     field term if LDOPER = .TRUE.
C     EGLLJ    delivered with the value of dV/d(L<LAMBDA>) at fix L<MU>
C     for the Lennard-Jones term if LDOPER = .TRUE.
C     EG34EL   delivered with the 3D to 4D free energy contribution
C     for the Coulomb term if LDOPER = .TRUE. .AND. NCALCD = 4
C     EG34RF   delivered with the 3D to 4D free energy contribution
C     for the distant dependent Poisson-Boltzmann reaction field
C     term if LDOPER = .TRUE. .AND. NCALCD = 4
C     EG34RC   delivered with the 3D to 4D free energy contribution
C     for the distant independent Poisson-Boltzmann reaction field
C     term if LDOPER = .TRUE. .AND. NCALCD = 4
C     EG34LJ   delivered with the 3D to 4D free energy contribution
C     for the Lennard-Jones term if LDOPER = .TRUE. .AND. NCALCD = 4
C     NUNRE2   number of energy groups
C     EEL(NUNRE2)
C     delivered with the Coulomb energy components as
C     defined by L<NRE>
C     ERF(NUNRE2)
C     delivered with the distant dependent
C     Poisson-Boltzmann energy components as defined by L<NRE>
C     ERC(NUNRE2)
C     delivered with the distant independent
C     Poisson-Boltzmann energy components as defined by L<NRE>
C     ELJ(NUNRE2)
C     delivered with the Lennard-Jones energy components
C     as defined by L<NRE>
C     LDOPER   perturbation is applied, lambda derivative computed
C     LDO34D   perturbation is applied, mu derivative computed
C     RLAM     actual value of the lambda parameter
C     RMU      actual value of the mu parameter
C     ALPHLJ   soft-core alpha parameter for the Lennard-Jones component
C     ALPHC    soft-core alpha parameter for the electrostatics component
C     NLAM     exponent of lambda in the Hamiltonian coupling scheme
C     MMU      exponent of mu in the Hamiltonian coupling scheme
C     LPION    a path-integral calculation is performed 
C     LEVERY
C     .FALSE.  a number of quantities depending on BOX, BETA etc.
C     are only calculated at the first SUBR. call.
C     .TRUE.   they are calculated every SUBR. call.
C     LFREE
C     .TRUE.   the pair list and interaction is calculated solely 
C     for perturbed atoms. This is used in free energy 
C     extrapolation.
C     
C     NOT passed as arguments, but accessed in common blocks are
C     the following:
C     
C     If L<LDOPER> .EQ. .TRUE., then L<RLAM>, L<RMU>, L<IPERT>, L<CGB>
C     and L<ISCC>
C     are accessed in order to calculate a perturbation dependent
C     potential.
C     IPERT(MAXNRP)
C     If IPERT(II) =0, the charge of atom II is not perturbed.
C     otherwise it denotes where the new charge can be found
C     in L<CGB>.
C     CGB(MAXPAT)
C     charges determining the perturbation potential for
C     state B.
COMMEND
C     ==--------------------------------------------------------------==      
      SUBROUTINE NONBML(NATTOT,NPM,NRAGT,NCALCD,XCOORD,F,XR,
C     ==--------------------------------------------------------------==      
#ifdef EWALD
     $     ALPHA,
#endif
#ifdef EWATCUT
     $     RCUTP,
#endif  
     $     NSPM,NSP,VIR,
     $     INB,NSZJNB,JNB,
     $     MAXNRE,NRELKP,IAGRP,
     $     RCRF,EPSRF,APPAK,
     $     EGLEL,EGLRF,EGLRC,EGLLJ,
     $     EG34EL,EG34RF,EG34RC,EG34LJ,
     $     NUNRE2,EEL,ERF,ERC,ELJ,
     $     LDOPER,L3D4D,RLAM,RMU,ALPHLJ,ALPHC,NLAM,MMU,
     $     LPION,LEVERY,LFREE)
      
C     
C     INCLUDES
C     
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'pertsz.h'
      INCLUDE 'pertar.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'vector.h'
      INCLUDE '../mm_input.inc'
C     
C     ARGUMENTS
C     
      LOGICAL LEVERY,LDOPER,L3D4D,LPION,LFREE
      INTEGER NATTOT,NPM,NRAGT,NCALCD
      INTEGER NSZJNB,INB(NRAGT),JNB(NSZJNB)
      INTEGER NLAM,MMU
      INTEGER MAXNRE,NRELKP(MAXNRE,*)
      INTEGER NSPM,NSP(NSPM+1)
      INTEGER NUNRE2
      real*8 XCOORD(*),F(*),XR(*)
      real*8 VIR(NDIM)
      real*8 RCRF,EPSRF,APPAK
      real*8 RLAM,RMU,ALPHLJ,ALPHC
      real*8 EGLEL,EGLRF,EGLRC,EGLLJ
      real*8 EG34EL,EG34RF,EG34RC,EG34LJ
      real*8 EEL(*),ERF(*),ERC(*),ELJ(*)
      CHARACTER IAGRP(NATTOT)
C     
C     LOCAL PARAMS
C     
      real*8 ONE
      PARAMETER (ONE = 1.0D0)
      real*8 ZEROC6, ZEROAL
      PARAMETER (ZEROC6 = 1.0D-20 , ZEROAL = 1.0D-20)
C     
C     SAVED VARIABLES      
C     
      LOGICAL LMONO,LOCTO,LVAC,LDOVIR,LDOTRA,LFIRST,L4D
      SAVE    LMONO,LOCTO,LVAC,LDOVIR,LDOTRA,LFIRST,L4D
      LOGICAL LPIDOP
      SAVE    LPIDOP
      real*8  PININV
      SAVE    PININV
      real*8 COSB,COSB2,BOXOH,BOXOQ
      SAVE COSB,COSB2,BOXOH,BOXOQ
      real*8 RFF, RFE, RFC, C1, RCRF2
      SAVE RFF, RFE, RFC, C1, RCRF2
C     
C     LOCAL VARIABLES
C     
      LOGICAL LIZERO, LJZERO, LPERTL, LIPERT

      INTEGER NRAGP,NI,NRPT
      real*8 BETAR

      CHARACTER PRGSTR*(6)
C     
C     PERTURBATION VARIABLES
C     
      real*8 RIJ3D2, VELA4D, VLJA4D, VRFA4D, DFA4D, DFRFA4,
     $     DLLJA4, DLELA4, DLRFA4, VELB4D,
     $     VLJB4D, VRFB4D, DFB4D, DFRFB4, DLLJB4, DLELB4,
     $     DLRFB4,
     $     VELA3D, VLJA3D, VRFA3D, DFA3D, DFRFA3, DLLJA3,
     $     DLELA3, DLRFA3, VELB3D, VLJB3D, VRFB3D ,
     $     DFB3D, DFRFB3, DLLJB3, DLELB3, DLRFB3,
     $     VEL4D, VRF4D, VLJ4D, DF4D4D, DF4D, DLLJ4D, DLEL4D,
     $     DLRF4D, 
     $     VEL3D, VRF3D, VLJ3D, DF3D, DLLJ3D, DLEL3D, DLRF3D,
     $     DDM3EL, DDM3RF, DDM3LJ,
     $     VRCA,VRCB,VRC,DLRC,
     $     QRFR,
     $     SOFTLJ, SOFTEL, CA126, CB126, SFEL, SFLJ,
     $     RSOFT, RSOFT2, RIJ4, RIJ6, RRF, RRF2, RFSQRT,
     $     RLMA, RL2A, RLMB, RL2B, RLLA, RLLB, RLB, RLA, RDLB, RDLA,
     $     RLARMU, RLBRMU,
     $     RM4D, RM3D, RDM3D, RDM4D,
     $     QIQJA, QIQJB ,
     $     RIJ2,RIJINV,RIJIN2,RIJIN6,
     $     CA12,CA6,CB12,CB6,
     $     CRA1,CRA12, RMULOC
      
      LOGICAL LERROR, LVIR, LIEX14
      
      INTEGER MAXLEN
      SAVE MAXLEN
      
      INTEGER ISRCHEQ, ISRCHFLE
      
      INTEGER JTOP, JRC, IPQ,
     $     I, J, JRCNUM, JRCL1X, NJNB, II,
     $     I1X, K1X, ILAST, ICE, IVIRL, ITMPA, ITMPB,
     $     JRCTOT, JRCSRT, JRCLEN, JRC1, JRC2,
     $     NRAG, NPAIR, K, JLEN, NR, JTYPE, INTIJA, INTIJB, M,
     $     NITOT, NN, IPC, JRCLOC, JRCEXL,
     $     IX, JX, NIMAX, ILEX, IPCB,
     $     IBASE, J1X, ILJ, I9X, NJIND, NJINDB, JN, I20X, I44X, I32X,
     $     JLOW, JLEN2, ILEN, KFORCE
      
      real*8 BOX1, BOX2, BOX3, BOX4, BOXH1, BOXH2, BOXH3, BOXH4,
     $     BOXIV1, BOXIV2, BOXIV3, BOXIV4, XIJ, YIJ, ZIJ,
     $     DX, DY, DZ, DW, DRIJ, DF, 
     $     XH, YH, ZH, WH, PIFAC, DF3T, DF4T, WIJ

      SAVE NR, NRAG, NRAGP, NRPT, KFORCE
      
C     
C     DATA
C     
      DATA LFIRST /.TRUE./
      DATA PRGSTR /'NONBML'/

#include "inc_nonbml2.h"
#ifndef EWALD
      CALL STOPGM('NONBML','vector+reaction field: debug it yourself')
#endif
C     
C------------BEGIN NONBML 

      IF (LFREE) THEN
         WRITE(6,*) 'FREE ENERGY EXTRAPOLATION (LFREE .EQ. TRUE)',
     $        ' NOT IMPLEMENTED IN VECTOR VERSION OF NONBML'
         WRITE(6,*) 'USE THE STANDARD VERSION OF NONBML'
         CALL STOPGM('NONBML','FATAL ERROR')
      ENDIF 

C     
C------------DO FIRST CALL SETUP (SWITCHES,"MONOCLINIC" VARIABLES
C     REACTION FIELD VARIABLES, SIMPLE CONSISTENCY CHECKS)     
C     
      IF (LEVERY .OR. LFIRST) THEN

         LMONO = (NTB .GT. 0)
         LOCTO = (NTB. LT. 0)
         LVAC  = (NTB .EQ. 0)
         LDOVIR= (ABS(NTB) .EQ. NTBVIR)
C     
         LDOTRA = .FALSE.
         IF (LMONO) THEN
            BETAR = BETA*DATAN(ONE)/45.0D0
            COSB = DCOS(BETAR)
            LDOTRA = (DABS(COSB).GE.1.D-04)
            COSB2 = 2.0D0*COSB
         ENDIF
         IF (LOCTO) THEN
            BOXOH = BOXH(1)
            BOXOQ = BOX(1)*0.75D0
         ELSE
            BOXOH = 0.0D0
            BOXOQ = 0.0D0
         ENDIF
      ENDIF

      LPERTL = (LDOPER.OR.L3D4D)

      IF (LFIRST) THEN
         LFIRST = .FALSE.
         L4D = (NCALCD .EQ. 4)
         LPIDOP = .TRUE.
         IF (NPID .NE. 0) THEN
            PININV = ONE / NPID
         ENDIF
C     
C     
C---- PERFORM A FEW FIRST CALL CHECKS
C     
C--   IF LPERTL AND ALPHLJ.NE.0.0 CHECK THAT C6.EQ.0 => C12.EQ.0 FOR ATOMS II 
C     WITH ISCLJ(II).NE.0 AND ANY OTHER ATOM ( STATE A AND B )
C     
         IF (LPERTL .AND. DABS(ALPHLJ).GT.ZEROAL) THEN
            DO 26 II=1,NRP
               IF (IPERT(II) .NE. 0) THEN
                  IF (ISCLJ(IPERT(II)) .NE. 0) THEN
                     ITMPA = IAC(II)
                     ITMPB = IACB(IPERT(II))
                     DO JTYPE=1,NRATT
                        INTIJA = MPAC(ITMPA,JTYPE)
                        INTIJB = MPAC(ITMPB,JTYPE)
                        IF (C6(INTIJA) .LT. ZEROC6.AND.
     $                       C12(INTIJA) .GT. ZEROC6) THEN
                           PRINT 900,PRGSTR,II
                           PRINT 910,'A',JTYPE
                        ENDIF
                        IF (C6(INTIJB) .LT. ZEROC6.AND.
     $                       C12(INTIJB) .GT. ZEROC6) THEN
                           PRINT 900,PRGSTR,II
                           PRINT 910,'B',JTYPE
                        ENDIF
                        IF (CS6(INTIJA) .LT. ZEROC6.AND.
     $                       CS12(INTIJA) .GT. ZEROC6) THEN
                           PRINT 900,PRGSTR,II
                           PRINT 915,'A',JTYPE
                        ENDIF
                        IF (CS6(INTIJB) .LT. ZEROC6.AND.
     $                       CS12(INTIJB) .GT. ZEROC6) THEN
                           PRINT 900,PRGSTR,II
                           PRINT 915,'B',JTYPE
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
 26         CONTINUE
         ENDIF

C     
C--   CALCULATING IN 4D WITH 3D COORDINATES IS WRONG
C     
            IF ( NCALCD.GT.NDIM ) THEN
               PRINT *,PRGSTR,' ERROR: NCALCD.GT.NDIM'
               CALL STOPGM('NONBML','FATAL ERROR')
            ENDIF
         
C     
C---- DONE FIRST CALL CHECKS
C     NOW REACTION FIELD VARIABLES
C     
         IF (EPSRF .LE. 0.0D0) THEN
            C1 = -1.0D0
            PRINT *,PRGSTR, ' : WARNING, EPSRF.LE.0 , PERMITTIVITY ',
     $           'SET TO INFINITY'
         ELSE
            C1 = ((2.0D0-2.0D0*EPSRF)*(1.0D0+APPAK*RCRF) 
     $           - EPSRF*(APPAK*RCRF)**2)
     $           /((1.0D0+2.0D0*EPSRF)*(1.0D0+APPAK*RCRF) 
     $           + EPSRF*(APPAK*RCRF)**2)
         ENDIF
         RFF = C1 / RCRF**3
         RFE = -0.5D0 * RFF
#ifdef EWALD
         RFF = 0.0D0
         RFE = 0.0D0
#endif
         RFC = -(1.0D0-0.5D0*C1)/RCRF
         RCRF2 = RCRF*RCRF
         
         NRAG = NPM*NCAG
         NRAGP = NPM*NCAG
         NRPT  = NPM*NRP
C     
         KFORCE = (NRATT*(NRATT+1))/2
      
         NR = NATTOT
   
      ENDIF
C     
C------------DONE FIRST CALL SETUP
C     CALCULATE RMU/RLAMBDA DEPENDANT VARIABLES AND A FEW MORE CHECKS
C
cmb-> Now fix some defaults and initialize variables
      LIZERO=.FALSE.
      LJZERO=.FALSE.
C      LPERTL=.FALSE.   ! this is a bug. MB. LPERTL is defined elsewhere...
      LIPERT=.FALSE.
      RL2A = 0.0D0
      RLMA = 0.0D0
      RL2B = 0.0D0
      RLMB = 0.0D0 
      RLB  = 0.0D0
      RDLB = 0.0D0
      RLA  = 0.0D0
      RDLA = 0.0D0
      RLLB = 0.0D0
      RLLA = 0.0D0
      RLARMU = 0.0D0
      RLBRMU = 0.0D0
      RM4D  = 0.0D0
      RDM4D = 0.0D0
      RM3D  = 0.0D0
      RDM3D = 0.0D0
cmb-> end defaults

      IF (LPERTL) THEN  ! ...and this statement would be always false
C     
C--   CALCULATING IN 3D DOESN'T ALLOW THE CALCULATION OF MU DERIVATIVE
C     
         IF (.NOT. L4D .AND. L3D4D) THEN
            PRINT *,PRGSTR,': NCALCD.EQ.3.AND.',
     $           '(NTG.EQ.NTGMU.OR.NTG.EQ.NTGBOT) NOT ALLOWED'
            CALL STOPGM('NONBML','FATAL ERROR')
         ENDIF
C     
C--   CALCULATING IN 3D WITH RMU.NE.0 SHOULDN'T BE DONE
C     
         IF (.NOT. L4D .AND. RMU .NE. 0.0D0) THEN
            PRINT *,PRGSTR,' : WARNING, NCALCD.EQ.3.AND.RMU.NE.0.0'
            PRINT *,'         LOCAL RMU WILL BE SET TO ZERO'
            RMULOC = 0.0D0
         ELSE
            RMULOC = RMU
         ENDIF
C     
C--   LOOKS LIKE SOME MACHINES DON'T LIKE 0**0 = 1
C     ... WHICH WE UNDERSTAND ...
C     
         RL2A = RLAM**2
         RLMA = RL2A + RMULOC**2
         RL2B = (1.0D0 - RLAM)**2
         RLMB = RL2B + RMULOC**2
C     
         IF (RLAM.EQ.0.0D0 .AND. NLAM.EQ.1) THEN
            RLB = 0.0D0
            RDLB = DBLE(NLAM)
         ELSE
            RLB = RLAM**NLAM
            RDLB = DBLE(NLAM)*RLB/RLAM
cmb            RDLB = RLAM**(NLAM - 1) * NLAM
         ENDIF
         IF (RLAM.EQ.1.0D0 .AND. NLAM.EQ.1) THEN
            RLA = 0.0D0
            RDLA = -DBLE(NLAM)
         ELSE
            RLA = (1.0D0 - RLAM)**NLAM          
            RDLA = -DBLE(NLAM)*RLA/(1.0D0-RLAM)
cmb            RDLA = - (1.0D0 - RLAM)**(NLAM - 1) * NLAM
         ENDIF
C     
         RLLB = -(1.0D0 - RLAM) * RLB
         RLLA = RLAM * RLA
         RLARMU = RLA * RMULOC
         RLBRMU = RLB * RMULOC
C     
         IF (RMULOC.EQ.0.0D0 .AND. MMU.EQ.1) THEN
            RM4D = 0.0D0
            RDM4D = DBLE(MMU)
         ELSE
            RM4D = RMULOC**MMU
            RDM4D = DBLE(MMU)*RM4D/RMULOC
cmb            RDM4D = RMULOC**(MMU - 1) * MMU
         ENDIF
         IF (RMULOC.EQ.1.0D0 .AND. MMU.EQ.1) THEN
            RM3D = 0.0D0
            RDM3D = -DBLE(MMU)
         ELSE
            RM3D = (1.0D0 - RMULOC)**MMU
            RDM3D = -DBLE(MMU)*RM3D/(1.0D0 - RMULOC)
cmb            RDM3D = - (1.0D0 - RMULOC)**(MMU - 1) * MMU
         ENDIF
C     
         IF (.NOT. L4D) THEN
            RDM3D = 0.0D0
         ENDIF     
C     
      ENDIF

C
C------------DONE RMU/RLAMBDA DEPENDANT VARIABLES

C     
C------------INITIALIZE ARRAYS AND GENERAL VARIABLES
C     

      BOX1 = BOX(1)
      BOX2 = BOX(2)
      BOX3 = BOX(3)
      BOX4 = BOX(4)
      IF (NCALCD .EQ. 3) THEN
         BOX4 = 1.0D10
      ENDIF
      IF (LVAC) THEN
         BOX1 = 1.0D10
         BOX2 = 1.0D10
         BOX3 = 1.0D10
         BOX4 = 1.0D10
      ENDIF
      IF (.NOT. LDOTRA) THEN
         COSB2 = 0.0D0
      ENDIF
      
      BOXH1 = 0.5D0 * BOX1
      BOXH2 = 0.5D0 * BOX2
      BOXH3 = 0.5D0 * BOX3
      BOXH4 = 0.5D0 * BOX4 

      BOXIV1 = 1.0D0 / BOX1
      BOXIV2 = 1.0D0 / BOX2
      BOXIV3 = 1.0D0 / BOX3
      BOXIV4 = 1.0D0 / BOX4

      BOXOH = BOXH1
      BOXOQ = BOX1 * 0.75D0
      
      EGLEL = 0.0D0
      EGLRF = 0.0D0
      EGLRC = 0.0D0
      EGLLJ = 0.0D0
C     
      EG34EL = 0.0D0
      EG34RF = 0.0D0
      EG34RC = 0.0D0
      EG34LJ = 0.0D0
C     
!$OMP parallel do private(M)
#ifdef _vpp_      
!OCL NOALIAS 
#endif
      DO M=1,NUNRE2
         EEL(M) = 0.0D0
#ifndef EWALD
         ERF(M) = 0.0D0
         ERC(M) = 0.0D0
#endif
         ELJ(M) = 0.0D0
      ENDDO
C
      CALL mm_AZZERO(VIR,NDIM)
      CALL mm_AZZERO(TVIR,3*MAXP)
      
!$OMP parallel do private(I,J) __COLLAPSE2
#ifdef _vpp_
!OCL NOALIAS 
#endif
      DO J = 1, IVECT
         DO I = 1, IEPAIR
            TEEL(I,J) = 0.0D0
            TERF(I,J) = 0.0D0
            TERC(I,J) = 0.0D0
            TELJ(I,J) = 0.0D0
         ENDDO
      ENDDO
      
C     FILL VECTOR ARRAYS
      IF (NDIM .EQ. 3) THEN
         IF (LDOVIR) THEN
!$OMP parallel do private(I)
#ifdef _vpp_
!OCL NOALIAS
#endif
            DO I = 1, NR
               XC(I) = XCOORD((I - 1) * 3 + 1)
               YC(I) = XCOORD((I - 1) * 3 + 2)
               ZC(I) = XCOORD((I - 1) * 3 + 3)
               FLRX(I) = 0.0D0
               FLRY(I) = 0.0D0
               FLRZ(I) = 0.0D0
               XRC(I) = XR((I - 1) * 3 + 1)
               YRC(I) = XR((I - 1) * 3 + 2)
               ZRC(I) = XR((I - 1) * 3 + 3)
            ENDDO
         ELSE
!$OMP parallel do private(I)
#ifdef _vpp_
!OCL NOALIAS 
#endif
            DO I = 1, NR
               XC(I) = XCOORD((I - 1) * 3 + 1)
               YC(I) = XCOORD((I - 1) * 3 + 2)
               ZC(I) = XCOORD((I - 1) * 3 + 3)
               FLRX(I) = 0.0D0
               FLRY(I) = 0.0D0
               FLRZ(I) = 0.0D0
            ENDDO
         ENDIF
      ELSEIF (NDIM .EQ. 4 .AND. NCALCD .EQ. 3) THEN
         IF (LDOVIR) THEN
!$OMP parallel do private(I)
#ifdef _vpp_
!OCL NOALIAS
#endif
            DO I = 1, NR
               XC(I) = XCOORD((I - 1) * 4 + 1)
               YC(I) = XCOORD((I - 1) * 4 + 2)
               ZC(I) = XCOORD((I - 1) * 4 + 3)
               FLRX(I) = 0.0D0
               FLRY(I) = 0.0D0
               FLRZ(I) = 0.0D0
               XRC(I) = XR((I - 1) * 4 + 1)
               YRC(I) = XR((I - 1) * 4 + 2)
               ZRC(I) = XR((I - 1) * 4 + 3)
            ENDDO
         ELSE
!$OMP parallel do private(I)
#ifdef _vpp_
!OCL NOALIAS
#endif
            DO I = 1, NR
               XC(I) = XCOORD((I - 1) * 4 + 1)
               YC(I) = XCOORD((I - 1) * 4 + 2)
               ZC(I) = XCOORD((I - 1) * 4 + 3)
               FLRX(I) = 0.0D0
               FLRY(I) = 0.0D0
               FLRZ(I) = 0.0D0
            ENDDO
         ENDIF
      ELSEIF (NDIM .EQ. 4 .AND. NCALCD .EQ. 4) THEN
         IF (LDOVIR) THEN
!$OMP parallel do private(I)
#ifdef _vpp_
!OCL NOALIAS
#endif
            DO I = 1, NR
               XC(I) = XCOORD((I - 1) * 4 + 1)
               YC(I) = XCOORD((I - 1) * 4 + 2)
               ZC(I) = XCOORD((I - 1) * 4 + 3)
               WC(I) = XCOORD((I - 1) * 4 + 4)
               FLRX(I) = 0.0D0
               FLRY(I) = 0.0D0
               FLRZ(I) = 0.0D0
               FLRW(I) = 0.0D0
               XRC(I) = XR((I - 1) * 4 + 1)
               YRC(I) = XR((I - 1) * 4 + 2)
               ZRC(I) = XR((I - 1) * 4 + 3)
            ENDDO
         ELSE
!$OMP parallel do private(I)
#ifdef _vpp_
!OCL NOALIAS
#endif
            DO I = 1, NR
               XC(I) = XCOORD((I - 1) * 4 + 1)
               YC(I) = XCOORD((I - 1) * 4 + 2)
               ZC(I) = XCOORD((I - 1) * 4 + 3)
               WC(I) = XCOORD((I - 1) * 4 + 4)
               FLRX(I) = 0.0D0
               FLRY(I) = 0.0D0
               FLRZ(I) = 0.0D0
               FLRW(I) = 0.0D0
            ENDDO
         ENDIF
      ENDIF
      
C     
C*****CHECK WHETHER THE VECTORLENGTH MAXP WILL BE TOO SHORT.
C     
      IF (LNEWP) THEN
         MAXLEN = 0
         NIMAX = 0
         JTOP = 1
         DO 2300 NI = 1, NRAGT - 1
            NPAIR = INB(NI)
            IF (NPAIR .GT. 0) THEN
               JLEN = 0
!$OMP parallel do private(IPQ) reduction(+:JLEN)
#ifdef _vpp_
!OCL NOALIAS 
#endif
               DO IPQ = 1, NPAIR
                  NJL(IPQ) = INAG(JNB(IPQ+JTOP-1))
                  JLEN = JLEN + NJL(IPQ)
               ENDDO
            ELSE
               JLEN = 0
            ENDIF
C     
            JTOP = JTOP + NPAIR
            IF (JLEN .GT. MAXLEN) THEN
               MAXLEN = JLEN
               NIMAX = NI
            ENDIF
C     
 2300    CONTINUE
C     
         IF (MAXLEN .GT. MAXP) THEN
            PRINT 9200, MAXLEN, MAXP
 9200       FORMAT(/' THE LENGTH OF THE SCRATCH ARRAYS (',I5,') IS ',
     $           ' LARGER THAN THE WORK SPACE PROVIDED (',I5,') ',
     $           'IN NONBML'/)
            LERROR = .TRUE.
            GOTO 5750
         ENDIF
         PRINT 9220, MAXLEN, NIMAX
 9220    FORMAT(/' THE LARGEST NUMBER OF ATOMS DENOTED',
     $        ' BY THE PAIRLIST IS ',I5,', AROUND CHARGE GROUP ',I5/)
C     
         LNEWP = .FALSE.
      ENDIF
C*****TRANSFORM CARTESIAN COORDINATES TO OBLIQUE
C*****CONTRAVARIANT COORDINATES, WHEN REQUIRED
C     
      IF (LDOTRA) CALL TRACO (NR, 0, XCOORD, BETA, 1, LEVERY)
      IF (LDOTRA) CALL TRACO (NR, 0, F, BETA, 1, LEVERY)
      
      IF (NCALCD .EQ. 3) THEN

C*******************************************************************************
C 3D
C*******************************************************************************

C     
C*****CALCULATE THE NON-BONDED INTERACTIONS IN 3D.
C     
C     
C*****TREAT THE INTRA-CHARGE GROUP INTERACTIONS.
C*****IT IS ASSUMED THAT CHARGE GROUPS CONTAINING LESS THAN 4 ATOMS
C*****CONTAIN NO INTRA CHARGE GROUP INTERACTIONS.
C*****SOLVENT MOLECULES ARE ASSUMED TO CONTAIN NO INTRA MOLECULE
C*****INTERACTIONS
C     
         CALL WHENIGT (NCAG, INAG(1),1,3,NILIST,NITOT)
C     
         DO 2700 NN = 1, NITOT
            NI = NILIST(NN)
            ILAST = IFIRST(NI) + INAG(NI) - 1
C     
C     
            DO 2650 I = IFIRST(NI), ILAST - 1
               DO 2600 J = I + 1, ILAST
C     
C*****CHECK WHETHER ATOM J IS EXCLUDED FROM ATOM I SOMEHOW.
C     
                  ILEX = KNE(I) + INE(I)
                  JRCLOC = KNE(I) + 1
                  JRCLEN = ILEX - JRCLOC + 1
                  IF(JRCLEN.GT.0) THEN
                     JRCEXL=ISRCHEQ(JRCLEN,JSNE(JRCLOC),1,J)-JRCLEN
                  ELSE
                     JRCEXL=1
                  ENDIF

                  IF (LPION) THEN
                     LIZERO = (IPICT(I) .EQ. 0)
                     LJZERO = (IPICT(J) .EQ. 0)
                     IF (LIZERO .AND. LJZERO) THEN
                        PIFAC = 1.0D0
                     ELSEIF (LIZERO .OR. LJZERO) THEN
                        PIFAC = PININV
                     ELSEIF (IPICT(I) .EQ. IPICT(J)) THEN
                        PIFAC = PININV
                     ELSE
                        PIFAC = 0.0D0
                        JRCEXL = 0
                     ENDIF
                  ELSE
                     PIFAC = 1.0D0
                  ENDIF 

                  IF(JRCEXL.GT.0) THEN
C     
C*****CHECK 1-4 EXCLUSIONS.
C     
                     ILEX = KNE14(I) + INE14(I)
                     JRCLOC = KNE14(I) + 1
                     JRCLEN = ILEX - JRCLOC + 1
                     IF(JRCLEN.GT.0) THEN
                        LIEX14 = ISRCHEQ(JRCLEN,JSNE14(JRCLOC),1,J)
     $                       .LE. JRCLEN
                     ELSE
                        LIEX14 = .FALSE.
                     ENDIF
C     
                     IPC = MPAC(IAC(I),IAC(J))
                     IF (LIEX14) THEN
                        CA12 = CS12(IPC)
                        CA6 = CS6(IPC)
                     ELSE
                        CA12 = C12(IPC)
                        CA6 = C6(IPC)
                     ENDIF
C     
                     IF (LPERTL) THEN
                        IPCB = MPAC(IACBT(I),IACBT(J))
                        IF (LIEX14) THEN
                           CB12 = CS12(IPCB)
                           CB6 = CS6(IPCB)
                        ELSE
                           CB12 = C12(IPCB)
                           CB6 = C6(IPCB)
                        ENDIF
                        IF (DABS(CA6) .LT. ZEROC6) THEN
                           CA126 = 0.0D0
                        ELSE
                           CA126 = CA12 / CA6
                        ENDIF
                        IF (DABS(CB6) .LT. ZEROC6) THEN
                           CB126 = 0.0D0
                        ELSE
                           CB126 = CB12 / CB6
                        ENDIF
                     ENDIF

C
C*****CALCULATE THE LJ-ENERGY AND FORCE.
C     
                     IBASE = 0
                     IF (LOCTO) THEN
                        DO K = 1, NPM
                           IX = I + IBASE
                           JX = J + IBASE
                           DX = (XC(IX)-XC(JX))*BOXIV1
                           DX = (DX-ANINT(DX))*BOX1
                           DY = (YC(IX)-YC(JX))*BOXIV2
                           DY = (DY-ANINT(DY))*BOX2
                           DZ = (ZC(IX)-ZC(JX))*BOXIV3
                           DZ = (DZ-ANINT(DZ))*BOX3
                           RIJ2 = DX*DX + DY*DY + DZ*(DZ+COSB2*DX)
                           DRIJ = BOXOQ - DABS(DX)-DABS(DY)-DABS(DZ)
                           IF (DRIJ .LE. 0.D0) THEN
                              RIJ2 = RIJ2 + DRIJ*BOX1
                              DX = DX - SIGN(BOXOH,DX)
                              DY = DY - SIGN(BOXOH,DY)
                              DZ = DZ - SIGN(BOXOH,DZ)
                           ENDIF
                           ILJ = NRELKP(ICET(IX),ICET(JX))
                           QIQJA = CGAT(IX) * CGAT(JX)
                           IF (LLPERT(IX) .OR. LLPERT(JX)) THEN
                             QIQJB = CGBT(IX) * CGBT(JX)
                             SOFTLJ = ALPHLJ * 
     $                          (1 - ISCLJT(IX) * ISCLJT(JX))
                             SOFTEL = ALPHC * 
     $                          (1 - ISCCT(IX) * ISCCT(JX))
                           
                             RIJ3D2 = RIJ2
                             RIJ4 = RIJ2 * RIJ2
                             RIJ6 = RIJ4 * RIJ2
C     
C--   3D, STATE A
C     
C     COULOMB
                             SFEL   = SOFTEL * RLMA
                             RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                             RSOFT  = DSQRT(RSOFT2)
                             VELA3D = QIQJA * RSOFT
C     REACTION FIELD
                             RRF2   = 1.0D0 / (SFEL + RCRF2)
                             RFSQRT = DSQRT(RRF2)
                             RRF    = RRF2 * RFSQRT
                             DFRFA3 = QIQJA * C1 * RRF
                             QRFR   = DFRFA3 * RIJ3D2
                             VRFA3D = -0.5D0 * QRFR
C     LENNARD-JONES
                             SFLJ   = SOFTLJ * CA126
                             CRA1   = 1.0D0 / (SFLJ * RLMA + RIJ6)
                             CRA12  = CA12 * CRA1
                             VLJA3D = CRA1 * (CRA12 - CA6)
                             DFA3D  = CRA1 * ( CRA1 * CRA12 + VLJA3D)
C     DV/DL
                             DLLJA3 = -2.0D0 * SFLJ * DFA3D
                             DLELA3 = -VELA3D * SOFTEL * RSOFT2 
                             DLRFA3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                             DFA3D =RSOFT2*VELA3D+6.0D0*RIJ4*DFA3D
C     
C--   3D, STATE B
C     
C     COULOMB
                             SFEL   = SOFTEL * RLMB
                             RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                             RSOFT  = DSQRT(RSOFT2)
                             VELB3D = QIQJB * RSOFT
C     REACTION FIELD
                             RRF2   = 1.0D0 / (SFEL + RCRF2)
                             RFSQRT = DSQRT(RRF2)
                             RRF    = RRF2 * RFSQRT
                             DFRFB3 = QIQJB * C1 * RRF
                             QRFR   = DFRFB3 * RIJ3D2
                             VRFB3D = -0.5D0 * QRFR
C     LENNARD-JONES
                             SFLJ   = SOFTLJ * CB126 
                             CRA1   = 1.0D0 / (SFLJ * RLMB + RIJ6)
                             CRA12  = CB12 * CRA1
                             VLJB3D = CRA1 * (CRA12 - CB6)
                             DFB3D  = CRA1 * ( CRA1 * CRA12 + VLJB3D)
C     DV/DL
                             DLLJB3 = -2.0D0 * SFLJ * DFB3D
                             DLELB3 = -VELB3D * SOFTEL * RSOFT2 
                             DLRFB3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                             DFB3D = RSOFT2*VELB3D+6.0D0*RIJ4*DFB3D

C     
C--   COMBINE 3D, STATE A AND B
C     
                             VEL3D = RLB * VELB3D + RLA * VELA3D
                             VRF3D = RLB * VRFB3D + RLA * VRFA3D
                             VLJ3D = RLB * VLJB3D + RLA * VLJA3D
                             DF3D  = RLB * DFB3D  + RLA * DFA3D
     $                             + RLB * DFRFB3 + RLA * DFRFA3
                             DLLJ3D = RDLB * VLJB3D + RLLB * DLLJB3
     $                            + RDLA * VLJA3D + RLLA * DLLJA3
                             DLEL3D = RDLB * VELB3D + RLLB * DLELB3 
     $                            + RDLA * VELA3D + RLLA * DLELA3
                             DLRF3D = RDLB * VRFB3D + RLLB * DLRFB3
     $                            + RDLA * VRFA3D + RLLA * DLRFA3

                             VRCA = RFC*QIQJA
                             VRCB = RFC*QIQJB
                             VRC  = RLB*VRCB  + RLA*VRCA
                             DLRC = RDLB*VRCB + RDLA*VRCA
C     
                             EEL(ILJ) = EEL(ILJ) + VEL3D
#ifndef EWALD
                             ERF(ILJ) = ERF(ILJ) + VRF3D
                             ERC(ILJ) = ERC(ILJ) + VRC
#endif
                             ELJ(ILJ) = ELJ(ILJ) + VLJ3D
                           
                             EGLEL = EGLEL + DLEL3D
                             EGLRF = EGLRF + DLRF3D
                             EGLRC = EGLRC + DLRC 
                             EGLLJ = EGLLJ + DLLJ3D
C     
                           ELSE
                             RIJIN2 = 1.0D0 / RIJ2
                             RIJINV = DSQRT(RIJIN2)
                             RIJIN6 = RIJIN2 * RIJIN2 * RIJIN2
C     
                             VELA3D = QIQJA * RIJINV
C     
                             CRA1   = CA12 * RIJIN6
                             CRA12  = CRA1 - CA6
                             VLJA3D = CRA12 * RIJIN6
C     
                             DF3D  =  RIJIN2*(VELA3D+ 
     .                          6.0D0*(CRA1+CRA12)*RIJIN6)+RFF*QIQJA
                             DF3D = DF3D * PIFAC
C     
#ifndef EWALD
                             ERF(ILJ) = ERF(ILJ) +
     $                          QIQJA * RFE * RIJ2 * PIFAC
                             ERC(ILJ) = ERC(ILJ) +
     $                          QIQJA * RFC * PIFAC
#endif
                             ELJ(ILJ) = ELJ(ILJ) + VLJA3D * PIFAC
                             EEL(ILJ) = EEL(ILJ) + VELA3D * PIFAC
                          ENDIF
                           
                          XH = DX*DF3D
                          YH = DY*DF3D
                          ZH = DZ*DF3D
                          FLRX(IX) = FLRX(IX) + XH
                          FLRY(IX) = FLRY(IX) + YH
                          FLRZ(IX) = FLRZ(IX) + ZH
                          FLRX(JX) = FLRX(JX) - XH
                          FLRY(JX) = FLRY(JX) - YH
                          FLRZ(JX) = FLRZ(JX) - ZH
                          IBASE = IBASE + NRP
                        END DO
                  ELSE
                     DO K = 1, NPM
                           IX = I + IBASE
                           JX = J + IBASE
                           DX = (XC(IX)-XC(JX))*BOXIV1
                           DX = (DX-ANINT(DX))*BOX1
                           DY = (YC(IX)-YC(JX))*BOXIV2
                           DY = (DY-ANINT(DY))*BOX2
                           DZ = (ZC(IX)-ZC(JX))*BOXIV3
                           DZ = (DZ-ANINT(DZ))*BOX3
                           RIJ2 = DX*DX + DY*DY + DZ*(DZ+COSB2*DX)
                           ILJ = NRELKP(ICET(IX),ICET(JX))
                           QIQJA = CGAT(IX) * CGAT(JX)
                           IF (LLPERT(IX) .OR. LLPERT(JX)) THEN
                           QIQJB = CGBT(IX) * CGBT(JX)
                           SOFTLJ = ALPHLJ * 
     $                          (1 - ISCLJT(IX) * ISCLJT(JX))
                           SOFTEL = ALPHC * 
     $                          (1 - ISCCT(IX) * ISCCT(JX))
                           
                           RIJ3D2 = RIJ2
                           RIJ4 = RIJ2 * RIJ2
                           RIJ6 = RIJ4 * RIJ2
C     
C--   3D, STATE A
C     
C     COULOMB
                           SFEL   = SOFTEL * RLMA
                           RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELA3D = QIQJA * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFA3 = QIQJA * C1 * RRF
                           QRFR   = DFRFA3 * RIJ3D2
                           VRFA3D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CA126
                           CRA1   = 1.0D0 / (SFLJ * RLMA + RIJ6)
                           CRA12  = CA12 * CRA1
                           VLJA3D = CRA1 * (CRA12 - CA6)
                           DFA3D  = CRA1 * ( CRA1 * CRA12 + VLJA3D)
C     DV/DL
                           DLLJA3 = -2.0D0 * SFLJ * DFA3D
                           DLELA3 = -VELA3D * SOFTEL * RSOFT2 
                           DLRFA3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFA3D = RSOFT2*VELA3D+6.0D0*RIJ4*DFA3D
C     
C--   3D, STATE B
C     
C     COULOMB
                           SFEL   = SOFTEL * RLMB
                           RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELB3D = QIQJB * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFB3 = QIQJB * C1 * RRF
                           QRFR   = DFRFB3 * RIJ3D2
                           VRFB3D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CB126 
                           CRA1   = 1.0D0 / (SFLJ * RLMB + RIJ6)
                           CRA12  = CB12 * CRA1
                           VLJB3D = CRA1 * (CRA12 - CB6)
                           DFB3D  = CRA1 * ( CRA1 * CRA12 + VLJB3D)
C     DV/DL
                           DLLJB3 = -2.0D0 * SFLJ * DFB3D
                           DLELB3 = -VELB3D * SOFTEL * RSOFT2 
                           DLRFB3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFB3D = RSOFT2*VELB3D+6.0D0*RIJ4*DFB3D

C     
C--   COMBINE 3D, STATE A AND B
C     
                           VEL3D = RLB * VELB3D + RLA * VELA3D
                           VRF3D = RLB * VRFB3D + RLA * VRFA3D
                           VLJ3D = RLB * VLJB3D + RLA * VLJA3D
                           DF3D  = RLB * DFB3D  + RLA * DFA3D
     $                           + RLB * DFRFB3 + RLA * DFRFA3
                           DLLJ3D = RDLB * VLJB3D + RLLB * DLLJB3
     $                            + RDLA * VLJA3D + RLLA * DLLJA3
                           DLEL3D = RDLB * VELB3D + RLLB * DLELB3 
     $                            + RDLA * VELA3D + RLLA * DLELA3
                           DLRF3D = RDLB * VRFB3D + RLLB * DLRFB3
     $                            + RDLA * VRFA3D + RLLA * DLRFA3

                           VRCA = RFC*QIQJA
                           VRCB = RFC*QIQJB
                           VRC  = RLB*VRCB  + RLA*VRCA
                           DLRC = RDLB*VRCB + RDLA*VRCA
C
                           EEL(ILJ) = EEL(ILJ) + VEL3D 
#ifndef EWALD
                           ERF(ILJ) = ERF(ILJ) + VRF3D 
                           ERC(ILJ) = ERC(ILJ) + VRC 
#endif
                           ELJ(ILJ) = ELJ(ILJ) + VLJ3D
                           
                           EGLEL = EGLEL + DLEL3D 
                           EGLRF = EGLRF + DLRF3D 
                           EGLRC = EGLRC + DLRC
                           EGLLJ = EGLLJ + DLLJ3D 
C       
                        ELSE
                           RIJIN2 = 1.0D0 / RIJ2
                           RIJINV = DSQRT(RIJIN2)
                           RIJIN6 = RIJIN2 * RIJIN2 * RIJIN2
                           CRA1   = CA12 * RIJIN6
                           CRA12  = CRA1 - CA6
                           VLJA3D = CRA12 * RIJIN6
C
#ifndef EWALD
                           VELA3D = QIQJA * RIJINV
                           DF3D   = RIJIN2*(VELA3D + 
     .                          6.0D0*(CRA1+CRA12)*RIJIN6)+RFF*QIQJA
#else
                           IF (RIJ2 .LE. RCUTP2) THEN
                             VELA3D = QIQJA
     $                         *( RIJINV
     $                          + FACG3*RIJ2/RIJINV
     $                          + FACG2*RIJ2
     $                          + FACG0 )
                             DF3D  = QIQJA*RIJINV
     $                         *(  RIJIN2
     $                          + FACGF1/RIJINV
     $                          + FACGF2*RIJ2 )
     $                          + RIJIN2
     $                          * (6.0D0*(CRA1+CRA12)
     $                          * RIJIN6)
                           ELSE
                             VELA3D = 0.0D0
                             DF3D = 0.0D0
                           ENDIF
#endif
C     
                           DF3D = DF3D  + RFF * QIQJA
                           DF3D = DF3D * PIFAC
C     
#ifndef EWALD
                           ERF(ILJ) = ERF(ILJ) +
     $                          QIQJA * RFE * RIJ2 * PIFAC
                           ERC(ILJ) = ERC(ILJ) +
     $                          QIQJA * RFC * PIFAC
#endif
                           ELJ(ILJ) = ELJ(ILJ) + VLJA3D * PIFAC

                           EEL(ILJ) = EEL(ILJ) + VELA3D * PIFAC
                        ENDIF
                           
                        XH = DX*DF3D 
                        YH = DY*DF3D 
                        ZH = DZ*DF3D 
                        FLRX(IX) = FLRX(IX) + XH
                        FLRY(IX) = FLRY(IX) + YH
                        FLRZ(IX) = FLRZ(IX) + ZH
                        FLRX(JX) = FLRX(JX) - XH
                        FLRY(JX) = FLRY(JX) - YH
                        FLRZ(JX) = FLRZ(JX) - ZH
                        IBASE = IBASE + NRP
                     END DO
                  ENDIF
               ENDIF
 2600       CONTINUE
 2650    CONTINUE
 2700 CONTINUE
C     
C     
C*****TREAT THE FORCES OF ALL PAIRS DENOTED BY THE PAIRLIST.
C*****FIRST, ALL CHARGE GROUPS NI BELONGING TO THE SOLUTES WILL BE PROC
C     
C     
      JN = 1
      NJNB = 1
      DO 4300 NI = 1, NRAG
         NPAIR = INB(NI)
C     
         IF (NPAIR .GT. 0) THEN
!$OMP parallel do private(IPQ) shared(NJNB)
            DO IPQ = 1, NPAIR
               JFIRST(IPQ) = IFIRST(JNB(IPQ+NJNB-1))
               NJL(IPQ) = INAG(JNB(IPQ+NJNB-1))
            ENDDO
C*****CODE EXPANDED FROM ROUTINE:  VEXPLD
C     
               JRCNUM = 1
               JRCL1X = NPAIR
               J1X = 0
77001          CONTINUE
               K1X = 0

               DO I1X = 1, JRCL1X
                  IF (NJL(I1X) .GE. JRCNUM) THEN
                     J1X = J1X + 1
                     JLIST(J1X) = JFIRST(I1X)
                     K1X = K1X + 1
                     JFIRST(K1X) = JFIRST(I1X) + 1
                     NJL(K1X) = NJL(I1X)
                  ENDIF
               END DO
               JRCL1X = K1X
               JRCNUM = JRCNUM + 1
               IF (JRCL1X .GT. 0) GO TO 77001
               JLEN = J1X
C*****END OF CODE EXPANDED FROM ROUTINE:  VEXPLD
               NJNB = NJNB + NPAIR
C     
               ILAST = IFIRST(NI) + INAG(NI) - 1
C     
C*****GATHER THE COORDINATES AND THE FORCES.
C     
               IF (LPERTL) THEN
!$OMP parallel do private(IPQ)
                  DO IPQ = 1, JLEN
                     XJ(IPQ) = XC(JLIST(IPQ))
                     YJ(IPQ) = YC(JLIST(IPQ))
                     ZJ(IPQ) = ZC(JLIST(IPQ))
C     
                     FXTEMP(IPQ) = FLRX(JLIST(IPQ))
                     FYTEMP(IPQ) = FLRY(JLIST(IPQ))
                     FZTEMP(IPQ) = FLRZ(JLIST(IPQ))
C     
C*****DETERMINE THE INTEGER ATOM CODE FOR EACH ATOM J.
C*****GATHER THE CHARGES AND DETERMINE TO WHICH GROUP EACH ATOM BELO
C     
                     IACAJ(IPQ) = IACT(JLIST(IPQ)) + JRCBK1(IPQ)
                     CGAJ(IPQ) = CGAT(JLIST(IPQ))
                     JRCILJ(IPQ) = ICET(JLIST(IPQ))
                     AMULT(IPQ) = 1.0D0
                     IFORCE(IPQ) = 0
                     IACBJ(IPQ) = IACBT(JLIST(IPQ))
                     CGBJ(IPQ) = CGBT(JLIST(IPQ))
                     LLPERJ(IPQ) = LLPERT(JLIST(IPQ))
                     ISCLJJ(IPQ) = ISCLJT(JLIST(IPQ))
                     ISCCJ(IPQ) = ISCCT(JLIST(IPQ))
                     IF (LDOVIR) NSPJ(IPQ) = NSPTT(JLIST(IPQ))
                  ENDDO
               ELSE
                  IF (LPION) THEN
!$OMP parallel do private(IPQ)
                     DO IPQ = 1, JLEN
                        IPICJ(IPQ) = IPICT(JLIST(IPQ))
                     ENDDO
                  ENDIF
                  IF (LDOVIR) THEN
!$OMP parallel do private(IPQ)
                     DO IPQ = 1, JLEN
                        XJ(IPQ) = XC(JLIST(IPQ))
                        YJ(IPQ) = YC(JLIST(IPQ))
                        ZJ(IPQ) = ZC(JLIST(IPQ))
C     
                        FXTEMP(IPQ) = FLRX(JLIST(IPQ))
                        FYTEMP(IPQ) = FLRY(JLIST(IPQ))
                        FZTEMP(IPQ) = FLRZ(JLIST(IPQ))
C     
C*****DETERMINE THE INTEGER ATOM CODE FOR EACH ATOM J.
C*****GATHER THE CHARGES AND DETERMINE TO WHICH GROUP EACH ATOM BELO
C     
                        IACAJ(IPQ) = IACT(JLIST(IPQ)) + JRCBK1(IPQ)
                        CGAJ(IPQ) = CGAT(JLIST(IPQ))
                        JRCILJ(IPQ) = ICET(JLIST(IPQ))
                        AMULT(IPQ) = 1.0D0
                        IFORCE(IPQ) = 0
                        NSPJ(IPQ) = NSPTT(JLIST(IPQ))
                     ENDDO
                  ELSE
!$OMP parallel do private(IPQ)
#ifdef _vpp_
!OCL NOALIAS 
#endif
                     DO IPQ = 1, JLEN
                        XJ(IPQ) = XC(JLIST(IPQ))
                        YJ(IPQ) = YC(JLIST(IPQ))
                        ZJ(IPQ) = ZC(JLIST(IPQ))
C     
                        FXTEMP(IPQ) = FLRX(JLIST(IPQ))
                        FYTEMP(IPQ) = FLRY(JLIST(IPQ))
                        FZTEMP(IPQ) = FLRZ(JLIST(IPQ))
C     
C*****DETERMINE THE INTEGER ATOM CODE FOR EACH ATOM J.
C*****GATHER THE CHARGES AND DETERMINE TO WHICH GROUP EACH ATOM BELO
C     
                        IACAJ(IPQ) = IACT(JLIST(IPQ)) + JRCBK1(IPQ)
                        CGAJ(IPQ) = CGAT(JLIST(IPQ))
                        JRCILJ(IPQ) = ICET(JLIST(IPQ))
                        AMULT(IPQ) = 1.0D0
                        SCALE(IPQ)=1.0D0
                        IFORCE(IPQ) = 0
                     END DO
                  ENDIF
               ENDIF
C     
C*****FOR EACH ATOM I MUST BE DETERMINED WHETHER IT MUST BE EXCLUDED
C*****SOMEHOW FROM ATOM I. SOME DATA FOR THIS DETERMINATION CAN BE
C*****GATHERED FOR ALL ATOMS J.
C     
               IBASE = INT((NI-1)/NCAG)*NRP
               JLEN2 = 0
               DO IPQ = 1, JLEN
                  IF (JLIST(IPQ).LE.IBASE+NRP .AND. JLIST(IPQ).GT.
     .                 IBASE) THEN
                     JLEN2 = JLEN2 + 1
                     JRCIND(JLEN2) = IPQ
                     JFIRST(JLEN2) = JLIST(IPQ)
                  ENDIF
               ENDDO
               
               IF (JLEN2 .GT. 0) II = IFIRST(NI) - IBASE
C     
               IF (LDOVIR) THEN
                  IVIRL = 0
                  DO IPQ = 1, JLEN
                     IF (NSPJ(IPQ) .NE. NSPTT(ILAST)) THEN
                        IVIRL = IVIRL + 1
                        IOCVIR(IVIRL) = IPQ
                     ENDIF
                  ENDDO
                  LVIR = .FALSE.
                  IF (IVIRL .GT. 0) THEN
                     LVIR = .TRUE.
!$OMP parallel do private(I9X)
                     DO I9X = 1, JLEN
                        XRJ(I9X) = XRC(JLIST(I9X))
                        YRJ(I9X) = YRC(JLIST(I9X))
                        ZRJ(I9X) = ZRC(JLIST(I9X))
                     ENDDO
                  ENDIF
               ENDIF
C     
C     USE 'IFORCE' TO REPLACE THE EXCLUDED PAIR-PARAMETERS BY 1-4
C     PARAMETERS.
C     
C     USE 'AMULT' TO CONTROL THE ASSIGNMENTS TO THE INVERSE DISTANCE.
C     
               DO 4200 I = IFIRST(NI), ILAST

                  ICE = ICET(I)

C--   CHECK FOR EXCLUDED ATOMS IN THIS PASS
                  IF (JLEN2 .GT. 0) THEN
C--   SET UP THE DEFAULT ARRAYS FOR EXCLUSIONS AND 1-4 INTERACTIONS
!$OMP parallel do private(K) shared(IBASE)
                     DO K = 1, JLEN2
                        JFORCE(JFIRST(K)-IBASE) = 0
                        BMULT(JFIRST(K)-IBASE) = 1.0D0
                        BSCALE(JFIRST(K)-IBASE) = 1.0D0
                     ENDDO
C     
                     ILEN = INE14(II)
                     DO K = 1, ILEN
                        JFORCE(JSNE14(K+KNE14(II))) = KFORCE
                        BSCALE(JSNE14(K+KNE14(II))) = scale_14_amber
                     ENDDO
                     ILEN = INE(II)
!$OMP parallel do private(K) shared(II)
                     DO K = 1, ILEN
                        BMULT(JSNE(K+KNE(II))) = 0.0D0
                     ENDDO
C     
!$OMP parallel do private(K) shared(IBASE)
                     DO K = 1, JLEN2
                        IFORCE(JRCIND(K)) = JFORCE(JFIRST(K)-IBASE)
                        AMULT(JRCIND(K)) = BMULT(JFIRST(K)-IBASE)
                        SCALE(JRCIND(K)) = BSCALE(JFIRST(K)-IBASE)
                     ENDDO
                     
                     II = II + 1
                  ENDIF
C     
C--   DETERMINE LJ-PARAMETERS.
C     
                  IF (LPERTL) THEN
!$OMP parallel do private(IPQ,NJIND,NJINDB)
                     DO IPQ = 1, JLEN
                        NJIND = JRCPAC(IACAJ(IPQ)+(IACT(I)-1)*NRATT) + 
     $                       IFORCE(IPQ) + JRCBK2(IPQ)
                        NJINDB = JRCPAC(IACBJ(IPQ)+(IACBT(I)-1)*NRATT) 
     $                       + IFORCE(IPQ) + JRCBK2(IPQ)
                        VC6A(IPQ) = FORC6(NJIND)
                        VC12A(IPQ) = FORC12(NJIND)
                        VC6B(IPQ) = FORC6(NJINDB)
                        VC12B(IPQ) = FORC12(NJINDB)
                        AC12C6(IPQ) = FC12C6(NJIND)
                        BC12C6(IPQ) = FC12C6(NJINDB)
                     ENDDO
                  ELSE
!$OMP parallel do private(IPQ,NJIND)
#ifdef _vpp_
!OCL NOALIAS 
#endif
                     DO IPQ = 1, JLEN
                        NJIND = JRCPAC(IACAJ(IPQ)+(IACT(I)-1)*NRATT) + 
     $                       IFORCE(IPQ) + JRCBK2(IPQ)
                        VC6A(IPQ) = FORC6(NJIND)
                        VC12A(IPQ) = FORC12(NJIND)
                     ENDDO
                  ENDIF
C     

                  IF (LVAC) THEN
!$OMP parallel do private(IPQ) shared(I)
                     DO IPQ = 1, JLEN
                        VXIJ(IPQ) = XC(I) - XJ(IPQ)
                        VYIJ(IPQ) = YC(I) - YJ(IPQ)
                        VZIJ(IPQ) = ZC(I) - ZJ(IPQ)
                        VRIJ2(IPQ) = VXIJ(IPQ)*VXIJ(IPQ) + 
     .                               VYIJ(IPQ)*VYIJ(IPQ) + 
     .                   VZIJ(IPQ)*(VZIJ(IPQ)+COSB2*VXIJ(IPQ))
                     ENDDO
                  ELSE
!$OMP parallel do private(IPQ,XIJ,YIJ,ZIJ)
#ifdef _vpp_
!OCL NOALIAS 
#endif
                     DO IPQ = 1, JLEN
                        XIJ = (XC(I)-XJ(IPQ))*BOXIV1
                        YIJ = (YC(I)-YJ(IPQ))*BOXIV2
                        ZIJ = (ZC(I)-ZJ(IPQ))*BOXIV3
C     
                        VXIJ(IPQ) = (XIJ-ANINT(XIJ))*BOX1
                        VYIJ(IPQ) = (YIJ-ANINT(YIJ))*BOX2
                        VZIJ(IPQ) = (ZIJ-ANINT(ZIJ))*BOX3
C     
                        VRIJ2(IPQ) = VXIJ(IPQ)*VXIJ(IPQ) + 
     $                               VYIJ(IPQ)*VYIJ(IPQ) +
     $                      VZIJ(IPQ)*(VZIJ(IPQ)+COSB2*VXIJ(IPQ))
                     ENDDO
                  ENDIF
                  IF (LOCTO) THEN
                     DO 3350 IPQ = 1, JLEN
                        VEELJ(IPQ) = BOXOQ - DABS(VXIJ(IPQ)) -
     .                             DABS(VYIJ(IPQ))-DABS(VZIJ(IPQ))
                        IF (VEELJ(IPQ) .LE. 0.D0) THEN
                           VRIJ2(IPQ)=VRIJ2(IPQ)+VEELJ(IPQ)*BOX1
                           VXIJ(IPQ) =VXIJ(IPQ)-SIGN(BOXOH,VXIJ(IPQ))
                           VYIJ(IPQ) =VYIJ(IPQ)-SIGN(BOXOH,VYIJ(IPQ))
                           VZIJ(IPQ) =VZIJ(IPQ)-SIGN(BOXOH,VZIJ(IPQ))
                        ENDIF
 3350                CONTINUE
                  ENDIF
C     
                  IPQ = ISRCHFLE(JLEN,VRIJ2,1,1.0D-20)
                  IF (IPQ .LE. JLEN) THEN
                     PRINT 9240, NI, I
 9240                FORMAT(/' A DISTANCE THAT IS TOO SMALL HAS',
     $                    ' BEEN FOUND',
     $                    ' FOR CHARGE GROUP ',I4,'  AND ATOM ',I5/)
                     LERROR = .TRUE.
                     GO TO 5750
                  ENDIF
C     
!$OMP parallel do private(IPQ)
#ifdef _vpp_
!OCL NOALIAS 
#endif
                  DO IPQ = 1, JLEN
                     VRIJ(IPQ) = DSQRT(VRIJ2(IPQ))
                     VRINV(IPQ) = AMULT(IPQ)/VRIJ(IPQ)
                     VRINV2(IPQ) = VRINV(IPQ)*VRINV(IPQ)
                     VRINV6(IPQ) = VRINV2(IPQ)*VRINV2(IPQ)*
     $                             VRINV2(IPQ)
                  ENDDO
C     
C     
                  IF (LPION) THEN
                     LIZERO = (IPICT(I) .EQ. 0)   
                     DO IPQ = 1, JLEN
                        LJZERO = (IPICJ(IPQ) .EQ. 0)
                        IF (LIZERO .AND. LJZERO) THEN
                           PIFACT(IPQ) = 1.0D0
                        ELSEIF (LIZERO .OR. LJZERO) THEN
                           PIFACT(IPQ) = PININV
                        ELSEIF (IPICT(I) .EQ. IPICJ(IPQ)) THEN
                           PIFACT(IPQ) = PININV
                        ELSE
                           PIFACT(IPQ) = 0.0D0
                        ENDIF
                     ENDDO
                  ENDIF 
C     
C     
                  IF (LPERTL) THEN
                     LIPERT = LLPERT(I)
                     IF (LIPERT) THEN
!$OMP parallel do private(IPQ,CA6,CA12,CB6,CB12,QIQJA,QIQJB,CA126,CB126)
!$OMP+ private(SOFTLJ,SOFTEL,RIJ3D2,RIJ4,RIJ6,SFEL,RSOFT2,RSOFT,VELA3D)
!$OMP+ private(RRF2,RFSQRT,RRF,DFRFA3,QRFR,VRFA3D,SFLJ,CRA1,CRA12)
!$OMP+ private(VLJA3D,DFA3D,DLLJA3,DLELA3,DLRFA3,VEL3D,VRF3D)
!$OMP+ private(VLJ3D,DF3D,DLLJ3D,DLEL3D,DLRF3D,VRCA,VRCB,VRC,DLRC)
                        DO IPQ = 1, JLEN
C     
                           CA6 = VC6A(IPQ)
                           CA12 = VC12A(IPQ)
                           CB6 = VC6B(IPQ)
                           CB12 = VC12B(IPQ)

                           QIQJA = CGAT(I) * CGAJ(IPQ) * AMULT(IPQ)
                           QIQJB = CGBT(I) * CGBJ(IPQ) * AMULT(IPQ)

                           CA126 = AC12C6(IPQ)
                           CB126 = BC12C6(IPQ)

                           SOFTLJ = ALPHLJ * 
     $                          (1 - ISCLJT(I) * ISCLJJ(IPQ))
                           SOFTEL = ALPHC * 
     $                          (1 - ISCCT(I) * ISCCJ(IPQ))
C     
C     CALCULATE FULL 3D CONTRIBUTION
C     
                           RIJ3D2 = VRIJ2(IPQ)
                           RIJ4 = RIJ3D2 * RIJ3D2
                           RIJ6 = RIJ4 * RIJ3D2
C     
C--   3D, STATE A
C     
C     COULOMB
                           SFEL   = SOFTEL * RLMA
                           RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELA3D = QIQJA * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFA3 = QIQJA * C1 * RRF
                           QRFR   = DFRFA3 * RIJ3D2
                           VRFA3D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CA126
                           CRA1   = 1.0D0 / (SFLJ * RLMA + RIJ6)
                           CRA12  = CA12 * CRA1
                           VLJA3D = CRA1 * (CRA12 - CA6)
                           DFA3D  = CRA1 * ( CRA1 * CRA12 + VLJA3D)
C     DV/DL
                           DLLJA3 = -2.0D0 * SFLJ * DFA3D
                           DLELA3 = -VELA3D * SOFTEL * RSOFT2 
                           DLRFA3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFA3D = RSOFT2*VELA3D+6.0D0*RIJ4*DFA3D
C     
C--   3D, STATE B
C     
C     COULOMB
                           SFEL   = SOFTEL * RLMB
                           RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELB3D = QIQJB * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFB3 = QIQJB * C1 * RRF
                           QRFR   = DFRFB3 * RIJ3D2
                           VRFB3D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CB126 
                           CRA1   = 1.0D0 / (SFLJ * RLMB + RIJ6)
                           CRA12  = CB12 * CRA1
                           VLJB3D = CRA1 * (CRA12 - CB6)
                           DFB3D  = CRA1 * ( CRA1 * CRA12 + VLJB3D)
C     DV/DL
                           DLLJB3 = -2.0D0 * SFLJ * DFB3D
                           DLELB3 = -VELB3D * SOFTEL * RSOFT2 
                           DLRFB3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFB3D = RSOFT2*VELB3D+6.0D0*RIJ4*DFB3D
C     
C--   COMBINE 3D, STATE A AND B
C     
                           VEL3D =  RLB  * VELB3D + RLA  * VELA3D
                           VRF3D =  RLB  * VRFB3D + RLA  * VRFA3D
                           VLJ3D =  RLB  * VLJB3D + RLA  * VLJA3D
                           DF3D  =  RLB  * DFB3D  + RLA  * DFA3D
     $                            + RLB  * DFRFB3 + RLA  * DFRFA3
                           DLLJ3D = RDLB * VLJB3D + RLLB * DLLJB3
     $                            + RDLA * VLJA3D + RLLA * DLLJA3
                           DLEL3D = RDLB * VELB3D + RLLB * DLELB3 
     $                            + RDLA * VELA3D + RLLA * DLELA3
                           DLRF3D = RDLB * VRFB3D + RLLB * DLRFB3
     $                            + RDLA * VRFA3D + RLLA * DLRFA3

                           VRCA = RFC*QIQJA
                           VRCB = RFC*QIQJB
                           VRC  = RLB*VRCB  + RLA*VRCA
                           DLRC = RDLB*VRCB + RDLA*VRCA
C     
                           VEELJ(IPQ) = VEL3D 
                           VERFJ(IPQ) = VRF3D 
                           VERCJ(IPQ) = VRC * AMULT(IPQ)
                           VLJJ(IPQ)  = VLJ3D * AMULT(IPQ)
                           
                           EGLELJ(IPQ) = DLEL3D * AMULT(IPQ)
                           EGLRFJ(IPQ) = DLRF3D * AMULT(IPQ)
                           EGLRCJ(IPQ) = DLRC * AMULT(IPQ)
                           EGLLJJ(IPQ) = DLLJ3D * AMULT(IPQ)
C     
                           DF3(IPQ) = DF3D * AMULT(IPQ)
                           
                        ENDDO

                     ELSE
                        DO IPQ = 1, JLEN
                           IF (LLPERJ(IPQ)) THEN

                           CA6 = VC6A(IPQ)
                           CA12 = VC12A(IPQ)
                           CB6 = VC6B(IPQ)
                           CB12 = VC12B(IPQ)

                           QIQJA = CGAT(I) * CGAJ(IPQ) * AMULT(IPQ)
                           QIQJB = CGBT(I) * CGBJ(IPQ) * AMULT(IPQ)

                           CA126 = AC12C6(IPQ)
                           CB126 = BC12C6(IPQ)

                           SOFTLJ = ALPHLJ * 
     $                          (1 - ISCLJT(I) * ISCLJJ(IPQ))
                           SOFTEL = ALPHC * 
     $                          (1 - ISCCT(I) * ISCCJ(IPQ))
C     
C     CALCULATE FULL 3D CONTRIBUTION
C     
                           RIJ3D2 = VRIJ2(IPQ)
                           RIJ4 = RIJ3D2 * RIJ3D2
                           RIJ6 = RIJ4 * RIJ3D2
C     
C--   3D, STATE A
C     
C     COULOMB
                           SFEL   = SOFTEL * RLMA
                           RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELA3D = QIQJA * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFA3 = QIQJA * C1 * RRF
                           QRFR   = DFRFA3 * RIJ3D2
                           VRFA3D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CA126
                           CRA1   = 1.0D0 / (SFLJ * RLMA + RIJ6)
                           CRA12  = CA12 * CRA1
                           VLJA3D = CRA1 * (CRA12 - CA6)
                           DFA3D  = CRA1 * ( CRA1 * CRA12 + VLJA3D)
C     DV/DL
                           DLLJA3 = -2.0D0 * SFLJ * DFA3D
                           DLELA3 = -VELA3D * SOFTEL * RSOFT2 
                           DLRFA3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFA3D = RSOFT2*VELA3D+6.0D0*RIJ4*DFA3D
C     
C--   3D, STATE B
C     
C     COULOMB
                           SFEL   = SOFTEL * RLMB
                           RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELB3D = QIQJB * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFB3 = QIQJB * C1 * RRF
                           QRFR   = DFRFB3 * RIJ3D2
                           VRFB3D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CB126 
                           CRA1   = 1.0D0 / (SFLJ * RLMB + RIJ6)
                           CRA12  = CB12 * CRA1
                           VLJB3D = CRA1 * (CRA12 - CB6)
                           DFB3D  = CRA1 * ( CRA1 * CRA12 + VLJB3D)
C     DV/DL
                           DLLJB3 = -2.0D0 * SFLJ * DFB3D
                           DLELB3 = -VELB3D * SOFTEL * RSOFT2 
                           DLRFB3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFB3D = RSOFT2*VELB3D+6.0D0*RIJ4*DFB3D
C     
C--   COMBINE 3D, STATE A AND B
C     
                           VEL3D =  RLB  * VELB3D + RLA  * VELA3D
                           VRF3D =  RLB  * VRFB3D + RLA  * VRFA3D
                           VLJ3D =  RLB  * VLJB3D + RLA  * VLJA3D
                           DF3D  =  RLB  * DFB3D  + RLA  * DFA3D
     $                            + RLB  * DFRFB3 + RLA  * DFRFA3
                           DLLJ3D = RDLB * VLJB3D + RLLB * DLLJB3
     $                            + RDLA * VLJA3D + RLLA * DLLJA3
                           DLEL3D = RDLB * VELB3D + RLLB * DLELB3 
     $                            + RDLA * VELA3D + RLLA * DLELA3
                           DLRF3D = RDLB * VRFB3D + RLLB * DLRFB3
     $                            + RDLA * VRFA3D + RLLA * DLRFA3

                           VRCA = RFC*QIQJA
                           VRCB = RFC*QIQJB
                           VRC  = RLB*VRCB  + RLA*VRCA
                           DLRC = RDLB*VRCB + RDLA*VRCA
C     
                           VEELJ(IPQ) = VEL3D 
                           VERFJ(IPQ) = VRF3D 
                           VERCJ(IPQ) = VRC * AMULT(IPQ)
                           VLJJ(IPQ)  = VLJ3D * AMULT(IPQ)
                           
                           EGLELJ(IPQ) = DLEL3D * AMULT(IPQ)
                           EGLRFJ(IPQ) = DLRF3D * AMULT(IPQ)
                           EGLRCJ(IPQ) = DLRC * AMULT(IPQ)
                           EGLLJJ(IPQ) = DLLJ3D * AMULT(IPQ)
C     
                           DF3(IPQ) = DF3D * AMULT(IPQ)
                           
                           ELSE
                              QIQJA = CGAT(I) * CGAJ(IPQ) * AMULT(IPQ)
                              VERFJ(IPQ) = QIQJA * RFE * VRIJ2(IPQ)
                              VERCJ(IPQ) = QIQJA * RFC
                              VEELJ(IPQ) = QIQJA * VRINV(IPQ)
                              VLJJ(IPQ) = VRINV6(IPQ)*(VC12A(IPQ)
     $                             *VRINV6(IPQ)-VC6A(IPQ))

                              EGLELJ(IPQ) = 0.0D0
                              EGLRFJ(IPQ) = 0.0D0
                              EGLRCJ(IPQ) = 0.0D0
                              EGLLJJ(IPQ) = 0.0D0

                              DF = VRINV6(IPQ)*(2.0D0*VRINV6(IPQ) *
     $                             VC12A(IPQ)-VC6A(IPQ))
                              DF3(IPQ) = (6.0D0 * DF + 
     $                             VEELJ(IPQ))*VRINV2(IPQ) 
     $                             +  QIQJA * RFF
                              DF3(IPQ) = DF3(IPQ)
                           ENDIF

                        ENDDO

                     ENDIF
!$OMP parallel do private(IPQ) reduction(+:EGLEL,EGLRF,EGLRC,EGLLJ)
                     DO  IPQ = 1, JLEN
                        EGLEL = EGLEL + EGLELJ(IPQ)
                        EGLRF = EGLRF + EGLRFJ(IPQ)
                        EGLRC = EGLRC + EGLRCJ(IPQ)
                        EGLLJ = EGLLJ + EGLLJJ(IPQ)
                     ENDDO

                  ELSE
                     
                     IF (LPION) THEN
!$OMP parallel do private(IPQ,QIQJA,DF)
                        DO IPQ = 1, JLEN
                           QIQJA = CGAT(I) * CGAJ(IPQ) * AMULT(IPQ) 
     $                          * PIFACT(IPQ)
                           VERFJ(IPQ) = QIQJA * RFE * VRIJ2(IPQ)
                           VERCJ(IPQ) = QIQJA * RFC
                           VEELJ(IPQ) = QIQJA * VRINV(IPQ)
                           VLJJ(IPQ) = VRINV6(IPQ)*(VC12A(IPQ)*VRINV6(
     $                          IPQ)-VC6A(IPQ)) * PIFACT(IPQ)
                           DF = VRINV6(IPQ)*(2.0D0*VRINV6(IPQ)*
     $                          VC12A(IPQ)-VC6A(IPQ)) * PIFACT(IPQ)
                           DF3(IPQ) = (6.0D0 * DF + VEELJ(IPQ))*
     $                          VRINV2(IPQ) +  QIQJA * RFF
                        ENDDO
                     ELSE
!$OMP parallel do private(IPQ,QIQJA,DF)
#ifdef _vpp_
!OCL NOALIAS 
#endif
                        DO IPQ = 1, JLEN
                           QIQJA = CGAT(I)*CGAJ(IPQ)*
     $                     AMULT(IPQ)*SCALE(IPQ)
                           VERFJ(IPQ) = QIQJA * RFE * VRIJ2(IPQ)
                           VERCJ(IPQ) = QIQJA * RFC
#ifndef EWALD
                           VEELJ(IPQ) = QIQJA*VRINV(IPQ)
                           DVEELJ(IPQ)= VEELJ(IPQ) * VRINV2(IPQ)+
     $                             6.D0*VRINV6(IPQ)*(2.0D0*VRINV6(IPQ)*
     $                             VC12A(IPQ)-VC6A(IPQ)) * VRINV2(IPQ)
#else
                           IF (VRIJ2(IPQ) .LE. RCUTP2) THEN 
                             VEELJ(IPQ) = QIQJA* (  VRINV(IPQ)
     $                          + FACG3* VRIJ(IPQ)*VRIJ2(IPQ)
     $                          + FACG2*VRIJ2(IPQ)
     $                          + FACG0 )
                             DVEELJ(IPQ)=  QIQJA* VRINV2(IPQ)*
     $                            (VRINV(IPQ)
     $                          +  FACGF2 * VRIJ(IPQ)*VRIJ2(IPQ)
     $                          +  FACGF1 * VRIJ2(IPQ)) +
     $                             6.D0*VRINV6(IPQ)*(2.0D0*VRINV6(IPQ)*
     $                             VC12A(IPQ)-VC6A(IPQ)) * VRINV2(IPQ) 
                           ELSE
                             VEELJ(IPQ) = 0.D0
                             DVEELJ(IPQ)= 0.D0
                             DF=0.D0
                           ENDIF
#endif 
                           VLJJ(IPQ) = VRINV6(IPQ)*(VC12A(IPQ)*VRINV6(
     $                          IPQ)-VC6A(IPQ))
                           DF3(IPQ) =  DVEELJ(IPQ) +  QIQJA * RFF  
                        ENDDO
                     ENDIF
                  ENDIF
C     
C     
                  JRCTOT = JLEN
                  JRCSRT = 0
 3950             CONTINUE
                  JRCLEN = MIN(IVECT,JRCTOT)
C     DIR$ SHORTLOOP
#ifdef _vpp_
!OCL NOALIAS 
#endif
                  DO JRC = 1, JRCLEN
                     TEEL(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) = 
     $                    TEEL(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) + 
     $                    VEELJ(JRC+JRCSRT)
                     TELJ(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) = 
     $                    TELJ(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) + 
     $                    VLJJ(JRC+JRCSRT)
                     
                     TERF(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TERF(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VERFJ(JRC+JRCSRT)
                     
                     TERC(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TERC(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VERCJ(JRC+JRCSRT)
                  ENDDO

                  JRCSRT = JRCSRT + JRCLEN
                  JRCTOT = JRCTOT - JRCLEN
                  IF (JRCTOT .GT. 0) GO TO 3950
C     
C*****CALCULATE THE FORCES.
C     
!$OMP parallel do private(IPQ)
#ifdef _vpp_
!OCL NOALIAS 
#endif
                  DO IPQ = 1, JLEN
                     VXH(IPQ) = VXIJ(IPQ)*DF3(IPQ)
                     VYH(IPQ) = VYIJ(IPQ)*DF3(IPQ)
                     VZH(IPQ) = VZIJ(IPQ)*DF3(IPQ)
C     
                     FXTEMP(IPQ) = FXTEMP(IPQ) - VXH(IPQ)
                     FYTEMP(IPQ) = FYTEMP(IPQ) - VYH(IPQ)
                     FZTEMP(IPQ) = FZTEMP(IPQ) - VZH(IPQ)
                  ENDDO
C
                  IF (LDOVIR) THEN
                     IF (LVIR) THEN
!$OMP parallel do private(IPQ) shared(I)
                        DO IPQ = 1, IVIRL
                           TVIR(IPQ,1) = TVIR(IPQ,1) - VXH(IOCVIR(
     $                          IPQ))*(VXIJ(IOCVIR(IPQ))+XRJ(IOCVIR(
     $                          IPQ))-XRC(I))
                           TVIR(IPQ,2) = TVIR(IPQ,2) - VYH(IOCVIR(
     $                          IPQ))*(VYIJ(IOCVIR(IPQ))+YRJ(IOCVIR(
     $                          IPQ))-YRC(I))
                           TVIR(IPQ,3) = TVIR(IPQ,3) - VZH(IOCVIR(
     $                          IPQ))*(VZIJ(IOCVIR(IPQ))+ZRJ(IOCVIR(
     $                          IPQ))-ZRC(I))
                        ENDDO
                     ENDIF
                  ENDIF
C     
!$OMP parallel do private(IPQ) shared(I)
#ifdef _vpp_
!OCL NOALIAS 
#endif
                  DO IPQ = 1, JLEN
                     FLRX(I) = FLRX(I) + VXH(IPQ)
                     FLRY(I) = FLRY(I) + VYH(IPQ)
                     FLRZ(I) = FLRZ(I) + VZH(IPQ)
                  ENDDO
C     
                  JN = JN + JLEN
 4200          CONTINUE
C     
C     
!$OMP parallel do private(I20X)
               DO I20X = 1, JLEN
                  FLRX(JLIST(I20X)) = FXTEMP(I20X)
                  FLRY(JLIST(I20X)) = FYTEMP(I20X)
                  FLRZ(JLIST(I20X)) = FZTEMP(I20X)
               ENDDO
C     
            ENDIF
 4300    CONTINUE
C     
C*****CALCULATE THE INTERACTION BETWEEN THE SOLVENT MOLECULES.
C     
         NIMAX = NRAGT - 1
         DO 5500 NI = NRAG + 1, NIMAX
            NPAIR = INB(NI)
            IF (NPAIR .GT. 0) THEN
!$OMP parallel do private(I44X)
#ifdef _vpp_
!OCL NOALIAS 
#endif
               DO I44X = 1, NPAIR
                  JLIST(I44X) = IFIRST(JNB(I44X+NJNB-1))
               ENDDO
               JLOW = 1
               JLEN = NPAIR + 1
               DO K = 2, NRAM
!$OMP parallel do private(IPQ)
#ifdef _vpp_
!OCL NOVREC
#endif
                  DO IPQ = 1, NPAIR
                     JLIST(JLEN+IPQ-1) = JLIST(JLOW+IPQ-1) + 1
                  ENDDO
                  JLOW = JLEN
                  JLEN = JLEN + NPAIR
               ENDDO
C     
               JLEN = JLEN - 1
               NJNB = NJNB + NPAIR
               
               ILAST = IFIRST(NI) + NRAM - 1

C     
C*****GATHER THE COORDINATES AND THE FORCES.
C     
               IF (LDOVIR) THEN
!$OMP parallel do private(IPQ)
#ifdef _vpp_
!OCL NOALIAS
#endif
                  DO IPQ = 1, JLEN
                     XJ(IPQ) = XC(JLIST(IPQ))
                     YJ(IPQ) = YC(JLIST(IPQ))
                     ZJ(IPQ) = ZC(JLIST(IPQ))
C     
                     FXTEMP(IPQ) = FLRX(JLIST(IPQ))
                     FYTEMP(IPQ) = FLRY(JLIST(IPQ))
                     FZTEMP(IPQ) = FLRZ(JLIST(IPQ))
                     XRJ(IPQ) = XRC(JLIST(IPQ))
                     YRJ(IPQ) = YRC(JLIST(IPQ))
                     ZRJ(IPQ) = ZRC(JLIST(IPQ))
                     JRCILJ(IPQ) = ICET(JLIST(IPQ))
                  ENDDO
               ELSE
!$OMP parallel do private(IPQ)
#ifdef _vpp_
!OCL NOALIAS 
#endif
                  DO IPQ = 1, JLEN
                     XJ(IPQ) = XC(JLIST(IPQ))
                     YJ(IPQ) = YC(JLIST(IPQ))
                     ZJ(IPQ) = ZC(JLIST(IPQ))
C     
                     FXTEMP(IPQ) = FLRX(JLIST(IPQ))
                     FYTEMP(IPQ) = FLRY(JLIST(IPQ))
                     FZTEMP(IPQ) = FLRZ(JLIST(IPQ))
                     JRCILJ(IPQ) = ICET(JLIST(IPQ))
                  ENDDO
               ENDIF
C     
C*****DETERMINE THE INTEGER ATOM CODE FOR EACH ATOM J.
C*****GATHER THE CHARGES AND DETERMINE TO WHICH GROUP EACH ATOM BELO
C     
cmb             JTOP = 1
!$OMP parallel do private(K,IPQ,JTOP)
               DO K = 1, NRAM
                 JTOP = 1 + (K-1)*NPAIR 
#ifdef _vpp_
!OCL NOALIAS 
#endif
                  DO IPQ = 1, NPAIR
                     IACAJ(JTOP+IPQ-1) = IACS(K)
                     CGAJ(JTOP+IPQ-1) = CGS(K)
                  ENDDO
cmb                JTOP = JTOP + NPAIR
               ENDDO
C     
               DO 5400 I = IFIRST(NI), ILAST
                  ICE = ICET(I)
C     
C*****DETERMINE LJ-PARAMETERS.
C     
cmb                JTOP = 1
!$OMP parallel do private(K,IPQ,IPC,JTOP) shared(I)
                  DO K = 1, NRAM
                     IPC = MPAC(IACS(K),IACT(I))
                     JTOP = 1 + (K-1)*NPAIR
#ifdef _vpp_
!OCL NOALIAS 
#endif
                     DO IPQ = 1, NPAIR
                        VC6A(JTOP+IPQ-1) = C6(IPC)
                        VC12A(JTOP+IPQ-1) = C12(IPC)
                     ENDDO
cmb                  JTOP = JTOP + NPAIR
                  ENDDO
C     
!$OMP parallel do private(IPQ,XIJ,YIJ,ZIJ) shared(I)
#ifdef _vpp_
!OCL NOALIAS 
#endif
                  DO IPQ = 1, JLEN
                     XIJ = (XC(I)-XJ(IPQ))*BOXIV1
                     YIJ = (YC(I)-YJ(IPQ))*BOXIV2
                     ZIJ = (ZC(I)-ZJ(IPQ))*BOXIV3
C     
                     VXIJ(IPQ) = (XIJ-ANINT(XIJ))*BOX1
                     VYIJ(IPQ) = (YIJ-ANINT(YIJ))*BOX2
                     VZIJ(IPQ) = (ZIJ-ANINT(ZIJ))*BOX3
C     
                     VRIJ2(IPQ) = VXIJ(IPQ)*VXIJ(IPQ) + VYIJ(
     $                    IPQ)*VYIJ(IPQ) + VZIJ(IPQ)*(VZIJ(IPQ)
     $                    +COSB2*VXIJ(IPQ))
                  ENDDO
                  IF (LOCTO) THEN
!$OMP parallel do private(IPQ)
                     DO IPQ = 1, JLEN
                        VEELJ(IPQ) = BOXOQ - DABS(VXIJ(IPQ)) - 
     .                      DABS(VYIJ(IPQ)) - DABS(VZIJ(IPQ))
                        IF (VEELJ(IPQ) .LE. 0.D0) THEN
                           VRIJ2(IPQ)=VRIJ2(IPQ)+VEELJ(IPQ)*BOX1
                           VXIJ(IPQ) =VXIJ(IPQ)-SIGN(BOXOH,VXIJ(IPQ))
                           VYIJ(IPQ) =VYIJ(IPQ)-SIGN(BOXOH,VYIJ(IPQ))
                           VZIJ(IPQ) =VZIJ(IPQ)-SIGN(BOXOH,VZIJ(IPQ))
                        ENDIF
                     ENDDO
                  ENDIF
C     
C--   CHECK FOR CLOSE CONTACTS
                  IPQ = ISRCHFLE(JLEN,VRIJ2,1,1.0D-20)
                  IF (IPQ .LE. JLEN) THEN
                     PRINT 9240, NI, I
                     LERROR = .TRUE.
                     GO TO 5750
                  ENDIF
C     
!$OMP parallel do private(IPQ)
#ifdef _vpp_
!OCL NOALIAS 
#endif
                  DO IPQ = 1, JLEN
                     VRIJ(IPQ) = DSQRT(VRIJ2(IPQ))
                     VRINV(IPQ) = 1.0D0/VRIJ(IPQ)
                     VRINV2(IPQ) = VRINV(IPQ)*VRINV(IPQ)
                     VRINV6(IPQ) = VRINV2(IPQ)*VRINV2(IPQ)*
     $                             VRINV2(IPQ)
                  ENDDO
C     
!$OMP parallel do private(IPQ,QIQJA)
#ifdef _vpp_
!OCL NOALIAS 
#endif
                  DO IPQ = 1, JLEN
                     QIQJA = CGAT(I) * CGAJ(IPQ)
                     VERFJ(IPQ) = QIQJA * RFE * VRIJ2(IPQ)
                     VERCJ(IPQ) = QIQJA * RFC
#ifndef EWALD
                     VEELJ(IPQ) = QIQJA*VRINV(IPQ)
                     DVEELJ(IPQ)= VEELJ(IPQ) * VRINV2(IPQ)+
     $                             6.D0*VRINV6(IPQ)*(2.0D0*VRINV6(IPQ)*
     $                             VC12A(IPQ)-VC6A(IPQ)) * VRINV2(IPQ)
#else
                     IF (VRIJ2(IPQ) .LE. RCUTP2) THEN
                       VEELJ(IPQ) = QIQJA* (  VRINV(IPQ)
     $                          + FACG3* VRIJ(IPQ)*VRIJ2(IPQ)
     $                          + FACG2*VRIJ2(IPQ)
     $                          + FACG0 )
                       DVEELJ(IPQ)=  QIQJA* VRINV2(IPQ)*
     $                            (VRINV(IPQ)
     $                          +  FACGF2 * VRIJ(IPQ)*VRIJ2(IPQ)
     $                          +  FACGF1 * VRIJ2(IPQ)) +
     $                             6.D0*VRINV6(IPQ)*(2.0D0*VRINV6(IPQ)*
     $                             VC12A(IPQ)-VC6A(IPQ)) * VRINV2(IPQ) 
                     ELSE
                       VEELJ(IPQ) = 0.D0
                       DVEELJ(IPQ)= 0.D0
                     ENDIF  
#endif
                     VLJJ(IPQ) = VRINV6(IPQ)*(VC12A(IPQ)*VRINV6(IPQ)
     .                         - VC6A(IPQ))
                     DF3(IPQ) =  DVEELJ(IPQ)+ QIQJA * RFF
                  ENDDO
C     
                  
                  JRCTOT = JLEN
                  JRCSRT = 0
 5150             CONTINUE
                  JRCLEN = MIN(IVECT,JRCTOT)

CDIR$ SHORTLOOP                  
!$OMP parallel do private(JRC) shared(JRCSRT,ICE)
#ifdef _vpp_
!OCL NOALIAS 
#endif
                  DO JRC = 1, JRCLEN
                     TEEL(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) = 
     $                    TEEL(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) + 
     $                    VEELJ(JRC+JRCSRT)
                     TELJ(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) = 
     $                    TELJ(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) + 
     $                    VLJJ(JRC+JRCSRT)
                     TERF(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TERF(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VERFJ(JRC+JRCSRT)
                     TERC(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TERC(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VERCJ(JRC+JRCSRT)
                  ENDDO
                  JRCSRT = JRCSRT + JRCLEN
                  JRCTOT = JRCTOT - JRCLEN
                  IF (JRCTOT .GT. 0) GO TO 5150
                  
C     
C*****CALCULATE THE FORCES.
C     
!$OMP parallel do private(IPQ)
#ifdef _vpp_
!OCL NOALIAS      
#endif
                  DO IPQ = 1, JLEN
                     VXH(IPQ) = VXIJ(IPQ)*DF3(IPQ)
                     VYH(IPQ) = VYIJ(IPQ)*DF3(IPQ)
                     VZH(IPQ) = VZIJ(IPQ)*DF3(IPQ)
                     FXTEMP(IPQ) = FXTEMP(IPQ) - VXH(IPQ)
                     FYTEMP(IPQ) = FYTEMP(IPQ) - VYH(IPQ)
                     FZTEMP(IPQ) = FZTEMP(IPQ) - VZH(IPQ)
                  ENDDO
C     
                  IF (LDOVIR) THEN
!$OMP parallel do private(IPQ) shared(I)
                     DO IPQ = 1, JLEN
                        TVIR(IPQ,1) = TVIR(IPQ,1) - VXH(IPQ)*
     $                       (VXIJ(IPQ)+XRJ(IPQ)-XRC(I))
                        TVIR(IPQ,2) = TVIR(IPQ,2) - VYH(IPQ)*
     $                       (VYIJ(IPQ)+YRJ(IPQ)-YRC(I))
                        TVIR(IPQ,3) = TVIR(IPQ,3) - VZH(IPQ)*
     $                       (VZIJ(IPQ)+ZRJ(IPQ)-ZRC(I))
                     ENDDO
                  ENDIF
C     
#ifdef _vpp_
!OCL NOALIAS 
#endif
                  DO IPQ = 1, JLEN
                     FLRX(I) = FLRX(I) + VXH(IPQ)
                     FLRY(I) = FLRY(I) + VYH(IPQ)
                     FLRZ(I) = FLRZ(I) + VZH(IPQ)
                  ENDDO
C     
                  JN = JN + JLEN
 5400          CONTINUE
C     
C     
!$OMP parallel do private(I32X)
#ifdef _vpp_
!OCL NOALIAS 
#endif
               DO I32X = 1, JLEN
                  FLRX(JLIST(I32X)) = FXTEMP(I32X)
                  FLRY(JLIST(I32X)) = FYTEMP(I32X)
                  FLRZ(JLIST(I32X)) = FZTEMP(I32X)
               ENDDO
C     
            ENDIF
 5500    CONTINUE
C     
C*****SCATTER THE FORCES BACK TO THEIR ORIGINAL POSITIONS.
C     
!$OMP parallel do private(JRC)
#ifdef _vpp_
!OCL NOVREC(F)
#endif
         DO JRC = 1, NR
            F((JRC-1)*3+1) = F((JRC-1)*3+1) + FLRX(JRC)
            F((JRC-1)*3+2) = F((JRC-1)*3+2) + FLRY(JRC)
            F((JRC-1)*3+3) = F((JRC-1)*3+3) + FLRZ(JRC)
         ENDDO
C     
         IF (LDOVIR) THEN
            DO JRC = 1, MAXLEN
               VIR(1) = VIR(1) + TVIR(JRC,1)
               VIR(2) = VIR(2) + TVIR(JRC,2)
               VIR(3) = VIR(3) + TVIR(JRC,3)
            ENDDO
            VIR(1) = VIR(1)*0.5D0
            VIR(2) = VIR(2)*0.5D0
            VIR(3) = VIR(3)*0.5D0
         ENDIF
C     
C--   UPDATE THE 'EEL' AND 'ELJ' VALUES
         DO JRC1 = 1, NUNRE2
            DO JRC2 = 1, IVECT
               EEL(JRC1) = EEL(JRC1) + TEEL(JRC1, JRC2)
#ifndef EWALD
               ERF(JRC1) = ERF(JRC1) + TERF(JRC1, JRC2)
               ERC(JRC1) = ERC(JRC1) + TERC(JRC1, JRC2)
#endif
               ELJ(JRC1) = ELJ(JRC1) + TELJ(JRC1, JRC2)
               
            ENDDO
         ENDDO

      ELSE

C*******************************************************************************
C 4D
C*******************************************************************************
         
C     
C*****CALCULATE THE NON-BONDED INTERACTIONS IN 4D.
C     
C     
C*****TREAT THE INTRA-CHARGE GROUP INTERACTIONS.
C*****IT IS ASSUMED THAT CHARGE GROUPS CONTAINING LESS THAN 4 ATOMS
C*****CONTAIN NO INTRA CHARGE GROUP INTERACTIONS.
C*****SOLVENT MOLECULES ARE ASSUMED TO CONTAIN NO INTRA MOLECULE
C*****INTERACTIONS
C     
         CALL WHENIGT (NCAG, INAG(1),1,3,NILIST,NITOT)
C     
         DO NN = 1, NITOT
            NI = NILIST(NN)
            ILAST = IFIRST(NI) + INAG(NI) - 1
C     
C     
            DO I = IFIRST(NI), ILAST - 1
               DO J = I + 1, ILAST
C     
C*****CHECK WHETHER ATOM J IS EXCLUDED FROM ATOM I SOMEHOW.
C     
                  ILEX = KNE(I) + INE(I)
                  JRCLOC = KNE(I) + 1
                  JRCLEN = ILEX - JRCLOC + 1
                  IF(JRCLEN.GT.0) THEN
                     JRCEXL=ISRCHEQ(JRCLEN,JSNE(JRCLOC),1,J)-JRCLEN
                  ELSE
                     JRCEXL=1
                  ENDIF

C                  IF (LPION) THEN
C                     LIZERO = (IPICT(I) .EQ. 0)
C                     LJZERO = (IPICT(J) .EQ. 0)
C                     IF (LIZERO .AND. LJZERO) THEN
C                        PIFAC = 1.0
C                     ELSEIF (LIZERO .OR. LJZERO) THEN
C                        PIFAC = PININV
C                     ELSEIF (IPICT(I) .EQ. IPICT(J)) THEN
C                        PIFAC = PININV
C                     ELSE
C                        PIFAC = 0.0
C                        JRCEXL = 0
C                     ENDIF
C                  ELSE
C                     PIFAC = 1.0
C                  ENDIF 

                  IF(JRCEXL.GT.0) THEN
C     
C*****CHECK 1-4 EXCLUSIONS.
C     
                     ILEX = KNE14(I) + INE14(I)
                     JRCLOC = KNE14(I) + 1
                     JRCLEN = ILEX - JRCLOC + 1
                     IF(JRCLEN.GT.0) THEN
                        LIEX14 = ISRCHEQ(JRCLEN,JSNE14(JRCLOC),1,J)
     $                       .LE. JRCLEN
                     ELSE
                        LIEX14 = .FALSE.
                     ENDIF
C     
                     IPC = MPAC(IAC(I),IAC(J))
                     IF (LIEX14) THEN
                        CA12 = CS12(IPC)
                        CA6 = CS6(IPC)
                     ELSE
                        CA12 = C12(IPC)
                        CA6 = C6(IPC)
                     ENDIF
C     
                     IF (LPERTL) THEN
                        IPCB = MPAC(IACBT(I),IACBT(J))
                        IF (LIEX14) THEN
                           CB12 = CS12(IPCB)
                           CB6 = CS6(IPCB)
                        ELSE
                           CB12 = C12(IPCB)
                           CB6 = C6(IPCB)
                        ENDIF
                        IF (DABS(CA6) .LT. ZEROC6) THEN
                           CA126 = 0.0D0
                        ELSE
                           CA126 = CA12 / CA6
                        ENDIF
                        IF (DABS(CB6) .LT. ZEROC6) THEN
                           CB126 = 0.0D0
                        ELSE
                           CB126 = CB12 / CB6
                        ENDIF
                     ENDIF

C     
C*****CALCULATE THE LJ-ENERGY AND FORCE.
C     

                     IBASE = 0
                     IF (LOCTO) THEN
                        DO K = 1, NPM
                           IX = I + IBASE
                           JX = J + IBASE
                           DX = (XC(IX)-XC(JX))*BOXIV1
                           DX = (DX-ANINT(DX))*BOX1
                           DY = (YC(IX)-YC(JX))*BOXIV2
                           DY = (DY-ANINT(DY))*BOX2
                           DZ = (ZC(IX)-ZC(JX))*BOXIV3
                           DZ = (DZ-ANINT(DZ))*BOX3
                           DW = (WC(IX)-WC(JX))*BOXIV4
                           DW = (DW-ANINT(DW))*BOX4
                           RIJ3D2 = DX*DX + DY*DY + DZ*(DZ+COSB2*DX)
                           RIJ2 = RIJ3D2 + DW * DW
                           DRIJ = BOXOQ - DABS(DX)-DABS(DY)-DABS(DZ)
                           IF (DRIJ .LE. 0.D0) THEN
                              RIJ3D2 = RIJ3D2 + DRIJ*BOX1
                              RIJ2 = RIJ2 + DRIJ*BOX1
                              DX = DX - SIGN(BOXOH,DX)
                              DY = DY - SIGN(BOXOH,DY)
                              DZ = DZ - SIGN(BOXOH,DZ)
                           ENDIF
                           ILJ = NRELKP(ICET(IX),ICET(JX))
                           QIQJA = CGAT(IX) * CGAT(JX)
!                          if(LIEX14) QIQJA=QIQJA* scale_14_amber
                           IF (LLPERT(IX) .OR. LLPERT(JX)) THEN
                           QIQJB = CGBT(IX) * CGBT(JX)
                           SOFTLJ = ALPHLJ * 
     $                          (1 - ISCLJT(IX) * ISCLJT(JX))
                           SOFTEL = ALPHC * 
     $                          (1 - ISCCT(IX) * ISCCT(JX))

C*******************************************************************************
C     4D
C*******************************************************************************
C
C---- SOFT-CORE PERTURBATION
C     CALCULATE FULL 4D CONTRIBUTION
C     
                           RIJ4 = RIJ2 * RIJ2
C     
C--   4D, STATE A
C     
C     COULOMB
                           SFEL   = SOFTEL * RL2A
                           RSOFT2 = 1.0D0 / (SFEL + RIJ2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELA4D = QIQJA * RSOFT                    
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFA4 = QIQJA * C1 * RRF
                           QRFR   = DFRFA4 * RIJ3D2
                           VRFA4D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CA126
                           CRA1   = 1.0D0 / (SFLJ * RL2A + RIJ6)
                           CRA12  = CA12 * CRA1
                           VLJA4D = CRA1 * (CRA12 - CA6)
                           DFA4D  = CRA1 * ( CRA1 * CRA12 + VLJA4D)
C     DV/DL
                           DLLJA4 = -2.0D0 * SFLJ * DFA4D
                           DLELA4 = -VELA4D * SOFTEL * RSOFT2
                           DLRFA4 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFA4D = RSOFT2*VELA4D+6.0D0*RIJ4*DFA4D
C     
C--   4D, STATE B
C     
C     COULOMB
                           SFEL   = SOFTEL * RL2B 
                           RSOFT2 = 1.0D0 / (SFEL + RIJ2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELB4D = QIQJB * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFB4 = QIQJB * C1 * RRF
                           QRFR   = DFRFB4 * RIJ3D2
                           VRFB4D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CB126
                           CRA1   = 1.0D0 / (SFLJ * RL2B + RIJ6)
                           CRA12  = CB12 * CRA1
                           VLJB4D = CRA1 * (CRA12 - CB6)
                           DFB4D  = CRA1 * ( CRA1 * CRA12 + VLJB4D)
C     DV/DL
                           DLLJB4 = -2.0D0 * SFLJ * DFB4D
                           DLELB4 = -VELB4D * SOFTEL * RSOFT2
                           DLRFB4 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFB4D = RSOFT2*VELB4D+6.0D0*RIJ4*DFB4D
C     
C--   COMBINE 4D STATE A AND B
C     
                           VEL4D  = RLB * VELB4D + RLA * VELA4D
                           VRF4D  = RLB * VRFB4D + RLA * VRFA4D
                           VLJ4D  = RLB * VLJB4D + RLA * VLJA4D
                           DF4D4D = RLB * DFB4D  + RLA * DFA4D
                           DF4D   = DF4D4D 
     .                            + RLB  * DFRFB4 + RLA  * DFRFA4
                           DLLJ4D = RDLB * VLJB4D + RLLB * DLLJB4
     .                            + RDLA * VLJA4D + RLLA * DLLJA4
                           DLEL4D = RDLB * VELB4D + RLLB * DLELB4
     .                            + RDLA * VELA4D + RLLA * DLELA4
                           DLRF4D = RDLB * VRFB4D + RLLB * DLRFB4
     .                            + RDLA * VRFA4D + RLLA * DLRFA4

C*******************************************************************************
C     3D
C*******************************************************************************
                           RIJ4 = RIJ3D2 * RIJ3D2
                           RIJ6 = RIJ4 * RIJ3D2
C     
C--   3D, STATE A
C     
C     COULOMB
                           SFEL   = SOFTEL * RLMA
                           RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELA3D = QIQJA * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFA3 = QIQJA * C1 * RRF
                           QRFR   = DFRFA3 * RIJ3D2
                           VRFA3D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CA126
                           CRA1   = 1.0D0 / (SFLJ * RLMA + RIJ6)
                           CRA12  = CA12 * CRA1
                           VLJA3D = CRA1 * (CRA12 - CA6)
                           DFA3D  = CRA1 * ( CRA1 * CRA12 + VLJA3D)
C     DV/DL
                           DLLJA3 = -2.0D0 * SFLJ * DFA3D
                           DLELA3 = -VELA3D * SOFTEL * RSOFT2 
                           DLRFA3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFA3D = RSOFT2*VELA3D+6.0D0*RIJ4*DFA3D
C     
C--   3D, STATE B
C     
C     COULOMB
                           SFEL   = SOFTEL * RLMB
                           RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELB3D = QIQJB * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFB3 = QIQJB * C1 * RRF
                           QRFR   = DFRFB3 * RIJ3D2
                           VRFB3D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CB126 
                           CRA1   = 1.0D0 / (SFLJ * RLMB + RIJ6)
                           CRA12  = CB12 * CRA1
                           VLJB3D = CRA1 * (CRA12 - CB6)
                           DFB3D  = CRA1 * ( CRA1 * CRA12 + VLJB3D)
C     DV/DL
                           DLLJB3 = -2.0D0 * SFLJ * DFB3D
                           DLELB3 = -VELB3D * SOFTEL * RSOFT2 
                           DLRFB3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFB3D = RSOFT2*VELB3D+6.0D0*RIJ4*DFB3D
C     
C--   COMBINE 3D, STATE A AND B
C     
                           VEL3D =  RLB  * VELB3D + RLA  * VELA3D
                           VRF3D =  RLB  * VRFB3D + RLA  * VRFA3D
                           VLJ3D =  RLB  * VLJB3D + RLA  * VLJA3D
                           DF3D  =  RLB  * DFB3D  + RLA  * DFA3D
     $                            + RLB  * DFRFB3 + RLA  * DFRFA3
                           DLLJ3D = RDLB * VLJB3D + RLLB * DLLJB3
     $                            + RDLA * VLJA3D + RLLA * DLLJA3
                           DLEL3D = RDLB * VELB3D + RLLB * DLELB3 
     $                            + RDLA * VELA3D + RLLA * DLELA3
                           DLRF3D = RDLB * VRFB3D + RLLB * DLRFB3
     $                            + RDLA * VRFA3D + RLLA * DLRFA3

                           DDM3EL =  RLBRMU * DLELB3 + RLARMU * DLELA3
                           DDM3RF =  RLBRMU * DLRFB3 + RLARMU * DLRFA3
                           DDM3LJ =  RLBRMU * DLLJB3 + RLARMU * DLLJA3

C     REACTION FIELD CONSTANT
                           VRCA = RFC*QIQJA
                           VRCB = RFC*QIQJB
                           VRC  = RLB*VRCB  + RLA*VRCA
                           DLRC = RDLB*VRCB + RDLA*VRCA
C     
C---- COMBINE 4D AND 3D COMPONENTS
                           DF3T = RM3D * DF3D + RM4D * DF4D
                           DF4T = RM4D * DF4D4D

                           EEL(ILJ) = EEL(ILJ) +  
     $                          RM3D * VEL3D + RM4D * VEL4D
#ifndef EWALD
                           ERF(ILJ) = ERF(ILJ) + 
     $                          RM3D * VRF3D + RM4D * VRF4D
                           ERC(ILJ) = ERC(ILJ) + 
     $                          RM3D * VRC   + RM4D * VRC
#endif
                           ELJ(ILJ) = ELJ(ILJ) + 
     $                          RM3D * VLJ3D + RM4D * VLJ4D
                           
                           EGLEL = EGLEL + 
     $                          RM3D * DLEL3D + RM4D * DLEL4D
                           EGLRF = EGLRF + 
     $                          RM3D * DLRF3D + RM4D * DLRF4D
                           EGLRC = EGLRC + 
     $                          RM3D * DLRC   + RM4D * DLRC
                           EGLLJ = EGLLJ +
     $                          RM3D * DLLJ3D + RM4D * DLLJ4D

                           EG34EL = EG34EL + 
     $                          RDM3D*VEL3D + RM3D*DDM3EL + RDM4D*VEL4D
                           EG34RF = EG34RF +
     $                          RDM3D*VRF3D + RM3D*DDM3RF + RDM4D*VRF4D
                           EG34RC = EG34RC + 
     $                          RDM3D*VRC + RDM4D*VRC
                           EG34LJ = EG34LJ +
     $                          RDM3D*VLJ3D + RM3D*DDM3LJ + RDM4D*VLJ4D
C     
                        ELSE
                           RIJIN2 = 1.0D0 / RIJ2
                           RIJINV = DSQRT(RIJIN2)
                           RIJIN6 = RIJIN2 * RIJIN2 * RIJIN2
C     
                           VELA3D = QIQJA * RIJINV
C     
                           CRA1   = CA12 * RIJIN6
                           CRA12  = CRA1 - CA6
                           VLJA3D = CRA12 * RIJIN6
C     
                           DF4T  =  RIJIN2 * (VELA3D + 
     $                          6.0D0 * (CRA1 + CRA12) * RIJIN6)
                           DF3T = DF4T + RFF * QIQJA
C     
#ifndef EWALD
                           ERF(ILJ) = ERF(ILJ) +
     $                          QIQJA * RFE * RIJ3D2
                           ERC(ILJ) = ERC(ILJ) +
     $                          QIQJA * RFC 
#endif
                           ELJ(ILJ) = ELJ(ILJ) + VLJA3D
                           EEL(ILJ) = EEL(ILJ) + VELA3D
                        ENDIF
                           
                        XH = DX*DF3T
                        YH = DY*DF3T
                        ZH = DZ*DF3T
                        WH = DW*DF4T
                        FLRX(IX) = FLRX(IX) + XH
                        FLRY(IX) = FLRY(IX) + YH
                        FLRZ(IX) = FLRZ(IX) + ZH
                        FLRW(IX) = FLRW(IX) + WH
                        FLRX(JX) = FLRX(JX) - XH
                        FLRY(JX) = FLRY(JX) - YH
                        FLRZ(JX) = FLRZ(JX) - ZH
                        FLRW(JX) = FLRW(JX) - WH
                        
                        IBASE = IBASE + NRP
                     END DO
                  ELSE
                     DO K = 1, NPM
                           IX = I + IBASE
                           JX = J + IBASE
                           DX = (XC(IX)-XC(JX))*BOXIV1
                           DX = (DX-ANINT(DX))*BOX1
                           DY = (YC(IX)-YC(JX))*BOXIV2
                           DY = (DY-ANINT(DY))*BOX2
                           DZ = (ZC(IX)-ZC(JX))*BOXIV3
                           DZ = (DZ-ANINT(DZ))*BOX3
                           DW = (WC(IX)-WC(JX))*BOXIV4
                           DW = (DW-ANINT(DW))*BOX4
                           RIJ3D2 = DX*DX + DY*DY + DZ*(DZ+COSB2*DX)
                           RIJ2 = RIJ3D2 + DW*DW
                           ILJ = NRELKP(ICET(IX),ICET(JX))
                           QIQJA = CGAT(IX) * CGAT(JX)
                           IF (LLPERT(IX) .OR. LLPERT(JX)) THEN
                           QIQJB = CGBT(IX) * CGBT(JX)
                           SOFTLJ = ALPHLJ * 
     $                          (1 - ISCLJT(IX) * ISCLJT(JX))
                           SOFTEL = ALPHC * 
     $                          (1 - ISCCT(IX) * ISCCT(JX))

C*******************************************************************************
C     4D
C*******************************************************************************
                           
C     
C---- SOFT-CORE PERTURBATION
C     CALCULATE FULL 4D CONTRIBUTION
C     
                           RIJ4 = RIJ2 * RIJ2
                           RIJ6 = RIJ4 * RIJ2
C     
C--   4D, STATE A
C     
C     COULOMB
                           SFEL   = SOFTEL * RL2A
                           RSOFT2 = 1.0D0 / (SFEL + RIJ2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELA4D = QIQJA * RSOFT                    
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFA4 = QIQJA * C1 * RRF
                           QRFR   = DFRFA4 * RIJ3D2
                           VRFA4D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CA126
                           CRA1   = 1.0D0 / (SFLJ * RL2A + RIJ6)
                           CRA12  = CA12 * CRA1
                           VLJA4D = CRA1 * (CRA12 - CA6)
                           DFA4D  = CRA1 * ( CRA1 * CRA12 + VLJA4D)
C     DV/DL
                           DLLJA4 = -2.0D0 * SFLJ * DFA4D
                           DLELA4 = -VELA4D * SOFTEL * RSOFT2
                           DLRFA4 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFA4D = RSOFT2*VELA4D+6.0D0*RIJ4*DFA4D
C     
C--   4D, STATE B
C     
C     COULOMB
                           SFEL   = SOFTEL * RL2B 
                           RSOFT2 = 1.0D0 / (SFEL + RIJ2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELB4D = QIQJB * RSOFT
C     REACTION FIELD
                           RRF2   =  1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFB4 = QIQJB * C1 * RRF
                           QRFR   = DFRFB4 * RIJ3D2
                           VRFB4D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CB126
                           CRA1   = 1.0D0 / (SFLJ * RL2B + RIJ6)
                           CRA12  = CB12 * CRA1
                           VLJB4D = CRA1 * (CRA12 - CB6)
                           DFB4D  = CRA1 * ( CRA1 * CRA12 + VLJB4D)
C     DV/DL
                           DLLJB4 = -2.0D0 * SFLJ * DFB4D
                           DLELB4 = -VELB4D * SOFTEL * RSOFT2
                           DLRFB4 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFB4D = RSOFT2*VELB4D+6.0D0*RIJ4*DFB4D
C     
C--   COMBINE 4D STATE A AND B
C     
                           VEL4D  =    RLB  * VELB4D + RLA  * VELA4D
                           VRF4D  =    RLB  * VRFB4D + RLA  * VRFA4D
                           VLJ4D  =    RLB  * VLJB4D + RLA  * VLJA4D
                           DF4D4D =    RLB  * DFB4D  + RLA  * DFA4D
                           DF4D   =    DF4D4D 
     $                               + RLB  * DFRFB4 + RLA  * DFRFA4
                           DLLJ4D =    RDLB * VLJB4D + RLLB * DLLJB4
     $                               + RDLA * VLJA4D + RLLA * DLLJA4
                           DLEL4D =    RDLB * VELB4D + RLLB * DLELB4
     $                               + RDLA * VELA4D + RLLA * DLELA4
                           DLRF4D =    RDLB * VRFB4D + RLLB * DLRFB4
     $                               + RDLA * VRFA4D + RLLA * DLRFA4

C*******************************************************************************
C     3D
C*******************************************************************************

                           RIJ4 = RIJ3D2 * RIJ3D2
                           RIJ6 = RIJ4 * RIJ3D2
C
C--   3D, STATE A
C     
C     COULOMB
                           SFEL   = SOFTEL * RLMA
                           RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELA3D = QIQJA * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFA3 = QIQJA * C1 * RRF
                           QRFR   = DFRFA3 * RIJ3D2
                           VRFA3D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CA126
                           CRA1   = 1.0D0 / (SFLJ * RLMA + RIJ6)
                           CRA12  = CA12 * CRA1
                           VLJA3D = CRA1 * (CRA12 - CA6)
                           DFA3D  = CRA1 * ( CRA1 * CRA12 + VLJA3D)
C     DV/DL
                           DLLJA3 = -2.0D0 * SFLJ * DFA3D
                           DLELA3 = -VELA3D * SOFTEL * RSOFT2 
                           DLRFA3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFA3D = RSOFT2*VELA3D+6.0D0*RIJ4*DFA3D
C     
C--   3D, STATE B
C     
C     COULOMB
                           SFEL   = SOFTEL * RLMB
                           RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELB3D = QIQJB * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFB3 = QIQJB * C1 * RRF
                           QRFR   = DFRFB3 * RIJ3D2
                           VRFB3D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CB126 
                           CRA1   = 1.0D0 / (SFLJ * RLMB + RIJ6)
                           CRA12  = CB12 * CRA1
                           VLJB3D = CRA1 * (CRA12 - CB6)
                           DFB3D  = CRA1 * ( CRA1 * CRA12 + VLJB3D)
C     DV/DL
                           DLLJB3 = -2.0D0 * SFLJ * DFB3D
                           DLELB3 = -VELB3D * SOFTEL * RSOFT2 
                           DLRFB3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFB3D = RSOFT2*VELB3D+6.0D0*RIJ4*DFB3D
C     
C--   COMBINE 3D, STATE A AND B
C     
                           VEL3D =  RLB  * VELB3D + RLA  * VELA3D
                           VRF3D =  RLB  * VRFB3D + RLA  * VRFA3D
                           VLJ3D =  RLB  * VLJB3D + RLA  * VLJA3D
                           DF3D  =  RLB  * DFB3D  + RLA  * DFA3D
     $                            + RLB  * DFRFB3 + RLA  * DFRFA3
                           DLLJ3D = RDLB * VLJB3D + RLLB * DLLJB3
     $                            + RDLA * VLJA3D + RLLA * DLLJA3
                           DLEL3D = RDLB * VELB3D + RLLB * DLELB3 
     $                            + RDLA * VELA3D + RLLA * DLELA3
                           DLRF3D = RDLB * VRFB3D + RLLB * DLRFB3
     $                            + RDLA * VRFA3D + RLLA * DLRFA3

                           DDM3EL =  RLBRMU * DLELB3 + RLARMU * DLELA3
                           DDM3RF =  RLBRMU * DLRFB3 + RLARMU * DLRFA3
                           DDM3LJ =  RLBRMU * DLLJB3 + RLARMU * DLLJA3  

C     REACTION FIELD CONSTANT
                           VRCA = RFC*QIQJA
                           VRCB = RFC*QIQJB
                           VRC  = RLB*VRCB  + RLA*VRCA
                           DLRC = RDLB*VRCB + RDLA*VRCA
C     
C---- COMBINE 4D AND 3D COMPONENTS
                           DF3T = RM3D * DF3D + RM4D * DF4D
                           DF4T = RM4D * DF4D4D
                           EEL(ILJ) = EEL(ILJ) +  
     $                          RM3D * VEL3D + RM4D * VEL4D
#ifndef EWALD
                           ERF(ILJ) = ERF(ILJ) + 
     $                          RM3D * VRF3D + RM4D * VRF4D
                           ERC(ILJ) = ERC(ILJ) + 
     $                          RM3D * VRC   + RM4D * VRC
#endif
                           ELJ(ILJ) = ELJ(ILJ) + 
     $                          RM3D * VLJ3D + RM4D * VLJ4D
                           
                           EGLEL = EGLEL + 
     $                          RM3D * DLEL3D + RM4D * DLEL4D
                           EGLRF = EGLRF + 
     $                          RM3D * DLRF3D + RM4D * DLRF4D
                           EGLRC = EGLRC + 
     $                          RM3D * DLRC + RM4D * DLRC
                           EGLLJ = EGLLJ +
     $                          RM3D * DLLJ3D + RM4D * DLLJ4D

                           EG34EL = EG34EL + 
     $                          RDM3D*VEL3D + RM3D*DDM3EL + RDM4D*VEL4D
                           EG34RF = EG34RF +
     $                          RDM3D*VRF3D + RM3D*DDM3RF + RDM4D*VRF4D
                           EG34RC = EG34RC + 
     $                          RDM3D*VRC + RDM4D*VRC
                           EG34LJ = EG34LJ +
     $                          RDM3D*VLJ3D + RM3D*DDM3LJ + RDM4D*VLJ4D
C       
                        ELSE
                           RIJIN2 = 1.0D0 / RIJ2
                           RIJINV = DSQRT(RIJIN2)
                           RIJIN6 = RIJIN2 * RIJIN2 * RIJIN2
C     
                           VELA3D = QIQJA * RIJINV
C     
                           CRA1   = CA12 * RIJIN6
                           CRA12  = CRA1 - CA6
                           VLJA3D = CRA12 * RIJIN6
C     
                           DF4T  =  RIJIN2 * (VELA3D + 
     $                          6.0D0*(CRA1 + CRA12) * RIJIN6)
                           DF3T = DF4T + RFF * QIQJA
C     
#ifndef EWALD
                           ERF(ILJ) = ERF(ILJ) +
     $                          QIQJA * RFE * RIJ3D2
                           ERC(ILJ) = ERC(ILJ) +
     $                          QIQJA * RFC
#endif
                           ELJ(ILJ) = ELJ(ILJ) + VLJA3D
                           EEL(ILJ) = EEL(ILJ) + VELA3D
                        ENDIF
                           
                        XH = DX*DF3T
                        YH = DY*DF3T 
                        ZH = DZ*DF3T
                        WH = DW*DF4T
                        FLRX(IX) = FLRX(IX) + XH
                        FLRY(IX) = FLRY(IX) + YH
                        FLRZ(IX) = FLRZ(IX) + ZH
                        FLRW(IX) = FLRW(IX) + WH
                        FLRX(JX) = FLRX(JX) - XH
                        FLRY(JX) = FLRY(JX) - YH
                        FLRZ(JX) = FLRZ(JX) - ZH
                        FLRW(JX) = FLRW(JX) - WH
                        IBASE = IBASE + NRP
                     END DO
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C     
C     
C*****TREAT THE FORCES OF ALL PAIRS DENOTED BY THE PAIRLIST.
C*****FIRST, ALL CHARGE GROUPS NI BELONGING TO THE SOLUTES WILL BE PROC
C     
C     
      JN = 1
      NJNB = 1
      DO NI = 1, NRAG
         NPAIR = INB(NI)
C     
         IF (NPAIR .GT. 0) THEN
            DO IPQ = 1, NPAIR
               JFIRST(IPQ) = IFIRST(JNB(IPQ+NJNB-1))
               NJL(IPQ) = INAG(JNB(IPQ+NJNB-1))
            ENDDO
C*****CODE EXPANDED FROM ROUTINE:  VEXPLD
C     
               JRCNUM = 1
               JRCL1X = NPAIR
               J1X = 0
88001          CONTINUE
               K1X = 0

               DO I1X = 1, JRCL1X
                  IF (NJL(I1X) .GE. JRCNUM) THEN
                     J1X = J1X + 1
                     JLIST(J1X) = JFIRST(I1X)
                     K1X = K1X + 1
                     JFIRST(K1X) = JFIRST(I1X) + 1
                     NJL(K1X) = NJL(I1X)
                  ENDIF
               END DO
               JRCL1X = K1X
               JRCNUM = JRCNUM + 1
               IF (JRCL1X .GT. 0) GO TO 88001
               JLEN = J1X
C*****END OF CODE EXPANDED FROM ROUTINE:  VEXPLD
               NJNB = NJNB + NPAIR
C     
               ILAST = IFIRST(NI) + INAG(NI) - 1
C     
C*****GATHER THE COORDINATES AND THE FORCES.
C     
               IF (LPERTL) THEN
!$OMP parallel do private(IPQ)
                  DO IPQ = 1, JLEN
                     XJ(IPQ) = XC(JLIST(IPQ))
                     YJ(IPQ) = YC(JLIST(IPQ))
                     ZJ(IPQ) = ZC(JLIST(IPQ))
                     WJ(IPQ) = WC(JLIST(IPQ))
C     
                     FXTEMP(IPQ) = FLRX(JLIST(IPQ))
                     FYTEMP(IPQ) = FLRY(JLIST(IPQ))
                     FZTEMP(IPQ) = FLRZ(JLIST(IPQ))
                     FWTEMP(IPQ) = FLRW(JLIST(IPQ))
C     
C*****DETERMINE THE INTEGER ATOM CODE FOR EACH ATOM J.
C*****GATHER THE CHARGES AND DETERMINE TO WHICH GROUP EACH ATOM BELO
C     
                     IACAJ(IPQ) = IACT(JLIST(IPQ)) + JRCBK1(IPQ)
                     CGAJ(IPQ) = CGAT(JLIST(IPQ))
                     JRCILJ(IPQ) = ICET(JLIST(IPQ))
                     AMULT(IPQ) = 1.D0
                     IFORCE(IPQ) = 0
                     IACBJ(IPQ) = IACBT(JLIST(IPQ))
                     CGBJ(IPQ) = CGBT(JLIST(IPQ))
                     LLPERJ(IPQ) = LLPERT(JLIST(IPQ))
                     ISCLJJ(IPQ) = ISCLJT(JLIST(IPQ))
                     ISCCJ(IPQ) = ISCCT(JLIST(IPQ))
C                     IF (LDOVIR) NSPJ(IPQ) = NSPTT(JLIST(IPQ))
                  ENDDO
               ELSE
C                  IF (LPION) THEN
C                     DO IPQ = 1, JLEN
C                        IPICJ(IPQ) = IPICT(JLIST(IPQ))
C                     ENDDO
C                  ENDIF
C                  IF (LDOVIR) THEN
C                     DO IPQ = 1, JLEN
C                        XJ(IPQ) = XC(JLIST(IPQ))
C                        YJ(IPQ) = YC(JLIST(IPQ))
C                        ZJ(IPQ) = ZC(JLIST(IPQ))
C                        WJ(IPQ) = WC(JLIST(IPQ))
C     
C                        FXTEMP(IPQ) = FLRX(JLIST(IPQ))
C                        FYTEMP(IPQ) = FLRY(JLIST(IPQ))
C                        FZTEMP(IPQ) = FLRZ(JLIST(IPQ))
C                        FWTEMP(IPQ) = FLRW(JLIST(IPQ))
C     
C*****DETERMINE THE INTEGER ATOM CODE FOR EACH ATOM J.
C*****GATHER THE CHARGES AND DETERMINE TO WHICH GROUP EACH ATOM BELO
C     
C                        IACAJ(IPQ) = IACT(JLIST(IPQ)) + JRCBK1(IPQ)
C                        CGAJ(IPQ) = CGAT(JLIST(IPQ))
C                        JRCILJ(IPQ) = ICET(JLIST(IPQ))
C                        AMULT(IPQ) = 1.
C                        IFORCE(IPQ) = 0
C                        NSPJ(IPQ) = NSPTT(JLIST(IPQ))
C                     END DO
C                  ELSE
!$OMP parallel do private(IPQ)
                     DO IPQ = 1, JLEN
                        XJ(IPQ) = XC(JLIST(IPQ))
                        YJ(IPQ) = YC(JLIST(IPQ))
                        ZJ(IPQ) = ZC(JLIST(IPQ))
                        WJ(IPQ) = WC(JLIST(IPQ))
C     
                        FXTEMP(IPQ) = FLRX(JLIST(IPQ))
                        FYTEMP(IPQ) = FLRY(JLIST(IPQ))
                        FZTEMP(IPQ) = FLRZ(JLIST(IPQ))
                        FWTEMP(IPQ) = FLRW(JLIST(IPQ))
C     
C*****DETERMINE THE INTEGER ATOM CODE FOR EACH ATOM J.
C*****GATHER THE CHARGES AND DETERMINE TO WHICH GROUP EACH ATOM BELO
C     
                        IACAJ(IPQ) = IACT(JLIST(IPQ)) + JRCBK1(IPQ)
                        CGAJ(IPQ) = CGAT(JLIST(IPQ))
                        JRCILJ(IPQ) = ICET(JLIST(IPQ))
                        AMULT(IPQ) = 1.0D0
                        IFORCE(IPQ) = 0
                     ENDDO
C                  ENDIF
               ENDIF
C     
C*****FOR EACH ATOM I MUST BE DETERMINED WHETHER IT MUST BE EXCLUDED
C*****SOMEHOW FROM ATOM I. SOME DATA FOR THIS DETERMINATION CAN BE
C*****GATHERED FOR ALL ATOMS J.
C     
               IBASE = INT((NI-1)/NCAG)*NRP
               JLEN2 = 0
               DO IPQ = 1, JLEN
                  IF (JLIST(IPQ).LE.IBASE+NRP .AND. JLIST(IPQ).GT.
     $                 IBASE) THEN
                     JLEN2 = JLEN2 + 1
                     JRCIND(JLEN2) = IPQ
                     JFIRST(JLEN2) = JLIST(IPQ)
                  ENDIF
               ENDDO
               
               IF (JLEN2 .GT. 0) II = IFIRST(NI) - IBASE
C     
C               IF (LDOVIR) THEN
C                  IVIRL = 0
C                  DO IPQ = 1, JLEN
C                     IF (NSPJ(IPQ) .NE. NSPTT(ILAST)) THEN
C                        IVIRL = IVIRL + 1
C                        IOCVIR(IVIRL) = IPQ
C                     ENDIF
C                  ENDDO
C                  LVIR = .FALSE.
C                  IF (IVIRL .GT. 0) THEN
C                     LVIR = .TRUE.
C                     DO I9X = 1, JLEN
C                        XRJ(I9X) = XRC(JLIST(I9X))
C                        YRJ(I9X) = YRC(JLIST(I9X))
C                        ZRJ(I9X) = ZRC(JLIST(I9X))
C                     ENDDO
C                  ENDIF
C               ENDIF
C     
C     USE 'IFORCE' TO REPLACE THE EXCLUDED PAIR-PARAMETERS BY 1-4
C     PARAMETERS.
C     
C     USE 'AMULT' TO CONTROL THE ASSIGNMENTS TO THE INVERSE DISTANCE.
C     
               DO I = IFIRST(NI), ILAST

                  ICE = ICET(I)

C--   CHECK FOR EXCLUDED ATOMS IN THIS PASS
                  IF (JLEN2 .GT. 0) THEN
C--   SET UP THE DEFAULT ARRAYS FOR EXCLUSIONS AND 1-4 INTERACTIONS
!$OMP parallel do private(K) shared(IBASE)
                     DO K = 1, JLEN2
                        JFORCE(JFIRST(K)-IBASE) = 0
                        BMULT(JFIRST(K)-IBASE) = 1.0D0
                     ENDDO
C     
                     ILEN = INE14(II)
!$OMP parallel do private(K) shared(II)
                     DO K = 1, ILEN
                        JFORCE(JSNE14(K+KNE14(II))) = KFORCE
                     ENDDO
                     ILEN = INE(II)
!$OMP parallel do private(K) shared(II)
                     DO K = 1, ILEN
                        BMULT(JSNE(K+KNE(II))) = 0.D0
                     ENDDO
C     
!$OMP parallel do private(K) shared(IBASE)
                     DO K = 1, JLEN2
                        IFORCE(JRCIND(K)) = JFORCE(JFIRST(K)-IBASE)
                        AMULT(JRCIND(K)) = BMULT(JFIRST(K)-IBASE)
                     ENDDO
                     
                     II = II + 1
                  ENDIF
C     
C--   DETERMINE LJ-PARAMETERS.
C     
                  IF (LPERTL) THEN
!$OMP parallel do private(IPQ,NJIND,NJINDB) shared(I)
                     DO IPQ = 1, JLEN
                        NJIND = JRCPAC(IACAJ(IPQ)+(IACT(I)-1)*NRATT) + 
     $                       IFORCE(IPQ) + JRCBK2(IPQ)
                        NJINDB = JRCPAC(IACBJ(IPQ)+(IACBT(I)-1)*NRATT) 
     $                       + IFORCE(IPQ) + JRCBK2(IPQ)
                        VC6A(IPQ) = FORC6(NJIND)
                        VC12A(IPQ) = FORC12(NJIND)
                        VC6B(IPQ) = FORC6(NJINDB)
                        VC12B(IPQ) = FORC12(NJINDB)
                        AC12C6(IPQ) = FC12C6(NJIND)
                        BC12C6(IPQ) = FC12C6(NJINDB)
                     ENDDO
                  ELSE
!$OMP parallel do private(IPQ,NJIND) shared(I)
                     DO IPQ = 1, JLEN
                        NJIND = JRCPAC(IACAJ(IPQ)+(IACT(I)-1)*NRATT) + 
     $                       IFORCE(IPQ) + JRCBK2(IPQ)
                        VC6A(IPQ) = FORC6(NJIND)
                        VC12A(IPQ) = FORC12(NJIND)
                        
                     ENDDO
                  ENDIF
C     

                  IF (LVAC) THEN
!$OMP parallel do private(IPQ) shared(I)
                     DO IPQ = 1, JLEN
                        VXIJ(IPQ) = XC(I) - XJ(IPQ)
                        VYIJ(IPQ) = YC(I) - YJ(IPQ)
                        VZIJ(IPQ) = ZC(I) - ZJ(IPQ)
                        VWIJ(IPQ) = WC(I) - WJ(IPQ)
                        R2IJ3D(IPQ) = VXIJ(IPQ)*VXIJ(IPQ) + 
     .                                VYIJ(IPQ)*VYIJ(IPQ) + 
     .                     VZIJ(IPQ)*(VZIJ(IPQ)+COSB2*VXIJ(IPQ))
                        VRIJ2(IPQ) = R2IJ3D(IPQ)+VWIJ(IPQ)*VWIJ(IPQ)
                     ENDDO
                  ELSE
!$OMP parallel do private(IPQ,XIJ,YIJ,ZIJ,WIJ) shared(I)
                     DO IPQ = 1, JLEN
                        XIJ = (XC(I)-XJ(IPQ))*BOXIV1
                        YIJ = (YC(I)-YJ(IPQ))*BOXIV2
                        ZIJ = (ZC(I)-ZJ(IPQ))*BOXIV3
                        WIJ = (WC(I)-WJ(IPQ))*BOXIV4
C     
                        VXIJ(IPQ) = (XIJ-ANINT(XIJ))*BOX1
                        VYIJ(IPQ) = (YIJ-ANINT(YIJ))*BOX2
                        VZIJ(IPQ) = (ZIJ-ANINT(ZIJ))*BOX3
                        VWIJ(IPQ) = (WIJ-ANINT(WIJ))*BOX4
C     
                        R2IJ3D(IPQ) = VXIJ(IPQ)*VXIJ(IPQ)+VYIJ(IPQ)
     .                              * VYIJ(IPQ)+VZIJ(IPQ)*(VZIJ(IPQ)
     .                              + COSB2*VXIJ(IPQ))
                        VRIJ2(IPQ) = R2IJ3D(IPQ)+VWIJ(IPQ)*VWIJ(IPQ)
                     ENDDO
                  ENDIF
                  IF (LOCTO) THEN
!$OMP parallel do private(IPQ)
                     DO IPQ = 1, JLEN
                        VEELJ(IPQ) = BOXOQ-DABS(VXIJ(IPQ))-
     .                       DABS(VYIJ(IPQ))-ABS(VZIJ(IPQ))
                        IF (VEELJ(IPQ) .LE. 0.D0) THEN
                           VRIJ2(IPQ) = VRIJ2(IPQ) +VEELJ(IPQ)*BOX1
                           R2IJ3D(IPQ)= R2IJ3D(IPQ)+VEELJ(IPQ)*BOX1
                           VXIJ(IPQ) = VXIJ(IPQ)-SIGN(BOXOH,VXIJ(IPQ))
                           VYIJ(IPQ) = VYIJ(IPQ)-SIGN(BOXOH,VYIJ(IPQ))
                           VZIJ(IPQ) = VZIJ(IPQ)-SIGN(BOXOH,VZIJ(IPQ))
                        ENDIF
                     ENDDO
                  ENDIF
C     
                  IPQ = ISRCHFLE(JLEN,VRIJ2,1,1.0D-20)
                  IF (IPQ .LE. JLEN) THEN
                     PRINT 8240, NI, I
 8240                FORMAT(/' A DISTANCE THAT IS TOO SMALL HAS',
     $                    ' BEEN FOUND',
     $                    ' FOR CHARGE GROUP ',I4,'  AND ATOM ',I5/)
                     LERROR = .TRUE.
                     GO TO 5750
                  ENDIF
C     
!$OMP parallel do private(IPQ)
                  DO IPQ = 1, JLEN
                     VRIJ(IPQ) = DSQRT(VRIJ2(IPQ))
                     VRINV(IPQ) = AMULT(IPQ)/VRIJ(IPQ)
                     VRINV2(IPQ) = VRINV(IPQ)*VRINV(IPQ)
                     VRINV6(IPQ) = VRINV2(IPQ)*VRINV2(IPQ)*
     .                             VRINV2(IPQ)
                  ENDDO
C     
C     
C                  IF (LPION) THEN
C                     LIZERO = (IPICT(I) .EQ. 0)   
C                     DO IPQ = 1, JLEN
C                        LJZERO = (IPICJ(IPQ) .EQ. 0)
C                        IF (LIZERO .AND. LJZERO) THEN
C                           PIFACT(IPQ) = 1.0
C                        ELSEIF (LIZERO .OR. LJZERO) THEN
C                           PIFACT(IPQ) = PININV
C                        ELSEIF (IPICT(I) .EQ. IPICJ(IPQ)) THEN
C                           PIFACT(IPQ) = PININV
C                        ELSE
C                           PIFACT(IPQ) = 0.0
C                        ENDIF
C                     ENDDO
C                  ENDIF 
C     
C     
                  IF (LPERTL) THEN
                     LIPERT = LLPERT(I)
                     IF (LIPERT) THEN
                        DO IPQ = 1, JLEN
C     
                           CA6 = VC6A(IPQ)
                           CA12 = VC12A(IPQ)
                           CB6 = VC6B(IPQ)
                           CB12 = VC12B(IPQ)

                           QIQJA = CGAT(I) * CGAJ(IPQ) * AMULT(IPQ)
                           QIQJB = CGBT(I) * CGBJ(IPQ) * AMULT(IPQ)

                           CA126 = AC12C6(IPQ)
                           CB126 = BC12C6(IPQ)

                           SOFTLJ = ALPHLJ * 
     $                          (1 - ISCLJT(I) * ISCLJJ(IPQ))
                           SOFTEL = ALPHC * 
     $                          (1 - ISCCT(I) * ISCCJ(IPQ))

C*******************************************************************************
C     4D
C*******************************************************************************
C
C---- SOFT-CORE PERTURBATION
C     CALCULATE FULL 4D CONTRIBUTION
C     
                           RIJ2 = VRIJ2(IPQ)
                           RIJ3D2 = R2IJ3D(IPQ)
                           RIJ4 = RIJ2 * RIJ2
                           RIJ6 = RIJ4 * RIJ2
C     
C--   4D, STATE A
C     
C     COULOMB
                           SFEL   = SOFTEL * RL2A
                           RSOFT2 = 1.0D0 / (SFEL + RIJ2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELA4D = QIQJA * RSOFT                    
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFA4 = QIQJA * C1 * RRF
                           QRFR   = DFRFA4 * RIJ3D2
                           VRFA4D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CA126
                           CRA1   = 1.0D0 / (SFLJ * RL2A + RIJ6)
                           CRA12  = CA12 * CRA1
                           VLJA4D = CRA1 * (CRA12 - CA6)
                           DFA4D  = CRA1 * ( CRA1 * CRA12 + VLJA4D)
C     DV/DL
                           DLLJA4 = -2.0D0 * SFLJ * DFA4D
                           DLELA4 = -VELA4D * SOFTEL * RSOFT2
                           DLRFA4 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFA4D = RSOFT2*VELA4D+6.0D0*RIJ4*DFA4D
C     
C--   4D, STATE B
C     
C     COULOMB
                           SFEL   = SOFTEL * RL2B 
                           RSOFT2 = 1.0D0 / (SFEL + RIJ2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELB4D = QIQJB * RSOFT
C     REACTION FIELD
                           RRF2   =  1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFB4 = QIQJB * C1 * RRF
                           QRFR   = DFRFB4 * RIJ3D2
                           VRFB4D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CB126
                           CRA1   = 1.0D0 / (SFLJ * RL2B + RIJ6)
                           CRA12  = CB12 * CRA1
                           VLJB4D = CRA1 * (CRA12 - CB6)
                           DFB4D  = CRA1 * ( CRA1 * CRA12 + VLJB4D)
C     DV/DL
                           DLLJB4 = -2.0D0 * SFLJ * DFB4D
                           DLELB4 = -VELB4D * SOFTEL * RSOFT2
                           DLRFB4 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFB4D = RSOFT2*VELB4D+6.0D0*RIJ4*DFB4D
C     
C--   COMBINE 4D STATE A AND B
C     
                           VEL4D  = RLB  * VELB4D + RLA  * VELA4D
                           VRF4D  = RLB  * VRFB4D + RLA  * VRFA4D
                           VLJ4D  = RLB  * VLJB4D + RLA  * VLJA4D
                           DF4D4D = RLB  * DFB4D  + RLA  * DFA4D
                           DF4D   = DF4D4D 
     .                            + RLB  * DFRFB4 + RLA  * DFRFA4
                           DLLJ4D = RDLB * VLJB4D + RLLB * DLLJB4
     .                            + RDLA * VLJA4D + RLLA * DLLJA4
                           DLEL4D = RDLB * VELB4D + RLLB * DLELB4
     .                            + RDLA * VELA4D + RLLA * DLELA4
                           DLRF4D = RDLB * VRFB4D + RLLB * DLRFB4
     .                            + RDLA * VRFA4D + RLLA * DLRFA4

C*******************************************************************************
C     3D
C*******************************************************************************
                           RIJ4 = RIJ3D2 * RIJ3D2
                           RIJ6 = RIJ4 * RIJ3D2
C     
C--   3D, STATE A
C     
C     COULOMB
                           SFEL   = SOFTEL * RLMA
                           RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELA3D = QIQJA * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFA3 = QIQJA * C1 * RRF
                           QRFR   = DFRFA3 * RIJ3D2
                           VRFA3D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CA126
                           CRA1   = 1.0D0 / (SFLJ * RLMA + RIJ6)
                           CRA12  = CA12 * CRA1
                           VLJA3D = CRA1 * (CRA12 - CA6)
                           DFA3D  = CRA1 * ( CRA1 * CRA12 + VLJA3D)
C     DV/DL
                           DLLJA3 = -2.0D0 * SFLJ * DFA3D
                           DLELA3 = -VELA3D * SOFTEL * RSOFT2 
                           DLRFA3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFA3D = RSOFT2*VELA3D+6.0D0*RIJ4*DFA3D
C     
C--   3D, STATE B
C     
C     COULOMB
                           SFEL   = SOFTEL * RLMB
                           RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELB3D = QIQJB * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFB3 = QIQJB * C1 * RRF
                           QRFR   = DFRFB3 * RIJ3D2
                           VRFB3D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CB126 
                           CRA1   = 1.0D0 / (SFLJ * RLMB + RIJ6)
                           CRA12  = CB12 * CRA1
                           VLJB3D = CRA1 * (CRA12 - CB6)
                           DFB3D  = CRA1 * ( CRA1 * CRA12 + VLJB3D)
C     DV/DL
                           DLLJB3 = -2.0D0 * SFLJ * DFB3D
                           DLELB3 = -VELB3D * SOFTEL * RSOFT2 
                           DLRFB3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFB3D = RSOFT2*VELB3D+6.0D0*RIJ4*DFB3D

C     
C--   COMBINE 3D, STATE A AND B
C     
                           VEL3D =  RLB  * VELB3D + RLA  * VELA3D
                           VRF3D =  RLB  * VRFB3D + RLA  * VRFA3D
                           VLJ3D =  RLB  * VLJB3D + RLA  * VLJA3D
                           DF3D  =  RLB  * DFB3D  + RLA  * DFA3D
     $                            + RLB  * DFRFB3 + RLA  * DFRFA3
                           DLLJ3D = RDLB * VLJB3D + RLLB * DLLJB3
     $                            + RDLA * VLJA3D + RLLA * DLLJA3
                           DLEL3D = RDLB * VELB3D + RLLB * DLELB3 
     $                            + RDLA * VELA3D + RLLA * DLELA3
                           DLRF3D = RDLB * VRFB3D + RLLB * DLRFB3
     $                            + RDLA * VRFA3D + RLLA * DLRFA3

                           DDM3EL =  RLBRMU * DLELB3 + RLARMU * DLELA3
                           DDM3RF =  RLBRMU * DLRFB3 + RLARMU * DLRFA3
                           DDM3LJ =  RLBRMU * DLLJB3 + RLARMU * DLLJA3

C     REACTION FIELD CONSTANT
                           VRCA = RFC*QIQJA
                           VRCB = RFC*QIQJB
                           VRC  = RLB*VRCB  + RLA*VRCA
                           DLRC = RDLB*VRCB + RDLA*VRCA
C     
C---- COMBINE 4D AND 3D COMPONENTS
                           DF3T = RM3D * DF3D + RM4D * DF4D
                           DF4T = RM4D * DF4D4D
                           DF3(IPQ) = DF3T * AMULT(IPQ)
                           DF4TH(IPQ) = DF4T * AMULT(IPQ)

                           VEELJ(IPQ) = 
     $                          (RM3D * VEL3D + RM4D * VEL4D)
     $                          * AMULT(IPQ)
                           VERFJ(IPQ) =
     $                          (RM3D * VRF3D + RM4D * VRF4D)
     $                          * AMULT(IPQ)
                           VERCJ(IPQ) =
     $                          (RM3D * VRC   + RM4D * VRC) 
     $                          * AMULT(IPQ)
                           VLJJ(IPQ)  = 
     $                          (RM3D * VLJ3D + RM4D * VLJ4D) 
     $                          * AMULT(IPQ)
C     
                           EGLELJ(IPQ) =
     $                          (RM3D * DLEL3D + RM4D * DLEL4D)
     $                          * AMULT(IPQ)
                           EGLRFJ(IPQ) =
     $                          (RM3D * DLRF3D + RM4D * DLRF4D)
     $                          * AMULT(IPQ)
                           EGLRCJ(IPQ) =
     $                          (RM3D * DLRC + RM4D * DLRC)
     $                          * AMULT(IPQ)
                           EGLLJJ(IPQ) =
     $                          (RM3D * DLLJ3D + RM4D * DLLJ4D)
     $                          * AMULT(IPQ)
C     
                           E34ELJ(IPQ) = 
     $                          (RDM3D*VEL3D + RM3D*DDM3EL + 
     $                          RDM4D*VEL4D) * AMULT(IPQ)
                           E34RFJ(IPQ) =
     $                          (RDM3D*VRF3D + RM3D*DDM3RF +
     $                          RDM4D*VRF4D) * AMULT(IPQ)
                           E34RCJ(IPQ) =
     $                          (RDM3D*VRC + RDM4D*VRC) * AMULT(IPQ) 
                           E34LJJ(IPQ) =
     $                          (RDM3D*VLJ3D + RM3D*DDM3LJ + 
     $                          RDM4D*VLJ4D) * AMULT(IPQ) 
                           
                        ENDDO

                     ELSE
                        DO IPQ = 1, JLEN
                           IF (LLPERJ(IPQ)) THEN

                           CA6 = VC6A(IPQ)
                           CA12 = VC12A(IPQ)
                           CB6 = VC6B(IPQ)
                           CB12 = VC12B(IPQ)

                           QIQJA = CGAT(I) * CGAJ(IPQ) * AMULT(IPQ)
                           QIQJB = CGBT(I) * CGBJ(IPQ) * AMULT(IPQ)

                           CA126 = AC12C6(IPQ)
                           CB126 = BC12C6(IPQ)

                           SOFTLJ = ALPHLJ * 
     $                          (1 - ISCLJT(I) * ISCLJJ(IPQ))
                           SOFTEL = ALPHC * 
     $                          (1 - ISCCT(I) * ISCCJ(IPQ))

C*******************************************************************************
C     4D
C*******************************************************************************
C     
C---- SOFT-CORE PERTURBATION
C     CALCULATE FULL 4D CONTRIBUTION
C     
                           RIJ2 = VRIJ2(IPQ)
                           RIJ3D2 = R2IJ3D(IPQ)
                           RIJ4 = RIJ2 * RIJ2
                           RIJ6 = RIJ4 * RIJ2
C     
C--   4D, STATE A
C     
C     COULOMB
                           SFEL   = SOFTEL * RL2A
                           RSOFT2 = 1.0D0 / (SFEL + RIJ2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELA4D = QIQJA * RSOFT                    
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFA4 = QIQJA * C1 * RRF
                           QRFR   = DFRFA4 * RIJ3D2
                           VRFA4D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CA126
                           CRA1   = 1.0D0 / (SFLJ * RL2A + RIJ6)
                           CRA12  = CA12 * CRA1
                           VLJA4D = CRA1 * (CRA12 - CA6)
                           DFA4D  = CRA1 * ( CRA1 * CRA12 + VLJA4D)
C     DV/DL
                           DLLJA4 = -2.0D0 * SFLJ * DFA4D
                           DLELA4 = -VELA4D * SOFTEL * RSOFT2
                           DLRFA4 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFA4D = RSOFT2*VELA4D+6.0D0*RIJ4*DFA4D
C     
C--   4D, STATE B
C     
C     COULOMB
                           SFEL   = SOFTEL * RL2B 
                           RSOFT2 = 1.0D0 / (SFEL + RIJ2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELB4D = QIQJB * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFB4 = QIQJB * C1 * RRF
                           QRFR   = DFRFB4 * RIJ3D2
                           VRFB4D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CB126
                           CRA1   = 1.0D0 / (SFLJ * RL2B + RIJ6)
                           CRA12  = CB12 * CRA1
                           VLJB4D = CRA1 * (CRA12 - CB6)
                           DFB4D  = CRA1 * ( CRA1 * CRA12 + VLJB4D)
C     DV/DL
                           DLLJB4 = -2.0D0 * SFLJ * DFB4D
                           DLELB4 = -VELB4D * SOFTEL * RSOFT2
                           DLRFB4 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFB4D = RSOFT2*VELB4D+6.0D0*RIJ4*DFB4D
C     
C--   COMBINE 4D STATE A AND B
C     
                           VEL4D  =    RLB  * VELB4D + RLA  * VELA4D
                           VRF4D  =    RLB  * VRFB4D + RLA  * VRFA4D
                           VLJ4D  =    RLB  * VLJB4D + RLA  * VLJA4D
                           DF4D4D =    RLB  * DFB4D  + RLA  * DFA4D
                           DF4D   =    DF4D4D 
     $                               + RLB  * DFRFB4 + RLA  * DFRFA4
                           DLLJ4D =    RDLB * VLJB4D + RLLB * DLLJB4
     $                               + RDLA * VLJA4D + RLLA * DLLJA4
                           DLEL4D =    RDLB * VELB4D + RLLB * DLELB4
     $                               + RDLA * VELA4D + RLLA * DLELA4
                           DLRF4D =    RDLB * VRFB4D + RLLB * DLRFB4
     $                               + RDLA * VRFA4D + RLLA * DLRFA4

C*******************************************************************************
C     3D
C*******************************************************************************
C     
C     CALCULATE FULL 3D CONTRIBUTION
C     
                           RIJ4 = RIJ3D2 * RIJ3D2
                           RIJ6 = RIJ4 * RIJ3D2
C     
C--   3D, STATE A
C     
C     COULOMB
                           SFEL   = SOFTEL * RLMA
                           RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELA3D = QIQJA * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFA3 = QIQJA * C1 * RRF
                           QRFR   = DFRFA3 * RIJ3D2
                           VRFA3D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CA126
                           CRA1   = 1.0D0 / (SFLJ * RLMA + RIJ6)
                           CRA12  = CA12 * CRA1
                           VLJA3D = CRA1 * (CRA12 - CA6)
                           DFA3D  = CRA1 * ( CRA1 * CRA12 + VLJA3D)
C     DV/DL
                           DLLJA3 = -2.0D0 * SFLJ * DFA3D
                           DLELA3 = -VELA3D * SOFTEL * RSOFT2 
                           DLRFA3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFA3D = RSOFT2*VELA3D+6.0D0*RIJ4*DFA3D
C     
C--   3D, STATE B
C     
C     COULOMB
                           SFEL   = SOFTEL * RLMB
                           RSOFT2 = 1.0D0 / (SFEL + RIJ3D2)
                           RSOFT  = DSQRT(RSOFT2)
                           VELB3D = QIQJB * RSOFT
C     REACTION FIELD
                           RRF2   = 1.0D0 / (SFEL + RCRF2)
                           RFSQRT = DSQRT(RRF2)
                           RRF    = RRF2 * RFSQRT
                           DFRFB3 = QIQJB * C1 * RRF
                           QRFR   = DFRFB3 * RIJ3D2
                           VRFB3D = -0.5D0 * QRFR
C     LENNARD-JONES
                           SFLJ   = SOFTLJ * CB126 
                           CRA1   = 1.0D0 / (SFLJ * RLMB + RIJ6)
                           CRA12  = CB12 * CRA1
                           VLJB3D = CRA1 * (CRA12 - CB6)
                           DFB3D  = CRA1 * ( CRA1 * CRA12 + VLJB3D)
C     DV/DL
                           DLLJB3 = -2.0D0 * SFLJ * DFB3D
                           DLELB3 = -VELB3D * SOFTEL * RSOFT2 
                           DLRFB3 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                           DFB3D = RSOFT2*VELB3D+6.0D0*RIJ4*DFB3D

C     
C--   COMBINE 3D, STATE A AND B
C     
                           VEL3D =  RLB  * VELB3D + RLA  * VELA3D
                           VRF3D =  RLB  * VRFB3D + RLA  * VRFA3D
                           VLJ3D =  RLB  * VLJB3D + RLA  * VLJA3D
                           DF3D  =  RLB  * DFB3D  + RLA  * DFA3D
     $                            + RLB  * DFRFB3 + RLA  * DFRFA3
                           DLLJ3D = RDLB * VLJB3D + RLLB * DLLJB3
     $                            + RDLA * VLJA3D + RLLA * DLLJA3
                           DLEL3D = RDLB * VELB3D + RLLB * DLELB3 
     $                            + RDLA * VELA3D + RLLA * DLELA3
                           DLRF3D = RDLB * VRFB3D + RLLB * DLRFB3
     $                            + RDLA * VRFA3D + RLLA * DLRFA3

                           DDM3EL =  RLBRMU * DLELB3 + RLARMU * DLELA3
                           DDM3RF =  RLBRMU * DLRFB3 + RLARMU * DLRFA3
                           DDM3LJ =  RLBRMU * DLLJB3 + RLARMU * DLLJA3

C     REACTION FILED CONSTANT

                           VRCA = RFC*QIQJA
                           VRCB = RFC*QIQJB
                           VRC  = RLB*VRCB  + RLA*VRCA
                           DLRC = RDLB*VRCB + RDLA*VRCA
C     
C---- COMBINE 4D AND 3D COMPONENTS
                           DF3T = RM3D * DF3D + RM4D * DF4D
                           DF4T = RM4D * DF4D4D
                           DF3(IPQ) = DF3T * AMULT(IPQ)
                           DF4TH(IPQ) = DF4T * AMULT(IPQ)

                           VEELJ(IPQ) = 
     $                          (RM3D * VEL3D + RM4D * VEL4D)
     $                          * AMULT(IPQ)
                           VERFJ(IPQ) =
     $                          (RM3D * VRF3D + RM4D * VRF4D)
     $                          * AMULT(IPQ)
                           VERCJ(IPQ) =
     $                          (RM3D * VRC   + RM4D * VRC) 
     $                          * AMULT(IPQ)
                           VLJJ(IPQ)  = 
     $                          (RM3D * VLJ3D + RM4D * VLJ4D) 
     $                          * AMULT(IPQ)
C     
                           EGLELJ(IPQ) =
     $                          (RM3D * DLEL3D + RM4D * DLEL4D)
     $                          * AMULT(IPQ)
                           EGLRFJ(IPQ) =
     $                          (RM3D * DLRF3D + RM4D * DLRF4D)
     $                          * AMULT(IPQ)
                           EGLRCJ(IPQ) =
     $                          (RM3D * DLRC   + RM4D * DLRC)
     $                          * AMULT(IPQ)
                           EGLLJJ(IPQ) =
     $                          (RM3D * DLLJ3D + RM4D * DLLJ4D) *
     $                          AMULT(IPQ) 
C     
                           E34ELJ(IPQ) = 
     $                          (RDM3D*VEL3D + RM3D*DDM3EL + 
     $                          RDM4D*VEL4D) * AMULT(IPQ)
                           E34RFJ(IPQ) =
     $                          (RDM3D*VRF3D + RM3D*DDM3RF +
     $                          RDM4D*VRF4D) * AMULT(IPQ)
                           E34RCJ(IPQ) =
     $                          (RDM3D*VRC + RDM4D*VRC) * AMULT(IPQ) 
                           E34LJJ(IPQ) =
     $                          (RDM3D*VLJ3D + RM3D*DDM3LJ + 
     $                          RDM4D*VLJ4D) * AMULT(IPQ) 
                           
                        ELSE
                           QIQJA = CGAT(I) * CGAJ(IPQ) * AMULT(IPQ)
                           VERFJ(IPQ) = QIQJA * RFE * VRIJ2(IPQ)
                           VERCJ(IPQ) = QIQJA * RFC
                           VEELJ(IPQ) = QIQJA * VRINV(IPQ)
                           VLJJ(IPQ) = VRINV6(IPQ)*(VC12A(IPQ)
     $                          *VRINV6(IPQ)-VC6A(IPQ))
                           
                           EGLELJ(IPQ) = 0.0D0
                           EGLRFJ(IPQ) = 0.0D0
                           EGLRCJ(IPQ) = 0.0D0
                           EGLLJJ(IPQ) = 0.0D0

                           E34ELJ(IPQ) = 0.0D0
                           E34RFJ(IPQ) = 0.0D0
                           E34RCJ(IPQ) = 0.0D0
                           E34LJJ(IPQ) = 0.0D0

                           DF = VRINV6(IPQ)*(2.0D0*VRINV6(IPQ) *
     $                          VC12A(IPQ)-VC6A(IPQ))
                           DF3T = (6.0D0 * DF + 
     $                          VEELJ(IPQ))*VRINV2(IPQ) * AMULT(IPQ)

                           DF3(IPQ) = DF3T +  QIQJA * RFF
                           DF4TH(IPQ) = DF3T
                        ENDIF
                        
                     ENDDO

                  ENDIF
!$OMP parallel do private(IPQ) 
!$OMP+ reduction(+:EGLEL,EGLRF,EGLRC,EGLLJ,EG34EL,EG34RF,EG34RC,EG34LJ)
                  DO  IPQ = 1, JLEN
                     EGLEL = EGLEL + EGLELJ(IPQ)
                     EGLRF = EGLRF + EGLRFJ(IPQ)
                     EGLRC = EGLRC + EGLRCJ(IPQ)
                     EGLLJ = EGLLJ + EGLLJJ(IPQ)
                     EG34EL =  EG34EL + E34ELJ(IPQ) 
                     EG34RF =  EG34RF + E34RFJ(IPQ)  
                     EG34RC =  EG34RC + E34RCJ(IPQ)  
                     EG34LJ =  EG34LJ + E34LJJ(IPQ)  
                  ENDDO
                  
               ELSE
                     
C                     IF (LPION) THEN
C                        DO IPQ = 1, JLEN
C                           QIQJA = CGAT(I) * CGAJ(IPQ) * AMULT(IPQ) 
C     $                          * PIFACT(IPQ)
C                           VERFJ(IPQ) = QIQJA * RFE * VRIJ2(IPQ)
C                           VERCJ(IPQ) = QIQJA * RFC
C                           VEELJ(IPQ) = QIQJA * VRINV(IPQ)
C                           VLJJ(IPQ) = VRINV6(IPQ)*(VC12A(IPQ)*VRINV6(
C     $                          IPQ)-VC6A(IPQ)) * PIFACT(IPQ)
C                           DF = VRINV6(IPQ)*(2.0E0*VRINV6(IPQ)*
C     $                          VC12A(IPQ)-VC6A(IPQ)) * PIFACT(IPQ)
C                           DF3(IPQ) = (6.0E0 * DF + VEELJ(IPQ))*
C     $                          VRINV2(IPQ) +  QIQJA * RFF
C                        ENDDO
C                     ELSE
!$OMP parallel do private(IPQ,QIQJA,DF,DF3T) shared(RFC,RFE,RFF)
                        DO IPQ = 1, JLEN
                           QIQJA = CGAT(I) * CGAJ(IPQ) * AMULT(IPQ)
                           VERFJ(IPQ) = QIQJA * RFE * R2IJ3D(IPQ)
                           VERCJ(IPQ) = QIQJA * RFC
                           VEELJ(IPQ) = QIQJA * VRINV(IPQ)
                           VLJJ(IPQ) = VRINV6(IPQ)*(VC12A(IPQ)*
     .                                 VRINV6(IPQ)-VC6A(IPQ))
                           DF = VRINV6(IPQ)*(2.0D0*VRINV6(IPQ)*
     .                          VC12A(IPQ)-VC6A(IPQ))
                           DF3T = (6.0D0 * DF + VEELJ(IPQ))*
     .                          VRINV2(IPQ) 
                           DF3(IPQ) = DF3T + QIQJA * RFF
                           DF4TH(IPQ) = DF3T
                        ENDDO
C                     ENDIF
                  ENDIF
C     
C     
                  JRCTOT = JLEN
                  JRCSRT = 0
 8950             CONTINUE
                  JRCLEN = MIN(IVECT,JRCTOT)
CDIR$ SHORTLOOP
                  DO JRC = 1, JRCLEN
                     TEEL(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) = 
     $                    TEEL(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) + 
     $                    VEELJ(JRC+JRCSRT)
                     TELJ(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) = 
     $                    TELJ(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) + 
     $                    VLJJ(JRC+JRCSRT)
                     
                     TERF(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TERF(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VERFJ(JRC+JRCSRT)
                     
                     TERC(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TERC(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VERCJ(JRC+JRCSRT)
                     
                  ENDDO

                  JRCSRT = JRCSRT + JRCLEN
                  JRCTOT = JRCTOT - JRCLEN
                  IF (JRCTOT .GT. 0) GO TO 8950
C     
C*****CALCULATE THE FORCES.
C     
!$OMP parallel do private(IPQ)
                  DO IPQ = 1, JLEN
                     VXH(IPQ) = VXIJ(IPQ)*DF3(IPQ)
                     VYH(IPQ) = VYIJ(IPQ)*DF3(IPQ)
                     VZH(IPQ) = VZIJ(IPQ)*DF3(IPQ)
                     VWH(IPQ) = VWIJ(IPQ)*DF4TH(IPQ)
C     
                     FXTEMP(IPQ) = FXTEMP(IPQ) - VXH(IPQ)
                     FYTEMP(IPQ) = FYTEMP(IPQ) - VYH(IPQ)
                     FZTEMP(IPQ) = FZTEMP(IPQ) - VZH(IPQ)
                     FWTEMP(IPQ) = FWTEMP(IPQ) - VWH(IPQ)
                  END DO
C
C                  IF (LDOVIR) THEN
C                     IF (LVIR) THEN
C                        DO IPQ = 1, IVIRL
C                           TVIR(IPQ,1) = TVIR(IPQ,1) - VXH(IOCVIR(
C     $                          IPQ))*(VXIJ(IOCVIR(IPQ))+XRJ(IOCVIR(
C     $                          IPQ))-XRC(I))
C                           TVIR(IPQ,2) = TVIR(IPQ,2) - VYH(IOCVIR(
C     $                          IPQ))*(VYIJ(IOCVIR(IPQ))+YRJ(IOCVIR(
C     $                          IPQ))-YRC(I))
C                           TVIR(IPQ,3) = TVIR(IPQ,3) - VZH(IOCVIR(
C     $                          IPQ))*(VZIJ(IOCVIR(IPQ))+ZRJ(IOCVIR(
C     $                          IPQ))-ZRC(I))
C                        ENDDO                           
C                     ENDIF
C                  ENDIF
C     
                  DO IPQ = 1, JLEN
                     FLRX(I) = FLRX(I) + VXH(IPQ)
                     FLRY(I) = FLRY(I) + VYH(IPQ)
                     FLRZ(I) = FLRZ(I) + VZH(IPQ)
                     FLRW(I) = FLRW(I) + VWH(IPQ)
                  ENDDO
C     
                  JN = JN + JLEN
               ENDDO
C     
C     
!$OMP parallel do private(I20X)
               DO I20X = 1, JLEN
                  FLRX(JLIST(I20X)) = FXTEMP(I20X)
                  FLRY(JLIST(I20X)) = FYTEMP(I20X)
                  FLRZ(JLIST(I20X)) = FZTEMP(I20X)
                  FLRW(JLIST(I20X)) = FWTEMP(I20X)
               ENDDO
C     
            ENDIF
         ENDDO
C     
C*****CALCULATE THE INTERACTION BETWEEN THE SOLVENT MOLECULES.
C     
         NIMAX = NRAGT - 1
         DO NI = NRAG + 1, NIMAX
            NPAIR = INB(NI)
            IF (NPAIR .GT. 0) THEN
               DO I44X = 1, NPAIR
                  JLIST(I44X) = IFIRST(JNB(I44X+NJNB-1))
               ENDDO
               JLOW = 1
               JLEN = NPAIR + 1
               DO K = 2, NRAM
                  DO IPQ = 1, NPAIR
                     JLIST(JLEN+IPQ-1) = JLIST(JLOW+IPQ-1) + 1
                  ENDDO
                  
                  JLOW = JLEN
                  JLEN = JLEN + NPAIR
               ENDDO
C     
               JLEN = JLEN - 1
               NJNB = NJNB + NPAIR
               
               ILAST = IFIRST(NI) + NRAM - 1

C     
C*****GATHER THE COORDINATES AND THE FORCES.
C     
C               IF (LDOVIR) THEN
C                  DO IPQ = 1, JLEN
C                     XJ(IPQ) = XC(JLIST(IPQ))
C                     YJ(IPQ) = YC(JLIST(IPQ))
C                     ZJ(IPQ) = ZC(JLIST(IPQ))
C                     WJ(IPQ) = WC(JLIST(IPQ))
C     
C                     FXTEMP(IPQ) = FLRX(JLIST(IPQ))
C                     FYTEMP(IPQ) = FLRY(JLIST(IPQ))
C                     FZTEMP(IPQ) = FLRZ(JLIST(IPQ))
C                     FWTEMP(IPQ) = FLRW(JLIST(IPQ)) 
C                     
C                     XRJ(IPQ) = XRC(JLIST(IPQ))
C                     YRJ(IPQ) = YRC(JLIST(IPQ))
C                     ZRJ(IPQ) = ZRC(JLIST(IPQ))
C                     JRCILJ(IPQ) = ICET(JLIST(IPQ))
C                  ENDDO
C               ELSE
!$OMP parallel do private(IPQ)
                  DO IPQ = 1, JLEN
                     XJ(IPQ) = XC(JLIST(IPQ))
                     YJ(IPQ) = YC(JLIST(IPQ))
                     ZJ(IPQ) = ZC(JLIST(IPQ))
                     WJ(IPQ) = WC(JLIST(IPQ))
C     
                     FXTEMP(IPQ) = FLRX(JLIST(IPQ))
                     FYTEMP(IPQ) = FLRY(JLIST(IPQ))
                     FZTEMP(IPQ) = FLRZ(JLIST(IPQ))
                     FWTEMP(IPQ) = FLRW(JLIST(IPQ))
                     JRCILJ(IPQ) = ICET(JLIST(IPQ))
                  END DO
C               ENDIF
C     
C*****DETERMINE THE INTEGER ATOM CODE FOR EACH ATOM J.
C*****GATHER THE CHARGES AND DETERMINE TO WHICH GROUP EACH ATOM BELO
C     
               JTOP = 1
               DO K = 1, NRAM
                  DO IPQ = 1, NPAIR
                     IACAJ(JTOP+IPQ-1) = IACS(K)
                     CGAJ(JTOP+IPQ-1) = CGS(K)
                  ENDDO
                  JTOP = JTOP + NPAIR
               ENDDO
C     
               
               DO I = IFIRST(NI), ILAST

                  ICE = ICET(I)

C     
C*****DETERMINE LJ-PARAMETERS.
C     
                  JTOP = 1
                  DO K = 1, NRAM
                     IPC = MPAC(IACS(K),IACT(I))
                     DO IPQ = 1, NPAIR
                        VC6A(JTOP+IPQ-1) = C6(IPC)
                        VC12A(JTOP+IPQ-1) = C12(IPC)
                     ENDDO
                     JTOP = JTOP + NPAIR
                  ENDDO
C     
!$OMP parallel do private(IPQ,XIJ,YIJ,ZIJ,WIJ) shared(I)
                  DO IPQ = 1, JLEN
                     XIJ = (XC(I)-XJ(IPQ))*BOXIV1
                     YIJ = (YC(I)-YJ(IPQ))*BOXIV2
                     ZIJ = (ZC(I)-ZJ(IPQ))*BOXIV3
                     WIJ = (WC(I)-WJ(IPQ))*BOXIV4 
C     
                     VXIJ(IPQ) = (XIJ-ANINT(XIJ))*BOX1
                     VYIJ(IPQ) = (YIJ-ANINT(YIJ))*BOX2
                     VZIJ(IPQ) = (ZIJ-ANINT(ZIJ))*BOX3
                     VWIJ(IPQ) = (WIJ-ANINT(WIJ))*BOX4
C     
                     R2IJ3D(IPQ) = VXIJ(IPQ)*VXIJ(IPQ) +
     $                             VYIJ(IPQ)*VYIJ(IPQ) +
     $                      VZIJ(IPQ)*(VZIJ(IPQ)+COSB2*VXIJ(IPQ))
                     VRIJ2(IPQ) = R2IJ3D(IPQ) +
     $                    VWIJ(IPQ) *  VWIJ(IPQ)
                  ENDDO
                  IF (LOCTO) THEN
!$OMP parallel do private(IPQ) shared(BOXOQ,BOX1)
                     DO IPQ = 1, JLEN
                        VEELJ(IPQ) = BOXOQ - DABS(VXIJ(IPQ)) - 
     .                       DABS(VYIJ(IPQ)) - DABS(VZIJ(IPQ))
                        IF (VEELJ(IPQ) .LE. 0.D0) THEN
                           VRIJ2(IPQ) = VRIJ2(IPQ)  + VEELJ(IPQ)*BOX1
                           R2IJ3D(IPQ)= R2IJ3D(IPQ) + VEELJ(IPQ)*BOX1
                           VXIJ(IPQ) = VXIJ(IPQ)-SIGN(BOXOH,VXIJ(IPQ))
                           VYIJ(IPQ) = VYIJ(IPQ)-SIGN(BOXOH,VYIJ(IPQ))
                           VZIJ(IPQ) = VZIJ(IPQ)-SIGN(BOXOH,VZIJ(IPQ))
                        ENDIF
                     ENDDO
                  ENDIF
C     
C--   CHECK FOR CLOSE CONTACTS
                  IPQ = ISRCHFLE(JLEN,VRIJ2,1,1.0D-20)
                  IF (IPQ .LE. JLEN) THEN
                     PRINT 9240, NI, I
                     LERROR = .TRUE.
                     GO TO 5750
                  ENDIF
C
!$OMP parallel do private(IPQ)
                  DO IPQ = 1, JLEN
                     VRIJ(IPQ) = DSQRT(VRIJ2(IPQ))
                     VRINV(IPQ) = 1.0D0/VRIJ(IPQ)
                     VRINV2(IPQ) = VRINV(IPQ)*VRINV(IPQ)
                     VRINV6(IPQ) = VRINV2(IPQ)*VRINV2(IPQ)*
     .                             VRINV2(IPQ)
                  ENDDO
C     
!$OMP parallel do private(IPQ,QIQJA,DF,DF3T) shared(I,RFC,RFE,RFF)
                  DO IPQ = 1, JLEN
                     QIQJA = CGAT(I) * CGAJ(IPQ)
                     VERFJ(IPQ) = QIQJA * RFE * VRIJ2(IPQ)
                     VERCJ(IPQ) = QIQJA * RFC
                     VEELJ(IPQ) = QIQJA*VRINV(IPQ)
                     VLJJ(IPQ) = VRINV6(IPQ)*(VC12A(IPQ)*
     .                           VRINV6(IPQ)-VC6A(IPQ))
                     DF = VRINV6(IPQ)*(2.0D0*VRINV6(IPQ)*VC12A(IPQ)
     .                    -VC6A(IPQ))
                     DF3T = (6.0D0 * DF + VEELJ(IPQ))*VRINV2(IPQ) 
                     DF3(IPQ) = DF3T +  QIQJA * RFF
                     DF4TH(IPQ) = DF3T
                  ENDDO
C     
                  
                  JRCTOT = JLEN
                  JRCSRT = 0
 8150             CONTINUE
                  JRCLEN = MIN(IVECT,JRCTOT)
                  
                  DO JRC = 1, JRCLEN

                     TEEL(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) = 
     $                    TEEL(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) + 
     $                    VEELJ(JRC+JRCSRT)
                     TELJ(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) = 
     $                    TELJ(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) + 
     $                    VLJJ(JRC+JRCSRT)
                     
                     TERF(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TERF(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VERFJ(JRC+JRCSRT)
                     
                     TERC(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TERC(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VERCJ(JRC+JRCSRT)
                     
                  ENDDO
                  JRCSRT = JRCSRT + JRCLEN
                  JRCTOT = JRCTOT - JRCLEN
                  IF (JRCTOT .GT. 0) GO TO 8150
                  
C     
C*****CALCULATE THE FORCES.
C     
!$OMP parallel do private(IPQ)
                  DO IPQ = 1, JLEN
                     VXH(IPQ) = VXIJ(IPQ)*DF3(IPQ)
                     VYH(IPQ) = VYIJ(IPQ)*DF3(IPQ)
                     VZH(IPQ) = VZIJ(IPQ)*DF3(IPQ)
                     VWH(IPQ) = VWIJ(IPQ)*DF4TH(IPQ)
C     
                     FXTEMP(IPQ) = FXTEMP(IPQ) - VXH(IPQ)
                     FYTEMP(IPQ) = FYTEMP(IPQ) - VYH(IPQ)
                     FZTEMP(IPQ) = FZTEMP(IPQ) - VZH(IPQ)
                     FWTEMP(IPQ) = FWTEMP(IPQ) - VWH(IPQ)
                  ENDDO
                  
C     
C                  IF (LDOVIR) THEN
C                     DO IPQ = 1, JLEN
C                        TVIR(IPQ,1) = TVIR(IPQ,1) - VXH(IPQ)*
C     $                       (VXIJ(IPQ)+XRJ(IPQ)-XRC(I))
C                        TVIR(IPQ,2) = TVIR(IPQ,2) - VYH(IPQ)*
C     $                       (VYIJ(IPQ)+YRJ(IPQ)-YRC(I))
C                        TVIR(IPQ,3) = TVIR(IPQ,3) - VZH(IPQ)*
C     $                       (VZIJ(IPQ)+ZRJ(IPQ)-ZRC(I))
C                     ENDDO
C                  ENDIF
C     
                  DO IPQ = 1, JLEN
                     FLRX(I) = FLRX(I) + VXH(IPQ)
                     FLRY(I) = FLRY(I) + VYH(IPQ)
                     FLRZ(I) = FLRZ(I) + VZH(IPQ)
                     FLRW(I) = FLRW(I) + VWH(IPQ)
                  ENDDO
C     
                  JN = JN + JLEN
               ENDDO
C     
C     
!$OMP parallel do private(I32X)
               DO I32X = 1, JLEN
                  FLRX(JLIST(I32X)) = FXTEMP(I32X)
                  FLRY(JLIST(I32X)) = FYTEMP(I32X)
                  FLRZ(JLIST(I32X)) = FZTEMP(I32X)
                  FLRW(JLIST(I20X)) = FWTEMP(I20X) 
               ENDDO
C     
            ENDIF
         ENDDO
C     
C*****SCATTER THE FORCES BACK TO THEIR ORIGINAL POSITIONS.
C     
         DO JRC = 1, NR
            F((JRC-1)*4+1) = F((JRC-1)*4+1) + FLRX(JRC)
            F((JRC-1)*4+2) = F((JRC-1)*4+2) + FLRY(JRC)
            F((JRC-1)*4+3) = F((JRC-1)*4+3) + FLRZ(JRC)
            F((JRC-1)*4+4) = F((JRC-1)*4+4) + FLRW(JRC)
         ENDDO
C     
C         IF (LDOVIR) THEN
C            DO JRC = 1, MAXLEN
C               VIR(1) = VIR(1) + TVIR(JRC,1)
C               VIR(2) = VIR(2) + TVIR(JRC,2)
C               VIR(3) = VIR(3) + TVIR(JRC,3)
C            ENDDO
C            VIR(1) = VIR(1)*0.5E0
C            VIR(2) = VIR(2)*0.5E0
C            VIR(3) = VIR(3)*0.5E0
C         ENDIF
C     
C--   UPDATE THE 'EEL' AND 'ELJ' VALUES
         DO JRC1 = 1, NUNRE2
            DO JRC2 = 1, IVECT
               EEL(JRC1) = EEL(JRC1) + TEEL(JRC1,JRC2)
#ifndef EWALD
               ERF(JRC1) = ERF(JRC1) + TERF(JRC1, JRC2)
               ERC(JRC1) = ERC(JRC1) + TERC(JRC1, JRC2)
#endif
               ELJ(JRC1) = ELJ(JRC1) + TELJ(JRC1,JRC2)
               
            ENDDO
         ENDDO

      ENDIF

C     
C*****TRANSFORM THE COORDINATES AND FORCES BACK TO CARTESIAN COORDINATE
C*****WHEN REQUIRED.
C     
         IF (LDOTRA) THEN
            CALL TRACO (NR, 0, XCOORD, BETA, -1, LEVERY)
            CALL TRACO (NR, 0, F, BETA, -1, LEVERY)
         ENDIF
     
C         
 900     FORMAT (1X,A6,': ATOM',I5,' SOFT-CORE PERTURBED')
 910     FORMAT (1X,'AND C6 .EQ. 0 BUT C12 .NE. 0 WITH STATE ',
     $        A1,' ATOM TYPE',I5)
 915     FORMAT (1X,'AND CS6 .EQ. 0 BUT CS12 .NE. 0 WITH STATE ',
     $        A1,' ATOM TYPE',I5)
         
         RETURN
C     
C     
C*****ERROR EXIT
C     
C     
C*****DUMP FROM SUBROUTINE.
C     
C     
 5750 CONTINUE
C     
 9260 FORMAT (10(1X,I9))
 9280 FORMAT (10(1X,G11.4))
 9300 FORMAT (/,1X,20A)
 9320 FORMAT (1X,A,'(', I2, '):', 3(5X, 3(F8.3, 2X)))
 9340 FORMAT (1X, 10L10)
 9360 FORMAT (3(1X,I9),1X,G11.4)
C     
      PRINT '(A)', '1THE SUBROUTINE NONBML HAS GONE WRONG. ',
     $     ' THE VARIABLES ARE:'
      PRINT 9300, '      NPM       NRP      NCAG       NSM',
     $     '      NRAM       NTB      NSPM      NJNB'
      PRINT 9260, NPM, NRP, NCAG, NRAM, NTB, NSPM, NJNB
      PRINT 9300, '     BOX(1)      BOX(2)      BOX(3)',
     $     '        BETA'
      PRINT 9280, (BOX(K), K = 1, 3), BETA
      PRINT 9300, '   EEL(1..10)'
      PRINT 9280, (EEL(K), K = 1, NUNRE2)
      PRINT 9300, '   ELJ(1..10)'
      PRINT 9280, (ELJ(K), K = 1, NUNRE2)
      PRINT 9300, '   VIR(1..3)'
      PRINT 9280, (VIR(K), K = 1, 3)
      PRINT 9300, '    NSP(1..NSPM)'
      PRINT 9260, (NSP(K), K = 1, NSPM)
      PRINT 9300, '    INC(1..20)'
      PRINT 9260, (INC(K), K = 1, 20)
      PRINT 9300, '     CG(1..20)'
      PRINT 9280, (CG(K), K = 1, 20)
      PRINT 9300, '    CGS(1..NRAM)'
      PRINT 9280, (CGS(K), K = 1, NRAM)
      PRINT 9300, '      XCOORD(1..45)'
      DO 6550 I = 1, 45, 9
         PRINT 9320, 'X', I, (XCOORD(M), M = I, I + 8)
 6550 CONTINUE
      PRINT 9300, '      XC, YC, ZC(1..15)'
      DO 6555 I = 1, 15, 3
         PRINT 9320, 'COORDS', I, (XC(M),YC(M),ZC(M),
     $        M = I, I + 2)
 6555 CONTINUE
      PRINT 9300, '    INB(1..20)'
      PRINT 9260, (INB(K), K = 1, 20)
      PRINT 9300, '    JNB(1..200)'
      PRINT 9260, (JNB(K), K = 1, 200)
      PRINT 9300, '                  LOCAL VARIABLES:'
      PRINT 9300, '      LEVERY    LDOVIR     LERROR',
     $     '    LDOTRA    '
      PRINT 9340, LEVERY, LDOVIR, LERROR, LDOTRA
      PRINT 9300, '   COSB2'
      PRINT 9280, COSB2
      PRINT 9300, '     NRAG     NRAGT',
     $     '        NR        JN       JLEN'
      PRINT 9260, NRAG, NRAGT, NR, JN, JLEN
      PRINT 9300, '      ICE        NI         I        II',
     $     '        NN'
      PRINT 9260, ICE, NI, I, II, NN
      PRINT 9300, ' INAG(NI)  IFIRST(NI)'
      PRINT 9260, INAG(NI), IFIRST(NI)
      PRINT 9300, '  ICET(NI)  NSPTT(NI)  IACT(NI)  CGAT(NI)'
      PRINT 9360, ICET(NI), NSPTT(NI), IACT(NI), CGAT(NI)
      PRINT 9300, 'BOXH(1..3), BOXINV(1..3)'
      PRINT 9280, (BOXH(K), K = 1, 3), (BOXINV(K), K = 1, 3)
      CALL STOPGM('NONBML','FATAL ERROR')
C     ==--------------------------------------------------------------==
      RETURN      
      END
#endif
