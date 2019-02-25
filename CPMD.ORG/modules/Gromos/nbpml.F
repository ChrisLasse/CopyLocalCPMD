#ifndef __VECTOR
C Id: nbpml.f,v 1.36 1996/11/22 12:30:04 wscott Exp 
C
COMMSUBR NBPML
C Id: nbpml.f,v 1.36 1996/11/22 12:30:04 wscott Exp
C W.F. VAN GUNSTEREN, GRONINGEN, JUNE 1986
C 4D/soft-core/reaction-field by Ilario Tironi & Philippe Huenenberger
C Path-integral code by Salomon Billeter, Zuerich, Aug. 95
Cmb - Revised on 17 August 2006
C
C      SUBROUTINE NBPML(NATTOT,NPM,NSM,NRAGT,NCALCD,
C     $     XCOORD,XCG,XR,FLR,
C     $     NSPM,NSP,VIRLR,
C     $     INB,NSZPL,JNB,RCUTP,RCUTL,
C     $     MAXNRE,NRELKP,IAGRP,
C     $     RCRF,EPSRF,APPAK,
C     $     ELGLEL,ELGLRF,ELGLRC,ELGLLJ,
C     $     EL34EL,EL34RF,EL34RC,EL34LJ,
C     $     NUNRE2,ELREL,ELRRF,ELRRC,ELRLJ,
C     $     LDOPER,L3D4D,RLAM,RMU,ALPHLJ,ALPHC,NLAM,MMU,
C     $     LPION,LEVERY,LFREE)
C
C     NBPML makes a list of non-bonded pairs using a cut-off criterion
C     for the centres of the geometry of solute charge groups
C     and the first atoms of solvent molecules, by scanning all
C     possible pairs (using 3D distances only). The cut-off radius
C     for mating the list is L<RCUTP>.
C     The pair list contains ordered sequence numbers of solute charge atom
C     groups and solvent molecules (NI < NJ).
C     Periodic images can be included in the list. In this case, the
C     atoms of a charge group (solute) must lie within BOX/2 or
C     BOX*SQRT(3)/4 of each other.
C
C     If L<RCUTL> > L<RCUTP>, a potential is evaluated for the twin-range
C     region in subroutine L<NBWITH> during the pair list construction. 
C     Otherwise the pair list is constructed without long range forces 
C     in subroutine L<NBNONE>. The term long-range in contrast to GROMOS87
C     includes Lennard-Jones interactions.
C
C     The potential in L<NBWITH> reads:
C     V(r_ij) = - C6/(r_ij**6) + C12/(r_ij**12)
C               + q_i*q_j/r_ij - (0.5*C1/RCRF**3)*q_i*q_j*(r_ij**2)
C                              - ((1-0.5*C1)/RCRF)*q_i*q_j
C
C     where RCRF is the radius of the reaction field boundary (in principle
C     RCUTL) and 
C
C     C1 = ( (2-2*EPSRF)*(1+APPAK*RCRF) - EPSRF*(APPAK*RCRF)**2 ) /
C          ( (1+2*EPSRF)*(1+APPAK*RCRF) + EPSRF*(APPAK*RCRF)**2 )
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
C     L3D4D (NTG.EQ.L<NTGMU>.OR.NTG.EQ.L<NTGBOT>) will not be 
C     The cut-off criterion is always applied in 3D.
C
C     Periodic boundary conditions can be applied, depending on L<NTB>,
C     where the periodic box can be a truncated octahedron, rectangular
C     or monoclinic. 
C
C     The periodic box can be a truncated octahedron, rectangular or
C     monoclinic, depending on L<NTB>.
C     If the long range virial is to be calculated
C     (ABS(L<NTB> .EQ. L<NTBVIR>), it must be specified whether the
C     solute consists of separate submolecules or not.
C     This can be done by setting the NSPM > 1 and setting the corresponding
C     values into the NSP array. The virial is only calculated
C     in 3D.
C
C
C     NATTOT   total number of atoms
C     NPM      number of (identical) solute molecules
C     NRAGT    total number of charge groups
C     NCALCD   the number of dimensions to calculate distance and
C              forces in. Must be 3 or 4.
C              Note that the 3D distance between charge groups
C              I and J is _always_ used for deciding whether to put
C              the pair into the list or not, independently of NCALCD.
C              The distance used in the force calculation, however,
C              does depend on NCALCD.
C     XCOORD(NDIM*NATTOT)
C              the atom cartesian coordinates
C     XCG(NDIM*NRAGT)
C              the geometric centre of solute charge groups and the
C              first atoms of solvent molecules.
C     XR(NDIM*NATTOT)
C              atom cartesian coordinates relative to submolecular
C              centres of mass. This array is only used if
C              the virial is calculated, i.e. ABS(NTB) = L<NTBVIR>.
C              It is still accessed (read only) in the opposite case 
C              and thus should have the correct dimension (e.g. XCOORD)
C     FLR(NDIM*NATTOT)
C              delivered with the long range coulomb forces if
C              RCUTL > RCUTP.
C     NSPM     If the virial is to calculated, then L<NSPM> and L<NSP> are
C              accessed in order to calculate the molecular virial.
C              and thus the dummy array used in the call should have 
C              the correct dimension (e.g. XCOORD).
C     NSP(NSPM+1)
C              atom sequence numbers of the last atoms of the
C              submolecules.
C              only elements 1..NSPM are considered to be filled
C              with data. We need NSPM+1 elements in order to avoid
C              a double IF statement when checking for submolecules
C     
C     VIRLR(NDIM)
C              delivered with the X-, Y-, and Z- components of the
C              long range virial if L<NTB> = L<NTBVIR>.
C     INB(NRAGT)
C              delivered with the pointer list;
C              INB(I) specifies the number of charge groups or
C              solvent molecules J with J>I forming a pair with charge
C              group or solvent molecule I
C     NSZPL    delivered with the number of elements in JNB
C     JNB(MAXJNB)
C              delivered with the pointer list;
C              it specifies the solute charge groups or solvent
C              molecules J forming a pair with solute charge group
C              or solvent molecule I in ascending order; for each I the
C              pointer in JNB is incremented by INB(I)
C     MAXNRE   maximum number of energy groups
C     NRELKP(MAXNRE,MAXNRE)   
C              index of pairs of energy groups in energy arrays
C     IAGRP(MAXNAT)
C              for each atom, the energy type code its
C              contribution is added to.
C     RCRF     radius of the reaction field boundary (in principle
C              RCUTL) and 
C     EPSRF    dielectric permittivity of the reaction field 
C              continuum outside the boundary
C     APPAK    inverse Debye screening length
C     NUNRE2   number of energy froups
C     ELREL(MXNRE2)
C              delivered with long range coulomb energies.
C     ELRRF(MXNRE2)
C              delivered with the long range distance dependent
C              Poisson-Boltzmann reaction field energies
C     ELRRC(MXNRE2)
C              delivered with the long range distance independent
C              Poisson-Boltzmann reaction field energies
C     ELGLEL   delivered with dV/d(L<LAMBDA>) at fix L<MU> 
C              for the long range Coulomb term if LDOPER = .TRUE. 
C     ELGLRF   delivered with with dV/d(L<LAMBDA>) at fix L<MU> 
C              for the distance dependent long range 
C              Poisson-Boltzmann term if LDOPER = .TRUE. 
C     ELGLRC   delivered with with dV/d(L<LAMBDA>) at fix L<MU> 
C              for the distance independent long range 
C              Poisson-Boltzmann term if LDOPER = .TRUE. 
C     EL34EL   delivered with the 3D to 4D long range free energy contribution
C              for the Coulomb term if LDOPER = .TRUE. .AND. NCALCD = 4
C     EL34RF   delivered with the 3D to 4D long range free energy contribution
C              for the long range distance dependent 
C              Poisson-Boltzmann reaction field
C              term if LDOPER = .TRUE. .AND. NCALCD = 4
C     EL34RC   delivered with the 3D to 4D long range free energy contribution
C              for the long range distance independent 
C              Poisson-Boltzmann reaction field
C              term if LDOPER = .TRUE. .AND. NCALCD = 4
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
C     .FALSE.  a number of quantities depending on L<NPM>, L<NRP>,
C              L<NSM>, L<NRAM>, L<NCAG>, L<NTB>, L<BOX>, L<BETA>,
C              L<RCUTP> and L<RCUTL> are only calculated at the 
C              first call
C     .TRUE.   they are calculated at every subr. call.
C     LFREE    
C     .TRUE.   the pair list and interaction is calculated solely 
C              for perturbed atoms. This is used in free energy 
C              extrapolation.
C
C     NOT passed as arguments, but accessed in common blocks are
C     the following:
C
C     RCUTP    cut-off radius for the pair list; since no check on
C              neighbour exclusions is performed, and since only charge
C              group pairs I<J are scanned, L<RCUTP> should be chosen
C              large enough when L<RCUTL> > L<RCUTP>.
C     RCUTL    cut-off radius for the long range coulomb force;
C              The following conditions must hold, depending on the
C              boundary conditions specified.
C              If L<NTB> <0, L<RCUTP>, L<RCUTL> < BOX(M)*SQRT(3)/4
C              If L<NTB> >0, L<RCUTP>, L<RCUTL> < BOX(M)/2
C              If L<RCUTL> is larger, only the nearest image is calculated.
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
C
COMMEND
C     ==--------------------------------------------------------------==
      SUBROUTINE NBPML(NATTOT,NPM,NSM,NRAGT,NCALCD,
     $     XCOORD,XCG,XR,FLR,
     $     NSPM,NSP,VIRLR,
     $     INB,NSZPL,JNB,RCUTP,RCUTL,
     $     MAXNRE,NRELKP,IAGRP,
     $     RCRF,EPSRF,APPAK,
     $     ELGLEL,ELGLRF,ELGLRC,ELGLLJ,
     $     EL34EL,EL34RF,EL34RC,EL34LJ,
     $     NUNRE2,ELREL,ELRRF,ELRRC,ELRLJ,
     $     LDOPER,L3D4D,RLAM,RMU,ALPHLJ,ALPHC,NLAM,MMU,
     $     LPION,LEVERY,LFREE)
C     ==--------------------------------------------------------------==
C     INCLUDES
C     
      IMPLICIT NONE 

      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'pertsz.h'
      INCLUDE 'pertar.h'
      INCLUDE 'nbpml.h'
C     
C     ARGUMENTS
C     
      LOGICAL LDOPER,L3D4D,LPION,LEVERY,LFREE
      INTEGER NATTOT,NRAGT,NPM,NSM,NCALCD,NSZPL
      INTEGER INB(NRAGT),JNB(MAXJNB)
      INTEGER NSPM,NSP(NSPM+1)
      INTEGER MAXNRE,NRELKP(MAXNRE,MAXNRE)
      INTEGER NLAM,MMU
      INTEGER NUNRE2
      real*8 RCUTP,RCUTL,RCRF,EPSRF,APPAK
      real*8 XCOORD(NDIM*NATTOT),XCG(NDIM*NRAGT)
      real*8 XR(NDIM*NATTOT),FLR(NDIM*NATTOT)
      real*8 VIRLR(NDIM)
      real*8 ELREL(NUNRE2),ELGLEL,EL34EL
      real*8 ELRRF(NUNRE2),ELGLRF,EL34RF
      real*8 ELRRC(NUNRE2),ELGLRC,EL34RC
      real*8 ELRLJ(NUNRE2),ELGLLJ,EL34LJ
      real*8 RLAM,RMU,ALPHLJ,ALPHC
      CHARACTER IAGRP(NATTOT)
C     
C     LOCAL PARAMS
C     
      real*8 ENEPS
      PARAMETER (ENEPS = 1.D-20)
      real*8 ZEROC6, ZEROAL
      PARAMETER (ZEROC6 = 1.0D-20, ZEROAL = 1.0D-20)

C     
C     SAVED VARIABLES      
C     
      LOGICAL LFIRST,LDOLNG
      SAVE LFIRST,LDOLNG
C     
C     LOCAL VARIABLES
C     
      INTEGER NATTO3,I,II,IFCG,ILCG,ILAST,J
      INTEGER ITMPA,ITMPB,JTYPE
      INTEGER INTIJA,INTIJB
      REAL*8 CONV,BETAR,PAPPAK
      CHARACTER PRGSTR*(5)

C     
C     DATA
C     
      DATA LFIRST/.TRUE./
      DATA PRGSTR  /'NBPML'/
C     
C------------BEGIN NBPML
C     
C------------DO FIRST CALL SETUP (SWITCHES,"MONOCLINIC" VARIABLES
C     REACTION FIELD VARIABLES, SIMPLE CONSISTENCY CHECKS)
C     

C
C---- IF L3D4D IS ON, WE HAVE TO GO THROUGH THE PERTURBATION CODE
C     ANYWAY
C     

C
      IF (LFIRST .OR. LEVERY) THEN

        LMONO = (NTB .GT. 0)
        LOCTO = (NTB. LT. 0)
        LVAC  = (NTB .EQ. 0)
        LDOVIR= (ABS(NTB) .EQ. NTBVIR)
C
        LDOTRA = .FALSE.
C
        IF (LMONO)THEN
          CONV = 45.0D0/DATAN(1.0D0)
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
C
        LDOLNG = (RCUTL .GT. RCUTP)        
C
      ENDIF
      
      IF (LFIRST.OR.LEVERY) THEN
        LFIRST = .FALSE.
        L4D = (NCALCD .EQ. 4)
        LPIDOP = .TRUE.
        IF (NPID .NE. 0) THEN
           PININV = 1.0D0 / NPID
        ENDIF
        LPERTL = (LDOPER.OR.L3D4D)
C     
C---- PERFORM A FEW FIRST CALL CHECKS
C     
C   
C--   IF LPERTL AND ALPHLJ.NE.0.0 CHECK THAT C6.EQ.0 => C12.EQ.0 FOR ATOMS II 
C     WITH ISCLJ(II).NE.0 AND ANY OTHER ATOM ( STATE A AND B )
C     
        IF (LPERTL .AND. DABS(ALPHLJ) .GT. ZEROAL) THEN
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
       ENDIF

C     NOW REACTION FIELD VARIABLES
C     
       IF (EPSRF .LE. 0.0D0) THEN
          C1 = -1.0D0
          PRINT *,PRGSTR,': WARNING, EPSRF.LE.0 , PERMITTIVITY SET ',
     $         'TO INFINITY'
       ELSE
          PAPPAK=APPAK*RCRF
          C1 = ((2.0D0-2.0D0*EPSRF)*(1.0D0+PAPPAK) 
     .         - EPSRF* PAPPAK**2)
     .        /((1.0D0+2.0D0*EPSRF)*(1.0D0+PAPPAK) 
     .         + EPSRF* PAPPAK**2)
       ENDIF
       RFF = C1 / RCRF**3
       RFE = -0.5D0*RFF
       RFC = -(1.0D0-0.5D0*C1)/RCRF
       RCRF2 = RCRF*RCRF
C
C--   CALCULATING IN 4D WITH 3D COORDINATES IS WRONG
C     
       IF ( NCALCD.GT.NDIM ) THEN
          PRINT *,PRGSTR,' ERROR: NCALCD.GT.NDIM'
          CALL STOPGM('NBPML','FATAL ERROR')
       ENDIF
C     
C--   CALCULATING IN 3D DOESNT ALLOW THE CALCULATION OF MU DERIVATIVE
C     
       IF (.NOT. L4D .AND. L3D4D) THEN
          PRINT *,PRGSTR,': NCALCD.EQ.3.AND.',
     $         '(NTG.EQ.NTGMU.OR.NTG.EQ.NTGBOT) NOT ALLOWED'
       ENDIF

       IF (LDOVIR) THEN
C     SET UP GLOBAL SUB-MOLECULE ARRAY
C     DENOTE SOLVENT MOLECULES BY -1
          IF (NSPM .GT. 0) THEN
             IFCG = 1
             ILCG = 1
             DO I = 1, NSPM
                ILAST = NSP(I)
 1880           IF (ILAST .GT. INC(ILCG)) THEN
                   ILCG = ILCG + 1
                   GOTO 1880
                ENDIF
                IF (ILAST .NE. INC(ILCG)) THEN
                   WRITE(6,*) 'NSP(I) ARRAY NOT CONSISTENT',
     $                  ' WITH CHARGE GROUP DEFINITION !'
                   WRITE(6,*) 'ERROR IN NBPML.F'
                   CALL STOPGM('NBPML','FATAL ERROR')
                ENDIF
!$OMP parallel do private(J)
                DO J = IFCG, ILCG
                   NSPT(J) = NSP(I)
                ENDDO
                IFCG = ILCG + 1
                ILCG = IFCG
             ENDDO
          ELSE
!$OMP parallel do private(I)
             DO I = 1, NCAG
                NSPT(I) = 1
             ENDDO
          ENDIF
          DO J = 2, NPM
             DO I = 1, NCAG
                NSPT((J - 1) * NCAG + I) = NSPT(I) + 
     $               (NSPT(NCAG) * (J - 1))
             ENDDO
          ENDDO
C     SOLVENT MOLECULES
          ILCG = NPM * NCAG
!$OMP parallel do private(I)
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
C---- DONE FIRST CALL CHECKS
C     
C--   CALCULATING IN 3D WITH RMU.NE.0 SHOULDNT BE DONE
C     
      IF ( .NOT.L4D.AND.RMU.NE.0.0D0 ) THEN
         PRINT *,PRGSTR,': WARNING, NCALCD .EQ. 3 .AND. RMU .NE. 0.0'
         PRINT *,'         LOCAL RMU WILL BE SET TO ZERO'
         RMULOC = 0.0D0
      ELSE
         RMULOC = RMU
      ENDIF

      RCUTP2 = RCUTP*RCUTP
      RCUTL2 = RCUTL*RCUTL
      NATTO3 = NDIM *NATTOT

C     
C------------DONE FIRST CALL SETUP
C     TRANSFORM COORDINATES IF REQUIRED
C     
      IF (LDOTRA) THEN
        CALL TRACO(NATTOT,0,XCOORD,BETA,1,LEVERY)
      ENDIF
C     
C------------INITIALIZE ARRAYS AND GENERAL VARIABLES
C     

C     
      ELGLEL = 0.0D0
      ELGLRF = 0.0D0
      ELGLRC = 0.0D0
      ELGLLJ = 0.0D0
         
      EL34EL = 0.0D0
      EL34RF = 0.0D0
      EL34RC = 0.0D0
      EL34LJ = 0.0D0
C     
      CALL mm_AZZERO(FLR,NATTO3)
      CALL mm_AZZERO(VIRLR,NDIM)
C     
!$OMP parallel do private(I)
      DO I=1,NUNRE2
         ELREL(I) = 0.0D0
         ELRRF(I) = 0.0D0
         ELRRC(I) = 0.0D0
         ELRLJ(I) = 0.0D0
      ENDDO
C     
C------------DONE INITIALIZE
C     GENERATE ARRAY XCG CONTAINING THE REFERENCE POINT OF 
C     EACH CHARGE GROUP FOR DISTANCE MEASUREMENT.
C     SOLUTE : CENTER OF GEOMETRY, SOLVENT : FIRST ATOM
C     
      CALL CLCGEO(NATTOT,NRAGT,NPM,NSM,NDIM,XCOORD,XCG)
C     
C------------DONE INITIALIZE
C     START CALCULATION ( WITH OR WITHOUT LONG RANGE FORCES )
C     
      IF (LDOLNG .OR. LFREE) THEN
        CALL NBWITH(NATTOT,NPM,NRAGT,NCALCD,
     $       XCG,XCOORD,XR,FLR,
     $       VIRLR,
     $       INB,JNB,NSZPL,
     $       MAXNRE,NRELKP,IAGRP,
     $       NUNRE2,
     $       ELREL,ELGLEL,EL34EL,
     $       ELRRF,ELGLRF,EL34RF,
     $       ELRRC,ELGLRC,EL34RC,
     $       ELRLJ,ELGLLJ,EL34LJ,
     $       LDOPER,L3D4D,RLAM,
     $       ALPHLJ,ALPHC,NLAM,MMU,
     $       LPION,LEVERY,LFREE)
      ELSE
        CALL NBNONE(NRAGT,3,XCG,INB,JNB,NSZPL)
      ENDIF
C     
C------------DONE CALCULATION
C     TRANSFORM BACK COORDINATES IF REQUIRED
C     
      IF (LDOTRA) THEN
        CALL TRACO(NATTOT,0,XCOORD,BETA,-1,LEVERY)
        CALL TRACO(NATTOT,0,FLR,BETA,-1,LEVERY)
      ENDIF
C     
 900  FORMAT (1X,A6,': ATOM',I5,' SOFT-CORE PERTURBED')
 910  FORMAT (1X,'AND C6 .EQ. 0 BUT C12 .NE. 0 WITH STATE ',
     $     A1,' ATOM TYPE',I5)
 915  FORMAT (1X,'AND CS6 .EQ. 0 BUT CS12 .NE. 0 WITH STATE ',
     $     A1,' ATOM TYPE',I5)
C------------END NBPML
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE CLCGEO(NATTOT,NRAGT,NPM,NSM,NDIM,XSRC,XDST)
C     ==--------------------------------------------------------------==
COMMSUBR CLCGEO
C     SUBROUTINE CLCGEO(NATTOT,NRAGT,NPM,NSM,NDIM,XSRC,XDST)
C
C     CLCGEO calulates the geometric centres of
C     all CG of the solute and the first atom coordinate of
C     solute molecules.
C     This routine is called by L<NBPML>
COMMEND
C
C     INCLUDES
C
      INCLUDE 'coordsz.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
C
C     ARGUMENTS
C
      INTEGER NATTOT,NRAGT,NPM,NSM,NDIM
      real*8 XSRC(NDIM*NATTOT),XDST(NDIM*NRAGT)
C
C     LOCAL VARIABLES
C
      INTEGER M,J,NSRC3,NDST3,NCG,NPP
      INTEGER NINC,ILACG,IFICG,NUMCG
      real*8 XTEMP(MAXDIM)
C
C------------ BEGIN CLCGEO
C
      NSRC3 = 0
      NDST3 = 0
C
C------------ FIRST DO THE SOLUTE ( CENTRE OF GEOMETRY )
C
      DO 10 NPP=1,NPM
         IFICG = 1
         DO 20 NCG=1,NCAG
            ILACG = INC(NCG)
            NUMCG = ILACG - IFICG + 1
            IF (NUMCG .GT. 1) THEN
!$OMP parallel do private(M)
               DO M=1,NDIM
                  XTEMP(M) = 0.0D0
               ENDDO
C     
               DO J=1,NUMCG
                  DO M=1,NDIM
                     XTEMP(M) = XTEMP(M) + XSRC(NSRC3+M)
                  ENDDO
                  NSRC3 = NSRC3 + NDIM
               ENDDO
C     
               DO M=1,NDIM
                  XDST(NDST3+M) = XTEMP(M)/NUMCG
               ENDDO
            ELSE
               DO M=1,NDIM
                  XDST(NDST3+M) = XSRC(NSRC3+M)
               ENDDO
               NSRC3 = NSRC3 + NDIM
            ENDIF
            NDST3 = NDST3 + NDIM
            IFICG = ILACG + 1
 20      CONTINUE
 10   CONTINUE
C
C------------ NOW DO THE SOLVENT (COORDINATE OF THE FIRST ATOM)
C
      NINC = NDIM*NRAM
      DO NPP=1,NSM
         DO M=1,NDIM
            XDST(NDST3+M) = XSRC(NSRC3+M)
         ENDDO
         NDST3 = NDST3 + NDIM
         NSRC3 = NSRC3 + NINC
      ENDDO
C
C------------ END CLCGEO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE NBNONE(NRAGT,NCALCD,X,INB,JNB,NJPTR)
C     ==--------------------------------------------------------------==
COMMSUBR NBNONE
C     SUBROUTINE NBNONE(NRAGT,NCALCD,X,INB,JNB,NJPTR)
C
C     NBNONE calculates a pairlist without calculating
C     a long range force.
C     This routine is called by L<NBPML>
COMMEND
C
C     INCLUDES
C
      INCLUDE 'coordsz.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'box.h'
      INCLUDE 'nbpml.h'
C
C     ARGUMENTS
C
      INTEGER NRAGT,NCALCD,INB(NRAGT),JNB(MAXJNB),NJPTR
      real*8 X(NRAGT*NDIM)
C
C     LOCAL VARIABLES
C
      INTEGER IX3,JX3,NI,NJ,M,JSAVE
      real*8 RIJCG2,DSTTMP,XIJ(MAXDIM),XCGI(MAXDIM)
      CHARACTER PRGSTR*(6)
C
C     DATA
C
      DATA PRGSTR /'NBNONE'/
C
C------------ BEGIN NBNONE
C
      NJPTR = 0
      IX3 = 0
C
C------------ LOOP OVER PRIMARY CHARGE GROUPS      
C
      DO 100 NI=1, NRAGT - 1
         JSAVE = NJPTR
         DO M=1,NDIM
            XCGI(M) = X(IX3+M)
         ENDDO
         IX3 = IX3 + NDIM

         JX3 = IX3
C
C-------- LOOP OVER SECONDARY CHARGE GROUPS (EITHER SOLUTE OR SOLVENT)
C
         DO 200 NJ=NI+1, NRAGT
C     
C---- CALCULATE DISTANCE AND CONNECTING VECTOR BETWEEN XCGI AND XCG(NJCOR)
C     CORRECT FOR BOUNDARY CONDITIONS IF NECESSARY
C
            DO M=1,NCALCD
               XIJ(M) = XCGI(M) - X(JX3 + M)
            ENDDO
            JX3 = JX3 + NDIM
            RIJCG2 = 0.0D0
            IF (LVAC) THEN
               DO M=1,NCALCD
                  RIJCG2 = RIJCG2 + XIJ(M)**2
               ENDDO
            ELSE
               DO M=1,NCALCD
                  DSTTMP = XIJ(M)
                  IF (DSTTMP .GE. BOXH(M)) THEN
                     DSTTMP = DSTTMP - BOX(M)
                  ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
                     DSTTMP = DSTTMP + BOX(M)
                  ENDIF
                  XIJ(M) = DSTTMP
                  RIJCG2 = RIJCG2 + DSTTMP*DSTTMP
               ENDDO
C
               IF (LOCTO) THEN
                  DSTTMP = BOXOQ - DABS(XIJ(1))-
     .                DABS(XIJ(2))-DABS(XIJ(3))
                  IF (DSTTMP .LT. 0.0D0) THEN
                     RIJCG2 = RIJCG2 + DSTTMP*BOX(1)
                     DO M=1,NCALCD
                        XIJ(M) = 
     .                       XIJ(M) - SIGN(BOXOH,XIJ(M))
                     ENDDO
                  ENDIF
               ELSEIF (LMONO .AND. LDOTRA) THEN
                  RIJCG2 = RIJCG2 + COSB2*XIJ(1)*XIJ(3)
               ENDIF
            ENDIF
C
C---- DONE DISTANCE CALCULATION
C     CHECK DISTANCE
C
            IF (RIJCG2 .LE. RCUTP2) THEN
C
C---- INSIDE SHORT-RANGE CUT-OFF, STICK INTO PAIRLIST
C
               NJPTR = NJPTR + 1
               IF (NJPTR .GT. MAXJNB) THEN
                  PRINT *,PRGSTR,' PAIRLIST OVERFLOW:'
                  PRINT *,'MAXIMUM SIZE: MAXJNB =',MAXJNB
                  PRINT *,'              NJPTR  =',NJPTR ! cmb
                  CALL STOPGM('NBPML','FATAL ERROR')
               ENDIF
               JNB(NJPTR) = NJ
            ENDIF
C     
C---------- END LOOP OVER NJ SECONDARY CHARGE GROUPS
C
 200     CONTINUE

         INB(NI) = NJPTR - JSAVE
C
C------------ END LOOP OVER PRIMARY CHARGE GROUPS      
C
 100  CONTINUE
C
      INB(NRAGT) = 0

      PRINT 500,NJPTR
C
C------------ END NBNONE
C
 500  FORMAT (1X,'NBPML:',I10,' ELEMENTS IN THE PAIRLIST')
C nbnone
C     ==--------------------------------------------------------------==
      RETURN
      END 
C     ==================================================================
      SUBROUTINE NBWITH(NATTOT,NPM,NRAGT,NCALCD,
     $     XCG,XCOORD,XR,FLR,
     $     VIRLR,
     $     INB,JNB,NJPTR,
     $     MAXNRE,NRELKP,IAGRP,
     $     NUNRE2,
     $     ELREL,ELGLEL,EL34EL,
     $     ELRRF,ELGLRF,EL34RF,
     $     ELRRC,ELGLRC,EL34RC,
     $     ELRLJ,ELGLLJ,EL34LJ,
     $     LDOPER,L3D4D,RLAM,
     $     ALPHLJ,ALPHC,NLAM,MMU,
     $     LPION,LEVERY,LFREE)
C     ==--------------------------------------------------------------==
C     NBWITH constructs a charge group pairlist and calculates a
C     long range force.
C     NBWITH is called by L<NBPML>.
C     ==--------------------------------------------------------------==
C     
C     INCLUDES
C     
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'pertsz.h'
      INCLUDE 'pertar.h'
      INCLUDE 'nbpml.h'
C     
C     LOCAL PARAMS
C     
      real*8 ENEPS
      PARAMETER (ENEPS = 1.D-20)      
C     
C     ARGUMENTS
C     
      LOGICAL LEVERY,LDOPER,L3D4D,LPION,LFREE
      INTEGER NRAGT,NCALCD,INB(NRAGT)
      INTEGER JNB(MAXJNB),NJPTR,NATTOT,NPM
      INTEGER NUNRE2
      INTEGER MAXNRE,NRELKP(MAXNRE,MAXNRE)
      INTEGER NLAM,MMU
      real*8  RLAM,ALPHLJ,ALPHC
      real*8 XCG(NDIM*NRAGT)
      real*8 XCOORD(NDIM*NATTOT),FLR(NDIM*NATTOT),XR(NDIM*NATTOT)
      real*8 VIRLR(NDIM)
      real*8 ELREL(NUNRE2),ELGLEL,EL34EL
      real*8 ELRRF(NUNRE2),ELGLRF,EL34RF
      real*8 ELRRC(NUNRE2),ELGLRC,EL34RC
      real*8 ELRLJ(NUNRE2),ELGLLJ,EL34LJ
      CHARACTER IAGRP(NATTOT)

C     
C     LOCAL VARIABLES
C     
      LOGICAL LJSOLV,LDOIT
      LOGICAL LIZERO, LJZERO
      LOGICAL LFULL4,LPART4,LFULL3,LPART3
      LOGICAL LFCGI,LGENPA
      INTEGER NRAGP,NRPT,NIPROT,M,IAUX,ITMP,JSAVE,JTOPCG,NJPROT
      INTEGER NITOP,NICOR,II
      INTEGER IFTOPA,ILTOPA,IFCORA,ILCORA,INUMAT,INUMA3
      INTEGER JFTOPA,JLTOPA,JFCORA,JLCORA,JNUMAT,JNUMA3
      INTEGER IATOFF,ICOPTR,ICGPTR
      INTEGER JATOFF,JCOPTR,JCGPTR
      INTEGER NJCOR,JJ
      INTEGER ILPICI
      real*8 RIJ2,RIJCG2,DSTTMP,XIJ(MAXDIM),XCGI(MAXDIM)
      real*8 PIFACT
C     
      real*8 ZEROC6
      PARAMETER (ZEROC6 = 1.0D-20)
C     
C     ISUBM, JSUBM ARE ONLY USED IF WE HAVE TO CALCULATE
C     THE VIRIAL AND WE HAVE A SOLUTE WITH SUBMOLECULES
C     
      CHARACTER PRGSTR*(6)

C     
C     LOCAL VARIABLES
C     
      INTEGER IM,ILOC3,JLOC3
      INTEGER ILGRP,JLGRP,NDXGRP
      INTEGER ILATYP,ILBTYP,JLATYP,JLBTYP
      INTEGER INTIJA,INTIJB
      real*8 DF3,DF4TH
      real*8 SOFTLJ, CA126, CB126, SFLJ
      real*8 RIJ4, RIJ6
C     
C     LOCAL VARIABLES FOR IMPROVED PERFORMANCE
C     
      LOGICAL LLIPER(MXATCG),LLJPER(MXATCG),LIPERT,LJPERT
      INTEGER ITYPEA(MXATCG),ITYPEB(MXATCG)
      INTEGER ISCLJL(MXATCG),JSCLJL(MXATCG),ISCLJI, ISCLJJ
      INTEGER JTYPEA(MXATCG),JTYPEB(MXATCG)
C     
C     ENERGY GROUPS OF THE ATOMS OF THE CG
C     
      INTEGER IGRP(MXATCG),JGRP(MXATCG)
      INTEGER ISCELL(MXATCG), JSCELL(MXATCG), ISCELI, ISCELJ
      INTEGER ILPIC(MXATCG), ILPJC(MXATCG)
      real*8 CGIA(MXATCG),CGIB(MXATCG)
      real*8 CGJA(MXATCG),CGJB(MXATCG)
      real*8 XIAT(MAXLAT),XJAT(MAXLAT)
      real*8 FIAT(MAXLAT),FJAT(MAXLAT)
      real*8 XRIAT(MAXLAT),XRJAT(MAXLAT)
C     
C     PERTURBATION VARIABLES
C     
      real*8 CGILA,CGJLA,CGILB,CGJLB
      real*8 XIJDST(MAXDIM),XH
      real*8 RIJ3D2, VELA4D, VRFA4D, DFA4D, DFRFA4
      real*8 VLJA4D
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
      real*8 QIQJA,QIQJB,QRFR
      real*8 SOFTEL, SFEL
      real*8 RSOFT, RSOFT2, RRF, RRF2, RFSQRT
      real*8 RLMA, RL2A, RLMB, RL2B, RLLA, RLLB, RLB, RLA, RDLB, RDLA
      real*8 RLARMU, RLBRMU
      real*8 RM4D, RM3D, RDM3D, RDM4D
      real*8 RMCST, RMDER
      real*8 RIJINV, RIJIN2, RIJIN6
      real*8 CA12,CA6,CB12,CB6
      real*8 CRA1,CRA12
C     
C     DATA
C     
      DATA PRGSTR /'NBWITH'/
C     
C------------BEGIN NBWITH
C     
C     
C------------CALCULATE RMU/RLAMBDA DEPENDANT VARIABLES
C     
      IF (LPERTL) THEN
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
          RDLB = RLAM**(NLAM - 1) * NLAM
        ENDIF
        IF (RLAM.EQ.1.0D0 .AND. NLAM.EQ.1) THEN
          RLA = 0.0D0
          RDLA = -DBLE(NLAM)
        ELSE
          RLA = (1.0D0 - RLAM)**NLAM          
          RDLA = - (1.0D0 - RLAM)**(NLAM - 1) * NLAM
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
          RDM4D = RMULOC**(MMU - 1) * MMU
        ENDIF
        IF (RMULOC.EQ.1.0D0 .AND. MMU.EQ.1) THEN
          RM3D = 0.0D0
          RDM3D = -DBLE(MMU)
        ELSE
          RM3D = (1.0D0 - RMULOC)**MMU
          RDM3D = - (1.0D0 - RMULOC)**(MMU - 1) * MMU
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
C--   YET A FEW MORE LOGICALS FOR SPEED
C     
        LFULL4 = (L4D .AND. (RMULOC .NE. 0.0D0))
        LPART4 = (RMULOC .EQ. 0.0D0 .AND. L3D4D .AND. MMU.EQ.1)
        LFULL3 = ((L4D .AND. (RMULOC .NE. 1.0D0)).OR. (.NOT. L4D))
        LPART3 = (RMULOC .EQ. 1.0D0 .AND. L3D4D .AND. MMU.EQ.1)
C
      ENDIF
C     
C------------INITIALIZE A FEW VARIABLES AND POINTERS
C     NRAGP: THE TOTAL NUMBER OF SOLUTE CHARGE GROUPS
C     NRPT : THE TOTAL NUMBER OF SOLUTE ATOMS
C     NJPTR IS USED TO ACCESS THE PAIRLIST
C     IATOFF IS THE ATOM SEQUENCE NUMBER OFFSET FOR 'I ATOMS'
C     ICOPTR IS THE COORDINATE ARRAY OFFSET GOT 'I ATOMS'
C     ICGPTR IS THE OFFSET FOR THE COORDS OF CG I
C     NICOR IS USED FOR INDEXING INB
C     
      NRAGP = NCAG*NPM
      NRPT = NRP*NPM
C     
      ICOPTR = 0
      IATOFF = 0
      NJPTR = 0
      ICGPTR = 0
      NICOR  = 1
C
C-- BY DEFAULT, WE CONSIDER ALL PAIRS
C
      LGENPA = .TRUE.
C     
C------------DONE INITIALIZE
C     START CALCULATION
C     
C------------LOOP OVER PROTEIN MOLECULES AND PRIMARY (PROTEIN ONLY)
C     CHARGE GROUPS      
C     
      DO 17 NIPROT=1,NPM
        IFTOPA = 1
        IFCORA = IFTOPA + IATOFF
        DO 100 NITOP=1,NCAG
          ILTOPA = INC(NITOP)
          ILCORA = ILTOPA + IATOFF
          INUMAT = ILTOPA - IFTOPA + 1
          INUMA3 = INUMAT*NDIM
C     
C---- PER PRIMARY CHARGE GROUP USE LOCAL ARRAYS
C     (TOPOLOGY INFORMATION, COORDINATES, FORCES , CENTER OF GEOMETRY)
C     DETERMINE SOLUTE SUBMOLECULES IF NECESSARY
C     
          IM = 1
          LFCGI = .FALSE.
          DO 11 II = IFTOPA,ILTOPA
            CGIA(IM)   = CG(II)
            ITYPEA(IM) = IAC(II)     
            ILPIC(IM)  = IPIC(II)
            IAUX       = IPERT(II)
            IF (LPERTL .AND. IAUX .NE. NOPERT) THEN
              LLIPER(IM) = .TRUE.
              LFCGI      = .TRUE.
              CGIB(IM)   = CGB(IAUX)
              ITYPEB(IM) = IACB(IAUX)
              ISCLJL(IM) = 1 - ISCLJ(IAUX)
              ISCELL(IM) = 1 - ISCC(IAUX)
            ELSE
              LLIPER(IM) = .FALSE.
              CGIB(IM)   = CGIA(IM)
              ITYPEB(IM) = ITYPEA(IM)
              ISCLJL(IM) = 1
              ISCELL(IM) = 1
            ENDIF
            IM = IM+1
 11       CONTINUE
C     
          IM = 1
          DO II=IFCORA,ILCORA
            IGRP(IM) = ICHAR(IAGRP(II))
            IM = IM + 1
          ENDDO
!$OMP parallel do private(II) shared(ICOPTR)
          DO II=1,INUMA3
            XIAT(II)  = XCOORD(ICOPTR+II)
            XRIAT(II) = XR(ICOPTR+II)
            FIAT(II)= 0.0D0
          ENDDO
C     
!$OMP parallel do private(M) shared(ICGPTR)
          DO M=1,NDIM
            XCGI(M) = XCG(ICGPTR+M)
          ENDDO
C     
C---- DONE WITH SETTING UP THE LOCAL ARRAYS
C     
          JSAVE = NJPTR
          JCGPTR = ICGPTR+NDIM
C     
C--------LOOP OVER SECONDARY CHARGE GROUPS (EITHER SOLUTE OR SOLVENT)
C     
          DO 201 NJCOR=NICOR+1, NRAGT
C     
C-----IF WE ARE DOING FREE ENERGY PERTURBATION ONLY CONTINUE IF 
C     CHARGE GROUP I OR CHARGE GROUP J CONTAIN PERTURBED ATOMS 
C     (AVOIDS DUPLICATION OF ROUTINES FOR USE WITH PROFEE) 
C     
            IF(LFREE) THEN 
              LGENPA = .FALSE.
              LJSOLV = (NJCOR .GT. NRAGP)
              IF (LFCGI) THEN
                LGENPA = .TRUE.
              ELSEIF (.NOT.LJSOLV)THEN
                JTOPCG = MOD(NJCOR-1,NCAG) + 1
                IF (JTOPCG .EQ. 1) THEN
                  JFTOPA = 1
                ELSE
                  JFTOPA = INC(JTOPCG-1) + 1
                ENDIF                
                JLTOPA = INC(JTOPCG)
                DO JJ=JFTOPA,JLTOPA
                  IAUX    = IPERT(JJ)
                  IF (IAUX .NE. NOPERT) THEN
                    LGENPA = .TRUE.
                  ENDIF
                ENDDO
              ENDIF
            ENDIF
C
            IF(LGENPA) THEN            

               LDOIT = (LDOVIR .AND. 
     $              (NSPT(NICOR) .NE. NSPT(NJCOR)))

C     
C---- THE CHARGE GROUP PAIR IS CONSIDERED.
C     CALCULATE DISTANCE AND CONNECTING VECTOR BETWEEN XCGI AND XCG(NJCOR)
C     CORRECT FOR BOUNDARY CONDITIONS IF NECESSARY
C     
              XIJ(1) = XCGI(1) - XCG(JCGPTR + 1)
              XIJ(2) = XCGI(2) - XCG(JCGPTR + 2)
              XIJ(3) = XCGI(3) - XCG(JCGPTR + 3)
              RIJCG2 = 0.0D0
              IF (LVAC) THEN
                RIJCG2=XIJ(1)*XIJ(1)+XIJ(2)*XIJ(2)+XIJ(3)*XIJ(3)
              ELSE
                DO M=1,3
                  DSTTMP = XIJ(M)
                  IF (DSTTMP .GE. BOXH(M)) THEN
                    DSTTMP = DSTTMP - BOX(M)
                  ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
                    DSTTMP = DSTTMP + BOX(M)
                  ENDIF
                  XIJ(M) = DSTTMP
                  RIJCG2 = RIJCG2 + DSTTMP*DSTTMP
                ENDDO
C     
                IF (LOCTO) THEN
                  DSTTMP = BOXOQ -DABS(XIJ(1))-
     .               DABS(XIJ(2))-DABS(XIJ(3))
                  IF (DSTTMP .LT. 0.0D0) THEN
                    RIJCG2 = RIJCG2 + DSTTMP*BOX(1)
                    XIJ(1) = XIJ(1) - SIGN(BOXOH,XIJ(1))
                    XIJ(2) = XIJ(2) - SIGN(BOXOH,XIJ(2))
                    XIJ(3) = XIJ(3) - SIGN(BOXOH,XIJ(3))
                  ENDIF
                ELSEIF (LMONO .AND. LDOTRA) THEN
                  RIJCG2 = RIJCG2 + COSB2*XIJ(1)*XIJ(3)
                ENDIF
              ENDIF
C
C---- DONE DISTANCE CALCULATION
C     CHECK DISTANCE
C     
              IF (RIJCG2 .LE. RCUTP2) THEN
C     
C---- INSIDE SHORT-RANGE CUT-OFF, STICK INTO PAIRLIST
C     
                NJPTR = NJPTR + 1
                IF (NJPTR .GT. MAXJNB) THEN
                  PRINT *,PRGSTR,' PAIRLIST OVERFLOW:'
                  PRINT *,'MAXIMUM SIZE: MAXJNB =',MAXJNB
                  PRINT *,'              NJPTR  =',NJPTR ! cmb
                  PRINT *,NICOR,NJCOR
                  CALL STOPGM('NBPML','FATAL ERROR')
                ENDIF
                JNB(NJPTR) = NJCOR
              ELSEIF (RIJCG2 .LE. RCUTL2) THEN
C     
C---- WITHIN THE TWIN RANGE, CALCULATE ELECTROSTATIC INTERACTION
C     
                LJSOLV = (NJCOR .GT. NRAGP)
C     
C---- PER SECONDARY CHARGE GROUP USE LOCAL ARRAYS
C     (TOPOLOGY INFORMATION), DISTINGUISH NJ SOLVENT OR SOLUTE
C     
                IF (LJSOLV) THEN
C     
C--   NJ IS A SOLVENT CG
C     
                  JFTOPA = 1
                  JLTOPA = NRAM
                  JFCORA = NRPT + NRAM*(NJCOR-NRAGP-1) + 1
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
                    JSCELL(IM) = 1
                    JSCLJL(IM) = 1
                    ILPJC(IM)  = 0
                    IM = IM+1
                  ENDDO
                ELSE
C     
C--   NJ IS A PROTEIN CHARGE GROUP
C     
                  NJPROT = ((NJCOR-1) / NCAG) + 1
                  JTOPCG = MOD(NJCOR-1,NCAG) + 1
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
                  IM = 1
                  ITMP = JFCORA-1
                  DO JJ=JFTOPA,JLTOPA
                    CGJA(IM)   = CG(JJ)
                    JTYPEA(IM) = IAC(JJ)
                    JGRP(IM)   = ICHAR(IAGRP(ITMP+IM))
                    ILPJC(IM)  = IPIC(JJ)
                    IAUX    = IPERT(JJ)
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
                JCOPTR = NDIM*(JFCORA-1)
!$OMP parallel do private(JJ)
                DO JJ=1,JNUMA3
                  XJAT(JJ)  = XCOORD(JCOPTR+JJ)
                  XRJAT(JJ) = XR(JCOPTR+JJ)
                  FJAT(JJ)= 0.0D0
                ENDDO
C     
C---- DONE WITH SETTING UP THE LOCAL ARRAYS
C     
                ILOC3 = 0
C     
C--------LOOP OVER ATOMS II OF CHARGE GROUPS NI
C     
                DO 301 II=1,INUMAT
                  CGILA = CGIA(II)
                  CGILB = CGIB(II)
                  ILGRP = IGRP(II)
                  ILATYP = ITYPEA(II)
                  ILBTYP = ITYPEB(II)
                  LIPERT = LLIPER(II)
                  ISCLJI = ISCLJL(II)
                  ISCELI = ISCELL(II)
                  ILPICI = ILPIC(II)
C     
                  JLOC3 = 0
C     
C------LOOP OVER ATOMS JJ OF CHARGE GROUP NJ
C     
                  DO 400 JJ=1,JNUMAT
                    CGJLA = CGJA(JJ)
                    CGJLB = CGJB(JJ)
                    JLGRP = JGRP(JJ)
                    LJPERT = LLJPER(JJ)
                    JLATYP = JTYPEA(JJ)
                    INTIJA = MPAC(ILATYP,JLATYP)
                    ISCELJ = JSCELL(JJ)
                    ISCLJJ = JSCLJL(JJ)
                    CA12  =  C12(INTIJA)
                    CA6  =  C6(INTIJA)
C     
C--------CHECK FOR PATH INTEGRAL EXCLUSIONS
C     
                    IF (LPION) THEN
                      LIZERO = (ILPICI .EQ. 0)
                      LJZERO = (ILPJC(JJ) .EQ. 0)
                      IF (LIZERO .AND. LJZERO) THEN
                        PIFACT = 1.0D0
                        LPIDOP = .TRUE.
                      ELSEIF (LIZERO .OR. LJZERO) THEN
                        PIFACT = PININV
                        LPIDOP = .TRUE.
                      ELSEIF (ILPICI .EQ. ILPJC(JJ)) THEN
                        PIFACT = PININV
                        LPIDOP = .TRUE.
                      ELSE
                        PIFACT = 0.0D0
                        LPIDOP = .FALSE.
                      ENDIF
                    ENDIF
C     
                    IF (LPIDOP) THEN
C     
C---- CALCULATE DISTANCE AND CONNECTING VECTOR BETWEEN II AND JJ
C     CORRECT FOR BOUNDARY CONDITIONS IF NECESSARY
C     
!$OMP parallel do private(M)
                      DO M=1,NCALCD
                        XIJDST(M) = XIAT(ILOC3+M) - XJAT(JLOC3+M)
                      ENDDO
                      RIJ2 = 0.0D0
                      IF (LVAC) THEN
                        DO M=1,NCALCD
                          RIJ2 = RIJ2 + XIJDST(M)**2
                        ENDDO
                        IF (L4D ) THEN
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
                          DSTTMP = BOXOQ -DABS(XIJDST(1))-
     .                         DABS(XIJDST(2))-DABS(XIJDST(3))
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
C     
                      NDXGRP =  NRELKP(ILGRP,JLGRP)
                      QIQJA = CGILA * CGJLA
C     
                      IF (LIPERT .OR. LJPERT) THEN
C     
C---- IF EITHER II OR JJ IS PERTURBED, A SOFT-CORE PERTURBATION
C     CALCULATION IS DONE
C     GET STATE B PARAMETERS AND USE COMBINATION RULE FOR SOFT-CORE
C     ALPHA PARAMETERS
C     
                        JLBTYP = JTYPEB(JJ)
                        INTIJB = MPAC(ILBTYP,JLBTYP)
                        QIQJB = CGILB * CGJLB
C     
                        CB12 = C12(INTIJB)
                        CB6  =  C6(INTIJB)
C     
                        SOFTLJ = ALPHLJ * (1 - ISCLJI * ISCLJJ)
                        SOFTEL = ALPHC * (1 - ISCELI * ISCELJ)
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
                          VRFA4D = -0.5D0*QRFR
C     LENNARD-JONES
                          SFLJ   = SOFTLJ * CA126
                          CRA1   = 1.0D0 / (SFLJ * RL2A + RIJ6)
                          CRA12  = CA12 * CRA1
                          VLJA4D = CRA1 * (CRA12 - CA6)
                          DFA4D  = CRA1 * (CRA1*CRA12 + VLJA4D)
C     DV/DL
                          DLLJA4 = -2.0D0 * SFLJ * DFA4D
                          DLELA4 = -VELA4D * SOFTEL * RSOFT2
                          DLRFA4 = 1.5D0 * QRFR * SOFTEL * RRF2
C     FORCES
                          DFA4D = RSOFT2 * VELA4D + 6.0D0*RIJ4*DFA4D
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
                          DFB4D = RSOFT2 * VELB4D + 6.0D0*RIJ4*DFB4D
C     
C--   COMBINE 4D STATE A AND B
C     
                          VEL4D  = RLB  * VELB4D + RLA  * VELA4D
                          VRF4D  = RLB  * VRFB4D + RLA  * VRFA4D
                          VLJ4D  = RLB  * VLJB4D + RLA  * VLJA4D
                          DF4D4D = RLB  * DFB4D  + RLA  * DFA4D
                          DF4D   = DF4D4D + RLB*DFRFB4 + RLA*DFRFA4
                          DLLJ4D = RDLB * VLJB4D + RLLB * DLLJB4
     .                           + RDLA * VLJA4D + RLLA * DLLJA4
                          DLEL4D = RDLB * VELB4D + RLLB * DLELB4
     .                           + RDLA * VELA4D + RLLA * DLELA4
                          DLRF4D = RDLB * VRFB4D + RLLB * DLRFB4
     .                           + RDLA * VRFA4D + RLLA * DLRFA4
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
                          CRA1   = 1.0D0/(SOFTLJ*CA126*RL2A + RIJ6)
                          VLJA4D = CRA1 * (CA12 * CRA1 - CA6)
C     
C--   4D, STATE B
C     
C     COULOMB
                          SFEL   = SOFTEL * RL2B 
                          VELB4D = QIQJB / DSQRT(SFEL + RIJ2)
C     REACTION FIELD
                          RRF2   = 1.0D0 / (SFEL + RCRF2)
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
C---- SOFT-CORE PERTURBATION, IF 0.0<=RMULOC<1.0 
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
C--   COMBINE 3D STATE A AND B
C     
                          VEL3D  = RLB * VELB3D + RLA * VELA3D
                          VRF3D  = RLB * VRFB3D + RLA * VRFA3D
                          VLJ3D  = RLB * VLJB3D + RLA * VLJA3D
                          DF3D   = RLB * (DFB3D + DFRFB3)
     .                           + RLA * (DFA3D + DFRFA3)
                          DLLJ3D = RDLB* VLJB3D + RLLB * DLLJB3
     .                           + RDLA* VLJA3D + RLLA * DLLJA3
                          DLEL3D = RDLB* VELB3D + RLLB * DLELB3 
     .                           + RDLA* VELA3D + RLLA * DLELA3
                          DLRF3D = RDLB* VRFB3D + RLLB * DLRFB3
     .                           + RDLA* VRFA3D + RLLA * DLRFA3
                          DDM3EL = RLBRMU * DLELB3 + RLARMU * DLELA3
                          DDM3RF = RLBRMU * DLRFB3 + RLARMU * DLRFA3
                          DDM3LJ = RLBRMU * DLLJB3 + RLARMU * DLLJA3
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
C     
                        ELREL(NDXGRP) = ELREL(NDXGRP) + 
     $                       RM3D * VEL3D + RM4D * VEL4D
                        ELRRF(NDXGRP) = ELRRF(NDXGRP) + 
     $                       RM3D * VRF3D + RM4D * VRF4D
                        ELRRC(NDXGRP) = ELRRC(NDXGRP) + 
     $                       RMCST * VRC
                        ELRLJ(NDXGRP) = ELRLJ(NDXGRP) + 
     $                       RM3D * VLJ3D + RM4D * VLJ4D
                        ELGLEL = ELGLEL + 
     $                       RM3D * DLEL3D + RM4D * DLEL4D
                        ELGLRF = ELGLRF + 
     $                       RM3D * DLRF3D + RM4D * DLRF4D
                        ELGLRC = ELGLRC + 
     $                       RMCST * DLRC
                        ELGLLJ = ELGLLJ + 
     $                       RM3D * DLLJ3D + RM4D * DLLJ4D
                        
                        EL34EL = EL34EL + 
     $                       RDM3D*VEL3D + RM3D*DDM3EL +
     $                       RDM4D*VEL4D
                        EL34RF = EL34RF +
     $                       RDM3D*VRF3D + RM3D*DDM3RF +
     $                       RDM4D*VRF4D
                        EL34RC = EL34RC + RMDER * VRC
                        EL34LJ = EL34LJ +
     $                       RDM3D*VLJ3D + RM3D*DDM3LJ + RDM4D*VLJ4D
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
                        VELA4D = QIQJA * RIJINV
C     
                        CRA1   = CA12 * RIJIN6
                        CRA12  = CRA1 - CA6
                        VLJA4D = CRA12 * RIJIN6
                        DFA4D  = RIJIN2 * (VELA4D + 
     $                       6.0D0 *(CRA1 + CRA12)*RIJIN6)
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
                        ELRRF(NDXGRP) = ELRRF(NDXGRP) + 
     $                       QIQJA * RFE * RIJ3D2
                        ELRRC(NDXGRP) = ELRRC(NDXGRP) + 
     $                       QIQJA * RFC
                        ELREL(NDXGRP) = ELREL(NDXGRP) + 
     $                       VELA4D
                        ELRLJ(NDXGRP) = ELRLJ(NDXGRP) + 
     $                       VLJA4D
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
                      FIAT(ILOC3 + 1) = FIAT(ILOC3 + 1) + XH
                      FJAT(JLOC3 + 1) = FJAT(JLOC3 + 1) - XH
                      XH = XIJDST(2) * DF3
                      FIAT(ILOC3 + 2) = FIAT(ILOC3 + 2) + XH
                      FJAT(JLOC3 + 2) = FJAT(JLOC3 + 2) - XH
                      XH = XIJDST(3) * DF3
                      FIAT(ILOC3 + 3) = FIAT(ILOC3 + 3) + XH
                      FJAT(JLOC3 + 3) = FJAT(JLOC3 + 3) - XH
C     
C--   THAT SHOULD WORK FOR 4D (PERT,NOPERT) AND >4D (NOPERT)
C     
                      DO M=4, NCALCD
                        XH = XIJDST(M) * DF4TH
                        FIAT(ILOC3 + M) = FIAT(ILOC3 + M) + XH
                        FJAT(JLOC3 + M) = FJAT(JLOC3 + M) - XH
                      ENDDO
C     
C---- CALCULATE THE VIRIAL IF REQUIRED
C     
                      IF (LDOIT) THEN
                        VIRLR(1) = VIRLR(1) - XIJDST(1)*DF3*
     .                    (XIJDST(1)-XRIAT(ILOC3+1)+XRJAT(JLOC3+1))
                        VIRLR(2) = VIRLR(2) - XIJDST(2)*DF3*
     .                    (XIJDST(2)-XRIAT(ILOC3+2)+XRJAT(JLOC3+2))
                        VIRLR(3) = VIRLR(3) - XIJDST(3)*DF3*
     .                    (XIJDST(3)-XRIAT(ILOC3+3)+XRJAT(JLOC3+3))
                      ENDIF
C     
C---- DONE WITH THE VIRIAL
C     
                    ENDIF
C     
                    JLOC3 = JLOC3 + NDIM
C     
C------END LOOP OVER JJ ATOMS
C     
 400              CONTINUE
C     
                  ILOC3 = ILOC3 + NDIM
C     
C--------END LOOP OVER II ATOMS
C     
 301            CONTINUE
C     
C---- ADD FORCE OF NI TO GLOBAL FORCE ARRAY
C
!$OMP parallel do private(JJ) shared(JCOPTR)
                DO JJ=1,JNUMA3
                  FLR(JCOPTR+JJ) = FLR(JCOPTR+JJ) + FJAT(JJ)
                ENDDO                        
C     
C---- ENDIF CUT-OFF CHECK
C     
              ENDIF
              
              JCGPTR = JCGPTR+NDIM
C     
C---- ENDIF FREE ENERGY PERTURBATION SWITCH
C     
            ENDIF
C     
C----------END LOOP OVER NJ SECONDARY CHARGE GROUPS
C     
 201      CONTINUE
C     
!$OMP parallel do private(II) shared(ICOPTR)
          DO II=1,INUMA3
            FLR(ICOPTR+II) = FLR(ICOPTR+II)+FIAT(II)
          ENDDO
C     
          INB(NICOR) = NJPTR - JSAVE                  
          NICOR  = NICOR + 1
          IFTOPA = ILTOPA + 1
          IFCORA = IFTOPA + IATOFF
          ICOPTR = ICOPTR+INUMA3
          ICGPTR = ICGPTR+NDIM
C     
C---- ADD FORCE OF NI TO GLOBAL FORCE ARRAY
C     
 100    CONTINUE                
        IATOFF = IATOFF+ NRP
 17   CONTINUE
C     
C     
C------------END LOOP OVER PROTEIN MOLECULES AND PRIMARY (PROTEIN ONLY)
C     CHARGE GROUPS      
C     
C--   IF WE ARE DOING FREE ENERGY PERTURBATION (INTERACTION/PAIR LIST
C     LIMITED TO PERTURBED ATOMS), WE ARE DONE
C     
      IF (LFREE) THEN
        PRINT 500,NJPTR
        RETURN
      ENDIF
C     
C------------NOW LOOP OVER PRIMARY (SOLVENT ONLY) CHARGE GROUPS
C     SECONDARY CHARGE GROUP WILL ALSO BE SOLVENT I.E.
C     NO PERTURBATION, NO INTERNAL CHARGE GROUP INTERACTION,
C     NO EXCLUDED AND THIRD NEIGHBOURS
C     
      INUMA3 = NDIM*NRAM
      IFCORA = NRPT+1
      ICOPTR = NRPT*NDIM
      ICGPTR = NRAGP*NDIM
C     
      DO 700 NICOR=NRAGP+1,NRAGT-1
        IAUX = IFCORA-1
        DO II=1,NRAM
          IGRP(II)   = ICHAR(IAGRP(IAUX+II))
        ENDDO
C     
!$OMP parallel do private(II) shared(ICOPTR)
        DO II=1,INUMA3
          XIAT(II)  = XCOORD(ICOPTR+II)
          XRIAT(II) = XR(ICOPTR+II)
          FIAT(II)= 0.0D0
        ENDDO
C     
!$OMP parallel do private(M) shared(ICGPTR)
        DO M=1,NDIM
          XCGI(M) = XCG(ICGPTR+M)
        ENDDO
C     
        JSAVE  = NJPTR
        JCGPTR = ICGPTR+NDIM
        JCOPTR = ICOPTR+INUMA3
        DO 800 NJCOR = NICOR+1,NRAGT
C     
C--------LOOP OVER SECONDARY CHARGE GROUPS (SOLVENT)
C     
C---- CALCULATE DISTANCE AND CONNECTING VECTOR BETWEEN XCGI AND XCG(NJCOR)
C     CORRECT FOR BOUNDARY CONDITIONS IF NECESSARY
C     
          DO M=1,NCALCD
            XIJ(M) = XCGI(M) - XCG(JCGPTR + M)
          ENDDO
          RIJCG2 = 0.0D0
          IF (LVAC) THEN
              RIJCG2 = XIJ(1)*XIJ(1)+XIJ(2)*XIJ(2)+XIJ(3)*XIJ(3)
          ELSE
            DO M=1,3
              DSTTMP = XIJ(M)
              IF (DSTTMP .GE. BOXH(M)) THEN
                DSTTMP = DSTTMP - BOX(M)
              ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
                DSTTMP = DSTTMP + BOX(M)
              ENDIF
              XIJ(M) = DSTTMP
              RIJCG2 = RIJCG2 + DSTTMP*DSTTMP
            ENDDO
C     
            IF (LOCTO) THEN
              DSTTMP = BOXOQ -DABS(XIJ(1))-
     .           DABS(XIJ(2))-DABS(XIJ(3))
              IF (DSTTMP .LT. 0.0D0) THEN
                RIJCG2 = RIJCG2 + DSTTMP*BOX(1)
                XIJ(1) = XIJ(1) - SIGN(BOXOH,XIJ(1))
                XIJ(2) = XIJ(2) - SIGN(BOXOH,XIJ(2))
                XIJ(3) = XIJ(3) - SIGN(BOXOH,XIJ(3))
              ENDIF
            ELSEIF (LMONO .AND. LDOTRA) THEN
              RIJCG2 = RIJCG2 + COSB2*XIJ(1)*XIJ(3)
            ENDIF
          ENDIF
C     
C---- DONE DISTANCE CALCULATION
C     CHECK DISTANCE
C     
          IF (RIJCG2 .LE. RCUTP2) THEN
C     
C---- INSIDE SHORT-RANGE CUT-OFF, STICK INTO PAIRLIST
C     
            NJPTR = NJPTR + 1
            IF (NJPTR .GT. MAXJNB) THEN
              PRINT *,PRGSTR,' PAIRLIST OVERFLOW:'
              PRINT *,'MAXIMUM SIZE: MAXJNB =',MAXJNB
              PRINT *,'              NJPTR  =',NJPTR ! cmb
              PRINT *,NICOR,NJCOR
              CALL STOPGM('NBPML','FATAL ERROR')
            ENDIF
            JNB(NJPTR) = NJCOR
          ELSEIF (RIJCG2.LE.RCUTL2 ) THEN
C     
C---- WITHIN THE TWIN RANGE, CALCULATE ELECTROSTATIC INTERACTION
C     IF THE SOLVENT IS CHARGED
C     
            JFCORA = NRPT + NRAM*(NJCOR-NRAGP-1) + 1
            IAUX = JFCORA-1
            DO JJ=1,NRAM
              JGRP(JJ)   = ICHAR(IAGRP(IAUX+JJ))
            ENDDO
C     
C---- PER SECONDARY CHARGE GROUP USE LOCAL ARRAYS
C     (TOPOLOGY INFORMATION) - SOLVENT
C     
!$OMP parallel do private(JJ) shared(JCOPTR)
            DO JJ=1,INUMA3
              XJAT(JJ)  = XCOORD(JCOPTR+JJ)
              XRJAT(JJ) = XR(JCOPTR+JJ)
              FJAT(JJ)= 0.0D0
            ENDDO
C     
C---- DONE WITH SETTING UP THE LOCAL ARRAYS
C     
C--------LOOP OVER ATOMS II OF CHARGE GROUPS NI
C     
            ILOC3 = 0
            DO 850 II=1,NRAM
              CGILA = CGS(II)
              ILATYP = IACS(II)
              ILGRP = IGRP(II)
C     
C------LOOP OVER ATOMS JJ OF CHARGE GROUP NJ
C     
              JLOC3 = 0
              DO 900 JJ=1,NRAM
                CGJLA = CGS(JJ)
                JLATYP = IACS(JJ)
                JLGRP = JGRP(JJ)                          
C     
C---- CALCULATE DISTANCE AND CONNECTING VECTOR BETWEEN II AND JJ
C     CORRECT FOR BOUNDARY CONDITIONS IF NECESSARY
C     
                DO M=1,NCALCD
                  XIJDST(M) = XIAT(ILOC3+M) - XJAT(JLOC3+M)
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
     .               DABS(XIJDST(2))-DABS(XIJDST(3))
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
C     

                NDXGRP =  NRELKP(ILGRP,JLGRP)
C     
                INTIJA = MPAC(ILATYP,JLATYP)
                CA12   =  C12(INTIJA)
                CA6    =  C6(INTIJA)
C     
                RIJIN2 = 1.0D0 / RIJ2
                RIJINV = DSQRT(RIJIN2)
                RIJIN6 = RIJIN2 * RIJIN2 * RIJIN2

                QIQJA  = CGILA * CGJLA
                VELA4D = QIQJA * RIJINV
C     
                CRA1   = CA12 * RIJIN6
                CRA12  = CRA1 - CA6
                VLJA4D = CRA12 * RIJIN6
                DFA4D  =  RIJIN2 * (VELA4D + 
     $               6.0D0 * (CRA1 + CRA12) * RIJIN6)
C     
                ELRRF(NDXGRP) = ELRRF(NDXGRP) + 
     $               QIQJA * RFE * RIJ3D2
                ELRRC(NDXGRP) = ELRRC(NDXGRP) + 
     $               QIQJA * RFC
                ELREL(NDXGRP) = ELREL(NDXGRP) + VELA4D
                ELRLJ(NDXGRP) = ELRLJ(NDXGRP) + 
     $               VLJA4D
C     
                DF3   = DFA4D + RFF * QIQJA
                DF4TH = DFA4D
C     
C     CALCULATE AND STORE THE FORCE INTO THE LOCAL FORCE ARRAY
C     
                XH = XIJDST(1) * DF3
                FIAT(ILOC3 + 1) = FIAT(ILOC3 + 1) + XH
                FJAT(JLOC3 + 1) = FJAT(JLOC3 + 1) - XH
                XH = XIJDST(2) * DF3
                FIAT(ILOC3 + 2) = FIAT(ILOC3 + 2) + XH
                FJAT(JLOC3 + 2) = FJAT(JLOC3 + 2) - XH
                XH = XIJDST(3) * DF3
                FIAT(ILOC3 + 3) = FIAT(ILOC3 + 3) + XH
                FJAT(JLOC3 + 3) = FJAT(JLOC3 + 3) - XH
C     
C--   THAT SHOULD WORK FOR 4D (PERT,NOPERT) AND >4D (NOPERT)
C     
                DO M=4, NCALCD
                  XH = XIJDST(M) * DF4TH
                  FIAT(ILOC3 + M) = FIAT(ILOC3 + M) + XH
                  FJAT(JLOC3 + M) = FJAT(JLOC3 + M) - XH
                ENDDO
C     
C---- CALCULATE THE VIRIAL IF REQUIRED
C     
                IF (LDOVIR) THEN
                  VIRLR(1) = VIRLR(1) - XIJDST(1)*DF3*
     .               (XIJDST(1)-XRIAT(ILOC3+1)+XRJAT(JLOC3+1))
                  VIRLR(2) = VIRLR(2) - XIJDST(2)*DF3*
     .               (XIJDST(2)-XRIAT(ILOC3+2)+XRJAT(JLOC3+2))
                  VIRLR(3) = VIRLR(3) - XIJDST(3)*DF3*
     .               (XIJDST(3)-XRIAT(ILOC3+3)+XRJAT(JLOC3+3))
                ENDIF
C     
C---- DONE WITH THE VIRIAL
C     
                JLOC3 = JLOC3 + NDIM
C     
C------END LOOP OVER JJ ATOMS
C     
 900          CONTINUE
C     
              ILOC3 = ILOC3 + NDIM
C     
C--------END LOOP OVER II ATOMS
C     
 850        CONTINUE
C     
C---- ADD FORCE OF NI TO GLOBAL FORCE ARRAY
C     
!$OMP parallel do private(JJ) shared(JCOPTR)
            DO JJ=1,INUMA3
              FLR(JCOPTR+JJ) = FLR(JCOPTR+JJ) + FJAT(JJ)
            ENDDO
C     
          ENDIF
C     
          JCOPTR = JCOPTR+INUMA3
          JCGPTR = JCGPTR+NDIM
 800    CONTINUE
C     
C---- ADD FORCE OF NI TO GLOBAL FORCE ARRAY
C     
!$OMP parallel do private(II) shared(ICOPTR)
        DO II=1,INUMA3
          FLR(ICOPTR+II) = FLR(ICOPTR+II)+FIAT(II)
        ENDDO
C
        INB(NICOR) = NJPTR - JSAVE
        IFCORA = IFCORA+NRAM
        ICOPTR = ICOPTR+INUMA3
        ICGPTR = ICGPTR+NDIM

C     
C------------END LOOP OVER SOLVENT MOLECULES
C     
 700  CONTINUE

      INB(NRAGT) = 0

      PRINT 500,NJPTR
C     
C------------CORRECT VIRIAL IF NECESSARY
C     
      IF (LDOVIR) THEN
!$OMP parallel do private(M)
        DO M=1,NDIM
          VIRLR(M) = VIRLR(M)*0.5D0
        ENDDO
      ENDIF
C     
C----------RETURN DERIVATIVE VARIABLES ACCORDING TO REQUEST
C     
      IF (.NOT.L3D4D) THEN
C     
C--   IN THIS CASE THESE WOULD BE *WRONG*
C     
        EL34EL = 0.0D0
        EL34RF = 0.0D0
        EL34RC = 0.0D0
        EL34LJ = 0.0D0
      ENDIF
      IF (.NOT.LDOPER) THEN
C     
C--   IN THIS CASE THESE WOULD BE OK, BUT SINCE BONDED ROUTINES
C     WONT RETURN THEM, RETURNING THEM MIGHT GENERATE CONFUSION... 
C     
        ELGLEL = 0.0D0
        ELGLRF = 0.0D0
        ELGLRC = 0.0D0
        ELGLLJ = 0.0D0
      ENDIF
C     
C------------END NBWITH
C     
 500  FORMAT (1X,'NBPML:',I10,' ELEMENTS IN THE PAIRLIST')
C     ==--------------------------------------------------------------==
      RETURN
      END

#endif
#ifdef __VECTOR
C     Id: nbpmlv.f,v 1.18 1996/07/23 11:45:37 ili Exp 
C     

COMMSUBR NBPML CCCCCC I.G. TIRONI, ZUERICH, JUNE 1996 CCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     SUBROUTINE NBPML(NATTOT,NPM,NSM,NRAGT,NCALCD,
C     $     XCOORD,XCG,XR,FLR,
C     $     NSPM,NSP,VIRLR,
C     $     INB,NSZPL,JNB,RCUTP,RCUTL,
C     $     MAXNRE,NRELKP,IAGRP,
C     $     RCRF,EPSRF,APPAK,
C     $     ELGLEL,ELGLRF,ELGLRC,ELGLLJ,
C     $     EL34EL,EL34RF,EL34RC,EL34LJ,
C     $     NUNRE2,ELREL,ELRRF,ELRRC,ELRLJ,
C     $     LDOPER,L3D4D,RLAM,RMU,ALPHLJ,ALPHC,NLAM,MMU,
C     $     LPION,LEVERY,LFREE)
C     
C     
C     This is a special version suited for vectorprocessor machines
C     like a CRAY. It is derived from the old cray vector routine.
C     Path-integral and free energy problems do not run as efficient
c     as a normal MD calculation.
C     Ilario Tironi June 1996.
C     
C     NBPML makes a list of non-bonded pairs using a cut-off criterion
C     for the centres of the geometry of solute charge groups
C     and the first atoms of solvent molecules, by scanning all
C     possible pairs (using 3D distances only). The cut-off radius
C     for mating the list is L<RCUTP>.
C     The pair list contains ordered sequence numbers of solute charge atom
C     groups and solvent molecules (NI<NJ).
C     Periodic images can be included in the list. In this case, the
C     atoms of a charge group (solute) must lie within BOX/2 or
C     BOX*SQRT(3)/4 of each other.
C     
C     If L<RCUTL> > L<RCUTP>, a potential is evaluated for the twin-range
C     region in subroutine L<NBWITH> during the pair list construction. 
C     Otherwise the pair list is constructed without long range forces 
C     in subroutine L<NBNONE>. The term long-range in contrast to GROMOS87
C     includes Lennard-Jones interactions.
C     
C     The potential in L<NBWITH> reads:
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
C     The cut-off criterion is always applied in 3D.
C     
C     Periodic boundary conditions can be applied, depending on L<NTB>,
C     where the periodic box can be a truncated octahedron, rectangular
C     or monoclinic. 
C     
C     The periodic box can be a truncated octahedron, rectangular or
C     monoclinic, depending on L<NTB>.
C     If the long range virial is to be calculated
C     (ABS(L<NTB> .EQ. L<NTBVIR>), it must be specified whether the
C     solute consists of separate submolecules or not.
C     This can be done by setting the NSPM > 1 and setting the corresponding
C     values into the NSP array. The virial is only calculated
C     in 3D.
C     
C     
C     NATTOT   total number of atoms
C     NPM      number of (identical) solute molecules
C     NRAGT    total number of charge groups
C     NCALCD   the number of dimensions to calculate distance and
C     forces in. Must be 3 or 4.
C     Note that the 3D distance between charge groups
C     I and J is _always_ used for deciding whether to put
C     the pair into the list or not, independently of NCALCD.
C     The distance used in the force calculation, however,
C     does depend on NCALCD.
C     XCOORD(NDIM*NATTOT)
C     the atom cartesian coordinates
C     XCG(NDIM*NRAGT)
C     the geometric centre of solute charge groups and the
C     first atoms of solvent molecules.
C     XR(NDIM*NATTOT)
C     atom cartesian coordinates relative to submolecular
C     centres of mass. This array is only used if
C     the virial is calculated, i.e. ABS(NTB) = L<NTBVIR>.
C     It is still accessed (read only) in the opposite case 
C     and thus should have the correct dimension (e.g. XCOORD)
C     FLR(NDIM*NATTOT)
C     delivered with the long range coulomb forces if
C     RCUTL > RCUTP.
C     NSPM     If the virial is to calculated, then L<NSPM> and L<NSP> are
C     accessed in order to calculate the molecular virial.
C     and thus the dummy array used in the call should have 
C     the correct dimension (e.g. XCOORD).
C     NSP(NSPM+1)
C     atom sequence numbers of the last atoms of the
C     submolecules.
C     only elements 1..NSPM are considered to be filled
C     with data. We need NSPM+1 elements in order to avoid
C     a double IF statement when checking for submolecules
C     
C     VIRLR(NDIM)
C     delivered with the X-, Y-, and Z- components of the
C     long range virial if L<NTB> = L<NTBVIR>.
C     INB(NRAGT)
C     delivered with the pointer list;
C     INB(I) specifies the number of charge groups or
C     solvent molecules J with J>I forming a pair with charge
C     group or solvent molecule I
C     NSZPL    delivered with the number of elements in JNB
C     JNB(MAXJNB)
C     delivered with the pointer list;
C     it specifies the solute charge groups or solvent
C     molecules J forming a pair with solute charge group
C     or solvent molecule I in ascending order; for each I the
C     pointer in JNB is incremented by INB(I)
C     MAXNRE   maximum number of energy groups
C     NRELKP(MAXNRE,MAXNRE)   
C     index of pairs of energy groups in energy arrays
C     IAGRP(MAXNAT)
C     for each atom, the energy type code its
C     contribution is added to.
C     RCRF     radius of the reaction field boundary (in principle
C     RCUTL) and 
C     EPSRF    dielectric permittivity of the reaction field 
C     continuum outside the boundary
C     APPAK    inverse Debye screening length
C     NUNRE2   number of energy froups
C     ELREL(MXNRE2)
C     delivered with long range coulomb energies.
C     ELRRF(MXNRE2)
C     delivered with the long range distance dependent
C     Poisson-Boltzmann reaction field energies
C     ELRRC(MXNRE2)
C     delivered with the long range distance independent
C     Poisson-Boltzmann reaction field energies
C     ELGLEL   delivered with dV/d(L<LAMBDA>) at fix L<MU> 
C     for the long range Coulomb term if LDOPER = .TRUE. 
C     ELGLRF   delivered with with dV/d(L<LAMBDA>) at fix L<MU> 
C     for the distance dependent long range 
C     Poisson-Boltzmann term if LDOPER = .TRUE. 
C     ELGLRC   delivered with with dV/d(L<LAMBDA>) at fix L<MU> 
C     for the distance independent long range 
C     Poisson-Boltzmann term if LDOPER = .TRUE. 
C     EL34EL   delivered with the 3D to 4D long range free energy contribution
C     for the Coulomb term if LDOPER = .TRUE. .AND. NCALCD = 4
C     EL34RF   delivered with the 3D to 4D long range free energy contribution
C     for the long range distance dependent 
C     Poisson-Boltzmann reaction field
C     term if LDOPER = .TRUE. .AND. NCALCD = 4
C     EL34RC   delivered with the 3D to 4D long range free energy contribution
C     for the long range distance independent 
C     Poisson-Boltzmann reaction field
C     term if LDOPER = .TRUE. .AND. NCALCD = 4
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
C     .FALSE.  a number of quantities depending on L<NPM>, L<NRP>,
C     L<NSM>, L<NRAM>, L<NCAG>, L<NTB>, L<BOX>, L<BETA>,
C     L<RCUTP> and L<RCUTL> are only calculated at the 
C     first call
C     .TRUE.   they are calculated at every subr. call.
C     LFREE    
C     .TRUE.   the pair list and interaction is calculated solely 
C     for perturbed atoms. This is used in free energy 
C     extrapolation.
C     
C     NOT passed as arguments, but accessed in common blocks are
C     the following:
C     
C     RCUTP    cut-off radius for the pair list; since no check on
C     neighbour exclusions is performed, and since only charge
C     group pairs I<J are scanned, L<RCUTP> should be chosen
C     large enough when L<RCUTL> > L<RCUTP>.
C     RCUTL    cut-off radius for the long range coulomb force;
C     The following conditions must hold, depending on the
C     boundary conditions specified.
C     If L<NTB> <0, L<RCUTP>, L<RCUTL> < BOX(M)*SQRT(3)/4
C     If L<NTB> >0, L<RCUTP>, L<RCUTL> < BOX(M)/2
C     If L<RCUTL> is larger, only the nearest image is calculated.
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
C     
COMMEND
C     ==--------------------------------------------------------------==      
      SUBROUTINE NBPML(NATTOT,NPM,NSM,NRAGT,NCALCD,
     $     XCOORD,XCG,XR,FLR,
     $     NSPM,NSP,VIRLR,
     $     INB,NSZPL,JNB,RCUTP,RCUTL,
     $     MAXNRE,NRELKP,IAGRP,
     $     RCRF,EPSRF,APPAK,
     $     ELGLEL,ELGLRF,ELGLRC,ELGLLJ,
     $     EL34EL,EL34RF,EL34RC,EL34LJ,
     $     NUNRE2,ELREL,ELRRF,ELRRC,ELRLJ,
     $     LDOPER,L3D4D,RLAM,RMU,ALPHLJ,ALPHC,NLAM,MMU,
     $     LPION,LEVERY,LFREE)
C     ==--------------------------------------------------------------==     
C     INCLUDES
C     
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'pertsz.h'
      INCLUDE 'pertar.h'
      INCLUDE 'vector.h'
      INCLUDE 'nbpml.h'
C     
C     ARGUMENTS
C     
      LOGICAL LDOPER,L3D4D,LPION,LEVERY,LFREE
      INTEGER NATTOT,NRAGT,NPM,NSM,NCALCD,NSZPL
      INTEGER INB(NRAGT),JNB(MAXJNB)
      INTEGER NSPM,NSP(NSPM+1)
      INTEGER MAXNRE,NRELKP(MAXNRE,MAXNRE)
      INTEGER NLAM,MMU
      INTEGER NUNRE2
      real*8 RCUTP,RCUTL,RCRF,EPSRF,APPAK
      real*8 XCOORD(NDIM*NATTOT),XCG(NDIM*NRAGT)
      real*8 XR(NDIM*NATTOT),FLR(NDIM*NATTOT)
      real*8 VIRLR(NDIM)
      real*8 ELREL(NUNRE2),ELGLEL,EL34EL
      real*8 ELRRF(NUNRE2),ELGLRF,EL34RF
      real*8 ELRRC(NUNRE2),ELGLRC,EL34RC
      real*8 ELRLJ(NUNRE2),ELGLLJ,EL34LJ
      real*8 RLAM,RMU,ALPHLJ,ALPHC
      CHARACTER IAGRP(NATTOT)
C     
C     LOCAL PARAMS
C     
      real*8 ENEPS
      PARAMETER (ENEPS = 1.D-20)
      real*8 ZEROC6, ZEROAL
      PARAMETER (ZEROC6 = 1.0D-20, ZEROAL = 1.0D-20)
C     
C     SAVED VARIABLES      
C     
      LOGICAL LFIRST,LDOLNG, LERROR
      SAVE LFIRST,LDOLNG
C     
C     LOCAL VARIABLES
C     
      INTEGER NATTO3,I,II
      INTEGER ITMPA,ITMPB,JTYPE
      INTEGER INTIJA,INTIJB
      real*8 CONV,BETAR,PAPPAK
      real*8 XTEMP(MAXDIM)
      CHARACTER PRGSTR*(5)
      
      INTEGER NI, NNOLD, NS, NN, NAG,
     $     LEN, JTOP, JRC, J1X, ICLO,
     $     JRCRAT, JRCFOR, NFORCE, JRCDUP,
     $     IPQ, JRCOFF, JRC1, JRC2,
     $     NVIR, NSPI, J, KFORCE, M, N, IPC

      INTEGER J2X, NSRC3, NDST3, NCG, NPP, ILACG, IFICG, NUMCG
C     
C     DATA
C     
      DATA LFIRST /.TRUE./
      DATA LNEWP /.FALSE./
      DATA PRGSTR  /'NBPML'/
C     
C------------BEGIN NBPML

      IF (LFREE) THEN
         WRITE(6,*) 'FREE ENERGY EXTRAPOLATION (LFREE .EQ. TRUE)',
     $        ' NOT IMPLEMENTED IN VECTOR VERSION OF NBPML'
         WRITE(6,*) 'USE THE STANDARD VERSION OF NBPML'
         CALL STOPGM('NBPML','FATAL ERROR')
      ENDIF
C     
C------------DO FIRST CALL SETUP (SWITCHES,"MONOCLINIC" VARIABLES
C     REACTION FIELD VARIABLES, SIMPLE CONSISTENCY CHECKS)
C     
      IF (LFIRST .OR. LEVERY) THEN

         LMONO = (NTB .GT. 0)
         LOCTO = (NTB. LT. 0)
         LVAC  = (NTB .EQ. 0)
         LDOVIR= (ABS(NTB) .EQ. NTBVIR)
C     
         LDOTRA = .FALSE.
         IF (LMONO)THEN
            CONV = 45.0D0/DATAN(1.0D0)
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


      IF (LFIRST) THEN
         LFIRST = .FALSE.
         L4D = (NCALCD .EQ. 4)
         LPIDOP = .TRUE.
         IF (NPID .NE. 0) THEN
            PININV = 1.0D0 / NPID
         ENDIF
         LPERTL = (LDOPER.OR.L3D4D)
         LDOLNG = (RCUTL .GT. RCUTP)       

         IF (MAXINB .LT. MAXP) THEN
            WRITE(6,*) 'WARNING IN NBPML'
            WRITE(6,*) 'MAXINB IS SMALLER THAN MAXP!'
            WRITE(6,*) 'VECTOR LENGTH COULD BE TO SMALL',
     .           ' INCREASE MAXINB OR SWAP MAXINB WITH MAXP IN',
     .           ' THE ARRAYSIZES SECTION OF QMMM'
cmb - see defaults in mm_read_qmmm_input.F - former vector.h
         ENDIF

         IF (MAXGP .LT. (NPM * NRP)) THEN
            WRITE(6,*) 'WARNING IN NBPML'
            WRITE(6,*) 'MAXGP IS SMALLER THAN NUMBER OF SOLUTE ATOMS!'
            WRITE(6,*) 'VECTOR LENGTH COULD BE TO SMALL INCREASE MAXGP',
     .           ' IN THE ARRAYSIZES SECTION OF QMMM'
cmb - see defaults in mm_read_qmmm_input.F - former vector.h
         ENDIF
            
C     
C---- PERFORM A FEW FIRST CALL CHECKS
C     
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
     $                       C12(INTIJA) .GT. ZEROC6) THEN
                           PRINT 900,PRGSTR,II
                           PRINT 910,'A',JTYPE
                           CALL STOPGM('NBPML','FATAL ERROR')
                        ENDIF
                        IF (C6(INTIJB) .LT. ZEROC6.AND.
     $                       C12(INTIJB) .GT. ZEROC6) THEN
                           PRINT 900,PRGSTR,II
                           PRINT 910,'B',JTYPE
                           CALL STOPGM('NBPML','FATAL ERROR')
                        ENDIF
                        IF (CS6(INTIJA) .LT. ZEROC6.AND.
     $                       CS12(INTIJA) .GT. ZEROC6) THEN
                           PRINT 900,PRGSTR,II
                           PRINT 915,'A',JTYPE
                           CALL STOPGM('NBPML','FATAL ERROR')
                        ENDIF
                        IF (CS6(INTIJB) .LT. ZEROC6.AND.
     $                       CS12(INTIJB) .GT. ZEROC6) THEN
                           PRINT 900,PRGSTR,II
                           PRINT 915,'B',JTYPE
                           CALL STOPGM('NBPML','FATAL ERROR')
                        ENDIF
 28                  CONTINUE
                  ENDIF
               ENDIF
 26         CONTINUE
         ENDIF

C     
C--   CALCULATING IN 4D WITH 3D COORDINATES IS WRONG
C     
         IF ( NCALCD.GT.NDIM ) THEN
            PRINT *,PRGSTR,' ERROR: NCALCD.GT.NDIM'
            CALL STOPGM('NBPML','FATAL ERROR')
         ENDIF
C     
C--   CALCULATING IN 3D DOESN'T ALLOW THE CALCULATION OF MU DERIVATIVE
C     
         IF (.NOT. L4D .AND. L3D4D) THEN
            PRINT *,PRGSTR,': NCALCD.EQ.3.AND.',
     $           '(NTG.EQ.NTGMU.OR.NTG.EQ.NTGBOT) NOT ALLOWED'
         ENDIF
C     
C--   CALCULATING IN 3D WITH RMU.NE.0 SHOULDN'T BE DONE
C     
         IF ( .NOT.L4D.AND.RMU.NE.0.0D0 ) THEN
            PRINT *,PRGSTR,': WARNING, NCALCD .EQ. 3 .AND. RMU .NE. 0.0'
            PRINT *,'         LOCAL RMU WILL BE SET TO ZERO'
            RMULOC = 0.0D0
         ELSE
            RMULOC = RMU
         ENDIF
C     
C---- DONE FIRST CALL CHECKS
C     NOW REACTION FIELD VARIABLES
C     
         IF (EPSRF .LE. 0.0D0) THEN
            C1 = -1.0D0
            PRINT *,PRGSTR,': WARNING, EPSRF.LE.0 , PERMITTIVITY SET ',
     $           'TO INFINITY'
         ELSE
            PAPPAK=APPAK*RCRF
            C1 = ((2.0D0-2.0D0*EPSRF)*(1.0D0+PAPPAK) 
     $           - EPSRF* PAPPAK**2)
     $          /((1.0D0+2.0D0*EPSRF)*(1.0D0+PAPPAK) 
     $           + EPSRF* PAPPAK**2)
         ENDIF
         RFF = C1 / RCRF**3
         RFE = -0.5D0*RFF
         RFC = -(1.0D0-0.5D0*C1)/RCRF
         RCRF2 = RCRF*RCRF

C      CHECK IF THERE IS ENOUGH SPACE FOR ENERGY GROUPS
         IF (NUNRE2 .GT. IEPAIR) THEN
            WRITE(6,*) 'NBPML: ERROR'
            WRITE(6,*) 'NUMBER OF ENERGY GROUPS BIGGER THAN'
            WRITE(6,*) 'WORKSPACE IN NBPML. INCREASE MAXGRP IN',
     .           ' THE ARRAYSIZES SECTION OF QMMM'
cmb - see defaults in mm_read_qmmm_input.F - former vector.h
            CALL STOPGM('NBPML','FATAL ERROR')
         ENDIF

C********************************************************************************
C***************************** FILL THE VECTOR ARRAYS ***************************
C********************************************************************************

C*****DETERMINE FOR EACH ATOM THE INTEGER ATOM CODE.
C     
         JTOP = 1
         DO NN = 1, NPM
            DO JRC = 1, NRP
               IACT(JTOP+JRC-1) = IAC(JRC)
            ENDDO
            JTOP = JTOP + NRP
         ENDDO
C     
         IF (NSM .GT. 0) THEN
            LEN = (NSM-1)*NRAM + 1
            DO NN = 1, NRAM
               DO JRC = 1, NSM
                  IACT((JRC-1)*NRAM+JTOP) = IACS(NN)
               ENDDO
              JTOP = JTOP + 1
            ENDDO
         ENDIF
         
C*****REARRANGE THE DATA STRUCTURE OF THE CHARGES
C     
         JTOP = 1
         DO NN = 1, NPM
            DO JRC = 1, NRP
               CGAT(JTOP+JRC-1) = CG(JRC)
            ENDDO
            JTOP = JTOP + NRP
         ENDDO
C     
         IF (NSM .GT. 0) THEN
            LEN = (NSM-1)*NRAM + 1
            DO NN = 1, NRAM
               DO JRC = 1, NSM
                  CGAT((JRC-1)*NRAM+JTOP) = CGS(NN)
               ENDDO
               JTOP = JTOP + 1
            ENDDO
         ENDIF
         
C***** STORE FOR EACH ATOM THE PATH INTEGRAL DATA
         IF (LPION) THEN
!$OMP parallel do private(I)
            DO I = 1, NRP
               JLIST(I) = IPIC(I)
            ENDDO
            
            ICLO = 0
            DO IPQ = 1, NPM
               DO I = 1, NRP
                  IPICT(ICLO + I) = JLIST(I)
               ENDDO
               ICLO = ICLO + NRP
            ENDDO

            ICLO = ICLO + 1
            IF (NSM .GT. 0) THEN
               LEN = (NSM - 1) * NRAM + 1
               DO I = 1, NRAM
                  DO IPQ = 1, NSM
                     IPICT((IPQ - 1) * NRAM + ICLO) = 0
                  ENDDO
                  ICLO = ICLO + 1
               ENDDO
            ENDIF
         ENDIF
                  
C*****REARRANGE THE DATA STRUCTURE OF THE PERTURBED ATOMS (CHARGES ETC.)
C*****USE LBVIR AND JLIST AS SCRATCH ARRAYS.
C     
         IF (LPERTL) THEN
!$OMP parallel do private(IPQ)
            DO IPQ = 1, NRP
               IF (IPERT(IPQ) .EQ. NOPERT) THEN
                  LBVIR(IPQ) = .TRUE.
                  JLIST(IPQ) = 1
                  LLPERT(IPQ) = .FALSE.
               ELSE
                  LBVIR(IPQ) = .FALSE.
                  JLIST(IPQ) = IPERT(IPQ)
                  LLPERT(IPQ) = .TRUE.
               ENDIF
               IACBT(IPQ) = IACB(JLIST(IPQ))
               CGBT(IPQ) = CGB(JLIST(IPQ))
               ISCCT(IPQ) = 1 - ISCC(JLIST(IPQ))
               ISCLJT(IPQ) = 1 - ISCLJ(JLIST(IPQ))
               IF (LBVIR(IPQ)) THEN
                  IACBT(IPQ) = IACT(IPQ)
                  CGBT(IPQ) = CGAT(IPQ)
                  ISCCT(IPQ) = 1
                  ISCLJT(IPQ) = 1
               ENDIF
            ENDDO   
            
            ICLO = NRP + 1
            DO NN = 2, NPM
!$OMP parallel do private(J1X) shared(ICLO)
               DO J1X = 1, NRP
                  IACBT(ICLO+J1X-1) = IACBT(J1X)
                  CGBT(ICLO+J1X-1) = CGBT(J1X)
                  ISCCT(ICLO+J1X-1) = ISCCT(J1X)
                  ISCLJT(ICLO+J1X-1) = ISCLJT(J1X)
                  LLPERT(ICLO+J1X-1) = LLPERT(J1X)
               ENDDO
               ICLO = ICLO + NRP
            ENDDO
            LEN = NSM*NRAM
            IF (LEN .GT. 0) THEN
!$OMP parallel do private(J1X) shared(ICLO)
#ifdef _vpp_
!OCL NOALIAS
#endif
               DO J1X = 1, LEN
                  IACBT(ICLO+J1X-1) = IACT(ICLO+J1X-1)
                  CGBT(J1X-1+ICLO) = CGAT(J1X-1+ICLO)
                  ISCCT(J1X-1+ICLO) = 1
                  ISCLJT(J1X-1+ICLO) = 1
                  LLPERT(J1X-1+ICLO) = .FALSE.
               ENDDO
            ENDIF
         ENDIF

C     
C--   SET UP THE FORCE CONSTANT TEMPROARIES - CHECK THE SPACE
         KFORCE = (NRATT*(NRATT+1))/2
         NFORCE = 2*KFORCE
C--   CHECK THAT THERE IS ENOUGH SPACE FOR THE CHARGE GROUP TEMPORARIES
         IF (NFORCE .GT. MFORCE) THEN
            PRINT 9140, NFORCE, MFORCE
 9140       FORMAT(/' THE SPACE NEEDED FOR FORCE CONSTANTS (',I5,
     $           ') IS GREATER',
     $           ' THAN THE WORK SPACE PROVIDED (',I5,') IN NONBML'/)
            LERROR = .TRUE.
C     GO TO 5750
            CALL STOPGM('NBPML','FATAL ERROR')
         ENDIF
C--   STORE THE FORCE CONSTANTS
!$OMP parallel do private(JRC) shared(KFORCE)
#ifdef _vpp_
!OCL NOALIAS
#endif
         DO JRC = 1, KFORCE
            FORC6(JRC) = C6(JRC)
            FORC6(JRC+KFORCE) = CS6(JRC)
            FORC12(JRC) = C12(JRC)
            FORC12(JRC+KFORCE) = CS12(JRC)
         ENDDO

C     BUILD UP THE SOFTCORE C12 / C6 COEFFICIENTS
         IF (LPERTL) THEN
!$OMP parallel do private(IPQ)
            DO IPQ = 1, NFORCE
               IF (FORC6(IPQ) .LT. ZEROC6) THEN
                  FC12C6(IPQ) = 0.0D0
               ELSE
                  FC12C6(IPQ) = FORC12(IPQ) / FORC6(IPQ)
               ENDIF
            ENDDO
         ENDIF

C     
C--   COMPUTE THE SPACE NEEDED FOR DUPLICATING THE 'MPAC' DATA
         JRCRAT = NRATT*NRATT + 1
C--   COMPUTE THE SPACE NEEDED TO DUPLICATE THE 'FORCE' DATA
         JRCFOR = NFORCE + 1
C--   COMPUTE THE NUMBER OF DUPLICATES POSSIBLE
         JRCDUP = MIN(JRCNPC/JRCRAT,MFORCE/JRCFOR)
C--   CHECK THIS VALUE         
         IF (JRCDUP .GT. 0) THEN
            PRINT 9160, JRCDUP
 9160       FORMAT(/I3,' COPIES OF THE PAIR CODES WILL BE USED'/)
            PRINT *,' JRCNPC = ',JRCNPC,'  JRCRAT = ',JRCRAT,
     $           '  QUOTIENT = ',JRCNPC/JRCRAT
            PRINT *,' MFORCE = ',MFORCE,'  JRCFOR = ',JRCFOR,
     $           '  QUOTIENT = ',MFORCE/JRCFOR
         ELSE
            PRINT 9180, JRCRAT, MAXP
 9180       FORMAT(/' THE NUMBER OF PAIR CODES (',I4,') IS LARGER',
     $           ' THAN THE WORK SPACE PROVIDED (',I5,') ',
     $           'IN NONBML AND NBPML'/)
            LERROR = .TRUE.
C     GO TO 5750
            CALL STOPGM('NBPML','FATAL ERROR')
         ENDIF

C     USE JRCILJ AS SCRATCH ARRAY TO PACK THE LJ CODES
         IF (MAXNAT .GT. JRCRAT) THEN
            DO I = 1, NRATT
               DO J = I, NRATT
                  M=NRATT*(I-1)+J
                  N=NRATT*(J-1)+I
                  IPC=I+J*(J-1)/2
                  JRCILJ(M) = MPAC(J, I)
                  JRCILJ(N) = MPAC(J, I)
               ENDDO
            ENDDO
         ELSE
            WRITE(6,*) 'SCRATCH ARRAY FOR MPAC DATA TOO SMALL'
            WRITE(6,*) 'INCREASE DIMENSION OF ARRAY JRCILJ IN',
     .           ' THE ARRAYSIZES SECTION OF QMMM'
cmb - see defaults in mm_read_qmmm_input.F - former vector.h
            CALL STOPGM('NBPML','FATAL ERROR')
         ENDIF

C--   DUPLICATE THIS TABLE NOW
         JRCOFF = 0
         DO JRC1 = 1, JRCDUP
            DO JRC2 = 1, JRCRAT
C     JRCPAC(JRC2+JRCOFF) = MPAC(JRC2,1)
               JRCPAC(JRC2+JRCOFF) = JRCILJ(JRC2)
            ENDDO
            JRCOFF = JRCOFF + JRCRAT
         ENDDO
C--   NOW GENERATE THE OFFSET VALUES TO USE THE MODIFIED TABLE
!$OMP parallel do private(JRC1)
#ifdef _vpp_
!OCL NOALIAS
#endif
         DO JRC1 = 1, MAXP
            JRCBK1(JRC1) = MOD(JRC1-1,JRCDUP)*JRCRAT
            JRCBK2(JRC1) = MOD(JRC1-1,JRCDUP)*JRCFOR
         ENDDO
C--   DUPLICATE THE FORCE CONSTANTS
         IF (LPERTL) THEN
            JRCOFF = 0
            DO JRC1 = 2, JRCDUP
               JRCOFF = JRCOFF + JRCFOR
!$OMP parallel do private(JRC2) shared(JRCOFF)
               DO JRC2 = 1, JRCFOR
                  FORC6(JRCOFF+JRC2) = FORC6(JRC2)
                  FORC12(JRCOFF+JRC2) = FORC12(JRC2)
                  FC12C6(JRCOFF+JRC2) = FC12C6(JRC2)
               ENDDO
            ENDDO
         ELSE
            JRCOFF = 0
            DO JRC1 = 2, JRCDUP
               JRCOFF = JRCOFF + JRCFOR
!$OMP parallel do private(JRC2) shared(JRCOFF)
               DO JRC2 = 1, JRCFOR
                  FORC6(JRCOFF+JRC2) = FORC6(JRC2)
                  FORC12(JRCOFF+JRC2) = FORC12(JRC2)
               ENDDO
            ENDDO
         ENDIF
         
C*****DETERMINE THE ENERGY GROUP TO WHICH EACH ATOM BELONGS
C*****AND STORE THE RESULTS.
         DO I = 1, NATTOT
            ICET(I) = ICHAR(IAGRP(I))
         ENDDO
         
C     
C*****IF THE VIRIAL MUST BE CALCULATED THEN DETERMINE FOR EACH
C*****ATOM TO WHICH SUBMOLECULE IT BELONGS.
C     
         IF (LDOVIR) THEN
            JTOP = 1
            NVIR = 0
            IF (NPM .GE. 1) THEN
               DO NSPI = 1, NSPM
                  NVIR = NVIR + 1
                  LEN = NSP(NSPI) - JTOP + 1
                  DO J1X = 1, LEN
                     NSPTT(J1X-1+JTOP) = NVIR
                  ENDDO
                  JTOP = JTOP + LEN
               ENDDO
C     
               DO NS = 2, NPM
                  DO J1X = 1, NRP
                     JRCIND(J1X) = NSPTT(J1X) + NVIR
                     NSPTT(J1X-1+JTOP) = JRCIND(J1X)
                  ENDDO
                  NVIR = NVIR + NSPM
                  JTOP = JTOP + NRP
               ENDDO
            ENDIF
C     
            IF (NSM .GT. 0) THEN
!$OMP parallel do private(IPQ) shared(NVIR)               
               DO IPQ = 1, NSM
                  JFIRST(IPQ) = NVIR + IPQ
               ENDDO
               DO NN = 1, NRAM
                  CALL SCOPY (NSM, JFIRST(1),1,NSPTT(JTOP),NRAM)
                  JTOP = JTOP + 1
               ENDDO
            ENDIF

         ENDIF
         
C*****DETERMINE FOR EACH CHARGE GROUP WHERE THE DATA IN THE ARRAYS
C*****FOR THE COORDINATES AND FORCES CAN BE FOUND.
C*****STORE THE POINTERS IN IFIRST.
C*****DETERMINE THE NUMBER OF ATOMS IN EACH CHARGE GROUP.
C*****STORE THE NUMBERS IN INAG.
C     
         I = 1
         NI = 1
         DO NS = 1, NPM
            NNOLD = 0
            DO NN = 1, NCAG
               NAG = INC(NN) - NNOLD
               NNOLD = INC(NN)
C     
               IFIRST(NI) = I
               INAG(NI) = NAG
C     
               I = I + NAG
               NI = NI + 1
            ENDDO
         ENDDO
C     
         DO NS = 1, NSM
            IFIRST(NI) = I
            INAG(NI) = NRAM
            I = I + NRAM
            NI = NI + 1
         ENDDO
!$OMP parallel do private(IPQ)         
         DO IPQ = 1, NRAGT
            INTROW(IPQ) = IPQ
         ENDDO

      ENDIF      
C
C------------DONE FIRST CALL SETUP


C     TRANSFORM COORDINATES IF REQUIRED
C     
      IF (LDOTRA) THEN
         CALL TRACO(NATTOT,0,XCOORD,BETA,1,LEVERY)
      ENDIF
C     
C------------INITIALIZE ARRAYS AND GENERAL VARIABLES
C     
      RCUTP2 = RCUTP**2
      RCUTL2 = RCUTL**2
      NATTO3 = NDIM * NATTOT
C     
      ELGLEL = 0.0D0
      ELGLRF = 0.0D0
      ELGLRC = 0.0D0
      ELGLLJ = 0.0D0
      
      EL34EL = 0.0D0
      EL34RF = 0.0D0
      EL34RC = 0.0D0
      EL34LJ = 0.0D0
C
      call mm_AZZERO(VIRLR,NDIM)
C     
!$OMP parallel do private(JRC)
      DO JRC=1,MAXP
         TVIR(JRC,1)=0.0D0
         TVIR(JRC,2)=0.0D0
         TVIR(JRC,3)=0.0D0
      ENDDO  

!$OMP parallel do private(I)
      DO I=1,NUNRE2
         ELREL(I) = 0.0D0
         ELRRF(I) = 0.0D0
         ELRRC(I) = 0.0D0
         ELRLJ(I) = 0.0D0
      ENDDO
         
!$OMP parallel do private(JI)
      DO JI = 1, IVECT*IEPAIR
         TEEL(JI, 1) = 0.0D0
         TERF(JI, 1) = 0.0D0
         TERC(JI, 1) = 0.0D0
         TELJ(JI, 1) = 0.0D0
      ENDDO

C     FILL COORD. ARRAYS
      IF (NDIM .EQ. 3) THEN
         IF (LDOVIR) THEN
!$OMP parallel do private(I)
#ifdef _vpp_
!OCL NOALIAS
#endif
            DO I = 1, NATTOT
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
            DO I = 1, NATTOT
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
            DO I = 1, NATTOT
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
            DO I = 1, NATTOT
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
            DO I = 1, NATTOT
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
            DO I = 1, NATTOT
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
      
      LNEWP = .TRUE.
C     
C------------DONE INITIALIZE
C     GENERATE ARRAY XCG CONTAINING THE REFERENCE POINT OF 
C     EACH CHARGE GROUP FOR DISTANCE MEASUREMENT.
C     SOLUTE : CENTER OF GEOMETRY, SOLVENT : FIRST ATOM
C     
C*****  Code Expanded From Routine:  CLCGEO
C
C------------BEGIN CLCGEO
C
 
C     THE PAIRLIST IS ALWAYS CALCULATED WITH 3D DISTANCES
      NSRC3 = 1
      NDST3 = 1
C
C------------FIRST DO THE SOLUTE ( CENTRE OF GEOMETRY )
C
      DO NPP = 1, NPM
         IFICG = 1
         DO NCG = 1, NCAG
            ILACG = INC(NCG)
            NUMCG = ILACG - IFICG + 1
            IF (NUMCG .GT. 1) THEN
               XTEMP(1) = 0.0D0
               XTEMP(2) = 0.0D0
               XTEMP(3) = 0.0D0
C
               DO J2X = 1, NUMCG
                  XTEMP(1) = XTEMP(1) + XC(NSRC3)
                  XTEMP(2) = XTEMP(2) + YC(NSRC3)
                  XTEMP(3) = XTEMP(3) + ZC(NSRC3)
                  NSRC3 = NSRC3 + 1
               ENDDO
C
               XDST(NDST3) = XTEMP(1)/NUMCG
               YDST(NDST3) = XTEMP(2)/NUMCG
               ZDST(NDST3) = XTEMP(3)/NUMCG
 
            ELSE
               XDST(NDST3) = XC(NSRC3)
               YDST(NDST3) = YC(NSRC3)
               ZDST(NDST3) = ZC(NSRC3)
               NSRC3 = NSRC3 + 1
            ENDIF
            NDST3 = NDST3 + 1
            IFICG = ILACG + 1
         END DO
      END DO
C
C
C------------NOW DO THE SOLVENT (COORDINATE OF THE FIRST ATOM)
C
!$OMP parallel do private(NPP) shared(NDST3,NSRC3,NRAM)
#ifdef _vpp_
!OCL NOALIAS
#endif
      DO NPP = 1, NSM
         XDST(NDST3+NPP-1) = XC(NSRC3+(NPP-1)*NRAM)
         YDST(NDST3+NPP-1) = YC(NSRC3+(NPP-1)*NRAM)
         ZDST(NDST3+NPP-1) = ZC(NSRC3+(NPP-1)*NRAM)
      END DO
C*****  End of Code Expanded From Routine:  CLCGEO
      
      IF (NCALCD .EQ. 3) THEN
C     
C------------DONE INITIALIZE
C     START CALCULATION ( WITH OR WITHOUT LONG RANGE FORCES )
C     
         IF (LDOLNG .OR. LFREE) THEN
            CALL NBWIT3(NATTOT,NPM,NSM,NRAGT,NCALCD,
     $           XCG,XCOORD,XR,FLR,
     $           NSPM,NSP,VIRLR,
     $           INB,JNB,NSZPL,
     $           MAXNRE,NRELKP,IAGRP,
     $           NUNRE2,
     $           ELREL,ELGLEL,EL34EL,
     $           ELRRF,ELGLRF,EL34RF,
     $           ELRRC,ELGLRC,EL34RC,
     $           ELRLJ,ELGLLJ,EL34LJ,
     $           LDOPER,L3D4D,RLAM,
     $           ALPHLJ,ALPHC,NLAM,MMU,
     $           LPION,LEVERY,LFREE)
         ELSE
            CALL NBNONE(NRAGT,3,XCG,INB,JNB,NSZPL)
         ENDIF
         
      ELSEIF (NCALCD .EQ. 4) THEN
         
         IF (LDOLNG .OR. LFREE) THEN
            CALL NBWIT4(NATTOT,NPM,NSM,NRAGT,NCALCD,
     $           XCG,XCOORD,XR,FLR,
     $           NSPM,NSP,VIRLR,
     $           INB,JNB,NSZPL,
     $           MAXNRE,NRELKP,IAGRP,
     $           NUNRE2,
     $           ELREL,ELGLEL,EL34EL,
     $           ELRRF,ELGLRF,EL34RF,
     $           ELRRC,ELGLRC,EL34RC,
     $           ELRLJ,ELGLLJ,EL34LJ,
     $           LDOPER,L3D4D,RLAM,
     $           ALPHLJ,ALPHC,NLAM,MMU,
     $           LPION,LEVERY,LFREE)
         ELSE
            CALL NBNONE(NRAGT,3,XCG,INB,JNB,NSZPL)
         ENDIF
      ENDIF

C     
C------------DONE CALCULATION
C     TRANSFORM BACK COORDINATES IF REQUIRED
C     
      IF (LDOTRA) THEN
         CALL TRACO(NATTOT,0,XCOORD,BETA,-1,LEVERY)
      ENDIF

C     
 900  FORMAT (1X,A6,': ATOM',I5,' SOFT-CORE PERTURBED')
 910  FORMAT (1X,'AND C6 .EQ. 0 BUT C12 .NE. 0 WITH STATE ',
     $     A1,' ATOM TYPE',I5)
 915  FORMAT (1X,'AND CS6 .EQ. 0 BUT CS12 .NE. 0 WITH STATE ',
     $     A1,' ATOM TYPE',I5)
C     
C------------END NBPML
C     ==--------------------------------------------------------------==
      RETURN     
      END
C     ==================================================================
      SUBROUTINE NBNONE(NRAGT,NCALCD,X,INB,JNB,NJNB)
C     ==--------------------------------------------------------------==
COMMSUBR NBNONE
C     SUBROUTINE NBNONE(NRAGT,NCALCD,X,INB,JNB,NJNB)
C     
C     NBNONE CALCULATES A PAIRLIST WITHOUT CALCULATING
C     A LONG RANGE FORCE IN 3D.
C     THIS ROUTINE IS CALLED BY L<NBPML>
COMMEND
C     
C     INCLUDES
C     
      INCLUDE 'coordsz.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'box.h'
      INCLUDE 'vector.h'
      INCLUDE 'nbpml.h'
C     
C     ARGUMENTS
C     
      INTEGER NRAGT,NCALCD,INB(NRAGT),JNB(MAXJNB),NJNB
      real*8 X(NRAGT*NDIM)
C     
C     LOCAL VARIABLES
C     
      INTEGER IX3,NI,NJ,IPQ,LEN,LRLEN
      CHARACTER PRGSTR*(6)
      REAL*8 BOX1, BOX2, BOX3, BOXH1, BOXH2, BOXH3,
     $     BOXIV1, BOXIV2, BOXIV3, XIJ, YIJ, ZIJ
C     
C     DATA
C     
      DATA PRGSTR /'NBNONE'/
C     
C------------BEGIN NBNONE
C     
      BOX1 = BOX(1)
      BOX2 = BOX(2)
      BOX3 = BOX(3)
      IF (LVAC) THEN
         BOX1 = 1.0D10
         BOX2 = 1.0D10
         BOX3 = 1.0D10
      ENDIF
      IF (.NOT. LDOTRA) THEN
         COSB2 = 0.0D0
      ENDIF
      
      BOXH1 = 0.5D0*BOX1
      BOXH2 = 0.5D0*BOX2
      BOXH3 = 0.5D0*BOX3
      BOXIV1 = 1.0D0/BOX1
      BOXIV2 = 1.0D0/BOX2
      BOXIV3 = 1.0D0/BOX3
      BOXOH = BOXH1
      BOXOQ = BOX1*0.75D0
      
      NJNB = 0
      IX3 = 0
C     
C------------LOOP OVER PRIMARY CHARGE GROUPS      
C     
      DO 100 NI = 1, NRAGT - 1
         IF (NJNB + NRAGT - NI .GT. MAXJNB) GO TO 4550
         NJ = NI + 1
         LEN = NRAGT - NI         
C     
C--------LOOP OVER SECONDARY CHARGE GROUPS (EITHER SOLUTE OR SOLVENT)         
C---- CALCULATE DISTANCE AND CONNECTING VECTOR BETWEEN XCGI AND XCG(NJCOR)
C     CORRECT FOR BOUNDARY CONDITIONS IF NECESSARY
C     
!$OMP parallel do private(IPQ,XIJ,YIJ,ZIJ)
#ifdef _vpp_
!OCL NOALIAS 
#endif
         DO IPQ = 1, LEN
            XIJ = (XDST(NI)-XDST(IPQ+NJ-1))*BOXIV1
            YIJ = (YDST(NI)-YDST(IPQ+NJ-1))*BOXIV2
            ZIJ = (ZDST(NI)-ZDST(IPQ+NJ-1))*BOXIV3
            
            VXIJ(IPQ) = (XIJ-ANINT(XIJ))*BOX1
            VYIJ(IPQ) = (YIJ-ANINT(YIJ))*BOX2
            VZIJ(IPQ) = (ZIJ-ANINT(ZIJ))*BOX3
            
            VRIJ2(IPQ) = VXIJ(IPQ) * VXIJ(IPQ) + VYIJ(IPQ)
     .                 * VYIJ(IPQ) + VZIJ(IPQ)*(VZIJ(IPQ)+
     .                   COSB2*VXIJ(IPQ))
         ENDDO
         
         IF (LOCTO) THEN
!$OMP parallel do private(IPQ)
            DO IPQ = 1, LEN
               VEELJ(IPQ) = BOXOQ - DABS(VXIJ(IPQ))-
     .            DABS(VYIJ(IPQ)) - DABS(VZIJ(IPQ))
               IF (VEELJ(IPQ) .LE. 0.0D0) VRIJ2(IPQ) = 
     .             VRIJ2(IPQ) + VEELJ(IPQ)*BOX1
            ENDDO
         ENDIF
C     
C---- DONE DISTANCE CALCULATION
C     CHECK DISTANCE
C     
         CALL WHENFLE (LEN, VRIJ2(1),1,RCUTP2,JRCIND,LRLEN)
C     
C---- INSIDE SHORT-RANGE CUT-OFF, STICK INTO PAIRLIST
C     
         IF (LRLEN .EQ. 0) THEN
            INB(NI) = 0
         ELSE
!$OMP parallel do private(IPQ)
            DO IPQ = 1, LRLEN
               JNB(NJNB + IPQ) = INTROW(NJ+JRCIND(IPQ)-1)
            ENDDO
            INB(NI) = LRLEN
            NJNB = NJNB + LRLEN
         ENDIF
C     
C------------END LOOP OVER PRIMARY CHARGE GROUPS      
C     
 100  CONTINUE
      
      INB(NRAGT) = 0
C     
      PRINT 500,NJNB
      GOTO 9999
      
C     ERRORS
 4550 WRITE(6,*) 'TOO MANY PAIRS WITHIN CUT-OFF DISTANCE'
      WRITE(6,*) 'NUMBER OF PAIRS',NJNB, 
     $     ' INCREASE MAXJNB IN forcesz.h'
      WRITE(6,*) 'ERROR IN ',PRGSTR
      CALL STOPGM('NBPML','FATAL ERROR')
C     
C------------END NBNONE
C     
 500  FORMAT (1X,'NBPML:',I10,' ELEMENTS IN THE PAIRLIST')
 9999 CONTINUE
C     ==--------------------------------------------------------------==
      RETURN     
      END 
C     ==================================================================
      SUBROUTINE NBWIT3(NR,NPM,NSM,NRAGT,NCALCD,
     $     XCG,XCOORD,XR,FLR,
     $     NSPM,NSP,VIRLR,
     $     INB,JNB,NJNB,
     $     MAXNRE,NRELKP,IAGRP,
     $     NUNRE2,
     $     ELREL,ELGLEL,EL34EL,
     $     ELRRF,ELGLRF,EL34RF,
     $     ELRRC,ELGLRC,EL34RC,
     $     ELRLJ,ELGLLJ,EL34LJ,
     $     LDOPER,L3D4D,RLAM,
     $     ALPHLJ,ALPHC,NLAM,MMU,
     $     LPION,LEVERY,LFREE)
C     ==--------------------------------------------------------------==
COMMSUBR NBWIT3 CCCCCC I.G. TIRONI, ZUERICH, JUNE 1996 CCCCCCCCCCCCCCCCCCCCCCCC
C     
C     SUBROUTINE NBWIT3(NR,NPM,NSM,NRAGT,NCALCD,
C     $     XCG,XCOORD,XR,FLR,
C     $     NSPM,NSP,VIRLR,
C     $     INB,JNB,NJNB,
C     $     MAXNRE,NRELKP,IAGRP,
C     $     NUNRE2,
C     $     ELREL,ELGLEL,EL34EL,
C     $     ELRRF,ELGLRF,EL34RF,
C     $     ELRRC,ELGLRC,EL34RC,
C     $     ELRLJ,ELGLLJ,EL34LJ,
C     $     LDOPER,L3D4D,RLAM,
C     $     ALPHLJ,ALPHC,NLAM,MMU,
C     $     LPION,LEVERY,LFREE)
C 
C     NBWIT3 constructs a charge group pairlist and calculates a
C     long range force in 3D.
C     NBWIT3 is called by L<NBPML>.
C     
COMMEND
C     
C     INCLUDES
C     
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'pertsz.h'
      INCLUDE 'pertar.h'
      INCLUDE 'vector.h'
      INCLUDE 'nbpml.h'
C     
C     LOCAL PARAMS
C     
      real*8 ENEPS
      PARAMETER (ENEPS = 1.0D-20)      
C     
C     ARGUMENTS
C     
      LOGICAL LEVERY,LDOPER,L3D4D,LPION,LFREE
      INTEGER NRAGT,NCALCD,INB(NRAGT)
      INTEGER JNB(MAXJNB),NJNB,NR,NPM,NSM
      INTEGER NSPM,NSP(NSPM+1),NUNRE2
      INTEGER MAXNRE,NRELKP(MAXNRE,MAXNRE)
      INTEGER NLAM,MMU
      real*8  RLAM,ALPHLJ,ALPHC
      real*8 XCG(NDIM*NRAGT)
      real*8 XCOORD(NDIM*NR),FLR(NDIM*NR),XR(NDIM*NR)
      real*8 VIRLR(NDIM)
      real*8 ELREL(NUNRE2),ELGLEL,EL34EL
      real*8 ELRRF(NUNRE2),ELGLRF,EL34RF
      real*8 ELRRC(NUNRE2),ELGLRC,EL34RC
      real*8 ELRLJ(NUNRE2),ELGLLJ,EL34LJ
      CHARACTER IAGRP(NR)
C     
C     LOCAL VARIABLES
C     
      real*8 ZEROC6
      PARAMETER (ZEROC6 = 1.0D-20)
      CHARACTER PRGSTR*(6)
      LOGICAL LIZERO, LJZERO, LERROR, LVIR, LIPERT
      INTEGER NRAGP, NRPT,
     $     JRC, IPQ, LRLEN, NSPAIR, LRLEN2,
     $     NI, NJ, JLEN, LEN, I, J, JRCNUM, JRCL1X,
     $     I1X, K1X, MAXLEN, ILAST, ICE, I41X, IVIRL,
     $     JRCTOT, JRCSRT, JRCLEN, JRC1, JRC2,
     $     NRAG, NIJ, NIJB, M
C     
C     PERTURBATION VARIABLES
C     
      real*8 RIJ3D2, SOFTLJ, CA126, CB126, SFLJ, RIJ2, RIJ4, RIJ6,
     $     VELA3D, VLJA3D, VRFA3D, DFA3D, DFRFA3, DLLJA3,
     $     DLELA3, DLRFA3, VELB3D, VLJB3D, VRFB3D ,
     $     DFB3D, DFRFB3, DLLJB3, DLELB3, DLRFB3,
     $     VEL3D, VRF3D, VLJ3D, DF3D, DLLJ3D, DLEL3D, DLRF3D,
     $     VRCA,VRCB,VRC,DLRC,
     $     QIQJA,QIQJB,QRFR,
     $     SOFTEL, SFEL,
     $     RSOFT, RSOFT2, RRF, RRF2, RFSQRT,
     $     RLMA, RL2A, RLMB, RL2B, RLLA, RLLB, RLB, RLA, RDLB, RDLA,
     $     RLARMU, RLBRMU,
     $     RM4D, RM3D, RDM3D, RDM4D,
     $     CA12,CA6,CB12,CB6,
     $     CRA1,CRA12
      
      real*8 BOX1, BOX2, BOX3, BOXH1, BOXH2, BOXH3,
     $     BOXIV1, BOXIV2, BOXIV3, XIJ, YIJ, ZIJ,
     $     TEMP1, DF
C     
C     DATA
C     
      DATA PRGSTR /'NBWIT3'/
C     
C------------BEGIN NBWIT3
C 
      NRAG = NPM*NCAG
      BOX1 = BOX(1)
      BOX2 = BOX(2)
      BOX3 = BOX(3)
      IF (LVAC) THEN
         BOX1 = 1.0D10
         BOX2 = 1.0D10
         BOX3 = 1.0D10
      ENDIF
      IF (.NOT. LDOTRA) THEN
         COSB2 = 0.0D0
      ENDIF
      
      BOXH1 = 0.5D0*BOX1
      BOXH2 = 0.5D0*BOX2
      BOXH3 = 0.5D0*BOX3
      BOXIV1 = 1.0D0/BOX1
      BOXIV2 = 1.0D0/BOX2
      BOXIV3 = 1.0D0/BOX3
      BOXOH = BOXH1
      BOXOQ = BOX1*0.75D0
C     
C------------CALCULATE RMU/RLAMBDA DEPENDANT VARIABLES
C     
      IF (LPERTL) THEN
C     
C--   LOOKS LIKE SOME MACHINES DON'T LIKE 0**0 = 1
C     ... WHICH WE UNDERSTAND ...
C     
         RL2A = RLAM**2
         RLMA = RL2A+RMULOC**2
         RL2B = (1.0D0-RLAM)**2
         RLMB = RL2B+RMULOC**2
C     
         IF (RLAM.EQ.0.0D0 .AND. NLAM.EQ.1) THEN
            RLB = 0.0D0
            RDLB = DBLE(NLAM)
         ELSE
            RLB = RLAM**NLAM
            RDLB = RLAM**(NLAM - 1) * NLAM
         ENDIF
         IF (RLAM.EQ.1.0D0 .AND. NLAM.EQ.1) THEN
            RLA = 0.0D0
            RDLA = -DBLE(NLAM)
         ELSE
            RLA = (1.0D0-RLAM)**NLAM          
            RDLA = - (1.0D0-RLAM)**(NLAM-1)*NLAM
         ENDIF
C     
         RLLB = -(1.0D0-RLAM)*RLB
         RLLA = RLAM*RLA
         RLARMU = RLA*RMULOC
         RLBRMU = RLB*RMULOC
C     
         IF (RMULOC.EQ.0.0D0 .AND. MMU.EQ.1) THEN
            RM4D = 0.0D0
            RDM4D = DBLE(MMU)
         ELSE
            RM4D = RMULOC**MMU
            RDM4D = RMULOC**(MMU-1)*MMU
         ENDIF
         IF (RMULOC.EQ.1.0D0 .AND. MMU.EQ.1) THEN
            RM3D = 0.0D0
            RDM3D = -DBLE(MMU)
         ELSE
            RM3D = (1.0D0-RMULOC)**MMU
            RDM3D = -(1.0D0-RMULOC)**(MMU-1)*MMU
         ENDIF
C     
         IF (.NOT. L4D) THEN
            RDM3D = 0.0D0
         ENDIF
C     
      ENDIF
C     
C     NJNB IS USED TO ACCESS THE PAIRLIST
      
      NRAGP = NCAG*NPM
      NRPT = NRP*NPM
C     
      NJNB = 0
C     
      DO 4450 NI = 1, NRAGT - 1
         IF (NJNB + NRAGT - NI .GT. MAXJNB) GO TO 4550
         NJ = NI + 1
         LEN = NRAGT - NI         
C     
C     CORRECT FOR BOUNDARY CONDITIONS IF NECESSARY
C     
!$OMP parallel do private(IPQ,XIJ,YIJ,ZIJ)
#ifdef _vpp_
!OCL NOALIAS
#endif
         DO IPQ = 1, LEN
            XIJ = (XDST(NI)-XDST(IPQ+NJ-1))*BOXIV1
            YIJ = (YDST(NI)-YDST(IPQ+NJ-1))*BOXIV2
            ZIJ = (ZDST(NI)-ZDST(IPQ+NJ-1))*BOXIV3
            
            VXIJ(IPQ) = (XIJ-ANINT(XIJ))*BOX1
            VYIJ(IPQ) = (YIJ-ANINT(YIJ))*BOX2
            VZIJ(IPQ) = (ZIJ-ANINT(ZIJ))*BOX3
            
            VRIJ2(IPQ) = VXIJ(IPQ)*VXIJ(IPQ) + VYIJ(IPQ)*VYIJ(IPQ)
     .                 + VZIJ(IPQ)*(VZIJ(IPQ)+COSB2*VXIJ(IPQ))
         ENDDO
         
         IF (LOCTO) THEN
!$OMP parallel do private(IPQ)
            DO IPQ = 1, LEN
               VEELJ(IPQ) = BOXOQ - DABS(VXIJ(IPQ))
     .          - DABS(VYIJ(IPQ)) - DABS(VZIJ(IPQ))
               IF (VEELJ(IPQ) .LE. 0.0D0) VRIJ2(IPQ) = 
     .              VRIJ2(IPQ) + VEELJ(IPQ)*BOX1
            ENDDO
         ENDIF
C     
C     CHECK DISTANCE
C     
         CALL WHENFLE (LEN, VRIJ2(1),1,RCUTL2,JRCIND,LRLEN)
C     
         NSPAIR=0
         LRLEN2=0
C     
         IF (LRLEN .EQ. 0) THEN
            INB(NI) = 0
         ELSE
            DO IPQ = 1, LRLEN
               R2D(IPQ) = VRIJ2(JRCIND(IPQ))
               NJL(IPQ) = INTROW(NJ+JRCIND(IPQ)-1)
C     
C*****PUT THE CHARGE GROUPS NJ WHERE R*R.LE.RCUTP2 IN THE NEIGHBOURLIST.
C     
               IF(R2D(IPQ).LE.RCUTP2) THEN
                  NSPAIR=NSPAIR+1
                  JNB(NJNB+NSPAIR)=NJL(IPQ)
               ELSE
                  LRLEN2=LRLEN2+1
                  JFIRST(LRLEN2) = IFIRST(NJL(IPQ))
                  NJL(LRLEN2) = INAG(NJL(IPQ))
               ENDIF
            ENDDO
            
            INB(NI) = NSPAIR
            NJNB = NJNB + NSPAIR

            IF (LRLEN2 .GT. 0) THEN
C     
C*****CALCULATE THE LONG RANGE FORCES AND ENERGY.
C     
               JRCNUM = 1
               JRCL1X = LRLEN2
               J = 0
77001          CONTINUE
               K1X = 0
               
               DO I1X = 1, JRCL1X
                  IF (NJL(I1X) .GE. JRCNUM) THEN
                     J = J + 1
                     JLIST(J) = JFIRST(I1X)
                     K1X = K1X + 1
                     JFIRST(K1X) = JFIRST(I1X) + 1
                     NJL(K1X) = NJL(I1X)
                  ENDIF
               END DO
               JRCL1X = K1X
               JRCNUM = JRCNUM + 1
               IF (JRCL1X .GT. 0) GO TO 77001
               JLEN = J
C     
               MAXLEN=MAX(MAXLEN,JLEN)
               IF (MAXLEN .GT. MAXNAT) THEN
                  PRINT 9023,MAXLEN,MAXNAT
 9023             FORMAT(/' THE LENGTH OF THE SCRATCH ARRAYS (',I5,
     $              ') IS ',' LARGER THAN THE WORK SPACE PROVIDED (',
     $              I5,') IN NBPML'/)
                  LERROR = .TRUE.
                  GO TO 4600
               ENDIF
C     
               ILAST = IFIRST(NI) + INAG(NI) - 1
C               ICE = ICET(ILAST)
C     
C*****GATHER THE COORDINATES FROM ALL ATOMS J.
C     
!$OMP parallel do private(IPQ)
#ifdef _vpp_
!OCL NOALIAS
#endif
               DO IPQ = 1, JLEN
                  XJ(IPQ) = XC(JLIST(IPQ))
                  YJ(IPQ) = YC(JLIST(IPQ))
                  ZJ(IPQ) = ZC(JLIST(IPQ))
C     
C*****GATHER THE FORCES CALCULATED SO FAR.
C     
                  FXTEMP(IPQ) = FLRX(JLIST(IPQ))
                  FYTEMP(IPQ) = FLRY(JLIST(IPQ))
                  FZTEMP(IPQ) = FLRZ(JLIST(IPQ))
C     
C*****CREATE VECTORS FOR THE CHARGES AND FOR GROUP-CODES.
C     
                  JRCILJ(IPQ)=ICET(JLIST(IPQ))
C     
                  IACAJ(IPQ) = IACT(JLIST(IPQ)) + JRCBK1(IPQ)
                  CGAJ(IPQ) = CGAT(JLIST(IPQ))
               ENDDO
C     
               IF (LPERTL) THEN
!$OMP parallel do private(IPQ)
#ifdef _vpp_
!OCL NOALIAS
#endif
                  DO IPQ = 1, JLEN
                     CGBJ(IPQ) = CGBT(JLIST(IPQ))
                     IACBJ(IPQ) = IACBT(JLIST(IPQ))
                     ISCLJJ(IPQ) = ISCLJT(JLIST(IPQ))
                     ISCCJ(IPQ) = ISCCT(JLIST(IPQ))
                     LLPERJ(IPQ) = LLPERT(JLIST(IPQ))
                  ENDDO
               ENDIF
C     
               IF (LPION) THEN
!$OMP parallel do private(IPQ)
                  DO IPQ = 1, JLEN
                     IPICJ(IPQ) = IPICT(JLIST(IPQ))
                  ENDDO
               ENDIF
C     
               IF (LDOVIR) THEN
!$OMP parallel do private(I41X)                  
                  DO I41X = 1, JLEN
                     NSPJ(I41X) = NSPTT(JLIST(I41X))
                  ENDDO
                  IVIRL=0
                  DO IPQ = 1, JLEN
                     IF(NSPJ(IPQ) .NE. NSPTT(ILAST)) THEN
                        IVIRL=IVIRL+1
                        IOCVIR(IVIRL)=IPQ
                     ENDIF
                  ENDDO
                  LVIR=.FALSE.
                  IF(IVIRL .GT. 0) THEN
                     LVIR=.TRUE.
!$OMP parallel do private(IPQ)
                     DO IPQ = 1, JLEN
                        XRJ(IPQ) = XRC(JLIST(IPQ))
                        YRJ(IPQ) = YRC(JLIST(IPQ))
                        ZRJ(IPQ) = ZRC(JLIST(IPQ))
                     ENDDO
                  ENDIF
               ENDIF
C     
C*****PROCES EACH I-ATOM, AND TREAT ALL J-ATOMS IN VECTORMODE.
C     
               DO 4350 I = IFIRST(NI), ILAST
                  ICE = ICET(I)
C     
C     GET THE LJ PARAMETERS
C     
                  IF (LPERTL) THEN
                     LIPERT = LLPERT(I)
!$OMP parallel do private(IPQ,NIJ,NIJB)
                     DO IPQ = 1, JLEN
                        NIJ = JRCPAC(IACAJ(IPQ) + (IACT(I) - 1) *
     $                       NRATT) + JRCBK2(IPQ)
                        NIJB = JRCPAC(IACBJ(IPQ) + (IACBT(I) - 1) *
     $                       NRATT) + JRCBK2(IPQ)
                        VC6A(IPQ) = FORC6(NIJ)
                        VC12A(IPQ) = FORC12(NIJ)
                        VC6B(IPQ) = FORC6(NIJB)
                        VC12B(IPQ) = FORC12(NIJB)
                        AC12C6(IPQ) = FC12C6(NIJ)
                        BC12C6(IPQ) = FC12C6(NIJB)
                     ENDDO
                  ELSE
!$OMP parallel do private(IPQ,NIJ)
                     DO IPQ = 1, JLEN
                        NIJ = JRCPAC(IACAJ(IPQ) + (IACT(I) - 1) *
     $                       NRATT) + JRCBK2(IPQ)
                        VC6A(IPQ) = FORC6(NIJ)
                        VC12A(IPQ) = FORC12(NIJ)
                     ENDDO
                  ENDIF

                  IF (LPION) THEN
                     LIZERO = (IPICT(I) .EQ. 0)
!$OMP parallel do private(IPQ,LJZERO)
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
C*****CALCULATE THE MINIMAL DISTANCE IN ACCORD WITH THE PERIODIC BOX.
C     
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
     $                            VYIJ(IPQ)*VYIJ(IPQ) + 
     $               VZIJ(IPQ)*(VZIJ(IPQ)+COSB2*VXIJ(IPQ))
                  ENDDO
                  IF (LOCTO) THEN
!$OMP parallel do private(IPQ,TEMP1)
                     DO IPQ = 1, JLEN
                        TEMP1 = BOXOQ - DABS(VXIJ(IPQ)) - 
     .                     DABS(VYIJ(IPQ)) - DABS(VZIJ(IPQ))
                        IF (TEMP1 .LE. 0.0D0) THEN
                           VRIJ2(IPQ) = VRIJ2(IPQ) + TEMP1*BOX1
                           VXIJ(IPQ) = VXIJ(IPQ) -
     $                          SIGN(BOXOH,VXIJ(IPQ))
                           VYIJ(IPQ) = VYIJ(IPQ) -
     $                          SIGN(BOXOH,VYIJ(IPQ))
                           VZIJ(IPQ) = VZIJ(IPQ) -
     $                          SIGN(BOXOH,VZIJ(IPQ))
                        ENDIF
                     ENDDO
                  ENDIF
C     
!$OMP parallel do private(IPQ)
                  DO IPQ = 1, JLEN
                     VRIJ(IPQ) = DSQRT(VRIJ2(IPQ))
                     VRINV(IPQ) = 1.0D0 / VRIJ(IPQ)
                     VRINV2(IPQ) = VRINV(IPQ)**2
                     VRINV6(IPQ) = VRINV2(IPQ)**3
                  ENDDO
                  
                  IF (LPERTL) THEN
                     IF (LIPERT) THEN
                        DO IPQ = 1, JLEN
C     
C---- A SOFT-CORE PERTURBATION
C     CALCULATION IS DONE
C     GET STATE B PARAMETERS AND USE COMBINATION RULE FOR SOFT-CORE
C     ALPHA PARAMETERS
C     
                           QIQJA = CGAT(I) * CGAJ(IPQ)
                           QIQJB = CGBT(I) * CGBJ(IPQ)
C     
                           CA12 = VC12A(IPQ)
                           CA6 = VC6A(IPQ)
                           CB12 = VC12B(IPQ)
                           CB6  = VC6B(IPQ)
C     
                           SOFTLJ = ALPHLJ * 
     $                          (1 - ISCLJT(I) * ISCLJJ(IPQ))
                           SOFTEL = ALPHC * 
     $                          (1 - ISCCT(I) * ISCCJ(IPQ))
C     
                           CA126 = AC12C6(IPQ)
                           CB126 = BC12C6(IPQ)
C     
C---- SOFT-CORE PERTURBATION
C     CALCULATE FULL 3D CONTRIBUTION
C     
                           RIJ3D2 = VRIJ2(IPQ)
                           RIJ2 = RIJ3D2
                           RIJ4 = VRIJ2(IPQ) * VRIJ2(IPQ)
                           RIJ6 = RIJ4 * VRIJ2(IPQ)
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
                           DFA3D  = CRA1 * (CRA1 * CRA12 + VLJA3D)
C     DV/DL
                           DLLJA3 = -2.0D0 * SFLJ * DFA3D
                           DLELA3 = -VELA3D * SOFTEL * RSOFT2
                           DLRFA3 = 1.5D0 * QRFR*SOFTEL*RRF2
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
                           DLRFB3 = 1.5D0 * QRFR*SOFTEL*RRF2
C     FORCES
                           DFB3D = RSOFT2*VELB3D + 6.0D0*RIJ4*DFB3D
C     
C--   COMBINE 3D STATE A AND B
C     
                           VEL3D  =  RLB  * VELB3D + RLA  * VELA3D
                           VRF3D  =  RLB  * VRFB3D + RLA  * VRFA3D
                           VLJ3D  =  RLB  * VLJB3D + RLA  * VLJA3D
                           DF3D   =  RLB  * (DFB3D + DFRFB3)
     $                             + RLA  * (DFA3D + DFRFA3)
                           DLLJ3D =  RDLB * VLJB3D + RLLB * DLLJB3
     $                             + RDLA * VLJA3D + RLLA * DLLJA3
                           DLEL3D =  RDLB * VELB3D + RLLB * DLELB3 
     $                             + RDLA * VELA3D + RLLA * DLELA3
                           DLRF3D =  RDLB * VRFB3D + RLLB * DLRFB3
     $                             + RDLA * VRFA3D + RLLA * DLRFA3

C     
C---- DONE SOFT-CORE 3D (FULL CONTRIBUTION)
C     
C     
C---- THE DISTANCE INDEPENDANT REACTION FIELD COMPONENT IS
C     DIMENSIONALITY INDEPENDENT. CALCULATE IT HERE.
C     
                           VRCA = RFC*QIQJA
                           VRCB = RFC*QIQJB
                           VRC  = RLB*VRCB  + RLA*VRCA
                           DLRC = RDLB*VRCB + RDLA*VRCA
C     
                           VEELJ(IPQ) = VEL3D
                           VERFJ(IPQ) = VRF3D 
                           VERCJ(IPQ) = VRC
                           VLJJ(IPQ) = VLJ3D
                           EGLELJ(IPQ) = DLEL3D 
                           EGLRFJ(IPQ) = DLRF3D
                           EGLRCJ(IPQ) = DLRC
                           EGLLJJ(IPQ) = DLLJ3D
                           
                           DF3(IPQ)   = DF3D 

                        ENDDO

                     ELSE
                        DO IPQ = 1, JLEN
                           IF (LLPERJ(IPQ)) THEN
                              QIQJA = CGAT(I) * CGAJ(IPQ)
                              QIQJB = CGBT(I) * CGBJ(IPQ)
C     
                              CA12 = VC12A(IPQ)
                              CA6 = VC6A(IPQ)
                              CB12 = VC12B(IPQ)
                              CB6  = VC6B(IPQ)
C     
                              SOFTLJ = ALPHLJ * 
     $                             (1 - ISCLJT(I) * ISCLJJ(IPQ))
                              SOFTEL = ALPHC * 
     $                             (1 - ISCCT(I) * ISCCJ(IPQ))
C     
                              CA126 = AC12C6(IPQ)
                              CB126 = BC12C6(IPQ)

C     
C---- SOFT-CORE PERTURBATION
C     CALCULATE FULL 3D CONTRIBUTION
C     
                              RIJ3D2 = VRIJ2(IPQ)
                              RIJ2 = RIJ3D2
                              RIJ4 = VRIJ2(IPQ) * VRIJ2(IPQ)
                              RIJ6 = RIJ4 * VRIJ2(IPQ)
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
                              DLRFA3 = 1.5D0 * QRFR*SOFTEL*RRF2
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
                              DLRFB3 = 1.5D0 * QRFR*SOFTEL*RRF2
C     FORCES
                              DFB3D = RSOFT2*VELB3D+6.0D0*RIJ4*DFB3D 
C     
C--   COMBINE 3D STATE A AND B
C     
                              VEL3D  =  RLB  * VELB3D + RLA  * VELA3D
                              VRF3D  =  RLB  * VRFB3D + RLA  * VRFA3D
                              VLJ3D  =  RLB  * VLJB3D + RLA  * VLJA3D
                              DF3D   =  RLB  * (DFB3D + DFRFB3)
     $                                + RLA  * (DFA3D + DFRFA3)
                              DLLJ3D =  RDLB * VLJB3D + RLLB * DLLJB3
     $                                + RDLA * VLJA3D + RLLA * DLLJA3
                              DLEL3D =  RDLB * VELB3D + RLLB * DLELB3 
     $                                + RDLA * VELA3D + RLLA * DLELA3
                              DLRF3D =  RDLB * VRFB3D + RLLB * DLRFB3
     $                                + RDLA * VRFA3D + RLLA * DLRFA3
C     
C---- DONE SOFT-CORE 3D (FULL CONTRIBUTION)
C     

C     
C---- THE DISTANCE INDEPENDANT REACTION FIELD COMPONENT IS
C     DIMENSIONALITY INDEPENDENT. CALCULATE IT HERE.
C     
                              VRCA = RFC*QIQJA
                              VRCB = RFC*QIQJB
                              VRC  = RLB*VRCB  + RLA*VRCA
                              DLRC = RDLB*VRCB + RDLA*VRCA
C     
                              VEELJ(IPQ) = VEL3D
                              VERFJ(IPQ) = VRF3D 
                              VERCJ(IPQ) = VRC
                              VLJJ(IPQ) = VLJ3D
                              EGLELJ(IPQ) = DLEL3D 
                              EGLRFJ(IPQ) = DLRF3D
                              EGLRCJ(IPQ) = DLRC
                              EGLLJJ(IPQ) = DLLJ3D
                              
                              DF3(IPQ)   = DF3D 
                              
                           ELSE
                              QIQJA = CGAT(I)*CGAJ(IPQ)
                              VERFJ(IPQ) = QIQJA * RFE * VRIJ2(IPQ)
                              VERCJ(IPQ) = QIQJA * RFC
                              VEELJ(IPQ) = QIQJA*VRINV(IPQ)
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
                           ENDIF

                        ENDDO
                        
                     ENDIF
!$OMP parallel do private(IPQ) reduction(+:ELGLEL,ELGLRF,ELGLRC,ELGLLJ)
                     DO  IPQ = 1, JLEN
                        ELGLEL = ELGLEL + EGLELJ(IPQ)
                        ELGLRF = ELGLRF + EGLRFJ(IPQ)
                        ELGLRC = ELGLRC + EGLRCJ(IPQ)
                        ELGLLJ = ELGLLJ + EGLLJJ(IPQ)
                     ENDDO
                     
                  ELSE
                     
                     IF (LPION) THEN
!$OMP parallel do private(IPQ,QIQJA,DF)
                        DO IPQ = 1, JLEN
                           QIQJA = CGAT(I)*CGAJ(IPQ) * PIFACT(IPQ)
                           VERFJ(IPQ) = QIQJA * RFE * VRIJ2(IPQ)
                           VERCJ(IPQ) = QIQJA * RFC
                           VEELJ(IPQ) = QIQJA*VRINV(IPQ)
                           VLJJ(IPQ) = VRINV6(IPQ)*(VC12A(IPQ)*VRINV6(
     $                          IPQ)-VC6A(IPQ)) * PIFACT(IPQ)
                           DF = VRINV6(IPQ)*(2.0D0*VRINV6(IPQ)*
     $                          VC12A(IPQ)-VC6A(IPQ)) * PIFACT(IPQ)
                           DF3(IPQ) = (6.0D0 * DF + VEELJ(IPQ))*
     $                          VRINV2(IPQ) +  QIQJA * RFF
                        ENDDO
                     ELSE
!$OMP parallel do private(IPQ,QIQJA,DF)
                        DO IPQ = 1, JLEN
                           QIQJA = CGAT(I)*CGAJ(IPQ)
                           VERFJ(IPQ) = QIQJA * RFE * VRIJ2(IPQ)
                           VERCJ(IPQ) = QIQJA * RFC
                           VEELJ(IPQ) = QIQJA*VRINV(IPQ)
                           VLJJ(IPQ) = VRINV6(IPQ)*(VC12A(IPQ)*VRINV6(
     $                          IPQ)-VC6A(IPQ))
                           DF = VRINV6(IPQ)*(2.0D0*VRINV6(IPQ)*
     $                          VC12A(IPQ)-VC6A(IPQ))
                           DF3(IPQ) = (6.0D0 * DF + VEELJ(IPQ))*
     $                          VRINV2(IPQ) +  QIQJA * RFF
                        ENDDO
                     ENDIF
                     
                  ENDIF
C     
                  JRCTOT=JLEN
                  JRCSRT=0
 4590             CONTINUE
                  JRCLEN=MIN(IVECT,JRCTOT)
CDIR$ SHORTLOOP                  
                  DO JRC=1,JRCLEN
                     
                     TEEL(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TEEL(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VEELJ(JRC+JRCSRT)
                     
                     TERF(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TERF(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VERFJ(JRC+JRCSRT)
                     
                     TERC(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TERC(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VERCJ(JRC+JRCSRT)
                     
                     TELJ(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TELJ(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VLJJ(JRC+JRCSRT)
                     
                  ENDDO

                  JRCSRT=JRCSRT+JRCLEN
                  JRCTOT=JRCTOT-JRCLEN
                  IF(JRCTOT.GT.0) GOTO 4590
C     
C*****CALCULATE AND STORE THE FORCES.
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
                        
                        DO IPQ=1,IVIRL
                           TVIR(IPQ,1) = TVIR(IPQ,1) -
     $                          VXH(IOCVIR(IPQ))*
     $                          (VXIJ(IOCVIR(IPQ)) +
     $                          XRJ(IOCVIR(IPQ))-XRC(I))
                           TVIR(IPQ,2) = TVIR(IPQ,2) -
     $                          VYH(IOCVIR(IPQ))*
     $                          (VYIJ(IOCVIR(IPQ)) +
     $                          YRJ(IOCVIR(IPQ))-YRC(I))
                           TVIR(IPQ,3) = TVIR(IPQ,3) -
     $                          VZH(IOCVIR(IPQ))*
     $                          (VZIJ(IOCVIR(IPQ)) +
     $                          ZRJ(IOCVIR(IPQ))-ZRC(I))
                        ENDDO
                     ENDIF
                  ENDIF
C     
                  DO IPQ = 1, JLEN
                     FLRX(I) = FLRX(I) + VXH(IPQ)
                     FLRY(I) = FLRY(I) + VYH(IPQ)
                     FLRZ(I) = FLRZ(I) + VZH(IPQ)
                  ENDDO
C     
 4350          CONTINUE
C     
!$OMP parallel do private(IPQ)
               DO IPQ = 1, JLEN
                  FLRX(JLIST(IPQ)) = FXTEMP(IPQ)
                  FLRY(JLIST(IPQ)) = FYTEMP(IPQ)
                  FLRZ(JLIST(IPQ)) = FZTEMP(IPQ)
               ENDDO
            ENDIF
         ENDIF
 4450 CONTINUE
C     
      INB(NRAGT) = 0
C--   UPDATE THE 'ELREL' VALUES
      DO JRC1=1,NUNRE2
         DO JRC2=1,IVECT
            ELREL(JRC1) = ELREL(JRC1) + TEEL(JRC1, JRC2)
            ELRRF(JRC1) = ELRRF(JRC1) + TERF(JRC1, JRC2)
            ELRRC(JRC1) = ELRRC(JRC1) + TERC(JRC1, JRC2)
            ELRLJ(JRC1) = ELRLJ(JRC1) + TELJ(JRC1, JRC2)
         ENDDO
      ENDDO
C
!$OMP parallel do private(JRC)
      DO JRC = 1, NR
         FLR((JRC-1)*3+1) = FLRX(JRC)
         FLR((JRC-1)*3+2) = FLRY(JRC)
         FLR((JRC-1)*3+3) = FLRZ(JRC)
      ENDDO
C     
      IF (LDOTRA) CALL TRACO (NR, 0, FLR, BETA, -1, LEVERY)
C     
      IF (LDOVIR) THEN
         
         DO JRC=1,MAXLEN
            VIRLR(1)=VIRLR(1)+TVIR(JRC,1)
            VIRLR(2)=VIRLR(2)+TVIR(JRC,2)
            VIRLR(3)=VIRLR(3)+TVIR(JRC,3)
         ENDDO
         VIRLR(1) = VIRLR(1)*0.5D0
         VIRLR(2) = VIRLR(2)*0.5D0
         VIRLR(3) = VIRLR(3)*0.5D0
      ENDIF
C     
      IF (LDOTRA) THEN
         CALL TRACO (NR, 0, XCOORD, BETA, -1, LEVERY)
      ENDIF
C     
C     
C----------RETURN DERIVATIVE VARIABLES ACCORDING TO REQUEST
C     
      IF (.NOT.L3D4D) THEN
C     
C--   IN THIS CASE THESE WOULD BE *WRONG*
C     
        EL34EL = 0.0D0
        EL34RF = 0.0D0
        EL34RC = 0.0D0
        EL34LJ = 0.0D0
      ENDIF
      IF (.NOT.LDOPER) THEN
C     
C--   IN THIS CASE THESE WOULD BE OK, BUT SINCE BONDED ROUTINES
C     WON'T RETURN THEM, RETURNING THEM MIGHT GENERATE CONFUSION... 
C     
        ELGLEL = 0.0D0
        ELGLRF = 0.0D0
        ELGLRC = 0.0D0
        ELGLLJ = 0.0D0
      ENDIF

      PRINT 9040, NJNB      
 9040 FORMAT (1X,'NBPML:',I10,' ELEMENTS IN THE PAIRLIST')

      RETURN
C     
C     
C*****ERROR EXIT
C     
C     
 9100 FORMAT (A, I7)
 4550 CONTINUE
      PRINT *,PRGSTR,' PAIRLIST OVERFLOW:'
      PRINT *,'MAXIMUM SIZE: MAXJNB =',MAXJNB
      PRINT *,'FOUND: ',NJNB
C     
C     
C*****DUMP FROM SUBROUTINE.
C     
C     
 4600 CONTINUE
C     
 9120 FORMAT (10(1X,I9))
 9140 FORMAT (10(1X,G11.4))
 9160 FORMAT (/,1X,20A)
 9180 FORMAT (1X,A,'(', I2, '):', 3(5X, 3(F8.3, 2X)))
 9200 FORMAT (1X, 10L10)
C     
 5750 CONTINUE
      PRINT '(A)', '1THE SUBROUTINE NBPML HAS GONE WRONG. ',
     1     ' THE VARIABLES ARE:'
      PRINT 9160, '      NPM       NRP      NCAG       NSM',
     1     '      NRAM       NTB      NSPM      NJNB'
      PRINT 9120, NPM, NRP, NCAG, NSM, NRAM, NTB, NSPM, NJNB
      PRINT 9160, '     BOX(1)      BOX(2)      BOX(3)',
     1     '        BETA       RCUTP       RCUTL'
      PRINT 9140, (BOX(I), I = 1, 3), BETA, DSQRT(RCUTP2), 
     $     DSQRT(RCUTL2)
C     PRINT 9160, '   NRE(1)    NRE(2)    NRE(3)    NRE(4)'
C     PRINT 9120, (NRE(K), K = 1, 4)
      PRINT 9160, '   ELREL(1..10)'
      PRINT 9140, (ELREL(I), I = 1, 10)
      PRINT 9160, '   VIRLR(1..3)'
      PRINT 9140, (VIRLR(I), I = 1, 3)
      PRINT 9160, '    NSP(1..NSPM)'
      PRINT 9120, (NSP(I), I = 1, NSPM)
      PRINT 9160, '    INC(1..20)'
      PRINT 9120, (INC(I), I = 1, 20)
      PRINT 9160, '     CG(1..20)'
      PRINT 9140, (CG(I), I = 1, 20)
      PRINT 9160, '    CGS(1..NRAM)'
      PRINT 9140, (CGS(I), I = 1, NRAM)
      PRINT 9160, '      X(1..45)'
      DO 4650 I = 1, 36, 9
         PRINT 9180, 'XCOORD', I, (XCOORD(M), M = I, I + 8)
 4650 CONTINUE
      PRINT 9160, '    FLR(1..30)'
      PRINT 9140, (FLR(I), I = 1, 30)
      PRINT 9160, '     XR(1..10)'
      PRINT 9140, (XR(I), I = 1, 10)
      PRINT 9160, '    INB(1..20)'
      PRINT 9120, (INB(I), I = 1, 20)
      PRINT 9160, '    JNB(1..200)'
      PRINT 9120, (JNB(I), I = 1, 200)
      PRINT 9160, '                  LOCAL VARIABLES:'
      PRINT 9160, '    LDOVIR     LVIR    LERROR',
     1     '    LDOTRA'
      PRINT 9200, LDOVIR, LVIR, LERROR, LDOTRA
      PRINT 9160, '       ZERO      RCUTP2      RCUTL2       COSB2'
      PRINT 9140, RCUTP2, RCUTL2, COSB2
      PRINT 9160, '     NRAG     NRAGT ',
     1     '        NR'
      PRINT 9120, NRAG, NRAGT, NR
      PRINT 9160, '        I       ICE        NI',
     1     '        NN     NNOLD', '        NS       NAG'
C     PRINT 9120, I, ICE, NI, NN, NNOLD, NS, NAG
      PRINT 9160, '        I     LRLEN    LRLEN2      JLEN'
      PRINT 9120, I, LRLEN, LRLEN2, JLEN
      PRINT 9160, 'BOXH(1..3), BOXINV(1..3)'
      PRINT 9140, (BOXH(I), I = 1, 3), (BOXINV(I), I = 1, 3)
      CALL STOPGM('NBPML','FATAL ERROR')

C     END NBWIT3
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE NBWIT4(NR,NPM,NSM,NRAGT,NCALCD,
     $     XCG,XCOORD,XR,FLR,
     $     NSPM,NSP,VIRLR,
     $     INB,JNB,NJNB,
     $     MAXNRE,NRELKP,IAGRP,
     $     NUNRE2,
     $     ELREL,ELGLEL,EL34EL,
     $     ELRRF,ELGLRF,EL34RF,
     $     ELRRC,ELGLRC,EL34RC,
     $     ELRLJ,ELGLLJ,EL34LJ,
     $     LDOPER,L3D4D,RLAM,
     $     ALPHLJ,ALPHC,NLAM,MMU,
     $     LPION,LEVERY,LFREE)
C     ==--------------------------------------------------------------==
COMMSUBR NBWIT4 CCCCCC I.G. TIRONI, ZUERICH, JUNE 1996 CCCCCCCCCCCCCCCCCCCCCCCC
C     
C     SUBROUTINE NBWIT4(NR,NPM,NSM,NRAGT,NCALCD,
C     $     XCG,XCOORD,XR,FLR,
C     $     NSPM,NSP,VIRLR,
C     $     INB,JNB,NJNB,
C     $     MAXNRE,NRELKP,IAGRP,
C     $     NUNRE2,
C     $     ELREL,ELGLEL,EL34EL,
C     $     ELRRF,ELGLRF,EL34RF,
C     $     ELRRC,ELGLRC,EL34RC,
C     $     ELRLJ,ELGLLJ,EL34LJ,
C     $     LDOPER,L3D4D,RLAM,
C     $     ALPHLJ,ALPHC,NLAM,MMU,
C     $     LPION,LEVERY,LFREE)
C 
C     NBWIT4 constructs a charge group pairlist and calculates a
C     long range force in 4D.
C     NBWIT4 is called by L<NBPML>.
C     
COMMEND
C     
C     INCLUDES
C     
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'pertsz.h'
      INCLUDE 'pertar.h'
      INCLUDE 'vector.h'
      INCLUDE 'nbpml.h'
C     
C     LOCAL PARAMS
C     
      real*8 ENEPS
      PARAMETER (ENEPS = 1.D-20)      
C     
C     ARGUMENTS
C     
      LOGICAL LEVERY,LDOPER,L3D4D,LPION,LFREE
      INTEGER NRAGT,NCALCD,INB(NRAGT)
      INTEGER JNB(MAXJNB),NJNB,NR,NPM,NSM
      INTEGER NSPM,NSP(NSPM+1),NUNRE2
      INTEGER MAXNRE,NRELKP(MAXNRE,MAXNRE)
      INTEGER NLAM,MMU
      real*8  RLAM,ALPHLJ,ALPHC
      real*8 XCG(NDIM*NRAGT)
      real*8 XCOORD(NDIM*NR),FLR(NDIM*NR),XR(NDIM*NR)
      real*8 VIRLR(NDIM)
      real*8 ELREL(NUNRE2),ELGLEL,EL34EL
      real*8 ELRRF(NUNRE2),ELGLRF,EL34RF
      real*8 ELRRC(NUNRE2),ELGLRC,EL34RC
      real*8 ELRLJ(NUNRE2),ELGLLJ,EL34LJ
      CHARACTER IAGRP(NR)
C     
C     LOCAL VARIABLES
C     
      real*8 ZEROC6
      PARAMETER (ZEROC6 = 1.0D-20)

      CHARACTER PRGSTR*(6)
      LOGICAL LERROR, LIPERT

      INTEGER NRAGP, NRPT,
     $     JRC, IPQ, LRLEN, NSPAIR, LRLEN2,
     $     NI, NJ, JLEN, LEN, I, J, JRCNUM, JRCL1X,
     $     I1X, K1X, MAXLEN, ILAST, ICE,
     $     JRCTOT, JRCSRT, JRCLEN, JRC1, JRC2,
     $     NRAG, NIJ, NIJB, M
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
     $     QIQJA, QIQJB, RIJ2, CA12,CA6,CB12,CB6,
     $     CRA1,CRA12
      
      real*8 BOX1, BOX2, BOX3, BOX4, BOXH1, BOXH2, BOXH3, BOXH4,
     $     BOXIV1, BOXIV2, BOXIV3, BOXIV4, XIJ, YIJ, ZIJ, WIJ,
     $     TEMP1, DF
C     
C     DATA
C     
      DATA PRGSTR /'NBWIT4'/
C     
C------------BEGIN NBWIT4
C 

      NRAG = NPM*NCAG
      BOX1 = BOX(1)
      BOX2 = BOX(2)
      BOX3 = BOX(3)
      BOX4 = BOX(4)
      IF (LVAC) THEN
         BOX1 = 1.0D10
         BOX2 = 1.0D10
         BOX3 = 1.0D10
         BOX4 = 1.0D10
      ENDIF
      IF (.NOT. LDOTRA) THEN
         COSB2 = 0.0D0
      ENDIF
      
      BOXH1 = 0.5D0*BOX1
      BOXH2 = 0.5D0*BOX2
      BOXH3 = 0.5D0*BOX3
      BOXH4 = 0.5D0*BOX4
      BOXIV1 = 1.0D0/BOX1
      BOXIV2 = 1.0D0/BOX2
      BOXIV3 = 1.0D0/BOX3
      BOXIV4 = 1.0D0/BOX4
      BOXOH = BOXH1
      BOXOQ = BOX1*0.75D0
C     
C------------CALCULATE RMU/RLAMBDA DEPENDANT VARIABLES
C     
      IF (LPERTL) THEN
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
            RDLB = RLAM**(NLAM-1)*NLAM
         ENDIF
         IF (RLAM.EQ.1.0D0 .AND. NLAM.EQ.1) THEN
            RLA = 0.0D0
            RDLA = -DBLE(NLAM)
         ELSE
            RLA = (1.0D0-RLAM)**NLAM          
            RDLA = -(1.0D0-RLAM)**(NLAM-1)*NLAM
         ENDIF
C     
         RLLB = -(1.0D0-RLAM)*RLB
         RLLA = RLAM*RLA
         RLARMU = RLA*RMULOC
         RLBRMU = RLB*RMULOC
C     
         IF (RMULOC.EQ.0.0D0 .AND. MMU.EQ.1) THEN
            RM4D = 0.0D0
            RDM4D = DBLE(MMU)
         ELSE
            RM4D = RMULOC**MMU
            RDM4D = RMULOC**(MMU - 1) * MMU
         ENDIF
         IF (RMULOC.EQ.1.0D0 .AND. MMU.EQ.1) THEN
            RM3D = 0.0D0
            RDM3D = -DBLE(MMU)
         ELSE
            RM3D = (1.0D0-RMULOC)**MMU
            RDM3D = -(1.0D0-RMULOC)**(MMU-1)*MMU
         ENDIF
C     
         IF (.NOT. L4D) THEN
            RDM3D = 0.0D0
         ENDIF
C     
      ENDIF
C     
C     NJNB IS USED TO ACCESS THE PAIRLIST
      
      NRAGP = NCAG*NPM
      NRPT = NRP*NPM
C     
      NJNB = 0
C     
      DO 4450 NI = 1, NRAGT - 1
         IF (NJNB + NRAGT - NI .GT. MAXJNB) GO TO 4550
         NJ = NI + 1
         LEN = NRAGT - NI         
         
C     
C     CORRECT FOR BOUNDARY CONDITIONS IF NECESSARY
C    jusq 
!$OMP parallel do private(IPQ,XIJ,YIJ,ZIJ)
#ifdef _vpp_
!OCL NOALIAS
#endif
         DO IPQ = 1, LEN
            XIJ = (XDST(NI)-XDST(IPQ+NJ-1))*BOXIV1
            YIJ = (YDST(NI)-YDST(IPQ+NJ-1))*BOXIV2
            ZIJ = (ZDST(NI)-ZDST(IPQ+NJ-1))*BOXIV3
            
            VXIJ(IPQ) = (XIJ-ANINT(XIJ))*BOX1
            VYIJ(IPQ) = (YIJ-ANINT(YIJ))*BOX2
            VZIJ(IPQ) = (ZIJ-ANINT(ZIJ))*BOX3
            
            VRIJ2(IPQ) = VXIJ(IPQ)*VXIJ(IPQ)+VYIJ(IPQ)*VYIJ(IPQ)
     .                  +VZIJ(IPQ)*(VZIJ(IPQ)+COSB2*VXIJ(IPQ))
         ENDDO
         
         IF (LOCTO) THEN
!$OMP parallel do private(IPQ)
            DO IPQ = 1, LEN
               VEELJ(IPQ) = BOXOQ - DABS(VXIJ(IPQ)) -
     .            DABS(VYIJ(IPQ)) - DABS(VZIJ(IPQ))
               IF (VEELJ(IPQ) .LE. 0.0D0) VRIJ2(IPQ) = 
     .                    VRIJ2(IPQ) + VEELJ(IPQ)*BOX1
            ENDDO
         ENDIF
C     
C     CHECK DISTANCE
C     
         CALL WHENFLE (LEN, VRIJ2(1),1,RCUTL2,JRCIND,LRLEN)
C     
         NSPAIR=0
         LRLEN2=0
C     
         IF (LRLEN .EQ. 0) THEN
            INB(NI) = 0
         ELSE
            DO IPQ = 1, LRLEN
               R2D(IPQ) = VRIJ2(JRCIND(IPQ))
               NJL(IPQ) = INTROW(NJ+JRCIND(IPQ)-1)
C     
C*****PUT THE CHARGE GROUPS NJ WHERE R*R.LE.RCUTP2 IN THE NEIGHBOURLIST.
C     
               IF(R2D(IPQ).LE.RCUTP2) THEN
                  NSPAIR=NSPAIR+1
                  JNB(NJNB+NSPAIR)=NJL(IPQ)
               ELSE
                  LRLEN2=LRLEN2+1
                  JFIRST(LRLEN2) = IFIRST(NJL(IPQ))
                  NJL(LRLEN2) = INAG(NJL(IPQ))
               ENDIF
            ENDDO
            
            INB(NI) = NSPAIR
            NJNB = NJNB + NSPAIR

            IF (LRLEN2 .GT. 0) THEN
C     
C*****CALCULATE THE LONG RANGE FORCES AND ENERGY.
C     
               JRCNUM = 1
               JRCL1X = LRLEN2
               J = 0
77001          CONTINUE
               K1X = 0
               
               DO I1X = 1, JRCL1X
                  IF (NJL(I1X) .GE. JRCNUM) THEN
                     J = J + 1
                     JLIST(J) = JFIRST(I1X)
                     K1X = K1X + 1
                     JFIRST(K1X) = JFIRST(I1X) + 1
                     NJL(K1X) = NJL(I1X)
                  ENDIF
               END DO
               JRCL1X = K1X
               JRCNUM = JRCNUM + 1
               IF (JRCL1X .GT. 0) GO TO 77001
               JLEN = J
C     
               MAXLEN=MAX(MAXLEN,JLEN)
               IF (MAXLEN .GT. MAXNAT) THEN
                  PRINT 9023,MAXLEN,MAXNAT
 9023             FORMAT(/' THE LENGTH OF THE SCRATCH ARRAYS (',I5,
     $                 ') IS ',' LARGER THAN THE WORK SPACE PROVIDED (',
     $                 I5,') IN NBPML'/)
                  LERROR = .TRUE.
                  GO TO 4600
               ENDIF
C     
               ILAST = IFIRST(NI) + INAG(NI) - 1
C               ICE = ICET(ILAST)
C     
C*****GATHER THE COORDINATES FROM ALL ATOMS J.
C
!$OMP parallel do private(IPQ)
#ifdef _vpp_
!OCL NOALIAS
#endif
               DO IPQ = 1, JLEN
                  XJ(IPQ) = XC(JLIST(IPQ))
                  YJ(IPQ) = YC(JLIST(IPQ))
                  ZJ(IPQ) = ZC(JLIST(IPQ))
                  WJ(IPQ) = WC(JLIST(IPQ))
C     
C*****GATHER THE FORCES CALCULATED SO FAR.
C     
                  FXTEMP(IPQ) = FLRX(JLIST(IPQ))
                  FYTEMP(IPQ) = FLRY(JLIST(IPQ))
                  FZTEMP(IPQ) = FLRZ(JLIST(IPQ))
                  FWTEMP(IPQ) = FLRW(JLIST(IPQ))
C     
C*****CREATE VECTORS FOR THE CHARGES AND FOR GROUP-CODES.
C     
                  JRCILJ(IPQ)=ICET(JLIST(IPQ))
C     
                  IACAJ(IPQ) = IACT(JLIST(IPQ)) + JRCBK1(IPQ)
                  CGAJ(IPQ) = CGAT(JLIST(IPQ))
               ENDDO
C     
               IF (LPERTL) THEN
!$OMP parallel do private(IPQ)
#ifdef _vpp_
!OCL NOALIAS
#endif
                  DO IPQ = 1, JLEN
                     CGBJ(IPQ) = CGBT(JLIST(IPQ))
                     IACBJ(IPQ) = IACBT(JLIST(IPQ))
                     ISCLJJ(IPQ) = ISCLJT(JLIST(IPQ))
                     ISCCJ(IPQ) = ISCCT(JLIST(IPQ))
                     LLPERJ(IPQ) = LLPERT(JLIST(IPQ))
                  ENDDO
               ENDIF
C     
C               IF (LPION) THEN
C                  DO IPQ = 1, JLEN
C                     IPICJ(IPQ) = IPICT(JLIST(IPQ))
C                  ENDDO
C               ENDIF
C     
C               IF (LDOVIR) THEN
C                  
C                  DO 3550 I41X = 1, JLEN
C                     NSPJ(I41X) = NSPTT(JLIST(I41X))
C 3550             CONTINUE
C                  IVIRL=0
C                  DO 3600 IPQ = 1, JLEN
C                     IF(NSPJ(IPQ) .NE. NSPTT(ILAST)) THEN
C                        IVIRL=IVIRL+1
C                        IOCVIR(IVIRL)=IPQ
C                     ENDIF
C 3600             CONTINUE
C                  LVIR=.FALSE.
C                  IF(IVIRL .GT. 0) THEN
C                     LVIR=.TRUE.
C                     
C                     DO 3700 IPQ = 1, JLEN
C                        XRJ(IPQ) = XRC(JLIST(IPQ))
C                        YRJ(IPQ) = YRC(JLIST(IPQ))
C                        ZRJ(IPQ) = ZRC(JLIST(IPQ))
C 3700                CONTINUE
C                  ENDIF
C               ENDIF
C     
C*****PROCES EACH I-ATOM, AND TREAT ALL J-ATOMS IN VECTORMODE.
C     
               DO 4350 I = IFIRST(NI), ILAST

                  ICE = ICET(I)
C     
C     GET THE LJ PARAMETERS
C     
                  IF (LPERTL) THEN
                     LIPERT = LLPERT(I)
!$OMP parallel do private(IPQ,NIJ,NIJB) shared(I)
                     DO IPQ = 1, JLEN
                        NIJ = JRCPAC(IACAJ(IPQ) + (IACT(I) - 1) *
     $                       NRATT) + JRCBK2(IPQ)
                        NIJB = JRCPAC(IACBJ(IPQ) + (IACBT(I) - 1) *
     $                       NRATT) + JRCBK2(IPQ)

                        VC6A(IPQ) = FORC6(NIJ)
                        VC12A(IPQ) = FORC12(NIJ)
                        VC6B(IPQ) = FORC6(NIJB)
                        VC12B(IPQ) = FORC12(NIJB)
                        AC12C6(IPQ) = FC12C6(NIJ)
                        BC12C6(IPQ) = FC12C6(NIJB)
                     ENDDO
                  ELSE
!$OMP parallel do private(IPQ,NIJ) shared(I)
                     DO IPQ = 1, JLEN
                        NIJ = JRCPAC(IACAJ(IPQ) + (IACT(I) - 1) *
     $                       NRATT) + JRCBK2(IPQ)
                        VC6A(IPQ) = FORC6(NIJ)
                        VC12A(IPQ) = FORC12(NIJ)
                     ENDDO
                  ENDIF

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
C*****CALCULATE THE MINIMAL DISTANCE IN ACCORD WITH THE PERIODIC BOX.
C     
!$OMP parallel do private(IPQ,XIJ,YIJ,ZIJ,WIJ)
#ifdef _vpp_
!OCL NOALIAS
#endif
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
                     R2IJ3D(IPQ) = VXIJ(IPQ)*VXIJ(IPQ) 
     $                           + VYIJ(IPQ)*VYIJ(IPQ)
     $                 +VZIJ(IPQ)*(VZIJ(IPQ)+COSB2*VXIJ(IPQ))

                     VRIJ2(IPQ) = R2IJ3D(IPQ)+VWIJ(IPQ)*VWIJ(IPQ)
                  ENDDO
                  IF (LOCTO) THEN
!$OMP parallel do private(IPQ,TEMP1)                     
                     DO IPQ = 1, JLEN
                        TEMP1 = BOXOQ - DABS(VXIJ(IPQ)) -
     .                      DABS(VYIJ(IPQ)) - DABS(VZIJ(IPQ))
                        IF (TEMP1 .LE. 0.0D0) THEN
                           VRIJ2(IPQ) = VRIJ2(IPQ) + TEMP1*BOX1
                           R2IJ3D(IPQ) =  R2IJ3D(IPQ) + TEMP1*BOX1
                           VXIJ(IPQ) = VXIJ(IPQ) -
     $                          SIGN(BOXOH,VXIJ(IPQ))
                           VYIJ(IPQ) = VYIJ(IPQ) -
     $                          SIGN(BOXOH,VYIJ(IPQ))
                           VZIJ(IPQ) = VZIJ(IPQ) -
     $                          SIGN(BOXOH,VZIJ(IPQ))
                        ENDIF
                     ENDDO
                  ENDIF
C     
!$OMP parallel do private(IPQ)
                  DO IPQ = 1, JLEN
                     VRIJ(IPQ) = DSQRT(VRIJ2(IPQ))
                     VRINV(IPQ) = 1.0D0/VRIJ(IPQ)
C     
                     VRINV2(IPQ) = VRINV(IPQ)*VRINV(IPQ)
                     VRINV6(IPQ) = VRINV2(IPQ)*VRINV2(IPQ)*
     .                             VRINV2(IPQ)
                  ENDDO
                  
                  IF (LPERTL) THEN
                     IF (LIPERT) THEN
                        DO IPQ = 1, JLEN
C     
C---- A SOFT-CORE PERTURBATION
C     CALCULATION IS DONE
C     GET STATE B PARAMETERS AND USE COMBINATION RULE FOR SOFT-CORE
C     ALPHA PARAMETERS
C     
                           QIQJA = CGAT(I) * CGAJ(IPQ)
                           QIQJB = CGBT(I) * CGBJ(IPQ)
C     
                           CA12 = VC12A(IPQ)
                           CA6 = VC6A(IPQ)
                           CB12 = VC12B(IPQ)
                           CB6  = VC6B(IPQ)
C     
                           SOFTLJ = ALPHLJ * 
     $                          (1 - ISCLJT(I) * ISCLJJ(IPQ))
                           SOFTEL = ALPHC * 
     $                          (1 - ISCCT(I) * ISCCJ(IPQ))
C     
                           CA126 = AC12C6(IPQ)
                           CB126 = BC12C6(IPQ)

                           RIJ3D2 = R2IJ3D(IPQ)
                           RIJ2 = VRIJ2(IPQ)

C*******************************************************************************
C     4D
C*******************************************************************************
C---- SOFT-CORE PERTURBATION
C     CALCULATE FULL 4D CONTRIBUTION
C     
                           RIJ4 = VRIJ2(IPQ) * VRIJ2(IPQ)
                           RIJ6 = RIJ4 * VRIJ2(IPQ)
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
                           VEL4D  =   RLB  * VELB4D + RLA  * VELA4D
                           VRF4D  =   RLB  * VRFB4D + RLA  * VRFA4D
                           VLJ4D  =   RLB  * VLJB4D + RLA  * VLJA4D
                           DF4D4D =   RLB  * DFB4D  + RLA  * DFA4D
                           DF4D   =   DF4D4D 
     $                              + RLB  * DFRFB4 + RLA  * DFRFA4
                           DLLJ4D =   RDLB * VLJB4D + RLLB * DLLJB4
     $                              + RDLA * VLJA4D + RLLA * DLLJA4
                           DLEL4D =   RDLB * VELB4D + RLLB * DLELB4
     $                              + RDLA * VELA4D + RLLA * DLELA4
                           DLRF4D =   RDLB * VRFB4D + RLLB * DLRFB4
     $                              + RDLA * VRFA4D + RLLA * DLRFA4


C*******************************************************************************
C  3D
C*******************************************************************************
C     
C---- SOFT-CORE PERTURBATION
C     CALCULATE FULL 3D CONTRIBUTION
C     
                           RIJ4 = R2IJ3D(IPQ) * R2IJ3D(IPQ)
                           RIJ6 = RIJ4 * R2IJ3D(IPQ)
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
C--   COMBINE 3D STATE A AND B
C     
                           VEL3D  =  RLB  * VELB3D + RLA  * VELA3D
                           VRF3D  =  RLB  * VRFB3D + RLA  * VRFA3D
                           VLJ3D  =  RLB  * VLJB3D + RLA  * VLJA3D
                           DF3D   =  RLB  * (DFB3D + DFRFB3)
     $                             + RLA  * (DFA3D + DFRFA3)
                           DLLJ3D =  RDLB * VLJB3D + RLLB * DLLJB3
     $                             + RDLA * VLJA3D + RLLA * DLLJA3
                           DLEL3D =  RDLB * VELB3D + RLLB * DLELB3 
     $                             + RDLA * VELA3D + RLLA * DLELA3
                           DLRF3D =  RDLB * VRFB3D + RLLB * DLRFB3
     $                             + RDLA * VRFA3D + RLLA * DLRFA3
                           
                           DDM3EL =  RLBRMU * DLELB3 + RLARMU * DLELA3
                           DDM3RF =  RLBRMU * DLRFB3 + RLARMU * DLRFA3
                           DDM3LJ =  RLBRMU * DLLJB3 + RLARMU * DLLJA3
C     
C---- DONE SOFT-CORE 3D (FULL CONTRIBUTION)
C     
C     
C---- THE DISTANCE INDEPENDANT REACTION FIELD COMPONENT IS
C     DIMENSIONALITY INDEPENDENT. CALCULATE IT HERE.
C     
                           VRCA = RFC*QIQJA
                           VRCB = RFC*QIQJB
                           VRC  = RLB*VRCB  + RLA*VRCA
                           DLRC = RDLB*VRCB + RDLA*VRCA
C     
                           DF3(IPQ)   = RM3D*DF3D + RM4D*DF4D
                           DF4TH(IPQ) = RM4D*DF4D4D

                           VEELJ(IPQ) = RM3D*VEL3D + RM4D*VEL4D
                           VERFJ(IPQ) = RM3D*VRF3D + RM4D*VRF4D
                           VERCJ(IPQ) = RM3D*VRC   + RM4D*VRC
                           VLJJ(IPQ)  = RM3D*VLJ3D + RM4D*VLJ4D
C
                           EGLELJ(IPQ) = RM3D*DLEL3D + RM4D*DLEL4D
                           EGLRFJ(IPQ) = RM3D*DLRF3D + RM4D*DLRF4D
                           EGLRCJ(IPQ) = RM3D*DLRC   + RM4D*DLRC
                           EGLLJJ(IPQ) = RM3D*DLLJ3D + RM4D*DLLJ4D
C
                           E34ELJ(IPQ) = 
     $                          RDM3D*VEL3D + RM3D*DDM3EL + RDM4D*VEL4D
                           E34RFJ(IPQ) =
     $                          RDM3D*VRF3D + RM3D*DDM3RF + RDM4D*VRF4D
                           E34RCJ(IPQ) =
     $                          RDM3D*VRC + RDM4D*VRC
                           E34LJJ(IPQ) =
     $                          RDM3D*VLJ3D + RM3D*DDM3LJ + RDM4D*VLJ4D
                           
                        ENDDO
                     ELSE
                        DO IPQ = 1, JLEN
                           IF (LLPERJ(IPQ)) THEN
C     
C---- A SOFT-CORE PERTURBATION
C     CALCULATION IS DONE
C     GET STATE B PARAMETERS AND USE COMBINATION RULE FOR SOFT-CORE
C     ALPHA PARAMETERS
C     
                           QIQJA = CGAT(I) * CGAJ(IPQ)
                           QIQJB = CGBT(I) * CGBJ(IPQ)
C     
                           CA12 = VC12A(IPQ)
                           CA6 = VC6A(IPQ)
                           CB12 = VC12B(IPQ)
                           CB6  = VC6B(IPQ)
C     
                           SOFTLJ = ALPHLJ*(1-ISCLJT(I)*ISCLJJ(IPQ))
                           SOFTEL = ALPHC *(1-ISCCT(I) *ISCCJ(IPQ)) 
C     
                           CA126 = AC12C6(IPQ)
                           CB126 = BC12C6(IPQ)

                           RIJ3D2 = R2IJ3D(IPQ)
                           RIJ2 = VRIJ2(IPQ)

C*******************************************************************************
C     4D
C*******************************************************************************
C---- SOFT-CORE PERTURBATION
C     CALCULATE FULL 4D CONTRIBUTION
C     
                           RIJ4 = VRIJ2(IPQ) * VRIJ2(IPQ)
                           RIJ6 = RIJ4 * VRIJ2(IPQ)
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
                           VEL4D  =   RLB  * VELB4D + RLA  * VELA4D
                           VRF4D  =   RLB  * VRFB4D + RLA  * VRFA4D
                           VLJ4D  =   RLB  * VLJB4D + RLA  * VLJA4D
                           DF4D4D =   RLB  * DFB4D  + RLA  * DFA4D
                           DF4D   =   DF4D4D 
     $                              + RLB  * DFRFB4 + RLA  * DFRFA4
                           DLLJ4D =   RDLB * VLJB4D + RLLB * DLLJB4
     $                              + RDLA * VLJA4D + RLLA * DLLJA4
                           DLEL4D =   RDLB * VELB4D + RLLB * DLELB4
     $                              + RDLA * VELA4D + RLLA * DLELA4
                           DLRF4D =   RDLB * VRFB4D + RLLB * DLRFB4
     $                              + RDLA * VRFA4D + RLLA * DLRFA4

C*******************************************************************************
C  3D
C*******************************************************************************
C     
C---- SOFT-CORE PERTURBATION
C     CALCULATE FULL 3D CONTRIBUTION
C     
                           RIJ4 = R2IJ3D(IPQ) * R2IJ3D(IPQ)
                           RIJ6 = RIJ4 * R2IJ3D(IPQ)
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
C--   COMBINE 3D STATE A AND B
C     
                           VEL3D  =  RLB  * VELB3D + RLA  * VELA3D
                           VRF3D  =  RLB  * VRFB3D + RLA  * VRFA3D
                           VLJ3D  =  RLB  * VLJB3D + RLA  * VLJA3D
                           DF3D   =  RLB  * (DFB3D + DFRFB3)
     $                             + RLA  * (DFA3D + DFRFA3)
                           DLLJ3D =  RDLB * VLJB3D + RLLB * DLLJB3
     $                             + RDLA * VLJA3D + RLLA * DLLJA3
                           DLEL3D =  RDLB * VELB3D + RLLB * DLELB3 
     $                             + RDLA * VELA3D + RLLA * DLELA3
                           DLRF3D =  RDLB * VRFB3D + RLLB * DLRFB3
     $                             + RDLA * VRFA3D + RLLA * DLRFA3

                           DDM3EL =  RLBRMU * DLELB3 + RLARMU * DLELA3
                           DDM3RF =  RLBRMU * DLRFB3 + RLARMU * DLRFA3
                           DDM3LJ =  RLBRMU * DLLJB3 + RLARMU * DLLJA3
C     
C---- DONE SOFT-CORE 3D (FULL CONTRIBUTION)
C     
C     
C---- THE DISTANCE INDEPENDANT REACTION FIELD COMPONENT IS
C     DIMENSIONALITY INDEPENDENT. CALCULATE IT HERE.
C     
                           VRCA = RFC*QIQJA
                           VRCB = RFC*QIQJB
                           VRC  = RLB*VRCB  + RLA*VRCA
                           DLRC = RDLB*VRCB + RDLA*VRCA
                           
                           DF3(IPQ)    = RM3D * DF3D + RM4D * DF4D
                           DF4TH(IPQ) = RM4D * DF4D4D
C     
                           VEELJ(IPQ) = 
     $                          RM3D * VEL3D + RM4D * VEL4D
                           VERFJ(IPQ) = 
     $                          RM3D * VRF3D + RM4D * VRF4D
                           VERCJ(IPQ) =
     $                          RM3D * VRC + RM4D * VRC
                           VLJJ(IPQ) = 
     $                          RM3D * VLJ3D + RM4D * VLJ4D
C
                           EGLELJ(IPQ) = 
     $                          RM3D * DLEL3D + RM4D * DLEL4D
                           EGLRFJ(IPQ) = 
     $                          RM3D * DLRF3D + RM4D * DLRF4D
                           EGLRCJ(IPQ) = 
     $                          RM3D * DLRC + RM4D * DLRC
                           EGLLJJ(IPQ) = 
     $                          RM3D * DLLJ3D + RM4D * DLLJ4D
C
                           E34ELJ(IPQ) = 
     $                          RDM3D*VEL3D + RM3D*DDM3EL + RDM4D*VEL4D
                           E34RFJ(IPQ) =
     $                          RDM3D*VRF3D + RM3D*DDM3RF + RDM4D*VRF4D
                           E34RCJ(IPQ) =
     $                          RDM3D*VRC + RDM4D*VRC
                           E34LJJ(IPQ) =
     $                          RDM3D*VLJ3D + RM3D*DDM3LJ + RDM4D*VLJ4D
                           
                           ELSE
                              QIQJA = CGAT(I)*CGAJ(IPQ)
                              VERFJ(IPQ) = QIQJA * RFE * R2IJ3D(IPQ)
                              VERCJ(IPQ) = QIQJA * RFC
                              VEELJ(IPQ) = QIQJA*VRINV(IPQ)
                              VLJJ(IPQ) = VRINV6(IPQ)*(VC12A(IPQ)
     $                             *VRINV6(IPQ)-VC6A(IPQ))

                              EGLELJ(IPQ) = 0.0D0
                              EGLRFJ(IPQ) = 0.0D0
                              EGLRCJ(IPQ) = 0.0D0
                              EGLLJJ(IPQ) = 0.0D0

                              E34ELJ(IPQ) = 0.0D0
                              E34RFJ(IPQ) = 0.0D0
                              E34RCJ(IPQ) = 0.0D0
                              E34LJJ(IPQ) = 0.0D0
                           
                              DF = VRINV6(IPQ)*(2.0D0*VRINV6(IPQ) *
     $                             VC12A(IPQ)-VC6A(IPQ))
                              DF = (6.0D0*DF+VEELJ(IPQ))*VRINV2(IPQ)
                              DF3(IPQ) = DF +  QIQJA * RFF
                              DF4TH(IPQ) = DF
                           ENDIF

                        ENDDO
                     
                     ENDIF
!$OMP parallel do private(IPQ) reduction(+:ELGLEL,ELGLRF,ELGLRC)
!$OMP+ reduction(+:ELGLLJ,EL34EL,EL34RF,EL34RC,EL34LJ)
                     DO  IPQ = 1, JLEN
                        ELGLEL = ELGLEL + EGLELJ(IPQ)
                        ELGLRF = ELGLRF + EGLRFJ(IPQ)
                        ELGLRC = ELGLRC + EGLRCJ(IPQ)
                        ELGLLJ = ELGLLJ + EGLLJJ(IPQ)

                        EL34EL = EL34EL + E34ELJ(IPQ)
                        EL34RF = EL34RF + E34RFJ(IPQ)
                        EL34RC = EL34RC + E34RCJ(IPQ)
                        EL34LJ = EL34LJ + E34LJJ(IPQ)
                     ENDDO
                     
                  ELSE
                     
C                     IF (LPION) THEN
C                        DO IPQ = 1, JLEN
C                           
C                           QIQJA = CGAT(I)*CGAJ(IPQ) * PIFACT(IPQ)
C                           VERFJ(IPQ) = QIQJA * RFE * VRIJ2(IPQ)
C                           VERCJ(IPQ) = QIQJA * RFC
C                           VEELJ(IPQ) = QIQJA*VRINV(IPQ)                    
C                           VLJJ(IPQ) = VRINV6(IPQ)*(VC12A(IPQ)*VRINV6(
C     $                          IPQ)-VC6A(IPQ)) * PIFACT(IPQ)
C                           
C                           DF = VRINV6(IPQ)*(2.0E0*VRINV6(IPQ)*
C     $                          VC12A(IPQ)-VC6A(IPQ)) * PIFACT(IPQ)
C                           DF3(IPQ) = (6.0E0 * DF + VEELJ(IPQ))*
C     $                          VRINV2(IPQ) +  QIQJA * RFF
C                        ENDDO

C                     ELSE
!$OMP parallel do private(IPQ,QIQJA,DF) shared(I)
                     DO IPQ = 1, JLEN
                        QIQJA = CGAT(I)*CGAJ(IPQ)
                        VERFJ(IPQ) = QIQJA * RFE * R2IJ3D(IPQ)
                        VERCJ(IPQ) = QIQJA * RFC
                        VEELJ(IPQ) = QIQJA*VRINV(IPQ)
                        VLJJ(IPQ) = VRINV6(IPQ)*(VC12A(IPQ)*
     $                              VRINV6(IPQ)-VC6A(IPQ))
                        DF = (6.0D0*(
     .                       VRINV6(IPQ)*(2.0D0*VRINV6(IPQ)*
     .                       VC12A(IPQ)-VC6A(IPQ))
     .                       )+VEELJ(IPQ))*VRINV2(IPQ) 
                        DF3(IPQ) = DF +  QIQJA * RFF
                        DF4TH(IPQ) = DF
                     ENDDO   
C                  ENDIF
                     
                  ENDIF
C     
                  JRCTOT=JLEN
                  JRCSRT=0
 4590             CONTINUE
                  JRCLEN=MIN(IVECT,JRCTOT)

CDIR$ SHORTLOOP                  
                  DO JRC=1,JRCLEN

                     TEEL(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TEEL(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VEELJ(JRC+JRCSRT)
                     
                     TERF(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TERF(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VERFJ(JRC+JRCSRT)
                     
                     TERC(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TERC(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VERCJ(JRC+JRCSRT)
                     
                     TELJ(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) =
     $                    TELJ(NRELKP(JRCILJ(JRC+JRCSRT),ICE),JRC) +
     $                    VLJJ(JRC+JRCSRT)
                     
                  ENDDO

                  JRCSRT=JRCSRT+JRCLEN
                  JRCTOT=JRCTOT-JRCLEN
                  IF(JRCTOT.GT.0) GOTO 4590
C     
C*****CALCULATE AND STORE THE FORCES.
C     
!$OMP parallel do private(IPQ)
#ifdef _vpp_
!OCL NOALIAS
#endif
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
C                     IF (LVIR) THEN
C                        
C                        DO 4300 IPQ=1,IVIRL
C                           TVIR(IPQ,1) = TVIR(IPQ,1) -
C     $                          VXH(IOCVIR(IPQ))*
C     $                          (VXIJ(IOCVIR(IPQ)) +
C     $                          XRJ(IOCVIR(IPQ))-XRC(I))
C                           TVIR(IPQ,2) = TVIR(IPQ,2) -
C     $                          VYH(IOCVIR(IPQ))*
C     $                          (VYIJ(IOCVIR(IPQ)) +
C     $                          YRJ(IOCVIR(IPQ))-YRC(I))
C                           TVIR(IPQ,3) = TVIR(IPQ,3) -
C     $                          VZH(IOCVIR(IPQ))*
C     $                          (VZIJ(IOCVIR(IPQ)) +
C     $                          ZRJ(IOCVIR(IPQ))-ZRC(I))
C 4300                   CONTINUE
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
 4350          CONTINUE
C     
!$OMP parallel do private(IPQ)
               DO IPQ = 1, JLEN
                  FLRX(JLIST(IPQ)) = FXTEMP(IPQ)
                  FLRY(JLIST(IPQ)) = FYTEMP(IPQ)
                  FLRZ(JLIST(IPQ)) = FZTEMP(IPQ)
                  FLRW(JLIST(IPQ)) = FWTEMP(IPQ)
               ENDDO
            ENDIF
         ENDIF
 4450 CONTINUE
C     
      INB(NRAGT) = 0
C--   UPDATE THE 'ELREL' VALUES
      DO JRC1=1,NUNRE2
         DO JRC2=1,IVECT
            ELREL(JRC1) = ELREL(JRC1) + TEEL(JRC1, JRC2)
            ELRRF(JRC1) = ELRRF(JRC1) + TERF(JRC1, JRC2)
            ELRRC(JRC1) = ELRRC(JRC1) + TERC(JRC1, JRC2)
            ELRLJ(JRC1) = ELRLJ(JRC1) + TELJ(JRC1, JRC2)
         ENDDO
      ENDDO
C     
!$OMP parallel do private(JRC)
      DO JRC = 1, NR
         FLR((JRC-1)*4+1) = FLRX(JRC)
         FLR((JRC-1)*4+2) = FLRY(JRC)
         FLR((JRC-1)*4+3) = FLRZ(JRC)
         FLR((JRC-1)*4+4) = FLRW(JRC)
      ENDDO
C     
      IF (LDOTRA) CALL TRACO (NR, 0, FLR, BETA, -1, LEVERY)
C     
C      IF (LDOVIR) THEN
C         
C         DO 6460 JRC=1,MAXLEN
C            VIRLR(1)=VIRLR(1)+TVIR(JRC,1)
C            VIRLR(2)=VIRLR(2)+TVIR(JRC,2)
C            VIRLR(3)=VIRLR(3)+TVIR(JRC,3)
C 6460    CONTINUE
C         VIRLR(1) = VIRLR(1)*0.5E0
C         VIRLR(2) = VIRLR(2)*0.5E0
C         VIRLR(3) = VIRLR(3)*0.5E0
C      ENDIF
C     
      IF (LDOTRA) THEN
         CALL TRACO (NR, 0, XCOORD, BETA, -1, LEVERY)
      ENDIF
C     
C     
C----------RETURN DERIVATIVE VARIABLES ACCORDING TO REQUEST
C     
      IF (.NOT.L3D4D) THEN
C     
C--   IN THIS CASE THESE WOULD BE *WRONG*
C     
        EL34EL = 0.0D0
        EL34RF = 0.0D0
        EL34RC = 0.0D0
        EL34LJ = 0.0D0
      ENDIF
      IF (.NOT.LDOPER) THEN
C     
C--   IN THIS CASE THESE WOULD BE OK, BUT SINCE BONDED ROUTINES
C     WON'T RETURN THEM, RETURNING THEM MIGHT GENERATE CONFUSION... 
C     
        ELGLEL = 0.0D0
        ELGLRF = 0.0D0
        ELGLRC = 0.0D0
        ELGLLJ = 0.0D0
      ENDIF

      PRINT 9040, NJNB      
 9040 FORMAT (1X,'NBPML:',I10,' ELEMENTS IN THE PAIRLIST')

      RETURN
C     
C     
C*****ERROR EXIT
C     
C     
 9100 FORMAT (A, I7)
 4550 CONTINUE
      PRINT *,PRGSTR,' PAIRLIST OVERFLOW:'
      PRINT *,'MAXIMUM SIZE: MAXJNB =',MAXJNB
      PRINT *,'FOUND: ',NJNB
C     
C     
C*****DUMP FROM SUBROUTINE.
C     
C     
 4600 CONTINUE
C     
 9120 FORMAT (10(1X,I9))
 9140 FORMAT (10(1X,G11.4))
 9160 FORMAT (/,1X,20A)
 9180 FORMAT (1X,A,'(', I2, '):', 3(5X, 3(F8.3, 2X)))
 9200 FORMAT (1X, 10L10)
C     
 5750 CONTINUE
      PRINT '(A)', '1THE SUBROUTINE NBPML HAS GONE WRONG. ',
     1     ' THE VARIABLES ARE:'
      PRINT 9160, '      NPM       NRP      NCAG       NSM',
     1     '      NRAM       NTB      NSPM      NJNB'
      PRINT 9120, NPM, NRP, NCAG, NSM, NRAM, NTB, NSPM, NJNB
      PRINT 9160, '     BOX(1)      BOX(2)      BOX(3)',
     1     '        BETA       RCUTP       RCUTL'
      PRINT 9140, (BOX(I), I = 1, 3), BETA,DSQRT(RCUTP2), 
     $     DSQRT(RCUTL2)
C     PRINT 9160, '   NRE(1)    NRE(2)    NRE(3)    NRE(4)'
C     PRINT 9120, (NRE(K), K = 1, 4)
      PRINT 9160, '   ELREL(1..10)'
      PRINT 9140, (ELREL(I), I = 1, 10)
      PRINT 9160, '   VIRLR(1..3)'
      PRINT 9140, (VIRLR(I), I = 1, 3)
      PRINT 9160, '    NSP(1..NSPM)'
      PRINT 9120, (NSP(I), I = 1, NSPM)
      PRINT 9160, '    INC(1..20)'
      PRINT 9120, (INC(I), I = 1, 20)
      PRINT 9160, '     CG(1..20)'
      PRINT 9140, (CG(I), I = 1, 20)
      PRINT 9160, '    CGS(1..NRAM)'
      PRINT 9140, (CGS(I), I = 1, NRAM)
      PRINT 9160, '      X(1..45)'
      DO 4650 I = 1, 36, 9
         PRINT 9180, 'XCOORD', I, (XCOORD(M), M = I, I + 8)
 4650 CONTINUE
      PRINT 9160, '    FLR(1..30)'
      PRINT 9140, (FLR(I), I = 1, 30)
      PRINT 9160, '     XR(1..10)'
      PRINT 9140, (XR(I), I = 1, 10)
      PRINT 9160, '    INB(1..20)'
      PRINT 9120, (INB(I), I = 1, 20)
      PRINT 9160, '    JNB(1..200)'
      PRINT 9120, (JNB(I), I = 1, 200)
      PRINT 9160, '                  LOCAL VARIABLES:'
      PRINT 9160, '     LERROR',
     1     '    LDOTRA'
      PRINT 9200, LERROR, LDOTRA
      PRINT 9160, '       ZERO      RCUTP2      RCUTL2       COSB2'
      PRINT 9140, RCUTP2, RCUTL2, COSB2
      PRINT 9160, '     NRAG     NRAGT ',
     1     '        NR'
      PRINT 9120, NRAG, NRAGT, NR
      PRINT 9160, '        I       ICE        NI',
     1     '        NN     NNOLD', '        NS       NAG'
C     PRINT 9120, I, ICE, NI, NN, NNOLD, NS, NAG
      PRINT 9160, '        I     LRLEN    LRLEN2      JLEN'
      PRINT 9120, I, LRLEN, LRLEN2, JLEN
      PRINT 9160, 'BOXH(1..3), BOXINV(1..3)'
      PRINT 9140, (BOXH(I), I = 1, 3), (BOXINV(I), I = 1, 3)
      CALL STOPGM('NBPML','FATAL ERROR')

C     END NBWIT4
C     ==--------------------------------------------------------------==
      RETURN
      END
C MIMICS SOME CRAY INTRINSIC ROUTINES. NOT FULLY TESTED.
C I.G. TIRONI MAY 1996
C
C     ==================================================================
      SUBROUTINE WHENFLE (LEN, VRIJ2, INCX, RCUTP2, JRCIND, LRLEN)
C     ==--------------------------------------------------------------==
COMMSUBR WHENFLE
C     SUBROUTINE WHENFLE (LEN, VRIJ2, INCX, RCUTP2, JRCIND, LRLEN)
C      
C     MIMICS THE CRAY WHENFLE ROUTINE.
C
COMMEND
C     ARGS
      INTEGER LEN, INCX, LRLEN, JRCIND(LEN)
      real*8 RCUTP2, VRIJ2(*)

C     LOCAL
      INTEGER I, ICOUNT

C     BEGIN WHENFLE
      IF (INCX .EQ. 0) THEN
         RETURN
      ENDIF
      ICOUNT = 0
      IF (INCX .GT. 0) THEN

         DO 10 I = 1, (LEN - 1)*INCX + 1, INCX
            IF (VRIJ2(I) .LE. RCUTP2) THEN
               ICOUNT = ICOUNT + 1
               JRCIND(ICOUNT) = I
            ENDIF
 10      CONTINUE

      ELSE
         
         DO 20 I = LEN - 1, 0, -1
              IF (VRIJ2(1 - INCX*I) .LE. RCUTP2) THEN
               ICOUNT = ICOUNT + 1
               JRCIND(ICOUNT) = 1 - INCX*I
            ENDIF
 20      CONTINUE

      ENDIF

      LRLEN = ICOUNT
C     END WHENFLE
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SCOPY(N, X, INCX, Y, INCY)
C     ==--------------------------------------------------------------==
COMMSUBR SCOPY
C      SUBROUTINE SCOPY(N, X, INCX, Y, INCY)
C      
C     MIMICS THE CRAY SCOPY ROUTINE
C
COMMEND
C     ARGS
      INTEGER N, INCX, INCY
      real*8 X(*), Y(*)

C     LOCAL
      INTEGER I

C     BEGIN

      IF (N .LE. 0 .OR. INCX .EQ. 0 .OR. INCY .EQ. 0) THEN
         RETURN
      ENDIF

      IF (INCX .GT. 0) THEN

         IF (INCY .GT. 0) THEN
            DO I = 1, N
               Y((I - 1)*INCY + 1) = X((I - 1)*INCX + 1)
            ENDDO
         ELSE
            DO I = N - 1, 0 , -1
               Y(1 - I * INCY) = X((N - I - 1) * INCX + 1)
            ENDDO
         ENDIF
         
      ELSE

         IF (INCY .GT. 0) THEN
            DO I =  N - 1, 0 , -1
               Y((N - I - 1) * INCX + 1) = X(1 - I * INCX)
            ENDDO
         ELSE
            DO I = N - 1, 0 , -1
               Y(1 - I * INCY) = X(1 - I * INCX)
            ENDDO
         ENDIF 
      ENDIF

C     END SCOPY
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      FUNCTION ISRCHEQ(N, X, INCX, TARGET)
C     ==--------------------------------------------------------------==
COMMSUBR ISRCHEQ
C     FUNCTION ISRCHEQ(N, X, INCX, TARGET)
C      
C     MIMICS THE CRAY ISRCHEQ ROUTINE.
C
COMMEND
C     ARGS
      INTEGER N, X(*), INCX, TARGET, ISRCHEQ

C     LOCAL
      INTEGER I, J

C     BEGIN ISRCHEQ
      J = 1
      ISRCHEQ = 0
      IF (N .LE. 0) THEN
         RETURN
      ENDIF

      IF (INCX .LT. 0) THEN
         J = 1 - (N - 1) * INCX
      ENDIF

      DO I = 1, N
         IF (X(J) .EQ. TARGET) GOTO 200
         J = J + INCX
      ENDDO

 200  ISRCHEQ = I
      RETURN
    
C     END ISRCHEQ
      END
C     ==================================================================
      SUBROUTINE WHENIGT(N, X, INCX, ITARGET, INDEX, NN)
C     ==--------------------------------------------------------------==
COMMSUBR SCOPY
C     SUBROUTINE WHENIGT(N, X, INCX, ITARGET, INDEX, NN)
C      
C     MIMICS THE CRAY WHENIGT ROUTINE
C
COMMEND
C
C     ARGS
      INTEGER N, X(*), INCX, ITARGET, INDEX(N), NN

C     LOCAL
      INTEGER I, J

C     BEGIN WHENIGT

      J = 1
      NN = 0
      
      IF (N .EQ. 0) THEN
         RETURN
      ENDIF

      IF (INCX .LT. 0) THEN
         J = 1 - (N - 1) * INCX
      ENDIF

      DO I = 1, N
         IF (X(J) .GT. ITARGET) THEN
            NN = NN + 1
            INDEX(NN) = J
         ENDIF
         J = J + INCX
      ENDDO

C     END WHENIGT
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      FUNCTION ISRCHFLE(N, X, INCX, FTARGET)
C     ==--------------------------------------------------------------==
COMMSUBR ISRCHFLE
C     SUBROUTINE ISRCHFLE(N, X, INCX, FTARGET)
C      
C     MIMICS THE CRAY ISRCHFLE ROUTINE
C
COMMEND
C
C     ARGS
      INTEGER N, INCX, ISRCHFLE
      real*8 X(*), FTARGET

C     LOCAL
      INTEGER I, J

C     BEGIN ISRCHFLE

      J = 1
      ISRCHFLE = 0
      IF (N .EQ. 0) THEN
         RETURN
      ENDIF
      
      IF (INCX .LT. 0) THEN
         J = 1 - (N - 1) * INCX
      ENDIF

      DO I = 1, N
         IF (X(J) .LE. FTARGET) GOTO 200
         J = J + INCX
      ENDDO

 200  ISRCHFLE = I
      RETURN
    
C     END ISRCHEQ
      END
       
#endif
