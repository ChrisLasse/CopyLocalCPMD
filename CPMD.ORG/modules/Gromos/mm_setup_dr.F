C     ==--------------------------------------------------------------==
      SUBROUTINE mm_setup_dr(NRPTp,NSOLVp,NATTOTp,NUM4DAp,TIMEp,
C     ==--------------------------------------------------------------==
Cmb - Revised on 24 March 2006
C IN   : NRPT,NSOLV,NATTOT,NUM4DA,IUTRJX,IUTRJV,IUTRJE,
C        IUTRJG,IFNCFG
C OUT  :
C INOUT: TIME
C
C     RUNMD will carry out a molecular dynamics run for a member
C     of solute (protein) molecules plus solvent molecules.
C     Periodic boundary conditions can be applied.
C
C     It starts at time T, runs with time step L<DT> and stops when
C     the number of steps reaches L<NSTLIM>.
C
C     If L<NTT> > 0, the temperature T of the system is coupled to 
C     that of a bath of temperature L<TEMP0> with coupling time
C     constant L<TAUTP> or L<TAUTS>.
C
C     If NTP >0, the pressure P of the system is coupled to that of
C     a bath of pressure PRES0 with coupling time constant TAUP
C     (see Berendsen et al., J.Chem.Phys 81 (1984) 3684).
C
C     The differential equations are integrated using a leap-frog
C     scheme.
C     At each step, subroutine L<SHAKE> is called to satisfy actual
C     constraints, such as covalent bond-lengths or angles.
C
C     Subroutine L<FORCE> supplies the potential energies and forces.
C     Subroutine L<CENMAS> calculates the center of mass coordinates
C     and motion, when required.
C     Subroutine L<STOPCM> stops the translation of and rotation
C     about the centre of mass when required.
C     Subroutine L<LSQFIT> performs a least square fit of the total
C     energy to a straight line and prints the results.
C     The subroutines L<WRPOS>, L<WRVEL> and L<WRNRG> write
C     the position, velocity and energy trajectories, respectively.
C
C     NRPT   : number of solute (protein) atoms
C     NSOLV  : number of solvent atoms
C     NATTOT : total number of atoms
C     NUM4DA : number of atoms which are in 4D
C     TIME   : time of the system. This variable is updated
C              within this routine.
C
C     All of these files have been opened before calling this
C     routine.
C     IUTRJX : coordinate trajectory
C     IUTRJV : velocity trajectory file unit number
C     IUTRJE : energies trajectory file unit number
C     IUTRJG : free energies trajectory file
C     IFNCFG : final coordinates unit number
C     ==--------------------------------------------------------------==
     $     DLSUMp,DMSUMp,
     $     IUTRJXp,IUTRJVp,IUTRJEp,IUTRJGp,
#ifdef WRITE_FORCES
     $     IUTRJFp,
#endif
#ifdef READ_COORD
     $     IURDCOp,
#endif
     $     IFNCFGp)
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'md.h'
      INCLUDE 'box.h'
      INCLUDE 'formats.h'
      INCLUDE 'cenmas.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'forcear.h'
      INCLUDE 'runmd.h'
      INCLUDE 'pertsz.h'
      INCLUDE 'pertar.h'
C args
      INTEGER NRPTp,NSOLVp,NATTOTp,NUM4DAp
      real*8 TIMEp, DLSUMp(MXETBL),DMSUMp(MXETBL)
      INTEGER IUTRJXp,IUTRJVp,IUTRJEp,IUTRJGp,IFNCFGp
      INTEGER maxnrsp,maxnrpp,maxnatp,nsmp,iii,nrpp
#ifdef WRITE_FORCES
      integer     IUTRJFp
#endif
#ifdef READ_COORD
      integer     IURDCOp
#endif
C Qmmm: this common is used to save argumants
      INTEGER NRPT,NSOLV,NATTOT,NUM4DA
      real*8 TIME, DLSUM(MXETBL),DMSUM(MXETBL)
      INTEGER IUTRJX,IUTRJV,IUTRJE,IUTRJG,IFNCFG
      common/arg_i/NRPT,NSOLV,NATTOT,NUM4DA,IUTRJX,IUTRJV,IUTRJE,
     c             IUTRJG,IFNCFG
#ifdef WRITE_FORCES
     $     ,IUTRJF
#endif
#ifdef READ_COORD
     $     ,IURDCO
#endif
      common/arg_r/TIME,DLSUM,DMSUM
C local params
C used for the LS fitting array FIT
      INTEGER NFIT,NFITI,NFIT2,MAXFIT
      PARAMETER (NFIT = 1,NFITI = 2,NFIT2 = 3,MAXFIT = NFIT2)
C LPRKIN: whether we print out kinetic information or not.
C in MD and SD this is true, in EM it is false
      LOGICAL LPRKIN
      PARAMETER (LPRKIN  = .TRUE.)

C local vars
      INTEGER NRPT3,NSOLV3,NATTO3,I3,K,N,NITP,NITS
      INTEGER NCONP,NCONBC
      INTEGER NSTEP,II
      LOGICAL LMONO,LOCTO,LVAC,LDOVIR,LDOPER,LPR4
      LOGICAL LDOPIN
      LOGICAL LDODLA,LDOVDL,LDODMU,LDOVDM
C     LWR4: whether to write all four coords to trajectory or not.
      LOGICAL LWR4
      LOGICAL LREMCM,LPRLSQ
C whether to shake or not
      LOGICAL LSHKPR,LSHKSO
C
C TOTMAS = total mass of the whole system
C TMASS= total mass of one solvent molecule
C SUBMAS: the mass of the solute submolecules 1..NSP
C SUBMIN: the inverse mass of the submolecules 1..NSP
      real*8 TOTMAS,TMASS,SUBMAS(MAXNSP),SUBMIN(MAXNSP)
      real*8 DTHALF,BOHALF
      real*8 XCM(MAXDIM),VCM(MAXDIM),ACM(MAXDIM),OCM(MAXDIM)
      real*8 EKCM(MAXDIM),EKCMTO,EKROT
      real*8 XMIN(MAXDIM),TMP,EDUM
      DOUBLE PRECISION DTMP

C     NOTE: for all those arrays for which averages and r.m.s fluctuations
C     are collected and printed at the end of the run, five arrays
C     are kept. For example, for the temperature we have:
C
C     TEMP: contains the temperatures calculated at every step.
C     
C     TMPSUM, TMPSQ: contain the sum and the sum of the squares
C     of all TEMP values over the run. These variables are always
C     in DOUBLE PRECISION for better statistics.
C
C     TMPAVE,TMPFLC: after calling CLCAVE, these arrays contain
C     the average and mean square fluctuations of the TEMP
C     values over the run. These variables are always in declared
C     as REAL. 
C
C temperature arrays
C TFACPR: the factors to get from EKIN to TEMP for printing
      INTEGER NDOF(NFTMAX)
      real*8 TFACPR(NFTMAX),SCAL(NFTMAX)
      real*8 ESCL(NFTMAX),EKNOW(NFTMAX)
      real*8 TEMP(NFTMAX),TMPAVE(NFTMAX),TMPFLC(NFTMAX)
      DOUBLE PRECISION TMPSUM(NFTMAX),TMPSQ(NFTMAX)

C aux arrays for temp coupling
      INTEGER NBATH
      INTEGER NBNUM(MXBATH),NBNDX(MXBATH,MXBATH)
      real*8 EKREF(MXBATH),EKBATH(MXBATH),SCLBAT(MXBATH)
      real*8 DTBATH(MXBATH),TFACBT(MXBATH)
      real*8 EKTMP,SCLTMP

C used in scaling of 4th D reference temperature 
      real*8 EKIN4B,EKIN40,CW4D

C velocity of centres of mass of the solute submolecules
      real*8 VELCM(MXSBCO)

C pressure scaling vars
      real*8 DTCP

C used for LSQ fitting
      INTEGER IFSTEP
      real*8 FITAR(MAXFIT)

C virial, pressure stuff
      real*8 VIR(MAXDIM),PRES(MAXDIM),VIRTOT,PRESTO
      LOGICAL LERR,LEVERY,LFORM

C nonbonded atom group energy components
C force calculates these energy arrays
      real*8 EPLJ(MXNRE2),EPLJAV(MXNRE2),EPLJFL(MXNRE2)
      DOUBLE PRECISION EPLJSU(MXNRE2),EPLJSQ(MXNRE2)
C
      real*8 EPEL(MXNRE2),EPELAV(MXNRE2),EPELFL(MXNRE2)
      DOUBLE PRECISION EPELSU(MXNRE2),EPELSQ(MXNRE2)
C
      real*8 EPRF(MXNRE2),EPRFAV(MXNRE2),EPRFFL(MXNRE2)
      DOUBLE PRECISION EPRFSU(MXNRE2),EPRFSQ(MXNRE2)
C
      real*8 EPRC(MXNRE2),EPRCAV(MXNRE2),EPRCFL(MXNRE2)
      DOUBLE PRECISION EPRCSU(MXNRE2),EPRCSQ(MXNRE2)

C energy , dE/dLamda and dE/dmu
      real*8 ENER(MXETBL),EPAVE(MXETBL),EPFLC(MXETBL)
      DOUBLE PRECISION EPSUM(MXETBL),EPSQ(MXETBL)
C
      real*8 DEDLAM(MXETBL),DEDLAV(MXETBL),DEDLFL(MXETBL)
      DOUBLE PRECISION DEDLSU(MXETBL),DEDLSQ(MXETBL)
C
      real*8 DEDMU(MXETBL),DEDMAV(MXETBL),DEDMFL(MXETBL)
      DOUBLE PRECISION DEDMSU(MXETBL),DEDMSQ(MXETBL)


C (dE/dLamda)*DLAM and (dE/dMu)*DMU
C we do NOT calculate the averages of these values.
C instead, we add up their contributions to DLSUM and DMSUM respectively.
C DLSUM and DMSUM are passed into RUNMD as arguments and saved to the
C final configuration by PROMD.
      real*8 DVDL(MXETBL),DVDM(MXETBL)
C     we do not calculate the averages of these values as they
C     are meaningless.
C     DVDMAV(MXETBL),DVDMSQ(MXETBL)
C     DVDLAV(MXETBL),DVDLSQ(MXETBL)

C arrays for additional force field terms..and averages thereof
      real*8 ENERES(MXCTBL),EREAVE(MXCTBL),EREFLC(MXCTBL)
      DOUBLE PRECISION ERESUM(MXCTBL),ERESQ(MXCTBL)

C VOLPRT: some sytem parameters saved to file.
C NOTE that the values in this array do NOT necessarily reflect
C the current status of the system at all times;
C it is only updated when needed for printing or writing to file.
C Hence, e.g., when using the pressure of the system in a calculation
C within the code, the variable PRESTO should be used, and not the
C corresponding element in VOLPRT.
      real*8 VOLPRT(MXVTBL),VPAVE(MXVTBL),VPFLC(MXVTBL)
      DOUBLE PRECISION VPSUM(MXVTBL),VPSQ(MXVTBL)
C variables used for saving configs (NTWSE .NE. 0)
      INTEGER NLIM
      LOGICAL LSEWRT,LSECPY,LSEFST

C auxiliary array for leap with perturbations
      INTEGER NSPG(MAXPAT)
      real*8 DLAM,DMU
C
      CHARACTER PRGSTR*(5)

C external function definition
      LOGICAL LTENER
      EXTERNAL LTENER

C data statements
      DATA PRGSTR /'RUNMD'/
      DATA LSEFST /.TRUE./

#include "inc_runmd5.h"

      INCLUDE 'mm_save_config.h'  
C begin

c   saves argument in the common arg for qmmm
      NRPT=NRPTp
      NSOLV=NSOLVp
      NATTOT=NATTOTp
      NUM4DA=NUM4DAp
      IUTRJX=IUTRJXp
      IUTRJV=IUTRJVp
      IUTRJE=IUTRJEp
      IUTRJG=IUTRJGp
      IFNCFG=IFNCFGp
      TIME=TIMEp
!$OMP parallel do private(iii)
      DO iii=1,MXETBL
        DLSUM(iii)=DLSUMp(iii)
        DMSUM(iii)=DMSUMp(iii)
      ENDDO
c
      NRPp=NRP
      NSMp=NSM
      MAXNRSp= MAXNRS
      MAXNRPp= MAXNRP
      MAXNATp= MAXNAT
#ifdef WRITE_FORCES
      IUTRJFp = IUTRJF
#endif
#ifdef READ_COORD
      IURDCOp = IURDCO
#endif
C init my local vars
      CW4D = CW4DA
      LPR4 = (NT4DIM .NE. NT4OFF)
      LWR4 = (NT4DIM .NE. NT4OFF .AND. NT4XO .EQ. NT4XOY)

      LDOPER = (NTG .NE. NTGOFF)
      LDODLA = (NTG .EQ. NTGLAM .OR. NTG .EQ. NTGBOT)
      LDODMU = (NTG .EQ. NTGMU .OR. NTG .EQ. NTGBOT)
      LDOVDL = (LDODLA .AND. DLAMT .NE. 0.0D0)
      LDOVDM = (LDODMU .AND. DMUT .NE. 0.0D0)

      LEVERY = (NTP .NE. NTPOFF)
      LFORM  = (NTPW .EQ. NTPWFO)

      LDOPIN = (NTPI .EQ. NTPION)

      LSEWRT = .FALSE.
      LSECPY = .FALSE.

      DTHALF = DT*0.5D0
      BOHALF = BOLTZ*0.5D0

      NRPT3  = NRPT*NDIM
      NSOLV3 = NSOLV*NDIM
      NATTO3 = NATTOT*NDIM
      NITP = 0
      NITS = 0

C do we shake the solute ?
      LSHKPR = (NTC .NE. NTCNON .AND. NPM .GT. 0 .AND. NRP .GT. 1)
C do we shake the solvent ?
      LSHKSO = (NSM .GT. 0 .AND. NRAM .GT. 1)

!$OMP parallel do private(II)
      DO II=1,MXETBL
          EPSUM(II) = 0.0D0
           EPSQ(II) = 0.0D0
          EPAVE(II) = 0.0D0
          EPFLC(II) = 0.0D0

         DEDLSU(II) = 0.0D0
         DEDLSQ(II) = 0.0D0
         DEDLAV(II) = 0.0D0
         DEDLFL(II) = 0.0D0

         DEDMSU(II) = 0.0D0
         DEDMSQ(II) = 0.0D0
         DEDMAV(II) = 0.0D0
         DEDMFL(II) = 0.0D0
      ENDDO

!$OMP parallel do private(II)
      DO II=1,MXCTBL
         ERESUM(II) = 0.0D0
         ERESQ(II)  = 0.0D0
         EREAVE(II) = 0.0D0
         EREFLC(II) = 0.0D0
      ENDDO

!$OMP parallel do private(II)
      DO II=1,MXVTBL
         VPSUM(II)  = 0.0D0
         VPSQ(II)   = 0.0D0
         VPAVE(II)  = 0.0D0
         VPFLC(II)  = 0.0D0
      ENDDO

!$OMP parallel do private(II)
      DO II=1,NFTMAX
         TMPSUM(II) = 0.0D0
         TMPSQ(II)  = 0.0D0
         TMPAVE(II) = 0.0D0
         TMPFLC(II) = 0.0D0
      ENDDO

!$OMP parallel do private(II)
      DO II=1,MXNRE2
         EPLJ(II)   = 0.0D0
         EPLJSU(II) = 0.0D0
         EPLJSQ(II) = 0.0D0
         EPLJAV(II) = 0.0D0
         EPLJFL(II) = 0.0D0

         EPEL(II)   = 0.0D0
         EPELSU(II) = 0.0D0
         EPELSQ(II) = 0.0D0
         EPELAV(II) = 0.0D0
         EPELFL(II) = 0.0D0

         EPRF(II)   = 0.0D0
         EPRFSU(II) = 0.0D0
         EPRFSQ(II) = 0.0D0
         EPRFAV(II) = 0.0D0
         EPRFFL(II) = 0.0D0

         EPRC(II)   = 0.0D0
         EPRCSU(II) = 0.0D0
         EPRCSQ(II) = 0.0D0
         EPRCAV(II) = 0.0D0
         EPRCFL(II) = 0.0D0
      ENDDO

      VIRTOT = 0.0D0
      PRESTO = 0.0D0
      EKCMTO = 0.0D0
!$OMP parallel do private(II)
      DO II=1,MAXDIM
         VIR(II) = 0.0D0
         PRES(II)= 0.0D0
         EKCM(II)= 0.0D0
         XMIN(II)= 0.0D0
      ENDDO

C     calc degrees of freedom for the solute.
C     calc number of constraints of the solute.
C     If we have position constraints, the WMAS and WINV
C     arrays are also modified in this routine.
      CALL CLCDOF(NRP,NRAM,NCONS,NCONG,NSPM,NRCON,JRC,V,
     $     WINV,WMAS,NUM4DA,
     $     NBON,IB,JB,
     $     NBONH,IBH,JBH,
     $     NSP,
     $     NCONP,NCONBC,NDOF)

      PRINT *,'NDEGBLOCK'
      PRINT 598,NRCON
      PRINT 599,NCONBC
      PRINT 600,NCONP
      PRINT 605,NDOF(NFTPIR)
      PRINT 606,NDOF(NFTPCM)
      PRINT 610,NCONS
      PRINT 615,NDOF(NFTSLV)
      PRINT 620,NDFMIN
      PRINT 625,NDOF(NFTTO)

      IF (LPR4) THEN
         PRINT *
         PRINT 650,NUM4DA
         PRINT 660,NDFMI4
         PRINT 665,NDOF(NFT4D)
         PRINT *
         IF (NDOF(NFT4D) .LE. 0) THEN
            PRINT *,PRGSTR,': illegal number of degrees',
     $           ' of freedom NDF4 =',NDOF(NFT4D)
            CALL STOPGM('MM_SETUP_DR','FATAL ERROR')
         ENDIF
      ENDIF
      PRINT *,'END NDEGBLOCK'
      PRINT *

      IF (NDOF(NFTTO) .LE. 0) THEN
         PRINT *,PRGSTR,': illegal number of degrees',
     $        ' of freedom NDF = ',NDOF(NFTTO)
         CALL STOPGM('MM_SETUP_DR','FATAL ERROR')
      ENDIF

C init temp coupling arrays, and also write out info
      CALL CLCBTH(NDOF,NBATH,NBNUM,NBNDX,EKREF,DTBATH,TFACBT)

C init TFAC array for printing
!$OMP parallel do private(II)
      DO II=1,NFTMAX
         SCAL(II)   = 1.0D0
         TFACPR(II) = BOHALF*NDOF(II)
      ENDDO

C     EKIN40      : reference ekin at beginning of run
C     EKIN4B      : reference ekin at end of run
      EKIN40 = TFACPR(NFT4D)*TEMP04
      EKIN4B = TFACPR(NFT4D)*TEMP0B

C pressure coupling
      IF (NTP .EQ. NTPOFF) THEN
         DTCP = 0.0D0
      ELSE
         DTCP = COMP*DT/TAUP
      ENDIF

C boundary conditions
      LMONO = (NTB .GT. 0)
      LOCTO = (NTB. LT. 0)
      LVAC  = (NTB .EQ. 0)
      LDOVIR= (ABS(NTB) .EQ. NTBVIR)

C init NSPG if required
      IF (LDOPER) THEN
         DLAM = DLAMT*DT
         DMU  = DMUT*DT
         DO 30 K = 1,NJLA
            II = JLA(K)
            N = 1
 40         IF (N .LE. NSPM) THEN
               IF (II .GT. NSP(N)) THEN
                  N = N + 1
                  GOTO 40
               ENDIF
            ENDIF
            NSPG(K) = N
 30      CONTINUE
      ELSE
         DLAM = 0.0D0
         DMU  = 0.0D0
      ENDIF

C     calc total mass of the system.
C     this must be done after CLCDOF, as the WMAS arrays could
C     have been changed in that routine (position constraining)
C     Note that the WMAS and WINV arrays are also modified
C     where the mass of an atom is perturbed.
C     This is done in L<RDPERT>.
      CALL CLCMAS(NPM,NSM,NSPM,NSP,TOTMAS,TMASS,SUBMAS,SUBMIN)

C apply initial shake if required
      IF (INIT .NE. INITCO .AND. INIT .NE. INNOSH) THEN
         CALL INICRD(IFNCFG,NATTOT,NRPT,NDIM,PRGSTR,X,V,F,
     $        C4D,NCONG,ICOG,JCOG,CONP,FCON,NRCON,JRC)
      ENDIF

C stop COM motion if required, else just print it out
!$OMP parallel do private(II)
      DO II=1,NATTO3
         F(II) = X(II) - V(II)*DTHALF
      ENDDO

C Skip over centre of mass calculations if atoms are positionally
C constrained
      IF (NTR .NE. NTRCON .AND. LTCM .AND. INIT .NE. INITCO) THEN
         LTCM = .FALSE.
!        CALL CENMAS(NATTOT,NPM,NRP,NSM,NRAM,0,NDIM,NDRMAX,F,V,
!    $        TOTMAS,0,WMAS,WMASS,
!    $        EKCMTO,XCM,VCM,ACM,EKROT,OCM,ICMROT)
!        CALL STOPCM(NATTOT,NDIM,F,V,XCM,VCM,OCM,ISCROT)
      ENDIF

C (re)calc COM motion
!     CALL CENMAS(NATTOT,NPM,NRP,NSM,NRAM,0,NDIM,NDRMAX,F,V,
!    $     TOTMAS,0,WMAS,WMASS,
!    $     EKCMTO,XCM,VCM,ACM,EKROT,OCM,ICMROT)

      IF (.NOT. LVAC .AND. NPM .GT. 0) THEN
         CALL CONCG(NRPT,NPM,X,LEVERY)
      ENDIF

C-- calc initial EKIN
      CALL CLCTEM(NATTOT,NPM,NSM,NSPM,NSP,SUBMIN,V,ESCL,VELCM)

C print out a measure of 80 chars for checking output.
C we should not exceed 80 chars in width
CCC
C     PRINT '(8(I10))',(NSTEP, NSTEP=1,8)
CCC
!$OMP parallel do private(II)
      DO II=1,NFTMAX
         IF (TFACPR(II) .GT. EPS) THEN
            TEMP(II) = ESCL(II)/TFACPR(II)
         ELSE
            TEMP(II) = 0.0D0
         ENDIF
      ENDDO

      PRINT *,'INITIAL TEMPERATURES ARE :'
      CALL PRTEMP(TEMP)

      DO II=1,NBATH
         EKTMP = 0.0D0
         DO K=1,NBNUM(II)
            EKTMP = EKTMP + ESCL(NBNDX(II,K))
         ENDDO

C divide by zero measure:
C this is a counter divide-by-zero measure for the first step:
C if we have T coupling and zero initial velocities,
C we would get a divide by zero error in
C calculating SCALP, SCALS or SCAL4 ==> set energy to the
C reference temperature

         IF (EKTMP .LE. EPS) THEN
            EKBATH(II) = EKREF(II)
         ELSE
            EKBATH(II) = EKTMP
         ENDIF
      ENDDO

      IFSTEP = 0
      CALL mm_AZZERO(FITAR,MAXFIT)
c
c    QMMM : end setup 
c
      RETURN

C if we write a minimum energy trajectory, then we should write
C the saved configuration and energy of the last time frame here.
c     IF (NTWSE .NE. 0) THEN
c        CALL WRSAVE(IUTRJX,LFORM,
c    $        NTWX .GT. 0,NTB .NE. NTBVAC,LWR4,
c    $        LPRKIN,LDOVIR,LDOPIN,
c    $        NUSNRE,NUNRE2,NATTOT,NRPT,NDIM,XSAVE)
c     ENDIF

c calc averages and r.m.s. fluctuations and print them
c     CALL CLCAVE(NSTLIM,MXETBL,EPSUM,EPSQ,EPAVE,EPFLC)
c     CALL CLCAVE(NSTLIM,MXCTBL,ERESUM,ERESQ,EREAVE,EREFLC)
c     CALL CLCAVE(NSTLIM,MXVTBL,VPSUM,VPSQ,VPAVE,VPFLC)
c     CALL CLCAVE(NSTLIM,NFTMAX,TMPSUM,TMPSQ,TMPAVE,TMPFLC)

c     CALL CLCAVE(NSTLIM,MXNRE2,EPLJSU,EPLJSQ,EPLJAV,EPLJFL)
c     CALL CLCAVE(NSTLIM,MXNRE2,EPELSU,EPELSQ,EPELAV,EPELFL)
c     CALL CLCAVE(NSTLIM,MXNRE2,EPRFSU,EPRFSQ,EPRFAV,EPRFFL)
c     CALL CLCAVE(NSTLIM,MXNRE2,EPRCSU,EPRCSQ,EPRCAV,EPRCFL)

c     CALL CLCAVE(NSTLIM,MXETBL,DEDLSU,DEDLSQ,DEDLAV,DEDLFL)
c     CALL CLCAVE(NSTLIM,MXETBL,DEDMSU,DEDMSQ,DEDMAV,DEDMFL)

c     PRINT 5010
c     PRINT 5020
c     CALL PRARR(LPRKIN,LDOVIR,LDOPIN,
c    $     .TRUE.,EPAVE,EREAVE,VPAVE,TMPAVE,
c    $     LDODLA,DEDLAV,
c    $     .FALSE.,EDUM,
c    $     LDODMU,DEDMAV,
c    $     .FALSE.,EDUM,
c    $     NUSNRE,NUNRE2,EPELAV,EPRFAV,EPRCAV,EPLJAV)

c     PRINT 5030
c     CALL PRARR(LPRKIN,LDOVIR,LDOPIN,
c    $     .TRUE.,EPFLC,EREFLC,VPFLC,TMPFLC,
c    $     LDODLA,DEDLFL,
c    $     .FALSE.,EDUM,
c    $     LDODMU,DEDMFL,
c    $     .FALSE.,EDUM,
c    $     NUSNRE,NUNRE2,EPELFL,EPRFFL,EPRCFL,EPLJFL)


c print out the sum of dE/dLamda*DLamda and corresponding for Mu
c     IF (NTG .NE. NTGOFF) THEN
c        PRINT 5040
c        PRINT *,'FREEENERGYCOMPS (dE/dLamda)*DLamda'
c        CALL PRENER(.TRUE.,DLSUM,.FALSE.,EDUM,LDOPIN,.FALSE.,
c    $        .FALSE.,EDUM,EDUM)

c
c        PRINT *,'FREEENERGYCOMPS (dE/dMu)*DMu'
c        CALL PRENER(.TRUE.,DMSUM,.FALSE.,EDUM,LDOPIN,.FALSE.,
c    $        .FALSE.,EDUM,EDUM)
c     ENDIF

c     PRINT *
C end of code

C formats
 500  FORMAT (1X,A6,2A5,(1X,A12  ))
 501  FORMAT (1X,I6,2I5,(1X,E12.5))

 598  FORMAT (4X,'NUMBER OF POSITIONALLY CONSTRAINED ATOMS   :',I8)
 599  FORMAT (4X,'NUMBER OF BONDS BETWEEN POS. CONSTR. ATOMS :',I8)

 600  FORMAT (4X,'NUMBER OF SOLUTE CONSTRAINTS               :',I8)
 605  FORMAT (4X,'NUMBER OF 3D SOLUTE (INT+ROT) D.O.F        :',I8)
 606  FORMAT (4X,'NUMBER OF 3D SOLUTE C.O.M. TRANSL. D.O.F   :',I8)

 610  FORMAT (4X,'NUMBER OF SOLVENT CONSTRAINTS              :',I8)
 615  FORMAT (4X,'NUMBER OF 3D SOLVENT D.O.F.                :',I8)
 620  FORMAT (4X,'NUMBER OF 3D D.O.F. TO SUBTRACT (NDFMIN)   :',I8)
 625  FORMAT (4X,'NUMBER OF TOTAL D.O.F. IN 3D               :',I8)

 650  FORMAT (4X,'NUMBER OF ATOMS IN 4D                      :',I8)
 660  FORMAT (4X,'NUMBER OF 4th D D.O.F. TO SUBTRACT (NDFMI4):',I8)
 665  FORMAT (4X,'NUMBER OF TOTAL D.O.F in 4th D             :',I8)

C title start from 5000
 5000 FORMAT (//,1X,'6. D A T A   P E R   S T E P',//)
 5010 FORMAT (//,1X,'7. A V E R A G E S   O V E R   R U N')
 5020 FORMAT (//,1X,'A: A V E R A G E S')
 5030 FORMAT (//,1X,'B: R O O T   M E A N   S Q U A R E   ',
     $     'F L U C T U A T I O N S')
 5040 FORMAT (//,1X,'C: C U M U L A T I V E   V A L U E S')
C end runmd
      END
C     ==================================================================
      SUBROUTINE SAVTIM(NL,X,NDIM,BOX,TIME,NSTEP,XSAVE)
C     ==--------------------------------------------------------------==
C     SUBROUTINE SAVTIM(NL,X,NDIM,BOX,TIME,NSTEP,XSAVE)
C
C     SAVTIM saves the time and coordinates for
C     later writing to a minimum energy trajectory
C     ( L<NTWSE> .NE. 0) by subroutine L<WRSAVE>.
C
C     SAVTIM is called by L<RUNMD>.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'md.h'
      INCLUDE 'runmd.h'
C args
      INTEGER NL,NSTEP,NDIM
      real*8 X(NL),XSAVE(NL),BOX(NDIM),TIME
COMMON
      real*8 BOXSAV
      real*8 TIMSAV,ENRGSE,ERESSE,VOPSAV,TEMPSV
      real*8 EPLJSE,EPELSE,EPRFSE,EPRCSE
      INTEGER NSTSAV,NLIM
      COMMON /RMDSRE/BOXSAV(MAXDIM),TIMSAV,ENRGSE(MXETBL),
     $     ERESSE(MXCTBL),VOPSAV(MXVTBL),TEMPSV(NFTMAX),
     $     EPLJSE(MXNRE2),EPELSE(MXNRE2),EPRFSE(MXNRE2),
     $     EPRCSE(MXNRE2)
      COMMON /RMDSIN/NSTSAV,NLIM
COMMONEND
C local vars
      INTEGER II
C begin
      NLIM = NL
!$OMP parallel do private(II)
      DO II=1,NLIM
         XSAVE(II) = X(II)
      ENDDO
!$OMP parallel do private(II)
      DO II=1,NDIM
         BOXSAV(II) = BOX(II)
      ENDDO
      TIMSAV = TIME
      NSTSAV = NSTEP
C savtim
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SAVNRG(ENER,ENERES,VOLPRT,TEMP,
     $     EPLJ,EPEL,EPRF,EPRC)
C     ==--------------------------------------------------------------==
C     SUBROUTINE SAVNRG(ENER,ENERES,VOLPRT,TEMP,
C     $     EPLJ,EPEL,EPRF,EPRC)
C
C     SAVNRG saves the values of the arrays passed to it
C     for minimum energy trajectory writing (L<NTWSE> .NE. 0)
C     The thus saved energies can be written to 
C     trajecotory by L<WRSAVE>.
C     SAVNRG is called by L<RUNMD>
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'md.h'
      INCLUDE 'runmd.h'
C args
      real*8 ENER(MXETBL),ENERES(MXCTBL)
      real*8 VOLPRT(MXVTBL),TEMP(NFTMAX)
      real*8 EPLJ(MXNRE2),EPEL(MXNRE2),EPRF(MXNRE2),EPRC(MXNRE2)
COMMON
      real*8 BOXSAV
      real*8 TIMSAV,ENRGSE,ERESSE,VOPSAV,TEMPSV
      real*8 EPLJSE,EPELSE,EPRFSE,EPRCSE
      INTEGER NSTSAV,NLIM
      COMMON /RMDSRE/BOXSAV(MAXDIM),TIMSAV,ENRGSE(MXETBL),
     $     ERESSE(MXCTBL),VOPSAV(MXVTBL),TEMPSV(NFTMAX),
     $     EPLJSE(MXNRE2),EPELSE(MXNRE2),EPRFSE(MXNRE2),
     $     EPRCSE(MXNRE2)
      COMMON /RMDSIN/NSTSAV,NLIM
COMMONEND
C local vars
      INTEGER II
C begin
!$OMP parallel do private(II)
      DO II=1,MXETBL
         ENRGSE(II) = ENER(II)
      ENDDO

!$OMP parallel do private(II)
      DO II=1,MXCTBL
         ERESSE(II) = ENERES(II)
      ENDDO

!$OMP parallel do private(II)
      DO II=1,MXNRE2
         EPLJSE(II) = EPLJ(II)
         EPELSE(II) = EPEL(II)
         EPRFSE(II) = EPRF(II)
         EPRCSE(II) = EPRC(II)
      ENDDO

!$OMP parallel do private(II)
      DO II=1,MXVTBL
         VOPSAV(II) = VOLPRT(II)
      ENDDO

!$OMP parallel do private(II)
      DO II=1,NFTMAX
         TEMPSV(II) = TEMP(II)
      ENDDO
C savnrg
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      LOGICAL FUNCTION LTENER(ECMP)
C     ----------------------------------------------------------------==
C     Logical function LTENER returns .TRUE. if
C     ECMP is smaller than the currently saved
C     energy value. LTENER is used in minimum 
C     energy trajectory writing (L<NTWSE> .NE. 0)
C     and is called by L<RUNMD>.
C     ----------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'md.h'
      INCLUDE 'runmd.h'
C args
      real*8 ECMP
COMMON
      real*8 BOXSAV
      real*8 TIMSAV,ENRGSE,ERESSE,VOPSAV,TEMPSV
      real*8 EPLJSE,EPELSE,EPRFSE,EPRCSE
      INTEGER NSTSAV,NLIM
      COMMON /RMDSRE/BOXSAV(MAXDIM),TIMSAV,ENRGSE(MXETBL),
     $     ERESSE(MXCTBL),VOPSAV(MXVTBL),TEMPSV(NFTMAX),
     $     EPLJSE(MXNRE2),EPELSE(MXNRE2),EPRFSE(MXNRE2),
     $     EPRCSE(MXNRE2)
      COMMON /RMDSIN/NSTSAV,NLIM
COMMONEND
C begin
      LTENER = (ECMP .LT. ENRGSE(NTWSE))
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE WRSAVE(IUTRJX,LFORM,LWRALL,LWRBOX,LWR4,
     $     LPRKIN,LDOVIR,LDOPIN,
     $     NUSNRE,NUNRE2,NATTOT,NRPT,NDIM,XSAVE)
C     ----------------------------------------------------------------==
C     WRSAVE writes out the data needed in a minimum energy
C     trajectory ( L<NTWSE> .NE. 0).
C     The data is saved to a common block by routines
C     L<SAVNRG> and L<SAVTIM>.
C     It is called by L<RUNMD>.
C     ----------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'md.h'
      INCLUDE 'runmd.h'
C args
      INTEGER IUTRJX,NATTOT,NRPT,NDIM,NUSNRE,NUNRE2
      LOGICAL LFORM,LWRALL,LWRBOX,LWR4,LPRKIN,LDOVIR,LDOPIN
      real*8 XSAVE(NDIM*NATTOT)
COMMON
      real*8 BOXSAV
      real*8 TIMSAV,ENRGSE,ERESSE,VOPSAV,TEMPSV
      real*8 EPLJSE,EPELSE,EPRFSE,EPRCSE
      INTEGER NSTSAV,NLIM
      COMMON /RMDSRE/BOXSAV(MAXDIM),TIMSAV,ENRGSE(MXETBL),
     $     ERESSE(MXCTBL),VOPSAV(MXVTBL),TEMPSV(NFTMAX),
     $     EPLJSE(MXNRE2),EPELSE(MXNRE2),EPRFSE(MXNRE2),
     $     EPRCSE(MXNRE2)
      COMMON /RMDSIN/NSTSAV,NLIM
COMMONEND
C local vars
      real*8 XDUM
C begin
      CALL WRTIME(IUTRJX,LFORM,NSTSAV,TIMSAV)
      CALL WRNRG(IUTRJX,LFORM,
     $     MXEWRT,ENRGSE,
     $     MXCTBL,ERESSE,
     $     NUSNRE,EPLJSE,EPELSE,EPRFSE,EPRCSE)
      IF (LWRALL) THEN
         CALL WRXRED(IUTRJX,LFORM,LWR4,NATTOT,NDIM,XSAVE)
      ELSE
C     only write solute coords
         CALL WRXRED(IUTRJX,LFORM,LWR4,NRPT,NDIM,XSAVE)
      ENDIF
      IF (LWRBOX) THEN
         CALL WRBOX(IUTRJX,LFORM,BOXSAV)
      ENDIF
C also print the energies
      PRINT *,'MINIMUM ENERGY BLOCK'
      PRINT 500,'STEP','TIME'
      PRINT 501,NSTSAV,TIMSAV
      CALL PRARR(LPRKIN,LDOVIR,LDOPIN,
     $     .TRUE.,ENRGSE,ERESSE,VOPSAV,TEMPSV,
     $     .FALSE.,XDUM,
     $     .FALSE.,XDUM,
     $     .FALSE.,XDUM,
     $     .FALSE.,XDUM,
     $     NUSNRE,NUNRE2,EPELSE,EPRFSE,EPRCSE,EPLJSE)
      PRINT *,'END MINIMUM ENERGY BLOCK'
      PRINT *

 500  FORMAT (1X,A6,(1X,A12  ))
 501  FORMAT (1X,I6,(1X,E12.5))
C wrsave
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE CLCBTH(NDOF,NBATH,NBNUM,NBNDX,EKREF,DTBATH,TFACBT)
C     ----------------------------------------------------------------==
C     Subroutine CLCBTH calculates some auxiliary variables involved in
C     temperature coupling from the input (L<NTT>, L<TEMP0> and L<TAUT>)
C     in order to handle the many different cases more easily.
C     It is called in subroutine L<RUNMD>.
C     ----------------------------------------------------------------==
C includes
      INCLUDE 'md.h'
      INCLUDE 'runmd.h'
C args
      INTEGER NDOF(NFTMAX)
      INTEGER NBATH,NBNUM(MXBATH),NBNDX(MXBATH,MXBATH)
      real*8 EKREF(MXBATH),DTBATH(MXBATH),TFACBT(MXBATH)
C local vars
      INTEGER I,J,K,IABSVL
      LOGICAL LINEG,LJNEG,LKNEG
      LOGICAL LCHK(NTTNUM)
      INTEGER NDOFBT(MXBATH)
      real*8 BOHALF
C translate from NTT array to NFT indices
      INTEGER ITRANS(NTTNUM)
C data
      DATA ITRANS(NTTPIR) /NFTPIR/
      DATA ITRANS(NTTPCM) /NFTPCM/
      DATA ITRANS(NTTSLV) /NFTSLV/
C begin
      BOHALF = BOLTZ*0.5D0
      NBATH = 0
!$OMP parallel do private(I)
      DO I=1,NTTNUM
         LCHK(I) = .TRUE.
      ENDDO

C work through the NTT switches
      DO 30 I=1,NTTNUM
         IF (LCHK(I) .AND. NTT(I) .NE. NTTOFF) THEN
C do something depending on value of switch
            IABSVL = ABS(NTT(I))
            LINEG   = (NTT(I) .LT. NTTOFF)
            IF (IABSVL .EQ. NTTONE) THEN
C coupling to one bath
               NBATH = NBATH + 1
               NBNUM(NBATH) = 1
               NBNDX(NBATH,1) = ITRANS(I)
               EKREF(NBATH)   = TEMP0(I)
               DTBATH(NBATH)  = TAUT(I)
               NDOFBT(NBATH)  = NDOF(I)
            ELSEIF (IABSVL .EQ. NTTTWO) THEN
C look for another
               J=I+1
C while loop
 50            IF (J .LE. NTTNUM) THEN
                  IF (ABS(NTT(J)) .NE. IABSVL) THEN
                     J = J+1
                     GOTO 50
                  ENDIF
               ENDIF
C end while
               IF (J .LE. NTTNUM) THEN
                  LJNEG = (NTT(J) .LT. NTTOFF)
C found, do something
C we couple to baths together
                  NBATH = NBATH + 1
                  NBNUM(NBATH) = 2
                  NBNDX(NBATH,1) = ITRANS(I)
                  NBNDX(NBATH,2) = ITRANS(J)
                  IF (LJNEG) THEN
                     EKREF(NBATH)   = TEMP0(I)
                     DTBATH(NBATH)  = TAUT(I)
                  ELSE
                     EKREF(NBATH)   = TEMP0(J)
                     DTBATH(NBATH)  = TAUT(J)
                  ENDIF
                  NDOFBT(NBATH)  = NDOF(I) + NDOF(J)
C everything is OK: mark J
                  LCHK(J) = .FALSE.
               ENDIF
            ELSEIF (IABSVL .EQ. NTTHRE) THEN
               IF (I .EQ. 1) THEN
                  J = I+1
                  K = J+1
                  IF (ABS(NTT(J)) .EQ. NTTHRE .AND.
     $                 ABS(NTT(K)) .EQ. NTTHRE) THEN
                     LJNEG = (NTT(J) .LT. NTTOFF)
                     LKNEG = (NTT(K) .LT. NTTOFF)
C we have three baths together --> couple to the TOTAL
C system (NDFMIN comes into play)
                  NBATH = NBATH + 1
CCC
C     NBNUM(NBATH) = 3
C     NBNDX(NBATH,1) = ITRANS(I)
C     NBNDX(NBATH,2) = ITRANS(J)
C     NBNDX(NBATH,3) = ITRANS(K)
C     NDOFBT(NBATH) = NDOF(I) + NDOF(J) + NDOF(K)
CCC
                  NBNUM(NBATH) = 1
                  NBNDX(NBATH,1) = NFTTO
                  NDOFBT(NBATH)  = NDOF(NFTTO)

                  IF (.NOT. LINEG) THEN
                     EKREF(NBATH)   = TEMP0(I)
                     DTBATH(NBATH)  = TAUT(I)
                  ELSEIF (.NOT. LJNEG) THEN
                     EKREF(NBATH)   = TEMP0(J)
                     DTBATH(NBATH)  = TAUT(J)
                  ELSEIF (.NOT. LKNEG) THEN
                     EKREF(NBATH)   = TEMP0(K)
                     DTBATH(NBATH)  = TAUT(K)
                  ENDIF
C mark the J and K items
                     LCHK(J) = .FALSE.
                     LCHK(K) = .FALSE.
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
 30   CONTINUE

C now handle the 4D bath switch
      IF (NTT4 .NE. NTT4OF) THEN
         NBATH = NBATH + 1
         NBNUM(NBATH) = 1
         NBNDX(NBATH,1) = NFT4D
         EKREF(NBATH)   = TEMP04
         DTBATH(NBATH)  = TAUT4
         NDOFBT(NBATH)  = NDOF(NFT4D)
      ENDIF

C write out information about the temperatur baths
      PRINT *,'BATHBLOCK'
      PRINT 700,NBATH
      IF (NBATH .GT. 0) THEN
         PRINT 705,'TEMP0','TAUT','DEG-OF-FR.'
         DO I=1,NBATH
            PRINT 710,EKREF(I),DTBATH(I),NDOFBT(I)
         ENDDO
      ENDIF
      PRINT *,'END BATHBLOCK'
      PRINT *

C now calc DTBATH: DT/TTAUT
C and TFACBT: o.5* kb* degrees of freedom
C and the ekin of the reference bath (from temperature to Ekin)
!$OMP parallel do private(I)
      DO I=1,NBATH
         DTBATH(I) = DT/DTBATH(I)
         TFACBT(I) = BOHALF*NDOFBT(I)
         EKREF(I)  = EKREF(I)*TFACBT(I)
      ENDDO

 700  FORMAT (4X,'TEMPERATURE COUPLING TO',I3,' BATH(S)')
 705  FORMAT (3A12)
 710  FORMAT (2(1X,F11.3),I12)
C end
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE MDLEAP(IFNCFG,NATTOT,NPM,NRP,NSM,NRAM,N4DBON,
     $     NSPM,NSP,SUBMAS,SUBMIN,
     $     WINV,WMAS,WINVS,WMASS,DT,
     $     X,V,F,C4D,VELCM,
     $     LSHKPR,NCONG,ICOG,JCOG,CONP,FCON,NRCON,JRC,SHKTOL,
     $     LSHKSO,NCONS,ICONS,JCONS,CONS,
     $     SCAL,EKNOW,EKNEXT,
     $     LEVERY,NITP,NITS)
C     ==--------------------------------------------------------------==
C     MDLEAP performs a leap-frog integration with L<SHAKE> if
C     specified.
C     If pressure coupling is specified, the resulting coordinates
C     are scaled appropriately.
C
C     NOTE: in order to save space, the F array is used as
C     temporary storage in this routine.
C     THUS, THE FORCES ARE NO LONGER IN THE ARRAY F AFTER THIS ROUTINE.
C
C     IFNCFG  unit number onto which the configuration will be written
C             in case of a SHAKE error
C     NATTOT  the total number of atoms in the system
C     X,V,F   coordinate, velocity, force 
C              arrays, respectively.
C     NCONG   number of constraints in ICOG,JCOG and CONP
C     ICOG,JCOG atoms forming constraints
C     CONP    constraints for shake
C     SCALP   scaling factor for solute temperature bath.
C     SCALS   sclaing factor for solvent temperature bath.
C     SCAL4   scaling factor in the fourth dim.
C
C     EKPR3   delivered with kinetic energy of solute in 3D at time T
C     EKPR4   delivered with kinetic energy of solute in 4D at time T
C     EKSO3   delivered with kinetic energy of solvent in 3D at time T
C     EKSO4   delivered with kinetic energy of solvent in 3D at time T
C     NITP    delivered with number of L<SHAKE> iterations for solute
C     NITS    delivered with number of L<SHAKE> iterations fo solvent
C
C     THE FOLLOWING ENERGIES ARE AT TIME T+DT/2 !!
C     The above kinetic energies are used for calculating the total
C     energy at time T.
C     The following kinetic energies are used for scaling the
C     temperature in the next step, i.e. T+DT.
C     ESCLP   delivered with kinetic energy of solute
C             in 3D at time T+DT/2
C     ESCLS   delivered with kinetic energy of solvent
C             in 3D at time T+DT/2
C     ESCL4   delivered with kinetic energy in the fourth dim.
C             at time T+DT/2
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'toposz.h'
      INCLUDE 'formats.h'
      INCLUDE 'runmd.h'
C args
      INTEGER NATTOT,NPM,NRP,NSM,NRAM,NSPM,NCONS
      INTEGER IFNCFG,NSP(NSPM)
      INTEGER NCONG,ICOG(NCONG),JCOG(NCONG),NRCON,JRC(NRCON)
      INTEGER ICONS(NCONS),JCONS(NCONS),N4DBON
      LOGICAL LSHKPR,LSHKSO,LEVERY
      real*8 C4D(NATTOT),SUBMAS(NSPM),SUBMIN(NSPM),DT,SHKTOL
      real*8 WINV(NRP),WMAS(NRP),WINVS(NRAM),WMASS(NRAM)
      real*8 X(NATTOT*NDIM),V(NATTOT*NDIM),F(NATTOT*NDIM)
      real*8 VELCM(NPM*NSPM*NDIM)
      real*8 SCAL(NFTMAX),EKNOW(NFTMAX),EKNEXT(NFTMAX)

      real*8 CONP(NCONG),FCON(NCONG)
      real*8 CONS(NCONS)

      INTEGER NITP,NITS
C local params
      real*8 ONE8TH
      PARAMETER (ONE8TH = 1.0D0/8.0D0)
C local vars
      INTEGER I3,NN,J,M,K,IFIRST,ILAST,ICM3,I3FRST
      INTEGER NPROT,NPROT3,NSOLV,NSOLV3
      real*8 RMINV,RMASS,DTMINV,VOLD,VNEW,XOLD,VMNEW

      real*8 EKPR3,EKSO3,EKS4,EKP4,EKPCM
      real*8 ESCLP,ESCLS,ESCL4S,ESCL4P,ESPCM

      real*8 DTINV
      real*8 FCONS(MXCONS)
      real*8 SCLP(MAXDIM),SCLCM(MAXDIM),SCLS(MAXDIM)
      real*8 EKP(MAXDIM),ESP(MAXDIM),EKS(MAXDIM),ESS(MAXDIM)
      real*8 ESPSUB(MAXDIM),EKPSUB(MAXDIM),VELSUB(MAXDIM),VELCMN(MAXDIM)
      CHARACTER PRGSTR*(6)
C data
      DATA PRGSTR /'MDLEAP'/
C begin
C
      NPROT = NPM*NRP
      NPROT3= NDIM*NPROT

      NSOLV = NSM*NRAM
      NSOLV3= NDIM*NSOLV
      
      IF (NATTOT .NE. NPROT+NSOLV) THEN
         PRINT *,PRGSTR,': NATTOT inconsistency!'
         CALL STOPGM('MM_SETUP_DR','FATAL ERROR')
      ENDIF

      DTINV = 1.0D0/DT
!$OMP parallel do private(M)
      DO M=1,NDRMAX
         SCLP(M)  = SCAL(NFTPIR)
         SCLCM(M) = SCAL(NFTPCM) - SCAL(NFTPIR)
         SCLS(M)  = SCAL(NFTSLV)
      ENDDO
      SCLP(MAXDIM) = SCAL(NFT4D)
      SCLCM(MAXDIM)= 0.0D0
      SCLS(MAXDIM) = SCAL(NFT4D)

!$OMP parallel do private(M)
      DO M=1,MAXDIM
         EKP(M)    = 0.0D0
         ESP(M)    = 0.0D0
         ESPSUB(M) = 0.0D0
         EKPSUB(M) = 0.0D0
         EKS(M)    = 0.0D0
         ESS(M)    = 0.0D0
      ENDDO

      CALL mm_AZZERO(EKNOW,NFTMAX)

      IF (LSHKPR) THEN
C solute is shaken, first calc new free flight pos.
C Here we DONT calc Ekin, V or VCM but copy old X to F for shaking
C the new free flight X.
C We have to overwrite F with the unscaled free flight vel in the first
C loop over sub-molecules for scaling in the second one.

         I3 = 0
         DO 110 NN=1,NPM
            IFIRST = 1
            DO 115 K=1,NSPM
!$OMP parallel do private(M)
               DO M=1,NDIM
                  VELSUB(M) = 0.0D0
               ENDDO
               ILAST = NSP(K)
               I3FRST = I3
               DO J=IFIRST,ILAST
                  DTMINV = WINV(J)*DT
                  RMASS  = WMAS(J)
!$OMP parallel do private(M,VOLD,VNEW) shared(DTMINV,RMASS)
                  DO M=1,NDIM
                     VOLD = V(I3+M)
                     VNEW = VOLD + F(I3+M)*DTMINV
                     F(I3+M) = VNEW
                     VELSUB(M) = VELSUB(M) + RMASS*VNEW
                  ENDDO
                  I3 = I3+NDIM
               ENDDO

               RMINV = SUBMIN(K)
!$OMP parallel do private(M) shared(RMINV)
               DO M=1,NDIM
                  VELSUB(M) = VELSUB(M)*RMINV
               ENDDO

C now scale new velocity and calc new pos
               I3 = I3FRST
               DO J=IFIRST,ILAST
!$OMP parallel do private(M,VOLD,VNEW,XOLD) shared(DT)
                  DO M=1,NDIM
                     VOLD = F(I3+M)
                     VNEW = VOLD*SCLP(M) + VELSUB(M)*SCLCM(M)
                     XOLD = X(I3+M)
                     F(I3+M)= XOLD
                     X(I3+M)= XOLD + VNEW*DT
                  ENDDO
                  I3 = I3+NDIM
               ENDDO

               IFIRST = ILAST + 1
 115        CONTINUE
 110     CONTINUE

C F : old pos
C V : old vel
C VCM: old com velocity
C X : new free flight pos

C now shake solute and calc new V, VCM and ekin
         CALL mm_SHAKE(NPROT,NRP,NPM,N4DBON,WINV,0,
     $        NCONG,ICOG,JCOG,CONP,FCON,
     $        SHKTOL,F,X,C4D,
     $        NRCON,JRC,LEVERY,NITP,.FALSE.)

         IF (NITP .EQ. 0)THEN
            PRINT FMSHKP,PRGSTR
            CALL EMGSTP(IFNCFG)
         ENDIF

C the arrays now contain:
C F  : old pos (unchanged)
C V  : old vel (unchanged)
C VCM: old com vels (unchanged)
C X  : new shaken pos
         I3 = 0
         ICM3 = 0
         DO 240 NN=1, NPM
            IFIRST = 1
            DO 245 K=1,NSPM
!$OMP parallel do private(M)
               DO M=1,NDIM
                  VELSUB(M) = 0.0D0
               ENDDO
               ILAST = NSP(K)

               DO J=IFIRST,ILAST
                  RMASS = WMAS(J)
!$OMP parallel do private(M,VOLD,VNEW,VMNEW) shared(DTINV,RMASS)
                  DO M=1,NDIM
                     VOLD = V(I3+M)
                     VNEW = (X(I3+M) - F(I3+M))*DTINV
                     V(I3+M)= VNEW
                     VMNEW     = RMASS*VNEW
                     VELSUB(M) = VELSUB(M) + VMNEW
                     ESP(M)    = ESP(M) + VMNEW*VNEW
                     EKP(M)    = EKP(M) + RMASS * (VNEW+VOLD)**2
                  ENDDO
                  I3 = I3 + NDIM
               ENDDO
C
               RMASS = SUBMAS(K)
               RMINV = SUBMIN(K)
!$OMP parallel do private(M,VOLD,VNEW) shared(RMINV,RMASS)
               DO M=1,NDIM
                  VNEW = VELSUB(M)*RMINV
                  VOLD = VELCM(ICM3+M)
                  EKPSUB(M) = EKPSUB(M) + RMASS*(VOLD+VNEW)**2
                  ESPSUB(M) = ESPSUB(M) + RMASS*VNEW*VNEW
                  VELCM(ICM3+M) = VNEW
               ENDDO

               ICM3 = ICM3 + NDIM
               IFIRST = ILAST + 1
 245        CONTINUE
 240     CONTINUE
      ELSE
C just do free flight leap frog step
         I3 = 0
         ICM3 = 0
         DO 10 NN=1,NPM
            IFIRST = 1
            DO 20 K=1,NSPM
!$OMP parallel do private(M)
               DO M=1,NDIM
                  VELSUB(M) = 0.0D0
                  VELCMN(M) = 0.0D0
               ENDDO
               ILAST = NSP(K)
               I3FRST = I3
C calculate unscaled V 
               DO J=IFIRST,ILAST
                  RMASS  = WMAS(J)
                  DTMINV = WINV(J)*DT
!$OMP parallel do private(M,VOLD,VNEW) shared(DTMINV,RMASS)
                  DO M=1,NDIM
                     VOLD = V(I3+M)
                     VNEW = (VOLD + F(I3+M)*DTMINV)
                     F(I3+M) = VNEW
                     VELSUB(M) = VELSUB(M) + RMASS*VNEW
                  ENDDO
                  I3 = I3+NDIM
               ENDDO

               RMINV = SUBMIN(K)
!$OMP parallel do private(M)
               DO M=1,NDIM
                  VELSUB(M) = VELSUB(M)*RMINV
               ENDDO

               I3 = I3FRST
               DO J=IFIRST,ILAST
                  RMINV  = WINV(J)
                  DTMINV = RMINV*DT
                  RMASS  = WMAS(J)
!$OMP parallel do private(M,VOLD,VNEW,VMNEW) shared(DT,RMASS)
                  DO M=1,NDIM
                     VOLD = F(I3+M)
                     VNEW = VOLD*SCLP(M) + VELSUB(M)*SCLCM(M)
                     VOLD = V(I3+M)
                     V(I3+M) = VNEW
                     VMNEW     = RMASS*VNEW
                     VELCMN(M) = VELCMN(M) + VMNEW
                     ESP(M)    = ESP(M) + VMNEW*VNEW
                     EKP(M)    = EKP(M) + RMASS * (VNEW+VOLD)**2
                     X(I3+M)= X(I3+M) + VNEW*DT
                  ENDDO
                  I3 = I3+NDIM
               ENDDO
CC
               RMASS = SUBMAS(K)
               RMINV = SUBMIN(K)
!$OMP parallel do private(M,VOLD,VNEW) shared(RMINV,RMASS)
               DO M=1,NDIM
                  VNEW = VELCMN(M)*RMINV
                  VOLD = VELCM(ICM3+M)
                  EKPSUB(M) = EKPSUB(M) + RMASS*(VNEW+VOLD)**2
                  ESPSUB(M) = ESPSUB(M) + RMASS* VNEW**2
                  VELCM(ICM3+M) = VNEW
               ENDDO

               ICM3 = ICM3 + NDIM
               IFIRST = ILAST + 1
 20         CONTINUE
 10      CONTINUE
      ENDIF

      EKPR3 = 0.0D0
      ESCLP = 0.0D0
      ESPCM = 0.0D0
      EKPCM = 0.0D0
!$OMP parallel do private(M) reduction(+:EKPR3,ESCLP,ESPCM,EKPCM)
      DO M=NDRMIN,NDRMAX
         EKPR3 = EKPR3 + EKP(M)
         ESCLP = ESCLP + ESP(M)
         ESPCM = ESPCM + ESPSUB(M)
         EKPCM = EKPCM + EKPSUB(M)
      ENDDO
      EKPR3 = EKPR3*ONE8TH
      EKPCM = EKPCM*ONE8TH
      ESCLP = ESCLP*0.5D0
      ESPCM = ESPCM*0.5D0

      EKP4   = EKP(MAXDIM)*ONE8TH
      ESCL4P = ESP(MAXDIM)*0.5D0

C do solvent
      IF (LSHKSO) THEN
C solvent is shaken, first calc new free flight pos.
C Here we DONT overwrite V, and dont calc Ekin
C but copy old X to F for shaking
         I3 = NPROT3
         DO NN=1,NSM
            DO J=1,NRAM
               DTMINV = WINVS(J)*DT
!$OMP parallel do private(M,VOLD,VNEW,XOLD) shared(DTMINV,DT)
               DO M=1,NDIM
                  VOLD = V(I3+M)
                  VNEW = SCLS(M)*(VOLD + F(I3+M)*DTMINV)
                  XOLD    = X(I3+M)
                  F(I3+M) = XOLD
                  X(I3+M) = XOLD + VNEW*DT
               ENDDO
               I3 = I3+NDIM
            ENDDO
         ENDDO

C F : old pos
C V : old vel
C X : new free flight pos

         CALL mm_SHAKE(NATTOT,NRAM,NSM,N4DBON,WINVS,NPROT,
     $        NCONS,ICONS,JCONS,CONS,FCONS,
     $        SHKTOL,F,X,C4D,
     $        0,JRC,.TRUE.,NITS,.TRUE.)

         IF (NITS .EQ. 0)THEN
            PRINT FMSHKS,PRGSTR
            CALL EMGSTP(IFNCFG)
         ENDIF

         I3 = NPROT3
         DO NN=1,NSM
            DO J=1,NRAM
               RMASS = WMASS(J)
!$OMP parallel do private(M,VOLD,VNEW) shared(DTINV,RMASS)
               DO M=1,NDIM
                  VOLD = V(I3+M)
                  VNEW = (X(I3+M)-F(I3+M))*DTINV
                  V(I3+M) = VNEW
                  ESS(M) = ESS(M) + RMASS*VNEW*VNEW
                  EKS(M) = EKS(M) + RMASS*(VOLD+VNEW)**2
               ENDDO
               I3 = I3 + NDIM
            ENDDO
         ENDDO
      ELSE
         I3 = NPROT3
         DO NN=1,NSM
            DO J=1,NRAM
               RMINV = WINVS(J)
               DTMINV = RMINV*DT
               RMASS = WMASS(J)
!$OMP parallel do private(M,VOLD,VNEW) shared(DT,RMASS)
               DO M=1,NDIM
                  VOLD = V(I3+M)
                  VNEW = SCLS(M)*(VOLD + F(I3+M)*DTMINV)
                  V(I3+M) = VNEW
                  ESS(M) = ESS(M) + RMASS*VNEW*VNEW
                  EKS(M) = EKS(M) + RMASS*(VOLD + VNEW)**2
                  X(I3+M )= X(I3+M) + VNEW*DT
               ENDDO
               I3 = I3+NDIM
            ENDDO
         ENDDO
      ENDIF

      EKSO3  = 0.0D0
      ESCLS  = 0.0D0
!$OMP parallel do private(M) reduction(+:EKSO3,ESCLS)
      DO M=NDRMIN,NDRMAX
         EKSO3 = EKSO3 + EKS(M)
         ESCLS = ESCLS + ESS(M)
      ENDDO
      ESCLS = ESCLS*0.5D0
      EKSO3 = EKSO3*ONE8TH

      ESCL4S = ESS(MAXDIM)*0.5D0
      EKS4   = EKS(MAXDIM)*ONE8TH

C now fill my arrays to return
      EKNOW(NFTTO)  = EKPR3 + EKSO3 + EKS4 + EKP4
      EKNOW(NFTPIR) = EKPR3 - EKPCM
      EKNOW(NFTPCM) = EKPCM
      EKNOW(NFTSLV) = EKSO3
      EKNOW(NFT4D)  = EKS4 + EKP4
      EKNOW(NFTS4D) = EKS4
      EKNOW(NFTP4D) = EKP4

      EKNEXT(NFTTO) = ESCLP + ESCLS + ESCL4S + ESCL4P
      EKNEXT(NFTPIR)= ESCLP - ESPCM
      EKNEXT(NFTPCM)= ESPCM
      EKNEXT(NFTSLV)= ESCLS
      EKNEXT(NFT4D) = ESCL4S + ESCL4P
      EKNEXT(NFTS4D)= ESCL4S
      EKNEXT(NFTP4D)= ESCL4P

C reset the 3D atoms to 3D if we are running in 4D
      IF (NDIM .EQ. MAXDIM) THEN
         I3 = NDIM
         DO J=1,NATTOT
            IF (C4D(J) .LT. 0.0D0) THEN
               X(I3) = 0.0D0
               V(I3) = 0.0D0
            ENDIF
            I3 = I3 + NDIM
         ENDDO
      ENDIF
C end MDLEAP
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SDLEAP(IFNCFG,NATTOT,NPM,NRP,NSM,NRAM,N4DBON,
     $     NSPM,NSP,SUBMAS,SUBMIN,
     $     WINV,WMAS,WINVS,WMASS,DT,
     $     X,V,F,SX,C4D,VELCM,
     $     LSHKPR,NCONG,ICOG,JCOG,CONP,FCON,NRCON,JRC,SHKTOL,
     $     LSHKSO,NCONS,ICONS,JCONS,CONS,
     $     SCAL,EKNOW,EKNEXT,
     $     LEVERY,NITP,NITS,
     $     IG,CC1,CC2,CC3,CC4,CC5,CC6,CC7,CC8,CC9,
     $     SWINK,SWINKS)
C     ==--------------------------------------------------------------==
C     SDLEAP performs a stochastic dynamics integration step
C     with L<SHAKE> if specified.
C     
C     NOTE: in order to save space, the F array is used as
C     temporary storage in this routine.
C     THUS, THE FORCES ARE NO LONGER IN THE ARRAY F AFTER THIS ROUTINE.
C
C     IFNCFG  unit number onto which the configuration will be written
C             in case of a SHAKE error
C     NATTOT  the total number of atoms in the system
C     X,V,F   coordinate, velocity, force 
C     arrays, respectively.
C     NCONG   number of constraints in ICOG,JCOG and CONP
C     ICOG,JCOG atoms forming constraints
C     CONP    constraints for shake
C     SCALP   scaling factor for solute temperature bath.
C     SCALS   sclaing factor for solvent temperature bath.
C     SCAL4   scaling factor in the fourth dim.
C     
C     EKPR3   delivered with kinetic energy of solute in 3D at time T
C     EKPR4   delivered with kinetic energy of solute in 4D at time T
C     EKSO3   delivered with kinetic energy of solvent in 3D at time T
C     EKSO4   delivered with kinetic energy of solvent in 3D at time T
C     NITP    delivered with number of L<SHAKE> iterations for solute
C     NITS    delivered with number of L<SHAKE> iterations fo solvent
C     
C     THE FOLLOWING ENERGIES ARE AT TIME T+DT/2 !!
C     The above kinetic energies are used for calculating the total
C     energy at time T.
C     The following kinetic energies are used for scaling the
C     temperature in the next step, i.e. T+DT.
C     ESCLP   delivered with kinetic energy of solute
C     in 3D at time T+DT/2
C     ESCLS   delivered with kinetic energy of solvent
C     in 3D at time T+DT/2
C     ESCL4   delivered with kinetic energy in the fourth dim.
C     at time T+DT/2
C     
C     The following pseudo code might serve to elucidate the
C     code; familiarity with SD algorithms is assumed.
C     The solute and solvent parts are conceptually the same.
C     
C     In the following, VNEW, VOLD are scalar auxiliary variables,
C     not whole arrays.
C     Every "step" below corresponds to one set of loops in the
C     routine.
C     
C     initially, we have:
C     X: contains positions at T
C     F: contains forces at T
C     V: contains velocities at T-DT/2
C     
C     IF (SHAKE) THEN
C     step 1)
C        VNEW := new velocity
C        update SX array
C        F := X    (*copy old positions to F*)
C        X := X + VNEW*CVDT (* new free flight pos*)
C        (*note that V is not modified in this step*)
C     step 2)
C        CALL mm_SHAKE(F,X) --> modified X array
C     step 3)
C        VOLD := V
C        VNEW := (X - F)/CVINV
C        calculate ESCL?? (kinetic energy at time T+DT/2)
C        and EK?? (kinetic energy at time T)
C        now add random fluctuations to the new X values:
C        X := X perturbed
C        update SX array
C     step 4)
C        CALL mm_SHAKE(F,X) --> modified X array
C     ELSE
C        (*if we dont use shake, we can combine
C          step 1) and 3) into one; neither do we copy old X to F*)
C     step 4)
C        VOLD := old velocity
C        VNEW := new velocity
C        V := VNEW
C        calculate ESCL?? (kinetic energy at time T+DT/2
C        and EK?? (kinetic energy at time T)
C        X := X + VNEW*CVDT (* new free flight pos*)
C        now add random fluctuations to the new X values:
C        X := X perturbed
C        update SX array
C     ENDIF
C     ==--------------------------------------------------------------==
C     includes
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'toposz.h'
      INCLUDE 'formats.h'
      INCLUDE 'runmd.h'
C     args
      INTEGER NATTOT,NPM,NRP,NSM,NRAM,NSPM,NCONS,N4DBON
      INTEGER NITP,NITS
      INTEGER IFNCFG,IG,NSP(NSPM)
      INTEGER NCONG,ICOG(NCONG),JCOG(NCONG),NRCON,JRC(NRCON)
      INTEGER ICONS(NCONS),JCONS(NCONS)
      LOGICAL LSHKPR,LSHKSO,LEVERY
      real*8 SUBMAS(NSPM),SUBMIN(NSPM),DT,SHKTOL
      real*8 X(NATTOT*NDIM),V(NATTOT*NDIM),F(NATTOT*NDIM)
      real*8 SX(NATTOT*NDIM),C4D(NATTOT)
      real*8 CONP(NCONG),FCON(NCONG),CONS(NCONS)
      real*8 WINV(NRP),WMAS(NRP),WINVS(NRAM),WMASS(NRAM)
      real*8 VELCM(NPM*NSPM*NDIM)
      real*8 SCAL(NFTMAX),EKNOW(NFTMAX),EKNEXT(NFTMAX)
C
      real*8 CC1(NATTOT),CC2(NATTOT),CC3(NATTOT),CC4(NATTOT)
      real*8 CC5(NATTOT),CC6(NATTOT),CC7(NATTOT),CC8(NATTOT)
      real*8 CC9(NATTOT)
      real*8 SWINK(NRP),SWINKS(NRAM)
C     
C     local params
      real*8 ONE8TH
      PARAMETER (ONE8TH = 1.0D0/8.0D0)
C     local vars
      INTEGER I,J,K,M,I3,NN,IFIRST,ILAST,I3FRST,ICM3,ISAVE
      INTEGER NPROT,NPROT3,NSOLV,NSOLV3
      real*8 RMASS,RMINV,VOLD,VNEW,XOLD,XNEW,VMNEW
      real*8 CV,CVDT,CVINV,SD1,SD2,SD3,SD4,CF,CSX,SVH,CSV,SXH,SWTMP

      real*8 EKPR3,EKSO3,EKS4,EKP4,EKPCM
      real*8 ESCLP,ESCLS,ESCL4S,ESCL4P,ESPCM

      real*8 RV1(MAXDIM),RV2(MAXDIM)
      real*8 RX1(MAXDIM),RX2(MAXDIM)
      real*8 FCONS(MXCONS)
      real*8 SCLP(MAXDIM),SCLCM(MAXDIM),SCLS(MAXDIM)
      real*8 EKP(MAXDIM),ESP(MAXDIM),EKS(MAXDIM),ESS(MAXDIM)
      real*8 ESPSUB(MAXDIM),EKPSUB(MAXDIM),VELSUB(MAXDIM),VELCMN(MAXDIM)
      CHARACTER PRGSTR*(6)
C data
      DATA PRGSTR /'SDLEAP'/
C begin
      NPROT = NPM*NRP
      NPROT3= NDIM*NPROT

      NSOLV = NSM*NRAM
      NSOLV3= NDIM*NSOLV

      IF (NATTOT .NE. NPROT+NSOLV) THEN
         PRINT *,PRGSTR,': NATTOT inconsistency!'
         CALL STOPGM('MM_SETUP_DR','FATAL ERROR')
      ENDIF

!$OMP parallel do private(M)
      DO M=1,NDRMAX
         SCLP(M) = SCAL(NFTPIR)
         SCLCM(M)= SCAL(NFTPCM) - SCAL(NFTPIR)
         SCLS(M) = SCAL(NFTSLV)
      ENDDO
      SCLP(MAXDIM) = SCAL(NFT4D)
      SCLCM(MAXDIM)= 0.0D0
      SCLS(MAXDIM) = SCAL(NFT4D)
!$OMP parallel do private(M)
      DO M=1,MAXDIM
         EKP(M) = 0.0D0
         ESP(M) = 0.0D0
         ESPSUB(M) = 0.0D0
         EKPSUB(M) = 0.0D0
         EKS(M) = 0.0D0
         ESS(M) = 0.0D0
      ENDDO

      IF (LSHKPR) THEN
C---  FIRST: calc new free flight position into X
C     and calc values for SX.
C Here we DONT calc Ekin, V or VCM but copy old X to F for shaking
C the new free flight X.
C We have to overwrite F with the unscaled free flight vel in the first
C loop over sub-molecules for scaling in the second one.

         I = 0
         I3 = 0
         DO 10 NN=1,NPM
            IFIRST = 1
            DO 15 K=1,NSPM
               CALL mm_AZZERO(VELSUB,NDIM)
               ILAST = NSP(K)
               I3FRST = I3
               ISAVE = I
C calc new unscaled velocity --> F
               DO 20 J=IFIRST,ILAST
                  I = I+1
                  RMASS = WMAS(J)

                  CV = CC1(I)
                  CF = WINV(J)*DT*CC2(I)
                  SD1 = SWINK(J)*CC3(I)
                  SD2 = SWINK(J)*CC4(I)
                  CSX = CC5(I)
               
                  CALL GAUSS(0.0D0,SD1,NDIM,RV1,IG)
                  CALL GAUSS(0.0D0,SD2,NDIM,RV2,IG)
!$OMP parallel do private(M,SVH,VNEW) shared(CSX,CV,CF,RMASS,RV1,RV2,V)               
                  DO M=1,NDIM
                     SVH = SX(I3+M)*CSX + RV2(M)
                     SX(I3+M)= RV1(M)
                     VNEW = ((V(I3+M)-SVH)*CV+F(I3+M)*CF+RV1(M))
                     F(I3+M) = VNEW
                     VELSUB(M) = VELSUB(M) + RMASS*VNEW
                  ENDDO
                  I3 = I3 + NDIM
 20            CONTINUE

               RMINV = SUBMIN(K)
!$OMP parallel do private(M) shared(RMINV)               
               DO M=1,NDIM
                  VELSUB(M) = VELSUB(M)*RMINV
               ENDDO

C now scale new velocity and calc new pos
               I3 = I3FRST
               I = ISAVE
               DO 720 J=IFIRST,ILAST
                  I = I+1
                  CVDT = CC6(I)*DT
!$OMP parallel do private(M,VOLD,VNEW,XOLD) shared(CVDT)               
                  DO M=1,NDIM
                     VOLD = F(I3+M)
                     VNEW = VOLD*SCLP(M) + VELSUB(M)*SCLCM(M)
                     XOLD = X(I3+M)
                     F(I3+M)= XOLD
                     X(I3+M)= XOLD + VNEW*CVDT
                  ENDDO
                  I3 = I3+NDIM
 720           CONTINUE

               IFIRST = ILAST + 1
 15         CONTINUE
 10      CONTINUE
         
         CALL mm_SHAKE(NPROT,NRP,NPM,N4DBON,WINV,0,
     $        NCONG,ICOG,JCOG,CONP,FCON,
     $        SHKTOL,F,X,C4D,
     $        NRCON,JRC,LEVERY,NITP,.FALSE.)

         IF (NITP .EQ. 0)THEN
            PRINT FMSHKP,PRGSTR
            CALL EMGSTP(IFNCFG)
         ENDIF

C the arrays now contain:
C F  : old pos (unchanged)
C V  : old vel (unchanged)
C VCM: old com vels (unchanged)
C X  : new shaken pos
         I3 = 0
         ICM3 = 0
         I = 0
         DO 30 NN = 1,NPM
            IFIRST = 1
            DO 245 K=1,NSPM
               CALL mm_AZZERO(VELSUB,NDIM)
               ILAST = NSP(K)
               DO 40 J=IFIRST,ILAST
                  I = I + 1
                  RMASS = WMAS(J)
                  CVINV = 1.0D0/(CC6(I)*DT)
                  SD1   = SWINK(J)*CC7(I)
                  SD2   = SWINK(J)*CC8(I)
                  CSV   = CC9(I)
                  CALL GAUSS(0.0D0,SD1,NDIM,RX1,IG)
                  CALL GAUSS(0.0D0,SD2,NDIM,RX2,IG)
!$OMP parallel do private(M,VOLD,VNEW,VMNEW,SXH) shared(CVINV,CSV,RMASS)
                  DO M = 1,NDIM
                     VOLD  = V(I3+M)
                     VNEW  = (X(I3+M)-F(I3+M))*CVINV
                     V(I3+M) = VNEW
                     VMNEW     = RMASS*VNEW
                     VELSUB(M) = VELSUB(M) + VMNEW
                     ESP(M)    = ESP(M) + VMNEW*VNEW
                     EKP(M)    = EKP(M) + RMASS*(VNEW+VOLD)**2
                     SXH   = SX(I3+M)*CSV + RX2(M)
                     X(I3+M) = X(I3+M) + RX1(M) - SXH
                     SX(I3+M)= RX1(M)
                  ENDDO
                  I3 = I3+NDIM
 40            CONTINUE
C
               RMASS = SUBMAS(K)
               RMINV = SUBMIN(K)
!$OMP parallel do private(M,VOLD,VNEW) shared(RMINV,RMASS)
               DO M=1,NDIM
                  VNEW = VELSUB(M)*RMINV
                  VOLD = VELCM(ICM3+M)
                  EKPSUB(M) = EKPSUB(M) + RMASS*(VOLD+VNEW)**2
                  ESPSUB(M) = ESPSUB(M) + RMASS*VNEW*VNEW
                  VELCM(ICM3+M) = VNEW
               ENDDO

               ICM3 = ICM3 + NDIM
               IFIRST = ILAST + 1
 245        CONTINUE
 30      CONTINUE

C     shake the perturbed new coordinates
         CALL mm_SHAKE(NPROT,NRP,NPM,N4DBON,WINV,0,
     $        NCONG,ICOG,JCOG,CONP,FCON,
     $        SHKTOL,F,X,C4D,
     $        NRCON,JRC,LEVERY,NITP,.FALSE.)

         IF (NITP .EQ. 0)THEN
            PRINT FMSHKP,PRGSTR
            CALL EMGSTP(IFNCFG)
         ENDIF
      ELSE
C     no shaking of solute
         I = 0
         I3 = 0
         ICM3 = 0
         DO 50 NN=1,NPM
            IFIRST = 1
            DO 220 K=1,NSPM
!$OMP parallel do private(M)
               DO M=1,NDIM
                  VELSUB(M) = 0.0D0
                  VELCMN(M) = 0.0D0
               ENDDO
               ILAST = NSP(K)
               I3FRST = I3
               ISAVE = I
C calculate unscaled V 
               DO 60 J=IFIRST,ILAST
                  RMASS = WMAS(J)
                  I = I+1
                  CV   = CC1(I)
                  CF   = WINV(J)*DT*CC2(I)
                  SWTMP= SWINK(J)
                  SD1  = SWTMP*CC3(I)
                  SD2  = SWTMP*CC4(I)
                  CSX  = CC5(I)

                  CALL GAUSS(0.0D0,SD1,NDIM,RV1,IG)
                  CALL GAUSS(0.0D0,SD2,NDIM,RV2,IG)

!$OMP parallel do private(M,SVH,VOLD,VNEW) shared(CV,CF,RMASS)
                  DO M=1,NDIM
                     SVH = SX(I3+M)*CSX + RV2(M)
                     SX(I3+M) = RV1(M)
                     VOLD = V(I3+M)
                     VNEW = ((VOLD-SVH)*CV+F(I3+M)*CF+RV1(M))
                     F(I3+M)= VNEW
                     VELSUB(M) = VELSUB(M) + RMASS*VNEW
                  ENDDO
                  I3=I3+NDIM
 60            CONTINUE

               RMINV = SUBMIN(K)
!$OMP parallel do private(M) shared(RMINV)
               DO M=1,NDIM
                  VELSUB(M) = VELSUB(M)*RMINV
               ENDDO

               I3 = I3FRST
               I  = ISAVE
               DO 240 J=IFIRST,ILAST
                  RMASS  = WMAS(J)
                  I = I+1
                  CVDT = CC6(I)*DT
                  SWTMP= SWINK(J)
                  SD3  = SWTMP*CC7(I)
                  SD4  = SWTMP*CC8(I)
                  CSV  = CC9(I)
                  CALL GAUSS(0.0D0,SD3,NDIM,RX1,IG)
                  CALL GAUSS(0.0D0,SD4,NDIM,RX2,IG)
!$OMP parallel do private(M,VOLD,VNEW,VMNEW,XNEW,SXH) 
!$OMP+ shared(CVDT,CSV,RMASS)
                  DO M=1,NDIM
                     VOLD = F(I3+M)
                     VNEW = VOLD*SCLP(M) + VELSUB(M)*SCLCM(M)
                     VOLD = V(I3+M)
                     V(I3+M) = VNEW
                     VMNEW     = RMASS*VNEW
                     VELCMN(M) = VELCMN(M) + VMNEW
                     ESP(M)    = ESP(M) + VMNEW*VNEW
                     EKP(M)    = EKP(M) + RMASS * (VNEW+VOLD)**2
                     XNEW  = X(I3+M) +  VNEW*CVDT
                     SXH   = SX(I3+M)*CSV + RX2(M)
                     X(I3+M) = XNEW + RX1(M) - SXH
                     SX(I3+M)= RX1(M)
                  ENDDO
                  I3 = I3 + NDIM
 240           CONTINUE
CC
               RMASS = SUBMAS(K)
               RMINV = SUBMIN(K)
!$OMP parallel do private(M,VOLD,VNEW) shared(RMINV,RMASS)
               DO M=1,NDIM
                  VNEW = VELCMN(M)*RMINV
                  VOLD = VELCM(ICM3+M)
                  EKPSUB(M) = EKPSUB(M) + RMASS*(VNEW+VOLD)**2
                  ESPSUB(M) = ESPSUB(M) + RMASS* VNEW*VNEW
                  VELCM(ICM3+M) = VNEW
               ENDDO

               ICM3 = ICM3 + NDIM
               IFIRST = ILAST + 1
 220        CONTINUE
 50      CONTINUE
      ENDIF

      EKPR3 = 0.0D0
      ESCLP = 0.0D0
      ESPCM = 0.0D0
      EKPCM = 0.0D0
!$OMP parallel do private(M) reduction(+:EKPR3,ESCLP,ESPCM,EKPCM)
      DO M=NDRMIN,NDRMAX
         EKPR3 = EKPR3 + EKP(M)
         ESCLP = ESCLP + ESP(M)
         ESPCM = ESPCM + ESPSUB(M)
         EKPCM = EKPCM + EKPSUB(M)
      ENDDO
      EKPR3 = EKPR3*ONE8TH
      EKPCM = EKPCM*ONE8TH
      ESCLP = ESCLP*0.5D0
      ESPCM = ESPCM*0.5D0

      EKP4   = EKP(MAXDIM)*ONE8TH
      ESCL4P = ESP(MAXDIM)*0.5D0

C------now do solvent --------------
C     
      IF (LSHKSO) THEN
C---  FIRST: calc new free flight position into X
C     and calc values for SX.
C     copy old X position to F for shaking.
C     do NOT calculate kinetic energies, do NOT update V

         I = NPROT
         I3 = NPROT3
         DO NN=1,NSM
            DO J=1,NRAM
               I = I+1
               CV = CC1(I)
               CF = WINVS(J)*DT*CC2(I)
               SD1 = SWINKS(J)*CC3(I)
               SD2 = SWINKS(J)*CC4(I)
               CSX = CC5(I)
               CVDT = CC6(I)*DT
               CALL GAUSS(0.0D0,SD1,NDIM,RV1,IG)
               CALL GAUSS(0.0D0,SD2,NDIM,RV2,IG)
               DO M=1,NDIM
                  I3 = I3 + 1
                  SVH = SX(I3)*CSX + RV2(M)
                  SX(I3)= RV1(M)
                  VNEW = ((V(I3)-SVH)*CV+F(I3)*CF+RV1(M))*SCLS(M)
                  XOLD = X(I3)
                  F(I3)= XOLD
                  X(I3)= XOLD +  VNEW*CVDT
               ENDDO
            ENDDO
         ENDDO

         CALL mm_SHAKE(NATTOT,NRAM,NSM,N4DBON,WINVS,NPROT,
     $        NCONS,ICONS,JCONS,CONS,FCONS,
     $        SHKTOL,F,X,C4D,
     $        0,JRC,.TRUE.,NITS,.TRUE.)

         IF (NITS .EQ. 0)THEN
            PRINT FMSHKS,PRGSTR
            CALL EMGSTP(IFNCFG)
         ENDIF

         I3 = NPROT3
         I = NPROT
         DO NN = 1,NSM
            DO J = 1,NRAM
               I = I + 1
               RMASS = WMASS(J)
               CVINV = 1.0D0/(CC6(I)*DT)
               SD1   = SWINKS(J)*CC7(I)
               SD2   = SWINKS(J)*CC8(I)
               CSV   = CC9(I)
               CALL GAUSS(0.0D0,SD1,NDIM,RX1,IG)
               CALL GAUSS(0.0D0,SD2,NDIM,RX2,IG)
               DO M = 1,NDIM
                  I3 = I3 + 1
                  VOLD  = V(I3)
                  VNEW  = (X(I3)-F(I3))*CVINV
                  V(I3) = VNEW
                  ESS(M) = ESS(M) + RMASS*VNEW*VNEW
                  EKS(M) = EKS(M) + RMASS*(VOLD+VNEW)**2
                  SXH   = SX(I3)*CSV + RX2(M)
                  X(I3) = X(I3) + RX1(M) - SXH
                  SX(I3)= RX1(M)
               ENDDO
            ENDDO
         ENDDO

C     shake the perturbed new coordinates
         CALL mm_SHAKE(NATTOT,NRAM,NSM,N4DBON,WINVS,NPROT,
     $        NCONS,ICONS,JCONS,CONS,FCONS,
     $        SHKTOL,F,X,C4D,
     $        0,JRC,.TRUE.,NITS,.TRUE.)
         IF (NITS .EQ. 0)THEN
            PRINT FMSHKS,PRGSTR
            CALL EMGSTP(IFNCFG)
         ENDIF
      ELSE
C     no shaking of solvent
         I = NPROT
         I3 = NPROT3
         DO 150 NN=1,NSM
            DO 160 J=1,NRAM
               RMASS = WMASS(J)
               I = I+1
               CV   = CC1(I)
               SWTMP= SWINKS(J)
               CF   = WINVS(J)*DT*CC2(I)
               SD1  = SWTMP*CC3(I)
               SD2  = SWTMP*CC4(I)
               CSX  = CC5(I)
               CVDT = CC6(I)*DT
               SD3  = SWTMP*CC7(I)
               SD4  = SWTMP*CC8(I)
               CSV  = CC9(I)

               CALL GAUSS(0.0D0,SD1,NDIM,RV1,IG)
               CALL GAUSS(0.0D0,SD2,NDIM,RV2,IG)
               CALL GAUSS(0.0D0,SD3,NDIM,RX1,IG)
               CALL GAUSS(0.0D0,SD4,NDIM,RX2,IG)

               DO M=1,NDIM
                  I3 = I3 + 1
                  SVH = SX(I3)*CSX + RV2(M)
                  VOLD = V(I3)
                  VNEW = ((VOLD-SVH)*CV+F(I3)*CF+RV1(M))*SCLS(M)
                  V(I3)= VNEW
                  ESS(M) = ESS(M) + RMASS*VNEW*VNEW
                  EKS(M) = EKS(M) + RMASS*(VOLD+VNEW)**2
                  XNEW  = X(I3) +  VNEW*CVDT
                  SXH   = RV1(M)*CSV + RX2(M)
                  X(I3) = XNEW + RX1(M) - SXH
                  SX(I3)= RX1(M)
               ENDDO
 160        CONTINUE
 150     CONTINUE
      ENDIF

      ESCLS  = 0.0D0
      EKSO3  = 0.0D0
!$OMP parallel do private(M) reduction(+:ESCLS,EKSO3)
      DO M=NDRMIN,NDRMAX
         ESCLS = ESCLS + ESS(M)
         EKSO3 = EKSO3 + EKS(M)
      ENDDO
      ESCLS  = ESCLS*0.5D0
      EKSO3  = EKSO3*ONE8TH

      ESCL4S = ESS(MAXDIM)*0.5D0
      EKS4   = EKS(MAXDIM)*ONE8TH

C now fill my arrays to return
      EKNOW(NFTTO)  = EKPR3 + EKSO3 + EKS4 + EKP4
      EKNOW(NFTPIR) = EKPR3 - EKPCM
      EKNOW(NFTPCM) = EKPCM
      EKNOW(NFTSLV) = EKSO3
      EKNOW(NFT4D)  = EKS4 + EKP4
      EKNOW(NFTS4D) = EKS4
      EKNOW(NFTP4D) = EKP4

      EKNEXT(NFTTO) = ESCLP + ESCLS + ESCL4S + ESCL4P
      EKNEXT(NFTPIR)= ESCLP - ESPCM
      EKNEXT(NFTPCM)= ESPCM
      EKNEXT(NFTSLV)= ESCLS
      EKNEXT(NFT4D) = ESCL4S + ESCL4P
      EKNEXT(NFTS4D)= ESCL4S
      EKNEXT(NFTP4D)= ESCL4P

C     reset the atoms to 3D if necessary
      IF (NDIM .EQ. MAXDIM) THEN
         I3 = NDIM
         DO J=1,NATTOT
            IF (C4D(J) .LT. 0.0D0) THEN
               X(I3) = 0.0D0
               V(I3) = 0.0D0
               SX(I3)= 0.0D0
            ENDIF
            I3 = I3 + NDIM
         ENDDO
      ENDIF
C     end SDLEAP
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SCLCRD(NATTOT,X,XC,PRES,PRESTO)
C     ==--------------------------------------------------------------==
C     Perform pressure coupling (rescale box and coords)
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'md.h'
C args
      INTEGER NATTOT
      real*8 X(NDIM*NATTOT),XC(NDIM*NATTOT)
C PRESTO is total pressure
      real*8 PRES(NDIM),PRESTO
C local params
      real*8 ONETHI
      PARAMETER (ONETHI = 1.0D0/3.0D0)
C local vars
      INTEGER M,J,I3
      real*8 DTCP,RMUFAC(MAXDIM)
C begin
      DTCP = COMP*DT/TAUP
      IF (NTP .EQ. NTPISO) THEN
         RMUFAC(1) = (1.0D0-DTCP* (PRES0-PRESTO))**ONETHI
         RMUFAC(2) = RMUFAC(1)
         RMUFAC(3) = RMUFAC(1)
      ELSEIF (NTP .EQ. NTPANI) THEN
!$OMP parallel do private(M)
         DO M = 1,NDIM
            RMUFAC(M) = (1.0D0-DTCP* (PRES0-PRES(M)))**ONETHI
         ENDDO
      ENDIF 

!$OMP parallel do private(M)
      DO M = 1,NDIM
         BOX(M)   = BOX(M)*RMUFAC(M)
         BOXH(M)  = BOX(M)*0.5D0
         BOXINV(M)= 1.0D0/BOX(M)
      ENDDO

      I3 = 0
      DO J = 1,NATTOT
!$OMP parallel do private(M)
         DO M = 1,NDIM
            X(I3+M) = X(I3+M)*RMUFAC(M)
         ENDDO
         I3 = I3 + NDIM
      ENDDO

C if we have position re(con)straints: rescale reference positions
      IF (NTR .NE. NTROFF) THEN
         I3 = 0
         DO J=1,NATTOT
!$OMP parallel do private(M)
            DO M=1,NDIM
               XC(I3+M) = XC(I3+M)*RMUFAC(M)
            ENDDO
            I3 = I3 + NDIM
         ENDDO
      ENDIF
C end SCLCRD
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SETPRT(NRP,NPM,NSPM,RLAM,DLAM,
     $     WMAS,WINV,TOTMAS,NSPG,SUBMAS,SUBMIN,
     $     LDOSD,BOLTZT,SWINK,
     $     LDOSHK,LSHDOB,NSIZE,CONP)
C     ----------------------------------------------------------------== 
C     (1)    adjust the masses according to RLAM
C            If SD is performed, the SWINK array is adjusted
C            accordingly.
C
C     (2)    Update the constraint length according to RLAM
C
C     TOTMAS   delivered with total mass of the system
C     (all solute and solvent atoms)
C     SUBMAS(NSPM) delivered with total mass of each submolecule
C     SUBMIN(NSPM) delivered with inverse total mass of each submolecule
C
C     If all atoms in a submolecule I have mass zero (due to position
C     constraining), then SUBMAS (I) and SUBMIN(I) are both set to 0.0.
C     ----------------------------------------------------------------== 
C includes
      INCLUDE 'pertsz.h'
      INCLUDE 'pertar.h'
C args
      INTEGER NRP,NPM,NSPM,NSPG(NJLA),NSIZE,ISUB
      LOGICAL LDOSD,LDOSHK,LSHDOB
      real*8 RLAM,DLAM,TOTMAS,BOLTZT
      real*8 WMAS(NRP),WINV(NRP),SWINK(NRP)
      real*8 SUBMAS(NSPM),SUBMIN(NSPM)
      real*8 CONP(NSIZE)
C local vars
      INTEGER JL,J,NN,N
      real*8 DW
C update the mass arrays to reflect the new lamda value
C also update SWINK if we are doing SD.
      DO JL = 1,NJLA
         J = JLA(JL)
         WMAS(J) = (1.0D0-RLAM)*WMA(JL)+RLAM*WMB(JL)
         WINV(J) = 1.0D0/WMAS(J)
         IF (LDOSD) THEN
            SWINK(J) = DSQRT(BOLTZT*WINV(J))
         ENDIF
         DW     = DLAM* (WMB(JL)-WMA(JL))
         TOTMAS = TOTMAS + DW*NPM
         ISUB = NSPG(JL)
         SUBMAS(ISUB) = SUBMAS(ISUB) + DW
         IF (SUBMAS(ISUB) .GT. 0.0D0) THEN
            SUBMIN(ISUB) = 1.0D0/SUBMAS(ISUB)
         ELSE
            SUBMIN(ISUB) = 0.0D0
         ENDIF
      ENDDO

C     set distance constraint (bond) lengths
      IF (LDOSHK) THEN
         NN = 0
         DO N = 1,NBONHG
            NN = NN + 1
C     Note this is a linear dependence of the constraint distance on lambda 
            CONP(NN) = ((1.0D0-RLAM)*BHA0(N)+RLAM*BHB0(N))**2
         ENDDO

         IF (LSHDOB) THEN
            DO N = 1,NBONG
               NN = NN + 1
C     Note this is a linear dependence of the constraint distance on lambda
               CONP(NN) = ((1.0D0-RLAM)*BA0(N)+RLAM*BB0(N))**2
            ENDDO
         ENDIF
      ENDIF
C setprt
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE CLCPRT(NATTOT,NRP,NPM,NSPM,NSP,
     $     DT,RLAM,DEDLAM,
     $     V,VCM,
     $     LDOSHK,LSHDOB,NSIZE,FCON)
C     ----------------------------------------------------------------== 
C     (1) Calculate the kinetic energy contribution to the free energy
C         derivative.
C
C     (3) Calculate the constraint contribution to the free energy.
C         using formula from eq 4.5.9 in (MD 93.29)
C
C     DEDLAM  is delivered with the derivatives of the various
C     terms of the Hamiltionian with respect to L<LAMDA>
C     ----------------------------------------------------------------== 
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'pertsz.h'
      INCLUDE 'pertar.h'
C args
      INTEGER NATTOT,NRP,NPM,NSPM,NSIZE,NSP(NSPM)
      LOGICAL LDOSHK,LSHDOB
      real*8 DT,RLAM
      real*8 DEDLAM(MXETBL)
      real*8 V(NDIM*NATTOT),VCM(NDIM*NPM*NSPM)
      real*8 FCON(NSIZE)
C local vars
      INTEGER NN,M,JNN,JL,J3,N,ICM3,ITOPA,IFIRST,ILAST,ISUB
      real*8 DT2INV,DMASS,CURLEN
      real*8 EKP3TO,EKP4TH,EKPCM
C
C     EKP3TO : total (internal + rotational + com translational)
C     kinetic energy in 3 dimensions
C     EKPCM  : com translational kinetic energy in 3 dimensions
C     EKP4TH : kinetic energy in the fourth dimension
C     
C begin
C     perturbation contributions:
C     calculate the derivative of the kinetic energy
C     with respect to lamda for solute internal-rotational and
C     translational centre of mass.
      EKP3TO = 0.0D0
      EKPCM = 0.0D0
      EKP4TH = 0.0D0
      IF (NJLA .GT. 0) THEN
         ICM3 = 0
         DO NN = 1,NPM
C     centre of mass contribution
            IFIRST = 1
            DO ISUB=1,NSPM
               ILAST = NSP(ISUB)
               DO JL = 1,NJLA
                  ITOPA = JLA(JL)
                  IF (ITOPA .GE. IFIRST .AND. ITOPA .LE. ILAST) THEN
                     DMASS = (WMA(JL)-WMB(JL))
                     DO M = NDRMIN,NDRMAX
                        EKPCM = EKPCM + DMASS* VCM(ICM3+M)**2
                     ENDDO
                  ENDIF
               ENDDO
               IFIRST = ILAST + 1
               ICM3 = ICM3 + NDIM
            ENDDO

C     no calculate the total atomic contribution
            JNN = (NN-1)*NRP
            DO JL = 1,NJLA
               DMASS = (WMA(JL)-WMB(JL))
               J3 = NDIM* (JLA(JL)+JNN-1)
               DO M = NDRMIN,NDRMAX
                  EKP3TO = EKP3TO + DMASS*V(J3+M)**2
               ENDDO
               DO M = NDHMIN,NDHMAX
                  EKP4TH = EKP4TH + DMASS*V(J3+M)**2
               ENDDO
            ENDDO
         ENDDO

         EKPCM = EKPCM*0.5D0
         EKP3TO = EKP3TO*0.5D0
         EKP4TH = EKP4TH*0.5D0
      ENDIF

      DEDLAM(IKSLU3) = EKP3TO
      DEDLAM(IKSUCM) = EKPCM
      DEDLAM(IKSLU4) = EKP3TO + EKP4TH

C bond contributions
      IF (LDOSHK) THEN
C inverse of the square of DT
         DT2INV = 1.0D0/(DT*DT)
         DEDLAM(IPBH) = 0.0D0 
         NN = 0
         DO N = 1,NBONHG
            NN = NN + 1
C  linear dependence from eq 4.5.9 MD 93.29
C  NOTE FCON = 2 * DT**2 * <l>
            CURLEN = ((1.0D0-RLAM)*BHA0(N)+RLAM*BHB0(N))
            DEDLAM(IPBH) = DEDLAM(IPBH)
     $           +  DT2INV * FCON(NN)*(BHB0(N)-BHA0(N))
     $           * CURLEN
         ENDDO

C add constraint contrib. DEDLAM(IPBH) to totals 
C     DEDLAM(IPTOT),DEDLAM(IPBTOT),DEDLAM(IPHTOT),DEDLAM(IPCOVT)

         DEDLAM(IPTOT) = DEDLAM(IPTOT)  + DEDLAM(IPBH) 
         DEDLAM(IPBTOT)= DEDLAM(IPBTOT) + DEDLAM(IPBH) 
         DEDLAM(IPHTOT)= DEDLAM(IPHTOT) + DEDLAM(IPBH) 
         DEDLAM(IPCOVT)= DEDLAM(IPCOVT) + DEDLAM(IPBH) 

         IF (LSHDOB) THEN
            DEDLAM(IPBN) = 0.0D0
            DO N = 1,NBONG
               NN = NN + 1
C  linear dependence from eq 4.5.9 MD 93.29
C  NOTE FCON = 2 * DT**2 * <l>
               CURLEN = ((1.0D0-RLAM)*BA0(N)+RLAM*BB0(N))
               DEDLAM(IPBN) = DEDLAM(IPBN)
     $              +  DT2INV * FCON(NN)*(BB0(N)-BA0(N))
     $              * CURLEN
            ENDDO

C     add constraint contrib. DEDLAM(IPBN) to totals 
C     DEDLAM(IPTOT),DEDLAM(IPBTOT),DEDLAM(IPNHTO),DEDLAM(IPCOVT)
            DEDLAM(IPTOT) = DEDLAM(IPTOT)  + DEDLAM(IPBN)
            DEDLAM(IPBTOT)= DEDLAM(IPBTOT) + DEDLAM(IPBN)
            DEDLAM(IPNHTO)= DEDLAM(IPNHTO) + DEDLAM(IPBN)
            DEDLAM(IPCOVT)= DEDLAM(IPCOVT) + DEDLAM(IPBN)
         ENDIF
      ENDIF
C     end pertst
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PRNRG(NSIZE,ERG)
C     ----------------------------------------------------------------==
C     PRNRG prints out the energies in ERG to screen
C     real*8 ERG(NSIZE) the energies to print out.
C     this routine uses format strings that are
C     calculated in routine L<CHKNRE>.
C     ----------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'md.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'forcear.h'
C args
      INTEGER NSIZE
      real*8 ERG(NSIZE)
C local vars
      INTEGER ICOLNM,IOFF,II,IR
C begin
C write title string along the top
      PRINT '(1X,A12,20(1X,A12))',
     $     'ATOM NRS',(EPTSTR(II), II=1,NUSNRE)

      ICOLNM = NUSNRE
      IOFF = 0
      DO IR = 1, NUSNRE
         PRINT EPTFMT(IR), EPTSTR(IR),(ERG(IOFF+II),II=1,ICOLNM)
         IOFF = IOFF+ICOLNM
         ICOLNM = ICOLNM-1
      ENDDO
C end prnrg
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE CLCMAS(NPM,NSM,NSPM,NSP,TOTMAS,TMASS,SUBMAS,SUBMIN)
C     ==--------------------------------------------------------------==
C     CLCMAS calculates masses from the topology (WMAS, WINV arrays
C     for the solute and WINVS and WMASS arrays for the solvent)
C
C     NPM    number of solute molecules
C     NSM    number of lolvent molecules
C     NSPM   number of solute submolecules
C     NSP(NSPM)
C            index number of last atom in each submolecule
C
C
C     TOTMAS   delivered with the total mass of the system
C            (all solute and solvent atoms)
C     TMASS  delivered with the total mass of ONE solvent
C            molecule
C     SUBMAS(NSPM)
C            delivered with the total mass of each solute submolecule
C     SUBMIN(NSPM)
C            delivered with the inverse of the total mass of
C            each solute submolecule
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
C local parsms
      real*8 EPS
      PARAMETER (EPS = 1.0D-5)
C args
      INTEGER NPM,NSM,NSPM,NSP(NSPM)
      real*8 TOTMAS,TMASS
      real*8 SUBMAS(NSPM),SUBMIN(NSPM)
C local vars
      real*8 SUB
      INTEGER I,J,NSUB,IFRST,ILAST
C begin
      TOTMAS = 0.0D0
C calc mass of each submolecule individually     
      IFRST = 1
      DO NSUB=1,NSPM
         SUB = 0.0D0
         ILAST = NSP(NSUB)
C while NSP(J .LE. NSP(NSUB) do
!$OMP parallel do private(J) reduction(+:SUB)
         DO J=IFRST,ILAST
            SUB = SUB + WMAS(J)
         ENDDO
         IFRST = ILAST + 1
         SUBMAS(NSUB) = SUB
         IF (SUB .GT. EPS) THEN
            SUBMIN(NSUB) = 1.0D0/SUB
         ELSE
            SUBMIN(NSUB) = 0.0D0
         ENDIF
         TOTMAS = TOTMAS + SUB
      ENDDO

      TOTMAS = TOTMAS*NPM

C calc mass of one solvent molecule
      TMASS = 0.0D0
!$OMP parallel do private(I) reduction(+:TMASS)
      DO I=1,NRAM
         TMASS = TMASS + WMASS(I)
      ENDDO

      TOTMAS  = TOTMAS + NSM*TMASS
C end clcmas
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE CLCTEM(NATTOT,NPM,NSM,NSPM,NSP,SUBMIN,V,EKRES,VELCM)
C     ==--------------------------------------------------------------==
C     CLCTEM calculates the kinetic energy of solute and solvent
C     in 3D, and the total 4th D contribution from the velocity 
C     array. The velocity of the centres of mass of the solute
C     submolecules is also calculated from V and returned in VELCM.
C
C     NATTOT   the total number of atoms
C     NPM      the number of solute molecules
C     NSM      the number of solvent molecules
C     NSPM     the number of submolecules per solute molecule
C     NSP      array defining the submolecules
C     SUBMIN   inverse mass of the submolecules
C     V        velocity array
C     EKRES(NFTMAX)
C              delivered with the kinetic energies.
C              for the size of EKRES see runmd.h
C     VELCM(NPM*NSPM*NDIM)
C              delivered with the velocity of the centres of mass 
C              of the solute submolecules.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'runmd.h'
C args
      INTEGER NATTOT,NPM,NSM,NSPM,NSP(NSPM)
      real*8 V(NATTOT*NDIM),EKRES(NFTMAX)
      real*8 SUBMIN(NSPM)
      real*8 VELCM(NPM*NSPM*NDIM)
C local vars
      INTEGER J,NN,I3,M,K,IFRST,ILAST,ICM3
      real*8 EKLOCP(MAXDIM),EKLOCS(MAXDIM),EKPSUB(MAXDIM)
      real*8 VELSUB(MAXDIM),VEL,VM,RMINV
      real*8 TMPMSS,EKP,EKS,EK4,EKPCM
C begin
!$OMP parallel do private(M)
      DO M=1,MAXDIM
         EKPSUB(M) = 0.0D0
         EKLOCP(M) = 0.0D0
         EKLOCS(M) = 0.0D0
      ENDDO

      CALL mm_AZZERO(EKRES,NFTMAX)

      I3 = 0
      ICM3 = 0
      DO 10 NN=1,NPM
         IFRST = 1
         DO 20 K=1,NSPM
!$OMP parallel do private(M)
            DO M=1,NDIM
               VELSUB(M) = 0.0D0
            ENDDO
            ILAST = NSP(K)
            DO J=IFRST,ILAST
               TMPMSS = WMAS(J)
!$OMP parallel do private(M,VEL,VM) shared(TMPMSS)
               DO M = 1,NDIM
                  VEL = V(I3+M)
                  VM = VEL*TMPMSS
                  VELSUB(M) = VELSUB(M) + VM
                  EKLOCP(M) = EKLOCP(M) + VM*VEL
               ENDDO
               I3 = I3 + NDIM
            ENDDO
            IFRST = ILAST+1
C
            RMINV = SUBMIN(K)
!$OMP parallel do private(M,VM) shared(RMINV)
            DO M=1,NDRMAX
               VM = RMINV*VELSUB(M)
               VELCM(ICM3+M) = VM
               EKPSUB(M)     = EKPSUB(M) + VM*VELSUB(M)
            ENDDO

            ICM3 = ICM3 + NDIM
 20      CONTINUE
 10   CONTINUE

      DO NN=1,NSM
         DO J=1,NRAM
            TMPMSS = WMASS(J)
!$OMP parallel do private(M) shared(TMPMSS)
            DO M=1,NDIM
               EKLOCS(M) = EKLOCS(M) + TMPMSS*V(I3+M)**2
            ENDDO
            I3 = I3+NDIM
         ENDDO
      ENDDO

      EKP   = 0.0D0
      EKPCM = 0.0D0
      EKS   = 0.0D0
!$OMP parallel do private(M) reduction(+:EKP,EKPCM,EKS)
      DO M=1,NDRMAX
         EKP   = EKP   + EKLOCP(M)
         EKPCM = EKPCM + EKPSUB(M)
         EKS   = EKS   + EKLOCS(M)
      ENDDO
      EKP   = EKP * 0.5D0
      EKPCM = EKPCM*0.5D0
      EKS   = EKS * 0.5D0
      EK4 = (EKLOCP(MAXDIM) + EKLOCS(MAXDIM))*0.5D0

      EKRES(NFTTO) = EKP + EKS + EK4
      EKRES(NFTPIR)= EKP-EKPCM
      EKRES(NFTPCM)= EKPCM
      EKRES(NFTSLV)= EKS
      EKRES(NFT4D) = EK4
C end clctem
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PRPVIR(NATTOT,NPM,NSM,X,V,XR,TMASS,NSPM,NSP,SUBMAS,
     $     EKCM,EKCMTO,LEVERY)
C     ==--------------------------------------------------------------==
C     PRPVIR does some preparatory calculations for calculating
C     the virial.
C     Essentially, what we do in this routine is calculate the
C     centre of masses for later use.
C     After the force routines, L<CLCVIR> is called
C     by L<RUNMD> to calculate the (molecular) virial.
C     ==--------------------------------------------------------------==
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
      INTEGER NSKIPA,NSKPA3,NSKIPT
      real*8 XREF(MAXDIM),TM
      real*8 XCM(MAXDIM),VCM(MAXDIM),ACMDUM(MAXDIM),
     c       OCMDUM(MAXDIM),EKRDUM
      real*8 ESCLCM
C begin
      LSCLCM = .FALSE.
      EKCMTO = 0.0D0
      ESCLCM = 0.0D0

      NATTO3 = NDIM*NATTOT
      CALL mm_AZZERO(EKCM,NDIM)
      CALL mm_AZZERO(XR,NATTO3)

C NSKIPA : number of atoms to skip in the coord array
C NSKIPT : number of atoms to skip in the topology(WINV array)
      NSKIPA = 0
      NSKPA3 = 0
      DO 30 NNP=1,NPM
         NFIRST = 1
         NSKIPT = 0
         DO 40 NSUB=1,NSPM
            NLAST = NSP(NSUB)
            NUMAT = NLAST - NFIRST + 1
C XREF is always the first atom in the current submolecule
!$OMP parallel do private(M)
            DO M=1,NDIM
               XREF(M) = XR(NSKPA3+M)
            ENDDO
            CALL CONAT(NATTOT,NUMAT,NSKIPA,XREF,XR,LEVERY)
            TM = SUBMAS(NSUB)
!           CALL CENMAS(NATTOT,1,NUMAT,0,0,NSKIPA,NDIM,NDRMAX,XR,V,TM,
!    $           NSKIPT,WMAS,WMASS,
!    $           EKCMTO,XCM,VCM,ACMDUM,EKRDUM,OCMDUM, -ICMVEL)

!$OMP parallel do private(M) shared(TM,VCM)
            DO M=1,NDIM
               EKCM(M) = EKCM(M) + TM*VCM(M)**2
            ENDDO
            DO N=1,NUMAT
!$OMP parallel do private(M) shared(XCM)
               DO M=1, NDIM
                  XR(NSKPA3+M) = XR(NSKPA3+M) - XCM(M)
               ENDDO
               NSKPA3 = NSKPA3 + NDIM
            ENDDO
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
      DO NNP=1,NSM
!        CALL CENMAS(NATTOT,0,0,1,NRAM,NSKIPA,NDIM,NDRMAX,XR,V,
!    $        TMASS,NSKIPT,WMAS,WMASS,
!    $        EKCMTO,XCM,VCM,ACMDUM,EKRDUM,OCMDUM, -ICMVEL)

!$OMP parallel do private(M) shared(TMASS,VCM)
         DO M=1,NDIM
            EKCM(M) = EKCM(M) + TMASS*VCM(M)**2
         ENDDO

         DO N=1,NRAM
!$OMP parallel do private(M) shared(XCM)
            DO M=1,NDIM
               XR(NSKPA3+M) = XR(NSKPA3+M) - XCM(M)
            ENDDO
            NSKPA3 = NSKPA3 + NDIM
         ENDDO
         NSKIPA = NSKIPA + NRAM
      ENDDO
C end PRPVIR
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE CLCVIR(EKCM,VIR,PRES,EKCMTO,VIRTOT,PRESTO)
C     ==--------------------------------------------------------------==
C     CLCVIR calculates the molecular virial and pressure.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
C args
      real*8 EKCM(NDIM),VIR(NDIM),PRES(NDIM)
      real*8 EKCMTO,VIRTOT,PRESTO
C local vars
      LOGICAL LOCTO
      INTEGER M
      real*8 VOL
C begin
      IF (NDIM .NE. 3) THEN
         PRINT *,'CLCVIR: implementation restriction!'
         PRINT *,'NDIM must be three !'
         CALL STOPGM('MM_SETUP_DR','FATAL ERROR')
      ENDIF
      LOCTO = (NTB .LT. 0)

      VOL = BOX(1)*BOX(2)*BOX(3)
      IF (LOCTO) THEN
         VOL = VOL*0.5D0
      ENDIF

      EKCMTO = 0.0D0
      VIRTOT = 0.0D0
      PRESTO = 0.0D0
!$OMP parallel do private(M) 
!$OMP+ reduction(+:EKCMTO,VIRTOT,PRESTO)
      DO M = 1,NDIM
         EKCM(M) = EKCM(M)*0.5D0
         PRES(M) = 2.0D0* (EKCM(M)-VIR(M))/VOL
         EKCMTO  = EKCMTO + EKCM(M)
         VIRTOT  = VIRTOT + VIR(M)
         PRESTO  = PRESTO + PRES(M)
      ENDDO
      PRESTO = PRESTO/3.0D0
C end CLCVIR
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PRENER(LDOKIN,EE,LDORES,EERES,LDOPIN,LDOVIR,
     $     LDOTMP,VOLPRT,TEMP) 
C     ----------------------------------------------------------------==
C     PRENER prints out the EE, EERES, VOLPRT and TEMP
C     arrays in a standard format according to the switches passed in.
C
C     LDOKIN       print out kinetic data in array EE if .TRUE.
C     EE(MXETBL)   the array to print out.
C     LDORES       print out the EERES array if .TRUE.
C     EERES(MXCTBL)additional array to print out if LDORES = .TRUE.
C     LDOPIN       print out path integral information
C     LDOVIR       print out virial and pressure information
C     LDOTMP       print out temperature information
C     ----------------------------------------------------------------==
      INCLUDE 'coordsz.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'runmd.h'
C args
      LOGICAL LDOKIN,LDORES,LDOPIN,LDOVIR,LDOTMP
      real*8 EE(MXETBL),EERES(MXCTBL)
      real*8 VOLPRT(MXVTBL),TEMP(NFTMAX)
C local vars
      INTEGER II
C begin
      IF (LDOKIN) THEN
         PRINT 500,'ENERGY','TOTAL','KINETIC','POTENTIAL'
         PRINT 550,'E-TOT',EE(IETOT),EE(IKTOT4),EE(IPTOT)
         PRINT *

         PRINT 500,'KINETIC','TOTAL','SOLUTE-IR-3D','SOLUTE-COM-3D',
     $        'SOLVENT-3D','4th D'
         PRINT 550,'E-KIN',EE(IKTOT4),EE(IKSLU3)-EE(IKSUCM),
     $        EE(IKSUCM),EE(IKSLV3),EE(IK4THD)
         IF (LDOTMP) THEN
            PRINT 550,'TEMP',TEMP(NFTTO),TEMP(NFTPIR),TEMP(NFTPCM),
     $           TEMP(NFTSLV),TEMP(NFT4D)
            PRINT 560,'SCALE',VOLPRT(IVSPIR),VOLPRT(IVSPCM),
     $           VOLPRT(IVSCLS),VOLPRT(IVSCL4)
            PRINT *
         ENDIF
      ELSE
         PRINT 500,'ENERGY'
         PRINT 550,'E-TOT',EE(IETOT)
      ENDIF

C now potential energy
      PRINT 500,'POTENTIAL','TOTAL','COVALENT','NON-BONDED','SPECIAL'
C     PRINT 550,'E-POT',EE(IPTOT),EE(IPCOVT),ENOBON,ESPEC
      PRINT 550,'E-POT',EE(IPTOT),EE(IPCOVT),EE(IPNBON),EE(IPSPEC)
      PRINT *

C components
      PRINT 500,'COVALENT','TOTAL','H-ATOMS','NO H-ATOMS'
      PRINT 550,'E-BOND',EE(IPBTOT),EE(IPBH),EE(IPBN)
      PRINT 550,'E-BOND-ANGLE',EE(IPBATO),EE(IPBAH),EE(IPBA)
      PRINT 550,'E-IMPR.DIH',EE(IPQTOT),EE(IPQH),EE(IPQE)
      PRINT 550,'E-DIHEDRAL',EE(IPPTOT),EE(IPPH),EE(IPPE)
      PRINT 550,'E-COVALENT',EE(IPCOVT),EE(IPHTOT),EE(IPNHTO)
      PRINT *
 
#ifdef EWALD
      PRINT 500,'NON-BONDED','TOTAL','R-SPACE','EXCL-CORR',
     $     'SELF-TERM','K-SPACE'
      PRINT 550,'E-ELECTRO',EE(IPELEC),EE(IPEL),EE(IPRF),EE(IPRC),
     $     EE(IPKW)
#else
      PRINT 500,'NON-BONDED','TOTAL','SHORT+LONG','RF-DIST-DEP',
     $     'RF-DST-INDEP'
      PRINT 550,'E-ELECTRO',EE(IPELEC),EE(IPEL),EE(IPRF),EE(IPRC)
#endif
      PRINT 550,'E-VDWAALS',EE(IPLJ),EE(IPLJ)
      PRINT 550,'E-NONBON',EE(IPNBON)
      PRINT *

      IF (LDORES) THEN
         PRINT 500,'SPECIAL','TOTAL','E-POSRE','E-DISRE','E-DIHRE',
     $        'E-JVAL'
         PRINT 550,'E-SPECIAL1',EE(IPSPEC),EERES(ICPOSR),EERES(ICDISR),
     $        EERES(ICDHRE),EERES(ICJVAL)

         PRINT 500,' ',' ','E-LOCELEV','E-4D-REST'
         PRINT 560,'E-SPECIAL2',EERES(ICLOCE),EERES(ICRE4D)
         PRINT *
      ENDIF

      IF (LDOPIN) THEN
         PRINT 500,'PATH-INT','RING'
         PRINT 550,'E-RING',EE(IPPISP)
         PRINT *
      ENDIF

      IF (LDOTMP .OR. LDOVIR) THEN
         PRINT 500,'VOLPRT','PROD/AVERAGE','X-AXIS','Y-AXIS','Z-AXIS'

         IF (LDOTMP) THEN
            PRINT 550,'VOL/LENGTH',VOLPRT(IVVOL),
     $           (VOLPRT(II),II=IVBOXX,IVBOXZ)
         ENDIF
         IF (LDOVIR) THEN
            PRINT 550,'E-KIN-TR',VOLPRT(IVEKCM),
     $           (VOLPRT(II),II=IVCMX,IVCMZ)
            PRINT 550,'VIRIAL',VOLPRT(IVVIR),
     $           (VOLPRT(II),II=IVVIRX,IVVIRZ)
            PRINT 550,'PRESSURE',VOLPRT(IVPRES),
     $           (VOLPRT(II),II=IVPRSX,IVPRSZ)
         ENDIF
         PRINT *
      ENDIF
C formats
 500  FORMAT (6(1X,A12))

#ifdef HIGH_PREC_OUTPUT
C Phil
 550  FORMAT (1X,A12,5(1X,E16.9))
#else
 550  FORMAT (1X,A12,5(1X,E12.5))
#endif
 560  FORMAT (1X,A12,13X,5(1X,E12.5))
 501  FORMAT (1X,A12,5(1X,E12.5))
 502  FORMAT (1X,A12,2(13X),1X,E12.5)
 510  FORMAT (9(1X,E12.5))
C end PRENER
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PRTEMP(TP)
C     ----------------------------------------------------------------==
C     PRTEMP prints out the temperature arrsy TP.
C     It is called by L<RUNMD>.
C     ----------------------------------------------------------------==
      INCLUDE 'runmd.h'
C args
      real*8 TP(NFTMAX)
C begin
      PRINT 500,'TOTAL',TP(NFTTO)
      PRINT 500,'SOLUTE INTERNAL AND ROTATIONAL',TP(NFTPIR)
      PRINT 500,'SOLUTE C.O.M. TRANSLATIONAL   ',TP(NFTPCM) 
      PRINT 500,'SOLVENT ',TP(NFTSLV)
      PRINT 500,'4th DIMENSION',TP(NFT4D)
      PRINT *

 500  FORMAT (A32,(1X,E14.8))
C end prtemp
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE CLCAVE(NVAL,MAX,ESUM,ESUMSQ,EAVE,EFLUCT)
C     ----------------------------------------------------------------==
C     CLCAVE calculates the averages and mean square
C     fluctuations of values in ESUM and ESUMSQ which are DOUBLE PRECISION
C     into EAVE and EFLUCT which are REAL.
C     ESUM,ESUMSQ,EAVE and EFLUCT are arrays of length MAX.
C
C     NVAL   the number of samples taken.
C     MAX    the size of the arrays
C     ESUM(MAX)
C            input with the sum of the samples of each value.
C            This array is not modified in this routine.
C            The array is in DOUBLE PRECISION.
C     ESUMSQ(MAX)
C            input with the sum of the square of each sample
C            for each value.
C            This array is not modified in this routine.
C            The array is in DOUBLE PRECISION.
C
C     EAVE(MAX)
C            delivered with the average of each value.
C            This array is REAL*8.
C     EFLUCT(MAX)
C            delivered with the mean square fluctuations
C            of each value.
C            This array is REAL*8.
C     ----------------------------------------------------------------==
C args
      INTEGER NVAL,MAX
      DOUBLE PRECISION ESUM(MAX),ESUMSQ(MAX)
      real*8 EAVE(MAX),EFLUCT(MAX)
C local vars
      INTEGER I
      DOUBLE PRECISION VALINV,AVE,AVESQ
C begin
      IF (NVAL .LE. 0) THEN
         PRINT *,'CLCAVE: ERROR ! averages of 0 steps!!'
         CALL STOPGM('MM_SETUP_DR','FATAL ERROR')
      ELSE
cmb         VALINV = NVAL
cmb         VALINV = 1.0D0/VALINV
         VALINV = 1.0D0/DBLE(NVAL)
      ENDIF

      DO I=1,MAX
         AVE   = ESUM(I)*VALINV
         AVESQ = ESUMSQ(I)*VALINV
         EAVE(I) = AVE
C due to numeric error, AVE**2 can become <= AVESQ
         AVE = AVE*AVE
         IF (AVE .GE. AVESQ) THEN
            EFLUCT(I) = 0.0D0
         ELSE
            EFLUCT(I) = DSQRT(AVESQ-AVE)
         ENDIF
      ENDDO
C end clcave
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PRARR(LDOKIN,LDOVIR,LDOPIN,
     $     LDOENE,ENER,ENERES,VOLPRT,TEMP,
     $     LDODLA,DEDLAM,
     $     LDOVDL,DVDL,
     $     LDODMU,DEDMU,
     $     LDOVDM,DVDM,
     $     NUSNRE,NSIZE,EPEL,EPRF,EPRC,EPLJ)
C     ----------------------------------------------------------------==
C     PRARR prints out a collection of data to screen in a standard
C     format.
C     Exactly what gets printed out is defined by the controlling 
C     logicals LDOENE,LDOKIN, LDOVIR, LDOPIN.
C
C     The individual arrays DEDLAM, DVDL, DEDMU and DVDM are
C     only printed if the corresponding logicals are .TRUE.
C
C     The arrays EPEL,EPRF,EPRC and EPLJ are printed out by calling
C     subroutine L<PRNRG> if NUSNRE is > 0.
C     ----------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'runmd.h'
      INCLUDE 'forcesz.h'
C args
      INTEGER NSIZE,NUSNRE
      LOGICAL LDOKIN,LDOVIR,LDOPIN,LDOENE
      LOGICAL LDODLA,LDOVDL,LDODMU,LDOVDM
      real*8 ENER(MXETBL),ENERES(MXCTBL)
      real*8 VOLPRT(MXVTBL),TEMP(NFTMAX)
      real*8 DEDLAM(MXETBL),DVDL(MXETBL)
      real*8 DEDMU(MXETBL),DVDM(MXETBL)
      real*8 EPLJ(NSIZE),EPEL(NSIZE),EPRF(NSIZE),EPRC(NSIZE)
C local vars
      real*8 EDUM
C begin

      IF (LDOENE) THEN
         CALL PRENER(LDOKIN,ENER,.TRUE.,ENERES,LDOPIN,LDOVIR,
     $        .TRUE.,VOLPRT,TEMP)
         PRINT *

C nonbonded energy components
         IF (NUSNRE .GT. 0) THEN
#ifdef EWALD
            PRINT *,'R-SPACE ELECTROSTATIC INTERACTION ENERGY'
#else
            PRINT *,'ELECTROSTATIC INTERACTION ENERGY'
#endif
            CALL PRNRG(NSIZE,EPEL)
         
#ifdef EWALD
            PRINT *,'EXCLUSION K-SPACE CORRECTION ENERGY'
#else
            PRINT *,'DISTANCE DEPENDENT REACTION FIELD ENERGY'
#endif
            CALL PRNRG(NSIZE,EPRF)
         
#ifdef EWALD
            PRINT *,'K-SPACE SELF TERM ENERGY'
#else
            PRINT *,'DISTANCE INDEPENDENT REACTION FIELD ENERGY'
#endif
            CALL PRNRG(NSIZE,EPRC)
         
            PRINT *,'LJ INTERACTION ENERGY'
            CALL PRNRG(NSIZE,EPLJ)
            PRINT *
         ENDIF
      ENDIF

      IF (LDODLA .OR. LDOVDL) THEN
         PRINT FMNDIM,'RLAM',VOLPRT(IVRLAM)
      ENDIF
      IF (LDODMU .OR. LDOVDM) THEN
         PRINT FMNDIM,'RMU',VOLPRT(IVRMU)
      ENDIF
#if 0
      IF (LDODLA) THEN
         PRINT *,'FREEENERGYDERIVS (dE/dLamda)'
         CALL PRENER(LDOKIN,DEDLAM,.FALSE.,EDUM,LDOPIN,.FALSE.,
     $        .FALSE.,EDUM,EDUM)
      ENDIF
      IF (LDOVDL) THEN
         PRINT *
         PRINT *,'FREEENERGYCOMPS (dE/dLamda)*DLamda'
         CALL PRENER(LDOKIN,DVDL,.FALSE.,EDUM,LDOPIN,.FALSE.,
     $        .FALSE.,EDUM,EDUM)
      ENDIF
      IF (LDODMU) THEN
         PRINT *
         PRINT *,'FREEENERGYDERIVS (dE/dMu)'
         CALL PRENER(LDOKIN,DEDMU,.FALSE.,EDUM,LDOPIN,.FALSE.,
     $        .FALSE.,EDUM,EDUM)
      ENDIF
      IF (LDOVDM) THEN
         PRINT *
         PRINT *,'FREEENERGYCOMPS (dE/dMu)*DMu'
         CALL PRENER(LDOKIN,DVDM,.FALSE.,EDUM,LDOPIN,.FALSE.,
     $        .FALSE.,EDUM,EDUM)
      ENDIF
#endif
C end prarr
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE LSQFIT (ISTEP,FIT,FITI,FIT2)
C     ----------------------------------------------------------------==
C     W.F. VAN GUNSTEREN, GRONINGEN, OCT. 1984
C
C     SUBROUTINE LSQFIT (ISTEP,FIT,FITI,FIT2)
C
C     LSQFIT PERFORMS A LEAST SQUARES FIT OF A FUNCTION Y(X) TO
C     A STRAIGHT LINE AND PRINTS THE RESULTS.
C
C     ISTEP = NUMBER OF EQUIDISTANT X-VALUES X1, X2, ..., XISTEP, THAT
C     ARE USED IN THE FIT (.GT.1)
C     FIT = SUM OVER I=1,ISTEP OF Y(XI)
C     FITI = SUM OVER I=1,ISTEP OF Y(XI)*I
C     FIT2 = SUM OVER I=1,ISTEP OF Y(XI)**2
C     ----------------------------------------------------------------==
C args
      INTEGER ISTEP
      real*8 FIT,FITI,FIT2
C local vars
      INTEGER IM1,IP1,IPM1
      real*8 RLINE,ALINE,CHILIN,BLINE
C begin
      FIT=FIT/ISTEP
      FITI=FITI/ISTEP
      FIT2=FIT2/ISTEP
      IF (ISTEP.LT.2 .OR. DABS(FIT2-FIT*FIT) .LT. 1.D-8) RETURN
      IM1=ISTEP-1
      IP1=ISTEP+1
      IPM1=IP1*IM1
      ALINE=(2.0D0*FITI-IP1*FIT)*6.0D0/IPM1
      BLINE=FIT-ALINE*IP1/2.0D0
      RLINE=IPM1/(FIT2-FIT*FIT)/12.0D0
Comment
C     some machines (e.g. SGI when using hardware sqrt)
C     dont like taking the SQRT of zero, so
C     we use .LE. and not .LT. below
C end comment
      IF (RLINE .LE. 0.0D0)THEN
         RLINE = 0.0D0
      ELSE
         RLINE=ALINE* DSQRT(RLINE)
      ENDIF
      CHILIN=FIT2-2.0D0*ALINE*FITI+BLINE*(BLINE-2.0D0*FIT)+ALINE*
     .    (ALINE*DBLE(IP1)*DBLE(IP1+IM1+1)/6.0D0+BLINE*DBLE(IP1))
cmb      CHILIN=FIT2-2.0D0*ALINE*FITI-2.0D0*BLINE*FIT+BLINE*BLINE+
cmb     $  ALINE*ALINE*IP1*DBLE(IP1+IM1+1)/6.0D0+ALINE*BLINE*IP1

      IF (CHILIN .LE. 0.0D0) THEN
         CHILIN = 0.0D0
      ELSE
         CHILIN= DSQRT(CHILIN)
      ENDIF

      PRINT *,'FITBLOCK'
      PRINT 500,ISTEP
      PRINT 510,ALINE
      PRINT 520,BLINE
      PRINT 530,RLINE
      PRINT 540,CHILIN
      PRINT *,'END FITBLOCK'
      PRINT *

 500  FORMAT (4X,'FIT OVER',I5,' STEPS')
 510  FORMAT (4X,'DRIFT PER STEP (DY/DX)            :',E16.8)
 520  FORMAT (4X,'FUNCTION VALUE AT ZERO (Y AT X=0) :',E16.8)
 530  FORMAT (4X,'CORRELATION COEFFICIENT R         :',E16.8)
 540  FORMAT (4X,'R.M.S. FLUCTUATION CHI            :',E16.8)
C end lsqfit
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE CLCDOF(NRP,NRAM,NCONS,NCONG,NSPM,NRCON,JRC,V,
     $     WINV,WMAS,NUM4DA,
     $     NBON,IB,JB,
     $     NBONH,IBH,JBH,
     $     NSP,
     $     NCONP,NCONBC,NDOF)
C     ----------------------------------------------------------------==
C     CLCDOF calculates the degrees of freedom for solute and
C     solvent in three and four dimensions. In addition, in the case
C     of position constraints, the WINV and WMAS arrays and the 
C     velocity array V for those atoms that are position constrained
C     are modified set to zero.
C
C     NRP     number of solute atoms per molecule
C     NRAM    number of solvent atoms per molecule
C     NCONS   number of distance constraints in a solvent molecule
C     NCONG   the number of distance constraints per solute molecule
C     NSPM    the number of submolecules per solute molecule
C     NRCON   the number of positionally constrained atoms per solute
C             molecule (Note that there can only be one solute molecule
C             when using position constraining. This, however, is not
C             checked for in this routine.)
C     JRC(NRCON)
C             contains the atom sequence numbers of positionally
C             constrained atoms.
C     V(NRP*NDIM)
C             contains the atomic velocities. Modified if there are
C             positionally constrained atoms.
C             Note that NPM is always 1 if we have position
C             constraining.
C     WINV(NRP) inverse mass of each solute atom. Modified in the case of
C             position constraints
C     WMAS(NRP) mass of each solute atom. Modified in the case of
C             position constraints.
C     NUM4DA  the number of atoms that are in 4D.
C
C     NBON    the number of bonds per solute molecule between atoms
C             that are not hydrogens.
C     IB(NBON),JB(NBON)
C             atom numbers forming a bond in the solute.
C     NBONH    the number of bonds per solute molecule between atoms
C             that are not hydrogens.
C     IBH(NBONH),JB(NBONH)
C             atom numbers forming a bond in the solute.
C
C     NSP(NSPM)
C             NSP(I) contains the number of the last atom in
C             submolecule I
C     The following values are delivered:
C     NCONP   delivered with the degrees of freedom to subtract from
C             a solute molecule due to shaken bonds between atoms
C             that are not both positionally constrained.
C     NCONBC  delivered with the number of shaken bonds between atoms
C             that are both positionally constrained.
C     NDOF    delivered with the calculated degrees of freedom.
C             The values for indexing NDOF are in runmd.h
C     ----------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'runmd.h'
      INCLUDE 'md.h'
      INCLUDE 'box.h'
C args
      INTEGER NRP,NRAM,NBON,NBONH,NCONS
      INTEGER NCONP,NCONG,NSPM
      INTEGER NUM4DA
      INTEGER NRCON,JRC(NRCON)
      INTEGER IBH(NBONH),JBH(NBONH)
      INTEGER IB(NBON),JB(NBON)
      INTEGER NSP(NSPM)
      INTEGER NDOF(NFTMAX)
      real*8 V(NRP*NDIM)
      real*8 WINV(NRP),WMAS(NRP)
C local vars
      INTEGER I,J,K,M,NCONBC
      INTEGER IFIRST,ILAST,NN,NDIR,NDCOM
      INTEGER IPCNUM(MAXNSP)
C begin
      CALL mm_IAZZERO(NDOF,NFTMAX)

      NCONBC = 0

      CALL mm_IAZZERO(IPCNUM,NSPM)

C do we have position constraints ?
      IF (NTR .EQ. NTRCON) THEN
C calculate the number of positionally constrained atoms per submolecule
         IFIRST = 1
         DO NN=1,NSPM
            ILAST = NSP(NN)
            DO I=1,NRCON
               K = JRC(I)
               IF (IFIRST .LE. K .AND. K .LE. ILAST) THEN
                  IPCNUM(NN) = IPCNUM(NN) + 1
               ENDIF
            ENDDO
            IFIRST = ILAST + 1
         ENDDO

C     NCONBC = Number of constrained bonds for which both atoms are 
C     also positionally constrained.
C     To prevent constrained atoms from responding to forces or having
C     a momentum: 
C     set mass of constrained atoms to zero AND
C     set inverse mass of constrained atoms to zero.

C     Set initial velocities of constrained atoms to zero
C     Centre of mass motion for other atoms is not conserved.
C     As at least one atom is fixed, use this to suppress centre of mass 
C     motion.

         DO I = 1, NRCON
            WINV(JRC(I)) = 0.0D0
            WMAS(JRC(I)) = 0.0D0
            DO M = 1, NDIM
               V(((JRC(I)-1)*NDIM)+M) = 0.0D0
            ENDDO
         ENDDO


C     Correct number of constraints for the internal/rotational
C     degrees of freedom.
C     Add back 1 degree of freedom for each bond for which
C     both atoms are positionally constrained in 3D.

         IF (NTC .NE. NTCNON) THEN
            DO 3000 I = 1, NBONH
               DO 3010 J = 1,NRCON
                  IF (JRC(J) .EQ. IBH(I)) THEN
                     DO 3020 K = 1,NRCON
                        IF (JRC(K) .EQ. JBH(I)) THEN
                           NCONBC = NCONBC + 1
                        ENDIF
 3020                CONTINUE
                  ENDIF
 3010          CONTINUE
 3000       CONTINUE
            
            IF (NTC .EQ. NTCDOB) THEN
               DO 3200 I = 1, NBON
                  DO 3210 J = 1,NRCON
                     IF (JRC(J) .EQ. IB(I)) THEN
                        DO 3220 K = 1,NRCON
                           IF (JRC(K) .EQ. JB(I)) THEN
                              NCONBC = NCONBC + 1
                           ENDIF
 3220                   CONTINUE
                     ENDIF
 3210             CONTINUE
 3200          CONTINUE
            ENDIF
         ENDIF
      ENDIF

      NCONP = NCONG - NCONBC
CC
C     now calculate the degrees of freedom to subtract from
C     a) internal and rotational degrees of freedom (NDIR)
C     b) c.o.m. translational (NDCOM)
C
C     A 4D atom may not be positionally constrained (this
C     is checked in subroutine L<GTCOOR> in promd.f), so we never
C     have to subtract any degrees of freedom from the 4th dimension
C     due to position constraining.
CC
      NDIR = 0
      NDCOM= 0
      DO NN=1,NSPM
         IF (IPCNUM(NN) .GE. 1) THEN
C for the first positionally constrained atom in a submolecule,
C subtract degrees of freedom from the translational c.o.m. d.o.f.
C for any further pos. cons. atoms, subtract from the
C internal/ rotational d.o.f.
            NDCOM = NDCOM + 1
            NDIR  = NDIR + IPCNUM(NN) - 1
         ENDIF
      ENDDO

      NDOF(NFTPIR) = NPM*(NDRMAX*(NRP-NSPM-NDIR)-NCONP)
      NDOF(NFTPCM) = NPM*NDRMAX*(NSPM - NDCOM)
      NDOF(NFTSLV) = NSM*(NDRMAX*NRAM-NCONS)
      NDOF(NFTTO)  = NDOF(NFTPIR) + NDOF(NFTPCM)+ NDOF(NFTSLV)- NDFMIN
      NDOF(NFT4D)  = NUM4DA - NDFMI4
C end clcdof
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE EMGSTP(IFNCFG)
C     ----------------------------------------------------------------==
C     EMGSTP is called in the event of a shake error during a
C     simulation.
C     The final configuration is written to file IFNCFG
C     which is typically the final configuration file in
C     L<PROMD> after which STOP is envoked.
C     A line is added to the title so that a configuration
C     that was produced as the result of a shake error can
C     be identifed.
C     EMGSTP can be called by any routine that calls L<SHAKE>,
C     e.g. L<MDLEAP>, L<SDLEAP>, L<INICRD>.
C     Typically, EMGSTP is called if the value for NITP or NITS
C     returned by L<SHAKE> is zero.
C     ----------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'md.h'
      INCLUDE 'box.h'
C args
      INTEGER IFNCFG
C local vars
      LOGICAL LFORMF
      INTEGER NATTOT,ILN
      real*8 XDUM
C begin
      LFORMF = (NTXO .EQ. NTXOFO)
      NATTOT = NPM*NRP + NSM*NRAM

      IF (MDTLNS .LT. MDTLIN) THEN
C     append a message to the title
         MDTLNS = MDTLNS + 1
         ILN = MDTLNS
      ELSE
C     overwrite the last line
         ILN = MDTLIN
      ENDIF

      MDTITL(ILN) = '*******FINAL CONFIGURATION ON SHAKE ERROR*******'

      CALL WRCFG(IFNCFG,LFORMF,
     $     MDTLNS,MDTITL,
     $     NATTOT,NPM,NSM,NDIM,NDR*NPM,NDJV*NPM,
     $     .TRUE.,X,
     $     .TRUE.,V,
     $     NTB .NE. NTBVAC,BOX,
     $     NTP .NE. NTPOFF .AND. NTR .NE. NTROFF,XC,
     $     IDOPRO .EQ. IDOSD, SX,IG,
     $     .FALSE.,0,XDUM,XDUM,XDUM,XDUM,
     $     NTDR .LT. 0,RIIAVE,
     $     NTJR .NE. NTJROF,COSQAV,COSIAV,
     $     NTLE .NE. NTLEOF,NDLE,NLECFG,NLEMEM,ILEMEM,NLEVST)

      CALL CLSFIL(IFNCFG)

      CALL STOPGM('MM_SETUP_DR','FATAL ERROR')
C end EMGSTP
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE INICRD(IFNCFG,NATTOT,NRPT,NDIM,PRGSTR,X,V,F,
     $     C4D,NCONG,ICOG,JCOG,CONP,FCON,NRCON,JRC)
C     ----------------------------------------------------------------==
C     INIRCD initialises the startup coordinates and velocities
C     by calling L<SHAKE> if the input switches L<INIT> and L<NTC>
C     requires it.
C     The shaken coordinates and velocities are returned in the
C     Arrays X anf V. F is used as an auxiliary array.
C
C     INICRD is called by subroutine L<RUNMD>.
C     ----------------------------------------------------------------==
C     includes
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'md.h'
      INCLUDE 'formats.h'
C args
      INTEGER NATTOT,NDIM,NRPT,NCONG,IFNCFG
      INTEGER ICOG(NCONG),JCOG(NCONG)
      INTEGER NRCON,JRC(NRCON)
      real*8 CONP(NCONG),FCON(NCONG),C4D(NATTOT)
      real*8 X(NDIM*NATTOT),V(NDIM*NATTOT),F(NDIM*NATTOT)
      CHARACTER PRGSTR*(*)
C local vars
      INTEGER NATTO3,I3,NITP,NITS
      real*8 FCONS(MXCONS)
C begin
      NATTO3 = NDIM*NATTOT
      NITP = 0
      NITS = 0

C solute shaking is controlled by NTC and INIT
C solvent shaking is only controlled by INIT

      IF (INIT .EQ. INSHVX) THEN
C     shake coordinates
!$OMP parallel do private(I3)
         DO I3 = 1,NATTO3
            F(I3) = X(I3)
         ENDDO

         IF (NTC .NE. NTCNON .AND. 
     $        NPM .GT. 0 .AND. NRP .GT. 1) THEN
            CALL mm_SHAKE(NRPT,NRP,NPM,NDO4D(N4DBON),WINV,0,
     $           NCONG,ICOG,JCOG,CONP,FCON,
     $           SHKTOL,F,X,C4D,
     $           NRCON,JRC,.FALSE.,NITP,.FALSE.)
            PRINT 400,NITP
            IF (NITP .EQ. 0) THEN
               PRINT FMSHKP,PRGSTR
               CALL EMGSTP(IFNCFG)
            ENDIF
         ENDIF

         IF (NSM .GT. 0 .AND. NRAM .GT. 1) THEN
            CALL mm_SHAKE(NATTOT,NRAM,NSM,NDO4D(N4DBON),WINVS,
     $           NRPT,NCONS,ICONS,JCONS,CONS,FCONS,
     $           SHKTOL,F,X,C4D,
     $           0,JRC,.TRUE.,NITS,.TRUE.)
            PRINT 405,NITS
            IF (NITS .EQ. 0) THEN
               PRINT FMSHKS,PRGSTR
               CALL EMGSTP(IFNCFG)
            ENDIF
         ENDIF
      ENDIF

      IF (INIT .EQ. INSHV .OR. INIT .EQ. INSHVX) THEN
C     shake velocities
C     F: calc the positions a time step before
!$OMP parallel do private(I3)
         DO I3 = 1,NATTO3
            F(I3) = X(I3) - DT*V(I3)
         ENDDO

         IF (NTC .NE. NTCNON .AND. 
     $        NPM .GT. 0 .AND. NRP .GT. 1) THEN
            CALL mm_SHAKE(NRPT,NRP,NPM,NDO4D(N4DBON),WINV,0,
     $           NCONG,ICOG,JCOG,CONP,FCON,
     $           SHKTOL,X,F,C4D,
     $           NRCON,JRC,.FALSE.,NITP,.FALSE.)
            PRINT 410,NITP
            IF (NITP .EQ. 0) THEN
               PRINT FMSHKP,PRGSTR
               CALL EMGSTP(IFNCFG)
            ENDIF
         ENDIF

         IF (NSM .GT. 0 .AND. NRAM .GT. 1) THEN
            CALL mm_SHAKE(NATTOT,NRAM,NSM,NDO4D(N4DBON),WINVS,
     $           NRPT,NCONS,ICONS,JCONS,CONS,FCONS,
     $           SHKTOL,X,F,C4D,
     $           0,JRC,.TRUE.,NITS,.TRUE.)
            PRINT 415,NITS
            IF (NITS .EQ. 0) THEN
               PRINT FMSHKS,PRGSTR
               CALL EMGSTP(IFNCFG)
            ENDIF
         ENDIF

!$OMP parallel do private(I3)
         DO I3=1,NATTO3
            V(I3) = (X(I3)-F(I3))/DT
         ENDDO
      ENDIF
      PRINT *

 400  FORMAT (1X,'mm_SHAKE Protein Coordinates, NITP =',I5)
 405  FORMAT (1X,'mm_SHAKE Solvent Coordinates, NITS =',I5)
 410  FORMAT (1X,'mm_SHAKE Protein Velocities,  NITP =',I5)
 415  FORMAT (1X,'mm_SHAKE Solvent Velocities,  NITS =',I5)
C end inicrd
C     ==--------------------------------------------------------------==
      RETURN
      END

#include "inc_runmd2.h"

