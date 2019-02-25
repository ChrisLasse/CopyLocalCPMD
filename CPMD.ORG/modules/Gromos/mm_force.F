C     ==================================================================
      SUBROUTINE mm_force(NAX,NSX,tau,mm_FION,epot_mm)
      IMPLICIT NONE
C     ==--------------------------------------------------------------==
c     == This routine is the second part of runmd.                    ==
C     == It is used to calculate MM forces                            ==
C     ==--------------------------------------------------------------==
cmb - Revised on 22 June 2006
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
C ES
      INCLUDE '../ropt.inc' ! INFI
c qmmm
      integer NAX,NSX
      real*8 epot_mm
      real*8 tau(3,NAX,NSX)
      real*8 mm_FION(3,NAX,NSX)
C Qmmm: the arguments of mm_setup_dr are passed through this commons
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
      real*8 FION_mm(3,NAX,NSX)
      INTEGER NRPT3,NSOLV3,NATTO3,I3,K,N,NITP,NITS
      INTEGER NCONP,NCONBC
#ifdef __VECTOR
      INTEGER NSTEP,II,ia,is
#else
      INTEGER NSTEP,II,i
#endif
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
      real*8 DTMP,DTMP2,DTMP3,DTMP4

C     NOTE: for all those arrays for which averages and r.m.s fluctuations
C     are collected and printed at the end of the run, five arrays
C     are kept. For example, for the temperature we have:
C
C     TEMP: contains the temperatures calculated at every step.
C     
C     TMPSUM, TMPSQ: contain the sum and the sum of the squares
C     of all TEMP values over the run. These variables are always
C     in real*8 for better statistics.
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
      real*8 TMPSUM(NFTMAX),TMPSQ(NFTMAX)

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
      real*8 EPLJSU(MXNRE2),EPLJSQ(MXNRE2)
C
      real*8 EPEL(MXNRE2),EPELAV(MXNRE2),EPELFL(MXNRE2)
      real*8 EPELSU(MXNRE2),EPELSQ(MXNRE2)
C
      real*8 EPRF(MXNRE2),EPRFAV(MXNRE2),EPRFFL(MXNRE2)
      real*8 EPRFSU(MXNRE2),EPRFSQ(MXNRE2)
C
      real*8 EPRC(MXNRE2),EPRCAV(MXNRE2),EPRCFL(MXNRE2)
      real*8 EPRCSU(MXNRE2),EPRCSQ(MXNRE2)

C energy , dE/dLamda and dE/dmu
      real*8 ENER(MXETBL),EPAVE(MXETBL),EPFLC(MXETBL)
      real*8 EPSUM(MXETBL),EPSQ(MXETBL)
C
      real*8 DEDLAM(MXETBL),DEDLAV(MXETBL),DEDLFL(MXETBL)
      real*8 DEDLSU(MXETBL),DEDLSQ(MXETBL)
C
      real*8 DEDMU(MXETBL),DEDMAV(MXETBL),DEDMFL(MXETBL)
      real*8 DEDMSU(MXETBL),DEDMSQ(MXETBL)


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
      real*8 ERESUM(MXCTBL),ERESQ(MXCTBL)

C VOLPRT: some sytem parameters saved to file.
C NOTE that the values in this array do NOT necessarily reflect
C the current status of the system at all times;
C it is only updated when needed for printing or writing to file.
C Hence, e.g., when using the pressure of the system in a calculation
C within the code, the variable PRESTO should be used, and not the
C corresponding element in VOLPRT.
      real*8 VOLPRT(MXVTBL),VPAVE(MXVTBL),VPFLC(MXVTBL)
      real*8 VPSUM(MXVTBL),VPSQ(MXVTBL)
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

#include "inc_runmd5.h"

      INCLUDE 'mm_save_config.h'
c
C data statements
      DATA PRGSTR /'RUNMD'/
      DATA LSEFST /.TRUE./
      DATA NSTEP /0/
C begin
      SAVE prgstr
      SAVE lsefst
      SAVE nstep

c qmmm: positions in cp format are converted in gromos format
      CALL mm_gcp(NAX,NSX,ener(iptot),epot_mm,'G','P',x,tau)

cmb    write (60,*) 'IN THE MM_FORCE ROUTINE'
cmb    do i=1,NRPT+NSOLV
cmb      write(60,*)i,'pos',(x(3*(i-1)+im),im=1,3),
cmb  $                    'vel',(v(3*(i-1)+im),im=1,3)
cmb    enddo

C     ------------- MAIN MD LOOP --------------------------
#ifdef READ_COORD
c     DO 1000 NSTEP=0,NSTLIM-1,NSTINC
#else
c     DO 1000 NSTEP=0,NSTLIM-1
#endif
C ES: take the timestep information from CPMD -> synchronous
      NSTEP=INFI
C     set the constraint lengths and masses according to RLAM
         IF (LDOPER) THEN
            CALL SETPRT(NRP,NPM,NSPM,RLAM,DLAM,
     $           WMAS,WINV,TOTMAS,NSPG,SUBMAS,SUBMIN,
     $           IDOPRO .EQ. IDOSD,TEMPSD*BOLTZ,SWINK,
     $           NTC .NE. NTCNON, NTC .EQ. NTCDOB,NCONG,CONP)
         ENDIF

C place particles back into box if necessary
         IF (.NOT. LVAC) THEN
            CALL mm_SHIAG(NATTOT,NPM,NSM,XMIN,X,LEVERY)
         ENDIF

C prepare virial calculation if necessary
         IF (LDOVIR) THEN
#ifdef EWALD
            CALL PRPVIR2(NATTOT,NPM,NSM,X,V,XR,XCEN,TMASS,
     $           NSPM,NSP,SUBMAS,
     $           EKCM,EKCMTO,LEVERY)
#else
            CALL PRPVIR(NATTOT,NPM,NSM,X,V,XR,TMASS,NSPM,NSP,SUBMAS,
     $           EKCM,EKCMTO,LEVERY)
#endif
         ENDIF

#include "inc_runmd3.h"

         CALL FORCE(NATTOT,X,F,XR,
     $        XCEN,
     $        C4D,CW4D,
     $        ENER,ENERES,DEDLAM,DEDMU,EPEL,EPRF,EPRC,EPLJ,
     $        VIR,
     $        NSPM,NSP,
     $        NDR,
     $        IDR1,JDR1,KDR1,LDR1,ICDR1,
     $        IDR2,JDR2,KDR2,LDR2,ICDR2,
     $        R0,W0,DISH,DISC,RIIAVE,
     $        NRRST,JRC,CXC,XC,
     $        NDJV,IPJV, JPJV, KPJV, LPJV,
     $        CPJV,PJR0,PSJR,AJV,BJV,CJV,
     $        COSQAV,COSIAV,
     $        NDLR,IPLR,JPLR,KPLR,LPLR,ICPLR,
     $        CPLR,PDLR,
     $        NDLE,IPLE,JPLE,KPLE,LPLE,
     $        NLECFG,NLEMEM,ILEMEM,NLEVST,
     $        IAGRP,LEVERY,LERR,.FALSE.)

c
c qmmm: gromos forces are converted into cp format and summed to the
c       forces of the quantum part
c
         CALL azzero(FION_mm,NAX*NSX*3)
         CALL mm_gcp(NAX,NSX,ener(iptot),epot_mm,'C','F',f,FION_mm)

         IF (LERR) THEN
            PRINT *,PRGSTR,': Error in subroutine FORCE'
            CALL STOPGM('MM_FORCE','FATAL ERROR')
         ENDIF

         CALL DAXPY(3*NSX*NAX,1.0d0,FION_mm(1,1,1),1,mm_FION(1,1,1),1)

#include "inc_runmd4.h"

C set this here to control pairlist construction in FORCE
C for next step
         LTNB = (MOD(NSTEP+1,NSNB) .EQ. 0)

C        QMMM QMMM QMMM QMMM time to return, we do NOT want gromos to do anything
C        more QMMM QMMM QMMM

C        ...hmmm. we could have gromos compute some energies
#if 0 
       GOTO 8000
       RETURN

cmb    write (6,*) 'IN THE MM_FORCE ROUTINE'
cmb    do i=1,NRPT+NSOLV
cmb      write(6,*)i,'pos',(x(3*(i-1)+im),im=1,3),
cmb   $        'vel',(v(3*(i-1)+im),im=1,3),
cmb   $        'force',(f(3*(i-1)+im),im=1,3)
cmb    enddo

         IF (LDOVIR) THEN
            CALL CLCVIR(EKCM,VIR,PRES,EKCMTO,VIRTOT,PRESTO)
         ENDIF

C write pos and box to trajectory
         IF (NTWX .NE. 0) THEN
            IF (NTWSE .EQ. 0) THEN
C     normal trajectory
               IF (MOD(NSTEP,NTWX) .EQ. 0) THEN
                  CALL WRTIME(IUTRJX,LFORM,NSTEP,TIME)
                  IF (NTWX .GT. 0) THEN
                     CALL WRXRED(IUTRJX,LFORM,LWR4,NATTOT,NDIM,X)
                  ELSE
C     only write solute coords
                     CALL WRXRED(IUTRJX,LFORM,LWR4,NRPT,NDIM,X)
                  ENDIF
                  IF (NTB .NE. NTBVAC) THEN
                     CALL WRBOX(IUTRJX,LFORM,BOX)
                  ENDIF
               ENDIF
            ELSE
C minimum energy trajectory
               LSEWRT  = .NOT. LSEFST .AND. (MOD(NSTEP,NTWX) .EQ. 0)
               LSECPY  = LSEFST .OR. (MOD(NSTEP,NTWX) .EQ. 0)
               LSEFST  = .FALSE.

               IF (LSEWRT) THEN
C     write what we have saved in a preceding time frame
                  CALL WRSAVE(IUTRJX,LFORM,
     $                 NTWX .GT. 0,NTB .NE. NTBVAC,LWR4,
     $                 LPRKIN,LDOVIR,LDOPIN,
     $                 NUSNRE,NUNRE2,NATTOT,NRPT,NDIM,XSAVE)
               ENDIF

C compare current energy with saved energy and store configuration
C if necessary. The updating of the ESAVE array occurs AFTER mdleap
C when we have completely filled in the ENER array (also the kinetic terms)

               IF (LSECPY .OR. LTENER(ENER(NTWSE))) THEN
C we use LSECPY for copying the energy arrays lower down
                  LSECPY = .TRUE.
                  IF (NTWX .GT. 0) THEN
                     NLIM = NATTO3
                  ELSE
                     NLIM = NRPT*NDIM
                  ENDIF
                  CALL SAVTIM(NLIM,X,NDIM,BOX,TIME,NSTEP,XSAVE)
               ENDIF
            ENDIF
         ENDIF

C now write velocities
         IF (NTWV .NE. 0) THEN
            IF (MOD(NSTEP,NTWV) .EQ. 0) THEN
               CALL WRTIME(IUTRJV,LFORM,NSTEP,TIME)
               IF (NTWV .LT. 0) THEN
C     only write solute
                  CALL WRVRED(IUTRJV,LFORM,LWR4,NRPT,NDIM,V)
               ELSE
                  CALL WRVRED(IUTRJV,LFORM,LWR4,NATTOT,NDIM,V)
               ENDIF
            ENDIF
         ENDIF


C scaling for Berendsen thermostat
         DO II=1,NBATH
            EKTMP = 0.0D0
            DO K=1,NBNUM(II)
              EKTMP = EKTMP + ESCL(NBNDX(II,K))
            ENDDO

            IF (EKTMP .GT. EPS) THEN
               SCLTMP=DSQRT(1.0D0+DTBATH(II)*(EKREF(II)/EKTMP-1.0D0))
               SCLBAT(II) = SCLTMP
               IF (NBNUM(II) .EQ. 1 .AND. NBNDX(II,1) .EQ. NFTTO) THEN
C special case: we couple to the WHOLE system in 3D: set all SCAL vals
                  SCAL(NFTPIR) = SCLTMP
                  SCAL(NFTPCM) = SCLTMP
                  SCAL(NFTSLV) = SCLTMP
C special case: 4D harmonic force constant scaling and simulataneous
C scaling of reference temperature in 4th dimension
               ELSEIF (NBNUM(II) .EQ. 1 .AND. 
     $                 NBNDX(II,1) .EQ. NFT4D) THEN
                  SCAL(NFT4D) = SCLTMP
                  IF (NTCW4D .EQ. NTC4YE) THEN
                     TMP = DBLE(NSTEP)/DBLE(NSTLIM)
                     CW4D = (CW4DB - CW4DA)*TMP + CW4DA
                     IF (NTT4 .NE. NTT4OF) THEN
                        EKREF(II) = (EKIN4B-EKIN40)*TMP+EKIN40
                     ENDIF
                  ENDIF
               ELSE
                  DO K=1,NBNUM(II)
                     SCAL(NBNDX(II,K)) = SCLTMP
                  ENDDO
               ENDIF
            ENDIF
         ENDDO

#ifdef READ_COORD
         IF (NRDCOR.EQ.0) THEN
#endif
C integration step
         IF (IDOPRO .EQ. IDOMD) THEN
            CALL MDLEAP(IFNCFG,NATTOT,NPM,NRP,NSM,NRAM,NDO4D(N4DBON),
     $           NSPM,NSP,SUBMAS,SUBMIN,
     $           WINV,WMAS,WINVS,WMASS,DT,
     $           X,V,F,C4D,VELCM,
     $           LSHKPR,NCONG,ICOG,JCOG,CONP,FCON,NRCON,JRC,SHKTOL,
     $           LSHKSO,NCONS,ICONS,JCONS,CONS,
     $           SCAL,EKNOW,ESCL,
     $           LEVERY,NITP,NITS)
         ELSEIF (IDOPRO .EQ. IDOSD) THEN
            CALL SDLEAP(IFNCFG,NATTOT,NPM,NRP,NSM,NRAM,NDO4D(N4DBON),
     $           NSPM,NSP,SUBMAS,SUBMIN,
     $           WINV,WMAS,WINVS,WMASS,DT,
     $           X,V,F,SX,C4D,VELCM,
     $           LSHKPR,NCONG,ICOG,JCOG,CONP,FCON,NRCON,JRC,SHKTOL,
     $           LSHKSO,NCONS,ICONS,JCONS,CONS,
     $           SCAL,EKNOW,ESCL,
     $           LEVERY,NITP,NITS,
     $           IG,CC1,CC2,CC3,CC4,CC5,CC6,CC7,CC8,CC9,
     $           SWINK,SWINKS)
         ELSE
            PRINT *,PRGSTR,': illegal value for IDOPRO = ',IDOPRO
            CALL STOPGM('MM_FORCE','FATAL ERROR')
         ENDIF

         ENER(IKSUCM) = EKNOW(NFTPCM)
C path integral energy is subtracted from the 3d internal kinetic
C energy, hence only from the total kinetic energy and not from the
C translational kinetic energy
         ENER(IKSLU3) = EKNOW(NFTPIR) + EKNOW(NFTPCM) - ENER(IPPISP)
         ENER(IKSLU4) = ENER(IKSLU3)  + EKNOW(NFTP4D) 
         ENER(IKSLV3) = EKNOW(NFTSLV)
         ENER(IKSLV4) = EKNOW(NFTSLV) + EKNOW(NFTS4D)

         ENER(IKTOT4) = ENER(IKSLU4)+ENER(IKSLV4)
         ENER(IKTOT3) = ENER(IKSLU3)+ENER(IKSLV3)
         ENER(IK4THD) = ENER(IKTOT4)-ENER(IKTOT3)

#include "inc_runmd6.h"

         ENER(IETOT) = ENER(IKTOT4)+ENER(IPTOT)

         IF (LDOPER) THEN
C     calculate kinetic and distance constraint contributions
            CALL CLCPRT(NATTOT,NRP,NPM,NSPM,NSP,
     $           DT,RLAM,DEDLAM,
     $           V,VCM,
     $           NTC .NE. NTCNON, NTC .EQ. NTCDOB,NCONG,FCON)

C path integral correction to kinetic contribution 
            DEDLAM(IKSLU3) = DEDLAM(IKSLU3)-DEDLAM(IPPISP)
            DEDLAM(IKSLU4) = DEDLAM(IKSLU4)-DEDLAM(IPPISP)
            DEDLAM(IKTOT4) = DEDLAM(IKSLU4)+DEDLAM(IKSLV4)
            DEDLAM(IKTOT3) = DEDLAM(IKSLU3)+DEDLAM(IKSLV3)
            DEDLAM(IK4THD) = DEDLAM(IKTOT4)-DEDLAM(IKTOT3)

            DEDLAM(IETOT)  = DEDLAM(IKTOT4)+DEDLAM(IPTOT)

            DEDMU(IKTOT4) = DEDMU(IKSLU4)+DEDMU(IKSLV4)
            DEDMU(IKTOT3) = DEDMU(IKSLU3)+DEDMU(IKSLV3)
            DEDMU(IK4THD) = DEDMU(IKTOT4)-DEDMU(IKTOT3)

            DEDMU(IETOT)  = DEDMU(IKTOT4)+DEDMU(IPTOT)

C calculate DLAM* (dV/dLamda)
C and DMU*(dV/dMu)
#if defined(__VECTOR)
!$OMP parallel do private(II) shared(DLAM,DMU)
#endif
            DO II=1,MXETBL
               DVDL(II) = DEDLAM(II)*DLAM
               DVDM(II) = DEDMU(II)*DMU
            ENDDO
         ENDIF

C calc temperatures
         DO II=1,NFTMAX
            IF (TFACPR(II) .GE. EPS) THEN
               TEMP(II) = EKNOW(II)/TFACPR(II)
            ELSE
               TEMP(II) = 0.0D0
            ENDIF
         ENDDO

C bring VOLPRT up to date for writing to file/printing to screen
 8000    CONTINUE
#endif
         VOLPRT(IVSPIR) = SCAL(NFTPIR)
         VOLPRT(IVSPCM) = SCAL(NFTPCM)
         VOLPRT(IVSCLS) = SCAL(NFTSLV)
         VOLPRT(IVSCL4) = SCAL(NFT4D)
         VOLPRT(IVBOXX) = BOX(1)
         VOLPRT(IVBOXY) = BOX(2)
         VOLPRT(IVBOXZ) = BOX(3)
         VOLPRT(IVVOL)  = BOX(1)*BOX(2)*BOX(3)
         IF (LOCTO) THEN
            VOLPRT(IVVOL) = VOLPRT(IVVOL)*0.5D0
         ENDIF

         VOLPRT(IVPRSX) = PRES(1)
         VOLPRT(IVPRSY) = PRES(2)
         VOLPRT(IVPRSZ) = PRES(3)
         VOLPRT(IVPRES) = PRESTO

         VOLPRT(IVVIRX) = VIR(1)
         VOLPRT(IVVIRY) = VIR(2)
         VOLPRT(IVVIRZ) = VIR(3)
         VOLPRT(IVVIR)  = VIRTOT

         VOLPRT(IVCMX)  = EKCM(1)
         VOLPRT(IVCMY)  = EKCM(2)
         VOLPRT(IVCMZ)  = EKCM(3)
         VOLPRT(IVEKCM) = EKCMTO

         VOLPRT(IVRLAM) = RLAM
         VOLPRT(IVRMU)  = RMU

C save the current energy arrays (for low energy trajectories)
!         IF (LSECPY) THEN
!            CALL SAVNRG(ENER,ENERES,VOLPRT,TEMP,
!     $           EPLJ,EPEL,EPRF,EPRC)
!         ENDIF

C writing of energies and volume,pressure etc. to energy trajectory
!
!         IF (NTWE .NE. 0) THEN
!            IF (MOD(NSTEP,NTWE) .EQ. 0) THEN
!               CALL WRTIME(IUTRJE,LFORM,NSTEP,TIME)
!               CALL WRNRG(IUTRJE,LFORM,
!     $              MXEWRT,ENER,
!     $              MXCTBL,ENERES,
!     $              NUSNRE,EPLJ,EPEL,EPRF,EPRC)
!               CALL WRVPRT(IUTRJE,LFORM,MXVWRT,VOLPRT)
!            ENDIF
!         ENDIF
!
C free energies
!         IF (NTG .NE. NTGOFF .AND. NTWG .NE. 0) THEN
!            IF (MOD(NSTEP,NTWG) .EQ. 0) THEN
!               CALL WRTIME(IUTRJG,LFORM,NSTEP,TIME)
!
!               IF (NTG .EQ. NTGLAM .OR. NTG .EQ. NTGBOT) THEN
!                  CALL WRRLAM(IUTRJG,LFORM,MXNRGF,ENER,
!     $                 RLAM,MXEWRT,DEDLAM)
!               ENDIF
!               IF (NTG .EQ. NTGMU .OR. NTG .EQ. NTGBOT) THEN
!                  CALL WRRMU(IUTRJG,LFORM,MXNRGF,ENER,
!     $                 RMU,MXEWRT,DEDMU)
!               ENDIF
!            ENDIF
!         ENDIF
!
C print values to screen 
         IF (NTPR .NE. 0) THEN
            IF (MOD(NSTEP,NTPR) .EQ. 0) THEN
               PRINT *,'STEPBLOCK'
               PRINT 500,'STEP','NITP','NITS','TIME'
               PRINT 501,NSTEP,NITP,NITS,TIME
               CALL PRARR(LPRKIN,LDOVIR,LDOPIN,
     $              .TRUE.,ENER,ENERES,VOLPRT,TEMP,
     $              LDODLA,DEDLAM,
     $              LDOVDL,DVDL,
     $              LDODMU,DEDMU,
     $              LDOVDM,DVDM,
     $              NUSNRE,NUNRE2,EPEL,EPRF,EPRC,EPLJ)
               PRINT *,'END STEPBLOCK'
               PRINT *
            ENDIF
         ENDIF
C (GM) hand the control back to the QMMM part
         RETURN

C use total energy for least squares fitting
         IF (NTPL .NE. 0) THEN
            IFSTEP = IFSTEP + 1
            FITAR(NFIT)  = FITAR(NFIT) + ENER(IETOT)
            FITAR(NFITI) = FITAR(NFITI)+ ENER(IETOT)*IFSTEP
            FITAR(NFIT2) = FITAR(NFIT2)+ ENER(IETOT)**2
         ENDIF

C add energies, volprt and temperatures to averages and average square
         DO II=1,MXETBL
            DTMP = ENER(II)
            EPSUM(II) = EPSUM(II) + DTMP
            EPSQ(II)  = EPSQ(II)  + DTMP*DTMP
         ENDDO

         DO II=1,MXCTBL
            DTMP = ENERES(II)
            ERESUM(II) = ERESUM(II) + DTMP
            ERESQ(II)  = ERESQ(II)  + DTMP*DTMP
         ENDDO

         DO II=1,MXVTBL
            DTMP = VOLPRT(II)
            VPSUM(II) = VPSUM(II) + DTMP
            VPSQ(II)  = VPSQ(II)  + DTMP*DTMP
         ENDDO

         DO II=1,MXNRE2
            DTMP = EPLJ(II)
            EPLJSU(II) = EPLJSU(II) + DTMP
            EPLJSQ(II) = EPLJSQ(II) + DTMP*DTMP

            DTMP2 = EPEL(II)
            EPELSU(II) = EPELSU(II) + DTMP2
            EPELSQ(II) = EPELSQ(II) + DTMP2*DTMP2

            DTMP3 = EPRF(II)
            EPRFSU(II) = EPRFSU(II) + DTMP3
            EPRFSQ(II) = EPRFSQ(II) + DTMP3*DTMP3

            DTMP4 = EPRC(II)
            EPRCSU(II) = EPRCSU(II) + DTMP4
            EPRCSQ(II) = EPRCSQ(II) + DTMP4*DTMP4
         ENDDO

         DO II=1,NFTMAX
            DTMP = TEMP(II)
            TMPSUM(II) = TMPSUM(II) + DTMP
            TMPSQ(II)  = TMPSQ(II)  + DTMP*DTMP
         ENDDO

         IF (LDOPER) THEN
            DO II=1,MXETBL
C dE/dLamda
               DTMP = DEDLAM(II)
               DEDLSU(II) = DEDLSU(II) + DTMP
               DEDLSQ(II) = DEDLSQ(II) + DTMP*DTMP

C (dE/dLamda)*DLAM
               DLSUM(II) = DLSUM(II) + DVDL(II)

C dE/dMu
               DTMP2 = DEDMU(II)
               DEDMSU(II) = DEDMSU(II) + DTMP2
               DEDMSQ(II) = DEDMSQ(II) + DTMP2*DTMP2

C (dE/dMu)*DMU
               DMSUM(II) = DMSUM(II) + DVDM(II)
            ENDDO
         ENDIF

C rescale coords if we have pressure coupling
         IF (NTP .NE. NTPOFF) THEN
            CALL SCLCRD(NATTOT,X,XC,PRES,PRESTO)
         ENDIF

C recalculate friction coeffs if necessary
         IF (IDOPRO .EQ. IDOSD .AND.
     $        NTFR .EQ. NTFRCA .AND. MOD(NSTEP+1,NSFR) .EQ. 0) THEN
!           CALL FRIC(NATTOT,NDIM,.TRUE.,NBREF,RCUTF,CFRIC,DT,
!    $           X,GAM,
!    $           CC1,CC2,CC3,CC4,CC5,CC6,CC7,CC8,CC9)
         ENDIF

C centre of mass printing (and removal if necessary)
         LREMCM = (NSCM .NE. 0)
         IF (LREMCM) THEN
            LREMCM = (MOD(NSTEP+1,NSCM) .EQ. 0)
         ENDIF

         LPRLSQ = (NTPL .NE. 0)
         IF (LPRLSQ) THEN
            LPRLSQ = (MOD(NSTEP+1,NTPL) .EQ. 0)
         ENDIF

         IF (LREMCM .OR. LPRLSQ) THEN
#if defined(__VECTOR)
!$OMP parallel do private(I3) shared(DTHALF)
#endif
#ifdef __SR11000
*poption parallel, tlocal(I3)
#endif
            DO I3 = 1,NATTO3
               F(I3) = X(I3) - V(I3)*DTHALF
            ENDDO
         ENDIF

         IF (LREMCM) THEN
!           CALL CENMAS(NATTOT,NPM,NRP,NSM,NRAM,0,NDIM,NDRMAX,F,V,
!    $           TOTMAS,0,WMAS,WMASS,
!    $           EKCMTO,XCM,VCM,ACM,EKROT,OCM,ICMROT)

!           CALL STOPCM(NATTOT,NDIM,F,V,XCM,VCM,OCM,ISCROT)
         ENDIF

         IF (LREMCM .OR. LPRLSQ) THEN
!           CALL CENMAS(NATTOT,NPM,NRP,NSM,NRAM,0,NDIM,NDRMAX,F,V,
!    $           TOTMAS,0,WMAS,WMASS,
!    $           EKCMTO,XCM,VCM,ACM,EKROT,OCM,ICMROT)
         ENDIF

C print out fit to total energy and reset every NTPL steps
         IF (LPRLSQ) THEN
            CALL LSQFIT(NTPL,FITAR(NFIT),FITAR(NFITI),FITAR(NFIT2))

            IFSTEP = 0
            CALL mm_AZZERO(FITAR,MAXFIT)
         ENDIF

C  end of main md loop
#ifdef READ_COORD
         TIME = TIME + DT * NSTINC
#else
         TIME = TIME + DT
#endif
         RLAM = RLAM + DLAM
         RMU  = RMU  + DMU

 1000 CONTINUE

C if we write a minimum energy trajectory, then we should write
C the saved configuration and energy of the last time frame here.
      IF (NTWSE .NE. 0) THEN
         CALL WRSAVE(IUTRJX,LFORM,
     $        NTWX .GT. 0,NTB .NE. NTBVAC,LWR4,
     $        LPRKIN,LDOVIR,LDOPIN,
     $        NUSNRE,NUNRE2,NATTOT,NRPT,NDIM,XSAVE)
      ENDIF

C calc averages and r.m.s. fluctuations and print them
      CALL CLCAVE(NSTLIM,MXETBL,EPSUM,EPSQ,EPAVE,EPFLC)
      CALL CLCAVE(NSTLIM,MXCTBL,ERESUM,ERESQ,EREAVE,EREFLC)
      CALL CLCAVE(NSTLIM,MXVTBL,VPSUM,VPSQ,VPAVE,VPFLC)
      CALL CLCAVE(NSTLIM,NFTMAX,TMPSUM,TMPSQ,TMPAVE,TMPFLC)

      CALL CLCAVE(NSTLIM,MXNRE2,EPLJSU,EPLJSQ,EPLJAV,EPLJFL)
      CALL CLCAVE(NSTLIM,MXNRE2,EPELSU,EPELSQ,EPELAV,EPELFL)
      CALL CLCAVE(NSTLIM,MXNRE2,EPRFSU,EPRFSQ,EPRFAV,EPRFFL)
      CALL CLCAVE(NSTLIM,MXNRE2,EPRCSU,EPRCSQ,EPRCAV,EPRCFL)

      CALL CLCAVE(NSTLIM,MXETBL,DEDLSU,DEDLSQ,DEDLAV,DEDLFL)
      CALL CLCAVE(NSTLIM,MXETBL,DEDMSU,DEDMSQ,DEDMAV,DEDMFL)

      PRINT 5010
      PRINT 5020
      CALL PRARR(LPRKIN,LDOVIR,LDOPIN,
     $     .TRUE.,EPAVE,EREAVE,VPAVE,TMPAVE,
     $     LDODLA,DEDLAV,
     $     .FALSE.,EDUM,
     $     LDODMU,DEDMAV,
     $     .FALSE.,EDUM,
     $     NUSNRE,NUNRE2,EPELAV,EPRFAV,EPRCAV,EPLJAV)

      PRINT 5030
      CALL PRARR(LPRKIN,LDOVIR,LDOPIN,
     $     .TRUE.,EPFLC,EREFLC,VPFLC,TMPFLC,
     $     LDODLA,DEDLFL,
     $     .FALSE.,EDUM,
     $     LDODMU,DEDMFL,
     $     .FALSE.,EDUM,
     $     NUSNRE,NUNRE2,EPELFL,EPRFFL,EPRCFL,EPLJFL)


C print out the sum of dE/dLamda*DLamda and corresponding for Mu
      IF (NTG .NE. NTGOFF) THEN
         PRINT 5040
         PRINT *,'FREEENERGYCOMPS (dE/dLamda)*DLamda'
         CALL PRENER(.TRUE.,DLSUM,.FALSE.,EDUM,LDOPIN,.FALSE.,
     $        .FALSE.,EDUM,EDUM)

         PRINT *,'FREEENERGYCOMPS (dE/dMu)*DMu'
         CALL PRENER(.TRUE.,DMSUM,.FALSE.,EDUM,LDOPIN,.FALSE.,
     $        .FALSE.,EDUM,EDUM)
      ENDIF

      PRINT *
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
C     ==--------------------------------------------------------------==
      RETURN
      END
