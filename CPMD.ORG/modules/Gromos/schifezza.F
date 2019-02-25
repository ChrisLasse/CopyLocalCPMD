C     multiple machines prefer blockdata here (actually all _should_ use it).
C     AK 2006/10/24
#if defined(__SGI) || defined(__Linux) || defined(__OSX) || defined(__IBM)
#define __USE_BLOCKDATA 1
#endif
C     ==================================================================
      SUBROUTINE schifezza
C     ==--------------------------------------------------------------==
cmb - Revised on 13 Oct. 2006
#if defined(__USE_BLOCKDATA)
      RETURN
      END
      BLOCKDATA SCHIFEZZADATA    
#endif       

      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
      INCLUDE 'box.h'
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'formats.h'
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
      INCLUDE 'pertsz.h'
      INCLUDE 'pertar.h'
      INCLUDE 'ptblock.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'units.h'

C     make sure this common block is changed
C     consistently within this file
C     --------------------------------------
      CHARACTER*(8)STR4D
      COMMON /PFOUR/STR4D(MX4DIN)
C     --------------------------------------
      DATA STR4D(N4DBON) /'BND/CONS'/
      DATA STR4D(N4DBAN) /'BND-ANGL'/
      DATA STR4D(N4DDIH) /'DIHEDRAL'/
      DATA STR4D(N4DNBD) /' NONBOND'/
      DATA STR4D(N4DRST) /'D-RESTR.'/
      DATA STR4D(N4DHAR) /' 4D HARM'/

      DATA BTNAME(NTITLE)   /       'TITLE'  /
      DATA LNEDMD(NTITLE)   /       .TRUE.   /
      DATA LNEDEM(NTITLE)   /       .TRUE.   /

      DATA BTNAME(NBOUND)   /    'BOUNDARY'  /
      DATA LNEDMD(NBOUND)   /       .TRUE.   /
      DATA LNEDEM(NBOUND)   /       .TRUE.   /

      DATA BTNAME(NTCOUP)   /     'TCOUPLE'  /
      DATA LNEDMD(NTCOUP)   /       .TRUE.   /
      DATA LNEDEM(NTCOUP)   /       .FALSE.  /

      DATA BTNAME(NPCOUP)   /     'PCOUPLE'  /
      DATA LNEDMD(NPCOUP)   /      .FALSE.   /
      DATA LNEDEM(NPCOUP)   /      .FALSE.   /

      DATA BTNAME(NNSP)     /'SUBMOLECULES'  /
      DATA LNEDMD(NNSP)     /       .TRUE.  /
      DATA LNEDEM(NNSP)     /       .FALSE.  /

      DATA BTNAME(NCOM)     /'CENTREOFMASS'  /
      DATA LNEDMD(NCOM)     /        .TRUE.  /
      DATA LNEDEM(NCOM)     /        .FALSE. /

      DATA BTNAME(NSTEP)    /        'STEP'  /
      DATA LNEDMD(NSTEP)    /        .TRUE.  /
      DATA LNEDEM(NSTEP)    /        .TRUE.  /

      DATA BTNAME(NSHAKE)   /       'SHAKE'  /
      DATA LNEDMD(NSHAKE)   /        .TRUE.  /
      DATA LNEDEM(NSHAKE)   /        .TRUE.  /

      DATA BTNAME(NFORCE)   /       'FORCE'  /
      DATA LNEDMD(NFORCE)   /        .TRUE.  /
      DATA LNEDEM(NFORCE)   /        .TRUE.  /

      DATA BTNAME(NPLIST)   /       'PLIST'  /
      DATA LNEDMD(NPLIST)   /        .TRUE.  /
      DATA LNEDEM(NPLIST)   /        .TRUE.  /

      DATA BTNAME(NPRINT)   /       'PRINT'  /
      DATA LNEDMD(NPRINT)   /        .TRUE.  /
      DATA LNEDEM(NPRINT)   /        .TRUE.  /

      DATA BTNAME(NWRITE)   /       'WRITE'  /
      DATA LNEDMD(NWRITE)   /        .FALSE. /
      DATA LNEDEM(NWRITE)   /        .FALSE. /

      DATA BTNAME(NPREST)   /     'POSREST'  /
      DATA LNEDMD(NPREST)   /       .FALSE.  /
      DATA LNEDEM(NPREST)   /       .FALSE.  /

      DATA BTNAME(NDREST)   /    'DISTREST'  /
      DATA LNEDMD(NDREST)   /       .FALSE.  /
      DATA LNEDEM(NDREST)   /       .FALSE.  /

      DATA BTNAME(NPERT)    /     'PERTURB'  /
      DATA LNEDMD(NPERT)    /       .FALSE.  /
      DATA LNEDEM(NPERT)    /       .FALSE.  /

      DATA BTNAME(NDHRES)   /    'DIHEREST'  /
      DATA LNEDMD(NDHRES)   /       .FALSE.  /
      DATA LNEDEM(NDHRES)   /       .FALSE.  /

      DATA BTNAME(NSTART)   /       'START'  /
      DATA LNEDMD(NSTART)   /        .TRUE.  /
      DATA LNEDEM(NSTART)   /        .TRUE.  /

      DATA BTNAME(NFOURD)   /      'FOURDIM' /
      DATA LNEDMD(NFOURD)   /        .FALSE. /
      DATA LNEDEM(NFOURD)   /        .FALSE. /

      DATA BTNAME(NLONGR)   /     'LONGRANGE'/
      DATA LNEDMD(NLONGR)   /         .TRUE. /
      DATA LNEDEM(NLONGR)   /         .TRUE. /

      DATA BTNAME(NEMBLK)   /      'MINIMISE'/
      DATA LNEDMD(NEMBLK)   /         .FALSE./
      DATA LNEDEM(NEMBLK)   /         .FALSE./

      DATA BTNAME(NSDBLK)   /    'STOCHASTIC'/
      DATA LNEDMD(NSDBLK)   /         .FALSE./
      DATA LNEDEM(NSDBLK)   /         .FALSE./

      DATA BTNAME(NPIBLK)   /       'PATHINT'/
      DATA LNEDMD(NPIBLK)   /         .FALSE./
      DATA LNEDEM(NPIBLK)   /         .FALSE./

      DATA BTNAME(NSYST)    /        'SYSTEM'/
      DATA LNEDMD(NSYST)    /          .TRUE./
      DATA LNEDEM(NSYST)    /          .TRUE./

      DATA BTNAME(NJVAL)    /         'J-VAL'/
      DATA LNEDMD(NJVAL)    /         .FALSE./
      DATA LNEDEM(NJVAL)    /         .FALSE./

      DATA BTNAME(NLEBLK)   /'LOCALELEVATION'/
      DATA LNEDMD(NLEBLK)   /         .FALSE./
      DATA LNEDEM(NLEBLK)   /         .FALSE./

      INCLUDE "inc_rdmd2.h"
C
C add other data statements here
C
      INTEGER I
C this is set to .TRUE. when the control file has been sucessfully
C read in
      LOGICAL LREDMD
      COMMON /MDLOC/LREDMD
      DATA LREDMD   /.FALSE. /

      DATA IDOPRO  /IDOMD/
C in block PCOUPLE
      DATA NTP /NTPOFF/
      DATA NRDBOX /NRDBXN/

C COM
      DATA NDFMIN /0/
      DATA NSCM   /0/
      DATA LTCM   /.FALSE./

C in block POSREST
      DATA NTR /NTROFF/

C in block DISTREST
      DATA NTDR /NTDROF/

C in block PERTURB
      DATA NTG /NTGOFF/

C in block DIHEREST
      DATA NTDLR /NTDLRN/

C in block FOURD
      DATA NT4DIM /NT4OFF/
      DATA NTCW4D /NTC4NO/
      DATA (NDO4D(I) ,I=1,MX4DIN)/MX4DIN* 3/

C in block LONGRANGE
      DATA EPSRF /1.0D0/
      DATA APPAK /1.0D0/
      DATA RCRF  /0.0D0/

C in path integral block
      DATA NTPI /NTPIOF/

C in local elevation
      DATA NTLE /NTLEOF/

C in j-vl restraining
      DATA NTJR /NTJROF/

C in write
      DATA NTWX  /0/
      DATA NTWSE /0/
      DATA NTWV  /0/
      DATA NTWE  /0/
      DATA NTWG  /0/
      DATA NTPW  /0/

C in SUBMOLECULES (not needed for EM)
C end blockdata mdinit
      DATA NSPM /0/

      DATA UNAME(JNMDCT)    /  'MDCT'/
      DATA IUNUM(JNMDCT)    /   IOMDCT/

      DATA UNAME(JNTRJX)    /  'TRJX' /
      DATA IUNUM(JNTRJX)    /   IOTRJX /

      DATA UNAME(JNTRJV)    /  'TRJV' /
      DATA IUNUM(JNTRJV)    /   IOTRJV /

      DATA UNAME(JNTRJE)    /  'TRJE' /
      DATA IUNUM(JNTRJE)    /   IOTRJE /

      DATA UNAME(JNTOPO)    /  'TOPO' /
      DATA IUNUM(JNTOPO)    /   IOTOPO/

      DATA UNAME(JNXVI)     /  'XVI' /
      DATA IUNUM(JNXVI)     /   IOXVI/

      DATA UNAME(JNGAM)     /  'GAMMA'/
      DATA IUNUM(JNGAM)     /   IOGAM/

      DATA UNAME(JNREST)    /  'RSTCOR'/
      DATA IUNUM(JNREST)    /   IOREST/

      DATA UNAME(JNRSTA)    /  'RSTNDX'/
      DATA IUNUM(JNRSTA)    /   IORSTA/

      DATA UNAME(JNRSTP)    /  'RSTPAIR'/
      DATA IUNUM(JNRSTP)    /   IORSTP/

      DATA UNAME(JNPERT)    /  'PERT'/
      DATA IUNUM(JNPERT)    /   IOPERT/

      DATA UNAME(JNDIHE)    /  'DIHESPEC'/
      DATA IUNUM(JNDIHE)    /   IODIHE/

      DATA UNAME(JNXVE)     /  'XVE'/
      DATA IUNUM(JNXVE)     /   IOXVE/

      DATA UNAME(JN4NDX)    /  '4DNDX'/
      DATA IUNUM(JN4NDX)    /   IO4NDX/

      DATA UNAME(JNCNV)     /  'CNV'/
      DATA IUNUM(JNCNV)     /   IOCNV/

      DATA UNAME(JNJVSP)    /  'JVALSPEC'/
      DATA IUNUM(JNJVSP)    /   IOJVSP/

      DATA UNAME(JNLESP)    /  'LESPEC'/
      DATA IUNUM(JNLESP)    /   IOLESP/

      DATA UNAME(JNTRJG)    /   'TRJG'/
      DATA IUNUM(JNTRJG)    /   IOTRJG/

C----
      DATA UNAME(JND00)     /  'IOD00'/
      DATA IUNUM(JND00)     /   IOD00/

      DATA UNAME(JND01)     /  'IOD01'/
      DATA IUNUM(JND01)     /   IOD01/

      DATA UNAME(JND02)     /  'IOD02'/
      DATA IUNUM(JND02)     /   IOD02/

      DATA UNAME(JND03)     /  'IOD03'/
      DATA IUNUM(JND03)     /   IOD03/

      DATA UNAME(JND04)     /  'IOD04'/
      DATA IUNUM(JND04)     /   IOD04/

      DATA UNAME(JND05)     /  'IOD05'/
      DATA IUNUM(JND05)     /   IOD05/

      DATA UNAME(JND06)     /  'IOD06'/
      DATA IUNUM(JND06)     /   IOD06/

      DATA UNAME(JND07)     /  'IOD07'/
      DATA IUNUM(JND07)     /   IOD07/

      DATA UNAME(JND08)     /  'IOD08' /
      DATA IUNUM(JND08)     /   IOD08/

      DATA UNAME(JND09)     /  'IOD09'/
      DATA IUNUM(JND09)     /   IOD09/

C     end

COMMVAR FMNII,FMNIR,FMGEI,FMGER,FMGTI,FMGTR,FMBTI,FMBTR,FMBOO,FMBTII,FMSHKP

C init my format statements defined in formats.h here
      DATA FMNII / '(1X,''illegal value for '',A,'' = '',I4)' /
      DATA FMNIR / '(1X,''illegal value for '',A,'' = '',F10.5)' /

      DATA FMGEI / '(1X,A,'' must be >= 0'')' /
      DATA FMGER / '(1X,A,'' must be >= 0.0'')' /

      DATA FMGTI / '(1X,A,'' must be > 0'')' /
      DATA FMGTR / '(1X,A,'' must be > 0.0'')' /


      DATA FMBTI  / '(1X,''must be '',I3,'' <= '',A,'' <= '',I3)' /
      DATA FMBTII / 
     $     '(1X,''must be '',I3,'' <= '',A,'' <= '',A,'' = '',I5)' /

      DATA FMBTR / '(1X,''must be '',F3.2,'' <= '',A,'' <= '',F3.2)' /

      DATA FMBOO / '(1X,A,''must be 0 or 1'')'/

      DATA FMSHKP/
     $'(1X,A,'': EMERGENCY STOP!'',/,''SHAKE FAILURE ON SOLUTE'',/)'/
      DATA FMSHKS/
     $'(1X,A,'': EMERGENCY STOP!'',/,''SHAKE FAILURE ON SOLVENT'',/)'/
COMMEND

      DATA FMRTAB(1)  /'(F1.1)'/
      DATA FMRTAB(2)  /'(F2.1)'/
      DATA FMRTAB(3)  /'(F3.1)'/
      DATA FMRTAB(4)  /'(F4.1)'/
      DATA FMRTAB(5)  /'(F5.1)'/
      DATA FMRTAB(6)  /'(F6.1)'/
      DATA FMRTAB(7)  /'(F7.1)'/
      DATA FMRTAB(8)  /'(F8.1)'/
      DATA FMRTAB(9)  /'(F9.1)'/
      DATA FMRTAB(10) /'(F10.1)'/
      DATA FMRTAB(11) /'(F11.1)'/
      DATA FMRTAB(12) /'(F12.1)'/
      DATA FMRTAB(13) /'(F13.1)'/
      DATA FMRTAB(14) /'(F14.1)'/
      DATA FMRTAB(15) /'(F15.1)'/
C
      DATA FMITAB(1)  /'(I1)'/
      DATA FMITAB(2)  /'(I2)'/
      DATA FMITAB(3)  /'(I3)'/
      DATA FMITAB(4)  /'(I4)'/
      DATA FMITAB(5)  /'(I5)'/
      DATA FMITAB(6)  /'(I6)'/
      DATA FMITAB(7)  /'(I7)'/
      DATA FMITAB(8)  /'(I8)'/
      DATA FMITAB(9)  /'(I9)'/
      DATA FMITAB(10) /'(I10)'/
C end 

C     LREDTO is initialized to .FALSE. here
C     and set to .TRUE. after a topology is
C     succesfully read in by RDTOPO.
      LOGICAL LREDTO
      COMMON /TOPOK/LREDTO

      DATA LREDTO /.FALSE./

C init the topology itself
      DATA    NRP /0/
      DATA  NRATT,NRATT2/0,0/
      DATA  NRAA2 /0/
      DATA   NBTY,NBON,NBONH /0,0,0/
      DATA   NTTY,NTHE,NTHEH /0,0,0/
      DATA   NQTY,NQHI,NQHIH /0,0,0/
      DATA   NPTY,NPHI,NPHIH /0,0,0/
      DATA   NCAG,NAEX,NAEX14/0,0,0/
      DATA   NRAM,NCONS      /0,0/

      DATA   FPEPSI,HBAR       /0.0D0,0.0D0/
      DATA   NPID /0/

C end blockdata
C begin
      DATA NRCON  /0/
      DATA NRRST  /0/
      DATA NDR    /0/
      DATA NDLR   /0/
      DATA NCONG  /0/
      DATA NLECFG /0/
      DATA NDJV   /0/
      DATA NLEMEM /0/
      DATA NDLE   /0/

C     The length of these strings must not be longer
C     than MAXIDL ( = 16) chars long.

      DATA NAMEID(IENDID)   / 'END' /
      DATA NAMEID(IFRMID)   / 'FORMAT' /
      DATA NAMEID(ITITID)   / 'TITLE' /
      DATA NAMEID(IBOXID)   / 'BOX' /

      DATA NAMEID(IPOSID)   / 'POSITION' /
      DATA NAMEID(IPORID)   / 'POSITIONRED' /
      DATA NAMEID(IP4DID)   / 'POSITION4THD'/
      DATA NAMEID(IP4RID)   / 'POSITION4THDRED'/
      DATA NAMEID(IPRCID)   / 'REFPOSITION' /
      DATA NAMEID(IPOBID)   / 'POSITIONOF' /

      DATA NAMEID(IVELID)   / 'VELOCITY' /
      DATA NAMEID(IVRDID)   / 'VELOCITYRED'/
      DATA NAMEID(IV4DID)   / 'VELOCITY4THD'/
      DATA NAMEID(IV4RID)   / 'VELOCITY4THDRED'/

      DATA NAMEID(ISDIID)   / 'STOCHINT'/
      DATA NAMEID(ISD4ID)   / 'STOCHINT4THD'/
      DATA NAMEID(IPRTID)   / 'PERTDATA'/
      DATA NAMEID(IRIIAV)   / 'DISRESEXPAVE'/
      DATA NAMEID(ICOSAV)   / 'JVALRESEXPAVE'/
      DATA NAMEID(ILEMID)   / 'LEMEMORY'/
      DATA NAMEID(ITIMID)   / 'TIMESTEP'/

      DATA NAMEID(IABFID)   / 'BFACTOR'/
      DATA NAMEID(IANBID)   / 'BFACTORANISO'/
      DATA NAMEID(IP2MID)   / 'POSITIONSECONDM'/
      DATA NAMEID(IP3MID)   / 'POSITIONTHIRDM'/
      DATA NAMEID(IP4MID)   / 'POSITIONFOURTHM'/
      DATA NAMEID(IP2TID)   / 'POSITIONSECONDMT'/

      DATA NAMEID(IQTAID)   / 'QUANTITYAVER' /
      DATA NAMEID(IQEAID)   / 'QUANENERAVER'/
      DATA NAMEID(IQESID)   / 'QUANSUMENERAVER'/
      DATA NAMEID(IQTSID)   / 'QUANTIMESERIES'/
      DATA NAMEID(IQDSID)   / 'QUANDISTRIB' /
      DATA NAMEID(IQTCID)   / 'QUANTIMECORR'/
      DATA NAMEID(IQSDID)   / 'QUANTIMECORRSPE'/

      DATA NAMEID(IXFRID)   / 'FRICTIONSPEC'/
      DATA NAMEID(IX4DID)   / 'FOURDIMATOMSPEC'/
      DATA NAMEID(IXPRID)   / 'POSRESSPEC' /

      DATA NAMEID(IXJVID)   / 'JVALRESSPEC'/
      DATA NAMEID(IXLEID)   / 'LEDIHSPEC'/
      DATA NAMEID(IDHBLK)   / 'DIHRESSPEC'/
      DATA NAMEID(IDRBLK)   / 'DISRESSPEC'/

      DATA NAMEID(INRGID)   / 'ENERGY'/
      DATA NAMEID(IDVPRT)   / 'VOLUMEPRESSURE'/
      DATA NAMEID(IDRLAM)   / 'FREEENERGYLAMBDA'/
      DATA NAMEID(IDRMU)    / 'FREEENERGY3D4'/

      DATA NAMEID(ISVSID)   / 'SOLVSTAT' /
      DATA NAMEID(IDFSID)   / 'DIFFSTAT' /
      DATA NAMEID(IDMSID)   / 'DIPMSTAT' /

C     format strings for writing out errors
      DATA FWFAIL / '(1X,A,'': write failed for '',A,'' on unit'',I3)'/
      DATA FRFAIL / '(1X,A,'': read failed for '',A,'' on unit'',I3)'/

      DATA FMEXP  / '(1X,A,'': '',A,'' expected on unit '',I3)'/
      DATA FMG2   / '(1X,A,'': '',A,'' encountered twice on unit'',I3)'/

      DATA ENDEXP /'END marker expected'/

      DATA EXPFIL /'(1X,''expected   :'',I5)'/
      DATA BONFIL /'(1X,''but on file:'',I5)'/
      DATA MAXVAL /'(1X,''max allowed:'',I5)'/

      DATA FLGBIN /'(1X,A6,'': Unit:'',I5)'/

      DATA MTYLIN
     $ /'(1X,''Read failed (skipped first'',I4,''CHARS on the line)'')'/

      DATA STREXP /'(1X,A,'': expected'')'/

C do not change this string lightly. Exactly 24 chars are written
C to file.
      DATA FMIGNO / '(''# first 24 chars ignored'')'/

C end blockdata

cmb      INTEGER I
C     --------------------------------------
C     make sure this common block is changed
C     consistently within this file
      real*8 EPLJ(MXNRE2),EPEL(MXNRE2),EPRF(MXNRE2),EPRC(MXNRE2)
      real*8 ENER(MXETBL),DEDLAM(MXETBL),DEDMU(MXETBL)
      real*8 ENERES(MXCTBL)
      real*8 XMIN(MAXDIM)
      real*8 VOLPRT(MXVTBL)
      COMMON /EMBLK/ EPLJ,EPEL,EPRF,EPRC,ENER,ENERES,
     $     DEDLAM,DEDMU,XMIN,VOLPRT
C     --------------------------------------

      DATA (XMIN(I), I=1,MAXDIM) /MAXDIM* 0.0D0/

      DATA PTBNAM (NPTTIT) / 'TITLE'/
      DATA PTBNAM (NPTATM) / 'PERTATOM'/
      DATA PTBNAM (NPTATP) / 'PERTATOMPAIR'/
      DATA PTBNAM (NPTBHG) / 'PERTBONDH'/
      DATA PTBNAM (NPTBNG) / 'PERTBOND'/
      DATA PTBNAM (NPTBAH) / 'PERTBANGLEH'/
      DATA PTBNAM (NPTBAG) / 'PERTBANGLE'/
      DATA PTBNAM (NPTQDH) / 'PERTIMPDIHEDRALH'/
      DATA PTBNAM (NPTQDG) / 'PERTIMPDIHEDRAL'/
      DATA PTBNAM (NPTPDH) / 'PERTDIHEDRALH'/
      DATA PTBNAM (NPTPDG) / 'PERTDIHEDRAL'/
      DATA PTBNAM (NPTPIN) / 'PERTPATHINT'/ 
      DATA PTBNAM (NPTEND) / 'END'/
C end blockdata ptbtin

      LOGICAL LREDPT
      COMMON /PTLOC/LREDPT
      DATA LREDPT /.FALSE./

      DATA   NJLA /0/
      DATA    NEB /0/
      DATA NBONHG /0/
      DATA  NBONG /0/
      DATA NTHEHG /0/
      DATA  NTHEG /0/
      DATA NQHIHG /0/
      DATA  NQHIG /0/
      DATA NPHIHG /0/
      DATA  NPHIG /0/

C end blockdata PTINIT

      DATA NDIM   /3/
      DATA NDHMIN /4/
      DATA NDHMAX /3/

#if defined(__USE_BLOCKDATA)
      END
#else
      RETURN 
      END
#endif
C     ==================================================================
      SUBROUTINE schifezza1
C     ==--------------------------------------------------------------==
cmb - Revised on 13 Oct. 2006
#if defined(__USE_BLOCKDATA)
      RETURN
      END
      BLOCKDATA SCHIFEZZA1DATA
#endif     
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
      INCLUDE 'box.h'
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'formats.h'
      INCLUDE 'md.h'
      INCLUDE 'pertsz.h'
      INCLUDE 'pertar.h'
      INCLUDE 'ptblock.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'units.h'
      INCLUDE 'topblock.h'

C init the block names
      DATA TPNAME(NPTITL) / 'TITLE'/
      DATA  LNEED(NPTITL) / .TRUE./

      DATA TPNAME(NATNAM) / 'ATOMTYPENAME'/
      DATA  LNEED(NATNAM) / .FALSE./

      DATA TPNAME(NRSNAM) / 'RESNAME'/
      DATA  LNEED(NRSNAM) / .FALSE./

      DATA TPNAME(NSOLAT) / 'SOLUTEATOM'/
      DATA  LNEED(NSOLAT) / .FALSE./

      DATA TPNAME(NBNDTY) / 'BONDTYPE'/
      DATA  LNEED(NBNDTY) / .FALSE. /

      DATA TPNAME(NBNDHY) / 'BONDH'/
      DATA  LNEED(NBNDHY) / .FALSE./

      DATA TPNAME(NBND)   / 'BOND'/
      DATA  LNEED(NBND)   / .FALSE./

      DATA TPNAME(NBATY)  / 'BONDANGLETYPE'/
      DATA  LNEED(NBATY)  / .FALSE./

      DATA TPNAME(NBAHY)  / 'BONDANGLEH'/
      DATA  LNEED(NBAHY)  / .FALSE./

      DATA TPNAME(NBAN)   / 'BONDANGLE'/
      DATA  LNEED(NBAN)   / .FALSE./

      DATA TPNAME(NIMDTY) / 'IMPDIHEDRALTYPE'/
      DATA  LNEED(NIMDTY) / .FALSE./

      DATA TPNAME(NIMDH)  / 'IMPDIHEDRALH'/
      DATA  LNEED(NIMDH)  / .FALSE./

      DATA TPNAME(NIMD)   / 'IMPDIHEDRAL'/
      DATA  LNEED(NIMD)   / .FALSE./

      DATA TPNAME(NDITY)  / 'DIHEDRALTYPE'/
      DATA  LNEED(NDITY)  / .FALSE./

      DATA TPNAME(NDIHY)  / 'DIHEDRALH'/
      DATA  LNEED(NDIHY)  / .FALSE./

      DATA TPNAME(NDIN)   / 'DIHEDRAL'/
      DATA  LNEED(NDIN)   / .FALSE./

      DATA TPNAME(NLJBLK) / 'LJPARAMETERS'/
      DATA  LNEED(NLJBLK) / .FALSE./

      DATA TPNAME(NTPEND) / 'END'/
      DATA  LNEED(NTPEND) / .FALSE./

      DATA TPNAME(NTPUNT) / 'TOPPHYSCON'/
      DATA  LNEED(NTPUNT) / .TRUE./

      DATA TPNAME(NSLVBK) / 'SOLVENTATOM'/
      DATA  LNEED(NSLVBK) / .FALSE./

      DATA TPNAME(NSLVCN) / 'SOLVENTCONSTR'/
      DATA  LNEED(NSLVCN) / .FALSE./

      DATA TPNAME(NTPVER) / 'TOPVERSION'/
      DATA  LNEED(NTPVER) / .TRUE./

      DATA TPNAME(NPITBL) / 'PATHINTSPEC'/
      DATA  LNEED(NPITBL) / .FALSE./

#if defined(__USE_BLOCKDATA)
      END
#else
      RETURN
      END
#endif   
C     ==--------------------------------------------------------------==
