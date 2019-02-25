C hello emacs, this is -*- fortran -*-
C     topoar.h
C     common block definitions for the GROMOS topology
C
C----------------------------------------------------------
C     Note: include the file toposz.h before including
C     this one.
C----------------------------------------------------------


COMMVAR IAC,IPERT,NRP,MPAC,NRATT,NRATT2
C*****atom types, atom interaction types
C     IAC   : integer atom code
C     IPERT : denotes whether an atom is perturbed or not
C     IPERT is NOT read in by L<RDTOPO>, but is set to L<NOPERT>.
C     The actual reading is performed (optionally) by
C     L<RDPERT> if perturbation is specified in the L<PROMD>
C     control file.
C
C     NRP   : number of atoms per solute molecule
C     MPAC  : interaction matrix
C     NRATT : number of atom types
C     NRATT2: is always NRATT*(NRATT+1)/2
C
      INTEGER NRP,NRATT,NRATT2
      integer, pointer :: IAC(:)
      integer, pointer :: IPERT(:)
      integer, pointer :: MPAC(:,:)
      COMMON/ATYP/
     $     NRP,
     $     NRATT,NRATT2
      COMMON/ATYP_p/
     $     IAC,IPERT,MPAC
COMMEND

COMMVAR NOPERT
C     a value used to denote that an atom is not perturbed.
C     See array L<IPERT>
      INTEGER NOPERT
      PARAMETER (NOPERT = 0)
COMMEND

COMMVAR FFTYPE,TOPTIT,NTPLNS
C*****title, residue names
C     FFTYPE is an array of length MAXATT of char(MAXTLE)
C     TOPTIT: the title read in from the molecular topology file.
C     NTPLNS: the number of lines in the title
C     FFTYPE: the names of the residues
C
      CHARACTER*(MAXTLE), pointer :: FFTYPE(:)
      CHARACTER*(MAXTIT), pointer :: TOPTIT(:)
      COMMON /ATYPCH/
     $     FFTYPE ,
     $     TOPTIT
      INTEGER NTPLNS
      COMMON /ATPINT/NTPLNS
COMMEND

COMMVAR C12,C6,CS12,CS6
C     C6, C12: normal interaction parameters
C     CS6,CS12: 1-4 interaction parameters
C
      real*8, pointer :: C12(:)
      real*8, pointer :: C6(:)
      real*8, pointer :: CS12(:)
      real*8, pointer :: CS6(:)
      COMMON /ATYPFP/ 
     $     C12,C6,
     $     CS12,CS6
COMMEND


COMMVAR PANM
C     PANM: the names of the solute atoms
C     PANM is an array of length MAXNRP of char(MAXNLE)
      CHARACTER*(MAXNLE), pointer :: PANM(:)
      COMMON /TOP1/
     $    PANM
COMMEND



COMMVAR NRAA2, AANM
C*****residue types
C     NRAA2: the number of residues in the topology
C     AANM : the names of the residues
C     AANM is an array of length MAXAA2 of char(MAXRLE)

      INTEGER NRAA2
      COMMON /RTYP/ NRAA2

      CHARACTER*(MAXRLE), pointer ::  AANM(:)
      COMMON /RTYPCH/AANM
COMMEND


COMMVAR NBTY, CB,B0,IB,JB,ICB,NBON,IBH,JBH,ICBH,NBONH
C*****bonds
C     NBTY: number of bond types
      INTEGER NBTY
      real*8, pointer ::  CB(:)
      real*8, pointer ::  B0(:)
      COMMON /BONDFP/
     $     CB,B0

      INTEGER NBON
      INTEGER, pointer :: IB(:)
      INTEGER, pointer :: JB(:)
      INTEGER, pointer :: ICB(:)
      COMMON /BOND_p/
     $     IB,JB,ICB
      COMMON /BOND/
     $     NBON,NBTY

      INTEGER NBONH
      INTEGER, pointer :: IBH(:)
      INTEGER, pointer :: JBH(:)
      INTEGER, pointer :: ICBH(:)
      COMMON /BONH_p/IBH,JBH,ICBH
      COMMON /BONH/NBONH
COMMEND


COMMVAR NTTY,CT,T0,ITH,JTH,KTH,ICTH,NTHEH,IT,JT,KT,ICT,NTHE
C*****angles
      INTEGER NTTY
      real*8, pointer ::  CT(:)
      real*8, pointer ::  T0(:) 
      COMMON /ANGLFP/CT,T0

      INTEGER NTHEH
      INTEGER NTHE
      INTEGER, pointer :: ITH(:)
      INTEGER, pointer :: JTH(:)
      INTEGER, pointer :: KTH(:)
      INTEGER, pointer :: ICTH(:)
      INTEGER, pointer :: IT(:)
      INTEGER, pointer :: JT(:)
      INTEGER, pointer :: KT(:)
      INTEGER, pointer :: ICT(:)
      COMMON /ANGL_p/
     $     ITH,JTH,KTH,ICTH,
     $     IT,JT,KT,ICT
      COMMON /ANGL/
     $     NTHEH,
     $     NTHE,NTTY
COMMEND


COMMVAR NQTY,CQ,Q0,IQ,JQ,KQ,LQ,ICQ,NQHI,IQH,JQH,KQH,LQH,ICQH,NQHIH
C*****improper dihedrals
      INTEGER NQTY
      real*8, pointer ::  CQ(:)
      real*8, pointer ::  Q0(:) 
      COMMON /IDHAFP/ CQ,Q0

      INTEGER NQHI
      INTEGER, pointer :: IQ(:)
      INTEGER, pointer :: JQ(:)
      INTEGER, pointer :: KQ(:)
      INTEGER, pointer :: LQ(:)
      INTEGER, pointer :: ICQ(:)
      COMMON/IDHA_p/
     $     IQ,JQ,KQ,LQ,
     $     ICQ
      COMMON/IDHA/
     $     NQHI,NQTY

      INTEGER NQHIH
      INTEGER, pointer :: IQH(:)
      INTEGER, pointer :: JQH(:)
      INTEGER, pointer :: KQH(:)
      INTEGER, pointer :: LQH(:)
      INTEGER, pointer :: ICQH(:)
      COMMON/IDHH_p/
     $     IQH,JQH,KQH,LQH,
     $     ICQH
      COMMON/IDHH/
     $     NQHIH
COMMEND


COMMVAR NP,NPTY,CP,PD,IP,JP,KP,LP,ICP,NPHI,IPH,JPH,KPH,LPH,ICPH,NPHIH
C*****dihedrals
      INTEGER NPTY
      real*8, pointer ::  CP(:)
      real*8, pointer ::  PD(:)
      INTEGER, pointer ::  NP(:)
      COMMON /DIHAFP/ CP,PD
      COMMON /DIHTYP_p/ NP
      COMMON /DIHTYP/ NPTY


      INTEGER NPHI
      INTEGER, pointer :: IP(:)
      INTEGER, pointer :: JP(:)
      INTEGER, pointer :: KP(:)
      INTEGER, pointer :: LP(:)
      INTEGER, pointer :: ICP(:)
      COMMON/DIHA_p/
     $     IP,JP,KP,LP,
     $     ICP
      COMMON/DIHA/
     $     NPHI


      INTEGER NPHIH
      INTEGER, pointer :: IPH(:)
      INTEGER, pointer :: JPH(:)
      INTEGER, pointer :: KPH(:)
      INTEGER, pointer :: LPH(:)
      INTEGER, pointer :: ICPH(:)
      COMMON/DIHH_p/
     $     IPH,JPH,KPH,LPH,
     $     ICPH
      COMMON/DIHH/
     $     NPHIH
COMMVAR


COMMVAR MRES,INC,NCAG,INE,KNE,INE14,KNE14,JSNE,NAEX,JSNE14,NAEX14
C*****residue numbers,mass, charge, exclusions, 1-4 interactions,
C     charge group definitions, interaction matrix (mpac)
C---
C     The variables INE,KNE and JSNE are organized as follows:
C     exclusions J of atom I are positioned at
C     JSNE(KNE(I)+1),...JSNE(KNE(I)+INE(I))
C     all J must be > I and in ascending order.
C
C     The variables INE14,KNE14 and JSNE14 are analogous.
C---

      INTEGER NCAG
      INTEGER NAEX,NAEX14
      INTEGER, pointer :: MRES(:)
      INTEGER, pointer :: INC(:)
      INTEGER, pointer :: INE(:)
      INTEGER, pointer :: KNE(:)
      INTEGER, pointer :: INE14(:)
      INTEGER, pointer :: KNE14(:)
      INTEGER, pointer :: JSNE(:)
      INTEGER, pointer :: JSNE14(:)
      COMMON/NON2_p/
     $     MRES, INC,INE,KNE,INE14,KNE14, JSNE,JSNE14
      COMMON/NON2/
     $     NCAG, NAEX,NAEX14
COMMEND

COMMVAR CG,WINV,WMAS
C     CG  : charge of solute atoms
C     WINV: inverse mass of solute atoms
C     WMAS: mass of solute atoms
      real*8, pointer :: CG(:)
      real*8, pointer :: WINV(:)
      real*8, pointer :: WMAS(:)
      COMMON /NON2FP/
     $     CG,WINV,WMAS
COMMEND


COMMVAR IACS,NRAM,CGS,WINVS,WMASS,ANMS
C*****solvent molecule atom data
C     IACS : integer atom code for solvent atoms
C     NRAM : number of atoms in a solvent molecule
C     CGS  : charge of solvent atoms
C     WINVS: inverse mass of solvent atoms
C     WMASS: mass of solvent atoms
C     ANMS : names of solvent atoms

      INTEGER   NRAM
      INTEGER, pointer :: IACS(:)
      COMMON/SOLV_p/ IACS
      COMMON/SOLV/ NRAM

      real*8, pointer :: CGS(:)
      real*8, pointer :: WINVS(:)
      real*8, pointer :: WMASS(:)
      COMMON /SOLVFP/ CGS,WINVS,WMASS

      CHARACTER*(MAXRLE), pointer :: ANMS(:)
      COMMON /SOLVCH/ ANMS
COMMEND


COMMVAR ICONS,JCONS,NCONS,CONS
C*****solvent molecule constraints data
      INTEGER NCONS
      INTEGER, pointer :: ICONS(:)
      INTEGER, pointer :: JCONS(:)
      COMMON/SOCON_p/
     $     ICONS,JCONS
      
      COMMON/SOCON/
     $     NCONS
      
      real*8, pointer :: CONS(:)
      COMMON /SCONFP/ CONS
COMMEND

COMMVAR FPEPSI,HBAR
C*****units defined in the topology file
C     EPS0 is the permittivity of vacuum
C     FPEPSI: we store 1/(4* PI* EPS0) 
C     HBAR  : Planck's constant HBAR = H/(2* PI)
      real*8 FPEPSI,HBAR
      COMMON /TOPNTS/FPEPSI,HBAR
COMMEND

COMMVAR NPIA,NPIA,NPID,IPIC,NPIT,NPIB,IPIB,JPIB,ICPIB,TPI,CPI,WMCL,BOLPI
C*****path integral topology data
C     NPIA:     number of path integral atoms
C     IPIA:     atom seq. numbers in PI atoms in the classical topology
C     NPID:     number of discretizations per atom
C     IPIC(I):  pseudoparticle number of the 'atom' I: if zero, the
C           'atom' is a classical atom
C     TPI:      temperature as basis for the harmonic 'spring' constants
C     NPIT:     number of harmonic 'spring' constants between
C           pseudoparticles
C     CPI(I):   harmonic 'spring' constant of code I
C     BOLPI:    Boltzmann constant for path integral
C     NPIB:     number of harmonic 'bonds' between pseudoparticles
C     IPIB(I):  'from' atom of 'bond' I
C     JPIB(I):  'to' atom of 'bond' I
C     ICPIB(I): 'bond' code of 'bond' I
 
      INTEGER NPIA, NPID, NPIT, NPIB
      INTEGER, pointer :: IPIA(:)
      INTEGER, pointer :: IPIC(:)
      INTEGER, pointer :: IPIB(:)
      INTEGER, pointer :: JPIB(:)
      INTEGER, pointer :: ICPIB(:)
      real*8 TPI, BOLPI
      real*8, pointer :: CPI(:)
      real*8, pointer :: WMCL(:)
      COMMON /PITINT_p/ IPIA, IPIC, IPIB, JPIB,
     $     ICPIB
      COMMON /PITINT/ NPIA, NPID, NPIT, NPIB
      COMMON /PITREL/ TPI,  BOLPI
      COMMON /PITREL_p/  CPI, WMCL
COMMEND

