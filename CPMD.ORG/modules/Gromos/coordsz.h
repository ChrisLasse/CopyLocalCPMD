C hello emacs, this is -*- fortran -*-
C     size limits of coordinate and related arrays
Cmb   enlarged due to insufficiency in several applications 


COMMVAR MAXDIM
C     the maximum number of dimensions we 
C     can calculate in.
C     the actual number of dimensions used is given by
C     L<NDIM> in md.h at run time.

      INTEGER MAXDIM
      PARAMETER (MAXDIM = 4)
COMMEND


COMMVAR MAXNSP
C     maximum number of submolecules per solute molecule
      INTEGER MAXNSP
cmb     PARAMETER (MAXNSP = 10)
      PARAMETER (MAXNSP = 80)
COMMEND

COMMVAR MXSBCO
C     maximum number of (solute) submolecule coordinates in total
C     I.e. NPM*NSPM*NDIM must be smaller than this value
C     L<MAXNSP>, L<MAXDIM>
C
      INTEGER MXSBCO
cmb      PARAMETER (MXSBCO = 20*MAXNSP*MAXDIM)
      PARAMETER (MXSBCO = 80*MAXNSP*MAXDIM)
COMMEND


C     max length of the coord title
      INTEGER MXCOTI
      PARAMETER (MXCOTI = 80)

COMMVAR MAXNAT
C     maximum total number of atoms in the system
      INTEGER MAXNAT

COMMVAR MAXXCO
C     This determines the array sizes for the
C     arrays holding atom coords such as X, F, V etc.
C     Choose this to be >= L<NDIM>*L<MAXNAT>, where NDIM is the
C     number of dimensions (3 or 4) you calculate in.

      INTEGER MAXXCO

COMMVAR MAXPER
C     maximum number of perturbed atoms
      INTEGER MAXPER
      PARAMETER (MAXPER = 5)

C     maximum number of atoms in the whole system that can be
C     in four dimensions.
C     This value can safely be set to 0 if one
C     only simulates in three dimensions and wants to save space.
C     OTHERWISE IT MUST BE SET TO L<MAXNAT>

      INTEGER MAX4DA

C     atom position re(con)straining
C     maximum number of atom position restraints
      INTEGER MAXXC

C     distance restraining
C     maximum number of distance restraints PER SOLUTE MOLECULE
      INTEGER MAXNDR

C     maximum TOTAL number of time averaged distance restraints.
C     In time averaging, each solute molecule must have its own
C     averages. Therefore this must be at least NPM*NDR.
C     Just set this to MAXNDR for now.
      INTEGER MAXTDR

C     maximum number of restrained dihedrals IN TOTAL
      INTEGER MAXNDL

C     maximum number of distance constraints PER SOLUTE MOLECULE
      INTEGER MAXCON

C     maximum number of J-coupling constant restraints PER SOLUTE MOLECULE.
C     If more than one solute molecule is specified (NPM > 1), the
C     averages are taken over all molecules.
      INTEGER MAXJ

C     dihedral angle monitoring
C     The maximum number of dihedrals that can be monitored IN TOTAL
C
      INTEGER MXDMON
cmb      PARAMETER (MXDMON  = 1000)
      PARAMETER (MXDMON  = 4000)
COMMEND



COMMVAR BITPDH,VALPDH,BITPBY,BYTPIN,BITPIN,GRDSZE,GRDSZ2,DIPINT
C*****local elevation
C     The number of bits needed to encode a dihedral angle.
C     With 4 bits (a nibble) we discretize each dihedral angle
C     into 16 values.
C     Don't change this value
      INTEGER BITPDH
      PARAMETER (BITPDH = 4)

      INTEGER VALPDH
      PARAMETER (VALPDH = 2**BITPDH)

C     the number of BITS per BYTE
C     Don't change this value.
      INTEGER BITPBY
      PARAMETER (BITPBY = 8)

C     the number of BYTES (= eight bits) per INTEGER
      INTEGER BYTPIN
C     assume we have 32 bit integers
C     THIS MIGHT HAVE TO BE CHANGED, DEPENDING ON MACHINE AND/OR
C     COMPILE OPTIONS
      PARAMETER (BYTPIN = 4)

C the number of BITS per INTEGER
      INTEGER  BITPIN
      PARAMETER (BITPIN = BYTPIN*BITPBY)

C the discretisation size of the complete circle.
      real*8 GRDSZE, GRDSZ2
      PARAMETER (GRDSZE=360.0d0/VALPDH,GRDSZ2=GRDSZE*0.5d0)


C     the number of dihedral angles that can be stored into one integer.
C     in order to avoid problems with signed integers in fortran, 
C     we don't use the BITPDH most significant bits.
C
      INTEGER DIPINT
      PARAMETER (DIPINT = (BYTPIN*BITPBY /BITPDH) -1)
COMMEND

COMMVAR MAXDLE,MXLECF,MXLECF
C     the maximum number of dihedral angles that can define
C     a configuration. Make this an integer multiple of L<DIPINT>
C     in order to avoid wasting space in the hash table
C     for the memory. (see the expression for L<NLECFG> below)
      INTEGER MAXDLE
cmb      PARAMETER (MAXDLE = 21)
      PARAMETER (MAXDLE = 121)

C     the maximum number of integers needed to store one configuration
      INTEGER MXLECF
      PARAMETER (MXLECF = (MAXDLE + DIPINT -1)/DIPINT)


C     the maximum number of configurations we want to be able
C     to store in our memory.
C     This can safely be set to zero when not using local elevation.
      INTEGER MLECFG
C      PARAMETER (MLECFG  = 1000000)
      PARAMETER (MLECFG = 1)
COMMEND



C-----these are sizes of arrays used in analysis tools and define
C     vectors of second, third and fourth moments

COMMVAR M2XX,M2XY,M2XZ,M2YY,M2YZ,M2ZZ,M2MAX
C     for second moment and anisotropic B -factors
      INTEGER M2XX,M2XY,M2XZ,M2YY,M2YZ,M2ZZ
      INTEGER M2MAX
      PARAMETER (M2XX = 1)
      PARAMETER (M2XY = 2)
      PARAMETER (M2XZ = 3)
      PARAMETER (M2YY = 4)
      PARAMETER (M2YZ = 5)
      PARAMETER (M2ZZ = 6)
      PARAMETER (M2MAX= M2ZZ)
COMMEND


COMMVAR M3XXX,M3XXY,M3XXZ,M3XYZ,M3XZZ,M3YYY,M3YYZ,M3YZY,M3YZZ,M3MAX
C     for third moment
      INTEGER M3XXX,M3XXY,M3XXZ,M3XYZ,M3XZZ
      INTEGER M3YYY,M3YYZ,M3YZY,M3YZZ
      INTEGER M3ZZZ
      INTEGER M3MAX

      PARAMETER (M3XXX = 1)
      PARAMETER (M3XXY = 2)
      PARAMETER (M3XXZ = 3)
      PARAMETER (M3XYZ = 4)
      PARAMETER (M3XZZ = 5)
      PARAMETER (M3YYY = 6)
      PARAMETER (M3YYZ = 7)
      PARAMETER (M3YZY = 8)
      PARAMETER (M3YZZ = 9)
      PARAMETER (M3ZZZ =10)
      PARAMETER (M3MAX = M3ZZZ)
COMMEND

COMMVAR M4XXXX,M4XXXY,M4XXXZ,M4XXYY,M4XXYZ,M4XXZZ,M4XYYY,M4MAX
C     for fourth moment
      INTEGER M4XXXX,M4XXXY,M4XXXZ,M4XXYY,M4XXYZ,M4XXZZ
      INTEGER M4XYYY,M4XYYZ,M4XYZZ,M4XZZZ
      INTEGER M4YYYY,M4YYYZ,M4YYZZ,M4YZZZ,M4ZZZZ
      INTEGER M4MAX

      PARAMETER (M4XXXX =  1)
      PARAMETER (M4XXXY =  2)
      PARAMETER (M4XXXZ =  3)
      PARAMETER (M4XXYY =  4)
      PARAMETER (M4XXYZ =  5)
      PARAMETER (M4XXZZ =  6)
      PARAMETER (M4XYYY =  7)
      PARAMETER (M4XYYZ =  8)
      PARAMETER (M4XYZZ =  9)
      PARAMETER (M4XZZZ = 10)
      PARAMETER (M4YYYY = 11)
      PARAMETER (M4YYYZ = 12)
      PARAMETER (M4YYZZ = 13)
      PARAMETER (M4YZZZ = 14)
      PARAMETER (M4ZZZZ = 15)
      PARAMETER (M4MAX = M4ZZZZ)
COMMEND
      INTEGER NRCOORD
      PARAMETER (NRCOORD = 9)
      COMMON /COORDSZ/ MAXNAT, MAXXCO, MAX4DA, MAXXC,
     &     MAXNDR, MAXTDR, MAXNDL, MAXCON, MAXJ

