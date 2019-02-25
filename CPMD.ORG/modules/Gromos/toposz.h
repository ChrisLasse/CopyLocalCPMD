C hello emacs, this is -*- fortran -*-
C Version not-to-be-changed by AK (2005)
C include file for GROMOS size limits on topology

C maximum number of atom types
      INTEGER MAXATT,MXATT2
C maximum number of residues
      INTEGER MAXAA2
C maximum number of atoms per solute
      INTEGER MAXNRP,MAXNP2
C maximum number of covalent bond types
      INTEGER MAXNBT
C maximum number of bonds involving H-atoms in the solute
      INTEGER MAXBNH
C maximum number of bonds NOT involving H-atoms in the solute
      INTEGER MAXBON
C maximum number of bond angle types
      INTEGER MAXTTY
C maximum number of bond angles involving H-atoms in the solute
      INTEGER MXQHEH
C maximum number of bond angles NOT involving H-atoms in the solute
      INTEGER MAXTHE
C maximum number of improper dihedral types
      INTEGER MAXQTY
C maximum number of improper dihedrals involving H-atoms in the solute
      INTEGER MAXHIH
C maximum number of improper dihedrals NOT involving H-atoms in the solute
      INTEGER MAXQHI
C maximum number of dihedral types
      INTEGER MAXPTY
C maximum number of dihedrals involving H-atoms in the solute
      INTEGER MXPHIH
C maximum number of dihedrals NOT involving H-atoms in the solute
      INTEGER MAXPHI
C maximum number of charge groups in a solute molecule
      INTEGER MAXCAG
C maximum total number of exclusions in a solute molecule
      INTEGER MAXAEX
C maximum number of third number atoms in a solute molecule
      INTEGER MXEX14
C maximum number of atoms per solvent molecule
      INTEGER MAXNRS
C maximum number of solvent constraints
      INTEGER MXCONS

COMMVAR MAXTIT,MAXLNS
C params defining the size of character arrays
C length of title string, and max number of lines allowed
      INTEGER MAXTIT, MAXLNS
      PARAMETER (MAXTIT = 80, MAXLNS = 10)
COMMEND

COMMVAR MAXTLE,MAXRLE,MAXNLE,MXNLE2,MXNLE3,MXNLE4
C----IT IS NOT ADVISED TO CHANGE THE VALUES OF THESE CONSTANTS
C length of atom type names
      INTEGER MAXTLE
      PARAMETER (MAXTLE = 5)

C length of residue type names
      INTEGER MAXRLE
      PARAMETER (MAXRLE = 5)

C length of atom name of solvent and solute atoms
      INTEGER MAXNLE
      PARAMETER (MAXNLE = 5)

C these used for pretty printing...
      INTEGER MXNLE2
      PARAMETER (MXNLE2 = 2*MAXNLE + 1)

      INTEGER MXNLE3
      PARAMETER (MXNLE3 = 3*MAXNLE + 2)

      INTEGER MXNLE4
      PARAMETER (MXNLE4 = 4*MAXNLE + 3)
COMMEND


COMMVAR MAXPIA,MAXPID,MAXPIT,MAXPIB,MAXPIW,MAXWR2
C-----PATH INTEGRAL
C maximum number of discretized atoms
      INTEGER MAXPIA
      PARAMETER (MAXPIA = 100)
 
C maximum number of discretizations
      INTEGER MAXPID
      PARAMETER (MAXPID = 10)
 
C maximum number of path integral 'bond' types, should be MAXPIA
      INTEGER MAXPIT
      PARAMETER (MAXPIT = MAXPIA)
 
C maximum number of path integral 'bonds': MAXPIA*MAXPID
      INTEGER MAXPIB
      PARAMETER (MAXPIB = 1000)
 
C maximum dimension of work arrays
      INTEGER MAXPIW
      PARAMETER (MAXPIW = 1000)
 
C maximum number of atoms forming a bond
      INTEGER MAXWR2
      PARAMETER (MAXWR2 = 4)
COMMEND

C number of parameters in the TOPOSZ common
      INTEGER NRTOPO
      PARAMETER (NRTOPO = 22)
      COMMON /TOPOSZ/ MAXATT,MXATT2,MAXAA2,MAXNRP,MAXNP2,MAXNBT,MAXBNH,
     & MAXBON,MAXTTY,MXQHEH,MAXTHE,MAXQTY,MAXHIH,MAXQHI,MAXPTY,MXPHIH,
     & MAXPHI,MAXCAG,MAXAEX,MXEX14,MAXNRS,MXCONS

