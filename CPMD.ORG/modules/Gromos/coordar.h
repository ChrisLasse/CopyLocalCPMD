C hello emacs, this is -*- fortran -*-
C GROMOS common blocks for coordinates, restraining and local elevation

C -----------------------------------------------------------
C always include the file coordsz.h before including this one
C -----------------------------------------------------------

C     The variables defining the extent of the
C     arrays used (e.g. NRCON) in these common blocks
C     are initialized to 0 in a BLOCKDATA statement in posio.f
C     These limits are then modified on reading data from
C     file.

COMMVAR X,V,F,XR,SX,XSAVE
C     Coordinate arrays, their use depends on whether
C     Molecular Dynamics, Stochastic Dynamics or Energy Minimisation
C     is performed.
C     X       atom coordinates                 MD,SD,EM
C     F       atom forces                      MD,SD,EM
C     V       atom velocities in MD,SD, scratch array in EM
C     XR      for virial calculation in MD,SD, scratch array in EM
C Phil
C     XCEN    for the submolecule and solvent molecules center of mass
C             (MAXXCO is overdimensionned)
C
C     SX      stochastic integral array        SD
C     XSAVE   atom coordinates for saving      
C             lowest energy configurations     MD,SD
      real*8, pointer ::  X(:),V(:),F(:),XR(:),XCEN(:),SX(:),XSAVE(:)
      COMMON /COORD/ X,F,V,
     $     XR,XCEN,SX,XSAVE
COMMEND

COMMVAR CORTIT
C     CORTIT: a title read from the coordinate file
      CHARACTER *(MXCOTI) CORTIT
      COMMON /CRDTIT/ CORTIT
COMMEND

COMMVAR IAGRP
C     IAGRP is set in L<PROMD> using the values specified in the
C     array L<NRE>.
C     IAGRP is used for a faster lookup of energy contributions in
C     the force calculation routines.
C     In order to save space, however, it is defined as a character
C     array instead of an integer array.
C     The values are written into the array using CHAR() and read
C     using ICHAR(), both of which are standard F77 functions.

      CHARACTER, pointer ::   IAGRP(:)
      COMMON /GROUPY/IAGRP
COMMEND

COMMVAR JRC,NRCON,NRRST
C*****position re(con)straining
C     We can have EITHER restraining OR constraining.
C     In the event of restraining : NRCON = 0 and NRRST > 0
C     In the event of constraining: NRCON > 0 and NRRST = 0
C
C     In both cases the array JRC holds the indices of re(con)strained
C     atoms. However, the allowed range is different:
C     for constraining JRC(I) 1..L<NRP>
C     for restraining  JRC(I) 1..L<NATTOT>
C
      INTEGER NRCON,NRRST
      INTEGER, pointer :: JRC(:)
      COMMON /POSINT_p/ JRC
      COMMON /POSINT/ NRCON,NRRST
COMMEND

COMMVAR XC,CXC,XC0,EC0
C*****position re(con)restraining
C     XC  : reference positions
C     CXC : force constants used in the position restraining
C     
      real*8, pointer :: XC(:),CXC(:)
      COMMON /POSR/ 
     $     XC,CXC
COMMEND


COMMVAR IDR1,JDR1,KDR1,LDR1,ICDR1,IDR2,JDR2,KDR2,LDR2,ICDR2,NDR
C*****atom distance restraining
      INTEGER, pointer ::  IDR1(:),JDR1(:),KDR1(:),LDR1(:),ICDR1(:)
      INTEGER, pointer ::  IDR2(:),JDR2(:)
      INTEGER, pointer ::  KDR2(:),LDR2(:),ICDR2(:)
      INTEGER NDR
      
      COMMON /DISRIN_p/
     $     IDR1,JDR1,KDR1,LDR1,
     $     ICDR1,
     $     IDR2,JDR2,KDR2,LDR2,
     $     ICDR2
      COMMON /DISRIN/
     $     NDR

COMMEND

COMMVAR R0,W0,RIIAVE,DISH,DISC
      real*8, pointer ::  R0(:),W0(:),RIIAVE(:)
      real*8 DISH,DISC
      COMMON /DISRFP_p/ R0,W0,
     $     RIIAVE
      COMMON /DISRFP/ DISH,DISC
C     RIIAVE is used in time averaged distant restraining
COMMEND


COMMVAR IPLR,JPLR,KPLR,LPLR,ICPLR,NDLR,CPLR,PDLR
C*****dihedral restraining
      INTEGER,  pointer ::   IPLR(:),JPLR(:),KPLR(:),LPLR(:),ICPLR(:)
      INTEGER NDLR 
      COMMON /DIHR_p/ IPLR,JPLR,KPLR,
     $     LPLR,ICPLR
      COMMON /DIHR/ 
     $     NDLR

      real*8,  pointer ::   CPLR(:),PDLR(:)
      COMMON /DIHRFP/ CPLR,PDLR
COMMEND

COMMVAR ICOG,JCOG,NCONG
C     distance constraining
      INTEGER, pointer :: ICOG(:),JCOG(:)
      INTEGER NCONG
      COMMON /DICO_p/ICOG,JCOG
      COMMON /DICO/NCONG
COMMEND

COMMVAR CONP,FCON
C     FCON is used to store and manipulate the
C     constraint forces
      real*8,  pointer ::  CONP(:),FCON(:)
      COMMON /DICOFP/ CONP,FCON
COMMEND

COMMVAR IPJV,JPJV,KPJV,LPJV,NDJV,CPJV,PJR0,PSJR,AJV,BJV,CJV,COSQAV,COSIAV
C     j-value restraining
      INTEGER NDJV
      INTEGER, pointer :: IPJV(:), JPJV(:), KPJV(:), LPJV(:)
      COMMON /JVINTY/NDJV
      COMMON /JVINTY_p/
     $     IPJV, JPJV, KPJV, LPJV

      real*8, pointer :: CPJV(:),PJR0(:),PSJR(:),AJV(:),BJV(:),CJV(:)
      real*8, pointer :: COSQAV(:),COSIAV(:)
      COMMON /JVREL/ CPJV,PJR0,PSJR,
     $     AJV,BJV,CJV,
     $     COSQAV,COSIAV
COMMEND


COMMVAR NDLE,IPLE,JPLE,KPLE,LPLE,NLECFG,NLEMEM
C*****local elevation

C     The local elevation dihedral angles.
      INTEGER NDLE
      INTEGER, pointer :: IPLE(:),JPLE(:),KPLE(:),LPLE(:)
      INTEGER, pointer :: ILEMEM(:,:)
      INTEGER, pointer :: NLEVST(:)

C     NLECFG: the number of INTS used to store one configuration
C     NLECFG <= L<MXLECF>
C     NLEMEM: the number of configurations in memory
C     NLEMEM <= L<MLECFG>

      INTEGER NLECFG,NLEMEM

      COMMON /LEINTY_p/
     $     IPLE,JPLE,KPLE,LPLE,
     $     ILEMEM,NLEVST
      COMMON /LEINTY/NDLE,NLECFG,NLEMEM
COMMEND



COMMVAR NSP,NSPM
C     NSP(I) contains the last atom of submolecule I in the solute
C     NSPM: the number of submolecules in the solute
C
C     NSP is defined as MAXNSP+1 in order to allow more efficient code
C     in the nonbonded force calculation routines.
C     (a double IF statement is avoided in the virial calculation)
C     However, only the elements 1..MAXNSP are ever actually filled with
C     useful data.
      INTEGER, pointer :: NSP(:)
      INTEGER NSPM
      COMMON /MNSP_p/ NSP
      COMMON /MNSP/ NSPM
COMMEND


COMMVAR C4D
C     Force consts of 4-th dim harmonic osc.
C     If the constants are positive or zero for an atom,
C     it is in 4D.
C     If it is strictly negative, the atom is in 3D.

      real*8, pointer :: C4D(:)
      COMMON /FORCFP/C4D
COMMEND

COMMVAR GAM,CC1,CC2,CC3,CC4,CC5,CC6,CC7,CC8,CC9,SWINK,SWINKS
C these used in SD
      real*8, pointer :: GAM(:),CC1(:),CC2(:),CC3(:),CC4(:)
      real*8, pointer :: CC5(:),CC6(:),CC7(:),CC8(:),CC9(:)
      real*8, pointer :: SWINK(:),SWINKS(:)
      COMMON /SDFRIC/GAM,
     $     CC1,CC2,CC3,CC4,
     $     CC5,CC6,CC7,CC8,
     $     CC9,SWINK,SWINKS
COMMEND



