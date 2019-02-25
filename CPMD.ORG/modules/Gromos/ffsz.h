C hello emacs, this is -*- fortran -*-
C sizes for force field arrays
C

C only define things here which are not defined in toposz.h

COMMVAR MAXIMC
C     maximum number of integer mass codes
      INTEGER MAXIMC
      PARAMETER (MAXIMC = 21)
COMMEND

COMMVAR MAXC12
C     for C12 interaction parameters
C     do not change this value
      INTEGER MAXC12
      PARAMETER (MAXC12 = 3)
COMMEND


COMMVAR MAXMIX
C     maximum number of mixed atom pairs
      INTEGER MAXMIX
      PARAMETER (MAXMIX = 10)
COMMEND



