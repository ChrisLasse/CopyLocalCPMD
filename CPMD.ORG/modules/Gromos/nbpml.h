C hello emacs, this is -*- fortran -*-
C common blocks for SUBR. NBPML

C-------------------------------------------
C INCLUDE forcesz.h BEFORE THIS FILE!
C-------------------------------------------


C     
C     COMMON BLOCK FOR NBPML
C     
      LOGICAL LPERTL
      LOGICAL LDOTRA,LMONO,LOCTO,LVAC,LDOVIR,L4D
      LOGICAL LPIDOP

      INTEGER, pointer :: NSPT(:)

      real*8 COSB,COSB2,BOXOH,BOXOQ,RCUTP2,RCUTL2
      real*8 RFF, RFE, RFC, C1, RCRF2, RMULOC
      real*8 PININV

      COMMON /PMLINT/
     $     LDOTRA,LMONO,LOCTO,LVAC,LDOVIR,L4D,LPIDOP,
     $     LPERTL

      COMMON /PMINT/NSPT

      COMMON /PMLREL/
     $     RFF, RFE, RFC, C1, RCRF2, 
     $     COSB,COSB2,BOXOH,BOXOQ,RCUTP2,RCUTL2,PININV,
     $     RMULOC

