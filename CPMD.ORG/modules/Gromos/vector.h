C hello emacs, this is -*- fortran -*-
C

C     INCLUDE forcesz.h AND coordsz.h BEFORE THIS ONE
C
      INTEGER MAXP, MFORCE, MAXANG,
     $        JRCNPC, JRCNJA, IVECT,
     $        MXGGRP, IEPAIR, MAXGP

      logical, pointer :: LBVIR(:), LLPERT(:), LLPERJ(:)
      logical LNEWP

      integer, pointer :: INAG(:), IFIRST(:), 
     $     ICET(:), NSPTT(:), ISCLJT(:),
     $     ISCCT(:), ISCLJJ(:), ISCCJ(:), INTROW(:), 
     $     JRCIND(:), IACT(:), IACBT(:), IACBJ(:),
     $     IACAJ(:), JLIST(:), NILIST(:),
     $     JRCPAC(:), JFIRST(:), JRCILJ(:),
     $     JRCBK1(:),JRCBK2(:), IPICJ(:), 
     $     IPICT(:), NJL(:), NSPJ(:),
     $     IOCVIR(:), IFORCE(:), JFORCE(:)

      real*8, pointer ::  XC(:), YC(:), ZC(:), WC(:),
     $     XRC(:), YRC(:), ZRC(:), DF4TH(:),
     $     XDST(:), YDST(:), ZDST(:), 
     $     XJ(:), YJ(:), ZJ(:), WJ(:),
     $     VXIJ(:), VYIJ(:), VZIJ(:), VWIJ(:),
     $     FLRX(:), FLRY(:), FLRZ(:), FLRW(:),
     $     FXTEMP(:), FYTEMP(:), FZTEMP(:), FWTEMP(:), 
     $     CGAT(:), CGBT(:), CGAJ(:), CGBJ(:),
     $     FORC6(:), FORC12(:), VC6A(:), VC12A(:), 
     $     VC6B(:), VC12B(:), FC12C6(:), AC12C6(:), 
     $     VXH(:), VYH(:), VZH(:), VWH(:),
     $     BC12C6(:), R2IJ3D(:),
     $     R2D(:), XRJ(:), YRJ(:), ZRJ(:), 
     $     VRIJ(:), VRINV(:), VRINV2(:), VRINV6(:),
     $     VEELJ(:),DVEELJ(:), 
     $     VERFJ(:),VERCJ(:),VLJJ(:),
     $     DFRF(:), VRIJ2(:), TVIR(:, :), 
     $     TEEL(:, :), TERF(:, :), PIFACT(:),
     $     TERC(:, :), TELJ(:, :), DF3(:),
     $     EGLELJ(:), EGLRFJ(:), EGLRCJ(:), BSCALE(:),
     $     EGLLJJ(:), AMULT(:), SCALE(:), BMULT(:), 
     $     E34ELJ(:), E34RFJ(:), E34RCJ(:), E34LJJ(:)


      COMMON /LVECT/ LBVIR, LLPERT, LLPERJ
      COMMON /LVECT_s/ LNEWP

      COMMON /KVECT/ INAG, IFIRST, 
     $     ICET, NSPTT, ISCLJT,
     $     ISCCT, ISCLJJ, ISCCJ, INTROW, 
     $     JRCIND, IACT, IACBT, IACBJ,
     $     IACAJ, JLIST, NILIST,
     $     JRCPAC, JFIRST, JRCILJ,
     $     JRCBK1,JRCBK2, IPICJ, 
     $     IPICT, NJL, NSPJ,
     $     IOCVIR, IFORCE, JFORCE

      COMMON /RVECT/ XC, YC, ZC, WC,
     $     XRC, YRC, ZRC, DF4TH,
     $     XDST, YDST, ZDST, 
     $     XJ, YJ, ZJ, WJ,
     $     VXIJ, VYIJ, VZIJ, VWIJ,
     $     FLRX, FLRY, FLRZ, FLRW,
     $     FXTEMP, FYTEMP, FZTEMP, FWTEMP, 
     $     CGAT, CGBT, CGAJ, CGBJ,
     $     FORC6, FORC12, VC6A, VC12A, 
     $     VC6B, VC12B, FC12C6, AC12C6, 
     $     VXH, VYH, VZH, VWH,
     $     BC12C6, R2IJ3D,
     $     R2D, XRJ, YRJ, ZRJ, 
     $     VRIJ, VRINV, VRINV2, VRINV6,
     $     VEELJ,DVEELJ, 
     $     VERFJ,VERCJ,VLJJ,
     $     DFRF, VRIJ2, TVIR, 
     $     TEEL, TERF, PIFACT,
     $     TERC, TELJ, DF3,
     $     EGLELJ, EGLRFJ, EGLRCJ, BSCALE,
     $     EGLLJJ, AMULT, SCALE, BMULT, 
     $     E34ELJ, E34RFJ, E34RCJ, E34LJJ


      INTEGER NRVECT
      PARAMETER (NRVECT=9)
      COMMON /VECTSZ/ MAXP, MFORCE, MAXANG, JRCNPC, JRCNJA,
     &     IVECT, MXGGRP, IEPAIR, MAXGP
