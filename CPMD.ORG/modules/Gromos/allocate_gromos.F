      subroutine allocate_gromos
   
      IMPLICIT NONE
      include 'toposz.h'
      include 'topoar.h'
      include 'coordsz.h'
      include 'coordar.h'
      include 'forcesz.h'
      include 'nbpml.h'
      include '../system.h'
      INCLUDE '../parac.inc'

      allocate (INC(MAXCAG))
      allocate (WMASS(MAXNRS))

      if(.not.mmnode) return
      allocate (IAC(MAXNRP), NSPT(MAXINB)) 
      allocate (IPERT(MAXNRP))
      allocate (MPAC(MAXATT,MAXATT))
      allocate (FFTYPE (MAXATT))
      allocate (TOPTIT(MAXLNS))
      allocate (C12(MXATT2),C6(MXATT2))
      allocate (CS12(MXATT2),CS6(MXATT2))
      allocate (PANM(MAXNRP))
      allocate (AANM(MAXAA2))
      allocate (CB(MAXNBT),B0(MAXNBT))
      allocate (IB(MAXBON),JB(MAXBON),ICB(MAXBON))
      allocate (IBH(MAXBNH),JBH(MAXBNH),ICBH(MAXBNH))
      allocate (CT(MAXTTY),T0(MAXTTY))
      allocate (ITH(MXQHEH),JTH(MXQHEH),KTH(MXQHEH),ICTH(MXQHEH))
      allocate (IT(MAXTHE),JT(MAXTHE),KT(MAXTHE),ICT(MAXTHE))
      allocate (CQ(MAXQTY),Q0(MAXQTY))
      allocate (IQ(MAXQHI),JQ(MAXQHI),KQ(MAXQHI),LQ(MAXQHI))
      allocate (ICQ(MAXQHI))
      allocate (IQH(MAXHIH),JQH(MAXHIH),KQH(MAXHIH),LQH(MAXHIH))
      allocate (ICQH(MAXHIH))
      allocate (CP(MAXPTY),PD(MAXPTY))
      allocate (NP(MAXPTY))
      allocate (IP(MAXPHI),JP(MAXPHI),KP(MAXPHI),LP(MAXPHI))
      allocate (ICP(MAXPHI))
      allocate (IPH(MXPHIH),JPH(MXPHIH),KPH(MXPHIH),LPH(MXPHIH))
      allocate (ICPH(MXPHIH))
      allocate (MRES(MAXNRP))
      allocate (INE(MAXNRP),KNE(MAXNRP),INE14(MAXNRP),KNE14(MAXNRP))
      allocate (JSNE(MAXAEX))
      allocate (JSNE14(MXEX14))
      allocate (CG(MAXNRP),WINV(MAXNRP),WMAS(MAXNRP))
      allocate (IACS(MAXNRS))
      allocate (CGS(MAXNRS),WINVS(MAXNRS))
      allocate (ANMS(MAXNRS))
      allocate (ICONS(MXCONS),JCONS(MXCONS))
      allocate (CONS(MXCONS))
      allocate (IPIA(MAXPIA),IPIC(MAXNRP))
      allocate (IPIB(MAXPIB), JPIB(MAXPIB), ICPIB(MAXPIB))
      allocate (CPI(MAXPIT), WMCL(MAXPIA))

      IPIC(:)=0  ! cmb-default

C     allocate coordar  

      allocate(   X(MAXXCO),F(MAXXCO),V(MAXXCO))
      allocate(   XR(MAXXCO),XCEN(MAXXCO),SX(MAXXCO),XSAVE(MAXXCO))
      allocate(   IAGRP(MAXNAT))
      allocate(   JRC(MAXXC))
      allocate(   XC(MAXXCO),CXC(MAXNAT))
      allocate(   IDR1(MAXNDR),JDR1(MAXNDR),KDR1(MAXNDR),LDR1(MAXNDR))
      allocate(   ICDR1(MAXNDR))
      allocate(   IDR2(MAXNDR),JDR2(MAXNDR),KDR2(MAXNDR),LDR2(MAXNDR))
      allocate(   ICDR2(MAXNDR))
      allocate(   R0(MAXNDR),W0(MAXNDR))
      allocate(   RIIAVE(MAXTDR))
      allocate(   IPLR(MAXNDL),JPLR(MAXNDL),KPLR(MAXNDL))
      allocate(   LPLR(MAXNDL),ICPLR(MAXNDL))
      allocate(   CPLR(MAXNDL),PDLR(MAXNDL))
      allocate(   ICOG(MAXCON),JCOG(MAXCON))
      allocate(   CONP(MAXCON),FCON(MAXCON))
      allocate(   IPJV(MAXJ), JPJV(MAXJ), KPJV(MAXJ), LPJV(MAXJ))
      allocate(   CPJV(MAXJ),PJR0(MAXJ),PSJR(MAXJ))
      allocate(   AJV(MAXJ),BJV(MAXJ),CJV(MAXJ))
      allocate(   COSQAV(MAXJ),COSIAV(MAXJ))
      allocate(   IPLE(MAXDLE),JPLE(MAXDLE),KPLE(MAXDLE),LPLE(MAXDLE))
      allocate(   ILEMEM(MXLECF,MLECFG),NLEVST(MLECFG))
      allocate(   NSP(MAXNSP+1))
      allocate(   C4D(MAX4DA))
      allocate(   GAM(MAXNAT))
      allocate(   CC1(MAXNAT),CC2(MAXNAT),CC3(MAXNAT),CC4(MAXNAT))
      allocate(   CC5(MAXNAT),CC6(MAXNAT),CC7(MAXNAT),CC8(MAXNAT))
      allocate(   CC9(MAXNAT),SWINK(MAXNAT),SWINKS(MAXNAT))

#ifdef __VECTOR 
      call allocate_gromos_vector
#endif

      CALL PRMEM('ALLOCATE_GROMOS')
      return
      end

#ifdef __VECTOR 
      subroutine allocate_gromos_vector

      IMPLICIT NONE

      include 'coordsz.h'
      include 'forcesz.h'
      include 'vector.h'

c stuff for the vector version

      ALLOCATE(LBVIR(MAXNAT), LLPERT(MAXNAT), LLPERJ(MAXP))

      ALLOCATE(INAG(MAXINB), IFIRST(MAXINB))
      ALLOCATE(ICET(MAXNAT), NSPTT(MAXNAT), ISCLJT(MAXNAT))
      ALLOCATE(ISCCT(MAXNAT),ISCLJJ(MAXP),ISCCJ(MAXP),INTROW(MAXINB)) 
      ALLOCATE(JRCIND(MAXINB),IACT(MAXNAT),IACBT(MAXNAT),IACBJ(MAXP))
      ALLOCATE(IACAJ(MAXP), JLIST(MAXNAT), NILIST(MAXP))
      ALLOCATE(JRCPAC(JRCNPC), JFIRST(MAXINB), JRCILJ(MAXNAT))
      ALLOCATE(JRCBK1(MAXP),JRCBK2(MAXP), IPICJ(MAXP)) 
      ALLOCATE(IPICT(MAXNAT), NJL(MAXINB), NSPJ(MAXP))
      ALLOCATE(IOCVIR(MAXNAT), IFORCE(MAXP), JFORCE(MAXGP))

      ALLOCATE(XC(MAXNAT), YC(MAXNAT), ZC(MAXNAT), WC(MAXNAT))
      ALLOCATE(XRC(MAXNAT), YRC(MAXNAT), ZRC(MAXNAT), DF4TH(MAXP))
      ALLOCATE(XDST(MAXINB), YDST(MAXINB), ZDST(MAXINB))
      ALLOCATE(XJ(MAXNAT), YJ(MAXNAT), ZJ(MAXNAT), WJ(MAXNAT))
      ALLOCATE(VXIJ(MAXINB), VYIJ(MAXINB), VZIJ(MAXINB), VWIJ(MAXINB))
      ALLOCATE(FLRX(MAXNAT), FLRY(MAXNAT), FLRZ(MAXNAT), FLRW(MAXNAT))
      ALLOCATE(FXTEMP(MAXP), FYTEMP(MAXP), FZTEMP(MAXP), FWTEMP(MAXP)) 
      ALLOCATE(CGAT(MAXNAT), CGBT(MAXNAT), CGAJ(MAXP), CGBJ(MAXP))
      ALLOCATE(FORC6(MFORCE), FORC12(MFORCE), VC6A(MAXP), VC12A(MAXP)) 
      ALLOCATE(VC6B(MAXP), VC12B(MAXP), FC12C6(MFORCE), AC12C6(MAXP)) 
      ALLOCATE(VXH(MAXP), VYH(MAXP), VZH(MAXP), VWH(MAXP))
      ALLOCATE(BC12C6(MAXP), R2IJ3D(MAXP))
      ALLOCATE(R2D(MAXINB), XRJ(MAXINB), YRJ(MAXINB), ZRJ(MAXINB)) 
      ALLOCATE(VRIJ(MAXP), VRINV(MAXP), VRINV2(MAXP), VRINV6(MAXP))
      ALLOCATE(VEELJ(MAXINB),DVEELJ(MAXINB)) 
      ALLOCATE(VERFJ(MAXP),VERCJ(MAXP),VLJJ(MAXP))
      ALLOCATE(DFRF(MAXP), VRIJ2(MAXINB), TVIR(MAXP, 3)) 
      ALLOCATE(TEEL(IEPAIR, IVECT), TERF(IEPAIR, IVECT), PIFACT(MAXP))
      ALLOCATE(TERC(IEPAIR, IVECT), TELJ(IEPAIR, IVECT), DF3(MAXP))
      ALLOCATE(EGLELJ(MAXP), EGLRFJ(MAXP), EGLRCJ(MAXP), BSCALE(MAXGP))
      ALLOCATE(EGLLJJ(MAXP), AMULT(MAXGP), SCALE(MAXGP), BMULT(MAXGP)) 
      ALLOCATE(E34ELJ(MAXP), E34RFJ(MAXP), E34RCJ(MAXP), E34LJJ(MAXP))

      return
      end

#endif

