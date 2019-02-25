C     ==================================================================
C      note: the virial is not yet properly in.
C     ----------------------------------------------------------------==
Cmb - Revised on 28 January 2013
C the default is to use the NAG routines in this file
#if defined(__WITH_OMP3)
#define __COLLAPSE3 collapse(3)
#else
#define __COLLAPSE3
#endif
#define NAG_FFT
C The NAG routines use two REAL*8 grids and need a packed grid.
C
C define SGI_FFT :
C to use ZFFT3D (SGI), compile with -lcomplib.sgimath.
C The performance is better when each of the dimensions 
C                (NGRDX,NGRDY,NGRDZ)
C are multiples of 2, 3, 4 and 5. The 3D FFT runs in parallel.
#ifdef SGI_FFT
C The SGI routines use a complex grid and dont need a packed grid.
#define COMPLEX_GRID
#define UNPACKED_GRID
#undef NAG_FFT
#endif
C
C define CRAY_FFT:
C to use CCFFT3D (CRAY)
#ifdef CRAY_FFT
C The CRAY routines use a complex grid and dont need a packed grid.
#define COMPLEX_GRID
#define UNPACKED_GRID
#undef NAG_FFT
#endif
C     ----------------------------------------------------------------==
      SUBROUTINE KPPPM(NATTOT,CGT,X,F,EEWK,ERCOR,ESELF,
     .                 ALPHA,NGRDX,NGRDY,NGRDZ,BOX)
C     ----------------------------------------------------------------==
C
C P.HUNENBERGER, UCSD, SEPT 97
C THIS ROUTINE WILL SUPPLY THE K-SPACE ATOMIC FORCES AND ENERGY USING 
C HOCKNEY & EASTWOOD PARTICLE-PARTICLE PARTICLE MESH (PPPM)
C ALGORITHM. THIS VERSION USES GAUSSIAN SHAPING CHARGES (WIDTH ALPHA), A TRIANGLE
C SHAPE CHARGE (TSC) CHARGE SAMPLING AND POTENTIAL INTERPOLATING FUNCTION, AND
C AN A LINEAR COMBINATION OF THE 1ST AND 2ND ORDER FINITE DIFFERENCE OPERATORS
C (WITH RESPECTIVE WEIGHTS PBETA AND 1-PBETA).
C
C IT IS ASSUMED THAT THE UNIT CELL IS RECTANGULAR AND THAT THE SYSTEM IS IN 3D.
C
C KPPPM USES SUBROUTINE GHATEV.
C
C IF THE ROUTINE IS COMPILED WITH THE MACRO "COMPLEX_GRID", A SINGLE COMPLEX GRID
C WILL BE USED INSTEAD OF TWO real*8 GRIDS FOR THE real*8 AND IMAGINARY PARTS.
C IF THE ROUTINE IS COMPILED WITH THE MACRO "UNPACKED_GRID", THE GRID WILL NOT
C BE PACKED ( I.E. NUMBERS WILL BE STORED IN GRID(MAXGRD,MAXGRD,MAXGRD) - 
C IN CONTRAST TO PACKED, GRID(NGRDX*NGRDY*NGRDZ) ).
C
C BY DEFAULT, THE NAG SUBROUTINE C06FJF IS USED FOR THE FFTS. THIS ROUTINE
C REQUIRES TWO PACKED real*8 GRIDS FOR THE real*8 AND IMAGINARY PARTS.
C
C NATTOT = TOTAL NUMBER OF ATOMS (CHARGES) IN THE SYSTEM
C CGT(1..NATTOT) = ATOMIC CHARGES - MULTIPLIED BY SQRT(1/(4 PI EPS0))
C X(1..3*NATTOT) = ATOMIC CARTESIAN COORDINATES, X- AND Y-AXES LIE ALONG 
C                  BOX EDGES (HEAD:SOLUTES, REAR:SOLVENT)
C F(1..3*NATTOT) = DELIVERED INCREASED BY THE ATOMIC MESH FORCES
C EEWK   = DELIVERED WITH THE RECIPROCAL SPACE MESH POTENTIAL ENERGY    
C ERCOR = CORRECTION TERM TO THE real*8 SPACE ENERGY
C ESELF = WIGNER SELF TERM
C         (ERCOR AND ESELF ARE ONLY CALCULATED ONCE IN THE FIRST 
C          CALL AND NOT UPDATED AFTERWARDS)
C ALPHA  = CONVERGENCE PARAMETER (WIDTH OF THE GAUSSIAN SHAPING FUNCTION)
C NGRDX  = NUMBER OF GRID POINTS ALONG THE X COORDINATE
C NGRDY  = NUMBER OF GRID POINTS ALONG THE Y COORDINATE
C NGRDZ  = NUMBER OF GRID POINTS ALONG THE Z COORDINATE
C BOX(1..3) = EDGES OF THE (RECTANGULAR) PERIODIC BOX
C     ----------------------------------------------------------------==
C PARAMETERS
C NUMBER OF ALIAS VECTORS TO USE FOR CALCULATING THE OPTIMAL 
C INFLUENCE FUNCTION (GHAT) 
      INTEGER NSUMX, NSUMY, NSUMZ
      PARAMETER (NSUMX = 3, NSUMY = 3, NSUMZ = 3)
C USE ANALITICAL FORM FOR THE WHAT SQUARED SUM ENTERING INTO GHAT 
      LOGICAL LWHSUM
      PARAMETER (LWHSUM=.TRUE.)
C MAXIMUM NUMBER OF GRID POINTS ALLOWED IN EACH DIMENSION - NGRDX,NGRDY
C AND NGRDZ WILL BE CHECKED AGAINST THIS VALUE AT FIRST CALL 
      INTEGER MAXGRD
C UNIT FOR READING THE GHAT FUNCTION (BETA<0)
      INTEGER IUNITG
      PARAMETER (IUNITG = 99)
C
      real*8 BETMIN,BETMAX,BETTOL
      PARAMETER (BETMIN=8.0D-1,BETMAX=3.0D0,BETTOL=1.0D-5)
C
C WEIGHTING BETWEEN FIRST AND SECOND ORDER DIFFERENTIAL OPERATORS
      real*8 PBETA
      PARAMETER (PBETA = 4.0D0/3.0D0)
C
      real*8 TOL
      PARAMETER (TOL=1.0D-6)
      INTEGER MAXIT
      PARAMETER (MAXIT=200) 
      real*8 XICUB
      PARAMETER (XICUB= -2.83729748D0)
C
C ARGS
      INTEGER NATTOT, NGRDX,NGRDY,NGRDZ
      real*8  F(3*NATTOT), X(3*NATTOT), BOX(3)
      real*8 EEWK, ERCOR, ESELF, ALPHA
      real*8 CGT(NATTOT)
C LOCAL
      INTEGER IIW(-1:1), JJW(-1:1), KKW(-1:1),
     $     IIZ(-3:3),JJZ(-3:3),KKZ(-3:3),
     $     IPX, IPY, IPZ, II1, JJ1, KK1, ICHG, ICHG3, LL ,
     $     IX,IY,IZ
C      INTEGER*4 NDIM, NGRD(3), NGRDXY, NGRDXP, NGRD3, NWRK
      INTEGER NDIM, NGRD(3), NGRDXY, NGRDXP, NGRD3, NWRK
#ifndef COMPLEX_GRID
C      INTEGER*4 IFAIL
      INTEGER IFAIL
#endif
#ifndef UNPACKED_GRID
C      INTEGER*4 ITMP,ITMP1,ITMP2,ITMP3
      INTEGER ITMP,ITMP1,ITMP2,ITMP3
#endif
      real*8 WX(-1:1),WY(-1:1),WZ(-1:1),
     $     HX, HY, HZ, VOL,
     $     RPX, RPY, RPZ,XP, YP, ZP, W1, W2, W3, RHS, FACT,
     $     CX1, CY1, CZ1, CX2, CY2, CZ2, Q, DFX, DFY, DFZ
      CHARACTER*5 PRGSTR
      DATA PRGSTR /'KPPPM'/
      CHARACTER*16 STR1,STR2,STR3
      real*8 ESELFL, ERCORL, EKCORL
      SAVE ESELFL, ERCORL, EKCORL
C SAVED
      real*8 PYE
      SAVE PYE
      LOGICAL LFIRST
      SAVE LFIRST
      DATA LFIRST /.TRUE./
C GRIDS
      real*8, pointer ::  GHAT(:,:,:)
      common/saveall1/GHAT
#ifdef COMPLEX_GRID
      real*8 grm_re,grm_im
      DOUBLE COMPLEX, pointer :: GRID(:,:,:)
      common/ss1/GRID
#else
      real*8, pointer ::  GRIDR(:)
      real*8, pointer ::  GRIDI(:)
      common/ss1/GRIDR,GRIDI
#endif
C
C WORK ARRAYS
#ifdef NAG_FFT
      real*8, pointer :: WRK(:)
      common/saveall2/WRK
#endif
#ifdef SGI_FFT
      complex*16, pointer ::   WRK(:)
      common/saveall2/WRK
#endif
#ifdef CRAY_FFT
      INTEGER ISYS(4)
      SAVE ISYS
C IF ISYS=0
      REAL*8, pointer ::  WRK(:)
      REAL*8, pointer :: TAB(:)
      common/saveall2/WRK,TAB
      logical checkfact
      external checkfact
#endif

      INTEGER I,II, JJ, KK, IIP, JJP, KKP
      real*8 TMP
      real*8 A2,A2H,A3,A1
      real*8 PI,PI2,PI2L1,PI2L2,PI2L3
      real*8 A2TMP
      INTEGER LMAX
      real*8 GAMHAT,ALPHAC,ALPHAK
      real*8 SUMQ,SUMQ2
      INTEGER KX,KY,KZ
      REAL*8 HX2,HY2,HZ2,K2
C
C------------------------------------------------------------------------------C
C***** BEGIN KPPPM *************************************************************
C------------------------------------------------------------------------------C
C
      MAXGRD=max(NGRDX,NGRDY)
      MAXGRD=max(MAXGRD,NGRDZ)
      MAXGRD=MAXGRD+1
      IF (LFIRST) THEN
C
C***** FIRST CALL SETUP
C
        LFIRST = .FALSE.
C     
!
        allocate(GHAT(MAXGRD,MAXGRD,MAXGRD)) 
#ifdef COMPLEX_GRID
        allocate(GRID(MAXGRD,MAXGRD,MAXGRD))
#else
        allocate(GRIDR(MAXGRD*MAXGRD*MAXGRD))
        allocate(GRIDI(MAXGRD*MAXGRD*MAXGRD))
#endif
#ifdef NAG_FFT
        allocate(WRK(MAXGRD*MAXGRD*MAXGRD + 100))
#endif
#ifdef NAG_FFT
        STR1 = 'NAG FFT ROUTINE'
#endif
#ifdef SGI_FFT
        allocate( WRK(MAXGRD+15 + 2*MAXGRD+15 + 2*MAXGRD+15) )
#endif
#ifdef SGI_FFT
        STR1 = 'SGI FFT ROUTINE'
#endif
#ifdef CRAY_FFT
        STR1 = 'CRAY FFT ROUTINE'
        isys(1)=3
        IF(checkfact(NGRDX)) THEN 
          isys(2)=0
        ELSE
          isys(2)=1
        ENDIF
        IF(checkfact(NGRDY)) THEN
          isys(3)=0
        ELSE
          isys(3)=1
        ENDIF
        IF(checkfact(NGRDZ)) THEN
          isys(4)=0
        ELSE
          isys(4)=1
        ENDIF
        if(isys(2).eq.0 .and. isys(3).eq.0  .and. isys(4).eq.0) then
          allocate (TAB(2*(NGRDX+NGRDY+NGRDZ)))
        else
          write(6,*)'NON FACTORIZABLE GRID FOR THE P3M FFT'
          write(6,*)'EFFICIENCY WILL BE AFFECTED'
          allocate (TAB(12*(NGRDX+NGRDY+NGRDZ)))
        endif
        allocate (WRK(2*(NGRDX*NGRDY*NGRDZ)))
       
#endif
C
#ifdef COMPLEX_GRID
        STR2 = 'COMPLEX GRID'
#else
        STR2 = 'TWO real*8 GRIDS'
#endif
C
#ifdef UNPACKED_GRID
        STR3 = 'UNPACKED GRID'
#else
        STR3 = 'PACKED GRID'
#endif
C
        PRINT *
        PRINT '(6A)', ' PPPM USING ',STR1,', ',STR2,', ',STR3
        PRINT '(1A)', ' AND A SPHERICAL HAT FUNCTION'
        PRINT *
C
        IF (NGRDX.GT.MAXGRD-1.OR.NGRDY.GT.MAXGRD-1
     $       .OR.NGRDZ.GT.MAXGRD-1) THEN
          PRINT *,PRGSTR,': TOO MANY GRID POINTS'
          PRINT *,'INCREASE MAXGRD IN ',PRGSTR,
     $         ' TO MAX GRID POINTS + 1'
          CALL STOPGM('LATSUM','FATAL ERROR')
        ENDIF
C
C
C***** THE INITIAL GHAT GRID WITH THE OPTIMAL INFLUENCE FUNCTION IS
C      CALCULATED FOR THE GIVEN BETA VALUE 
C
        CALL GHATEV(Q, GHAT, MAXGRD, ALPHA, PBETA, 
     $       NGRDX, NGRDY, NGRDZ, 
     $       NSUMX, NSUMY, NSUMZ, BOX, LWHSUM, .TRUE.)
C     
C**** NOW SYMMETRIZE
C     
#ifdef _vpp_
!OCL NOALIAS
!OCL NOVREC
#endif
      DO KK = 0, NGRDZ/2
        KKP = MOD(NGRDZ - KK, NGRDZ)
        DO JJ = 0, NGRDY/2
          JJP = MOD(NGRDY - JJ, NGRDY)
#ifdef __NEC
!CDIR NODEP(GHAT)
#endif
          DO II = 0, NGRDX/2
            IIP = MOD(NGRDX - II, NGRDX)
C     
C**   GHAT IS SYMMETRIC AND ONLY THE VALUES WITH KX,KY,KZ>=0 WERE STORED.
C     RESTORE FULL GRID
C     
            TMP = GHAT(II+1,JJ+1,KK+1)
            GHAT(IIP+1,JJ+1 ,KK+1 ) = TMP
            GHAT(II+1 ,JJP+1,KK+1 ) = TMP
            GHAT(IIP+1,JJP+1,KK+1 ) = TMP
            GHAT(II+1 ,JJ+1 ,KKP+1) = TMP
            GHAT(IIP+1,JJ+1 ,KKP+1) = TMP
            GHAT(II+1 ,JJP+1,KKP+1) = TMP
            GHAT(IIP+1,JJP+1,KKP+1) = TMP
C     
          ENDDO
        ENDDO
      ENDDO
C     
C**** END GHATLD
C     

        PRINT '(A,/,A,F8.3,A,E15.6)',
     $       ' GHAT FUNCTION HAS BEEN CALCULATED USING',' PBETA =  ',
     $       PBETA,' RESULTING IN Q = ',Q
        PRINT *

C
C CALCULATE THE SELF TERM - THIS TERM IS A CONSTANT AND 
C GENERATES NO FORCE.
C
C     THE CONSTANTS A1, A2 AND A3 ARE CALCULATED HERE:
C
C       A1 = -4pi/V INT{o->infty} dr r^2 eta(r)
C       A2 =  4pi/V SUM{l<>0} gammahat_p(k)/k^2
C       A3 =        LIM{r->0} [eta(r) - 1/r]
C       
C     WHERE eta(r) IS THE SHORT-RANGE REAL-SPACE FUNCTION AND
C     gammahat_p(k) IS THE FOURIER TRANSFORM OF THE SHAPING FUNCTION.
C     
C     THE QUANTITY
C
C       - A2/2 SUM q^2 + A1/2 [(SUM q)^2 - SUM q^2]
C
C     IS ADDED TO THE POTENTIAL ENERGY. THE A2 TERM CORRECTS FOR THE
C     FACT THAT THE K-SPACE ENERGY IS CALCULATED AS A SINGLE SUM (EWALD)
C     OR THROUGH PPPM - I.E. CONTAINS AN UNWANTED SELF TERM. THE A1 TERM
C     CORRECTS FOR THE FACT THAT eta(r) IS USED IN THE R-SPACE SUM,
C     INSTEAD OF eta(r)-<eta(r)>_V
C
C     THE EWALD SELF TERM IS CALCULATED 
C     AS (A1+A2+A3) SUM q^2. FOR A CUBIC UNIT CELL OF EDGE L,
C     (A1+A2+A3) = -2.837297/L
C
        VOL = BOX(1)*BOX(2)*BOX(3)
C
C GET A1 AND A3
C
        PI = 4.0D0*DATAN(1.0D0)
        A1 = -4.0D0*PI*ALPHA*ALPHA/(15.0D0*VOL)
        A3 = -2.0D0/ALPHA
C
C COMPUTE A2 = 4pi/V sum_{l<>0} gammahat(k)/k^2
C
        ALPHAC = 12.0D0/ALPHA**4
        PI2 = 2.0D0*PI
        PI2L1 = PI2/BOX(1)
        PI2L2 = PI2/BOX(2)
        PI2L3 = PI2/BOX(3)        
        IF(DABS(BOX(2)-BOX(1)).LT.1.0D-10 
     .    .AND.DABS(BOX(3)-BOX(1)).LT.1.0D-10
     .    .AND.1.EQ.1) THEN
C Use the numerical value for cubic box
          PRINT *,'USING NUMERICAL XI FOR CUBIC CELL'
          A2 = XICUB/BOX(1)-A1-A3
        ELSE
C Calculate explicitly by summation
          A2 = 0.0D0
          LMAX = 0
          
 599      CONTINUE
          LMAX = LMAX + 1
          A2TMP = 0.0D0
          
          DO KX = -LMAX, LMAX, 2*LMAX
            HX2 = (PI2L1*DBLE(KX))**2
            DO KY = -LMAX, LMAX
              HY2 = (PI2L2*DBLE(KY))**2
              DO KZ = -LMAX, LMAX
                K2 = HX2 + HY2 + (PI2L3*DBLE(KZ))**2
                ALPHAK = ALPHA*DSQRT(K2)
                GAMHAT = ALPHAC/(K2*K2*K2)
     .               * (2.0D0-2.0D0*DCOS(ALPHAK)
     .               - ALPHAK*DSIN(ALPHAK))
                A2TMP = A2TMP + GAMHAT
              ENDDO
            ENDDO
          ENDDO
C
          DO KX = -LMAX+1, LMAX-1
            HX2 = (PI2L1*DBLE(KX))**2
            DO KY = -LMAX, LMAX, 2*LMAX
              HY2 = (PI2L2*DBLE(KY))**2
              DO KZ = -LMAX, LMAX
                K2 = HX2 + HY2 + (PI2L3*DBLE(KZ))**2
                ALPHAK = ALPHA*DSQRT(K2)
                GAMHAT = ALPHAC/(K2*K2*K2)
     .               * (2.0D0-2.0D0*DCOS(ALPHAK)
     .               - ALPHAK*DSIN(ALPHAK))
                A2TMP = A2TMP + GAMHAT
              ENDDO
            ENDDO
          ENDDO
          
          DO KX = -LMAX+1, LMAX-1
            HX2 = (PI2L1*DBLE(KX))**2
            DO KY = -LMAX+1, LMAX-1
              HY2 = (PI2L2*DBLE(KY))**2
              DO KZ = -LMAX, LMAX, 2*LMAX
                K2 = HX2 + HY2 + (PI2L3*DBLE(KZ))**2
                ALPHAK = ALPHA*DSQRT(K2)
                GAMHAT = ALPHAC/(K2*K2*K2)
     .               * (2.0D0-2.0D0*DCOS(ALPHAK)
     .               - ALPHAK*DSIN(ALPHAK))
                A2TMP = A2TMP + GAMHAT
              ENDDO
            ENDDO
          ENDDO          
          A2 = A2 + A2TMP
          IF(DABS(A2TMP).GT.TOL.AND.LMAX.LT.MAXIT) GOTO 599
C
          IF (LMAX.EQ.MAXIT) THEN
            PRINT *,' LATSUM| WARNING: A2 EVALUATION DID NOT CONVERGE!!'
          ENDIF
          WRITE(6,'(X,A,I4,2X,A,F10.8)') 'A2 EVALUATED WITH LMAX = ',
     &         LMAX,' TOL = ',A2TMP
C
          A2 = A2*4.0D0*PI/VOL
        ENDIF
        A2H = 0.5D0*A2
C     
C Phil 
        SUMQ = 0.0D0
        SUMQ2 = 0.0D0
!$OMP parallel do private(I) reduction(+:SUMQ,SUMQ2)
#ifdef __SR11000
*poption parallel, tlocal(I)
*poption psum(SUMQ,SUMQ2)
#endif
        DO I = 1,NATTOT
          SUMQ = SUMQ + CGT(I)
          SUMQ2 = SUMQ2 + CGT(I)**2
        ENDDO
C
        ERCORL = +0.5D0*A1*(SUMQ*SUMQ-SUMQ2) 
        EKCORL = -A2H*SUMQ2
        ESELFL = +0.5D0*(A1+A2+A3)*SUMQ2
C
C report a few things...
C
        WRITE(6,'(3(A20,F14.8,/),A20,F14.8,A,/,3(A20,E14.6,/))')
     &       'CONSTANT A1       =',A1,
     &       'CONSTANT A2       =',A2,
     &       'CONSTANT A3       =',A3,
     &       'CONSTANT XI       =',BOX(1)*(A1+A2+A3),' = LX*(A1+A2+A3)',
     &       'EPSFAC*[SUM Q]**2 =',SUMQ*SUMQ,
     &       'EPSFAC*SUM Q**2   =',SUMQ2,
     &       'ESELFL            =',ESELFL
C
C
C***** GET READY FOR FFTS
C     
#ifdef SGI_FFT
      CALL ZFFT3DI (NGRDX,NGRDY,NGRDZ,WRK)
#endif
#ifdef CRAY_FFT
      CALL CCFFT3D(0,NGRDX,NGRDY,NGRDZ,1.0,DUM,1,1,DUM,1,1,TAB,WRK,ISYS)
#endif
C     
        PYE = 4.0D0*DATAN(1.0D0)
C
C***** DONE FIRST CALL SETUP
C
      ENDIF
      ERCOR = ERCORL
      ESELF = ESELFL
C     
C------------------------------------------------------------------------------C
C*****ASSIGN CHARGE DENSITY TO GRID POINTS ************************************
C------------------------------------------------------------------------------C
C     
C SAMPLE THE CHARGE DENSITY ONTO A 3D GRID, GRID (COMPLEX, IMAGINARY=0) OR 
C GRIDRE (REAL) USING A TRIANGLE SHAPE CHARGE (TSC)
C SINCE GROMOS CHARGES ARE MULTIPLIED BY 1/SQRT(4*PI*EPS0), THE CHARGE
C DENSITY WILL CONTAIN THIS MULTIPLICATIVE FACTOR.
C
C     
C*****A FEW DERIVED QUANTITIES
C     
      HX = BOX(1) / NGRDX
      HY = BOX(2) / NGRDY
      HZ = BOX(3) / NGRDZ
      VOL = BOX(1)*BOX(2)*BOX(3)
      NDIM = 3
      NGRD(1) = NGRDX
      NGRD(2) = NGRDY
      NGRD(3) = NGRDZ
      NWRK = MAXGRD*MAXGRD*MAXGRD + 100
      NGRDXY = NGRDX*NGRDY
      NGRDXP = NGRDXY + NGRDX
      NGRD3 = NGRDX*NGRDY*NGRDZ

C     
C*****ZERO CHARGE DENSITY GRID
C     
#ifdef COMPLEX_GRID
      CALL mm_ZAZZERO(GRID,NGRDX*NGRDY*NGRDZ)
cmb #ifdef __VECTOR
cmb !$OMP parallel do private(IZ,IY,IX)
cmb #ifdef _vpp_
cmb !OCL NOVREC
cmb #endif 
cmb      DO IZ = 1, NGRDX
cmb        DO IY = 1,NGRDY
cmb          DO IX = 1,NGRDZ
cmb            GRID(IX,IY,IZ) = (0.0D0,0.0D0)
cmb          ENDDO
cmb        ENDDO
cmb      ENDDO
cmb #else
cmb !$OMP parallel do private(IX)
cmb      DO IX = 1,NGRDX*NGRDY*NGRDZ
cmb        GRID(IX,1,1) = (0.0D0,0.0D0) 
cmb      ENDDO
cmb #endif
#else
!$OMP parallel do private(IZ,IY,IX,ITMP) __COLLAPSE3
#ifdef _vpp_
!OCL NOVREC
#endif 
      DO IZ = 1, NGRDX
        DO IY = 1,NGRDY
          DO IX = 1,NGRDZ
            ITMP = IX+(IY-1)*NGRDX+(IZ-1)*NGRDXY
            GRIDR(ITMP) = 0.0D0
            GRIDI(ITMP) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
#endif
C
C VOLUME OF A GRID CELL, THE FACTOR 1/8 COMES FROM THE WEIGHT FUNCTIONS - SEE BELOW
      FACT = 0.125D0/(HX*HY*HZ)
C
      DO 100 ICHG = 1, NATTOT
C
C**** LOOP OVER CHARGES
C
C GET THE CHARGE, SCALE IT
        RHS = FACT * CGT(ICHG)
C COORDINATES OF THE CHARGE IN MESH UNITS
        ICHG3 = 3*ICHG
        RPX = X(ICHG3-2)/HX
        RPY = X(ICHG3-1)/HY
        RPZ = X(ICHG3  )/HZ
C COORDINATES OF THE CLOSEST GRID POINT IN MESH UNITS
        IPX = NINT(RPX) 
        IPY = NINT(RPY) 
        IPZ = NINT(RPZ) 
C VECTOR FROM CLOSEST MESH POINT TO CHARGE IN MESH UNITS
        XP = RPX - DBLE(IPX)
        YP = RPY - DBLE(IPY)
        ZP = RPZ - DBLE(IPZ)
C THESE WEIGHTS EXCLUDE A FACTOR 1/2 WHICH IS ADDED IN "FACT" AS 1/8
        WX(-1) = (0.5d0-XP)**2
        WX(0)  =  1.5d0-2.0d0*XP*XP
        WX(1)  = (0.5d0+XP)**2
        WY(-1) = (0.5d0-YP)**2
        WY(0)  =  1.5d0-2.0d0*YP*YP
        WY(1)  = (0.5d0+YP)**2
        WZ(-1) = (0.5d0-ZP)**2
        WZ(0)  =  1.5d0-2.0d0*ZP*ZP
        WZ(1)  = (0.5d0+ZP)**2
C FIND THE PERIODIC IMAGE OF THE CORRESPONDING GRID POINTS IN THE CENTRAL CELL 
        IIW(-1) = MOD(IPX-1 + NGRDX, NGRDX) + 1
        IIW(0)  = MOD(IPX   + NGRDX, NGRDX) + 1
        IIW(1)  = MOD(IPX+1 + NGRDX, NGRDX) + 1
        JJW(-1) = MOD(IPY-1 + NGRDY, NGRDY) + 1
        JJW(0)  = MOD(IPY   + NGRDY, NGRDY) + 1
        JJW(1)  = MOD(IPY+1 + NGRDY, NGRDY) + 1
        KKW(-1) = MOD(IPZ-1 + NGRDZ, NGRDZ) + 1
        KKW(0)  = MOD(IPZ   + NGRDZ, NGRDZ) + 1
        KKW(1)  = MOD(IPZ+1 + NGRDZ, NGRDZ) + 1
C NOW DISTRIBUTE THE CHARGE DENSITY TO THESE POINTS
        DO II1 = -1, 1
          W1 = WX(II1)*RHS
          DO JJ1 = -1, 1
            W2 = W1*WY(JJ1)
            DO KK1 = -1, 1
              W3 = W2*WZ(KK1)
#ifdef COMPLEX_GRID
              GRID(IIW(II1),JJW(JJ1),KKW(KK1))
     $             = GRID(IIW(II1),JJW(JJ1),KKW(KK1)) + CMPLX(W3,0.0D0)
#else
              ITMP = IIW(II1)+(JJW(JJ1)-1)*NGRDX
     $             +(KKW(KK1)-1)*NGRDXY
              GRIDR(ITMP) = GRIDR(ITMP) + W3
#endif
            ENDDO
          ENDDO
        ENDDO
C
C**** DONE LOOP OVER CHARGES
C
 100  CONTINUE

C
C------------------------------------------------------------------------------C
C***** SOLVE FOR THE POTENTIAL BY FFT ******************************************
C------------------------------------------------------------------------------C
C
C FORWARD FFT THE CHARGE DENSITY ON THE GRID, MULTIPLY WITH THE 
C OPTIMIZED GREENS FUNCTION GHAT AND BACKWARDS FFT TO GET THE ELECTROSTATIC 
C POTENTIAL AT THE GRID POINTS. GET THE INTERACTION ENERGY EEWK
C (SHOULD BE MULTIPLIED BY 4*PI DUE TO THE DEFINITION OF GROMOS CHARGES)
C
C
C FORWARD FFT THE CHARGE DENSITY R-SPACE -> K-SPACE
C
#ifdef NAG_FFT
      IFAIL = 0
      CALL C06FJF(NDIM,NGRD,NGRD3,GRIDR,GRIDI,WRK,NWRK,IFAIL)
#endif
#ifdef SGI_FFT
      CALL ZFFT3D (-1,NGRDX,NGRDY,NGRDZ,GRID,MAXGRD,MAXGRD,WRK)
      CALL ZSCAL3D(NGRDX,NGRDY,NGRDZ,1.D0/DSQRT(DBLE(NGRD3)),
     $     GRID,MAXGRD,MAXGRD)
#endif
#ifdef CRAY_FFT
      CALL CCFFT3D(-1,NGRDX,NGRDY,NGRDZ,1.D0/DSQRT(DBLE(NGRD3)),
     $     GRID,MAXGRD,MAXGRD,GRID,MAXGRD,MAXGRD,TAB,WRK,ISYS)
#endif
C
C     MULTIPLY THE FT OF THE CHARGE DENSITY BY THE FT OF THE GREENS FUNCTION
      EEWK = 0.0D0
      DO IZ = 1, NGRDZ
        DO IY = 1, NGRDY
#ifdef _vpp_
!OCL NOALIAS 
#endif
          DO IX = 1,NGRDX
#ifdef COMPLEX_GRID
            grm_re=DBLE(GRID(IX,IY,IZ))  ! cmb
            grm_im=AIMAG(GRID(IX,IY,IZ)) ! cmb
            EEWK  =EEWK+GHAT(IX,IY,IZ)*
     .           (grm_re**2 + grm_im**2)
            GRID(IX,IY,IZ)=CMPLX(GHAT(IX,IY,IZ)*grm_re, 
     .                           GHAT(IX,IY,IZ)*grm_im)
#else 
            ITMP = IX+(IY-1)*NGRDX+(IZ-1)*NGRDXY
            EEWK = EEWK+GHAT(IX,IY,IZ) * 
     .           (GRIDR(ITMP)**2 + GRIDI(ITMP)**2)
            GRIDR(ITMP) =  GHAT(IX,IY,IZ)*GRIDR(ITMP)
C TO DO INVERSE TRANSFORM WITH NAG WE NEED TO FORM COMPLEX CONJUGATE OF 
C GHAT*GRID TO GET A BACKWARD FFT INSTEAD OF A FORWARD FFT.
            GRIDI(ITMP) = -GHAT(IX,IY,IZ)*GRIDI(ITMP)
#endif
          ENDDO
        ENDDO
      ENDDO

C CORRECT FOR GROMOS CHARGES 4PI and A2 term
      EEWK = 2.0D0*PYE*EEWK*HX*HY*HZ + EKCORL
C
C BACKWARD FFT THE POTENTIAL K-SPACE -> R-SPACE
C
#ifdef NAG_FFT
      IFAIL = 0
      CALL C06FJF(NDIM,NGRD,NGRD3,GRIDR,GRIDI,WRK,NWRK,IFAIL)
#endif
#ifdef SGI_FFT
      CALL ZFFT3D  (1,NGRDX,NGRDY,NGRDZ,GRID,MAXGRD,MAXGRD,WRK)
      CALL ZSCAL3D (NGRDX,NGRDY,NGRDZ,1.D0/DSQRT(DBLE(NGRD3)),
     $     GRID,MAXGRD,MAXGRD)
#endif
#ifdef CRAY_FFT
      CALL CCFFT3D(1,NGRDX,NGRDY,NGRDZ,1.D0/DSQRT(DBLE(NGRD3)),
     $     GRID,MAXGRD,MAXGRD,GRID,MAXGRD,MAXGRD,TAB,WRK,ISYS) 
#endif
C
C------------------------------------------------------------------------------C
C***** INTERPOLATE THE POTENTIAL TO GET THE MESH FORCES ************************
C------------------------------------------------------------------------------C
C
C INTERPOLATE THE ELECTRIC FIELD FROM THE POTENTIAL ON THE MESH POINTS
C BY APPLYING A CENTERED 4 POINT FINITE DIFFERENCE SCHEME TO THE POTENTIAL 
C VALUES AT THE MESH POINTS. THE FORCES ON THE PARTICLES ARE OBTAINED BY USING
C THE SAME SAMPLING FUNCTION WHICH WAS USED TO ASSIGNE THE CHARGES TO THE 
C GRID POINTS (TSC).
C (SHOULD BE MULTIPLIED BY 4*PI DUE TO THE DEFINITION OF GROMOS CHARGES)
C

C THE FACTOR 4PI COMES FROM GROMOS CHARGE DEFINITION
C THE FACTOR 1/8 COMES FROM THE WEIGHT FUNCTIONS - SEE BELOW
      FACT = PYE/2.0D0
C
      CX1 = 0.50d0*PBETA/HX
      CY1 = 0.50d0*PBETA/HY
      CZ1 = 0.50d0*PBETA/HZ
      CX2 = 0.25d0*(1.0D0-PBETA)/HX
      CY2 = 0.25d0*(1.0D0-PBETA)/HY
      CZ2 = 0.25d0*(1.0D0-PBETA)/HZ
C
      DO ICHG = 1, NATTOT
C
C**** LOOP OVER CHARGES
C
C GET THE CHARGE, SCALE IT
        RHS = CGT(ICHG) * FACT 
C COORDINATES OF THE CHARGE IN MESH UNITS
        ICHG3 = 3*ICHG
        IX = ICHG3-2
        IY = ICHG3-1
        IZ = ICHG3
        RPX = X(IX)/HX
        RPY = X(IY)/HY
        RPZ = X(IZ)/HZ
C COORDINATES OF THE CLOSEST GRID POINT IN MESH UNITS
        IPX = NINT(RPX) 
        IPY = NINT(RPY) 
        IPZ = NINT(RPZ) 
C VECTOR FROM CLOSEST MESH POINT TO CHARGE IN MESH UNITS
        XP = RPX - DBLE(IPX)
        YP = RPY - DBLE(IPY)
        ZP = RPZ - DBLE(IPZ)
C THESE WEIGHTS EXCLUDE A FACTOR 1/2 WHICH IS ADDED IN "FACT" AS 1/8
        WX(-1) = (0.5d0-XP)**2
        WX(0)  =  1.5d0-2.0d0*XP*XP
        WX(1)  = (0.5d0+XP)**2
        WY(-1) = (0.5d0-YP)**2
        WY(0)  =  1.5d0-2.0d0*YP*YP
        WY(1)  = (0.5d0+YP)**2
        WZ(-1) = (0.5d0-ZP)**2
        WZ(0)  =  1.5d0-2.0d0*ZP*ZP
        WZ(1)  = (0.5d0+ZP)**2
C
C FIND THE PERIODIC IMAGE OF THE CORRESPONDING GRID POINTS IN THE CENTRAL CELL 
        DO LL = -3, 3
#ifdef UNPACKED_GRID
          IIZ(LL) = MOD(IPX + LL + NGRDX, NGRDX) + 1
          JJZ(LL) = MOD(IPY + LL + NGRDY, NGRDY) + 1
          KKZ(LL) = MOD(IPZ + LL + NGRDZ, NGRDZ) + 1
#else
          IIZ(LL) = MOD(IPX + LL + NGRDX, NGRDX) + 1
          JJZ(LL) = MOD(IPY + LL + NGRDY, NGRDY) * NGRDX  + NGRDX
          KKZ(LL) = MOD(IPZ + LL + NGRDZ, NGRDZ) * NGRDXY + NGRDXY
#endif
        ENDDO
C NOW INTERPOLATE THE FORCES
        DO II1 = -1, 1
          W1 = RHS * WX(II1)
          DO JJ1 = -1, 1
            W2 = W1*WY(JJ1)
            DO KK1 = -1, 1
              W3 = W2*WZ(KK1)
#ifdef COMPLEX_GRID
              DFX = W3 * (
     $             CX1*( GRID(IIZ(II1-1),JJZ(JJ1),KKZ(KK1)) -
     $             GRID(IIZ(II1+1),JJZ(JJ1),KKZ(KK1)) )  +
     $             CX2*( GRID(IIZ(II1-2),JJZ(JJ1),KKZ(KK1)) -
     $             GRID(IIZ(II1+2),JJZ(JJ1),KKZ(KK1)) )  )
              DFY = W3 * (
     $             CY1*( GRID(IIZ(II1),JJZ(JJ1-1),KKZ(KK1)) -
     $             GRID(IIZ(II1),JJZ(JJ1+1),KKZ(KK1)) )  +
     $             CY2*( GRID(IIZ(II1),JJZ(JJ1-2),KKZ(KK1)) -
     $             GRID(IIZ(II1),JJZ(JJ1+2),KKZ(KK1)) )  )
              DFZ = W3 * (
     $             CZ1*( GRID(IIZ(II1),JJZ(JJ1),KKZ(KK1-1)) -
     $             GRID(IIZ(II1),JJZ(JJ1),KKZ(KK1+1)) )  +
     $             CZ2*( GRID(IIZ(II1),JJZ(JJ1),KKZ(KK1-2)) -
     $             GRID(IIZ(II1),JJZ(JJ1),KKZ(KK1+2)) )  )
#else
              ITMP1 = JJZ(JJ1)+KKZ(KK1)-NGRDXP
              ITMP2 = IIZ(II1)+KKZ(KK1)-NGRDXP
              ITMP3 = IIZ(II1)+JJZ(JJ1)-NGRDXP
              DFX = W3 * (
     $             CX1*( GRIDR(IIZ(II1-1)+ITMP1) -
     $                   GRIDR(IIZ(II1+1)+ITMP1) )  +
     $             CX2*( GRIDR(IIZ(II1-2)+ITMP1) -
     $                   GRIDR(IIZ(II1+2)+ITMP1) )  )
              DFY = W3 * (
     $             CY1*( GRIDR(JJZ(JJ1-1)+ITMP2) -
     $                   GRIDR(JJZ(JJ1+1)+ITMP2) )  +
     $             CY2*( GRIDR(JJZ(JJ1-2)+ITMP2) -
     $                   GRIDR(JJZ(JJ1+2)+ITMP2) )  )
              DFZ = W3 * (
     $             CZ1*( GRIDR(KKZ(KK1-1)+ITMP3) -
     $                   GRIDR(KKZ(KK1+1)+ITMP3) )  +
     $             CZ2*( GRIDR(KKZ(KK1-2)+ITMP3) -
     $                   GRIDR(KKZ(KK1+2)+ITMP3) )  )
#endif
              F(IX)   = F(IX) + DFX
              F(IY)   = F(IY) + DFY
              F(IZ)   = F(IZ) + DFZ
C     
            ENDDO
          ENDDO
        ENDDO
C     
      ENDDO
C
C**** DONE LOOP OVER CHARGES
C
C
C------------------------------------------------------------------------------C
C***** DONE KPPPM **************************************************************
C------------------------------------------------------------------------------C
      RETURN
      END
C     ==================================================================
      SUBROUTINE GHATEV(Q, GRID, MAXGRD, ALPHA, PBETA, 
     $     NGRDX, NGRDY, NGRDZ, 
     $     NSUMX, NSUMY, NSUMZ, BOX, LWHSUM, LGRID)
C     ----------------------------------------------------------------==
CCC   SUBROUTINE GETQ(Q, PBETA, LWHSUM, LGRID)
C
C P.HUNENBERGER, UCSD, SEPT 97
C
C CALCULATES THE OPTIMAL GHAT FUNCTION FOR K-VECTORS WITH POSITIVE OR 
C NULL COMPONENTS (LGRID=.TRUE.) AND THE ERROR BETWEEN MESH AND TRUE FIELD (Q) 
C AT A GIVEN VALUE OF ALPHA AND PBETA AND FOR A GIVEN MESH SPANNING A GIVEN
C BOX VOLUME - MAY HAVE TO BE UPDATED PERIODICALLY IF THE BOX VOLUME 
C CHANGES.
C     
C Q      = DELIVERED WITH THE Q OPTIMIZATION MEASURE
C GRID   = DELIVERED WITH THE GHAT FUNCTION FOR K-VECTORS WITH POSITIVE 
C          OR NULL COMPONENTS (LGRID=.TRUE.)
C MAXGRD = DIMENSIONS OF THE GRID ARRAY ALONG THE THREE INDICES
C ALPHA  = CONVERGENCE PARAMETER (WIDTH OF THE GAUSSIAN SHAPING FUNCTION)
C PBETA  = WEIGHTING BETWEEN THE 1ST AND 2ND ORDER FINITE DIFFERENCE OPERATOR
C NGRDX  = NUMBER OF GRID POINTS ALONG THE X COORDINATE
C NGRDY  = NUMBER OF GRID POINTS ALONG THE Y COORDINATE
C NGRDZ  = NUMBER OF GRID POINTS ALONG THE Z COORDINATE
C NSUMX  = NUMBER OF X-ALIAS VECTORS TO USE
C NSUMY  = NUMBER OF Y-ALIAS VECTORS TO USE
C NSUMZ  = NUMBER OF Z-ALIAS VECTORS TO USE
C BOX(1..3) = EDGES OF THE (RECTANGULAR) PERIODIC BOX
C LWHSUM = .TRUE. : USE ANALITICAL FORM FOR THE WHAT SQUARED SUM 
C          .FALSE.: DO THE SUM BY HAND
C LGRID  = .TRUE. : RETURN THE OPTIMAL GREENS FUNCTION IN GRID
C     ----------------------------------------------------------------==
C ARGS
      LOGICAL LWHSUM, LGRID
      INTEGER NGRDX,NGRDY,NGRDZ,NSUMX,NSUMY,NSUMZ,MAXGRD
      real*8 ALPHA,PBETA,Q,BOX(3),GRID(MAXGRD,MAXGRD,MAXGRD)
C     LOCAL
      INTEGER LX,LY,LZ,NALIAX,NALIAY,NALIAZ
      REAL*8 HX,HY,HZ,U1,U2,U3,K2,GAMHAT,KL(3),KN(3),
     .     DH(3), WHRH(3), RH(3),DENOM,NUMER,
     .     DH2,RHSUM2,RH2,WH,WH2,WHSUM2,PBETONE,
     .     SYMFAC,SINX,SINY,SINZ,PYE,DUPYE,VOL,
     .     kal_x,kal_y,kal_z,yuyu
C
C declarations
      real*8 ALPHAC,ALPHAK
      REAL*8 RH_x,RH_y,RH_z
C
C definitions
      ALPHAC = 12.0D0/ALPHA**4
C
C***** BEGIN GHATEV
C
      HX = BOX(1) / NGRDX
      HY = BOX(2) / NGRDY
      HZ = BOX(3) / NGRDZ
      VOL = BOX(1)*BOX(2)*BOX(3)
      PYE = 4.0D0*DATAN(1.0D0)
      DUPYE = 2.0D0*PYE
      PBETONE=1.0D0-PBETA
      kal_x=DUPYE/BOX(1)
      kal_y=DUPYE/BOX(2)
      kal_z=DUPYE/BOX(3)
      Q = 0.0D0
C
C**** LOOP OVER KL VECTORS IN POSITIVE OCTANT OF THE THE GRID
C
      DO 10, LX = 0, NGRDX / 2
        KL(1) = kal_x * LX
        DO 20, LY = 0, NGRDY / 2           
          KL(2) = kal_y * LY
          DO 30, LZ = 0, NGRDZ / 2
            KL(3) = kal_z * LZ
            IF ((LX + LY + LZ) .EQ. 0) THEN
C
C**** ORIGIN : SKIP AND SET GHAT GRID ELEMENT TO ZERO IF REQUIRED
C
              IF (LGRID) THEN
                GRID(1,1,1)=0.0D0
              ENDIF
            ELSE
C
C**** NOW, CALCULATE THE FOURIER COEFFICIENTS DHAT OF THE FINITE
C     DIFFERENCING VECTOR OPERATOR. 
C     DHAT IS real*8 (SEE EQ. 190 - DUE TO MULTIPLICATION BY I) AND ODD.
C
C     KL(3) = RECIPROCAL SPACE VECTOR FOR WHICH TO CALCULATE DHAT
C     DH(3) = CORRESPONDING VALUE OF DHAT
C     DH2   = SQUARE NORM OF DHAT
C
              DH(1) =  PBETA*DSIN(KL(1)*HX)/HX
     .             + 0.5d0*PBETONE*DSIN(2.0D0*KL(1)*HX)/HX
              DH(2) =  PBETA*DSIN(KL(2)*HY)/HY
     .             + 0.5d0*PBETONE*DSIN(2.0D0*KL(2)*HY)/HY
              DH(3) =  PBETA*DSIN(KL(3)*HZ)/HZ
     .             + 0.5d0*PBETONE*DSIN(2.0D0*KL(3)*HZ)/HZ
              DH2 = DH(1)**2+DH(2)**2+DH(3)**2
C
              WHRH(1) = 0.0D0
              WHRH(2) = 0.0D0
              WHRH(3) = 0.0D0
              RHSUM2 = 0.0D0
              WHSUM2 = 0.0D0
C              
C**** LOOP OVER THE ALIAS VECTORS N=(-NSUMX,-NSUMY,-NSUMZ)...
C     (NSUMX,NSUMY,NSUMZ) WITH  K_N = K + (2PI/H) N.
C     CALCULATE WHAT (REAL AND EVEN).
C     CALCULATE RHAT = SUM_KN OF KN/K2 GAMMAHAT(KN).
C     RHAT IS real*8 AND ODD.
C     
C     KN(3)   = RECIPROCAL SPACE ALIAS VECTOR FOR WHICH TO CALCULATE WHAT
C     K2      = SQUARE NORM OF KN
C     WH      = CORRESPONDING VALUE OF WHAT
C     WH2     = SQUARED VALUE OF WHAT
C     WHSUM2  = RUNNING SUM OF WH2
C     RH(3)   = CORRESPONDING VALUE OF RHAT
C     RH2     = SQUARE NORM OF RHAT
C     RHSUM2  = RUNNING SUM OF OF RH2
C     WHRH(3) = RUNNING SUM OF WH2*RH
C
*poption parallel, tlocal(KN), psum(WHSUM2,RHSUM2,WHRH) 
              DO NALIAX = -NSUMX, NSUMX              
                KN(1) = KL(1)+DUPYE*NALIAX/HX              
                IF (KN(1) .EQ. 0.0D0) THEN
                  U1=1.0D0
                ELSE
                  yuyu=0.5d0*KN(1)*HX
                  U1=DSIN(yuyu)/yuyu
                ENDIF
                DO NALIAY = -NSUMY, NSUMY                
                  KN(2) = KL(2)+DUPYE*NALIAY/HY                
                  IF (KN(2) .EQ. 0.0D0) THEN
                    U2=1.0D0
                  ELSE
                    yuyu=0.5d0*KN(2)*HY
                    U2=DSIN(yuyu)/yuyu
                  ENDIF
                  DO NALIAZ = -NSUMZ, NSUMZ                  
                    KN(3) = KL(3)+DUPYE*NALIAZ/HZ
                    K2 = KN(1)*KN(1)+KN(2)*KN(2)+KN(3)*KN(3)
C     
C**** CALCULATE WHAT
                    IF (KN(3) .EQ. 0.0D0) THEN
                      U3=1.0D0
                    ELSE
                      yuyu=0.5d0*KN(3)*HZ
                      U3=DSIN(yuyu)/yuyu
                    ENDIF
                    WH = (U1*U2*U3)**3
                    WH2 = WH*WH
                    WHSUM2 = WHSUM2+WH2
C     
C**** CALCULATE RHAT
C
C get GAMMAHAT
                    ALPHAK = ALPHA*DSQRT(K2)
                    GAMHAT = ALPHAC/(K2*K2*K2)
     .                       *(2.0D0-2.0D0*DCOS(ALPHAK)
     .                         - ALPHAK*DSIN(ALPHAK))
C
                    RH_x = KN(1)*GAMHAT
                    RH_y = KN(2)*GAMHAT
                    RH_z = KN(3)*GAMHAT
                    RH2 = RH_x*RH_x+RH_y*RH_y+RH_z*RH_z  
                    RHSUM2 = RHSUM2+RH2
                    WHRH(1) = WHRH(1)+WH2*RH_x
                    WHRH(2) = WHRH(2)+WH2*RH_y 
                    WHRH(3) = WHRH(3)+WH2*RH_z
C
                  ENDDO
                ENDDO
              ENDDO
C
C**** WHSUM2 CAN ALSO BE CALCULATED BY AN ANALYTICAL FORMULA
C
              IF (LWHSUM) THEN                
                SINX = (DSIN(0.5d0*KL(1)*HX))**2
                SINY = (DSIN(0.5d0*KL(2)*HY))**2
                SINZ = (DSIN(0.5d0*KL(3)*HZ))**2
                WHSUM2=(1.d0-SINX+2.0d0/1.5d+1*SINX*SINX)*
     .                 (1.d0-SINY+2.0d0/1.5d+1*SINY*SINY)*
     .                 (1.d0-SINZ+2.0d0/1.5d+1*SINZ*SINZ) 
              ENDIF
C              
C**** NOW COMBINE NUMERATOR AND DENOMINATOR
C     QOPT = 1/V * SUM [ RHSUM2**2 - (DH*WHRH)**2/(DH**2 WHSUM2**2)]
C          = 1/V * SUM [ RHSUM2**2 - NUMER**2/DENOM ]
C     GHATOPT = DH*WHRH/(DH**2 WHSUM2**2)
C             = NUMER/DENOM
C
              NUMER = DH(1)*WHRH(1)+DH(2)*WHRH(2) 
     .               +DH(3)*WHRH(3)
              DENOM = DH2*WHSUM2*WHSUM2
C              
              IF (DABS(DENOM).LE.1.0D-10) THEN
                IF (DABS(NUMER).GT.1.0D-10) THEN
                  PRINT *,'GHATEV: DENOM=0 BUT NUMER<>0'
                  CALL STOPGM('LATSUM','FATAL ERROR')
                ENDIF
                IF (LGRID) THEN
                  GRID(LX+1,LY+1,LZ+1)=0.0D0
                ENDIF
              ELSE
C
C**** CALCULATE THE SYMMETRY WEIGHTING
C
                SYMFAC = 8.0d+00
                IF ( (LX .EQ. 0) .OR. (LX .EQ. NGRDX/2) ) 
     .               SYMFAC = 0.5d0*SYMFAC
                IF ( (LY .EQ. 0) .OR. (LY .EQ. NGRDY/2) ) 
     .               SYMFAC = 0.5d0*SYMFAC
                IF ( (LZ .EQ. 0) .OR. (LZ .EQ. NGRDZ/2) ) 
     .               SYMFAC = 0.5d0*SYMFAC
C
                Q = Q + SYMFAC * (RHSUM2 - NUMER*NUMER / DENOM)
                IF (LGRID) THEN
                  GRID(LX+1, LY+1, LZ+1) = NUMER / DENOM
                ENDIF
              ENDIF
            ENDIF
 30       CONTINUE
 20     CONTINUE
 10   CONTINUE
C
      Q = 16.0D0*PYE*PYE*Q/(VOL*VOL)
C
C***** END GHATEV
C
      RETURN
      END
C     ==================================================================
      SUBROUTINE KWALD(NATTOT,NPM,NRP,NSM,NRAM,CG,CGS,X,F,
     .             EEWK,ALPHA,KXMAX,KYMAX,KZMAX,K2MAX,BOX)
C     ----------------------------------------------------------------==
C
C P.HUNENBERGER, UCSD, SEPT 97
C
C THIS ROUTINE WILL SUPPLY THE K-SPACE ATOMIC FORCES AND ENERGY
C USING THE EWALD SUM METHOD. THIS VERSION USES GAUSSIAN SHAPING 
C CHARGES (WIDTH ALPHA).
C
C IT IS ASSUMED THAT THE UNIT CELL IS RECTANGULAR AND THAT THE SYSTEM IS IN 3D.
C
C NATTOT = TOTAL NUMBER OF ATOMS (CHARGES) IN THE SYSTEM
C          (= NPM*NRP+NSM*NRAM)
C NPM    = NUMBER OF SOLUTE MOLECULES
C NRP    = NUMBER OF ATOMS PER SOLUTE MOLECULE
C NSM    = NUMBER OF SOLVENT MOLECULES
C NRAM   = NUMBER OF ATOMS PER SOLVENT MOLECULE
C CG(1..NRP) = CHARGES OF THE SOLUTE ATOMS
C CGS(1..NRAM) = CHARGES OF THE SOLVENT ATOMS
C X(1..3*NATTOT) = ATOMIC CARTESIAN COORDINATES, X- AND Y-AXES LIE ALONG 
C                  BOX EDGES (HEAD:SOLUTES, REAR:SOLVENT)
C F(1..3*NATTOT) = DELIVERED INCREASED BY THE ATOMIC MESH FORCES
C EEWK   = DELIVERED WITH THE RECIPROCAL SPACE MESH POTENTIAL ENERGY    
C ALPHA  = CONVERGENCE PARAMETER (WIDTH OF THE GAUSSIAN SHAPING FUNCTION)
C KXMAX  = MAXIMAL ABSOLUTE VALUE OF THE COMPONENT OF K-VECTORS INCLUDED ALONG X
C KYMAX  = MAXIMAL ABSOLUTE VALUE OF THE COMPONENT OF K-VECTORS INCLUDED ALONG Y
C KZMAX  = MAXIMAL ABSOLUTE VALUE OF THE COMPONENT OF K-VECTORS INCLUDED ALONG X
C K2MAX  = SQUARE NORM OF K-VECTORS INCLUDED (K-SPACE CUTOFF)
C BOX(1..3) = EDGES OF THE (RECTANGULAR) PERIODIC BOX
C
C TEMPORARY ARRAYS:
C COSKX,COXKY,COXKZ = real*8 PART OF EXP(2*I*PI*H+R)
C     SINKX,SINKY,SINKZ = IMAGINARY PART OF EXP(2*PI*I**H+R)
C     COSKXY, SINKXY = TEMPORARY STORGAE FOR real*8 AND IMAGINARY
C     PART OF EXP( 2*PI*I/L [ KX*RX + KY*RY])
C     COSKR, SINKR = STORAGE FOR real*8 AND IMAGINARY PART OF
C     EXP(2*PI*I/L [KX*RX + KY*RY + KZ*RZ])
C     ----------------------------------------------------------------==
C PARAMETERS
C MAXIMUM NUMBER OF ATOMS - USED FOR BUILDING THE CHARGE VECTOR OF THE 
C WHOLE SYSTEM - NATTOT WILL BE CHECKED AGAINST THIS VALUE AT FIRST CALL
      INTEGER MAXNAT
      PARAMETER (MAXNAT=5000)
C MAXIMUM ABSOLUTE VALUE OF K-VECTOR COMPONENTS INCLUDED - KXMAX,KYMAX AND
C KZMAX WILL BE CHECKED AGAINST THIS VALUE AT FIRST CALL
      INTEGER MAXK
      PARAMETER (MAXK=120)
C ARGS
      INTEGER NPM, NRP, NSM, NRAM, NATTOT,KXMAX,KYMAX,KZMAX,K2MAX
      real*8  F(3*NATTOT), X(3*NATTOT), BOX(3), 
     $     CG(NRP), CGS(NRAM)
      real*8 EEWK, ALPHA      
C     LOCAL
      INTEGER I, J, IX, IY, IZ, N,
     $     K, KX, KY, KZ, KYY, KZZ, KYMIN, KZMIN

      real*8 K2, SUMCOS, SUMSIN, DF, GAMHAT,
     $     HX, HY, HZ, HX2, HY2, HZ2,
     $     PI4V,PI2L1,PI2L2,PI2L3,VOL,H2MAX
      real*8 COSKX(MAXNAT,MAXK+1), COSKY(MAXNAT,MAXK+1), 
     $     COSKZ(MAXNAT,MAXK+1)
      real*8 SINKX(MAXNAT,MAXK+1), SINKY(MAXNAT,MAXK+1), 
     $     SINKZ(MAXNAT,MAXK+1)
      real*8 COSKXY(MAXNAT), SINKXY(MAXNAT)
      real*8 COSKR(MAXNAT),  SINKR(MAXNAT)
C
      CHARACTER*5 PRGSTR 
      DATA PRGSTR /'KWALD'/
C     
      real*8 CGT(MAXNAT)
      SAVE CGT
      real*8 PI,PI2
      SAVE PI,PI2
      INTEGER KCOMAX
      SAVE KCOMAX
C
      LOGICAL LFIRST
      SAVE LFIRST
      DATA LFIRST /.TRUE./
C declarations
      real*8 ALPHAC,ALPHAK
C definitions
      ALPHAC = 12.0D0/ALPHA**4
C
C------------------------------------------------------------------------------C
C***** BEGIN KWALD *************************************************************
C------------------------------------------------------------------------------C
C
      IF (LFIRST) THEN
C
C***** FIRST CALL SETUP
C
        LFIRST = .FALSE.
C
        PRINT *
        PRINT '(1A)', ' EWALD USING A SPHERICAL HAT FUNCTION'
        PRINT *
C
        IF (NATTOT.GT.MAXNAT) THEN
          PRINT *,PRGSTR,': TOO MANY ATOMS'
          PRINT *,'INCREASE MAXNAT IN ',PRGSTR,' TO AT LEAST ',NATTOT
          CALL STOPGM('KWALD','FATAL ERROR')
        ENDIF
        IF (KXMAX.GT.MAXK.OR.KYMAX.GT.MAXK.OR.KZMAX.GT.MAXK) THEN
          PRINT *,PRGSTR,': TOO MANY K-VECTOR COMPONENTS'
          PRINT *,'INCREASE MAXK IN ',PRGSTR
          CALL STOPGM('KWALD','FATAL ERROR')
        ENDIF
C
        PI    = 4.0d0*DATAN(1.0d0)
        PI2   = 2.0d0*PI
        KCOMAX = KXMAX
        IF (KYMAX.GT.KCOMAX) KCOMAX = KYMAX
        IF (KZMAX.GT.KCOMAX) KCOMAX = KZMAX
C     
C     MAKE CHARGE AND IAC VECTOR OF SYSTEM SIZE.
C     
        K = 0
        DO I = 1, NPM
          DO J = 1, NRP
            K = K + 1
            CGT(K) =  CG(J)
          ENDDO
        ENDDO        
!$OMP parallel do private(I,J) shared(K)
        DO I = 1, NRAM
          DO J = I, NSM*NRAM, NRAM
            CGT(K+J) = CGS(I)
          ENDDO
        ENDDO
C        
      ENDIF
C
      VOL   = BOX(1)*BOX(2)*BOX(3)
      PI4V   = 4.0D0*PI/VOL
      PI2L1 = PI2/BOX(1)
      PI2L2 = PI2/BOX(2)
      PI2L3 = PI2/BOX(3)
C     THE RECIPROCAL SPACE CUT-OFF IS TAKEN SUCH THAT FOR EVERY
C     COMPONENT AN EQUAL AMOUNT OF K-VECTORS IS TAKEN IF BOX-SIZES VARY.
C     THIS MEANS FOR THE SAME ACCURACY ON A LONG EDGE ONE NEEDS MORE KS.
C        
      H2MAX = (PI2/BOX(1) *DSQRT(DBLE(K2MAX)))**2
C
C     THIS LOOKS WEIRD, BUT GIVES NOT THE SAME AS (PI2/BOX(1))**2 * K2X
C
      EEWK = 0.0D0
C     
C*****FIRST WE FILL SIX ARRAYS COSKX(J,N), SINKX(J,N) [+ SIMILAR FOR Y,Z]
C     CONTAINING RE() AND IM() OF
C     
C     EXP(I (N-1) 2PI/L_X R_JX)
C     
C     THIS IS DONE BY RECURRENCE:
C     
C     COSKR(J,1) = 1    AND    SINKR(J,1) = 0   FOR N = 1,
C     COSKR(J,1) = COS(2PI/L_X R_J,X) AND SINKR(J,1) = SIN(2PI/L_X R_J,X)  FOR N = 2,
C     AND
C     COSKR(J,N) = COSKR(J,N-1)*COSKR(J,2) - SINKR(J,N-1)*SINKR(J,2)
C     COSKR(J,N) = COSKR(J,N-1)*SINKR(J,2) + SINKR(J,N-1)*COSKR(J,2) FOR N = 3,KCOMAX+1
C     
!$OMP parallel do private(I,IX,IY,IZ)
      DO I=1, NATTOT
cmb        IX=3*I-2
cmb        IY=3*I-1
        IZ=3*I
        IY=IZ-1
        IX=IY-1
        COSKX(I,1)=1.0d+00
        COSKX(I,2)=DCOS(PI2L1*X(IX))
        COSKY(I,1)=1.0d+00
        COSKY(I,2)=DCOS(PI2L2*X(IY))
        COSKZ(I,1)=1.0d+00      
        COSKZ(I,2)=DCOS(PI2L3*X(IZ))
        SINKX(I,1)=0.0d+00
        SINKX(I,2)=DSIN(PI2L1*X(IX))
        SINKY(I,1)=0.0d+00
        SINKY(I,2)=DSIN(PI2L2*X(IY))
        SINKZ(I,1)=0.0d+00
        SINKZ(I,2)=DSIN(PI2L3*X(IZ))
      ENDDO
C     
      DO N = 3, KCOMAX + 1
        DO I = 1, NATTOT
          COSKX(I, N) = COSKX(I, N-1) * COSKX(I, 2) - 
     $         SINKX(I, N-1) * SINKX(I, 2)
          COSKY(I, N) = COSKY(I, N-1) * COSKY(I, 2) - 
     $         SINKY(I, N-1) * SINKY(I, 2)
          COSKZ(I, N) = COSKZ(I, N-1) * COSKZ(I, 2) - 
     $         SINKZ(I, N-1) * SINKZ(I, 2)
          SINKX(I, N) = SINKX(I, N-1) * COSKX(I, 2) + 
     $         COSKX(I, N-1) * SINKX(I, 2)
          SINKY(I, N) = SINKY(I, N-1) * COSKY(I, 2) + 
     $         COSKY(I, N-1) * SINKY(I, 2)
          SINKZ(I, N) = SINKZ(I, N-1) * COSKZ(I, 2) + 
     $         COSKZ(I, N-1) * SINKZ(I, 2)          
        ENDDO
      ENDDO
C     
C*****DONE FILL THE SIX ARRAYS
C     
C*****NOW LOOP OVER K-VECTORS, USE CENTRAL SYMMETRY AND DOUBLE THE
C     FINAL RESULT
C     KX : 0 ... KXMAX
C     KY : -KYMAX ... KYMAX (KX<>0)
C     0 ... KYMAX      (KX=0)
C     KZ : -KZMAX ... KZMAX (KX<>0 OR KY<>0)
C     1 ... KZMAX      (KX=KY=0)
C     
C     FIRST LOOP:  KX=KY=0 -> KZMIN=1, KZMIN CHANGES TO -KZMAX AT THE END
C     SECOND LOOP: KX=0,KY<>0 -> KYMIN=0, KYMIN CHANGES TO -KYMAX AT THE END
      KYMIN=0
      KZMIN=1
C     
      DO 500 KX=0, KXMAX
        HX = PI2L1*DBLE(KX)
        HX2 = HX*HX
C     
        DO 400 KY=KYMIN , KYMAX
          HY = PI2L2*DBLE(KY)
          HY2 = HY*HY
C     
C     K-SPACE CUTOFF
          IF ((HX2+HY2).GT.H2MAX) GOTO 400
C     
C***  CALCULATE FOR ALL I, COSKXY(I,N) AND SINKXY(I,N) CONTAINING RE() AND IM() OF
C     
C     EXP(I N_X 2PI/L_X R_JX + I N_Y 2PI/L_Y R_JY )
C     
          IF ( KY .GE. 0 ) THEN
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
            DO I=1, NATTOT
              COSKXY(I) = COSKX(I, KX+1) * COSKY(I, KY+1) - 
     $             SINKX(I, KX+1) * SINKY(I, KY+1)
              SINKXY(I) = SINKX(I, KX+1) * COSKY(I, KY+1) + 
     $             COSKX(I, KX+1) * SINKY(I, KY+1)
            ENDDO
          ELSE
            KYY=IABS(KY)
!$OMP parallel do private(I)
#ifdef __SR11000
*poption parallel, tlocal(I)
#endif
            DO I=1, NATTOT
              COSKXY(I) = COSKX(I, KX+1) * COSKY(I, KYY+1) + 
     $             SINKX(I, KX+1) * SINKY(I, KYY+1)
              SINKXY(I) = SINKX(I, KX+1) * COSKY(I, KYY+1) - 
     $             COSKX(I, KX+1) * SINKY(I, KYY+1)
            ENDDO
          ENDIF
C     
          DO 300 KZ=KZMIN, KZMAX            
            HZ = PI2L3*DBLE(KZ)
            HZ2 = HZ*HZ
            K2 = HX2 + HY2 + HZ2
C     
C     K-SPACE CUTOFF
            IF (K2 .GT. H2MAX) GOTO 300
C     
C     K-SPACE PREFACTOR - THE 2.0 MULTIPLICATIVE FACTOR COMES FROM SYMMETRY

            ALPHAK = ALPHA*DSQRT(K2)
            GAMHAT = 2.0d0*PI4V*ALPHAC/(K2*K2*K2)
     .             * (2.0d0-2.0d0*DCOS(ALPHAK)
     .             - ALPHAK*DSIN(ALPHAK))
C     
            SUMCOS = 0.0d0
            SUMSIN = 0.0d0
C     
C***  CALCULATE FOR ALL I, COSKR(I) AND SINKR(I) CONTAINING RE() AND IM() OF
C     
C     EXP(I N_X 2PI/L_X R_JX + I N_Y 2PI/L_Y R_JY + I N_Z 2PI/L_Z R_JZ)
C     
C     AND SUMCOS, SUMSIN, THE SUM OF Q_J TIMES THIS FACTOR OVER ALL ATOMS.
C     
            IF (KZ .GE. 0 ) THEN
!$OMP parallel do private(I) reduction(+:SUMCOS,SUMSIN)
#ifdef __SR11000
*poption parallel, tlocal(I)
*poption psum(SUMCOS,SUMSIN)
#endif
              DO I=1, NATTOT
                COSKR(I) = COSKXY(I) * COSKZ(I, KZ+1) - 
     $               SINKXY(I) * SINKZ(I, KZ+1)
                SINKR(I) = SINKXY(I) * COSKZ(I, KZ+1) + 
     $               COSKXY(I) * SINKZ(I, KZ+1) 
                SUMCOS = SUMCOS + CGT(I)*COSKR(I)
                SUMSIN = SUMSIN + CGT(I)*SINKR(I)
              ENDDO
C     
            ELSE
              KZZ=IABS(KZ)
!$OMP parallel do private(I) reduction(+:SUMCOS,SUMSIN)
#ifdef __SR11000
*poption parallel, tlocal(I)
*poption psum(SUMCOS,SUMSIN)
#endif
              DO I=1, NATTOT
                COSKR(I) = COSKXY(I) * COSKZ(I, KZZ+1) + 
     $               SINKXY(I) * SINKZ(I, KZZ+1)
                SINKR(I) = SINKXY(I) * COSKZ(I, KZZ+1) - 
     $               COSKXY(I) * SINKZ(I, KZZ+1)
                SUMCOS = SUMCOS + CGT(I)*COSKR(I)
                SUMSIN = SUMSIN + CGT(I)*SINKR(I) 
              ENDDO
C     
            ENDIF
C     
C*****NOW CALCULATE AND SUM UP THE FORCES CORRESPONDING TO THIS K-VECTOR
C     
!$OMP parallel do private(I,IX,IY,IZ,DF)
            DO I=1, NATTOT
cmb              IX=3*I-2
cmb              IY=3*I-1
              IZ=3*I
              IY=IZ-1 
              IX=IY-1
              DF= CGT(I)*GAMHAT*(SINKR(I)*SUMCOS - COSKR(I)*SUMSIN)
              F(IX)=F(IX) + HX*DF
              F(IY)=F(IY) + HY*DF
              F(IZ)=F(IZ) + HZ*DF
            ENDDO
C     
C*****CALCULATE AND SUM UP ENERGY CORRESPONDING TO THIS K-VECTOR
C     
            EEWK=EEWK+0.5d0*GAMHAT*(SUMCOS*SUMCOS+SUMSIN*SUMSIN)
 300      CONTINUE
          KZMIN=-KZMAX
C     
 400    CONTINUE
        KYMIN=-KYMAX
C     
 500  CONTINUE
C     
C------------------------------------------------------------------------------C
C***** END KWALD ***************************************************************
C------------------------------------------------------------------------------C
      RETURN
      END
C     ==================================================================     
      SUBROUTINE KWCOR (NATTOT,NDIM,NPM,NRP,NRAM,NCAG,NRAGT,
     .              INE,JSNE,KNE,IAGRP,MAXNRE,NRELKP,CG,CGS,
     .              ALPHA,BOX,XCOORD,NUNRE2,ESLF,EXCL,F,scale)
C     ----------------------------------------------------------------==
C P.HUNENBERGER, UCSD, SEPT 97
C     
C     ALTHOUGH NDIM MAY NOT EQUAL 3, THE CORRECTION IS CALCULATED IN 3D ONLY.
C     NO OPTION FOR ATOMIC CUTOFF HERE, EXCLUDED ATOMS ARE ASSUMED TO ALWAYS
C     BE AT A SMALLER DISTANCE.
C     
C     EWALD CORRECTION FOR EXCLUDED ATOMS:
C     THE K-SPACE SUM INTRODUCES A CONTRIBUTION ERF(R)/R FOR ALL
C     EXCLUDED NEIGHBOURS (IN THE CENTRAL UNIT CELL), IRRESPECTIVE
C     WHETHER EXCLUDED NEIGHBOURS ARE OR ARE NOT IN THE PAIR LIST.
C     THIS CONTRIBUTION IS REMOVED HERE - THAT IS (ERFC(R)-1)/R IS ADDED.
C
C     EXCLUSION CORRECTION IS RETURNED IN EXCL
C     ALSO THE SELF TERM (INITIALIZED, AT FIRST CALL ONLY)
C
C
C------------------------------------------------------------------------------C
C***** BEGIN KWCOR *************************************************************
C------------------------------------------------------------------------------C
C     ARGS
      INTEGER NATTOT,NDIM,NPM,NRP,NRAM,NCAG,NRAGT,MAXNRE
      INTEGER INE(*),JSNE(*),KNE(*)
      INTEGER NRELKP(MAXNRE,MAXNRE)
      INTEGER NUNRE2
      real*8 XCOORD(NATTOT*NDIM),F(NATTOT*NDIM),BOX(NDIM)
      real*8 CG(*),CGS(*),ESLF(NUNRE2),EXCL(NUNRE2)
      CHARACTER IAGRP(NATTOT)
      real*8 ALPHA,scale
C PARAMS
      INTEGER NCALCD
      PARAMETER (NCALCD=3)
      real*8 TOL
      PARAMETER (TOL=1.0D-8)
      INTEGER MAXIT
      PARAMETER (MAXIT=200) 
      real*8 XICUB
      PARAMETER (XICUB= -2.83729748D0)
C LOCAL
      LOGICAL LFIRST
      INTEGER NRAGP,NRPT,
     $     NIPROT,NI,II,IILOC,IILOC3,JJ,JJLOC,JJLOC3,
     $     M,ILGRP,JLGRP,ITMP,
     $     NDXGRP
      real*8 RIJ2,RIJINV,RIJIN2,DSTTMP,QIQJA,CGILA,VELA4D,DFA4D,
     $     XH,BOXH(3),XIJDST(3)
C SAVE
      SAVE LFIRST,BOXH,NRAGP,NRPT
      DATA LFIRST /.TRUE./

      real*8 FACG0,FACG2,FACG3,FACGF1,FACGF2,ALPHA2
C
      FACG0  = -2.0d0/ALPHA
      FACG2  =  2.0d0/ALPHA**3
      FACG3  = -1.0d0/ALPHA**4
      FACGF1 = -2.0d0*FACG2
      FACGF2 = -3.0d0*FACG3
      ALPHA2 = ALPHA*ALPHA
c      ALPHAC = 12.0D0/ALPHA**4
C
C BEGIN
      IF (LFIRST) THEN
C
C***** FIRST CALL OPERATIONS
C
        LFIRST = .FALSE.
        BOXH(1) = 0.5d0*BOX(1)
        BOXH(2) = 0.5d0*BOX(2)
        BOXH(3) = 0.5d0*BOX(3)
C NRAGP: THE TOTAL NUMBER OF SOLUTE CHARGE GROUPS
        NRAGP = NPM*NCAG
C NRPT : THE TOTAL NUMBER OF SOLUTE ATOMS
        NRPT  = NPM*NRP
C
      ENDIF
C
C***** END FIRST CALL OPERATIONS
C      NOW CALCULATE THE EXCLUSION CORRECTION
C     
      CALL mm_AZZERO(EXCL,NUNRE2)
C
C FIRST SOLUTE-SOLUTE EXCLUSIONS
C     
      DO 7050 NIPROT = 1,NPM
        DO 7100 II = 1,NRP-1
          IILOC = (NIPROT-1)*NRP+II
          ILGRP = ICHAR(IAGRP(IILOC))
          CGILA = CG(II)
          IILOC3 = 3*IILOC-3
          DO 7150 ITMP = 1,INE(II)
            JJ = JSNE(KNE(II)+ITMP)
            JJLOC = (NIPROT-1)*NRP+JJ
            JLGRP = ICHAR(IAGRP(JJLOC))
            JJLOC3 = 3*JJLOC-3
C     
            DO M=1,NCALCD
              XIJDST(M) = XCOORD(IILOC3+M) - XCOORD(JJLOC3+M)
            ENDDO
            RIJ2 = 0.0d0
C     RECTANGULAR BOUNDARY ONLY, HERE
            DO M=1,NCALCD
              DSTTMP = XIJDST(M)
              IF (DSTTMP .GE. BOXH(M)) THEN
                DSTTMP = DSTTMP - BOX(M)
              ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
                DSTTMP = DSTTMP + BOX(M)
              ENDIF
              XIJDST(M) = DSTTMP
              RIJ2 = RIJ2 + DSTTMP*DSTTMP
            ENDDO
C     
            QIQJA  = scale*CGILA*CG(JJ)
            RIJIN2 = 1.0d0/RIJ2
            RIJINV = DSQRT(RIJIN2)
C     
            IF (RIJ2.LE.ALPHA2) THEN
C exclusion distance shorter than alpha, add eta(r)-1/r
              VELA4D = QIQJA
     .             *(RIJ2*(FACG3/RIJINV+FACG2)+FACG0)
              DFA4D  = QIQJA
     .             *(FACGF1+FACGF2*RIJ2*RIJINV)
            ELSE
C exclusion distance larger than alpha, add -1/r
              VELA4D = -QIQJA*RIJINV
              DFA4D  =  VELA4D*RIJIN2
            ENDIF
C     
            NDXGRP = NRELKP(ILGRP,JLGRP)
            EXCL(NDXGRP) = EXCL(NDXGRP)+VELA4D
C     
#ifdef __VECTOR
            XH = XIJDST(1)*DFA4D
            F(IILOC3+1) = F(IILOC3+1) + XH
            F(JJLOC3+1) = F(JJLOC3+1) - XH
            XH = XIJDST(2)*DFA4D
            F(IILOC3+2) = F(IILOC3+2) + XH
            F(JJLOC3+2) = F(JJLOC3+2) - XH
            XH = XIJDST(3)*DFA4D
            F(IILOC3+3) = F(IILOC3+3) + XH
            F(JJLOC3+3) = F(JJLOC3+3) - XH
#else
            DO M=1,3
              XH = XIJDST(M)*DFA4D
              F(IILOC3+M) = F(IILOC3+M) + XH
              F(JJLOC3+M) = F(JJLOC3+M) - XH
            ENDDO 
#endif
C     
 7150     CONTINUE
 7100   CONTINUE
 7050 CONTINUE
C     
C     THEN SOLVENT-SOLVENT EXCLUSIONS
C     
C ES: in order to avoid double counting of solvent solvent
C     interactions if one calculates scaled 1-4 interactions 
      if(nram.eq.0) return

      DO II = 1, NRAM - 1
        CGILA = CGS(II)
        DO JJ = II + 1, NRAM
          QIQJA  = scale * CGILA * CGS(JJ)
CCC   NRAGP+1,NRAGT
          DO NI = 1, NRAGT-NRAGP
            IILOC = NRPT + (NI - 1)*NRAM + II
            ILGRP = ICHAR(IAGRP(IILOC))
            JJLOC = NRPT + (NI - 1)*NRAM + JJ
            JLGRP = ICHAR(IAGRP(JJLOC))
            IILOC3 = 3*IILOC-3
            JJLOC3 = 3*JJLOC-3
C     
            DO M=1,NCALCD
              XIJDST(M) = XCOORD(IILOC3+M) - XCOORD(JJLOC3+M)
            ENDDO
C     RECTANGULAR BOUNDARY ONLY, HERE
            RIJ2 = 0.0D0
            DO M=1,NCALCD
              DSTTMP = XIJDST(M)
              IF (DSTTMP .GE. BOXH(M)) THEN
                DSTTMP = DSTTMP - BOX(M)
              ELSEIF (DSTTMP .LT. -BOXH(M)) THEN
                DSTTMP = DSTTMP + BOX(M)
              ENDIF
              XIJDST(M) = DSTTMP
              RIJ2 = RIJ2 + DSTTMP*DSTTMP
            ENDDO
C     
            RIJIN2 = 1.0d0/RIJ2
            RIJINV = DSQRT(RIJIN2)
C     
            IF (RIJ2.LE.ALPHA2) THEN
C exclusion distance shorter than alpha, add eta(r)-1/r
              VELA4D = QIQJA
     .             *((FACG3/RIJINV+FACG2)*RIJ2+FACG0)
              DFA4D  = QIQJA
     .             *(FACGF1+FACGF2*RIJ2*RIJINV)
            ELSE
C exclusion distance larger than alpha, add -1/r
              VELA4D = -QIQJA*RIJINV
              DFA4D  =  VELA4D*RIJIN2
            ENDIF
C     
            NDXGRP = NRELKP(ILGRP,JLGRP)
            EXCL(NDXGRP) = EXCL(NDXGRP)+VELA4D
C     
            DO M=1,NCALCD
              XH = XIJDST(M)*DFA4D
              F(IILOC3+M) = F(IILOC3+M) + XH
              F(JJLOC3+M) = F(JJLOC3+M) - XH
            ENDDO
C     
          ENDDO
        ENDDO
      ENDDO
C     
C------------------------------------------------------------------------------C
C***** DONE KWCOR *************************************************************
C------------------------------------------------------------------------------C
      RETURN
      END
C
C-----------------------------------------------------------------------------C
C-----------------------------------------------------------------------------C
C WHAT FOLLOWS ARE FFT ROUTINES FROM THE NAG LIBRARY (C06FFT.F)
C JUST CHANGED SO THAT INTEGERS ARE PASSED AFTER REALS (FOR SUNS)
C
C     ==================================================================
      SUBROUTINE C06FJF(NDIM,ND,N,X,Y,WORK,LWORK,IFAIL)
C     ----------------------------------------------------------------==
C     MARK 11 RELEASE. NAG COPYRIGHT 1984.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     MULTI-DIMENSIONAL DISCRETE FOURIER TRANSFORM OF A
C     MULTI-DIMENSIONAL SEQUENCE OF COMPLEX DATA VALUES
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FJF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LWORK, N, NDIM
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(LWORK), X(N), Y(N)
      INTEGER           ND(NDIM)
C     .. Local Scalars ..
      INTEGER           IFAIL1, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
c      EXTERNAL          P01ABF
C     .. External Subroutines ..
c      EXTERNAL          C06FFF
C     .. Executable Statements ..
      IF (NDIM.LT.1) GO TO 40
      DO 20 L = 1, NDIM
         IFAIL1 = 1
         CALL C06FFF(NDIM,L,ND,N,X,Y,WORK,LWORK,IFAIL1)
         IF (IFAIL1.NE.0) GO TO 60
   20 CONTINUE
      IFAIL = 0
      RETURN
   40 IFAIL1 = 1
   60 IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,P01REC)
      RETURN
      END
C     ==================================================================
      SUBROUTINE C06FFF(NDIM,L,ND,N,X,Y,WORK,LWORK,IFAIL)
C     ----------------------------------------------------------------==
C     MARK 11 RELEASE. NAG COPYRIGHT 1984.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     DISCRETE FOURIER TRANSFORM OF ONE VARIABLE IN A
C     MULTI-VARIABLE SEQUENCE OF COMPLEX DATA VALUES
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FFF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, L, LWORK, N, NDIM
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(LWORK), X(N), Y(N)
      INTEGER           ND(NDIM)
C     .. Local Scalars ..
      INTEGER           I, IFAIL1, NI, NK, NL
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
c      EXTERNAL          P01ABF
C     .. External Subroutines ..
c      EXTERNAL          C06FFZ
C     .. Executable Statements ..
      IF (NDIM.LT.1) GO TO 60
      IF (L.LT.1 .OR. L.GT.NDIM) GO TO 100
      NK = 1
      NI = 1
      DO 20 I = 1, NDIM
         IF (ND(I).LT.1) GO TO 120
         IF (I.GT.L) NK = NK*ND(I)
         IF (I.LT.L) NI = NI*ND(I)
   20 CONTINUE
      NL = ND(L)
      IF (NI*NL*NK.NE.N) GO TO 80
      IF (NL.EQ.1) GO TO 40
      IF (LWORK.LT.3*NL) GO TO 160
      CALL C06FFZ(X,Y,NI,NL,NK,WORK,WORK(NL+1),WORK(2*NL+1),IFAIL1)
      IF (IFAIL1.NE.0) GO TO 140
   40 IFAIL = 0
      RETURN
   60 IFAIL1 = 1
      GO TO 180
   80 IFAIL1 = 2
      GO TO 180
  100 IFAIL1 = 3
      GO TO 180
  120 IFAIL1 = 3 + 10*I
      GO TO 180
  140 IFAIL1 = IFAIL1 + 10*L
      GO TO 180
  160 IFAIL1 = 4 + 10*L
  180 IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,P01REC)
      RETURN
      END
C     ==================================================================
      SUBROUTINE C06FFZ(X,Y,NI,NJ,NK,W1,W2,W3,IFAIL)
C     ==--------------------------------------------------------------==
C     MARK 11 RELEASE. NAG COPYRIGHT 1984.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     DISCRETE FOURIER TRANSFORM OF THE 2ND VARIABLE IN A
C     3-DIMENSIONAL SEQUENCE OF COMPLEX DATA VALUES
C
C     .. Scalar Arguments ..
      INTEGER           IFAIL, NI, NJ, NK
C     .. Array Arguments ..
      DOUBLE PRECISION  W1(NJ), W2(NJ), W3(NJ), X(NI,NJ,NK), Y(NI,NJ,NK)
C     .. Local Scalars ..
      INTEGER           I, J, K
C     .. External Subroutines ..
c      EXTERNAL          C06FCF
      INTEGER           KI 
C     .. Executable Statements ..
#if defined(__SR11000) || defined(__ES) || defined(_vpp_)
!$OMP parallel do private(KI,K,I,J,W1,W2,W3)
      DO KI=1,NK*NI 
        K=(KI-1)/NI+1
        I=MOD(KI-1,NI)+1
        DO J=1,NJ 
          W1(J)=X(I,J,K) 
          W2(J)=Y(I,J,K) 
        ENDDO 
        IFAIL=1 
        CALL C06FCF(W1,W2,NJ,W3,IFAIL) 
        DO J=1,NJ 
          X(I,J,K)=W1(J) 
          Y(I,J,K)=W2(J) 
        ENDDO 
      ENDDO 
#else 
      DO K = 1, NK
        DO I = 1, NI
!$OMP parallel do private(J) shared(I,K)
          DO J = 1, NJ
            W1(J) = X(I,J,K)
            W2(J) = Y(I,J,K)
          ENDDO
          IFAIL = 1
          CALL C06FCF(W1,W2,NJ,W3,IFAIL)
!$OMP parallel do private(J) shared(I,K)
          DO J = 1, NJ
            X(I,J,K) = W1(J)
            Y(I,J,K) = W2(J)
          ENDDO
        ENDDO
      ENDDO
#endif 
      RETURN
      END
C     ==================================================================
      LOGICAL FUNCTION checkfact(N)
C     ==--------------------------------------------------------------==
      integer N,NP
      checkfact=.false.
      NP=N
10    continue
      if(mod(NP,2).EQ.0)then
        NP=NP/2
        goto 10
      endif
11    continue
      if(mod(NP,3).EQ.0)then
        NP=NP/3
        goto 11
      endif
12    continue
      if(mod(NP,5).EQ.0)then
        NP=NP/5
        goto 12
      endif
      if(NP.eq.1)checkfact=.true.
      return
      END
