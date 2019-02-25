C routines for reading PROMD control files
Cmb - Revised on 15 Sept. 2005
C     ==================================================================
      SUBROUTINE RDMD
C     ==--------------------------------------------------------------==
C     RDMD reads a PROMD control file and checks for errors.
C     The input file is opened using L<OPNFIL>.
C     RDMD writes to screen and terminates if an
C     error occurs.
C
C     Note that the correctness checking performed in RDMD is
C     limited to little more than checking ranges of the variables.
C     A consistency check between the topology and the control
C     file read here must be performed elsewhere, typically in
C     L<PROMD>.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
      INCLUDE 'fileio.h'
C define external
      INTEGER NDXBT
      EXTERNAL NDXBT
C local vars
      INTEGER IUNIT,NBLOCK,K
      LOGICAL LOK
      CHARACTER PRGSTR*(4)
C common
      LOGICAL LREDMD
      COMMON /MDLOC/LREDMD
C data
      DATA PRGSTR /'RDMD'/

C     AK: set some defaults. 
C     AK: no printing.
      NTPR=0
      NTPL=0
      NTPP=0
C     AK:

C begin
      CALL OPNFIL('MDCT','FORMATTED','OLD',IUNIT)
      IF (IUNIT .LT. 0)THEN
         PRINT *,PRGSTR,': failed to open PROMD control file!'
         PRINT *
         CALL STOPGM('RDMD','FATAL ERROR')
      ENDIF

      CALL INITLN(IUNIT)
!$OMP parallel do private(K)
      DO K=MINBT,MAXBT
         LGOT(K) = .FALSE.
      ENDDO

      CALL GETLN(IUNIT,LOK)
C     while OK do
 10   IF (LOK) THEN
         NBLOCK = NDXBT(FIOLIN)

         IF (NBLOCK .EQ. NBTERR) THEN
            CALL FLAGLN(PRGSTR)
            PRINT '(A,A10)','unknown mdblock type ',FIOLIN
            CALL STOPGM('RDMD','FATAL ERROR')
         ENDIF

         IF (LGOT(NBLOCK)) THEN
C we have got the same block twice
            CALL FLAGLN(PRGSTR)
            PRINT *,'Block ',BTNAME(NBLOCK), 'occurs twice!'
            CALL STOPGM('RDMD','FATAL ERROR')
         ENDIF
         LGOT(NBLOCK) = .TRUE.
C TITLE
         IF (NBLOCK .EQ. NTITLE) THEN
            CALL GTITLE(IUNIT)
C BOUND
         ELSEIF (NBLOCK .EQ. NBOUND) THEN
            CALL GBOUND(IUNIT,PRGSTR)
C TCOUP
         ELSEIF (NBLOCK .EQ. NTCOUP) THEN
            CALL GTCOUP(IUNIT,PRGSTR)
C NPCOUP
         ELSEIF (NBLOCK .EQ. NPCOUP) THEN
            IF (LGOT(NBOUND) .AND. LGOT(NTCOUP)) THEN
               CALL GPCOUP(IUNIT,PRGSTR)
            ELSE
               CALL FLAGLN(PRGSTR)
               PRINT 510
               IF (.NOT. LGOT(NBOUND)) THEN
                  PRINT 520,BTNAME(NBOUND),BTNAME(NPCOUP)
               ENDIF
               IF (.NOT. LGOT(NTCOUP)) THEN
                  PRINT 520,BTNAME(NTCOUP),BTNAME(NPCOUP)
               ENDIF
               CALL STOPGM('RDMD','FATAL ERROR')
            ENDIF
C NSP
         ELSEIF (NBLOCK .EQ. NNSP) THEN
            CALL GNSP(IUNIT,PRGSTR)
C COM
         ELSEIF (NBLOCK .EQ. NCOM) THEN
            CALL GCOM(IUNIT,PRGSTR)
C STEP
         ELSEIF (NBLOCK .EQ. NSTEP) THEN
            CALL GSTEP(IUNIT,PRGSTR)
C SHAKE
         ELSEIF (NBLOCK .EQ. NSHAKE) THEN
            CALL GSHAKE(IUNIT,PRGSTR)
C FORCE
         ELSEIF (NBLOCK .EQ. NFORCE) THEN
            CALL GFORCE(IUNIT,PRGSTR)
C PLIST
         ELSEIF (NBLOCK .EQ. NPLIST) THEN
            IF (LGOT(NBOUND)) THEN
               CALL GPLIST(IUNIT,PRGSTR)
            ELSE
               CALL FLAGLN(PRGSTR)
               PRINT 510
               PRINT 520,BTNAME(NBOUND),BTNAME(NPLIST)
               CALL STOPGM('RDMD','FATAL ERROR')
            ENDIF
C PRINT
         ELSEIF (NBLOCK .EQ. NPRINT) THEN
            CALL GPRINT(IUNIT,PRGSTR)
C WRITE
         ELSEIF (NBLOCK .EQ. NWRITE) THEN
            CALL GWRITE(IUNIT,PRGSTR)
C PREST
         ELSEIF (NBLOCK .EQ. NPREST) THEN
            CALL GPREST(IUNIT,PRGSTR)
C DREST
         ELSEIF (NBLOCK .EQ. NDREST) THEN
            CALL GDREST(IUNIT,PRGSTR)
C SYSTEM
         ELSEIF (NBLOCK .EQ. NSYST) THEN
            CALL GSYST(IUNIT,PRGSTR)
C PERT
         ELSEIF (NBLOCK .EQ. NPERT) THEN
            CALL GPERT(IUNIT,PRGSTR)
C DHRES
         ELSEIF (NBLOCK .EQ. NDHRES) THEN
            CALL GDHRES(IUNIT,PRGSTR)
C START
         ELSEIF (NBLOCK .EQ. NSTART) THEN
            CALL GSTART(IUNIT,PRGSTR)
C FOURD
         ELSEIF (NBLOCK .EQ. NFOURD) THEN
            IF (LGOT(NBOUND)) THEN
               CALL GFOURD(IUNIT,PRGSTR)
            ELSE
               CALL FLAGLN(PRGSTR)
               PRINT 510
               PRINT 520,BTNAME(NBOUND),BTNAME(NFOURD)
               CALL STOPGM('RDMD','FATAL ERROR')
            ENDIF
C LONGRANGE
         ELSEIF (NBLOCK .EQ. NLONGR) THEN
            CALL GLRAN(IUNIT,PRGSTR)
C PATH INTEGRAL
         ELSEIF (NBLOCK .EQ. NPIBLK) THEN
            CALL GPI(IUNIT,PRGSTR)
C J-VAL RESTRAINING
         ELSEIF (NBLOCK .EQ. NJVAL) THEN
            CALL GJVAL(IUNIT,PRGSTR)
C LOCAL ELEVATION
         ELSEIF (NBLOCK .EQ. NLEBLK) THEN
            CALL GLE(IUNIT,PRGSTR)
C MINIMISE
         ELSEIF (NBLOCK .EQ. NEMBLK) THEN
            IF (LGOT(NSDBLK)) THEN
               CALL FLAGLN(PRGSTR)
               PRINT 500
               CALL STOPGM('RDMD','FATAL ERROR')
            ELSE
               CALL GEM(IUNIT,PRGSTR)
            ENDIF
C STOCHASTIC
         ELSEIF (NBLOCK .EQ. NSDBLK) THEN
            IF (LGOT(NEMBLK)) THEN
               CALL FLAGLN(PRGSTR)
               PRINT 500
               CALL STOPGM('RDMD','FATAL ERROR')
            ELSE
               CALL GSD(IUNIT,PRGSTR)
            ENDIF
            INCLUDE "inc_rdmd1.h"
C
C add other block types here
C
C unknown block type
         ELSE
            CALL FLAGLN(PRGSTR)
            PRINT '(A18,A10)','unknown mdblock type',FIOLIN
            CALL STOPGM('RDMD','FATAL ERROR')
         ENDIF
C end of the mega case statement

         CALL GETLN(IUNIT,LOK)
C     end while loop
         GOTO 10
      ENDIF

C---------------------
C set the IDOPRO value
      IF (LGOT(NEMBLK)) THEN
         IDOPRO = IDOEM
      ELSEIF (LGOT(NSDBLK)) THEN
         IDOPRO = IDOSD
      ELSE
         IDOPRO = IDOMD
      ENDIF
C---------------------


C here we check we have read all the block types
C that we have to have read in, depending on whether we do
C MD/SD or EM.
      LOK = .TRUE.
      IF (IDOPRO .EQ. IDOEM) THEN
         DO 50 K=MINBT,MAXBT
            IF (.NOT. LGOT(K) .AND. LNEDEM(K)) THEN
               PRINT 505,PRGSTR,BTNAME(K)
               LOK = .FALSE.
            ENDIF
 50      CONTINUE
      ELSE
         DO 55 K=MINBT,MAXBT
            IF (.NOT. LGOT(K) .AND. LNEDMD(K)) THEN
               PRINT 505,PRGSTR,BTNAME(K)
               LOK = .FALSE.
            ENDIF
 55      CONTINUE
      ENDIF

      IF ( .NOT. LOK) THEN
         CALL STOPGM('RDMD','FATAL ERROR')
      ENDIF

      LREDMD = .TRUE.

 500  FORMAT (1X,'MINIMISE block may not occur with a STOCHASTIC block')
 505  FORMAT (1X,A4,': block type ',A,' missing!')
 510  FORMAT (1X,'Block dependence!')
 520  FORMAT (1X,'The block ',A,' must appear before the block ',A)
C end rdmd
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR GLBCHK
C     SUBROUTINE GLBCHK(PRGSTR)
C
C     Subroutine GLBCHK performs a 'global check' of the
C     input switches to program L<PROMD> and L<PROFEE>.
C     Consistency checks that are impossible to perform
C     in the individual blocks reading routines are
C     performed here.
C     This is called after all the control blocks have been read in
C     in the subroutine L<RDMD>.
C
C     If an inconsistency is detected, a message is printed to screen
C     and STOP is invoked.
COMMEND
      SUBROUTINE GLBCHK(PRGSTR)
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
      INCLUDE 'md.h'
      INCLUDE 'box.h'
      INCLUDE 'mdblock.h'
      INCLUDE 'formats.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'forcesz.h'
C args
      CHARACTER PRGSTR*(*)
C local vars
      LOGICAL LSHKBO,LSHKBH
      LOGICAL LCLCBO,LCLCBH
      LOGICAL LOK
C begin
C make sure I am not shaking AND calculating bonds (or neither)
C this is only a warning
      LCLCBH = LTF(NTFBNH)
      LCLCBO = LTF(NTFBND)

      LSHKBH = (NTC .EQ. NTCDOH .OR. NTC .EQ. NTCDOB)
      LSHKBO = (NTC .EQ. NTCDOB)


C check bonds with H atoms
      IF (LCLCBH .AND. LSHKBH) THEN
         PRINT 500
         PRINT 501
         PRINT 510,NTC
         PRINT 520,NTFBNH
      ENDIF
      IF (.NOT. LCLCBH .AND. .NOT. LSHKBH) THEN
         PRINT 500
         PRINT 501
         PRINT 511,NTC
         PRINT 521,NTFBNH
      ENDIF

C now bonds with no H atoms
      IF (LCLCBO .AND. LSHKBO) THEN
         PRINT 500
         PRINT 502
         PRINT 510,NTC
         PRINT 520,NTFBND
      ENDIF
      IF (.NOT. LCLCBO .AND. .NOT. LSHKBO) THEN
         PRINT 500
         PRINT 502
         PRINT 511,NTC
         PRINT 521,NTFBND
      ENDIF


C check I dont try to T couple non existent solvent
      IF (NTT(NTTSLV) .EQ. NTTONE .AND. NRAM*NSM .EQ. 0) THEN
         PRINT 500,PRGSTR
         PRINT *,'No solvent is present, but NTT(3) =',NTT(NTTSLV)
         PRINT *,'NRAM:',NRAM
         PRINT *,'NSM :',NSM
         PRINT *,'How can you couple nonexistent solvent '
         PRINT *,'to a separate temperature bath ?'
         CALL STOPGM('GLBCHK','FATAL ERROR')
      ENDIF

C check I dont try to T couple non existent solute
      IF  (NPM*NRP .EQ. 0 .AND. 
     $     (NTT(NTTPIR) .EQ. NTTONE .OR. NTT(NTTPCM) .EQ. NTTONE .OR.
     $     (ABS(NTTPIR) .EQ. NTTTWO .AND. ABS(NTTPCM) .EQ. NTTTWO))
     $    )THEN
         PRINT 500,PRGSTR
         PRINT *,'No solute is present, but'
         PRINT *,'NTT(1) =',NTT(NTTPIR)
         PRINT *,'NTT(2) =',NTT(NTTPCM)
         PRINT *,'NRP:',NRP
         PRINT *,'NPM:',NPM
         PRINT *,'How can you couple nonexistent solute '
         PRINT *,'to a separate temperature bath ?'
         CALL STOPGM('GLBCHK','FATAL ERROR')
      ENDIF

C check that for pressure coupling, the submolecule definitions
C blocks have been read in
      IF (NTP .NE. NTPOFF .AND.
     $     (.NOT. LGOT(NNSP))) THEN
         PRINT 500,PRGSTR
         PRINT 512,NTP
         PRINT *,'but submolecule definition for solute incomplete!'
         IF (.NOT. LGOT(NPCOUP)) THEN
            PRINT *,BTNAME(NPCOUP),' needed'
         ENDIF
         IF (.NOT. LGOT(NNSP)) THEN
            PRINT *,BTNAME(NNSP),' needed'
         ENDIF
         CALL STOPGM('GLBCHK','FATAL ERROR')
      ENDIF

C no anistropic P-coupling with truncated octahedron
      IF (NTP .EQ. NTPANI .AND. ABS(NTB) .LT. 0) THEN
         PRINT 500,PRGSTR
         PRINT *,'Anisotropic pressure scaling specified: NTP =',NTP
         PRINT *,'which is not possible with truncated octahedron'
         PRINT *,'periodic boundary conditions: NTB = ',NTB
         CALL STOPGM('GLBCHK','FATAL ERROR')
      ENDIF

C check for implementation restrictions for EM,SD,MD
      IF (IDOPRO .EQ. IDOEM) THEN
C check NTB
         IF (ABS(NTB) .EQ. NTBVIR) THEN
            PRINT 500,PRGSTR
            PRINT FMNII,'NTB',NTB
            PRINT *,'Energy minimisation requires ABS(NTB) <>',NTBVIR
            CALL STOPGM('GLBCHK','FATAL ERROR')
         ENDIF
      ELSEIF (IDOPRO .EQ. IDOSD) THEN
C we must be in vacuum
         IF (NTB .NE. NTBVAC) THEN
            PRINT 500,PRGSTR
            PRINT FMNII,'NTB',NTB
            PRINT *,'Stochastic dynamics requires vacuum',
     $           ' boundary conditions, i.e. NTB = ',NTBVAC
            CALL STOPGM('GLBCHK','FATAL ERROR')
         ENDIF
C we have restricted possibilities for the NTT switches
C     no coupling
         LOK =( NTT(NTTPIR) .EQ. NTTOFF .AND.
     $        NTT(NTTPCM) .EQ. NTTOFF .AND.
     $        NTT(NTTSLV) .EQ. NTTOFF)
         IF (.NOT. LOK) THEN
C     solute and solvent separately coupled
            LOK =( ABS(NTT(NTTPIR)) .EQ. NTTTWO .AND.
     $           ABS(NTT(NTTPCM)) .EQ. NTTTWO .AND.
     $           NTT(NTTSLV) .EQ. NTTONE)
         ENDIF
         IF (.NOT. LOK) THEN
C     solute and solvent jointly coupled
            LOK =( ABS(NTT(NTTPIR)) .EQ. NTTHRE .AND.
     $           ABS(NTT(NTTPCM)) .EQ. NTTHRE .AND.
     $           ABS(NTT(NTTSLV)) .EQ. NTTHRE)
         ENDIF
         IF (.NOT. LOK) THEN
            PRINT 500,PRGSTR
            PRINT *,'Stochastic dynamics allows only a reduced set',
     $           ' of NTT switch settings:'
            PRINT *,'EITHER no temperature coupling'
            PRINT 550,NTTPIR,NTTOFF
            PRINT 550,NTTPCM,NTTOFF
            PRINT 550,NTTSLV,NTTOFF
            PRINT *,'OR solute and solvent separately coupled'
            PRINT 551,NTTPIR,NTTTWO
            PRINT 551,NTTPCM,NTTTWO
            PRINT 550,NTTSLV,NTTONE
            PRINT *,'OR solute and solvent jointly coupled'
            PRINT 551,NTTPIR,NTTHRE
            PRINT 551,NTTPCM,NTTHRE
            PRINT 551,NTTSLV,NTTHRE
            CALL STOPGM('GLBCHK','FATAL ERROR')
         ENDIF

      ELSEIF (IDOPRO .EQ. IDOMD) THEN
C no special checks here for now...
      ENDIF

C I must run in 3D if NTG = NTGMU or NTGBOT
      IF ((NTG .EQ. NTGMU .OR. NTG .EQ. NTGBOT)
     $     .AND. NT4DIM .EQ. NT4OFF) THEN
         PRINT 500,PRGSTR
         PRINT *,'MU dependent perturbation specified: NTG =',NTG
         PRINT *,'But system running in 3D: NT4DIM =',NT4DIM
         CALL STOPGM('GLBCHK','FATAL ERROR')
      ENDIF

C I must run in 4D with path integral simulations
      IF (NT4DIM .NE. NT4OFF .AND. NTPI .NE. NTPIOF) THEN
         PRINT 500,PRGSTR
         PRINT *,'Path integral simulation specified: NTPI  =',NTPI
         PRINT *,'4D simulation specified           : NT4DIM=',NT4DIM
         PRINT *,'These are inconsistent switch combinations'
         CALL STOPGM('GLBCHK','FATAL ERROR')
      ENDIF

C we may only have one solute molecule with local elevation
      IF (NTLE .NE. NTLEOF .AND. NPM .GT. 1) THEN
         PRINT 500,PRGSTR
         PRINT *,'local elevation only possible with NPM = 1'
         CALL STOPGM('GLBCHK','FATAL ERROR')
      ENDIF


C NTWSE may only point to a potential energy
      IF (NTWSE .NE. 0) THEN
         IF ((NTWSE .LT. IPTOT) .OR.
     $        (NTWSE .GT. MXEWRT)) THEN
            PRINT *,PRGSTR,': implementation restriction!'
            PRINT *,'NTWSE must specify a potential energy'
            PRINT *,'NTWSE = ',NTWSE
            PRINT *
            CALL STOPGM('GLBCHK','FATAL ERROR')
         ENDIF
      ENDIF

#include "inc_rdmd8.h"

 500  FORMAT (1X,A,': inconsistent input!')
 501  FORMAT (1X,'WARNING: BONDS INVOLVING H ATOMS:')
 502  FORMAT (1X,'WARNING: BONDS NOT INVOLVING H ATOMS:')

 510  FORMAT (1X,'ARE SHAKEN    : NTC = ',I2)
 520  FORMAT (1X,'AND CALCULATED: NTF(',I1,') = 1',//)

 511  FORMAT (1X,'ARE NOT SHAKEN: NTC = ',I2)
 521  FORMAT (1X,'NOR CALCULATED: NTF(',I1,') = 0',//)

 512  FORMAT (1X,'Pressure coupling specified: NTP =',I3)

 550  FORMAT (5X,'NTT(',I1,')  = ',I1)
 551  FORMAT (1X,'ABS(NTT(',I1,')) = ',I1)
C  glbchk
C     ==--------------------------------------------------------------==
      RETURN
      END

C     Those block types that have a LNEDMD(K) = .FALSE. or
C     LNEDEM(K) = .FALSE. need not appear in the input file
C     in all cases.
C     For those block types, we define some control values here to
C     turn those options off.
C

COMMSUBR LMDOK
C     LOGICAL FUNCTION LMDOK()
C
C     LMDOK returns .TRUE. if the PROMD control file
C     has successfully been read in, .FALSE. otherwise.
COMMEND
      LOGICAL FUNCTION LMDOK()
      LOGICAL LREDMD
      COMMON /MDLOC/LREDMD

      LMDOK = LREDMD
      RETURN
      END

COMMSUBR NDXBT
C     INTEGER FUNCTION NDXBT(BTN)
C
C     NDXBT returns the token number of the string passed
C     in BTN. The string must represent a blockname for
C     the PROMD control input.
C     The values of the tokens are defined in 'mdblock.h'
COMMEND
      INTEGER FUNCTION NDXBT(BTN)
      CHARACTER *(*)BTN
      INCLUDE 'mdblock.h'
      INTEGER I
C begin
      I = MINBT
C look for BTN in table using while loop
 10   IF (I .LE. MAXBT) THEN
         IF (BTNAME(I) .NE. BTN) THEN
            I = I + 1
            GOTO 10
         ENDIF
      ENDIF
C     end of while loop
C
      IF (I .GT. MAXBT) THEN
         NDXBT = NBTERR
      ELSE
         NDXBT = I
      ENDIF
C end NDXBT
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GTITLE(IUNIT)
C     ==--------------------------------------------------------------==
C include
      INCLUDE 'coordsz.h'
      INCLUDE 'md.h'
C args
      INTEGER IUNIT
C local vars
      INTEGER NDUM
C begin
      CALL RDTIT(IUNIT,.TRUE.,NDUM,MDTLIN,MDTITL,MDTLNS)
C end gtitle
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PTITLE
C     ==--------------------------------------------------------------==
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C begin
      CALL PRTIT(BTNAME(NTITLE),MDTLNS,MDTITL)
C end ptitle
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GSYST(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(8)
      CHARACTER VSTR*(4)
C data
      DATA STR /'NPM,NSM'/
C begin
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NPM,PRGSTR,STR)
      CALL CHPINT(NSM,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)

      CALL GTEND(IUNIT,PRGSTR)

C check NPM
      IF (NPM .LT. 0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NPM'
         PRINT FMNII,VSTR,NPM
         PRINT FMGEI,VSTR
         CALL STOPGM('GSYST','FATAL ERROR')
      ENDIF

C check NSM
      IF (NSM .LT. 0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NSM'
         PRINT FMNII,VSTR,NSM
         PRINT FMGEI,VSTR
         CALL STOPGM('GSYST','FATAL ERROR')
      ENDIF
C end
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PSYST
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C begin
      PRINT *,BTNAME(NSYST)
      PRINT '(2A10)','NPM','NSM'
      PRINT '(2I10)',NPM,NSM
C end psyst
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GSTART(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'formats.h'
C local vars
#ifdef READ_COORD
      CHARACTER STR*(50)
#else
      CHARACTER STR*(36)
#endif
      CHARACTER VSTR*(6)
C data
#ifdef READ_COORD
      DATA STR /'NTX,INIT,NRDCOR,NSTCOR,IG,TEMPI,HEAT,NTXO,BOLTZ'/
#else
      DATA STR /'NTX,INIT,IG,TEMPI,HEAT,NTXO,BOLTZ'/
#endif
C begin
C NTX,NTCX,IG,TEMPI,HEAT,NTXO
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTX,PRGSTR,STR)
      CALL CHPINT(INIT,PRGSTR,STR)
#ifdef READ_COORD
      CALL CHPINT(NRDCOR,PRGSTR,STR)
      CALL CHPINT(NSTCOR,PRGSTR,STR)
#endif
      CALL CHPINT(IG,PRGSTR,STR)
      CALL CHPREL(TEMPI,PRGSTR,STR)
      CALL CHPREL(HEAT,PRGSTR,STR)
      CALL CHPINT(NTXO,PRGSTR,STR)
      CALL CHPREL(BOLTZ,PRGSTR,'BOLTZ')
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)
C check NTX
      IF (NTX .LT. NTXMIN .OR. NTX .GT. NTXMAX) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTX'
         PRINT FMNII,VSTR,NTX
         PRINT FMBTI,NTXMIN,VSTR,NTXMAX
         CALL STOPGM('GSTART','FATAL ERROR')
      ENDIF

C check INIT
      IF (INIT .LT. INITMI .OR. INIT .GT. INITMA) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'INIT'
         PRINT FMNII,VSTR,INIT
         PRINT FMBTI,INITMI,VSTR,INITMA
         CALL STOPGM('GSTART','FATAL ERROR')
      ENDIF

C check NTXO
      IF (NTXO .LT. NTXOMI .OR. NTXO .GT. NTXOMA) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTXO'
         PRINT FMNII,VSTR,NTXO
         PRINT FMBTI,NTXOMI,VSTR,NTXOMA
         CALL STOPGM('GSTART','FATAL ERROR')
      ENDIF

C check BOLTZ
      IF (BOLTZ .LE. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'BOLTZ'
         PRINT FMNIR,VSTR,BOLTZ
         PRINT FMGTR,VSTR
         CALL STOPGM('GSTART','FATAL ERROR')
      ENDIF
C end gstart
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PSTART
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C local params
C Boltzman constant for KJ/Mol
      real*8 BOLKJM
      PARAMETER (BOLKJM = 8.31441D-3)
C Boltzman constant for KCal/Mol
      real*8 BOLKCM
      PARAMETER (BOLKCM = BOLKJM/4.184D0)
C begin
      PRINT *,BTNAME(NSTART)
#ifdef READ_COORD
      PRINT '(7A10,A5,A12)',
     $     'NTX','INIT','NRDCOR','NSTCOR','IG',
     $     'TEMPI','HEAT','NTXO','BOLTZ'
      PRINT '(5I10,2F10.3,I5,E12.5)',NTX,INIT,NRDCOR,NSTCOR,IG,
     $     TEMPI,HEAT,NTXO,BOLTZ
#else
      PRINT '(5A10,A5,A12)',
     $     'NTX','INIT','IG','TEMPI','HEAT','NTXO','BOLTZ'
      PRINT '(3I10,2F10.3,I5,E12.5)',NTX,INIT,IG,
     $     TEMPI,HEAT,NTXO,BOLTZ
#endif
C
Cwrite out what we do, depending on NTX
C coords
      PRINT *,'   INITIAL COORDINATES ARE READ FROM FILE'
C
C generating/reading in of vel
      IF (TEMPI .GT. EPS) THEN
         PRINT 510
      ELSE
         IF (NTX .GE. NTXV) THEN
            PRINT *,'   INITIAL VELOCITIES ARE READ FROM FILE'
C     scaling of vel
            IF (HEAT .GT. EPS) THEN
               PRINT *,'   ALL (SOLUTE AND SOLVENT) VELOCITIES ARE ',
     $              'MULTIPLIED BY ''HEAT'''
            ELSEIF(HEAT .LT. -EPS) THEN
               PRINT *,'   ALL SOLUTE VELOCITIES ARE MULTIPLIED BY ',
     $              'ABS(HEAT). SOLVENT VELOCITIES ARE NOT CHANGED'
            ENDIF
         ELSE
            PRINT *,'   INITIAL VELOCITIES ARE SET TO ZERO'
         ENDIF
      ENDIF
C reading in of SX, IG
      IF (NTX .EQ. NTSX) THEN
         PRINT *,'   SX AND IG READ FROM FILE IF SD SPECIFIED'
      ENDIF
C INIT
      IF (INIT .EQ. INSHVX) THEN
         PRINT *,'   INITIAL COORDINATES ARE SHAKEN'
      ENDIF
      IF (INIT .GE. INSHVX .AND. INIT .LE. INSHV) THEN
         PRINT *,'   INITIAL VELOCITIES ARE SHAKEN'
      ENDIF
      IF (INIT .GE. INSHVX .AND. INIT .LE. INNOSH .AND. LTCM) THEN
         PRINT *,'   INITIAL CENTRE OF MASS MOTION IS REMOVED'
      ENDIF
      IF (INIT .EQ. INITCO) THEN
         PRINT *,'   NO SHAKE, NO C.O.M. MOTION REMOVAL ',
     $     'FROM INITIAL CONFIGURATION'
      ENDIF

#include "inc_rdmd7.h"

C NTXO
      IF (NTXO .EQ. NTXOUF)THEN
         PRINT *,'   UNFORMATTED OUTPUT OF:'
      ELSE
         PRINT *,'   FORMATTED OUTPUT OF:'
      ENDIF
      PRINT *,'      FINAL COORDINATES AND VELOCITIES'
      IF (NTP .NE. NTPOFF) THEN
         PRINT *,'      BOX DIMENSIONS'
         IF (NTR .NE. NTROFF) THEN
            PRINT *,
     $           '      SCALED REFERENCE POSITIONS FOR POS RESTRAINING'
         ENDIF
      ENDIF
      IF (NTG .NE. NTGOFF)THEN
         PRINT *,'      FINAL VALUES FOR RLAM,RMU,DSVDL AND DSKDL'
      ENDIF

C BOLTZ
      IF (DABS(BOLTZ-BOLKJM) .LE. EPS) THEN
         PRINT *,'   BOLTZ VALUE IS IN KJ/MOL'
      ELSEIF (DABS(BOLTZ-BOLKCM) .LE. EPS) THEN
         PRINT *,'   BOLTZ VALUE IS IN KCAL/MOL'
      ELSE
         PRINT *,'   BOLTZ is in unknown units'
         PRINT 500,BOLKJM,'KJ/MOL'
         PRINT 500,BOLKCM,'KCAL/MOL'
      ENDIF
      PRINT *

 500  FORMAT (4X,'(USE BOLTZ = ',E12.5,' FOR ',A,')')
 510  FORMAT (4X,'INITIAL VELOCITIES ARE TAKEN FROM A ',
     $     'MAXWELLIAN DISTRIBUTION',/,7X,
     $     'WITH TEMPERATURE TEMPI')

C end PSTART
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GTCOUP(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'formats.h'
C local vars
      INTEGER I,J,K,IABSVL,NNEG,II
      LOGICAL LINEG,LJNEG,LKNEG
      CHARACTER STR*(28)
      CHARACTER VSTR*(6)
      LOGICAL LCHK(NTTNUM)
C data
      DATA STR /'NTT,TEMP0,TAUT'/
C begin
C NTT,TEMP0,TAUT
      DO I=1,NTTNUM
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPINT(NTT(I),PRGSTR,STR)
         CALL CHPREL(TEMP0(I),PRGSTR,STR)
         CALL CHPREL(TAUT(I),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)
      ENDDO

      CALL GTEND(IUNIT,PRGSTR)

      DO I=1,NTTNUM
         LCHK(I) = .TRUE.
C check NTT
         IF (NTT(I) .LT. NTTMIN .OR. NTT(I) .GT. NTTMAX .OR. 
     $        NTT(I) .EQ. -NTTONE) THEN
            CALL FLAGLN(PRGSTR)
            VSTR = 'NTT'
            PRINT FMNII,VSTR,NTT(I)
            PRINT FMBTI,NTTMIN,VSTR,NTTMAX
            CALL STOPGM('GTCOUP','FATAL ERROR')
         ENDIF

C check TEMP0
         IF (TEMP0(I) .LE. 0.0D0) THEN
            CALL FLAGLN(PRGSTR)
            VSTR = 'TEMP0'
            PRINT FMNIR,VSTR, TEMP0(I)
            PRINT FMGER,VSTR
            CALL STOPGM('GTCOUP','FATAL ERROR')
         ENDIF

C check TAUT
         IF (TAUT(I) .LE. 0.0D0) THEN
            CALL FLAGLN(PRGSTR)
            VSTR = 'TAUT'
            PRINT FMNIR,VSTR,TAUT(I)
            PRINT FMGER,VSTR
            CALL STOPGM('GTCOUP','FATAL ERROR')
         ENDIF
      ENDDO

C check consistency of NTT switches
      DO 30 I=1,NTTNUM
         IF (LCHK(I) .AND. NTT(I) .NE. NTTOFF) THEN
C do something depending on value of switch
            IABSVL = ABS(NTT(I))
            LINEG   = (NTT(I) .LT. NTTOFF)
            IF (IABSVL .EQ. NTTTWO) THEN
C look for another
               J=I+1
C while loop
 50            IF (J .LE. NTTNUM) THEN
                  IF (ABS(NTT(J)) .NE. IABSVL) THEN
                     J = J+1
                     GOTO 50
                  ENDIF
               ENDIF
C end while
               IF (J .LE. NTTNUM) THEN
                  LJNEG = (NTT(J) .LT. NTTOFF)
C found, check sign
                  IF (LINEG .AND. LJNEG .OR. 
     $                 (.NOT. LINEG .AND. .NOT. LJNEG)) THEN
                     PRINT 500,PRGSTR,I,NTT(I)
                     PRINT 600
                     DO 200 II=1,NTTNUM
                        PRINT 601,NTT(II),TEMP0(II),TAUT(II)
 200                 CONTINUE
                     PRINT *,'one 2 and one -2 needed!'
                     CALL STOPGM('GTCOUP','FATAL ERROR')
                  ENDIF
C everything is OK: mark J
                  LCHK(J) = .FALSE.
               ELSE
C     no element found: error
                  PRINT 500,PRGSTR,I,NTT(I)
                  PRINT 600
                  DO 201 II=1,NTTNUM
                     PRINT 601,NTT(II),TEMP0(II),TAUT(II)
 201              CONTINUE
                  PRINT *,'no second value of 2 or -2 found!'
                  CALL STOPGM('GTCOUP','FATAL ERROR')
               ENDIF
            ELSEIF (IABSVL .EQ. NTTHRE) THEN
               IF (I .EQ. 1) THEN
                  J = I+1
                  K = J+1
                  IF (ABS(NTT(J)) .EQ. NTTHRE .AND.
     $                 ABS(NTT(K)) .EQ. NTTHRE) THEN
                     LJNEG = (NTT(J) .LT. NTTOFF)
                     LKNEG = (NTT(K) .LT. NTTOFF)
C check signs
                     NNEG = 0
                     IF (LINEG) THEN
                        NNEG = NNEG + 1
                     ENDIF
                     IF (LJNEG) THEN
                        NNEG = NNEG + 1
                     ENDIF
                     IF (LKNEG) THEN
                        NNEG = NNEG + 1
                     ENDIF
                     IF (NNEG .NE. 2) THEN
                        PRINT 500,PRGSTR,I,NTT(I)
                        PRINT 600
                        DO 205 II=1,NTTNUM
                           PRINT 601,NTT(II),TEMP0(II),TAUT(II)
 205                    CONTINUE
                        PRINT *,'one 3 and two -3 needed'
                        CALL STOPGM('GTCOUP','FATAL ERROR')
                     ENDIF
C we are happy with the input: mark the J and K items
                     LCHK(J) = .FALSE.
                     LCHK(K) = .FALSE.
                  ELSE
                     PRINT 500,PRGSTR,I,NTT(I)
                     PRINT 600
                     DO 207 II=1,NTTNUM
                        PRINT 601,NTT(II),TEMP0(II),TAUT(II)
 207                 CONTINUE
                     PRINT *,'all values must be 3 or -3'
                     CALL STOPGM('GTCOUP','FATAL ERROR')
                  ENDIF
               ELSE
                  PRINT 500,PRGSTR,I,NTT(I)
                  PRINT 600
                  DO 210 II=1,NTTNUM
                     PRINT 601,NTT(II),TEMP0(II),TAUT(II)
 210              CONTINUE
                  PRINT *,'all values must be 3 or -3'
                  CALL STOPGM('GTCOUP','FATAL ERROR')
               ENDIF
            ENDIF
         ENDIF
 30   CONTINUE

 500  FORMAT (1X,A,': inconsistent value for NTT(',I1,') =',I3)
 600  FORMAT ('       NTT     TEMP0      TAUT')
 601  FORMAT (I10,2F10.5)
C end GTCOUP
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PTCOUP
C     ==--------------------------------------------------------------==
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C local vars
      INTEGER I
C begin
      IF (LGOT(NTCOUP)) THEN
         PRINT *,BTNAME(NTCOUP)
         PRINT '(3A10)','NTT','TEMP0','TAUT'
         DO 10 I=1,NTTNUM
            PRINT '(I10,2F10.3)',NTT(I),TEMP0(I),TAUT(I)
 10      CONTINUE
C ntt
         PRINT *
      ENDIF
 500  FORMAT (4X,'NO TEMPERATURE COUPLING')
C end ptcoup
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GPCOUP(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
      INCLUDE 'box.h'
      INCLUDE 'md.h'
      INCLUDE 'formats.h'
C local vars
      INTEGER I
      LOGICAL LTC
      CHARACTER STR*(32)
      CHARACTER VSTR*(6)
C data
      DATA STR /'NTP,PRES0,COMP,TAUP'/
C begin
C NTP,PRES0,COMP,TAUP
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTP,PRGSTR,STR)
      CALL CHPREL(PRES0,PRGSTR,STR)
      CALL CHPREL(COMP,PRGSTR,STR)
      CALL CHPREL(TAUP,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C check NTP
      IF (NTP .LT. NTPMIN .OR. NTP .GT. NTPMAX ) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTP'
         PRINT FMNII,VSTR,NTP
         PRINT FMBTI,NTPMIN,VSTR,NTPMAX
         CALL STOPGM('GPCOUP','FATAL ERROR')
      ENDIF

      IF (NTP .NE. NTPOFF) THEN
C we have pressure scaling --
C check consistency with other vars
C we must have temp coupling too
         LTC = .FALSE.
         DO I=1,NTTNUM
            LTC = LTC .OR. (NTT(I) .NE. NTTOFF)
         ENDDO

         IF (.NOT. LTC) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'inconsistent coupling switches NTT(1..3) and NTP'
            PRINT *,'pressure coupling requires temperature coupling!'
            CALL STOPGM('GPCOUP','FATAL ERROR')
         ENDIF

C we must calc the virial,
C i.e. ABS(NTB) must be NTBVIR

         IF (ABS(NTB) .NE. NTBVIR) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'inconsistent switches NTB and NTP'
            PRINT 520,NTB
            PRINT 521,NTP
            PRINT *,'pressure coupling requires calculation of virial'
            CALL STOPGM('GPCOUP','FATAL ERROR')
         ENDIF

C check PRES0
         IF (PRES0 .LE. 0.0D0) THEN
            CALL FLAGLN(PRGSTR)
            VSTR = 'PRES0'
            PRINT FMNIR,VSTR,PRES0
            PRINT FMGER,VSTR
            CALL STOPGM('GPCOUP','FATAL ERROR')
         ENDIF
      ENDIF

C check COMP
      IF (COMP .LE. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'COMP'
         PRINT FMNIR,VSTR,COMP
         PRINT FMGTR,VSTR
         CALL STOPGM('GPCOUP','FATAL ERROR')
      ENDIF

C check TAUP
      VSTR = 'TAUP'
      IF (TAUP .LE. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNIR,VSTR, TAUP
         PRINT FMGTR,VSTR
         CALL STOPGM('GPCOUP','FATAL ERROR')
      ENDIF

C     only check TAUT(I) iff we use that value
      DO 20 I=1,NTTNUM
         IF (NTT(I) .GE. NTTONE .AND. TAUP .LE. TAUT(I)) THEN
            CALL FLAGLN(PRGSTR)
            PRINT FMNIR,VSTR,TAUP
            PRINT *,'pressure relaxation time TAUP must be larger'
            PRINT *,'than temperature relaxation time TAUT(I) =',TAUT(I)
            CALL STOPGM('GPCOUP','FATAL ERROR')
         ENDIF
 20   CONTINUE

 520  FORMAT (1X,'NTB =',I4)
 521  FORMAT (1X,'NTP =',I4)
C end GPCOUP
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PPCOUP
C     ==--------------------------------------------------------------==
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C begin
      IF (LGOT(NPCOUP)) THEN
         PRINT *,BTNAME(NPCOUP)
         PRINT '(4A10)','NTP','PRES0','COMP','TAUP'
         PRINT '(I10,2(1X,E9.4),F10.5)',NTP,PRES0,COMP,TAUP
         IF (NTP .EQ. NTPOFF) THEN
            PRINT 500
         ELSEIF (NTP .EQ. NTPISO) THEN
            PRINT *,'   PRESSURE COUPLING WITH ',
     $         'ISOTROPIC POSITION SCALING'
         ELSEIF (NTP .EQ. NTPANI) THEN
            PRINT *,'   PRESSURE COUPLING WITH ',
     $         'ANISOTROPIC POSITION SCALING'
         ELSE
            PRINT *,'illegal value for NTP = ',NTP
            CALL STOPGM('PPCOUP','FATAL ERROR')
         ENDIF
      ELSE
         PRINT 500
      ENDIF
      PRINT *
 500  FORMAT (4X,'NO PRESSURE COUPLING')
C end PPCOUP
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==--------------------------------------------------------------==
C-----MSDOS
C-----MSDOS The symbol table for the Microsoft compiler is too small
C-----MSDOS to compile this file.
C-----MSDOS Split the file here into rdmd1.f and rdmd2.f
C-----MSDOS and compile separately
C-----MSDOS - cmb - we do no longer care about the above
C     ==================================================================
      SUBROUTINE GNSP(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(18)
      CHARACTER VSTR*(8)
      INTEGER K
      DATA STR  /'NSPM,NSP(1..NSPM)'/
C begin
C NSPM
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NSPM,PRGSTR,STR)

C check NSPM
      IF (NSPM .LT. 1 .OR. NSPM .GT. MAXNSP) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NSPM'
         PRINT FMNII,VSTR,NSPM
         PRINT FMBTII,1,VSTR,'MAXNSP',MAXNSP
         CALL STOPGM('GNSP','FATAL ERROR')
      ENDIF

C read in NSP array
      DO K=1,NSPM
         CALL CHPINT(NSP(K),PRGSTR,STR)
C check NSP(K)
C must be 1 <= NSP(K) <= NRP, but cant check for
C <= NRP as topology hasnt been read in at this point in time.
C ==> do simple checking here, and proper checking in routine
C CHKNSP in file promd.f
         IF (NSP(K) .LE. 0) THEN
            CALL FLAGLN(PRGSTR)
            WRITE(6,*)'NSP(K) =',NSP(K),' K =',K
            WRITE(6,*)'FMNII =',FMNII,' FMGTI =',FMGTI,' STR =',STR
            CALL STOPGM('GNSP','FATAL ERROR')
         ENDIF
      ENDDO

      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C check NSP for monotonicity
      DO K=2,NSPM
         IF (NSP(K) .LE. NSP(K-1)) THEN
            CALL FLAGLN(PRGSTR)
            WRITE(6,*)'NSP(K) =',NSP(K),' NSP(K-1) =',NSP(K-1),' K =',K
            WRITE(6,*)'FMNII =',FMNII,' STR =',STR
            WRITE(6,*)'Values must increase monotonically'
            CALL STOPGM('GNSP','FATAL ERROR')
         ENDIF
      ENDDO
C end GNSP
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PNSP
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
      INCLUDE 'mdblock.h'
      INTEGER N
C begin
      IF (LGOT(NNSP)) THEN
         PRINT *,BTNAME(NNSP)
         PRINT 500,'NSPM'
         PRINT 510,NSPM
         PRINT 515,'NSP(1..NSPM)'
         PRINT 510,(NSP(N),N=1,NSPM)
         PRINT *
      ENDIF
 500  FORMAT (A10)
 510  FORMAT (8I10)
 515  FORMAT (A20)
C end PNSP
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GCOM(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(20)
      CHARACTER VSTR*(6)
      INTEGER NTCM
C data
      DATA STR /'NDFMIN,NTCM,NSCM'/
C begin
C NDFMIN,NTCM,NSCM
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NDFMIN,PRGSTR,STR)
      CALL CHPINT(NTCM,PRGSTR,STR)
      CALL CHPINT(NSCM,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C no checking of NDFMIN -- any value possible

C check NTCM
      IF (NTCM .NE. ITRUE .AND. NTCM .NE. IFALSE) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTCM'
         PRINT FMNII,VSTR,NTCM
         PRINT FMBOO,VSTR
         CALL STOPGM('GCOM','FATAL ERROR')
      ENDIF

      LTCM = (NTCM .EQ. ITRUE)

C check NSCM
      IF (NSCM .LT. 0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NSCM'
         PRINT FMNII,VSTR,NSCM
         PRINT FMGEI,VSTR
         CALL STOPGM('GCOM','FATAL ERROR')
      ENDIF
C end GCOM
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PCOM
C     ==--------------------------------------------------------------==
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C local vars
      INTEGER NTCM
C begin
      IF (LGOT(NCOM)) THEN
         IF (LTCM) THEN
            NTCM = ITRUE
         ELSE
            NTCM = IFALSE
         ENDIF
         PRINT *,BTNAME(NCOM)
         PRINT '(3A10)','NDFMIN','NTCM','NSCM'
         PRINT '(3I10)',NDFMIN,NTCM,NSCM
C ltcm
         IF (LTCM) THEN
            PRINT 500,'IS REMOVED IF INIT <>',INITCO
         ELSE
            PRINT 500,'IS NOT REMOVED'
         ENDIF
C NSCM
         IF (NSCM .EQ. 0) THEN
            PRINT *,'TRANSL. MOTION DURING RUN NOT REMOVED'
         ENDIF

         PRINT *
      ENDIF
 500  FORMAT (4X,'INITIAL C.O.M. MOTION ',A,I2)
C end PCOM
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GSTEP(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(16)
      CHARACTER VSTR*(6)
C data
      DATA STR /'NSTLIM,T,DT'/
C begin
C  NSTLIM,INIT,T,DT
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NSTLIM,PRGSTR,STR)
      CALL CHPREL(T,PRGSTR,STR)
      CALL CHPREL(DT,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C check NSTLIM
      IF (NSTLIM .LT. 0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NSTLIM'
         PRINT FMNII,VSTR,NSTLIM
         PRINT FMGEI,VSTR
         CALL STOPGM('GSTEP','FATAL ERROR')
      ENDIF

C check T
      IF (T .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'T'
         PRINT FMNIR,VSTR,T
         PRINT FMGER,VSTR
         CALL STOPGM('GSTEP','FATAL ERROR')
      ENDIF

C check DT
      IF (DT .LE. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'DT'
         PRINT FMNIR,VSTR,DT
         PRINT FMGTR,VSTR
         CALL STOPGM('GSTEP','FATAL ERROR')
      ENDIF
C end GSTEP
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PSTEP
C     ==--------------------------------------------------------------==
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C begin
      PRINT *,BTNAME(NSTEP)
      PRINT '(4A10)','NSTLIM','T','DT'
      PRINT '(I10,2F10.5)',NSTLIM,T,DT
      PRINT *
C end pstep
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GSHAKE(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(12)
      CHARACTER VSTR*(6)
C data
      DATA STR /'NTC,SHKTOL'/
C begin
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTC,PRGSTR,STR)
      CALL CHPREL(SHKTOL,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C check NTC
      IF (NTC .LT. NTCMIN .OR. NTC .GT. NTCMAX) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTC'
         PRINT FMNII,VSTR,NTC
         PRINT FMBTI,NTCMIN,VSTR,NTCMAX
         CALL STOPGM('GSHAKE','FATAL ERROR')
      ENDIF

C check SHKTOL
      IF (SHKTOL .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'SHKTOL'
         PRINT FMNIR,VSTR,SHKTOL
         PRINT FMGER,VSTR
         CALL STOPGM('GSHAKE','FATAL ERROR')
      ENDIF
C end GSHAKE
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PSHAKE
C     ==--------------------------------------------------------------==
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C begin
      PRINT *,BTNAME(NSHAKE)
      PRINT '(2A10)','NTC','TOL'
      PRINT '(I10,1X,E9.4)',NTC,SHKTOL
C ntc
      IF (NTC .EQ. NTCNON) THEN
         PRINT *,'   NO SHAKE IS PERFORMED'
      ELSE
         IF (NTC .EQ. NTCDOH) THEN
            PRINT *,'   SHAKE IS ONLY PERFORMED FOR BONDS ',
     $           'INVOLVING H-ATOMS'
         ELSEIF (NTC .EQ. NTCDOB) THEN
            PRINT *,'   SHAKE IS PERFORMED FOR ALL BONDS'
         ELSE
            PRINT *,'illegal value for NTC = ',NTC
            CALL STOPGM('PSHAKE','FATAL ERROR')
         ENDIF
      ENDIF
      PRINT *
C end pshake
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GFORCE(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(18)
      CHARACTER VSTR*(8)
      INTEGER NTMP,I
C data
      DATA STR /'NTF(1..10)'/
C begin
C NTF array
      CALL GIMME(IUNIT,PRGSTR,STR)
      DO 10 I=1,MAXNTF
         CALL CHPINT(NTMP,PRGSTR,STR)
         IF (NTMP .NE. IFALSE .AND. NTMP .NE. ITRUE) THEN
            CALL FLAGLN(PRGSTR)
            VSTR = 'NTF'
            PRINT FMNII,VSTR,NTMP
            PRINT FMBOO,VSTR
            CALL STOPGM('GFORCE','FATAL ERROR')
         ELSE
            LTF(I) = (NTMP .EQ. ITRUE)
         ENDIF
 10   CONTINUE
      CALL CHKMTY(PRGSTR,STR)

      IF (LTF(NTFCG) .AND. .NOT. LTF(NTFNBN)) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,'illegal input combination for NTF(9) and NTF(10)!'
         PRINT *,'non-bonded charges cannot be calculated without'
         PRINT *,'the non-bonded interaction'
         CALL STOPGM('GFORCE','FATAL ERROR')
      ENDIF
C next line: NEGR and NRE(..)
      CALL GIMME(IUNIT,PRGSTR,STR)
      STR = 'NEGR, NRE(1..NEGR)'
      CALL CHPINT(NEGR,PRGSTR,STR)
      IF (NEGR .LT. 0 .OR. NEGR .GT. MAXNRE) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NEGR'
         PRINT FMNII,VSTR,NEGR
         PRINT FMBTII,0,VSTR,'MAXNRE',MAXNRE
         CALL STOPGM('GFORCE','FATAL ERROR')
      ENDIF

      NTMP = 0
      DO 20 I=1,NEGR
         CALL CHPINT(NRE(I),PRGSTR,STR)

C do some checking of NRE
C NOTE: more thorough checking will have to be done
C when the topology has been read in.
C this is done in PROMD.f
         IF (NRE(I) .LT. 0) THEN
            CALL FLAGLN(PRGSTR)
            VSTR = 'NRE(..)'
            PRINT FMNII,VSTR,NRE(I)
            PRINT FMGEI,VSTR
            CALL STOPGM('GFORCE','FATAL ERROR')
         ENDIF


         IF (NTMP .GE. NRE(I)) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'inconsistent values for NRE(..)'
            PRINT *,'NRE(1) =',NRE(I)
            PRINT *,'NRE values must be in ascending order'
            CALL STOPGM('GFORCE','FATAL ERROR')
         ENDIF

         NTMP = NRE(I)
 20   CONTINUE

      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)
C end GFORCE
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PFORCE
C     ==--------------------------------------------------------------==
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C local params
      INTEGER MS
      PARAMETER (MS =16)
C local vars
      INTEGER NNTF(MAXNTF)
      CHARACTER STRNTF*(46)
      DIMENSION STRNTF(MAXNTF)
      INTEGER I
      CHARACTER STR*(MS)
      CHARACTER INCL*(MS)
      CHARACTER NOTINC*(MS)
C data
      DATA INCL   /'INCLUDES'/
      DATA NOTINC /'DOES NOT INCLUDE'/
C
      DATA STRNTF(NTFBNH) /'BOND FORCES INVOLVING H-ATOMS'/
      DATA STRNTF(NTFBND) /'BOND FORCES NOT INVOLVING H-ATOMS'/
      DATA STRNTF(NTFANH) /'BOND ANGLE FORCES INVOLVING H-ATOMS'/
      DATA STRNTF(NTFANG) /'BOND ANGLE FORCES NOT INVOLVING H-ATOMS'/
      DATA STRNTF(NTFIDH) /'IMPROPER DIHEDRAL FORCES INVOLVING H-ATOMS'/
      DATA STRNTF(NTFIDE) /
     $     'IMPROPER DIHEDRAL FORCES NOT INVOLVING H-ATOMS'/
      DATA STRNTF(NTFDHH) /'DIHEDRAL ANGLE FORCES INVOLVING H-ATOMS'/
      DATA STRNTF(NTFDIH) /
     $     'DIHEDRAL ANGLE FORCES NOT INVOLVING H-ATOMS'/
      DATA STRNTF(NTFCG)  /'CHARGES IN NON-BONDED INTERACTION'/
      DATA STRNTF(NTFNBN) /'NON-BONDED INTERACTION'/
C

C begin
      PRINT *,BTNAME(NFORCE)
      DO I=1,MAXNTF
         IF (LTF(I)) THEN
            NNTF(I) = ITRUE
         ELSE
            NNTF(I) = IFALSE
         ENDIF
      ENDDO

      PRINT '(A20)','NTF(1..10)'
      PRINT '(5A10)','BONDS','BND-ANGL','IMP-DIHE','DIHEDRAL','NONBOND'
      PRINT '(4(2X,A4,4X),(3X,A3,A4))','H','H','H','H','CHG','ALL'
      PRINT '(5(2X,2I4))',(NNTF(I),I=1,MAXNTF)
      
      PRINT '(A10,A20)','NEGR','NRE(1..NEGR)'
      PRINT '(20I10)',NEGR,(NRE(I),I=1,NEGR)

      PRINT 450
C do bonds involving H atoms
      DO I=1,MAXNTF
         IF (LTF(I)) THEN
            STR = INCL
         ELSE
            STR = NOTINC
         ENDIF
         PRINT 500,STR,STRNTF(I)
      ENDDO

      PRINT *
 450  FORMAT (4X,'FORCE CALCULATION')
 500  FORMAT (4X,A,1X,A)
C end pforce
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GPLIST(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'md.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(24)
      CHARACTER VSTR*(6)
      INTEGER NTNB
C data
      DATA STR /'NTNB,NSNB,RCUTP,RCUTL'/
C begin
C NTNB,NSNB,RCUTP,RCUTL
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTNB,PRGSTR,STR)
      CALL CHPINT(NSNB,PRGSTR,STR)
      CALL CHPREL(RCUTP,PRGSTR,STR)
      CALL CHPREL(RCUTL,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C check NTNB
      IF (NTNB .NE. ITRUE .AND. NTNB .NE. IFALSE) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTNB'
         PRINT FMNII,VSTR,NTNB
         PRINT FMBOO,VSTR
         CALL STOPGM('GPLIST','FATAL ERROR')
      ENDIF
      LTNB = (NTNB .EQ. 1)

C check NSNB
      IF (NSNB .LE. 0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NSNB'
         PRINT FMNII,VSTR,NSNB
         PRINT FMGTI,VSTR
         CALL STOPGM('GPLIST','FATAL ERROR')
      ENDIF

C check RCUTP
      IF (RCUTP .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'RCUTP'
         PRINT FMNIR,VSTR,RCUTP
         PRINT FMGER,VSTR
         CALL STOPGM('GPLIST','FATAL ERROR')
      ENDIF

C check RCUTL
      IF (RCUTL .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'RCUTL'
         PRINT FMNIR,VSTR,RCUTL
         PRINT FMGER,VSTR
         CALL STOPGM('GPLIST','FATAL ERROR')
      ENDIF

      IF (NRDBOX .EQ. NRDBXN) THEN
         CALL CHKRCT(PRGSTR,.TRUE.)
      ENDIF
C end GPLIST
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE CHKRCT(PRGSTR,LFLAG)
C     ==--------------------------------------------------------------==
C     Make sure the cut-offs are not too large for the box.
COMMEND
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'md.h'
C args
      CHARACTER PRGSTR*(*)
      LOGICAL LFLAG
C local vars
      INTEGER K
C begin
      IF (NTB .NE. NTBVAC) THEN
         DO 30 K=1,NDIM
            IF (2.0D0*RCUTL .GT. BOX(K)) THEN
               IF (LFLAG) THEN
                  CALL FLAGLN(PRGSTR)
               ENDIF
               PRINT *,'value for RCUTL =',RCUTL,' too large for BOX!'
               PRINT *,'BOX(',K,') = ',BOX(K)
               PRINT *,'2* RCUTL must be <= BOX(K) forall K'
#ifdef READ_COORD
               PRINT *,'READ_COORD is set, so I presume you recalculate'
               PRINT *,'a solvent trj without solvent... OK...'
#else               
               CALL STOPGM('CHKRCT','FATAL ERROR')
#endif
            ENDIF
 30      CONTINUE

         DO 20 K=1,NDIM
            IF (2.0D0*RCUTP .GT. BOX(K)) THEN
               IF (LFLAG) THEN
                  CALL FLAGLN(PRGSTR)
               ENDIF
               PRINT *,'value for RCUTP =',RCUTP,' too large for BOX!'
               PRINT *,'BOX(',K,') = ',BOX(K)
               PRINT *,'2* RCUTP must be <= BOX(K) forall K'
#ifdef READ_COORD
               PRINT *,'READ_COORD is set, so I presume you recalculate'
               PRINT *,'a solvent trj without solvent... OK...'
#else
               CALL STOPGM('CHKRCT','FATAL ERROR')
#endif
            ENDIF
 20      CONTINUE
      ENDIF
C end chkrct
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PPLIST
C     ==--------------------------------------------------------------==
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C local vars
      INTEGER NTNB
C begin
      IF (LTNB) THEN
         NTNB = ITRUE
      ELSE
         NTNB = IFALSE
      ENDIF
      PRINT *,BTNAME(NPLIST)
      PRINT '(4A10)','NTNB','NSNB','RCUTP','RCUTL'
      PRINT '(2I10,2F10.3)',NTNB,NSNB,RCUTP,RCUTL

      IF (LTNB) THEN
         PRINT 500,'A'
      ELSE
         PRINT 500,'NO'
      ENDIF
      IF (RCUTL .GT. RCUTP) THEN
         PRINT 510,'A'
      ELSE
         PRINT 510,'NO'
      ENDIF
      PRINT *

 500  FORMAT (4X,A,' PAIRLIST IS CONSTRUCTED IN THE FIRST STEP')
 510  FORMAT (4X,A,' LONGRANGE CONTRIBUTION IS CALCULATED ON',
     $     ' PL CONSTRUCTION')
C end pplist
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GPRINT(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(16)
C data
      DATA STR /'NTPR,NTPL,NTPP'/
C begin
C  NTPR,NTPP
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTPR,PRGSTR,STR)
      CALL CHPINT(NTPL,PRGSTR,STR)
      CALL CHPINT(NTPP,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C check NTPR
      IF (NTPR .LT. 0) THEN
         CALL FLAGLN(PRGSTR)
         STR = 'NTPR'
         PRINT FMNII,STR,NTPR
         PRINT FMGEI,STR
         CALL STOPGM('GPRINT','FATAL ERROR')
      ENDIF
C check NTPL
      IF (NTPL .LT. 0) THEN
         CALL FLAGLN(PRGSTR)
         STR = 'NTPL'
         PRINT FMNII,STR,NTPR
         PRINT FMGEI,STR
         CALL STOPGM('GPRINT','FATAL ERROR')
      ENDIF

C check NTPP
      IF (NTPP .NE. IFALSE .AND. NTPP .NE. ITRUE) THEN
         CALL FLAGLN(PRGSTR)
         STR = 'NTPP'
         PRINT FMNII,STR,NTPP
         PRINT FMBOO,STR
         CALL STOPGM('GPRINT','FATAL ERROR')
      ENDIF
C end GPRINT
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PPRINT
C     ==--------------------------------------------------------------==
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C begin
      PRINT *,BTNAME(NPRINT)
      PRINT '(3A10)','NTPR','NTPL','NTPP'
      PRINT '(3I10)',NTPR,NTPL,NTPP

      IF (NTPR .EQ. 0) THEN
         PRINT *,'   NO PRINTING OF ENERGIES'
      ENDIF
      IF (NTPL .EQ. 0) THEN
         PRINT *,'   NO TEMPERATURE LEAST SQUARES FITTING',
     $        ' AND C.O.M PRINTING'
      ENDIF
      IF (NTPP .EQ. IFALSE) THEN
         PRINT 500,'   NO'
      ELSE
         PRINT 500,'  '
      ENDIF
      PRINT *
C end pprint
 500  FORMAT (1X,A,' DIHEDRAL ANGLE MONITORING')
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GWRITE(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'forcesz.h'
      INCLUDE 'formats.h'
C local vars
#ifdef WRITE_FORCES
      CHARACTER STR*(42)
#else
      CHARACTER STR*(35)
#endif
      CHARACTER VSTR*(6)
C data
#ifdef WRITE_FORCES
      DATA STR /'NTWX,NTWSE,NTWV,NTWE,NTWF,NTWG,NTWFOR,NTPW'/
#else
      DATA STR /'NTWX,NTWSE,NTWV,NTWE,NTWF,NTWG,NTPW'/
#endif
C begin
C NTWX,NTWSE,NTWV,NTWE,NTWF,NTPW
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTWX,PRGSTR,STR)
      CALL CHPINT(NTWSE,PRGSTR,STR)
      CALL CHPINT(NTWV,PRGSTR,STR)
      CALL CHPINT(NTWE,PRGSTR,STR)
      CALL CHPINT(NTWG,PRGSTR,STR)
#ifdef WRITE_FORCES
      CALL CHPINT(NTWFOR,PRGSTR,STR)
#endif
      CALL CHPINT(NTPW,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

      IF (NTPW .LT. NTPWMI .OR. NTPW .GT. NTPWMA) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTPW'
         PRINT FMNII,VSTR,NTPW
         PRINT FMBTI,NTPWMI,VSTR,NTPWMA
         CALL STOPGM('GWRITE','FATAL ERROR')
      ENDIF

C check for ranges
C we check NTWSE properly in GLBCHK
      IF (NTWSE .LT. 0 .OR. NTWSE .GT. MXEWRT) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTWSE'
         PRINT FMNII,VSTR,NTWSE
         PRINT FMBTI,0,VSTR,MXEWRT
         CALL STOPGM('GWRITE','FATAL ERROR')
      ENDIF

      IF (NTWE .LT. 0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTWE'
         PRINT FMNII,VSTR,NTWE
         PRINT FMGEI,VSTR
         CALL STOPGM('GWRITE','FATAL ERROR')
      ENDIF

      IF (NTWX .EQ. 0 .AND. NTWSE .NE. 0)THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,'inconsistent values for NTWSE and NTWX'
         PRINT *,'NTWX =',NTWX
         PRINT *,'NTWSE=',NTWSE
         PRINT *,'NTWX must be > 0 if NTWSE > 0'
         CALL STOPGM('GWRITE','FATAL ERROR')
      ENDIF

      IF (NTWG .LT. 0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTWG'
         PRINT FMNII,VSTR,NTWG
         PRINT FMGEI,VSTR
         CALL STOPGM('GWRITE','FATAL ERROR')
      ENDIF

#include "inc_rdmd6.h"

C end GWRITE
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PWRITE
C     ==--------------------------------------------------------------==
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C local vars
      CHARACTER  STRWAT*(12)
      CHARACTER  STRSSS*(24)
      CHARACTER  STRSOL*(24)
      DATA STRSSS /'SOLUTE AND SOLVENT'/
      DATA STRSOL /'            SOLUTE'/
C begin
      IF (LGOT(NWRITE)) THEN
         PRINT *,BTNAME(NWRITE)
#ifdef WRITE_FORCES
         PRINT '(8A10)','NTWX','NTWSE','NTWV',
     $        'NTWE','NTWG','NTWFOR','NTPW'
         PRINT '(8I10)',NTWX,NTWSE,NTWV,NTWE,NTWG,NTWFOR,NTPW
#else
         PRINT '(7A10)','NTWX','NTWSE','NTWV',
     $        'NTWE','NTWG','NTPW'
         PRINT '(7I10)',NTWX,NTWSE,NTWV,NTWE,NTWG,NTPW
#endif
         STRWAT = 'COORDINATES'
         IF (NTWX .EQ. 0) THEN
            PRINT 500,STRWAT
         ELSEIF (NTWX .GT. 0) THEN
            PRINT 510,STRSSS,STRWAT
         ELSE
            PRINT 510,STRSOL,STRWAT
         ENDIF

         STRWAT = ' VELOCITIES'
         IF (NTWV .EQ. 0) THEN
            PRINT 500,STRWAT
         ELSEIF (NTWV .GT. 0) THEN
            PRINT 510,STRSSS,STRWAT
         ELSE
            PRINT 510,STRSOL,STRWAT
         ENDIF

         STRWAT = '   ENERGIES'
         IF (NTWE .EQ. 0) THEN
            PRINT 500,STRWAT
         ELSE
            PRINT 520,STRWAT
         ENDIF

         IF (NTPW .EQ. NTPWBI) THEN
            PRINT *,'   DATA WRITTEN IN BINARY FORM'
         ELSEIF (NTPW .EQ. NTPWFO) THEN
            PRINT *,'   DATA WRITTEN IN FORMATTED FORM'
         ELSE
            PRINT *,'unknown option for NTPW!'
            CALL STOPGM('PWRITE','FATAL ERROR')
         ENDIF
         PRINT *
      ENDIF

 500  FORMAT (4X,19X,   A12,' NOT WRITTEN TO FILE')
 510  FORMAT (4X,A18,1X,A12,'     WRITTEN TO FILE')
 520  FORMAT (4X,19X,   A12,'     WRITTEN TO FILE')
C end pwrite
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GPREST(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(14)
      CHARACTER VSTR*(6)
C data
      DATA STR /'NTR,CHO,NRDRX'/
C begin
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTR,PRGSTR,STR)
      CALL CHPREL(CHO,PRGSTR,STR)
      CALL CHPINT(NRDRX,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C check NTR
      IF (NTR .LT. NTRMIN .OR. NTR .GT. NTRMAX) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTR'
         PRINT FMNII,VSTR,NTR
         PRINT FMBTI,NTRMIN,VSTR,NTRMAX
         CALL STOPGM('GPREST','FATAL ERROR')
      ENDIF

C check NRDRX
      IF (NRDRX .NE. IFALSE .AND. NRDRX .NE. ITRUE) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NRDRX'
         PRINT FMNII,VSTR,NRDRX
         PRINT FMBOO,VSTR
         CALL STOPGM('GPREST','FATAL ERROR')
      ENDIF

C check CHO
      IF (CHO .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'CHO'
         PRINT FMNIR,VSTR,CHO
         PRINT FMGER,VSTR
         CALL STOPGM('GPREST','FATAL ERROR')
      ENDIF

      LGOT(NPREST) = .TRUE.
C end PREST
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PPREST
C     ==--------------------------------------------------------------==
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
      INCLUDE 'formats.h'
      CHARACTER STR*(24)
C begin
      IF (LGOT(NPREST)) THEN
         PRINT *,BTNAME(NPREST)
         PRINT '(3A10)','NTR','CHO','NRDRX'
         PRINT '(I10,1X,E9.4,I10)',NTR,CHO,NRDRX
         IF (NTR .EQ. NTROFF) THEN
            PRINT 500
         ELSE
            IF (NTR .EQ. NTRCHO) THEN
               STR = 'CHO'
            ELSEIF (NTR .EQ. NTRBFA) THEN
               STR = 'CHO AND ATOMIC B-FACTORS'
            ELSEIF (NTR .EQ. NTRCON) THEN
               STR = 'CONSTRAINTS'
            ELSE
               PRINT FMNII,'NTR',NTR
               CALL STOPGM('GPREST','FATAL ERROR')
            ENDIF
            PRINT *,'   POSITION RESTRAINING USING ',STR
         ENDIF

         IF (NRDRX .EQ. NRDRXY) THEN
            PRINT *,'   ATOM RESTRAINING REFERENCE POSITIONS ARE ',
     $           'READ FROM FILE'
         ENDIF
      ELSE
         PRINT 500
      ENDIF
      PRINT *
 500  FORMAT (4X,'NO POSITION RESTRAINING')
C end pprest
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GDREST(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(26)
      CHARACTER VSTR*(6)
      INTEGER NRDDR
C data
      DATA STR /'NTDR,CDIS,DR0,TAUDR,NRDDR'/
C begin
C  NTDR,CDIS,DR0,TAUDR,NRDDR
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTDR,PRGSTR,STR)
      CALL CHPREL(CDIS,PRGSTR,STR)
      CALL CHPREL(DR0,PRGSTR,STR)
      CALL CHPREL(TAUDR,PRGSTR,STR)
      CALL CHPINT(NRDDR,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C check NTDR
      IF (NTDR .LT. NTDRMI .OR. NTDR .GT. NTDRMA) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,'NTDR',NTDR
         PRINT FMBTI,NTDRMI,'NTDR',NTDRMA
         CALL STOPGM('GDREST','FATAL ERROR')
      ENDIF
C check DR0
      IF (DR0 .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'DR0'
         PRINT FMNIR,VSTR,DR0
         PRINT FMGER,VSTR
         CALL STOPGM('GDREST','FATAL ERROR')
      ENDIF

C check TAUDR
      IF (TAUDR .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'TAUDR'
         PRINT FMNIR,VSTR,TAUDR
         PRINT FMGER,VSTR
         CALL STOPGM('GDREST','FATAL ERROR')
      ENDIF

C check NRDDR
      IF (NRDDR .NE. 0 .AND.  NRDDR .NE. 1) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NRDDR'
         PRINT FMNII,VSTR,NRDDR
         PRINT FMBOO,VSTR
         CALL STOPGM('GDREST','FATAL ERROR')
      ENDIF
      LRDDR = (NRDDR .EQ. 1)

C check CDIS
      IF (CDIS .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'CDIS'
         PRINT FMNIR,VSTR,CDIS
         PRINT FMGER,VSTR
         CALL STOPGM('GDREST','FATAL ERROR')
      ENDIF

      LGOT(NDREST) = .TRUE.
C end GDREST
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PDREST
C     ==--------------------------------------------------------------==
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
      INCLUDE 'formats.h'
C local vars
      INTEGER NRDDR
C begin
      IF (LGOT(NDREST)) THEN
         IF (LRDDR) THEN
            NRDDR = ITRUE
         ELSE
            NRDDR = IFALSE
         ENDIF
         PRINT *,BTNAME(NDREST)
         PRINT '(5A10)','NTDR','CDIS','DR0',
     $        'TAUDR','NRDDR'
         PRINT '(I10,1X,E9.4,2F10.5,I10)',NTDR,CDIS,DR0,TAUDR,NRDDR

         IF (NTDR .EQ. NTDROF) THEN
            PRINT 500
         ELSE
            IF (NTDR .LT. 0) THEN
               PRINT 520
            ELSE
               PRINT 530
            ENDIF
            IF(ABS(NTDR) .EQ. NTDRCD) THEN
               PRINT 550,'CDIS'
            ELSEIF(ABS(NTDR) .EQ. NTDRW0) THEN
               PRINT 550,'CDIS * WEIGHT FACTORS READ FROM FILE'
            ELSE
               PRINT FMNII,'NTDR',NTDR
               CALL STOPGM('PDREST','FATAL ERROR')
            ENDIF
         ENDIF
         IF (NTDR .LT. NTDROF .AND. LRDDR) THEN
            PRINT 510
         ENDIF
      ELSE
         PRINT 500
      ENDIF
      PRINT *

 500  FORMAT (4X,'NO DISTANCE RESTRAINING')
 510  FORMAT (4X,
     $     'DISTANCE AVERAGES READ FROM INITIAL COORDINATES FILE')
 520  FORMAT (4X,'TIME AVERAGED DISTANCE RESTRAINING')
 530  FORMAT (4X,'INSTANTANEOUS DISTANCE RESTRAINING')
 550  FORMAT (4X,'DISTANCE RESTRAINING USING FORCE CONSTANT ',A)

C end pdrest
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GPERT(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(36)
      CHARACTER VSTR*(6)
C dat
      DATA STR /'NTG,NRDGL,RLAM,DLAMT,RMU,DMUT'/
C begin
C NTG,NTGL,RLAM,DLAMT
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTG,PRGSTR,STR)
      CALL CHPINT(NRDGL,PRGSTR,STR)
      CALL CHPREL(RLAM,PRGSTR,STR)
      CALL CHPREL(DLAMT,PRGSTR,STR)
      CALL CHPREL(RMU,PRGSTR,STR)
      CALL CHPREL(DMUT,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)

C check NTG
      IF (NTG .LT. NTGMIN .OR. NTG .GT. NTGMAX) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTG'
         PRINT FMNII,VSTR,NTG
         PRINT *,VSTR,' must be 0, 1 or 2'
         CALL STOPGM('GPERT','FATAL ERROR')
      ENDIF

C check NTGL
      IF (NRDGL .NE. IFALSE .AND. NRDGL .NE. ITRUE) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NRDGL'
         PRINT FMNII,VSTR,NRDGL
         PRINT FMBOO,VSTR
         CALL STOPGM('GPERT','FATAL ERROR')
      ENDIF

C check RLAM
      IF (RLAM .LT. 0.0D0 .OR. RLAM .GT. 1.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'RLAM'
         PRINT FMNIR,VSTR,RLAM
         PRINT FMBTR,0.0D0,VSTR,1.0D0
         CALL STOPGM('GPERT','FATAL ERROR')
      ENDIF

C check RMU
      IF (RMU .LT. 0.0D0 .OR. RMU .GT. 1.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'RMU'
         PRINT FMNIR,VSTR,RMU
         PRINT FMBTR,0.0D0,VSTR,1.0D0
         CALL STOPGM('GPERT','FATAL ERROR')
      ENDIF

C we cant have both DLAMT and DMUT <> 0.0
      IF (DLAMT .NE. 0.0D0 .AND. DMUT .NE. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'DLAMT'
         PRINT FMNIR,VSTR,DLAMT
         VSTR = 'DMUT'
         PRINT FMNIR,VSTR,DMUT
         PRINT *,'only ONE (not both values) may be non zero'
         CALL STOPGM('GPERT','FATAL ERROR')
      ENDIF

C read in a second line
      STR = 'ALPHLJ, ALPHC, NLAM, MMU'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPREL(ALPHLJ,PRGSTR,STR)
      CALL CHPREL(ALPHC,PRGSTR,STR)
      CALL CHPINT(NLAM,PRGSTR,STR)
      CALL CHPINT(MMU,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)

      CALL GTEND(IUNIT,PRGSTR)
C check ALPHLJ
      IF (ALPHLJ .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'ALPHLJ'
         PRINT FMNIR,VSTR,ALPHLJ
         PRINT FMGER,VSTR
         CALL STOPGM('GPERT','FATAL ERROR')
      ENDIF

C check ALPHC
      IF (ALPHC .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'ALPHC'
         PRINT FMNIR,VSTR,ALPHC
         PRINT FMGER,VSTR
         CALL STOPGM('GPERT','FATAL ERROR')
      ENDIF

C check NLAM
      IF (NLAM .LE. 0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NLAM'
         PRINT FMNII,VSTR,NLAM
         PRINT FMGTI,VSTR
      ENDIF

C check MMU
      IF (MMU .LE. 0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'MMU'
         PRINT FMNII,VSTR,MMU
         PRINT FMGTI,VSTR
      ENDIF

      LGOT(NPERT) = .TRUE.
C end GPERT
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PPERT
C     ==--------------------------------------------------------------==
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
      INCLUDE 'formats.h'
C begin
      IF (LGOT(NPERT)) THEN
         PRINT *,BTNAME(NPERT)
         PRINT '(6A10)','NTG','NRDGL','RLAM',
     $        'DLAMT','RMU','DMUT'
         PRINT '(2I10,4(1X,E9.4))',NTG,NRDGL,RLAM,DLAMT,RMU,DMUT
         PRINT '(4A10)','ALPHLJ','ALPHC','NLAM','MMU'
         PRINT '(2F10.5,2I10)',ALPHLJ,ALPHC,NLAM,MMU

         IF (NTG .EQ. NTGOFF) THEN
            PRINT 500
         ELSE
            IF (NTG .EQ. NTGLAM) THEN
               PRINT 550,'LAMDA'
               PRINT 600,'RMU'
            ELSEIF (NTG .EQ. NTGMU) THEN
               PRINT 550,'MU'
               PRINT 600,'RLAM'
            ELSEIF (NTG .EQ. NTGBOT) THEN
               PRINT 550,'LAMDA'
               PRINT 550,'MU'
            ELSE
               PRINT FMNII,'NTG',NTG
               CALL STOPGM('GPERT','FATAL ERROR')
            ENDIF

            IF (NRDGL .EQ. NRDGLY) THEN
               PRINT 610
            ENDIF
         ENDIF
      ELSE
         PRINT 500
      ENDIF

      PRINT *
 500  FORMAT (4X,'NO PERTURBATION PERFORMED')
 550  FORMAT (4X,'PERTURBATION IS PERFORMED FOR ',A)
 600  FORMAT (4X,'AT ',A,' = 0.0')
 610  FORMAT (4X,'RLAM, RMU, DSVDL, DSKDL, DSVDM AND DSKDM',
     $     ' ARE READ FROM FILE (TAPE 21)')

C end ppert
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GDHRES(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(10)
      CHARACTER VSTR*(5)
C data
      DATA STR /'NTDLR,CDLR'/
C begin
C NTDLR,CDLR
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTDLR,PRGSTR,STR)
      CALL CHPREL(CDLR,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C check NTDLR
      IF (NTDLR .LT. NTDLR1 .OR. NTDLR .GT. NTDLR2) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTDLR'
         PRINT FMNII,VSTR,NTDLR
         PRINT FMBTI,NTDLR1,VSTR,NTDLR2
         CALL STOPGM('GDHRES','FATAL ERROR')
      ENDIF

C check CDLR
      IF (ALPHLJ .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'CDLR'
         PRINT FMNIR,VSTR,CDLR
         PRINT FMGER,VSTR
         CALL STOPGM('GDHRES','FATAL ERROR')
      ENDIF

      LGOT(NDHRES) = .TRUE.
C end GDHRES
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PDHRES
C     ==--------------------------------------------------------------==
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C begin
      IF (LGOT(NDHRES)) THEN
         PRINT *,BTNAME(NDHRES)
         PRINT '(2A10)','NTDLR','CDLR'
         PRINT '(I10,1X,E9.4)',NTDLR,CDLR
         IF (NTDLR .EQ. NTDLRN) THEN
            PRINT 500
         ELSEIF  (NTDLR .EQ. NTDLRC) THEN
            PRINT 510,'CDLR'
         ELSEIF  (NTDLR .EQ. NTDLRF) THEN
            PRINT 510,'CDLR*CPLR'
         ELSE
            PRINT *,'invalid value for NTDLR = ',NTDLR
            CALL STOPGM('PDHRES','FATAL ERROR')
         ENDIF
      ELSE
         PRINT 500
      ENDIF
      PRINT *
 500  FORMAT (4X,'NO DIHEDRAL RESTRAINING')
 510  FORMAT (4X,'DIHEDRAL RESTRAINING USING ',A)
C end pdhres
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
C this is just a template for adding
C block types.
C      SUBROUTINE G(IUNIT,PRGSTR)

C args
C      INTEGER IUNIT
C      CHARACTER PRGSTR*(*)
C includes
C      INCLUDE 'md.h'
C      INCLUDE 'formats.h'
C
C local vars
C      CHARACTER STR*(20)

C begin

C end G
C      END
C     ==================================================================
      SUBROUTINE GEM(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(30)
      CHARACTER VSTR*(6)
C data
      DATA STR /'NTEM, NCYC, DELE, DX0, DXM'/
C begin
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTEM,PRGSTR,STR)
      CALL CHPINT(NCYC,PRGSTR,STR)
      CALL CHPREL(DELE,PRGSTR,STR)
      CALL CHPREL(DX0,PRGSTR,STR)
      CALL CHPREL(DXM,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C check NTEM
      IF (NTEM .LT. NTEMMI .OR. NTEM .GT. NTEMMA) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,'NTEM',NTEM
         PRINT FMBTI,NTEMMI,'NTEM',NTEMMA
         CALL STOPGM('GEM','FATAL ERROR')
      ENDIF

C check NCYC
      IF (NCYC .LE. 0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NCYC'
         PRINT FMNII,VSTR,NCYC
         PRINT FMGTI,VSTR
         CALL STOPGM('GEM','FATAL ERROR')
      ENDIF

C check DELE
      IF (DELE .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'DELE'
         PRINT FMNIR,VSTR,DELE
         PRINT FMGER,VSTR
         CALL STOPGM('GEM','FATAL ERROR')
      ENDIF

C check DX0
      IF (DX0 .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'DX0'
         PRINT FMNIR,VSTR,DX0
         PRINT FMGER,VSTR
         CALL STOPGM('GEM','FATAL ERROR')
      ENDIF

C check DXM
      IF (DXM .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'DXM'
         PRINT FMNIR,VSTR,DXM
         PRINT FMGER,VSTR
         CALL STOPGM('GEM','FATAL ERROR')
      ENDIF

      IF (DX0 .GT. DXM) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,'inconsistent input for DX0 and DXM'
         PRINT *,'DXM =',DXM
         PRINT *,'DX0 =',DX0
         PRINT *,'DX0 must be <= DXM'
         CALL STOPGM('GEM','FATAL ERROR')
      ENDIF
C end GEM
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PEM
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C begin
      IF (LGOT(NEMBLK)) THEN
         PRINT *,BTNAME(NEMBLK)
         PRINT '(5A10)','NTEM', 'NCYC', 'DELE', 'DX0', 'DXM'
         PRINT '(2I10,3F10.5)',NTEM, NCYC, DELE, DX0, DXM
         IF (NTEM .EQ. NTEMSD) THEN
            PRINT 500
         ELSEIF (NTEM .EQ. NTEMCG) THEN
C     a conjugate gradients
            PRINT 510
         ELSE
            PRINT *,'PEM: illegal value for NTEM =',NTEM
            CALL STOPGM('PEM','FATAL ERROR')
         ENDIF
      ENDIF
 500  FORMAT (4X,'PERFORMING A STEEPEST DESCENT MINIMISATION')
 510  FORMAT (4X,'PERFORMING A CONJUGATE GRADIENT MINIMISATION')
C end pem
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GSD(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(32)
      CHARACTER VSTR*(6)
C data
      DATA STR /'NTFR, NSFR, NBREF, RCUTF, CFRIC'/
C begin
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTFR,PRGSTR,STR)
      CALL CHPINT(NSFR,PRGSTR,STR)
      CALL CHPINT(NBREF,PRGSTR,STR)
      CALL CHPREL(RCUTF,PRGSTR,STR)
      CALL CHPREL(CFRIC,PRGSTR,STR)
      CALL CHPREL(TEMPSD,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C check NTFR
      IF (NTFR .LT. NTFRMI .OR. NTFR .GT. NTFRMA) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTFR'
         PRINT FMNII,VSTR,NTFR
         PRINT FMBTI,NTFRMI,VSTR,NTFRMA
         CALL STOPGM('GSD','FATAL ERROR')
      ENDIF

C check NSFR
      IF (NSFR .LE. 0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NSFR'
         PRINT FMNII,VSTR,NSFR
         PRINT FMGTI,VSTR
         CALL STOPGM('GSD','FATAL ERROR')
      ENDIF

C check NBREF
      IF (NBREF .LE. 0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NBREF'
         PRINT FMNII,VSTR,NBREF
         PRINT FMGTI,VSTR
         CALL STOPGM('GSD','FATAL ERROR')
      ENDIF

C check RCUTF
      IF (RCUTF .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'RCUTF'
         PRINT FMNIR,VSTR,RCUTF
         PRINT FMGER,VSTR
         CALL STOPGM('GSD','FATAL ERROR')
      ENDIF

C check CFRIC
      IF (CFRIC .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'CFRIC'
         PRINT FMNIR,VSTR,CFRIC
         PRINT FMGER,VSTR
         CALL STOPGM('GSD','FATAL ERROR')
      ENDIF

C check TEMPSD
      IF (TEMPSD .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'TEMPSD'
         PRINT FMNIR,VSTR,TEMPSD
         PRINT FMGER,VSTR
         CALL STOPGM('GSD','FATAL ERROR')
      ENDIF
C end GSD
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PSD
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C begin
      IF (LGOT(NSDBLK)) THEN
         PRINT *,BTNAME(NSDBLK)
         PRINT '(6A10)','NTFR','NSFR','NBREF','RCUTF',
     $        'CFRIC','TEMPSD'
         PRINT '(3I10,3F10.3)',NTFR,NSFR,NBREF,RCUTF,CFRIC,TEMPSD
C NTFR
         IF (NTFR .EQ. NTFR0) THEN
            PRINT 500,'SET TO 0.0'
         ELSEIF (NTFR .EQ. NTFRFR) THEN
            PRINT 500,'SET TO CFRIC'
         ELSEIF (NTFR .EQ. NTFRGA) THEN
            PRINT 500,'SET TO CFRIC*GAM'
         ELSEIF (NTFR .EQ. NTFRCA) THEN
            PRINT 500,'CALCULATED USING FRIC'
         ELSE
            PRINT *,'PSD: invalid value for NTFR = ',NTFR
            CALL STOPGM('PSD','FATAL ERROR')
         ENDIF
         PRINT *
      ENDIF

 500  FORMAT (4X,'ATOMIC FRICTION COEFFICIENTS ',A)
C end psd
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GLRAN(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'md.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(20)
      CHARACTER VSTR*(6)
C data
      DATA STR /'EPSRF, APPAK, RCRF'/
C begin
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPREL(EPSRF,PRGSTR,STR)
      CALL CHPREL(APPAK,PRGSTR,STR)
      CALL CHPREL(RCRF,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C check EPSRF
      IF ((EPSRF .LT. 0.0D0)
     $     .OR. (EPSRF.GT.0.0D0 .AND. EPSRF.LT.1.0D0)) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'EPSREF'
         PRINT FMNIR,VSTR,EPSRF
         PRINT *,VSTR,' must be = 0.0, or >= 1.0'
         CALL STOPGM('GLRAN','FATAL ERROR')
      ENDIF

C check APPAK
      IF (APPAK .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'APPAK'
         PRINT FMNIR,VSTR,APPAK
         PRINT FMGER,VSTR
         CALL STOPGM('GLRAN','FATAL ERROR')
      ENDIF

C check RCRF
      IF (RCRF .LE. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'RCRF'
         PRINT FMNIR,VSTR,RCRF
         PRINT FMGTR,VSTR
         CALL STOPGM('GLRAN','FATAL ERROR')
      ENDIF
C end GLRAN
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PLRAN
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C begin
C     this is necessary to avois a crash when reading OLD promd control
C     files
      IF (.NOT. LGOT(NLONGR)) THEN
         RCRF = RCUTL
      ENDIF
      PRINT *,BTNAME(NLONGR)
      PRINT '(3A10)','EPSRF','APPAK','RCRF'
      PRINT '(2F10.3,E10.3)',EPSRF,APPAK,RCRF
      IF (EPSRF .EQ. 1.0D0 .AND. APPAK .EQ. 0.0D0) THEN
         PRINT 500,'NO'
      ELSE
         PRINT 500,'A'
      ENDIF
      PRINT *
 500  FORMAT (4X,A,' REACTION FIELD CORRECTION IS APPLIED')
C end plran
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GJVAL(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      CHARACTER STR*(26)
      CHARACTER VSTR*(6)
C data
      DATA STR /'NTJR,NTJRH,CJR,TAUJR,NRDJR'/
C begin
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTJR,PRGSTR,STR)
      CALL CHPINT(NTJRH,PRGSTR,STR)
      CALL CHPREL(CJR,PRGSTR,STR)
      CALL CHPREL(TAUJR,PRGSTR,STR)
      CALL CHPINT(NRDJR,PRGSTR,STR)

      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C check NTJR
      IF (NTJR .LT. NTJRMI .OR. NTJR .GT. NTJRMA) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTJR'
         PRINT FMNII,VSTR,NTJR
         PRINT FMBTI,NTJRMI,VSTR,NTJRMA
         CALL STOPGM('GJVAL','FATAL ERROR')
      ENDIF

C check NTJRH
      IF (NTJRH .NE. NTJRHH .AND. NTJRH .NE. NTJRHF) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTJRH'
         PRINT FMNII,VSTR,NTJRH
         PRINT FMBOO,VSTR
         CALL STOPGM('GJVAL','FATAL ERROR')
      ENDIF

C check TAUJR
      IF (TAUJR .LE. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'TAUJR'
         PRINT FMNIR,VSTR,TAUJR
         PRINT FMGER,VSTR
         CALL STOPGM('GJVAL','FATAL ERROR')
      ENDIF

C check NRDJR
      IF (NRDJR .NE. IFALSE .AND. NRDJR .NE. ITRUE) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NRDJR'
         PRINT FMNII,VSTR,NRDJR
         PRINT FMBOO,VSTR
         CALL STOPGM('GJVAL','FATAL ERROR')
      ENDIF

      LGOT(NJVAL) = .TRUE.
C end gjval
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PJVAL
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C begin
      IF (LGOT(NJVAL)) THEN
         PRINT *,BTNAME(NJVAL)
         PRINT '(5A10)','NTJR','NTJRH','CJR','TAUJR','NRDJR'
         PRINT '(2I10,1X,E9.4,F10.3,I10)',NTJR,NTJRH,CJR,TAUJR,NRDJR
C
         IF (NTJR .EQ. NTJROF) THEN
            PRINT 500
         ELSE
            IF  (NTJR .LT. 0) THEN
               PRINT *,'   TIME AVERAGED J-VALUE RESTRAINING'
               IF (NRDJR .EQ. NRDJRY) THEN
                  PRINT 510
               ELSE
                  PRINT *,'   COSQAV,COSIAV NOT READ IN'
               ENDIF
            ELSEIF (NTJR .GT. 0) THEN
               PRINT *,'   INSTANTANEOUS J-VALUE RESTRAINING'
            ENDIF
            IF (ABS(NTJR) .EQ. NTJCJR) THEN
               PRINT *,'   USING CJR AS FORCE CONSTANT'
            ELSEIF (ABS(NTJR) .EQ. NTJRED) THEN
               PRINT *,'   USING CJR*CPJV AS FORCE CONSTANT'
            ELSE
               PRINT *,'UNKNOWN VALUE FOR NTJR =',NTJR
               CALL STOPGM('PJVAL','FATAL ERROR')
            ENDIF

            IF (NTJRH .EQ. NTJRHH) THEN
               PRINT *,'   A HALF HARMONIC POTENTIAL IS USED'
            ELSE
               PRINT *,'   A HARMONIC POTENTIAL IS USED'
            ENDIF
         ENDIF
      ELSE
         PRINT 500
      ENDIF
      PRINT *

 500  FORMAT (4X,'NO J-VALUE RESTRAINING')
 510  FORMAT (4X,'COSQAV,COSIAV READ IN FOR CONTINUATION RUN')
C end pjval
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GLE(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'md.h'
      INCLUDE 'formats.h'
      INCLUDE 'mdblock.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      CHARACTER STR*(18)
      CHARACTER VSTR*(5)
C data
      DATA STR /'NTLE, CWLE, NRDLE'/
C begin
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTLE,PRGSTR,STR)
      CALL CHPREL(CWLE,PRGSTR,STR)
      CALL CHPINT(NRDLE,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C check NTLE
      IF (NTLE .LT. NTLEMI .OR. NTLE .GT. NTLEMA) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTLE'
         PRINT FMNII,VSTR,NTLE
         PRINT FMBTI,NTLEMI,VSTR,NTLEMA
         CALL STOPGM('GLE','FATAL ERROR')
      ENDIF

C check CWLE
      IF (CWLE .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'CWLE'
         PRINT FMNIR,VSTR,CWLE
         PRINT FMGER,VSTR
         CALL STOPGM('GLE','FATAL ERROR')
      ENDIF
C check NRDLE
      IF (NRDLE .NE. NRDLEN .AND. NRDLE .NE. NRDLEY)THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NRDLE'
         PRINT FMNII,VSTR,NRDLE
         PRINT FMBOO,VSTR
         CALL STOPGM('GLE','FATAL ERROR')
      ENDIF

      LGOT(NLEBLK) = .TRUE.
C end gle
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PLE
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
C begin
      IF (LGOT(NLEBLK)) THEN
         PRINT *,BTNAME(NLEBLK)
         PRINT '(3A10)','NTLE','CWLE','NRDLE'
         PRINT '(I10,1X,E9.4,I10)',NTLE,CWLE,NRDLE
C ntle
         IF (NTLE .EQ. NTLEOF) THEN
            PRINT 500
         ELSEIF (NTLE .EQ. NTLEGS) THEN
            PRINT 501,'A GAUSSIAN'
         ELSEIF (NTLE .EQ. NTLEQU) THEN
            PRINT 501,'AN INVERSE QUADRATIC'
         ELSE
            PRINT *,'illegal value for NTLE =',NTLE
            CALL STOPGM('PLE','FATAL ERROR')
         ENDIF
C nrdle
         IF (NRDLE .EQ. NRDLEY) THEN
            PRINT 510
         ENDIF
      ELSE
         PRINT 500
      ENDIF
      PRINT *
 500  FORMAT (4X,'NO LOCAL ELEVATION INTERACTION')
 501  FORMAT (4X,'LOCAL ELEVATION USING ',A,' POTENTIAL')
 510  FORMAT (4X,'LOCAL ELEVATION MEMORY READ FROM STARTUP FILE')
C end ple
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GFOURD(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      CHARACTER STR*(40)
      CHARACTER VSTR*(6)
      INTEGER II,IDO4D(MX4DIN)
C data
      DATA STR /'NT4DIM,CW4DA,TEMP4I,NDFMI4,NT4XI,NT4XO'/
C begin
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NT4DIM,PRGSTR,STR)
      CALL CHPREL(CW4DA,PRGSTR,STR)
      CALL CHPREL(TEMP4I,PRGSTR,STR)
      CALL CHPINT(NDFMI4,PRGSTR,STR)
      CALL CHPINT(NT4XI,PRGSTR,STR)
      CALL CHPINT(NT4XO,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)

C test NT4DIM
      VSTR = 'NT4DIM'
      IF (NT4DIM .LT. NT4MIN .OR. NT4DIM .GT. NT4MAX) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,VSTR,NT4DIM
         PRINT FMBTI,NT4MIN,VSTR,NT4MAX
         CALL STOPGM('GFOURD','FATAL ERROR')
      ENDIF

      LGOT(NFOURD) = .TRUE.

      IF (NT4DIM .NE. NT4OFF) THEN
         NDIM = 4
         NDHMIN = NDIM
         NDHMAX = NDIM

C set the fourth dimension to a large number
         BOX(NDIM) = BOXW
         BOXH(NDIM)= BOX(NDIM)*0.5D0
         BOXINV(NDIM)= 1.0D0/BOX(NDIM)

C     test CW4DA
         IF (CW4DA .LT. 0.0D0) THEN
            CALL FLAGLN(PRGSTR)
            VSTR = 'CW4DA'
            PRINT FMNIR,VSTR,CW4DA
            PRINT FMGER,VSTR
            CALL STOPGM('GFOURD','FATAL ERROR')
         ENDIF

C     test TEMP4I
         IF (TEMP4I .LT. 0.0D0) THEN
            CALL FLAGLN(PRGSTR)
            VSTR = 'TEMP4I'
            PRINT FMNIR,VSTR,TEMP4I
            PRINT FMGER,VSTR
            CALL STOPGM('GFOURD','FATAL ERROR')
         ENDIF

C no checking of NDFMI4 -- any value possible
C
C     test for implementation restrictions
C     no virial in 4D
         IF (ABS(NTB) .EQ. NTBVIR) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'Virial cannot be calculated in 4D'
            PRINT *,'NTB    =',NTB
            PRINT *,'NT4DIM =',NT4DIM
            CALL STOPGM('GFOURD','FATAL ERROR')
         ENDIF
C     beta must be 90.0
         IF (NTB .GT. 0 .AND. BETA .NE. 9.0D+1) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'Monoclinic coordinate systems not possible in 4D'
            PRINT *,'NTB    =',NTB
            PRINT *,'BETA   =',BETA
            PRINT *,'NT4DIM =',NT4DIM
            CALL STOPGM('GFOURD','FATAL ERROR')
         ENDIF
      ENDIF
C     no pressure coupling in 4D.
C     test for this in PFOURD, as it is possible
C     that a PCOUPLE block is read after a FOURD one...
C     This is not possible for the BOUND case.


C check NT4XI
      IF (NT4XI .LT. NT4XMI .OR. NT4XI .GT. NT4XMA) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NT4XI'
         PRINT FMNII,VSTR,NT4XI
         PRINT FMBTI,NT4XI,VSTR,NT4XMA
      ENDIF

C check NT4XO
      IF (NT4XO .NE. NT4XON .AND. NT4XO .NE. NT4XOY) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NT4XO'
         PRINT FMNII,VSTR,NT4XO
         PRINT FMBOO,VSTR
         CALL STOPGM('GFOURD','FATAL ERROR')
      ENDIF

C next line: temp coupling in 4 dims
      STR = 'NTT4, TEMP04, TAUT4'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTT4,PRGSTR,STR)
      CALL CHPREL(TEMP04,PRGSTR,STR)
      CALL CHPREL(TAUT4,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)


C check NTT4
      IF (NTT4 .NE. NTT4OF .AND. NTT4 .NE. NTT4ON) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTT4'
         PRINT FMNII,VSTR,NTT4
         PRINT FMBOO,VSTR
         CALL STOPGM('GFOURD','FATAL ERROR')
      ENDIF
      IF (NT4DIM .EQ. NT4OFF) THEN
         NTT4 = NTT4OF
      ENDIF

C     test TEMP04
      IF (TEMP04 .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'TEMP04'
         PRINT FMNIR,VSTR,TEMP04
         PRINT FMGER,VSTR
         CALL STOPGM('GFOURD','FATAL ERROR')
      ENDIF

C     test TAUT4
      IF (TAUT4 .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'TAUT4'
         PRINT FMNIR,VSTR,TAUT4
         PRINT FMGER,VSTR
         CALL STOPGM('GFOURD','FATAL ERROR')
      ENDIF

C     next line: 4D force scaling
      STR = 'NTCW4D, CW4DB, TEMP0B'
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTCW4D,PRGSTR,STR)
      CALL CHPREL(CW4DB,PRGSTR,STR)
      CALL CHPREL(TEMP0B,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)

C     check NTCW4D
      IF (NTCW4D .NE. NTC4NO .AND. NTCW4D .NE. NTC4YE)THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTCW4D'
         PRINT FMNII,VSTR,NTCW4D
         PRINT FMBOO,VSTR
         CALL STOPGM('GFOURD','FATAL ERROR')
      ENDIF

C     check CW4DB
      IF (CW4DB .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'CW4DB'
         PRINT FMNIR,VSTR,CW4DB
         PRINT FMGER,VSTR
         CALL STOPGM('GFOURD','FATAL ERROR')
      ENDIF

C     test TEMP0B
      IF (TEMP0B .LT. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'TEMP0B'
         PRINT FMNIR,VSTR,TEMP0B
         PRINT FMGER,VSTR
         CALL STOPGM('GFOURD','FATAL ERROR')
      ENDIF

C     next line: the interaction on/off switches
      STR = '4D INTERACTION SWITCHES'
      CALL GIMME(IUNIT,PRGSTR,STR)
      DO II=1,MX4DIN
         CALL CHPINT(IDO4D(II),PRGSTR,STR)
      ENDDO
      CALL CHKMTY(PRGSTR,STR)

      CALL GTEND(IUNIT,PRGSTR)

C     check all the flags, and assign the vals to NDO4D
      DO II=1,MX4DIN
         IF (IDO4D(II) .EQ. ITRUE) THEN
            NDO4D(II) = NDIM
         ELSEIF (IDO4D(II) .EQ. IFALSE) THEN
            NDO4D(II) = NDRMAX
         ELSE
            CALL FLAGLN(PRGSTR)
            PRINT FMNII,'4D interaction switch',IDO4D(II)
            PRINT FMBOO,'NTF4'
            CALL STOPGM('GFOURD','FATAL ERROR')
         ENDIF
      ENDDO

C     check for illegal combinations
      IF (NTCW4D .EQ. NTC4YE .AND. IDO4D(N4DHAR) .EQ. IFALSE) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,'nonsensical combination of 4D flags!'
         PRINT *,'scaling of 4D force constant ON: NTCW4D =',NTC4YE
         PRINT *,'but 4D force turned OFF!'
         CALL STOPGM('GFOURD','FATAL ERROR')
      ENDIF
C end GFOURD
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PFOURD
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
C common blocks
C     make sure this common block is changed
C     consistently in this file
      CHARACTER*(8)STR4D
      COMMON /PFOUR/STR4D(MX4DIN)
C local vars
      INTEGER I
C begin
      IF (LGOT(NFOURD)) THEN
         PRINT *,BTNAME(NFOURD)
         PRINT '(6A10)','NT4DIM','CW4DA','TEMP4I','NDFMI4',
     $        'NT4XI','NT4XO'
         PRINT '(I10,1X,E9.4,F10.3,3I10)',NT4DIM,CW4DA,TEMP4I,NDFMI4,
     $        NT4XI,NT4XO
         PRINT '(4A10)','NTT4','TEMP04','TAUT4'
         PRINT '(I10,3F10.3)',NTT4,TEMP04,TAUT4
         PRINT '(2A10)','NTCW4D','CW4DB'
         PRINT '(I10,1X,E9.4)', NTCW4D, CW4DB
         PRINT '(20A10)',(STR4D(I),I=1,MX4DIN)
         PRINT '(20I10)',(NDO4D(I)-NDRMAX,I=1,MX4DIN)

         PRINT 500,NDIM

         IF (NT4DIM .NE. NT4OFF) THEN
C     check implementation restrictions
C     no pressure coupling in 4D
            IF (NTP .NE. NTPOFF) THEN
               PRINT *,'Pressure coupling not possible in 4D'
               PRINT *,'NTP    =',NTP
               PRINT *,'NT4DIM =',NT4DIM
            ENDIF

C     analyse NT4XI for positions
            IF (NT4XI .EQ. NT4X0V) THEN
               PRINT 560,'SET TO ZERO'
            ELSEIF (NT4XI .EQ. NT4XBM) THEN
               PRINT 560,'TAKEN FROM A BOLTZMANN DISTRIBUTION'
            ELSEIF (NT4XI .EQ. NT4XRM .OR. NT4XI .EQ. NT4XRR) THEN
               PRINT 560,'READ WITH INITIAL 3D CONFIGURATION'
            ELSE
               PRINT *,'UNKOWN NT4XI FLAG = ',NT4XI
               CALL STOPGM('PFOURD','FATAL ERROR')
            ENDIF

C     analyse NT4XI for velocities
            IF (NT4XI .EQ. NT4X0V .OR. NT4XI .EQ. NT4XBM .OR.
     $           NT4XI .EQ. NT4XRM) THEN
               PRINT 510
            ELSEIF (NT4XI .EQ. NT4XRR) THEN
               PRINT *,'   INITIAL 4th D VELOCITIES',
     $              ' READ WITH INITIAL 3D CONFIGURATION'
            ELSE
               PRINT *,'UNKOWN NT4XI FLAG = ',NT4XI
               CALL STOPGM('PFOURD','FATAL ERROR')
            ENDIF

C     analyse NT4XO
            IF (NT4DIM .NE. NT4OFF .AND. NT4XO .EQ. NT4XOY) THEN
               PRINT *,'   TRAJECTORY FILES WRITTEN IN FOUR DIMENSIONS'
            ENDIF
C     analyse NTT4
            IF (NTT4 .EQ. NTT4ON) THEN
               PRINT *,'   TEMPERATURE COUPLING IN FOURTH DIMENSION'
            ENDIF
            
C     analyse NTKC4D
            IF (NTCW4D .EQ. NTC4YE) THEN
               PRINT 550,CW4DA,CW4DB
            ENDIF
         ENDIF
      ELSE
         PRINT 500,NDIM
      ENDIF

      PRINT *

 500  FORMAT (4X,'SIMULATION PERFORMED IN',I2,' DIMENSIONS')
 510  FORMAT (4X,'INITIAL 4th D VELOCITIES TAKEN FROM A',
     .     ' MAXWELL DISTRIBUTION')
 550  FORMAT (4X,'SCALING OF 4th D HARM. FORCE CONST. FROM',
     .     E8.2,' TO',E8.2)
 560  FORMAT (4X,'INITIAL 4th D POSITIONS ',A)
C end pfourd
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GBOUND(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'md.h'
      INCLUDE 'box.h'
      INCLUDE 'formats.h'
C local vars
      CHARACTER STR*(36)
      CHARACTER VSTR*(6)
C data
      DATA STR /'NTB,BOX(1),BOX(2),BOX(3),BETA,NRDBOX'/
C begin
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTB,PRGSTR,STR)
      CALL CHPREL(BOX(1),PRGSTR,STR)
      CALL CHPREL(BOX(2),PRGSTR,STR)
      CALL CHPREL(BOX(3),PRGSTR,STR)
      CALL CHPREL(BETA,PRGSTR,STR)
      CALL CHPINT(NRDBOX,PRGSTR,STR)

      CALL CHKMTY(PRGSTR,STR)

      CALL GTEND(IUNIT,PRGSTR)

C     check range of NTB
      IF (NTB .LT. NTBMIN .OR. NTB .GT. NTBMAX) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTB'
         PRINT FMNII,VSTR,NTB
         PRINT FMBTI,NTBMIN,VSTR,NTBMAX
         CALL STOPGM('GBOUND','FATAL ERROR')
      ENDIF

      IF (NRDBOX .NE. NRDBXN .AND. NRDBOX .NE. NRDBXY) THEN
         VSTR = 'NRDBOX'
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,VSTR,NRDBOX
         PRINT FMBOO
      ENDIF

C only check the box if we are NOT reading it from startup file,
C I.e. we ARE reading it from MD input.
      IF (NRDBOX .EQ. NRDBXN) THEN
         CALL CHKBOX(PRGSTR,.TRUE.)
      ENDIF
C end gbound
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE CHKBOX(PRGSTR,LFLAG)
C     ==--------------------------------------------------------------==
      INCLUDE 'coordsz.h'
      INCLUDE 'md.h'
      INCLUDE 'box.h'
      INCLUDE 'formats.h'
C args
      CHARACTER PRGSTR*(*)
      LOGICAL LFLAG
C local vars
      INTEGER K
      CHARACTER VSTR*(6)
C begin
      IF (NTB .EQ. NTBVAC) THEN
C     check 'no periodicity' case

      ELSE
C NTB is not 0 --> we have periodicity
C     check range of BOX vals
         DO K=1,NDIM
            IF (BOX(K) .LE. 0.0D0) THEN
               IF (LFLAG) THEN
                  CALL FLAGLN(PRGSTR)
               ENDIF
               PRINT *,'illegal value for box dims'
               PRINT *,'BOX(',K,') = ',BOX(K)
               PRINT FMGTR,'BOX(K)'
               CALL STOPGM('CHKBOX','FATAL ERROR')
            ENDIF
         ENDDO
C if we have a truncated octahedron, the BOX vals must be equal
         IF (NTB .LT. 0) THEN
            IF (BOX(1) .NE. BOX(2) .OR. BOX(1) .NE. BOX(3)) THEN
               IF (LFLAG) THEN
                  CALL FLAGLN(PRGSTR)
               ENDIF
               PRINT *,'Inconsistent values for NTB and BOX!'
               PRINT *,'Truncated octahedron boundary conds (NTB <0)'
               PRINT *,'require a cubic box BOX(1) = BOX(2) = BOX(3)'
               CALL STOPGM('CHKBOX','FATAL ERROR')
            ENDIF
         ENDIF

C     calc BOXH and BOXINV
!$OMP parallel do private(K)
         DO K=1,NDIM
            BOXH(K) = BOX(K)*0.5d+0
            BOXINV(K) = 1.0d+0/BOX(K)
         ENDDO

C     check range of BETA
         VSTR = 'BETA'
         IF (BETA .LE. 0.0D0) THEN
            IF (LFLAG) THEN
               CALL FLAGLN(PRGSTR)
            ENDIF
            PRINT FMNIR,VSTR,BETA
            PRINT FMGTR,VSTR
            CALL STOPGM('CHKBOX','FATAL ERROR')
         ENDIF

         IF (MOD(BETA, 1.8D+2) .LE. EPS) THEN
            IF (LFLAG) THEN
               CALL FLAGLN(PRGSTR)
            ENDIF
            PRINT FMNIR,VSTR,BETA
            PRINT *,'You have specified a flat box!'
            CALL STOPGM('CHKBOX','FATAL ERROR')
         ENDIF

C     check virial calculation case with PBCs
         IF (ABS(NTB) .EQ. NTBVIR) THEN
C     BETA must be 90.0
            IF (BETA .NE. 9.0D+1) THEN
               IF (LFLAG) THEN
                  CALL FLAGLN(PRGSTR)
               ENDIF
               PRINT FMNIR,VSTR,BETA
               PRINT 501
               PRINT 510
               PRINT 520,NTB
               CALL STOPGM('CHKBOX','FATAL ERROR')
            ENDIF
         ENDIF

         IF (NTB .LT. 0) THEN
C     check truncated octahedron case
C     BETA must be 90
            IF (BETA .NE. 9.0D+1) THEN
               IF (LFLAG) THEN
                  CALL FLAGLN(PRGSTR)
               ENDIF 
               PRINT FMNIR,VSTR,BETA
               PRINT 500
               PRINT 510
               PRINT 520,NTB
               CALL STOPGM('CHKBOX','FATAL ERROR')
            ENDIF

C      box must be cubic, not rectangular
            IF ( BOX(1) .NE. BOX(2) .OR.
     $           BOX(2) .NE. BOX(3)) THEN
               IF (LFLAG) THEN
                  CALL FLAGLN(PRGSTR)
               ENDIF
               PRINT 500
               PRINT *,'This only works if a cube,',
     $              ' (not a rectangular prism) is specified.'
               PRINT 520,NTB
               PRINT '(A9,I6)','BOX(1) =',BOX(1)
               PRINT '(A9,I6)','BOX(2) =',BOX(2)
               PRINT '(A9,I6)','BOX(3) =',BOX(3)
               PRINT *,'These must all be the same for a cube'
               CALL STOPGM('CHKBOX','FATAL ERROR')
            ENDIF
         ELSE
C check rectangular (monoclinic) case

C no specific checks here...

C end check rectangular case
         ENDIF
C endif ntb .eq. 0
      ENDIF

 500  FORMAT (1X,'For NTB < 0 periodicity is applied with a ',
     $     'truncated octahedron.')
 501  FORMAT (1X,'For ABS(NTB) = 2 the virial is calculated.')
 510  FORMAT (1X,'BETA must be 90')
 520  FORMAT (1X,'NTB =',I4)
C enc chkbox
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PBOUND
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
      INCLUDE 'mdblock.h'
      INCLUDE 'md.h'
C
      CHARACTER PSTR*(33)
      DATA PSTR/'   PERIODICITY IS APPLIED WITH A '/
C begin
      PRINT *,BTNAME(NBOUND)
      PRINT '(6A10)','NTB','BOX(1)','BOX(2)','BOX(3)','BETA','NRDBOX'
      PRINT '(I10,3F10.6,F10.3,I10)',
     $     NTB,BOX(1),BOX(2),BOX(3),BETA,NRDBOX
      IF (NTB .GT. 0) THEN
         IF (DABS(BETA-9.0D+1) .GT. 1.0D-1) THEN
            PRINT *,PSTR,'MONOCLINIC BOX'
         ELSE
            PRINT *,PSTR,'RECTANGULAR BOX'
         ENDIF
      ELSEIF (NTB .EQ. 0) THEN
         PRINT *,'   NO PERIODICITY IS APPLIED'
      ELSEIF (NTB .LT. 0) THEN
         PRINT *,PSTR,'TRUNCATED OCTAHEDRON'
      ENDIF
      IF (ABS(NTB) .EQ. NTBVIR) THEN
         PRINT *,'   THE VIRIAL IS CALCULATED'
      ENDIF
C box
      IF (NTB .NE. NTBVAC .AND. NRDBOX .EQ. NRDBXY) THEN
         PRINT *,'   BOX DIMENSIONS ARE READ FROM COORDINATE FILE'
         PRINT *,'   THIS OVERWRITES THE BOX DIMS SPECIFIED IN',
     $        ' THE ''BOUNDARY'' BLOCK'
      ENDIF

      PRINT *
C end pbound
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GPI(IUNIT,PRGSTR)
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
C local vars
      CHARACTER STR*(12)
      CHARACTER VSTR*(6)
C data
      DATA STR /'NTPI'/
C begin
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPINT(NTPI,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
      CALL GTEND(IUNIT,PRGSTR)

C do not test NTT for NTTPI, but test for NPID
      IF (NTPI .EQ. NTPION .AND. NPID .EQ. 0) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTPI'
         PRINT FMNII,VSTR,NTPI
         VSTR = 'NPID'
         PRINT FMNII,VSTR,NPID
         PRINT *,'You cannot use path integral without beads!'
         CALL STOPGM('GPI','FATAL ERROR')
      ENDIF
         
C test for 4D
      IF (NTPI .EQ. NTPION .AND. NT4DIM .NE. NT4OFF) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTPI'
         PRINT FMNII,VSTR,NTPI
         VSTR = 'NT4DIM'
         PRINT FMNII,VSTR,NT4DIM
         PRINT *,'You cannot use path integral in 4d!'
         CALL STOPGM('GPI','FATAL ERROR')
      ENDIF

C test for perturbation
      IF (NTPI .EQ. NTPION .AND. NTG .NE. NTGOFF) THEN
         CALL FLAGLN(PRGSTR)
         VSTR = 'NTPI'
         PRINT FMNII,VSTR,NTPI
         VSTR = 'NTG'
         PRINT FMNII,VSTR,NTG
         PRINT *,'Path integral: perturbation not yet full there!'
         CALL STOPGM('GPI','FATAL ERROR')
      ENDIF

C test for temperature
C      IF (NTPI .EQ. NTPION .AND. TPI .NE. TEMP0) THEN
C         CALL FLAGLN(PRGSTR)
C         VSTR = 'TPI'
C         PRINT FMNIR,VSTR,TPI
C         VSTR = 'TEMP0'
C         PRINT FMNIR,VSTR,TEMP0
C         PRINT *,'Path integral warning:'
C         PRINT *,'Temperature of MD .NE. temp. of discretization'
C      ENDIF

      LGOT(NPIBLK) = (NTPI .EQ. NTPION)
C end GPI
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PPI
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'md.h'
      INCLUDE 'mdblock.h'

C begin
      IF (LGOT(NPIBLK)) THEN
         PRINT *,BTNAME(NPIBLK)
         PRINT '(A10)','NTPI'
         PRINT '(I10)',NTPI

C warning if the wrong heat bath is chosen
         IF (NTPI .EQ. NTPION .AND. NTT(NTTPCM) .EQ. NTTOFF) THEN
            PRINT *,'Warning! Center of mass heat bath is off while'
            PRINT *,'path integral is on.'
         ENDIF
         PRINT *
      ENDIF
C end ppi
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      INCLUDE 'inc_rdmd5.h'
C     ==================================================================
COMMSUBR PRMD
C      SUBROUTINE PRMD
C
C     PRMD prints out the annotated PROMD control file.
C
COMMEND
      SUBROUTINE PRMD
C     ==--------------------------------------------------------------==
C begin
      CALL PTITLE
      CALL PEM
      CALL PSD
      CALL PSYST
      CALL PSTART
      CALL PSTEP
      CALL PBOUND
      CALL PNSP
      CALL PTCOUP
      CALL PPCOUP
      CALL PCOM
      CALL PPRINT
      CALL PWRITE
      CALL PSHAKE
      CALL PFORCE
      CALL PPLIST
      CALL PLRAN
      CALL PPREST
      CALL PDREST
      CALL PDHRES
      CALL PJVAL
      CALL PLE
      CALL PFOURD
      CALL PPERT
      CALL PPI
      INCLUDE 'inc_rdmd4.h'
C add any other blocks here!

C end PRMD
C     ==--------------------------------------------------------------==
      RETURN
      END
