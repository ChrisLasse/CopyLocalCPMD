C Id: blockio.f,v 1.51 1997/04/04 08:56:11 wscott Exp
Cmb Revised 15 Sept. 2005
C GROMOS data block input/output routines.
C
C     There is a read and a write routine for each
C     block (see table below).
C
C     The routines in this file are organised along the
C     following lines:
C     a) low level IO routines
C        with names starting with 'DMP' for writing
C        and 'GET' for reading.
C        These routines read/write exactly one block
C        type as defined in 'dataid.h'.
C     b) higher level IO routines which can build on
C        the low level routines.
C        Names start with WR for writing, and RD for
C        reading.
C
C     Only the b) routines are 'exported', i.e. only
C     they should be called from routines outside of the
C     file. See e.g. the routines in posio.f for examples
C     of how to do this.
C
C     Currently, there are 4 class a) IO routine pairs:
C     DMPRED/GETRED  3D vector (binary or reduced formatted)
C     DMPARR/GETARR  3D vector in (verbose) formatted
C     There are two low level 3D vector routines with special behaviour,
C     GETRE1 and GETAR1. They have slightly different semantics 
C     than GETRED and GETARR.
C     DMP4RE/GET4RE  4D component of vector (binary or reduced formatted)
C     DMP4AR/GET4AR  4D component of vector in (verbose) formatted
C
C     Each GET routine reads a block written by the corresponding
C     DMP routine.
C
C     Higher level routines such as WRCFG call the low level routines
C     such as DMPARR and DMP4AR depending on NDIM.
C
C     For all the routines that write positions, velocities etc.
C     in formatted form and not in reduced form (e.g. a POSITION
C     block with LFORM = .TRUE.), a topology file and a PROMD control
C     file must previously have been read in in order to be able
C     to write the atom names etc. Reading such blocks
C     does not require this, as the atom names are ignored
C     on reading coordinate blocks.
C
C     The blocks are read/written by the following routines
C     in this file:
C
C     block name         write    read      comments
C     FORMAT             WRFMT    CHCKFR    only written to binary files
C     TITLE              RDTIT    WRTIT     use PRTIT for printing titles
C     BOX                WRBOX    RDBOX
C
C     POSITION           DMPARR   GETARR   |see NOTES 1 and 2 below.
C     POSITION4THD       DMP4AR   GET4AR   |
C NOTE 1:
C     For writing a POSITION and POSITION4THD block with one call,
C     use WRXPOS.
C     For reading these blocks in one call, use either RDXPOS or RDXPO1.
C
C NOTE 2:
C     There are two special routines for reading
C     POSITION blocks used by analysis programs: GETSPE and GETAR1
C     There is one special routines for writing POSITION blocks
C     used by analysis routines: DMPSPE
C
C     POSITIONRED        DMPRED   GETRED   |see NOTE 3 below
C     POSITION4THDRED    DMP4RE   GET4RE   |for these together
C
C NOTE 3:
C     For writing a  POSITIONRED and a POSITION4THDRED block with one call,
C     use WRXRED.
C     For reading these block in one call, use either RDXRED or RDXRE1.
C
C     REFPOSITION        WRXC     RDXC
C     POSITIONOF         WRXOF    RDXOF
C
C     VELOCITY           DMPARR   GETARR   |use WRXVEL and RDXVEL
C     VELOCITY4THD       DMP4AR   GET4AR   |for these together
C
C     VELOCITYRED        DMPRED   GETRED   |use WRVRED and RDVRED
C     VELOCITY4THDRED    DMP4RE   GET4RE   |for these together
C
C     STOCHINT           WRSX     RDSX     |use WRSX and RDSX
C     STOCHINT4THD       DMP4AR   GET4AR   |for these together
C
C     PERTDATA           WRPER    RDPER
C     DISRESEXPAVE       WRIIAV   RDIIAV
C     JVALRESEXPAVE      WRJVAV   RDJVAV
C     LEMEMORY           WRLEME   RDLEME
C     TIMESTEP           WRTIME   RDTIME
C
C     BFACTOR            WRABF    RDABF
C     BFACTORANISO       WRANBF   RDANBF
C Note: there are two special routines for writing
C     atomic B factors used by analysis routines.
C     WRSABF for writing a BFACTOR block and
C     WRSANB for writing a BFACTORANISO block
C
C     POSITIONSECONDM    WRP2M    RDP2M
C     POSITIONTHIRDM     WRP3M    RDP3M
C     POSITIONFOURTHM    WRP4M    RDP4M
C     POSITIONSECONDMT   WRP2MT   RDP2MT
C
C     QUANTITYAVER       WRAVQ    RDAVQ
C     QUANENERAVER       WRAVQ    RDAVQ
C     QUANSUMENERAVER    WRAVQ    RDAVQ
C
C     QUANTIMESERIES     WRTCF    RDTCF
C     QUANDISTRIB        WRTCF    RDTCF
C     QUANTIMECORR       WRTCF    RDTCF
C     QUANTIMECORRSPE    WRTCF    RDTCF
C
C 'SPEC' blocks are written by humans and only
C  have to be read in formatted form.
C  The exception is the POSRESSPEC block which is
C  used in PROCMT.
C     FRICTIONSPEC       none     RDGARR
C     FOURDIMATOMSPEC    none     RD4DX
C     POSRESSPEC         WRNDXC   RDNDXC
C
C     DIHRESSPEC         none     RDDHRS
C     DISRESSPEC         none     RDDRST
C     JVALRESSPEC        none     RDJVSP
C     LEDIHSPEC          none     RDLESP
C
C     ENERGY             WRNRG    RDNRG
C     VOLUMEPRESSURE     WRVPRT   RDVPRT
C     FREEENERGYLAMBDA   WRRLAM   RDRLAM
C     FREEENERGY3D4      WRRMU    RDRMU
C
C     SOLVSTAT           WRSVS    RDSVS
C     DIFFSTAT           WRDFS    RDDFS
C     DIPMSTAT           WRDMS    RDDMS
C-----------------------------------------------------


COMMSUBR NDXID
C     INTEGER FUNCTION NDXID(BTN)
C     CHARACTER *(*)BTN
C
C     Given a block name in string form in BTN the index
C     defined in dataid.h is returned.
C     If the string is unknown, L<IDUKN> is returned.
COMMEND
      INTEGER FUNCTION NDXID(BTN)
C args
      INCLUDE 'dataid.h'
      INTEGER I
      CHARACTER *(*)BTN
C begin
C a quick heuristic: if the string starts with a blank,
C it is an unknown.
      IF (BTN(1:1) .EQ. ' ') THEN
         NDXID = IDUKN
         RETURN
      ENDIF
C otherwise look through the table
      I = MINIDT
C look for BTN in table using while loop
 10   IF (I .LE. MAXIDT) THEN
         IF(NAMEID(I) .NE. BTN) THEN
            I = I + 1
            GOTO 10
         ENDIF
      ENDIF
C     end of while loop
C
      IF (I .GT. MAXIDT) THEN
         NDXID = IDUKN
      ELSE
         NDXID = I
      ENDIF
C end NDXID
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR GNDXID
C      SUBROUTINE GNDXID(IUNIT,PRGSTR,STR,NDX)
C
C     Calls L<GIMME> to get the next line
C     and then calls L<NDXID> to get the predefined
C     index which is returned in NDX
C     
C     This only works for formatted files!!
COMMEND

      SUBROUTINE GNDXID(IUNIT,PRGSTR,STR,NDX)
      INCLUDE 'fileio.h'
C
      INTEGER IUNIT,NDX
      CHARACTER PRGSTR*(*)
      CHARACTER STR*(*)
C externals
      INTEGER NDXID
      EXTERNAL NDXID
C begin
      CALL GIMME(IUNIT,PRGSTR,STR)
      NDX = NDXID(FIOLIN)
C end GNDXID
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR GTEND
C     SUBROUTINE GTEND(IUNIT,PRGSTR)
C
C     GTEND reads in a line using L<RDBHDR> and checks whether it
C     contains the END marker.
C     If not, L<FLAGLN> is called, an error is written out and STOP
C     is invoked.
C
C     Only used on reading formatted versions of blocks.
COMMEND
      SUBROUTINE GTEND(IUNIT,PRGSTR)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
      LOGICAL LFORM
      PARAMETER (LFORM = .TRUE.)
C local vars
      INTEGER ID,NDUMMY
C begin
      CALL RDBHDR(IUNIT,LFORM,ID,NDUMMY)
      IF (ID .NE. IENDID) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,ENDEXP
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF
C end gtend
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR  CHNREC
C     SUBROUTINE CHNREC(PRGSTR,NREC,NREXP)
C
C     CHMREC compares NREC with NREXP. If they are not the same
C     an error message is written out and STOP is invoked.
C     This routine is used to check the consistency of data blocks:
C     The number of FORTRAN records (I.e. the number of WRITE statements
C     used in writing the block) is the same as the one expected.
COMMEND

      SUBROUTINE CHNREC(PRGSTR,NREC,NREXP)
      INCLUDE 'dataid.h'
      INTEGER NREC,NREXP
      CHARACTER PRGSTR*(*)
C begin
      IF (NREC .NE. NREXP) THEN
         PRINT *,PRGSTR,': invalid value for NREC!'
         PRINT EXPFIL,NREXP
         PRINT BONFIL,NREC
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF
C end CHNREC
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR PRTIT
C     SUBROUTINE PRTIT(HSTR,NTIT,TITLE)
C
C     PRTIT prints the HSTR to one line 
C     (this is a header line describing which title follows)
C     and then NTIT lines contained in TITLE to standard
C     output.
C
C     Note that this routine CANNOT be used to write the
C     title to a file. Use L<WRTIT> for that.
C     It is used to print a title to standard output on reading
C     a title block from a file.
COMMEND
      SUBROUTINE PRTIT(HSTR,NTIT,TITLE)
C args
      CHARACTER HSTR*(*)
      INTEGER NTIT
      CHARACTER *(*) TITLE
      DIMENSION TITLE(NTIT)
C local vars
      INTEGER I,ILEN
C begin
      PRINT 500,HSTR
      DO 5 I=1,NTIT
C avoid trailing blanks
         ILEN = LEN(TITLE(I))
C decrement until we find the first nonblack char
 10      IF (TITLE(I)(ILEN:ILEN) .EQ. ' ' .AND. ILEN .GT. 1) THEN
            ILEN = ILEN -1
            GOTO 10
         ENDIF
         PRINT 500,TITLE(I)(1:ILEN)
 5    CONTINUE
      PRINT *
 500  FORMAT (1X,A)
C end prtit
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRTIT
C     SUBROUTINE WRTIT(IUNIT,LFORM,NLINE,TI)
C
C     WRTIT writes a string TI into a TITLE BLOCK
C     to the open file IUNIT. If LFORM = .TRUE.,
C     the string is written formatted, otherwise it
C     is written in binary form.
C     The title block can be read by L<RDTIT> and printed
C     to screen by L<PRTIT>.
COMMEND

      SUBROUTINE WRTIT(IUNIT,LFORM,NLINE,TI)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
C args
      INTEGER IUNIT,NLINE
      LOGICAL LFORM
      CHARACTER *(*) TI
      DIMENSION TI(NLINE)
C local param statements
      INTEGER MYID
      PARAMETER (MYID = ITITID)
      INTEGER NREC
      PARAMETER (NREC=2)
C local vars
      INTEGER IRES,I,NLEN
      CHARACTER PRGSTR*(5)
      DATA PRGSTR /'WRTIT'/
C begin
      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(MYID))
         DO I=1,NLINE
            FIOLIN = TI(I)
            CALL WRSTR(IUNIT,PRGSTR,FIOLIN)
         ENDDO
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
C unformatted output
         WRITE(UNIT=IUNIT, IOSTAT=IRES) NAMEID(MYID),NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

         NLEN = LEN(TI(1))

         WRITE(UNIT=IUNIT, IOSTAT=IRES) NLINE,NLEN
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'NLINE',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF


         WRITE(UNIT=IUNIT, IOSTAT=IRES)(TI(I),I=1,NLINE)
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'title string',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end wrtit
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDTIT
C     SUBROUTINE RDTIT(IUNIT,LFORM,NREC,MAXLIN,TITLE,NLRED)
C
C     Subroutine RDTIT reads a title block from
C     the open file IUNIT. If LFORM = .TRUE., a
C     formatted read is attempted, otherwise an
C     unformatted read is attempted.
C     The block may not contain more than MAXLIN lines.
C     The actual number of lines read are returned in NLRED.
C
C     The title can subsequently be printed to standard
C     output by using L<PRTIT>. A title block can
C     be writen to file by using L<WRTIT>.
COMMEND

      SUBROUTINE RDTIT(IUNIT,LFORM,NREC,MAXLIN,TITLE,NLRED)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
C args
      INTEGER IUNIT,NREC,MAXLIN,NLRED
      LOGICAL LFORM
      CHARACTER*(*) TITLE
      DIMENSION TITLE(MAXLIN)
C local vars
      INTEGER IRES,NLEN,I,NDX
      CHARACTER PRGSTR*(5)
      CHARACTER STR*(5)
      DATA PRGSTR /'RDTIT'/
      DATA STR    /'TITLE'/
C begin
      IF (LFORM) THEN
C get the title
         NLRED = 0
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
 10      IF (NDX .EQ. IDUKN .AND. NLRED .LT. MAXLIN) THEN
            NLRED = NLRED + 1
            TITLE(NLRED) = FIOLIN

            CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
            GOTO 10
         ENDIF
C check for the end marker
         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,ENDEXP
            PRINT *,'Maximum number of lines in TITLE block: ',MAXLIN
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ELSE
C unformatted input
         CALL CHNREC(PRGSTR,NREC,2)

         READ(UNIT = IUNIT, IOSTAT = IRES) NLRED,NLEN
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'NLRED,NLEN',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

         IF (NLRED .GT. MAXLIN) THEN
            PRINT 500,PRGSTR,IUNIT
            PRINT MAXVAL,MAXLIN
            PRINT BONFIL,NLRED
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         IF (LEN(TITLE(1)) .LT. NLEN) THEN
            PRINT 500,PRGSTR,IUNIT
            PRINT *,'string too short !'
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

         READ(UNIT = IUNIT, IOSTAT = IRES) (TITLE(I),I=1,NLRED)
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'TITLE string',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
 500  FORMAT (1X,A,': too many title lines encountered on unit ',I5)
C end rdtit
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR DMPARR
C     SUBROUTINE DMPARR(IUNIT,LFORM,PRGSTR,IDBLK,NAT,NPM,NSM,NDIM,X)
C
C     DMPARR is a generic coordinate dumping routine.
C     It writes the first three coords of every atom in the system
C     to file using data from the topology and promd input file
C     to annotate the coordinates with atom and residue name etc.
C     These verbose blocks only exist in formatted files.
C     If LFORM = .FALSE., L<DMPRED> is called.
C     This routine only writes 3D coords to file.
C     If NDIM > 3, then L<DMP4AR> must be called in addition, so that
C     higher dimensional coords are written to file as well.
C
C     This routine writes atom name information to file which
C     it gets from the topology (include files toposz.h and topoar.h).
C     Thus, a topology file must have been read in before calling this
C     routine.
COMMEND

      SUBROUTINE DMPARR(IUNIT,LFORM,PRGSTR,IDBLK,NAT,NPM,NSM,NDIM,X)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
C args
      INTEGER IUNIT,NAT,IDBLK,NDIM,NPM,NSM
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
      real*8 X(*)
C externals
      LOGICAL LTOPOK
      EXTERNAL LTOPOK
C local vars
      INTEGER I,J,IL,INA,IRES
      CHARACTER SOLCRD*(24)
      DATA SOLCRD /'solute atom coordinates'/
C begin
C check we have read the topology
      IF (.NOT. LTOPOK()) THEN
         PRINT *,PRGSTR,
     $      ': need a topology to write verbose coordinate format!'
         PRINT *
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF
C check the number of atoms is consistent
      IF (NRP*NPM + NSM*NRAM .NE. NAT) THEN
         PRINT *,PRGSTR,':inconsistent number of atoms!'
         PRINT *,'on topology/controlfile:',NRP*NPM,' solute atoms and',
     $      NSM*NRAM,' solvent atoms.'
         PRINT *,'total: ',NRP*NPM + NSM*NRAM
         PRINT *,'but  wanted to write NAT = ',NAT,' atoms'
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

      IF (.NOT. LFORM) THEN
         CALL DMPRED(IUNIT,LFORM,PRGSTR,IDBLK,NAT,NDIM,X)
         RETURN
      ENDIF

C LFORM is TRUE ==> write verbose block
      CALL WRSTR(IUNIT,PRGSTR,NAMEID(IDBLK))
C write a little comment
      WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=FMIGNO)
      IF (IRES .NE. 0)THEN
         PRINT FWFAIL,PRGSTR,SOLCRD,IUNIT
         PRINT *,FMIGNO
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF
C     write solute(s)
      IL = 1
      INA = 1
      DO 10 I=1, NPM
         DO 20 J=1,NRP
           IF (IAC(J).NE.0) THEN
             WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=12)
     $            MRES(J),AANM(MRES(J)),PANM(J),INA,
     $            X(IL),X(IL+1),X(IL+2)
             IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,SOLCRD,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
             ENDIF
           ENDIF
           IL = IL + NDIM
           INA = INA + 1
 20      CONTINUE
 10   CONTINUE

C write solvent atoms
      DO 30 I=1,NSM
         DO 40 J=1,NRAM
           IF (IACS(J).NE.0) THEN
             WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=12)
     $            I,'SOLV ',ANMS(J),INA,X(IL),X(IL+1),X(IL+2)
             IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,'solvent atom coordinates',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
             ENDIF
           ENDIF
           IL = IL + NDIM
           INA = INA + 1
 40      CONTINUE
 30   CONTINUE
C     write end marker
      CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
C
 12   FORMAT(I5,2(1X,A5),I7,3F15.9)
C end DMPARR
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR DMP4AR
C      SUBROUTINE DMP4AR(IUNIT,LFORM,PRGSTR,ID4BLK,NAT,NPM,NSM,NDIM,X)
C
C     Write the 4th coordinate in a coordinates array to file.
C     
C     A topology file must have been read in before calling this
C     routine.
COMMEND
      SUBROUTINE DMP4AR(IUNIT,LFORM,PRGSTR,ID4BLK,NAT,NPM,NSM,NDIM,X)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
C args
      INTEGER IUNIT,NAT,ID4BLK,NDIM,NPM,NSM
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
      real*8 X(NDIM*NAT)
C externals
      LOGICAL LTOPOK
      EXTERNAL LTOPOK
C local vars
      INTEGER I,J,IL,INA,IRES
      CHARACTER SOLAB*(20)
      DATA SOLAB /'solute atom block'/
C begin
C check we have read the topology
      IF (.NOT. LTOPOK()) THEN
         PRINT *,PRGSTR,
     $      ': need a topology to write in verbose format!'
         PRINT *
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF
C check the number of atoms is consistent
      IF (NRP*NPM + NSM*NRAM .NE. NAT) THEN
         PRINT *,PRGSTR,':inconsistent number of atoms!'
         PRINT *,'on topology/controlfile:',NRP*NPM,' solute atoms and',
     $      NSM*NRAM,' solvent atoms.'
         PRINT *,'total: ',NRP*NPM + NSM*NRAM
         PRINT *,'but  wanted to write NAT = ',NAT,' atoms'
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

      IF (.NOT. LFORM) THEN
         CALL DMP4RE(IUNIT,LFORM,PRGSTR,ID4BLK,NAT,NDIM,X)
         RETURN
      ENDIF


C LFORM is TRUE ==> write verbose block
      CALL WRSTR(IUNIT,PRGSTR,NAMEID(ID4BLK))
      WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=FMIGNO)
      IF (IRES .NE. 0)THEN
         PRINT FWFAIL,PRGSTR,SOLAB,IUNIT
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF
C     write solute(s)
      IL = 1
      INA = 1
      DO 10 I=1, NPM
         DO 20 J=1,NRP
            WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=12)
     $           MRES(J),AANM(MRES(J)),PANM(J),INA,X(IL+3)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,SOLAB,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            IL = IL + NDIM
            INA = INA + 1
 20      CONTINUE
 10   CONTINUE

C write solvent atoms
      DO 30 I=1,NSM
         DO 40 J=1,NRAM
            WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=12)
     $           I,'SOLV ',ANMS(J),INA,X(IL+3)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,'solvent atom block',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            IL = IL + NDIM
            INA = INA + 1
 40      CONTINUE
 30   CONTINUE
C     write end marker
      CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
C
 12   FORMAT(I5,2(1X,A5),I7,3F15.9)
C end DMP4AR
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR DMP4RE
C     SUBROUTINE DMP4RE(IUNIT,LFORM,PRGSTR,IDBLK,NAT,NDIM,X)
C
C     Write 4D coords in reduced format to file.
COMMEND
      SUBROUTINE DMP4RE(IUNIT,LFORM,PRGSTR,IDBLK,NAT,NDIM,X)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,IDBLK,NAT,NDIM
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
      real*8 X(NDIM*NAT)
C local vars
      INTEGER NREC,IRES,I,I3
C begin
      IF (LFORM) THEN
         NREC = NAT+1
C
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IDBLK))
C write the whole array
         DO 10 I=1,NAT
            WRITE (UNIT=IUNIT, IOSTAT=IRES,FMT='(20F12.4)')
     $           (X(NDIM*(I-1)+I3),I3=4,NDIM)

            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,'coordinate data',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
            IF (MOD(I,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,I)
            ENDIF
 10      CONTINUE
C
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
         NREC = 2
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAMEID(IDBLK), NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAT
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'NAT',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES)
     $        ((X(NDIM*I+I3),I3=4,NDIM),I=0,NAT-1)
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'4th D coordinate data',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end DMP4RE
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR GETARR
C     SUBROUTINE GETARR(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
C
C     GETARR is a generic coordinate block reading routine.
C     NAT three dimensional coordinate sets are read into array
C     X.
C     NOTE: NDIM is only used to step through the array.
C     If NDIM > 3, then still only three coords/atom are
C     read in.
C     Coordinates of any higher dimensions are always in a
C     separate block in GROMOS.
C     This routine ignores the first 24 character on each line.
C     It is used to read the more "verbose" versions of the
C     coordinates blocks, e.g. POSITION and VELOCITY blocks,
C     but NOT POSITIONRED blocks.
C
C     Alternative routines for reading this type of blocks
C     are L<GETAR1> and L<GETSPE>.
COMMEND

      SUBROUTINE GETARR(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NREC,NAT,NDIM
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
      real*8 X(NDIM*NAT)
C externals
      LOGICAL LISMTY
      EXTERNAL LISMTY
C local vars
      INTEGER I3,NATTMP,NDX
      CHARACTER STR*(16)
C begin
      IF (LFORM) THEN
C formatted read: read until we get an END marker or
C we have the right number of atoms
         NATTMP = 0
         I3 = 1
         STR = 'coordinate data'
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C     while no marker  AND still atoms to read
 10      IF (NDX .EQ. IDUKN .AND. NATTMP .LT. NAT) THEN
            NATTMP = NATTMP + 1
            CALL SKPCHR(24)
            IF (LISMTY()) THEN
               CALL FLAGLN(PRGSTR)
               PRINT MTYLIN,24
               PRINT STREXP,STR
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            CALL CHPNRE(3,X(I3),'X coord',PRGSTR)
C            print *, ndim, i3,x(i3),x(i3+1),x(i3+2)
            I3 = I3 + NDIM

            CALL CHKMTY(PRGSTR,STR)

            CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C end while
            GOTO 10
         ENDIF

         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 500,NAT
            PRINT *,ENDEXP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ELSEIF (NATTMP .NE. NAT) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 500,NAT
            PRINT BONFIL,NATTMP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ELSE
         CALL GETRED(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
      ENDIF

 500  FORMAT (1X,'attempting to read',I7,' atom coordinates')
C end GETARR
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR GETAR1
C     SUBROUTINE GETAR1(IUNIT,LFORM,NREC,PRGSTR,MAXNAT,NDIM,NATRED,X)
C
C     GETAR1 is a generic coordinate block reading routine.
C     NAT three dimensional coordinates are read into array
C     X.
C     NOTE: NDIM is only used to step through the array.
C     If NDIM > 3, then still only three coords/atom are
C     read in.
C     Coordinates of any higher dimensions are always in a
C     separate block in GROMOS.
C     This routine ignores the first 24 character on each line.
C     It is used to read the more "verbose" versions of the
C     coordinates blocks, e.g. POSITION and VELOCITY blocks,
C     but NOT POSITIONRED blocks.
C
C     GETAR1 has a slightly different behaviour than L<GETARR>:
C     The routine is called with MAXNAT definining the size of the
C     array (which is MAXNAT*NDOM). The routine returns the number
C     of atoms read from file in NATRED. This is not the same
C     as L<GETARR> where the caller states exactly how many
C     atoms are expected to be on file.
COMMEND

      SUBROUTINE GETAR1(IUNIT,LFORM,NREC,PRGSTR,MAXNAT,NDIM,NATRED,X)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NREC,MAXNAT,NDIM,NATRED
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
      real*8 X(NDIM*MAXNAT)
C externals
      LOGICAL LISMTY
      EXTERNAL LISMTY
C local vars
      INTEGER I3,NDX
      CHARACTER STR*(16)
C begin
      IF (LFORM) THEN
C formatted read: read until we get an END marker or
C we have the right number of atoms
         NATRED = 0
         I3 = 1
         STR = 'coordinate data'
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C     while no marker  AND still atoms to read
 10      IF (NDX .EQ. IDUKN .AND. NATRED .LT. MAXNAT) THEN
            NATRED = NATRED + 1
            CALL SKPCHR(24)
            IF (LISMTY()) THEN
               CALL FLAGLN(PRGSTR)
               PRINT MTYLIN,24
               PRINT STREXP,STR
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF

            CALL CHPNRE(3,X(I3),'X coord',PRGSTR)
            CALL CHKMTY(PRGSTR,STR)

            I3 = I3 + NDIM

            CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C end while
            GOTO 10
         ENDIF

         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,ENDEXP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ELSEIF (NATRED .GT. MAXNAT) THEN
            CALL FLAGLN(PRGSTR)
            PRINT MAXVAL,MAXNAT
            PRINT BONFIL,NATRED
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ELSE
         CALL GETRE1(IUNIT,LFORM,NREC,PRGSTR,MAXNAT,NDIM,NATRED,X)
      ENDIF
C end GETAR1
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR DMPMXN
C     SUBROUTINE DMPMXN(IUNIT,LFORM,MYID,PRGSTR,MPAT,NAT,NPM,NSM,XX)
C     
C     DMPMXN ("dump M times N") is a generic writing routine
C     for an MPATxNAT array of reals of dimension XX(MPAT,NAT).
C     The block can be read by calling L<GETMXN>.
C
C     In the formatted case, the routine ensures that the
C     individual lines do not exceed 80 characters in length.
C     In order to do this, not more than 6 reals are written
C     to one line if MPAT > 6, with the remaining values per atom
C     wrapped to subsequent lines.
COMMEND

      SUBROUTINE DMPMXN(IUNIT,LFORM,MYID,PRGSTR,MPAT,NAT,NPM,NSM,XX)
C     includes
      INCLUDE 'dataid.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
C     local param statements
      INTEGER NREC
      PARAMETER (NREC=2)
C     the maximum number of reals on one line
      INTEGER MAXREL
      PARAMETER (MAXREL = 6)
C     args
      INTEGER IUNIT,MYID,MPAT,NAT,NPM,NSM
      LOGICAL LFORM
      real*8 XX(MPAT,NAT)
      CHARACTER PRGSTR*(*)
C externals
      LOGICAL LTOPOK
      EXTERNAL LTOPOK
C     local vars
      INTEGER I,J,INA,IRES,K
      INTEGER NWRITE,NWROTE
C     begin
C     check we have read the topology
      IF (.NOT. LTOPOK()) THEN
         PRINT *,PRGSTR,
     $        ': need a topology to write verbose coordinate format!'
         PRINT *
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF
C     check the number of atoms is consistent
      IF (NRP*NPM + NSM*NRAM .NE. NAT) THEN
         PRINT *,PRGSTR,':inconsistent number of atoms!'
         PRINT *,'on topology/controlfile:',NRP*NPM,' solute atoms and',
     $        NSM*NRAM,' solvent atoms.'
         PRINT *,'total: ',NRP*NPM + NSM*NRAM
         PRINT *,'but  wanted to write NAT = ',NAT,' atoms'
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(MYID))
C     write a little comment
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=FMIGNO)
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'solute data',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

C     write solute(s)
         INA = 1
         DO 10 I=1, NPM
            DO 20 J=1,NRP
               IF (MPAT .GT. MAXREL) THEN
                  NWRITE = MAXREL
               ELSE
                  NWRITE = MPAT
               ENDIF
               WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=12)
     $              MRES(J),AANM(MRES(J)),PANM(J),INA,
     $              (XX(K,INA),K=1,NWRITE)
               IF (IRES .NE. 0) THEN
                  PRINT FWFAIL,PRGSTR,'solute data',IUNIT
                  CALL STOPGM('BLOCKIO','FATAL ERROR')
               ENDIF
               NWROTE = NWRITE
C     while loop
 25            IF (NWROTE .LT. MPAT) THEN
                  NWRITE = MPAT - NWROTE
                  IF (NWRITE .GT. MAXREL) THEN
                     NWRITE = MAXREL
                  ENDIF
                  WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=13)
     $                 (XX(K,INA),K=NWROTE+1,NWROTE+NWRITE)
                  IF (IRES .NE. 0) THEN
                     PRINT FWFAIL,PRGSTR,'solute data',IUNIT
                     CALL STOPGM('BLOCKIO','FATAL ERROR')
                  ENDIF
                  NWROTE = NWROTE + NWRITE
                  GOTO 25
               ENDIF
               INA = INA + 1
 20         CONTINUE
 10      CONTINUE
         
C     write solvent atoms
         DO 30 I=1,NSM
            DO 40 J=1,NRAM
               IF (MPAT .GT. MAXREL) THEN
                  NWRITE = MAXREL
               ELSE
                  NWRITE = MPAT
               ENDIF
               WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=12)
     $              I,'SOLV ',ANMS(J),INA,
     $              (XX(K,INA),K=1,NWRITE)
               IF (IRES .NE. 0) THEN
                  PRINT FWFAIL,PRGSTR,'solvent data',IUNIT
                  CALL STOPGM('BLOCKIO','FATAL ERROR')
               ENDIF
               NWROTE = NWRITE
C     while loop
 45            IF (NWROTE .LT. MPAT) THEN
                  NWRITE = MPAT - NWROTE
                  IF (NWRITE .GT. MAXREL) THEN
                     NWRITE = MAXREL
                  ENDIF
                  WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=13)
     $                 (XX(K,INA),K=NWROTE+1,NWROTE+NWRITE)
                  IF (IRES .NE. 0) THEN
                     PRINT FWFAIL,PRGSTR,'solvent data',IUNIT
                     CALL STOPGM('BLOCKIO','FATAL ERROR')
                  ENDIF
                  NWROTE = NWROTE + NWRITE
                  GOTO 45
               ENDIF
               INA = INA + 1
 40         CONTINUE
 30      CONTINUE
C     write end marker
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
C     binary write
         WRITE(UNIT=IUNIT, IOSTAT=IRES) NAMEID(MYID),NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C     
         WRITE(UNIT=IUNIT, IOSTAT=IRES) MPAT,NAT
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'NAT',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C     
         WRITE(UNIT=IUNIT, IOSTAT=IRES)XX
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'XX',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF

C     12   FORMAT(I5,2(1X,A5),I7,50F9.5)
C removed for testing precision, not for the definitive version.
 12   FORMAT(I5,2(1X,A5),I7,6F9.5)
 13   FORMAT(24X,6F9.5)
C     end DMPMXN
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR GETMXN
C     SUBROUTINE GETMXN(IUNIT,LFORM,NREC,PRGSTR,MPAT,NAT,DXY)
C
C     GETMXN ("Get M times N") reads in MPAT reals per atom
C     specified by NAT into DXY.
C     The array DXY is assumed to be a (MPAT,NAT) array.
C     The block can be written to file by L<DMPMXN>.
COMMEND

      SUBROUTINE GETMXN(IUNIT,LFORM,NREC,PRGSTR,MPAT,NAT,DXY)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NAT,NREC,MPAT
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
      real*8 DXY(MPAT,NAT)
C externals
      LOGICAL LISMTY
      EXTERNAL LISMTY
C local vars
      INTEGER IRES,NATTMP,MTMP,NDX,K
      CHARACTER STR*(16)
C begin
      IF (LFORM) THEN
C formatted read: read until we get an END marker or
C we have the right number of atoms
         NATTMP = 0
         STR = 'data (reals)'
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C     while no marker  AND still atoms to read
 10      IF (NDX .EQ. IDUKN .AND. NATTMP .LT. NAT) THEN
            NATTMP = NATTMP + 1
            CALL SKPCHR(24)
            IF (LISMTY()) THEN
               CALL FLAGLN(PRGSTR)
               PRINT MTYLIN,24
               PRINT STREXP,STR
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF

            DO 80 K=1,MPAT
               CALL CHPREL(DXY(K,NATTMP),'DXY',PRGSTR)
 80         CONTINUE
            CALL CHKMTY(PRGSTR,STR)

            CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C end while
            GOTO 10
         ENDIF

         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 530,NAT
            PRINT *,ENDEXP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ELSEIF (NATTMP .NE. NAT) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 530,NAT
            PRINT BONFIL,NATTMP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ELSE
C unformatted reading
         CALL CHNREC(PRGSTR,NREC,2)
         READ(UNIT = IUNIT, IOSTAT = IRES) MTMP,NATTMP
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'NAT',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         IF (NATTMP .NE. NAT) THEN
            PRINT 500,'NAT',IUNIT
            PRINT EXPFIL,NAT
            PRINT BONFIL,NATTMP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         IF (MTMP .NE. MPAT) THEN
            PRINT 500,'MPAT',IUNIT
            PRINT EXPFIL,MPAT
            PRINT BONFIL,MTMP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         READ(UNIT = IUNIT, IOSTAT = IRES)DXY
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'data',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF

 500  FORMAT (1X,'inconsistency in ',A,' on UNIT',I5)
 530  FORMAT (1X,'attempting to read',I7,' atom data lines')
C end GETMXN
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR DMPSPE
C     SUBROUTINE DMPSPE(IUNIT,LFORM,PRGSTR,IDBLK,NAT,
C     $     MRES,AANM,PANM,INA,
C     $     NDIM,X)
C
C     DMPSPE will dump all the arrays MRES,AANM,PANM, INA
C     and X to a block the name of which is specified
C     by IDBLK. This is only implemented for the formatted case.
C     This can be used to write a POSITION block with any atom and
C     residue names. In contrast L<DMPARR> gets the naming information
C     from the topology.
C
C     See L<GETSPE> for a discussion of the special routines.
COMMEND

      SUBROUTINE DMPSPE(IUNIT,LFORM,PRGSTR,IDBLK,NAT,
     $     MRES,AANM,PANM,INA,
     $     NDIM,X)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
      INCLUDE 'toposz.h'
C args
      INTEGER IUNIT,NAT,NDIM,IDBLK
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
      INTEGER MRES(NAT),INA(NAT)
      CHARACTER*(MAXRLE) AANM(NAT)
      CHARACTER*(MAXNLE) PANM(NAT)
      real*8 X(NDIM*NAT)
C local vars
      INTEGER IRES,I,IL
C begin
      IF (.NOT. LFORM) THEN
         PRINT *,PRGSTR,': implementation restriction'
         PRINT *,'file must be formatted !'
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF
      CALL WRSTR(IUNIT,PRGSTR,NAMEID(IDBLK))
      IL = 1
      DO 10 I=1, NAT
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=12)
     $        MRES(I),AANM(I),PANM(I),INA(I),X(IL),X(IL+1),X(IL+2)
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,'solute atom coordinates',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         IL = IL + NDIM
 10   CONTINUE

C     write end marker
      CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
C
 12   FORMAT(I5,2(1X,A5),I7,3F15.9)
C end dmpspe
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR GETSPE
C     SUBROUTINE GETSPE(IUNIT,LFORM,PRGSTR,MAXNAT,NATRED,
C     $     MRES,AANM,PANM,INA,
C     $     NDIM,X)
C     GETSPE is a generic coordinate block reading routine
C     used by certain analysis programs in order to read
C     a POSITION block and also reading the annotations
C     ignored by GETARR. It only works on formatted
C     files.
C     I.e. on finding the following line in a POSITION block:
C     1 MEBMT N          1    3.173000000    2.997000000    2.913000000
C
C     L<GETARR> reads only the three coordinates, whereas GETSPE
C     attempts to read an int, two strings, an int and the
C     three coordinates. The actual format of the data on
C     the line is not relevant,i.e. the data is still read in
C     free format.
C
C     In addition to reading atomic data, NATRED is returned
C     with the number of atoms successfully read from the
C     POSITION block. This may not exceed MAXNAT.
C     
C     NOTE: NDIM is only used to step through the array.
C     If NDIM > 3, then still only three coords/atom are
C     read in.
COMMEND

      SUBROUTINE GETSPE(IUNIT,LFORM,PRGSTR,MAXNAT,NATRED,
     $     MRES,AANM,PANM,INA,
     $     NDIM,X)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'toposz.h'
C args
      INTEGER IUNIT,MAXNAT,NATRED,NDIM
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
      INTEGER MRES(MAXNAT),INA(MAXNAT)
      CHARACTER*(MAXRLE) AANM(MAXNAT)
      CHARACTER*(MAXNLE) PANM(MAXNAT)
      real*8 X(NDIM*MAXNAT)
C local vars
      INTEGER I3,NATTMP,NDX
      CHARACTER STR*(16)
C begin
      IF (LFORM) THEN
C formatted read: read until we get an END marker or
C we have the right number of atoms
         NATTMP = 0
         I3 = 1
         STR = 'coordinate data'
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C     while no marker  AND still atoms to read
 10      IF (NDX .EQ. IDUKN .AND. NATTMP .LT. MAXNAT) THEN
            NATTMP = NATTMP + 1

            CALL CHPINT(MRES(NATTMP),PRGSTR,'MRES')
            CALL CHPSTR(AANM(NATTMP),PRGSTR,'AANM')
            CALL CHPSTR(PANM(NATTMP),PRGSTR,'PANM')
            CALL CHPINT(INA(NATTMP),PRGSTR,'INA')

            CALL CHPNRE(3,X(I3),PRGSTR,'X')
            CALL CHKMTY(PRGSTR,STR)

            I3 = I3 + NDIM

            CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C end while
            GOTO 10
         ENDIF

         IF (NATTMP .GT. MAXNAT) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'too many atoms on file'
            PRINT MAXVAL,MAXNAT
            PRINT BONFIL,NATTMP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,ENDEXP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         NATRED = NATTMP
      ELSE
         PRINT *,'GETSPE: implementation restriction'
         PRINT *,'file must be formatted'
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF
C end GETSPE
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR DMPASP
C     SUBROUTINE DMPASP(IUNIT,LFORM,PRGSTR,IDBLK,MAXNAT,
C     $     MRES,AANM,PANM,INA,
C     $     NDIM,X)
C     DMPASP is a generic array dumping routine used by certain
C     analysis programs in order to write a block containing one real
C     value per atom.
C     It is used to write atomic B factors by L<PROCS2> for example.
C
C     Data in the MRES,AANM,PANM arrays are used for atom and residue
C     names instead of using the data in the topology.
C
C     If LFORM = .FALSE., the data is written as specified
C     in the binary isotropic B factor block.
COMMEND

      SUBROUTINE DMPASP(IUNIT,LFORM,PRGSTR,IDBLK,NATTOT,
     $     MRES,AANM,PANM,INA,VAL)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
      INCLUDE 'toposz.h'
C args
      INTEGER IUNIT,NATTOT,IDBLK
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
      INTEGER MRES(NATTOT),INA(NATTOT)
      CHARACTER*(MAXRLE) AANM(NATTOT)
      CHARACTER*(MAXNLE) PANM(NATTOT)
      real*8 VAL(NATTOT)
C local vars
      INTEGER I,IRES,NREC
C begin
      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IDBLK))
         DO 10 I=1, NATTOT
            WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=12)
     $           MRES(I),AANM(I),PANM(I),INA(I),VAL(I)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,'data',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
 10      CONTINUE

C     write end marker
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
         NREC = 2
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAMEID(IDBLK), NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NATTOT
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'NAT',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES)VAL
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'data',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF

 12   FORMAT(I5,2(1X,A5),I7,3F15.9)
C end DMPASP
C     ==--------------------------------------------------------------==
      RETURN
      END 

COMMSUBR GET4AR
C     SUBROUTINE GET4AR(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
C
C     GET4AR is a generic 4th dimension coordinate reading
C     routine. It is used by L<RDXPOS> and L<RDXVEL>.
COMMEND

      SUBROUTINE GET4AR(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NREC,NAT,NDIM
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
      real*8 X(NDIM*NAT)
C externals
      LOGICAL LISMTY
      EXTERNAL LISMTY
C local vars
      INTEGER I3,NDX,NATTMP
      CHARACTER STR*(24)
C begin
      IF (LFORM) THEN
C formatted read: read until we get an END marker or
C we have the right number of atoms
         NATTMP = 0
         I3 = NDIM
         STR = '4D coordinate data'
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C     while no marker  AND still atoms to read
 10      IF (NDX .EQ. IDUKN .AND. NATTMP .LT. NAT) THEN
            NATTMP = NATTMP + 1

            CALL SKPCHR(24)
            IF (LISMTY()) THEN
               CALL FLAGLN(PRGSTR)
               PRINT MTYLIN,24
               PRINT STREXP,STR
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF

            CALL CHPREL(X(I3),'W-coord',PRGSTR)
            CALL CHKMTY(PRGSTR,STR)

            I3 = I3 + NDIM

            CALL GNDXID(IUNIT,PRGSTR,STR,NDX)  
C end while
            GOTO 10
         ENDIF

         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 500,NAT
            PRINT *,ENDEXP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ELSEIF (NATTMP .NE. NAT) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 500,NAT
            PRINT BONFIL,NATTMP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ELSE
         CALL GET4RE(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
      ENDIF
C 12   FORMAT (24X,3F15.9)
 500  FORMAT (1X,'attempting to read',I7,' atom coordinates')
C end GET4AR
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR DMPRED
C     SUBROUTINE DMPRED(IUNIT,LFORM,PRGSTR,IDBLK,NAT,NDIM,X)
C
C     DMPRED is a generic array dumping subroutine.
C     It dumps a coordinate array to file and is used for
C     pos, vel, force,...
C     DMPRED writes one block containing coords in 3D only.
C     Any additional blocks containing 4D coords can be written
C     by calling L<DMP4RED>.
C     The data can be read by subroutine L<GETRED> or L<GETRE1>.
COMMEND

      SUBROUTINE DMPRED(IUNIT,LFORM,PRGSTR,IDBLK,NAT,NDIM,X)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,IDBLK,NAT,NDIM
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
      real*8 X(NDIM*NAT)

C local vars
      INTEGER NREC,IRES,I,I3
C begin
      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IDBLK))
         I3 = 1
         DO 10 I=1,NAT
            WRITE (UNIT=IUNIT, IOSTAT=IRES,FMT='(3F15.9)')
     $           X(I3),X(I3+1),X(I3+2)
            IF (IRES .NE. 0) THEN
               PRINT '(2A)',PRGSTR,' Error writing data!'
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
            IF (MOD(I,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,I)
            ENDIF

            I3 = I3 + NDIM
 10      CONTINUE
C
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
         NREC = 2
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAMEID(IDBLK), NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAT
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'NAT',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES)
     $        ((X(NDIM*I+I3),I3=1,3),I=0,NAT-1)
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'coordinate data',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C     end DMPRED
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR GETRED
C     SUBROUTINE GETRED(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
C
C     GETRED is a generic array reading subroutine
C     which reads data written by L<DMPRED>.
C     It is used for coordinate arrays such as for positions,
C     velocities, etc.
C
C     NAT: the number of atom coordinates to read.
C     If the number of atoms on the block does not match,
C     an error messages is printed to screen and STOP
C     is envoked. If this behaviour is not wanted,
C     then use L<GETRE1>.
COMMEND
      SUBROUTINE GETRED(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
C args
      INTEGER IUNIT,NAT,NDIM,NREC
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
      real*8 X(NDIM*NAT)
C local params
C as this is a low level routine, we try to optimize for speed.
C if LFAST = .TRUE. we use the faster version of reading in reals
C in the formatted case.
C If this doesnt work on certain architectures, 
C set LFAST to .FALSE.
      LOGICAL LFAST
      PARAMETER (LFAST = .TRUE.)
C local vars
      INTEGER IRES,I,I3
      INTEGER NATTMP,NDX,NUMRED
      CHARACTER STR*(16)
C begin
      IF (LFORM) THEN
C try "fast" reading
         IF (LFAST) THEN
            NATTMP = 0
            I3 = 1
            CALL GTNREL(IUNIT,3,X(I3),NUMRED)
 50         IF (NUMRED .EQ. 3 .AND. NATTMP .LT. NAT) THEN
               NATTMP = NATTMP + 1
               I3 = I3 + NDIM
               IF (NATTMP .LT. NAT) THEN
                  CALL GTNREL(IUNIT,3,X(I3),NUMRED)
               ENDIF
               GOTO 50
            ENDIF
C
            CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
            IF (NDX .NE. IENDID) THEN
               CALL FLAGLN(PRGSTR)
               PRINT 530,NAT
               PRINT *,ENDEXP
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ELSEIF (NATTMP .NE. NAT) THEN
               CALL FLAGLN(PRGSTR)
               PRINT 530,NAT
               PRINT BONFIL,NATTMP
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
         ELSE
C----
C safe but slow reading. If the methods above dont work on
C a certain machine, we can always use this one
C----
            NATTMP = 0
            I3 = 1
            STR = 'coordinate data'
            CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C     while no marker  AND still atoms to read
 10         IF (NDX .EQ. IDUKN .AND. NATTMP .LT. NAT) THEN
               NATTMP = NATTMP + 1

               CALL CHPNRE(3,X(I3),PRGSTR,STR)
               CALL CHKMTY(PRGSTR,STR)

               I3 = I3 + NDIM

               CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C end while
               GOTO 10
            ENDIF

            IF (NDX .NE. IENDID) THEN
               CALL FLAGLN(PRGSTR)
               PRINT 530,NAT
               PRINT *,ENDEXP
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ELSEIF (NATTMP .NE. NAT) THEN
               CALL FLAGLN(PRGSTR)
               PRINT 530,NAT
               PRINT BONFIL,NATTMP
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
         ENDIF
      ELSE
C unformatted reading
         CALL CHNREC(PRGSTR,NREC,2)
         READ(UNIT = IUNIT, IOSTAT = IRES) NATTMP
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'NAT',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         IF (NATTMP .NE. NAT) THEN
            PRINT 500,IUNIT
            PRINT EXPFIL,NAT
            PRINT BONFIL,NATTMP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         READ(UNIT = IUNIT, IOSTAT = IRES)
     $        ((X(NDIM*I+I3),I3=1,3),I=0,NAT-1)

         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'data',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF

 500  FORMAT (1X,'inconsistency in NAT on UNIT',I5)
 530  FORMAT (1X,'attempting to read',I7,' atom coordinates')
C end GETRED
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR GETRE1
C     SUBROUTINE GETRE1(IUNIT,LFORM,NREC,PRGSTR,MAXNAT,NDIM,NATRED,X)
C
C     GETRE1 is similar to L<GETRED> bu has slightly different behaviour.
C     GETRE1 reads a block containig coordinate data
C     (positions, velocities, etc.).
C     NATRED is delivered with the number of atoms
C        of which the coordinates have been read in.
C        If the number of atoms on file is larger than MAXNAT,
C        the routine prints an error message and calls STOP.
C
C     X(NDIM*MAXNAT) delivered with NATRED*NDIM coordinates.
C     NDIM is only used to step through the array. Only three
C     coordinates are read. Thus if NDIM is 3, then the NATRED*NDIM
C     values of X are read. If NDIM = 4, then every fourth element
C     in X is not modified.
C
C     MAXNAT the maximum number of atoms that can be read into X.
C     This is not modified by the routine.
COMMEND

      SUBROUTINE GETRE1(IUNIT,LFORM,NREC,PRGSTR,MAXNAT,NDIM,NATRED,X)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,MAXNAT,NDIM,NREC,NATRED
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
      real*8 X(NDIM*MAXNAT)
C local vars
      INTEGER IRES,I,I3
      INTEGER NDX
      CHARACTER STR*(16)
C begin
      IF (LFORM) THEN
         NATRED = 0
         I3 = 1
         STR = 'coordinate data'
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C     while no marker  AND still atoms to read
 10      IF (NDX .EQ. IDUKN .AND. NATRED .LT. MAXNAT) THEN
            NATRED = NATRED + 1

            CALL CHPNRE(3,X(I3),PRGSTR,STR)
            CALL CHKMTY(PRGSTR,STR)

            I3 = I3 + NDIM

            CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C end while
            GOTO 10
         ENDIF

         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,ENDEXP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ELSEIF (NATRED .GT. MAXNAT) THEN
            CALL FLAGLN(PRGSTR)
            PRINT MAXVAL,MAXNAT
            PRINT BONFIL,NATRED
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ELSE
C unformatted reading
         CALL CHNREC(PRGSTR,NREC,2)
         READ(UNIT = IUNIT, IOSTAT = IRES) NATRED
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'NAT',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         IF (NATRED .GT. MAXNAT) THEN
            PRINT MAXVAL,MAXNAT
            PRINT BONFIL,NATRED
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         READ(UNIT = IUNIT, IOSTAT = IRES)
     $        ((X(NDIM*I+I3),I3=1,3),I=0,NATRED-1)

         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'data',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end GETRE1
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR GET4RE
C     SUBROUTINE GET4RE(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
C
C     Reads in a 4th dimension coordinates block.
COMMEND
      SUBROUTINE GET4RE(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
C args
      INTEGER IUNIT,NAT,NDIM,NREC
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
      real*8 X(NDIM*NAT)
C local vars
      INTEGER IRES,I,I3,NDX,NATTMP
      CHARACTER STR*(16)
C begin
      IF (LFORM) THEN
         NATTMP = 0
         I3 = 1
         STR = 'pos or vel data'
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C     while no marker  AND still atoms to read
 10      IF (NDX .EQ. IDUKN .AND. NATTMP .LT. NAT) THEN
            NATTMP = NATTMP + 1
            CALL CHPREL(X(I3+3),PRGSTR,STR)
            CALL CHKMTY(PRGSTR,STR)

            I3 = I3 + NDIM

            CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C     end while
            GOTO 10
         ENDIF
         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 550,NAT
            PRINT *,ENDEXP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ELSEIF (NATTMP .NE. NAT) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 550,NAT
            PRINT BONFIL,NATTMP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ELSE
C unformatted reading
         CALL CHNREC(PRGSTR,NREC,2)
         READ(UNIT = IUNIT, IOSTAT = IRES) NATTMP
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'NAT',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         IF (NATTMP .NE. NAT) THEN
            PRINT 500,IUNIT
            PRINT EXPFIL,NAT
            PRINT BONFIL,NATTMP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         READ(UNIT = IUNIT, IOSTAT = IRES)
     $        ((X(NDIM*I+I3),I3=4,NDIM),I=0,NAT-1)
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'4th D coordinate data',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF

 500  FORMAT (1X,'inconsistency in NAT on UNIT',I5)
 550  FORMAT (1X,'attempting to read',I7,' atom coordinates')
C end GET4RE
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR GNDX
C     SUBROUTINE GNDX(IUNIT,LFORM,PRGSTR,
C     $     NATTOT,MAXNDX,NDXARR,NUMRED)
C     
C     GNDX is a generic atom sequence number reading routine.
C     It is called directly by L<RDNDXC> (for sequence numbers
C     of position re(con) strained atoms) for example.
C
C     GNDX reads a block from IUNIT containing atom sequence
C     numbers of atoms which must be [1..NATTOT].
C     A maximum of MAXNDX sequence numbers are read.
C
C     The format on the file is the same as the POSITION block
C     format. 
C     However, ONLY THE SEQUENCE NUMBER COLUMN IS SIGNIFICANT,
C     I.e. the first 17 characters before and any remaining data
C     after the one integer read in on each line are ignored.
C
C     NATTOT and MAXNDX are not changed.
C     PRGSTR a string used to identify the callin subroutine
C        in error messages written to screen.
C     NDXARR(MAXNDX)
C        delivered with the sequence numbers of the atoms.
C     NUMRED
C        delivered with the number of seq. numbers read in.
COMMEND

      SUBROUTINE GNDX(IUNIT,LFORM,PRGSTR,
     $     NATTOT,MAXNDX,NDXARR,NUMRED)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NATTOT,MAXNDX,NUMRED
      INTEGER NDXARR(MAXNDX)
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
C externals
      LOGICAL LISMTY
      EXTERNAL LISMTY
C local vars
      INTEGER NDUMMY,NLAST,NDX
      CHARACTER STR*(28)
C data
      DATA STR /'atom index data (INTEGER)'/
C begin
      IF (LFORM) THEN
         NUMRED = 0
         NLAST = 0
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C     while no marker AND still atoms to read
 10      IF (NDX .EQ. IDUKN .AND. NUMRED .LT. MAXNDX) THEN

            CALL SKPCHR(17)
            IF (LISMTY()) THEN
               CALL FLAGLN(PRGSTR)
               PRINT MTYLIN,17
               PRINT STREXP,STR
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            CALL CHPINT(NDUMMY,PRGSTR,STR)
            IF (NDUMMY .LT. 1 .OR. NDUMMY .GT. NATTOT) THEN
               CALL FLAGLN(PRGSTR)
               PRINT 501,NDUMMY,IUNIT
               PRINT 502,NATTOT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            IF (NDUMMY .LT. NLAST) THEN
               CALL FLAGLN(PRGSTR)
               PRINT *,'Atom sequence numbers must be in ascending',
     $              ' order on UNIT ',IUNIT
               PRINT *,'Sequence number',NDUMMY,
     $              ' must be AFTER number ',NLAST
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ELSEIF (NDUMMY .EQ. NLAST) THEN
               CALL FLAGLN(PRGSTR)
               PRINT *,'Sequence number',NDUMMY,
     $              ' specified twice on UNIT = ',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF

            NUMRED = NUMRED+1
            NDXARR(NUMRED) = NDUMMY
            NLAST = NDUMMY

            CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C end while
            GOTO 10
         ENDIF

         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 500,NATTOT
            PRINT *,ENDEXP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

         IF (NUMRED .GT. MAXNDX) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 500,NUMRED,IUNIT
            PRINT *,'Too many atom sequence numbers!'
            PRINT *,'maximum allowed: MAXNDX = ',MAXNDX
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ELSE
C     unformatted file no supported
         PRINT *,PRGSTR,
     $        ' UNFORMATTED ATOM SEQUENCE BLOCK NOT IMPLEMENTED!'
         PRINT *,'UNIT = ',IUNIT
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

 500  FORMAT (1X,'NUMNDX =',I6,' out of range on UNIT ',I2)
 501  FORMAT (1X,'atom seq. num = ',I5,'out of range on UNIT ',I2)
 502  FORMAT (1X,'MUST be 1 <= seq. num.  <= NAT = ',I5)
C end gndx
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDNDXC
C     SUBROUTINE RDNDXC(IUNIT,LFORM,NLIM,MAXXC,JRC,NRC)
C
C     Reads a block containing the atom position con/restraining
C     indices written by L<WRNDXC>
C     This only works form LFORM = .TRUE.
C     This routine calls L<GNDX>.
COMMEND
      SUBROUTINE RDNDXC(IUNIT,LFORM,NLIM,MAXXC,JRC,NRC)
C args
      INTEGER IUNIT,NLIM,MAXXC,NRC
      INTEGER JRC(MAXXC)
      LOGICAL LFORM
C begin
      CALL GNDX(IUNIT,LFORM,'RDNDXC',NLIM,MAXXC,JRC,NRC)
C end RDNDXC
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRNDXC
C      SUBROUTINE WRNDXC(IUNIT,LFORM,NRC,MRES,AANM,PANM,JRC)
C
C     Writes a list of atom position con/restraining
C     indices to file. Information about atom names
C     are taken from the arrays MRES,AANM and PANM.
C     The block can be read by L<RDNDXC>.
COMMEND
      SUBROUTINE WRNDXC(IUNIT,LFORM,NRC,MRES,AANM,PANM,JRC)
C includes
      INCLUDE 'toposz.h'
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NRC
      INTEGER MRES(NRC)
      CHARACTER *(MAXRLE)AANM(NRC)
      CHARACTER *(MAXNLE)PANM(NRC)

      INTEGER JRC(NRC)
      LOGICAL LFORM
C local vars
      INTEGER IRES,I
      CHARACTER PRGSTR*(6)
C data
      DATA PRGSTR /'WRNDXC'/
C begin
      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IXPRID))
         DO 10 I=1,NRC
            WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=12)
     $           MRES(I),AANM(I),PANM(I),JRC(I)
 10      CONTINUE
C
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
         PRINT *,PRGSTR,': can only write FORMATTED!'
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

 12   FORMAT(I5,2(1X,A5),I7,3F15.9)
C end wrndxc
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RD4NDX
C     SUBROUTINE RD4DX(IUNIT,LFORM,NREC,NATTOT,C4D)
C
C     RD4NDX reads a block from IUNIT specifying the
C     atoms in three or four dimensions.
COMMEND
      SUBROUTINE RD4NDX(IUNIT,LFORM,NREC,NATTOT,C4D)
C args
      INTEGER IUNIT,NATTOT,NREC
      real*8 C4D(NATTOT)
      LOGICAL LFORM
C begin
      CALL GETANM(IUNIT,LFORM,NREC,'RD4NDX',NATTOT,C4D)
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDGARR
C     SUBROUTINE RDGARR(IUNIT,LFORM,NREC,NATTOT,GAM)
C
C     RDGARR reads a block from IUNIT specifying the
C     atomic friction coeffients of atoms for SD.
COMMEND
      SUBROUTINE RDGARR(IUNIT,LFORM,NREC,NATTOT,GAM)
C args
      INTEGER IUNIT,NATTOT,NREC
      real*8 GAM(NATTOT)
      LOGICAL LFORM
C begin
      CALL GETANM(IUNIT,LFORM,NREC,'RDGARR',NATTOT,GAM)
C     ==--------------------------------------------------------------==
      RETURN
      END
      
COMMSUBR WRSABF
C     SUBROUTINE WRSABF(IUNIT,LFORM,NATTOT,
C     $     MRES,AANM,PANM,INA,VAL)
C
C     WRSABF writes an isotropic atomic B factor
C     block to file which. Also see L<WRABF>.
C     It calls L<DMPASP>.
COMMEND
      SUBROUTINE WRSABF(IUNIT,LFORM,NATTOT,
     $     MRES,AANM,PANM,INA,VAL)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'toposz.h'
C args
      INTEGER IUNIT,NATTOT
      LOGICAL LFORM
      INTEGER MRES(NATTOT),INA(NATTOT)
      CHARACTER*(MAXRLE) AANM(NATTOT)
      CHARACTER*(MAXNLE) PANM(NATTOT)
      real*8 VAL(NATTOT)
C begin
      CALL DMPASP(IUNIT,LFORM,'WRSABF',IABFID,NATTOT,
     $     MRES,AANM,PANM,INA,VAL)
C end WRSABF
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRABF
C     SUBROUTINE WRABF(IUNIT,LFORM,NATTOT,NPM,NSM,CXC)
C
C     Writes an atomic isotropic B-factor block to file.
C     The CRC array is written as is WITHOUT MODIFICATION --
C     the 8*pi**2/3 factor should already be in CRC.
C
C     The block can be read using L<RDABF>.
COMMEND
      SUBROUTINE WRABF(IUNIT,LFORM,NATTOT,NPM,NSM,CXC)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NPM,NSM,NATTOT
      LOGICAL LFORM
      real*8 CXC(NATTOT)
C begin
      CALL DMPMXN(IUNIT,LFORM,IABFID,'WRABF',
     $     1,NATTOT,NPM,NSM,CXC)
C end wrabf
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDABF
C     SUBROUTINE RDABF(IUNIT,LFORM,NREC,NATTOT,CXC)
C
C     RDABF reads in a block specifying the isotropic atomic
C     B factors for NATTOT atoms.
C     The values are returned AS THEY ARE on file without any
C     changes (see L<WRABF>).
COMMEND

      SUBROUTINE RDABF(IUNIT,LFORM,NREC,NATTOT,CXC)
C args
      INTEGER IUNIT,NATTOT,NREC
      real*8 CXC(NATTOT)
      LOGICAL LFORM
C local vars
      INTEGER I
C begin
      CALL GETMXN(IUNIT,LFORM,NREC,'RDABF',1,NATTOT,CXC)

C      CALL GETANM(IUNIT,LFORM,NREC,'RDABF',NATTOT,CXC)
C make sure we dont have negative Bfactors
C we shouldnt really even have zero values, but some
C cristallographers do...
      DO 10 I=1,NATTOT
         IF (CXC(I) .LT. 0.0D0) THEN
            PRINT 500,'RDABF',I,IUNIT
            PRINT 501,CXC(I)
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
 10   CONTINUE

 500  FORMAT (1X,A,': illegal atomic B factor',
     $     'value for atom',I5,' on unit',I3)
 501  FORMAT (1X,'value read:', F15.6)
C rdabf
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRSANB
C     SUBROUTINE WRSANB(IUNIT,LFORM,NATTOT,
C     $     MRES,AANM,PANM,INA,DXY)
C
C     WRSANB writes an anisotropic atomic B factor
C     block to file.
C     Instead of taking the annotations from the topology
C     as L<WRANBF> does, they are taken from the arrays passed as
C     arguments.
C     This is only implemented for formatted writing.
COMMEND
      SUBROUTINE WRSANB(IUNIT,LFORM,NATTOT,
     $     MRES,AANM,PANM,INA,DXY)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'toposz.h'
      INCLUDE 'coordsz.h'
C args
      INTEGER IUNIT,NATTOT
      LOGICAL LFORM
      INTEGER MRES(NATTOT),INA(NATTOT)
      CHARACTER*(MAXRLE) AANM(NATTOT)
      CHARACTER*(MAXNLE) PANM(NATTOT)
      real*8 DXY(M2MAX,NATTOT)
C local vars
      INTEGER I,IRES,K
      CHARACTER PRGSTR*(6)
      CHARACTER STRDAT*(20)
      DATA PRGSTR /'WRSANB'/
      DATA STRDAT /'aniso B factor data'/
C begin
      IF (.NOT. LFORM) THEN
         PRINT *,PRGSTR,': implementation restriction'
         PRINT *,'file must be formatted !'
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF
      CALL WRSTR(IUNIT,PRGSTR,NAMEID(IANBID))
C write a little comment
      WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=FMIGNO)
      IF (IRES .NE. 0)THEN
         PRINT FWFAIL,PRGSTR,STRDAT,IUNIT
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

      DO 10 I=1, NATTOT
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=12)
     $        MRES(I),AANM(I),PANM(I),INA(I),
     $        (DXY(K,I),K=1,M2MAX)
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,STRDAT,IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
 10   CONTINUE

C     write end marker
      CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
C
 12   FORMAT(I5,2(1X,A5),I7,6F9.5)
C end WRSABF
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDANBF
C     SUBROUTINE RDANBF(IUNIT,LFORM,NREC,NATTOT,DXY)
C
C     Reads in an anisotropic B factor block from
C     file.
COMMEND
      SUBROUTINE RDANBF(IUNIT,LFORM,NREC,NATTOT,DXY)
      INCLUDE 'coordsz.h'
C args
      INTEGER IUNIT,NATTOT,NREC
      LOGICAL LFORM
      real*8 DXY(M2MAX,NATTOT)
C begin
      CALL GETMXN(IUNIT,LFORM,NREC,'RDANBF',M2MAX,NATTOT,DXY)
C end rdanbf
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRANBF
C     SUBROUTINE WRANBF(IUNIT,LFORM,NATTOT,NPM,NSM,DXY)
C
C     Writes an atomic anisotropic B-factor block to file.
C     The 8*pi**2 factor should already be in DXY.
C     The block can be read using L<RDANBF>.
COMMEND
      SUBROUTINE WRANBF(IUNIT,LFORM,NATTOT,NPM,NSM,DXY)
      INCLUDE 'coordsz.h'
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NATTOT,NPM,NSM
      LOGICAL LFORM
      real*8 DXY(M2MAX,NATTOT)
C begin
      CALL DMPMXN(IUNIT,LFORM,IANBID,'WRANBF',M2MAX,NATTOT,
     $     NPM,NSM,DXY)
C end wranbf
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR GETANM
C     SUBROUTINE GETANM(IUNIT,LFORM,NREC,PRGSTR,NATTOT,C4D)
C
C     GETANM is a generic subroutine that reads in 
C     NATTOT REAL*8 values into the C4D array,
C     (.i.e. one REAL*8 value for every atom)
C     ignoring 24 spaces at the head of every line.
C     Only formatted reading is supported in this release.
C     It is used in subroutines such as L<RDABF>,
C     L<RDGARR> and L<RD4NDX>.
C
C     In the event of an error, PRGSTR is used in order
C     identify the calling subroutine and STOP is invoked.
COMMEND
      SUBROUTINE GETANM(IUNIT,LFORM,NREC,PRGSTR,NATTOT,C4D)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NATTOT,NREC
      real*8 C4D(NATTOT)
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
C externals
      LOGICAL LISMTY
      EXTERNAL LISMTY
C local vars
      INTEGER NDX,NATTMP,IRES
      CHARACTER STR*(24)
C data
      DATA STR /'atom data (REAL)'/
C begin
      IF (LFORM) THEN
C formatted read: read until we get an END marker or
C we have the right number of atoms
         NATTMP = 0
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C     while no marker AND still atoms to read
 10      IF (NDX .EQ. IDUKN .AND. NATTMP .LT. NATTOT) THEN
            NATTMP = NATTMP + 1

            CALL SKPCHR(24)
            IF (LISMTY()) THEN
               CALL FLAGLN(PRGSTR)
               PRINT MTYLIN,24
               PRINT STREXP,STR
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            CALL CHPREL(C4D(NATTMP),PRGSTR,STR)
            CALL CHKMTY(PRGSTR,STR)

            CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C end while
            GOTO 10
         ENDIF

         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN(PRGSTR)
            PRINT EXPFIL,NATTOT
            PRINT *,ENDEXP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ELSEIF (NATTMP .NE. NATTOT) THEN
            CALL FLAGLN(PRGSTR)
            PRINT EXPFIL,NATTOT
            PRINT BONFIL,NATTMP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ELSE
C     unformatted file 
         CALL CHNREC(PRGSTR,NREC,2)
         READ(UNIT = IUNIT, IOSTAT = IRES) NATTMP
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'NAT',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         IF (NATTMP .NE. NATTOT) THEN
            PRINT FLGBIN,PRGSTR,IUNIT
            PRINT EXPFIL,NATTOT
            PRINT BONFIL,NATTMP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         READ(UNIT = IUNIT, IOSTAT = IRES) C4D
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'data',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C     end getanm
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDBOX
C     SUBROUTINE RDBOX(IUNIT,LFORM,NREC,BOX)
C
C     Reads in a box block from file written by L<WRBOX>.
C     Only three coordinates are ever present on file.
COMMEND RDBOX
      SUBROUTINE RDBOX(IUNIT,LFORM,NREC,BOX)
      INCLUDE 'coordsz.h'
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NREC
      LOGICAL LFORM
      real*8 BOX(MAXDIM)
C local vars
      INTEGER I,IRES
      CHARACTER PRGSTR*(5)
      CHARACTER STR*(5)
      DATA PRGSTR /'RDBOX'/
C begin
      IF (LFORM) THEN
C get box
         STR = 'BOX'
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPREL(BOX(1),PRGSTR,STR)
         CALL CHPREL(BOX(2),PRGSTR,STR)
         CALL CHPREL(BOX(3),PRGSTR,STR)
C get end marker
         CALL GTEND(IUNIT,PRGSTR)
      ELSE
C we have a binary file
         CALL CHNREC(PRGSTR,NREC,1)
C box
         READ(UNIT = IUNIT, IOSTAT = IRES)(BOX(I),I=1,3)
         write(6,*)'grooomos',box
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'BOX',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end rdbox
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRBOX
C     SUBROUTINE WRBOX(IUNIT,LFORM,BOX)
C
C     Writes a box block to file that can be read by L<RDBOX>.
C     Only three coordinates are ever written.
COMMEND
      SUBROUTINE WRBOX(IUNIT,LFORM,BOX)
      INCLUDE 'coordsz.h'
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT
      LOGICAL LFORM
      real*8 BOX(MAXDIM)
C local vars
      INTEGER I,IRES,NREC
      CHARACTER PRGSTR*(5)
      DATA PRGSTR /'WRBOX'/
C begin
      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IBOXID))
C box
         WRITE(UNIT = IUNIT, IOSTAT = IRES,FMT=500)(BOX(I),I=1,3)
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'BOX',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C end
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
         NREC = 1
         WRITE(UNIT=IUNIT, IOSTAT=IRES) NAMEID(IBOXID),NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C box
         WRITE(UNIT = IUNIT, IOSTAT = IRES)(BOX(I),I=1,3)
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'BOX',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF

 500  FORMAT (3F15.9)
C end wrbox
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRPER
C     SUBROUTINE WRPER(IUNIT,LFORM,NDSZE,RLAM,DLSUM,RMUE,DMSUM)
C
C     Write a perturbation continuation block to file.
C     The block can be read by L<RDPER>.
COMMEND
      SUBROUTINE WRPER(IUNIT,LFORM,NDSZE,RLAM,DLSUM,RMUE,DMSUM)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NDSZE
      LOGICAL LFORM
      real*8 RLAM,DLSUM(NDSZE)
      real*8 RMUE,DMSUM(NDSZE)
C local vars
      INTEGER IRES,NREC,I
      CHARACTER PRGSTR*(6)
      CHARACTER VARSTR*(24)
      DATA PRGSTR /'WRPER'/
      DATA VARSTR /'RLAM,DLSUM,RMU,DMSUM'/
C begin
      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IPRTID))
C RLAM
         CALL WRCOMM(IUNIT,PRGSTR,'RLAM')
         WRITE(UNIT=IUNIT, IOSTAT=IRES, FMT=500)RLAM
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,'RLAM',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C DLSUM
         CALL WRCOMM(IUNIT,PRGSTR,'DLSUM')
         DO 10 I=1,NDSZE
            WRITE(UNIT=IUNIT, IOSTAT=IRES, FMT=500)DLSUM(I)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,'DLSUM',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            IF (MOD(I,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,I)
            ENDIF
 10      CONTINUE
C RMUE
         CALL WRCOMM(IUNIT,PRGSTR,'RMU')
         WRITE(UNIT=IUNIT, IOSTAT=IRES, FMT=500)RMUE
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,'RMU',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C DMSUM
         CALL WRCOMM(IUNIT,PRGSTR,'DMSUM')
         DO 20 I=1,NDSZE
            WRITE(UNIT=IUNIT, IOSTAT=IRES, FMT=500)DMSUM(I)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,'DMSUM',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            IF (MOD(I,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,I)
            ENDIF
 20      CONTINUE

C write end marker
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
         NREC = 2
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAMEID(IPRTID),NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES)RLAM,DLSUM
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,VARSTR,IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES)RMUE,DMSUM
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,VARSTR,IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF

 500  FORMAT (3E15.7)
C end WRPER
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDPER
C      SUBROUTINE RDPER(IUNIT,LFORM,NREC,NDSZE,RLAM,DLSUM,RMUE,DMSUM)
C
C     Read a perturbation continuation block written by L<WRPER>
C     from file.
COMMEND
      SUBROUTINE RDPER(IUNIT,LFORM,NREC,NDSZE,RLAM,DLSUM,RMUE,DMSUM)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NREC,NDSZE
      LOGICAL LFORM
      real*8 RLAM,DLSUM(NDSZE)
      real*8 RMUE,DMSUM(NDSZE)
C local vars
      INTEGER IRES,I
      CHARACTER PRGSTR*(6)
      CHARACTER STR*(24)
      DATA PRGSTR /'RDPER'/
      DATA STR /'RLAM,DLSUM,RMU,DMSUM'/
C begin
      IF (LFORM) THEN
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPREL(RLAM,PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)   
         DO I=1,NDSZE
            CALL GIMME(IUNIT,PRGSTR,STR)
            CALL CHPREL(DLSUM(I),PRGSTR,STR)
            CALL CHKMTY(PRGSTR,STR)
         ENDDO

         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPREL(RMUE,PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)
         DO I=1,NDSZE
            CALL GIMME(IUNIT,PRGSTR,STR)
            CALL CHPREL(DMSUM(I),PRGSTR,STR)
            CALL CHKMTY(PRGSTR,STR)
         ENDDO

C get the end marker
         CALL GTEND(IUNIT,PRGSTR)
      ELSE
         CALL CHNREC(PRGSTR,NREC,2)
C
         READ(UNIT = IUNIT, IOSTAT = IRES)RLAM,DLSUM
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,STR,IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         READ(UNIT = IUNIT, IOSTAT = IRES)RMUE,DMSUM
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,STR,IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end RDPER
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRNRG
C     SUBROUTINE WRNRG(IUNIT,LFORM,
C     $     MXENER,ENER,
C     $     MXERES,ENERES,
C     $     NUMUSD,ENERLJ,ENERCL,ENERRF,ENERRC)
C
C     Writes an energy block to file.
C     The block can be read by L<RDNRG>.
COMMEND
      SUBROUTINE WRNRG(IUNIT,LFORM,
     $     MXENER,ENER,
     $     MXERES,ENERES,
     $     NUMUSD,ENERLJ,ENERCL,ENERRF,ENERRC)
Cincludes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,MXENER,MXERES,NUMUSD
      LOGICAL LFORM
      real*8 ENER(MXENER),ENERES(MXERES)
      real*8 ENERLJ(NUMUSD*(NUMUSD+1)/2),ENERCL(NUMUSD*(NUMUSD+1)/2)
      real*8 ENERRF(NUMUSD*(NUMUSD+1)/2),ENERRC(NUMUSD*(NUMUSD+1)/2)
C local vars
      INTEGER IRES,I,NREC,NWRITE
      CHARACTER PRGSTR*(5)
      DATA PRGSTR /'WRNRG'/
C begin
      NWRITE = NUMUSD*(NUMUSD+1)/2
      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(INRGID))
C write ENER
         CALL WRCOMM(IUNIT,PRGSTR,'ENER')
         DO 10 I=1,MXENER
            WRITE(UNIT=IUNIT, IOSTAT=IRES, FMT=500)ENER(I)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,'ENER',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            IF (MOD(I,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,I)
            ENDIF
 10      CONTINUE
         CALL WRCOMM(IUNIT,PRGSTR,'ENERES')
C write ENERES
         DO 20 I=1,MXERES
            WRITE(UNIT=IUNIT, IOSTAT=IRES, FMT=500)ENERES(I)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,'ENERES',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
            IF (MOD(I,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,I)
            ENDIF
 20      CONTINUE
C write NUMUSD
         CALL WRCOMM(IUNIT,PRGSTR,'NUMUSD')
         WRITE(UNIT=IUNIT, IOSTAT=IRES, FMT=510)NUMUSD
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,'NUMUSD',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C write ENERLJ, ENERCL,...
         CALL WRCOMM(IUNIT,PRGSTR,'ENERLJ,ENERCL,ENERRF,ENERRC')
         DO 30 I=1,NWRITE
            WRITE(UNIT=IUNIT, IOSTAT=IRES, FMT=500)
     $           ENERLJ(I),ENERCL(I),ENERRF(I),ENERRC(I)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,'ENERLJ,..',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
 30      CONTINUE
C write end marker
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
         NREC = 3
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAMEID(INRGID), NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES) ENER,ENERES
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'ENER,ENERES',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NUMUSD
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'NUMUSD',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT=IUNIT, IOSTAT = IRES)
     $        (ENERLJ(I),I=1,NWRITE),
     $        (ENERCL(I),I=1,NWRITE),
     $        (ENERRF(I),I=1,NWRITE),
     $        (ENERRC(I),I=1,NWRITE)
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'ENERLJ,...',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF

 500  FORMAT (4E15.7)
 510  FORMAT (I5)
C end wrnrg
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDNRG
C     SUBROUTINE RDNRG(IUNIT,LFORM,NREC,
C     $     MXENER,ENER,
C     $     MXERES,ENERES,
C     $     MXGGRP,NUMUSD,ENERLJ,ENERCL,ENERRF,ENERRC)
C
C     Read an energy block from file.
C     MXGGRP is is the size of the array ENERLJ,ENERCL,
C     ENERRF,ENERRC.
C     NUMUSD is returned with the number of energy groups
C     read from file. Thus, the number of elements valid in the arrays
C     ENERLJ,ENERCL,ENERRF,ENERRC will be NUMUSD*(NUMUSD+1)/2.
C
C     If NUMUSD*(NUMUSD+1)/2 is larger than MXGGRP, an error
C     is written out and STOP is envoked.
COMMEND
      SUBROUTINE RDNRG(IUNIT,LFORM,NREC,
     $     MXENER,ENER,
     $     MXERES,ENERES,
     $     MXGGRP,NUMUSD,ENERLJ,ENERCL,ENERRF,ENERRC)
Cincludes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NREC,MXENER,MXERES,MXGGRP,NUMUSD
      LOGICAL LFORM
      real*8 ENER(MXENER),ENERES(MXERES)
      real*8 ENERLJ(MXGGRP),ENERCL(MXGGRP)
      real*8 ENERRF(MXGGRP),ENERRC(MXGGRP)
C local vars
      INTEGER IRES,I,NREAD
      CHARACTER PRGSTR*(5)
      CHARACTER STR*(6)
      DATA PRGSTR /'RDNRG'/
C begin
      IF (LFORM) THEN
C read ENER
         STR = 'ENER'
         DO I=1,MXENER
            CALL GIMME(IUNIT,PRGSTR,STR)
            CALL CHPREL(ENER(I),PRGSTR,STR)
            CALL CHKMTY(PRGSTR,STR)
         ENDDO
C read ENERES
         STR = 'ENERES'
         DO I=1,MXERES
            CALL GIMME(IUNIT,PRGSTR,STR)
            CALL CHPREL(ENERES(I),PRGSTR,STR)
            CALL CHKMTY(PRGSTR,STR)
         ENDDO

C read NUMUSD
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPINT(NUMUSD,PRGSTR,STR)
         NREAD = NUMUSD*(NUMUSD+1)/2
         IF (NREAD .GT. MXGGRP) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'NUMUSD too big!'
            PRINT 500,NUMUSD
            PRINT 505,NREAD
            PRINT 510,MXGGRP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C read ENERLJ,...
         DO I=1,NREAD
            CALL GIMME(IUNIT,PRGSTR,STR)
            CALL CHPREL(ENERLJ(I),PRGSTR,STR)
            CALL CHPREL(ENERCL(I),PRGSTR,STR)
            CALL CHPREL(ENERRF(I),PRGSTR,STR)
            CALL CHPREL(ENERRC(I),PRGSTR,STR)
            CALL CHKMTY(PRGSTR,'Energies: 4 reals')
         ENDDO

C get the end marker
         CALL GTEND(IUNIT,PRGSTR)
      ELSE
         CALL CHNREC(PRGSTR,NREC,3)
C
         READ(UNIT = IUNIT, IOSTAT = IRES) ENER,ENERES
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'ENER,ENERES',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         READ(UNIT = IUNIT, IOSTAT = IRES) NUMUSD
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'NUMUSD',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

         NREAD = NUMUSD*(NUMUSD+1)/2
         IF (NREAD .GT. MXGGRP) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'NUMUSD too big!'
            PRINT 500,NUMUSD
            PRINT 505,NREAD
            PRINT 510,MXGGRP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         READ(UNIT=IUNIT, IOSTAT=IRES)
     $        (ENERLJ(I),I=1,NREAD),
     $        (ENERCL(I),I=1,NREAD),
     $        (ENERRF(I),I=1,NREAD),
     $        (ENERRC(I),I=1,NREAD)
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'ENERLJ,..',IUNIT
            PRINT 500,NUMUSD
            PRINT 505,NREAD
            PRINT 510,MXGGRP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
 500  FORMAT (1X,'NUMUSD =',I5)
 505  FORMAT (1X,'NREAD  =',I5)
 510  FORMAT (1X,'MXGGRP =',I5)
C end rdnrg
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRVPRT
C     SUBROUTINE WRVPRT(IUNIT,LFORM,MXVTBL,VOLPRT)
C     
COMMEND
      SUBROUTINE WRVPRT(IUNIT,LFORM,MXVTBL,VOLPRT)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,MXVTBL
      LOGICAL LFORM
      real*8 VOLPRT(MXVTBL)
C local vars
      INTEGER IRES,I,NREC
      CHARACTER PRGSTR*(6)
C data
      DATA PRGSTR /'WRVPRT'/
C begin
      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IDVPRT))
         DO 10 I=1,MXVTBL
            WRITE(IUNIT,IOSTAT=IRES,FMT=22)VOLPRT(I)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,'VOLPRT',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
            IF (MOD(I,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,I)
            ENDIF
 10      CONTINUE
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
C binary write
         NREC = 1
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAMEID(IDVPRT), NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES) VOLPRT
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'VOLPRT',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
 22   FORMAT (E15.7)
C end WRVPRT
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRRLAM
C     SUBROUTINE WRRLAM(IUNIT,LFORM,MXENER,ENER,RLAM,MAXFR,FRENER)
C
COMMEND
      SUBROUTINE WRRLAM(IUNIT,LFORM,MXENER,ENER,RLAM,MAXFR,FRENER)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,MXENER,MAXFR
      LOGICAL LFORM
      real*8 ENER(MXENER),RLAM,FRENER(MAXFR)
C local vars
      INTEGER I,IRES,NREC
      CHARACTER PRGSTR*(6)
      DATA PRGSTR /'WRRLAM'/
C begin
      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IDRLAM))
C ener
         CALL WRCOMM(IUNIT,PRGSTR,'ENER')
         DO 10 I=1,MXENER
            WRITE(IUNIT,IOSTAT=IRES,FMT=22)ENER(I)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,'ENER',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
            IF (MOD(I,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,I)
            ENDIF
 10      CONTINUE
C rlam
         CALL WRCOMM(IUNIT,PRGSTR,'RLAM')
         WRITE(IUNIT,IOSTAT=IRES,FMT=22)RLAM
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,'RLAM',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C fren
         CALL WRCOMM(IUNIT,PRGSTR,'FRENER')
         DO 20 I=1,MAXFR
            WRITE(IUNIT,IOSTAT=IRES,FMT=22)FRENER(I)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,'FRENER',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
            IF (MOD(I,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,I)
            ENDIF
 20      CONTINUE
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
         NREC = 1
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAMEID(IDRLAM), NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES) ENER,RLAM,FRENER
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'ENER,RLAM,FRENER',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF

 22   FORMAT (E15.7)
C end wrrlam
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDRLAM
C      SUBROUTINE RDRLAM(IUNIT,LFORM,NREC,MXENER,ENER,RLAM,MAXFR,FRENER)
C
COMMEND
      SUBROUTINE RDRLAM(IUNIT,LFORM,NREC,MXENER,ENER,RLAM,MAXFR,FRENER)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,MXENER,MAXFR,NREC
      LOGICAL LFORM
      real*8 ENER(MXENER),RLAM,FRENER(MAXFR)
C local vars
      INTEGER I,IRES
      CHARACTER PRGSTR*(6)
      CHARACTER STR*(12)
      DATA PRGSTR /'RDRLAM'/
C begin
      IF (LFORM) THEN
C ener
         STR = 'ENER(I)'
         DO I=1,MXENER
            CALL GIMME(IUNIT,PRGSTR,STR)
            CALL CHPREL(ENER(I),PRGSTR,STR)
            CALL CHKMTY(PRGSTR,STR)
         ENDDO
C rmu
         STR = 'RLAM'
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPREL(RLAM,PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)
C fren
         STR = 'FRENER(I)'
         DO I=1,MAXFR
            CALL GIMME(IUNIT,PRGSTR,STR)
            CALL CHPREL(FRENER(I),PRGSTR,STR)
            CALL CHKMTY(PRGSTR,STR)
         ENDDO

C get the end marker
         CALL GTEND(IUNIT,PRGSTR)
      ELSE
         CALL CHNREC(PRGSTR,NREC,1)
C
         READ(UNIT = IUNIT, IOSTAT = IRES) ENER,RLAM,FRENER
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'ENER,RLAM,FRENER',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end rdrlam
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRRMU
C     SUBROUTINE WRRMU(IUNIT,LFORM,MXENER,ENER,RMU,MAXFR,FRENER)
C
COMMEND
      SUBROUTINE WRRMU(IUNIT,LFORM,MXENER,ENER,RMU,MAXFR,FRENER)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,MXENER,MAXFR
      LOGICAL LFORM
      real*8 ENER(MXENER),RMU,FRENER(MAXFR)
C local vars
      INTEGER I,IRES,NREC
      CHARACTER PRGSTR*(6)
      DATA PRGSTR /'WRRMU'/
C begin
      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IDRMU))
C ener
         CALL WRCOMM(IUNIT,PRGSTR,'ENER')
         DO 10 I=1,MXENER
            WRITE(IUNIT,IOSTAT=IRES,FMT=22)ENER(I)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,'ENER',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
            IF (MOD(I,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,I)
            ENDIF
 10      CONTINUE
C rmu
         CALL WRCOMM(IUNIT,PRGSTR,'RMU')
         WRITE(IUNIT,IOSTAT=IRES,FMT=22)RMU
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,'RMU',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C fren
         CALL WRCOMM(IUNIT,PRGSTR,'FRENER')
         DO 20 I=1,MAXFR
            WRITE(IUNIT,IOSTAT=IRES,FMT=22)FRENER(I)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,'FRENER',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
            IF (MOD(I,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,I)
            ENDIF
 20      CONTINUE
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
         NREC = 1
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAMEID(IDRMU), NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES) ENER,RMU,FRENER
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'ENER,RMU,FRENER',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF

 22   FORMAT (E15.7)
C end wrrmu
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDRMU
C     SUBROUTINE RDRMU(IUNIT,LFORM,NREC,MXENER,ENER,RMU,MAXFR,FRENER)
C
COMMEND
      SUBROUTINE RDRMU(IUNIT,LFORM,NREC,MXENER,ENER,RMU,MAXFR,FRENER)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,MXENER,MAXFR,NREC
      LOGICAL LFORM
      real*8 ENER(MXENER),RMU,FRENER(MAXFR)
C local vars
      INTEGER I,IRES
      CHARACTER PRGSTR*(6)
      CHARACTER STR*(12)
      DATA PRGSTR /'RDRMU'/
C begin
      IF (LFORM) THEN
C ener
         STR = 'ENER(I)'
         DO I=1,MXENER
            CALL GIMME(IUNIT,PRGSTR,STR)
            CALL CHPREL(ENER(I),PRGSTR,STR)
            CALL CHKMTY(PRGSTR,STR)
         ENDDO
C rmu
         STR = 'RMU'
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPREL(RMU,PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)
C fren
         STR = 'FRENER(I)'
         DO I=1,MAXFR
            CALL GIMME(IUNIT,PRGSTR,STR)
            CALL CHPREL(FRENER(I),PRGSTR,STR)
            CALL CHKMTY(PRGSTR,STR)
         ENDDO

C get the end marker
         CALL GTEND(IUNIT,PRGSTR)
      ELSE
         CALL CHNREC(PRGSTR,NREC,1)
C
         READ(UNIT = IUNIT, IOSTAT = IRES) ENER,RMU,FRENER
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'ENER,RMU,FRENER',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end rdrmu
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRFMT
C     SUBROUTINE WRFMT(IUNIT,LFORM)
C
C     WRFMT writes a special FORMAT block to
C     the open file IUNIT which
C     allows GROMOS to verify that the binary format
C     is compatible when reading the file at a 
C     later date. The checking on reading is performed
C     by L<CHKFMT>. Thus, attempting to read a file 
C     written in single precision, when GROMOS is compiled
C     in double precision will produce an error on reading.
C
C     WRFMT only writes a block if LFORM = .FALSE.
COMMEND
      SUBROUTINE WRFMT(IUNIT,LFORM)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT
      LOGICAL LFORM
C local params
      INTEGER NREC
      PARAMETER (NREC = 1)
C local vars
      CHARACTER PRGSTR*(5)
      INTEGER IRES
      DATA PRGSTR /'WRFMT'/
C begin
C we dont write anything if we have formatted output
      IF (.NOT. LFORM) THEN
C block header
         WRITE (UNIT = IUNIT, IOSTAT = IRES) NAMEID(IFRMID),NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C format info
         WRITE (UNIT = IUNIT, IOSTAT = IRES) RMAGIC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'RMAGIC',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end WRFMT
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR CHCKFR
C     SUBROUTINE CHCKFR(IUNIT)
C
C     CHCKFR checks whether the format of data contained in
C     the open file IUNIT is compatible with the compiled program.
C     I.e. if data in single precision is read from a binary
C     file in double precision, the an error is generated
C     and STOP is envoked.
COMMEND
      SUBROUTINE CHCKFR(IUNIT)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT
C local vars
      INTEGER IRES
      real*8 RTEST
C begin
      READ(IUNIT, IOSTAT = IRES) RTEST
      IF (IRES .NE. 0)THEN
         PRINT *,'CHCKFR: failed to read data record!'
      ENDIF

      IF (IRES .NE. 0 .OR. RTEST .NE. RMAGIC) THEN
         PRINT *,'CHCKFR: incompatible binary format on unit',IUNIT
         PRINT *,'(wrong precision or byte swapping problem)'
         PRINT *
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF
C end CHCKFR
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDBHDR
C     SUBROUTINE RDBHDR(IUNIT,LFORM,IDBLK,NREC)
C
C     RDBHDR attempts to read a header from the open
C     file IUNIT.
C     The file may be formatted (use LFORM=.TRUE.) or
C     binary (use LFORM = .FALSE.).
C     IDBLK is returned with a value as defined
C     in dataid.h.
C     If a reading error occurs, IDERR is returned.
C     If a header is read which is not in the
C     NAMEID array defined in 'dataid.h' then IDUKN is
C     returned. This allows the user to either skip the
C     block using L<SKPBLK> or to call his/her own reading
C     routine for a block type which he/she has defined.
C
C     The actual block name in string form is read into
C     L<IDLN> defined in dataid.h.
C
C     In the binary case, NREC is returned with a value
C     specifying the number of FORTRAN records in the
C     following block. This is used as a consistency check
C     by the block reading routine on reading, and by
C     L<SKPBLK> when skipping over an unwanted block.
C
C     In the binary case, if this routine encounters a FORMAT
C     block, it is read and L<CHCKFR> is called which checks
C     that the format of the file and the binary representation
C     of numbers of the program are compatible.
C     On return from L<CHCKFR>, the next header is read which
C     is then returned. In other words, checking of correct
C     format is transparent to client programs.
C
C      IN: IUNIT,LFORM
C     OUT: IDBLK,NREC
COMMEND
      SUBROUTINE RDBHDR(IUNIT,LFORM,IDBLK,NREC)
C includes 
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
C args
      INTEGER IUNIT,IDBLK,NREC
      LOGICAL LFORM
C externals
      INTEGER NDXID
      EXTERNAL NDXID
C local vars
      INTEGER IRES
      LOGICAL LOK
C begin
      IF (LFORM) THEN
         CALL GETLN(IUNIT,LOK)
         IF (LOK) THEN
            IDBLK = NDXID(FIOLIN)
            IDLN = FIOLIN(1:MAXIDL)
         ELSE
C getln failed because we have reached
C the end of a file or because of a format error
            IDBLK = IDERR
         ENDIF
C in formatted files NREC is not used.
C just set this to an 'impossible' value
         NREC = -1
      ELSE
C read unformatted
         READ(UNIT = IUNIT, IOSTAT = IRES) IDLN,NREC
         IF (IRES .EQ. 0) THEN
            IDBLK = NDXID(IDLN)
         ELSE
            IDBLK = IDERR
            NREC  = -1
         ENDIF

         IF (IDBLK .EQ. IFRMID) THEN
C Check here if we have the correct precision and binary format
C by reading the FRM block. Stop if it isnt the
C same as the one we have been compiled in.
            CALL CHCKFR(IUNIT)
C read the next block for returning
            READ(UNIT = IUNIT, IOSTAT = IRES) IDLN,NREC
            IF (IRES .EQ. 0) THEN
               IDBLK = NDXID(IDLN)
            ELSE
               IDBLK = IDERR
               NREC  = -1
            ENDIF
         ENDIF
      ENDIF
C end RDBHDR
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR SKPBLK
C     SUBROUTINE SKPBLK(IUNIT,LFORM,NREC)
C
C     SKPBLK is used to skip over a block on the file.
C     The routine assumes that the header of the
C     block to be skipped has been read (by L<RDBHDR>).
C
C     In the formatted case (LFORM = .TRUE.), SKPBLK
C     skips over lines until an END marker is encountered.
C
C     In the binary case, NREC "empty" reads are performed
C     In order to skip over the block.
C     The correct value for NREC is returned by the
C     subroutine L<RDBHDR> which must then be passed
C     to this routine in order to skip the block.
COMMEND
      SUBROUTINE SKPBLK(IUNIT,LFORM,NREC)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
C args
      INTEGER IUNIT,NREC
      LOGICAL LFORM
C local vars
      INTEGER I,IRES,NDX
      LOGICAL LOK
C external
      INTEGER NDXID
      EXTERNAL NDXID
C begin
      IF (LFORM) THEN
C formatted read: read until we get an END marker
C NREC is ignored in this case
         CALL GETLN(IUNIT,LOK)
         NDX = NDXID(FIOLIN)
C while ID == IERRID AND LOK do
 10      IF (NDX .EQ. IDUKN .AND. LOK) THEN
            CALL GETLN(IUNIT,LOK)
            NDX = NDXID(FIOLIN)
C end while
            GOTO 10
         ENDIF
         IF (.NOT. LOK) THEN
            PRINT *,'SKPBLK: read error while skipping block!'
            PRINT *,'(probably END marker missing)'
            PRINT *,'LFORM = .TRUE.'
            PRINT *
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN('SKPBLK')
            PRINT *,ENDEXP
            PRINT *,'but got ',NAMEID(NDX)
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ELSE
C unformatted read
         I = 1
         LOK = .TRUE.
 20      IF (I .LE. NREC .AND. LOK) THEN
            READ(UNIT = IUNIT, IOSTAT = IRES)
            LOK = (IRES .EQ. 0)
            I = I+1
C end while
            GOTO 20
         ENDIF
      ENDIF
C end SKPBLK
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRXRED
C     SUBROUTINE WRXRED(IUNIT,LFORM,LDO4,NAT,NDIM,X)
C
C     WRXRED writes coordinates to a file in a "reduced
C     position format".
C     The data can be read by L<RDXRED>.
COMMEND

      SUBROUTINE WRXRED(IUNIT,LFORM,LDO4,NAT,NDIM,X)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NAT,NDIM
      LOGICAL LFORM,LDO4
      real*8 X(NDIM*NAT)
C begin
      CALL DMPRED(IUNIT,LFORM,'WRXRED',IPORID,NAT,NDIM,X)
      IF (LDO4 .AND. NDIM .GT. 3)THEN
         CALL DMP4RE(IUNIT,LFORM,'WRXRED',IP4RID,NAT,NDIM,X)
      ENDIF
C end WRXRED
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDXRED
C     SUBROUTINE RDXRED(IUNIT,LFORM,LDO4,NREC,NAT,NDIM,X)
C
C     RDXRED reads coordinates written by L<WRXRED>
C     from an open file.
COMMEND
      SUBROUTINE RDXRED(IUNIT,LFORM,LDO4,NREC,NAT,NDIM,X)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NAT,NDIM,NREC
      LOGICAL LFORM,LDO4
      real*8 X(NDIM*NAT)
C local vars
      INTEGER IDBLK
      CHARACTER PRGSTR*(6)
      DATA PRGSTR /'RDXRED'/
C begin
      CALL GETRED(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
      IF (LDO4 .AND. NDIM .GT. 3) THEN
         CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
         IF (IDBLK .EQ. IP4RID) THEN
            CALL GET4RE(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
         ELSE
            IF (LFORM) THEN
               CALL FLAGLN(PRGSTR)
            ENDIF
            PRINT *,NAMEID(IP4RID), 'expected!'
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end RDXRED
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDXRE1
C     SUBROUTINE RDXRE1(IUNIT,LFORM,LDO4,NREC,MAXNAT,NDIM,NATRED,X)
C
C     RDXRED reads coordinates written by L<WRXRED>
C     from an open file.
COMMEND
      SUBROUTINE RDXRE1(IUNIT,LFORM,LDO4,NREC,MAXNAT,NDIM,NATRED,X)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,MAXNAT,NATRED,NDIM,NREC
      LOGICAL LFORM,LDO4
      real*8 X(NDIM*MAXNAT)
C local vars
      INTEGER IDBLK
      CHARACTER PRGSTR*(6)
      DATA PRGSTR /'RDXRE1'/
C begin
      CALL GETRE1(IUNIT,LFORM,NREC,PRGSTR,MAXNAT,NDIM,NATRED,X)
      IF (LDO4 .AND. NDIM .GT. 3) THEN
         CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
         IF (IDBLK .EQ. IP4RID) THEN
            CALL GET4RE(IUNIT,LFORM,NREC,PRGSTR,NATRED,NDIM,X)
         ELSE
            IF (LFORM) THEN
               CALL FLAGLN(PRGSTR)
            ENDIF
            PRINT *,NAMEID(IP4RID), 'expected!'
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end RDXRE1
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRXPOS
C     SUBROUTINE WRXPOS(IUNIT,LFORM,LDO4,NAT,NPM,NSM,NDIM,X)
C
C     WRXPOS writes out atom coordinates in a 'verbose' form
C     with the atom and residue names before each atom.
C     This is the preferred format for single configuration files
C     such as startup and final configurations of an MD run
C     as it allows easier inspection of the coordinates
C     by humans and other animals.
C     For trajectories, the 'reduced form' is used (See L<WRXRED>).
C
C     WRXPOS will write out coordinates in three or
C     four dimensions by calling L<DMPARR> for three dims
C     and  L<DMPARR> and L<DMP4AR> for four dims
C     respectively.
C
C     WRXPOS requires that:
C     a topology has been read in so that the atom and residue
C     names can be written.
C     This is tested in the routine. If the condition is not
C     fulfilled, a message is printed and STOP is invoked.
COMMEND
      SUBROUTINE WRXPOS(IUNIT,LFORM,LDO4,NAT,NPM,NSM,NDIM,X)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NAT,NDIM,NPM,NSM
      LOGICAL LFORM,LDO4
      real*8 X(*)
C begin
      CALL DMPARR(IUNIT,LFORM,'WRXPOS',IPOSID,NAT,NPM,NSM,NDIM,X)
      IF (LDO4 .AND. NDIM .GT. 3) THEN
         CALL DMP4AR(IUNIT,LFORM,'WRXPOS',IP4DID,NAT,NPM,NSM,NDIM,X)
      ENDIF
C end WRXPOS
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDXPOS
C     SUBROUTINE RDXPOS(IUNIT,LFORM,LDO4,NREC,NAT,NDIM,X)
C
C     RDXPOS reads position coordinates written by L<WRXPOS>.
C     
C     In the formatted case, the first 24 characters on
C     each line are ignored.
COMMEND
      SUBROUTINE RDXPOS(IUNIT,LFORM,LDO4,NREC,NAT,NDIM,X)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NREC,NAT,NDIM
      LOGICAL LFORM,LDO4
      real*8 X(NDIM*NAT)
C local vars
      INTEGER IDBLK
      CHARACTER PRGSTR*(6)
      DATA PRGSTR /'RDXPOS'/
C begin
      CALL GETARR(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
      IF (LDO4 .AND. NDIM .GT. 3) THEN
         CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
         IF (IDBLK .EQ. IP4DID) THEN
            CALL GET4AR(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
         ELSE
            IF (LFORM) THEN
               CALL FLAGLN(PRGSTR)
            ENDIF
            PRINT *,NAMEID(IP4DID), 'expected!'
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end RDXPOS
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDXPO1
C     SUBROUTINE RDXPO1(IUNIT,LFORM,LDO4,NREC,MAXNAT,NDIM,NATRED,X)
C
C     RDXPO1 reads position coordinates written by L<WRXPOS>.
C     It has slightly different behaviour than to that of L<RDXPOS>.
C     
C     In the formatted case, the first 24 characters on
C     each line are ignored.
C
C     This routine returns the number of atoms read in NATRED.
C     It calls L<GETAR1> and L<GET4AR>.
COMMEND
      SUBROUTINE RDXPO1(IUNIT,LFORM,LDO4,NREC,MAXNAT,NDIM,NATRED,X)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NREC,MAXNAT,NDIM,NATRED
      LOGICAL LFORM,LDO4
      real*8 X(NDIM*MAXNAT)
C local vars
      INTEGER IDBLK
      CHARACTER PRGSTR*(6)
      DATA PRGSTR /'RDXPO1'/
C begin
      CALL GETAR1(IUNIT,LFORM,NREC,PRGSTR,MAXNAT,NDIM,NATRED,X)
      IF (LDO4 .AND. NDIM .GT. 3) THEN
         CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
         IF (IDBLK .EQ. IP4DID) THEN
            CALL GET4AR(IUNIT,LFORM,NREC,PRGSTR,NATRED,NDIM,X)
         ELSE
            IF (LFORM) THEN
               CALL FLAGLN(PRGSTR)
            ENDIF
            PRINT *,NAMEID(IP4DID), 'expected!'
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end RDXPO1
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRXOF
C     SUBROUTINE WRXOF(IUNIT,LFORM,NATTOT,NPM,NSM,NDIM,X,NEL)
C
C     WRXOF writes oblique coordinates to file. 
C     The file must have been opened formatted and a topology
C     must have been read in prior to calling this routine.
C     Only those solute atoms I for which L<IAC>(I) is nonzero
C     are written to file.
C     Only those solvent atoms I for which L<IACS>(I) is nonzero
C     are written to file.
C
C     The block can be read with L<RDXOF>.
COMMEND
      SUBROUTINE WRXOF(IUNIT,LFORM,NATTOT,NPM,NSM,NDIM,X,NEL)
      INCLUDE 'dataid.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
C args
      INTEGER IUNIT,NATTOT,NPM,NSM,NDIM
      INTEGER NEL(NATTOT)
      LOGICAL LFORM
      real*8 X(NATTOT*NDIM)
C local vars
      CHARACTER PRGSTR*(5)
      CHARACTER SOLCRD*(24)
      INTEGER I,II,J,I3,IRES
C externals
      LOGICAL LTOPOK
      EXTERNAL LTOPOK
C data
      DATA PRGSTR /'WRXOF'/
      DATA SOLCRD /'solute atom coordinates'/
C begin
      IF (.NOT. LTOPOK()) THEN
         PRINT *,PRGSTR,
     $        ': need a topology to write',
     $        ' oblique coordinates!'
         PRINT *
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

      IF (.NOT. LFORM) THEN
         PRINT *,PRGSTR,': implementation restriction'
         PRINT *,'file must be formatted !'
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

      IF (NPM*NRP + NSM*NRAM .NE. NATTOT) THEN
         PRINT *,PRGSTR,':inconsistent number of atoms!'
         PRINT *,'on topology/controlfile:',NRP*NPM,' solute atoms and',
     $        NSM*NRAM,' solvent atoms.'
         PRINT *,'total: ',NRP*NPM + NSM*NRAM
         PRINT *,'but  wanted to write NAT = ',NATTOT,' atoms'
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

C write block name
      CALL WRSTR(IUNIT,PRGSTR,NAMEID(IPOBID))

      I3 = 1
      II = 1
      DO 10 I=1,NPM
         DO 20 J=1,NRP
            IF (IAC(J).NE.0) THEN
               WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=500)
     $              X(I3),X(I3+1),X(I3+2),
     $              0.0D0, 1.0D0, NEL(IAC(J)), MRES(J), PANM(J),
     $              AANM(MRES(J))
               IF (IRES .NE. 0) THEN
                  PRINT FWFAIL,PRGSTR,SOLCRD,IUNIT
               ENDIF
            ENDIF
            II = II + 1
            I3 = I3 + NDIM
 20      CONTINUE
 10   CONTINUE

C write solvent atoms
      DO 30 I=1,NSM
         DO 40 J=1,NRAM
            IF (IACS(J).NE.0) THEN
               WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=500)
     $              X(I3),X(I3+1),X(I3+2),
     $              0.0D0,1.0D0,NEL(IACS(J)),I,ANMS(J),'SOLV '
               IF (IRES .NE. 0) THEN
                  PRINT FWFAIL,PRGSTR,'solvent atom coordinates',IUNIT
                  CALL STOPGM('BLOCKIO','FATAL ERROR')
               ENDIF
            ENDIF
            II = II + 1
            I3 = I3 + NDIM
 40      CONTINUE
 30   CONTINUE

C     write end marker
      CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
 500  FORMAT (5X,3F10.5,F10.5,F5.2,I5,7X,I3,1X,A5,1X,A5)
C wrxof
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDXOF
C     SUBROUTINE RDXOF(IUNIT,LFORM,NATTOT,NDIM,X,NEL)
C
C     RDXOF reads an block containing oblique coordinates
C     written by L<WRXOF>.
COMMEND
      SUBROUTINE RDXOF(IUNIT,LFORM,NATTOT,NDIM,X,NEL)
      INCLUDE 'dataid.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
C args
      INTEGER IUNIT,NATTOT,NDIM
      INTEGER NEL(NATTOT)
      LOGICAL LFORM
      real*8 X(NATTOT*NDIM)
C local vars
      CHARACTER PRGSTR*(5)
      CHARACTER STR*(24)
      INTEGER NATTMP,I3,NDX
      real*8 RDUMM
C data
      DATA PRGSTR /'RDXOF'/
      DATA STR /'coordinate data'/

C begin
      IF (.NOT. LFORM) THEN
         PRINT *,PRGSTR,': implementation restriction'
         PRINT *,'file must be formatted !'
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

      NATTMP = 0
      I3 = 1
      CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C     while no marker  AND still atoms to read
 10   IF (NDX .EQ. IDUKN .AND. NATTMP .LT. NATTOT) THEN
         NATTMP = NATTMP + 1

         CALL CHPNRE(3,X(I3),'X coord',PRGSTR)
         CALL CHPREL(RDUMM,'RDUMM',PRGSTR)
         CALL CHPREL(RDUMM,'RDUMM',PRGSTR)
         CALL CHPINT(NEL(NATTMP),'NEL',PRGSTR)

         CALL CHKMTY(PRGSTR,STR)

         I3 = I3 + NDIM

         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C end while
         GOTO 10
      ENDIF

      IF (NDX .NE. IENDID) THEN
         CALL FLAGLN(PRGSTR)
         PRINT 500,NATTOT
         PRINT *,ENDEXP
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ELSEIF (NATTMP .NE. NATTOT) THEN
         CALL FLAGLN(PRGSTR)
         PRINT 500,NATTOT
         PRINT BONFIL,NATTMP
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

 500  FORMAT (1X,'attempting to read',I7,' atom coordinates')
C rdxof
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRVRED
C     SUBROUTINE WRVRED(IUNIT,LFORM,LWR4,NAT,NDIM,X)
C
C     WRVRED writes velocities to a file.
C     The data can be read by L<RDVRED>
COMMEND
      SUBROUTINE WRVRED(IUNIT,LFORM,LWR4,NAT,NDIM,X)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NAT,NDIM
      LOGICAL LFORM,LWR4
      real*8 X(NDIM*NAT)
C begin
      CALL DMPRED(IUNIT,LFORM,'WRVRED',IVRDID,NAT,NDIM,X)
      IF (LWR4 .AND. NDIM .GT. 3) THEN
         CALL DMP4RE(IUNIT,LFORM,'WRVRED',IV4RID,NAT,NDIM,X)
      ENDIF
C end WRVRED
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDVRED
C     SUBROUTINE RDVRED(IUNIT,LFORM,LDO4,NREC,NAT,NDIM,X)
C
C     RDVRED reads velocities written by L<WRVRED>
C     from an open file.
COMMEND
      SUBROUTINE RDVRED(IUNIT,LFORM,LDO4,NREC,NAT,NDIM,X)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NAT,NDIM,NREC
      LOGICAL LFORM,LDO4
      real*8 X(NDIM*NAT)
C local vars
      INTEGER IDBLK
      CHARACTER PRGSTR*(6)
      DATA PRGSTR /'RDVRED'/
C begin
      CALL GETRED(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
      IF (LDO4 .AND. NDIM .GT. 3) THEN
         CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
         IF (IDBLK .EQ. IV4RID) THEN
            CALL GET4RE(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
         ELSE
            IF (LFORM) THEN
               CALL FLAGLN(PRGSTR)
            ENDIF
            PRINT *,NAMEID(IV4RID), 'expected!'
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end RDVRED
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRXVAL
C      SUBROUTINE WRXVEL(IUNIT,LFORM,LDO4,NAT,NPM,NSM,NDIM,X)
C
COMMEND
      SUBROUTINE WRXVEL(IUNIT,LFORM,LDO4,NAT,NPM,NSM,NDIM,X)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NAT,NDIM,NPM,NSM
      LOGICAL LFORM,LDO4
      real*8 X(NDIM*NAT)
C begin
      CALL DMPARR(IUNIT,LFORM,'WRXVEL',IVELID,NAT,NPM,NSM,NDIM,X)
      IF (LDO4 .AND. NDIM .GT. 3) THEN
         CALL DMP4AR(IUNIT,LFORM,'WRXVEL',IV4DID,NAT,NPM,NSM,NDIM,X)
      ENDIF
C end WRXVEL
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDXVEL
C     SUBROUTINE RDXVEL(IUNIT,LFORM,LDO4,NREC,NAT,NDIM,X)
C
C     RDXVEL reads a block written by L<WRXVEL>.
COMMEND
      SUBROUTINE RDXVEL(IUNIT,LFORM,LDO4,NREC,NAT,NDIM,X)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NREC,NAT,NDIM
      LOGICAL LFORM,LDO4
      real*8 X(NDIM*NAT)
C local vars
      INTEGER IDBLK
      CHARACTER PRGSTR*(6)
      DATA PRGSTR /'RDXVEL'/
C begin
      CALL GETARR(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
      IF (LDO4 .AND. NDIM .GT. 3) THEN
         CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
         IF (IDBLK .EQ. IV4DID) THEN
            CALL GET4AR(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
         ELSE
            IF (LFORM) THEN
               CALL FLAGLN(PRGSTR)
            ENDIF
            PRINT *,NAMEID(IV4DID), 'expected!'
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end RDXVEL
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRSX
C     SUBROUTINE WRSX(IUNIT,LFORM,LWR4,NAT,NPM,NSM,NDIM,SX,IG)
C
C     WRSX writes a stochastic integral block
C     to IUNIT. The block can be read by L<RDSX>.
COMMEND
      SUBROUTINE WRSX(IUNIT,LFORM,LWR4,NAT,NPM,NSM,NDIM,SX,IG)
      INCLUDE 'dataid.h'
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
C     args
      INTEGER IUNIT,NAT,NDIM,NPM,NSM
      LOGICAL LFORM,LWR4
      INTEGER IG
      real*8 SX(NDIM*NAT)
C externals
      LOGICAL LTOPOK
      EXTERNAL LTOPOK
C local vars
      INTEGER I,J,IL,INA,IRES,NREC,I3
      CHARACTER SOLCRD*(24)
      CHARACTER PRGSTR*(4)
C data
      DATA PRGSTR /'WRSX'/
      DATA SOLCRD /'solute atom coordinates'/
C begin
      IF (LFORM) THEN
C check we have read the topology
         IF (.NOT. LTOPOK()) THEN
            PRINT *,PRGSTR,
     $           ': need a topology to write verbose coordinate format!'
            PRINT *
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

C check the number of atoms is consistent
         IF (NRP*NPM + NSM*NRAM .NE. NAT) THEN
            PRINT *,PRGSTR,':inconsistent number of atoms!'
            PRINT *,'on topology/controlfile:',NRP*NPM,
     $           ' solute atoms and',
     $           NSM*NRAM,' solvent atoms.'
            PRINT *,'total: ',NRP*NPM + NSM*NRAM
            PRINT *,'but  wanted to write NAT = ',NAT,' atoms'
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

         CALL WRSTR(IUNIT,PRGSTR,NAMEID(ISDIID))
C write a little comment
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=FMIGNO)
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,SOLCRD,IUNIT
            PRINT *,FMIGNO
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

C     write solute(s)
         IL = 1
         INA = 1
         DO 10 I=1, NPM
            DO 20 J=1,NRP
               WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=12)
     $              MRES(J),AANM(MRES(J)),PANM(J),INA,
     $              SX(IL),SX(IL+1),SX(IL+2)
               IF (IRES .NE. 0) THEN
                  PRINT FWFAIL,PRGSTR,SOLCRD,IUNIT
                  CALL STOPGM('BLOCKIO','FATAL ERROR')
               ENDIF
               IL = IL + NDIM
               INA = INA + 1
 20         CONTINUE
 10      CONTINUE

C write solvent atoms
         DO 30 I=1,NSM
            DO 40 J=1,NRAM
               WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=12)
     $              I,'SOLV ',ANMS(J),INA,SX(IL),SX(IL+1),SX(IL+2)
               IF (IRES .NE. 0) THEN
                  PRINT FWFAIL,PRGSTR,'solvent atom coordinates',IUNIT
                  CALL STOPGM('BLOCKIO','FATAL ERROR')
               ENDIF
               IL = IL + NDIM
               INA = INA + 1
 40         CONTINUE
 30      CONTINUE

C now write IG
         CALL WRCOMM(IUNIT,PRGSTR,'IG')
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT='(I10)')IG

C     write end marker
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
C     binary case
         NREC = 3
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAMEID(ISDIID),NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAT
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'NAT',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES)
     $        ((SX(NDIM*I+I3),I3=1,3),I=0,NAT-1)
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'coordinate data',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES)IG
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'IG',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF

C now write 4D part if required
      IF (LWR4 .AND. NDIM .GT. 3) THEN
         IF (LFORM) THEN
C We cant use DMP4AR here, as we write  stochastic integrals and
C positions in a different format
            CALL WRSTR(IUNIT,PRGSTR,NAMEID(ISD4ID))
            WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=FMIGNO)
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,SOLCRD,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C     write solute(s)
            IL = 1
            INA = 1
            DO 110 I=1, NPM
               DO 120 J=1,NRP
                  WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=12)
     $                 MRES(J),AANM(MRES(J)),PANM(J),INA,SX(IL+3)
                  IF (IRES .NE. 0) THEN
                     PRINT FWFAIL,PRGSTR,SOLCRD,IUNIT
                     CALL STOPGM('BLOCKIO','FATAL ERROR')
                  ENDIF
                  IL = IL + NDIM
                  INA = INA + 1
 120           CONTINUE
 110        CONTINUE

C     write solvent atoms
            DO 130 I=1,NSM
               DO 140 J=1,NRAM
                  WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=12)
     $                 I,'SOLV ',ANMS(J),INA,SX(IL+3)
                  IF (IRES .NE. 0) THEN
                     PRINT FWFAIL,PRGSTR,'solvent atom block',IUNIT
                     CALL STOPGM('BLOCKIO','FATAL ERROR')
                  ENDIF
                  IL = IL + NDIM
                  INA = INA + 1
 140           CONTINUE
 130        CONTINUE
C     write end marker
            CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
         ELSE
C     binary 4D write
            CALL DMP4AR(IUNIT,LFORM,'WRSX',ISD4ID,NAT,NPM,NSM,NDIM,SX)
         ENDIF
      ENDIF

 12   FORMAT(I5,2(1X,A5),I7,3E15.7)
C end WRSX
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDSX
C     SUBROUTINE RDSX(IUNIT,LFORM,LDO4,NREC,NAT,NDIM,SX,IG)
C
C     RDSX reads a stochastic integral block
C     from IUNIT which has been written by L<WRSX>.
COMMEND 
      SUBROUTINE RDSX(IUNIT,LFORM,LDO4,NREC,NAT,NDIM,SX,IG)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NAT,NDIM,NREC
      LOGICAL LFORM,LDO4
      INTEGER IG
      real*8 SX(NDIM*NAT)
C externals
      LOGICAL LISMTY
      EXTERNAL LISMTY
C local vars
      INTEGER IDBLK,IRES,I,I3
      INTEGER NATTMP,NDX
      CHARACTER STR*(24)
      CHARACTER PRGSTR*(4)
      DATA PRGSTR /'RDSX'/
      DATA STR /'coordinate data'/
C begin
C      CALL GETARR(IUNIT,LFORM,NREC,'RDSX',NAT,NDIM,SX)
      IF (LFORM) THEN
C formatted read: read until we get an END marker or
C we have the right number of atoms
         NATTMP = 0
         I3 = 1
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C     while no marker  AND still atoms to read
 10      IF (NDX .EQ. IDUKN .AND. NATTMP .LT. NAT) THEN
            NATTMP = NATTMP + 1
            CALL SKPCHR(24)
            IF (LISMTY()) THEN
               CALL FLAGLN(PRGSTR)
               PRINT MTYLIN,24
               PRINT STREXP,STR
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF

            CALL CHPNRE(3,SX(I3),'SX coord',PRGSTR)
            CALL CHKMTY(PRGSTR,STR)

            I3 = I3 + NDIM

            CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C end while
            GOTO 10
         ENDIF

         IF (NATTMP .NE. NAT) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 500,NAT
            PRINT BONFIL,NATTMP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C get IG
         CALL CHPINT(IG,'IG',PRGSTR)
         CALL GTEND(IUNIT,PRGSTR)
      ELSE
C unformatted reading
         CALL CHNREC(PRGSTR,NREC,3)
         READ(UNIT = IUNIT, IOSTAT = IRES) NATTMP
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'NAT',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         IF (NATTMP .NE. NAT) THEN
            PRINT 500,IUNIT
            PRINT EXPFIL,NAT
            PRINT BONFIL,NATTMP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         READ(UNIT = IUNIT, IOSTAT = IRES)
     $        ((SX(NDIM*I+I3),I3=1,3),I=0,NAT-1)

         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'SX data',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         READ(UNIT = IUNIT, IOSTAT = IRES)IG
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'IG',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF

C now read 4D part if required
      IF (LDO4 .AND. NDIM .GT. 3) THEN
         CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
         IF (IDBLK .EQ. ISD4ID) THEN
            CALL GET4AR(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,SX)
         ELSE
            IF (LFORM) THEN
               CALL FLAGLN(PRGSTR)
            ENDIF
            PRINT *,NAMEID(ISD4ID), 'expected!'
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF

 500  FORMAT (1X,'inconsistency in NAT on UNIT',I5)
C end RDSX
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRXC
C     SUBROUTINE WRXC(IUNIT,LFORM,NAT,NPM,NSM,NDIM,XC)
C
C     WRXC writes a position re(con)straint reference
C     position coordinate block to IUNIT.
COMMEND
      SUBROUTINE WRXC(IUNIT,LFORM,NAT,NPM,NSM,NDIM,XC)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NAT,NPM,NSM,NDIM
      LOGICAL LFORM
      real*8 XC(NDIM*NAT)
C begin
      CALL DMPARR(IUNIT,LFORM,'WRXC',IPRCID,NAT,NPM,NSM,NDIM,XC)
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDXC
C     SUBROUTINE RDXC(IUNIT,LFORM,NREC,NAT,NDIM,XC)
C
C     RDXC reads a position re(con)straint reference
C     position coordinate block from IUNIT.
COMMEND 
      SUBROUTINE RDXC(IUNIT,LFORM,NREC,NAT,NDIM,XC)
C args
      INTEGER IUNIT,NAT,NDIM,NREC
      LOGICAL LFORM
      real*8 XC(NDIM*NAT)
C begin
      CALL GETARR(IUNIT,LFORM,NREC,'RDXC',NAT,NDIM,XC)
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDDRST
C      SUBROUTINE RDDRST(IUNIT,NRPT,MAXNDR,NDR,
C     $     DISH,DISC,
C     $     IDR1,JDR1,KDR1,LDR1,ICDR1,
C     $     IDR2,JDR2,KDR2,LDR2,ICDR2,
C     $     R0,W0)
C
C     RDDRST reads in a block of distance restraints
C     from an open formatted file.
C
C     If the number of distance restraints in the block
C     is larger than MAXNDR or any atom index larger
C     than NRPT, an error is written to screen and the
C     program is terminated.
C
C     NDR is returned with the number of distance restraints
C     read in.
C
C     This routine calls L<CHPDRS> twice per line
C     in order to read in the I,J,K,L, and IC variables
C     and check for consistency.
COMMEND
      SUBROUTINE RDDRST(IUNIT,NRPT,MAXNDR,NDR,
     $     DISH,DISC,
     $     IDR1,JDR1,KDR1,LDR1,ICDR1,
     $     IDR2,JDR2,KDR2,LDR2,ICDR2,
     $     R0,W0)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
C args
      INTEGER IUNIT,NRPT,MAXNDR,NDR
      INTEGER IDR1(MAXNDR),JDR1(MAXNDR),KDR1(MAXNDR)
      INTEGER LDR1(MAXNDR),ICDR1(MAXNDR)
      INTEGER IDR2(MAXNDR),JDR2(MAXNDR),KDR2(MAXNDR)
      INTEGER LDR2(MAXNDR),ICDR2(MAXNDR)
      real*8    DISH,DISC
      real*8    R0(MAXNDR),W0(MAXNDR)
C external
      INTEGER NDXID
      EXTERNAL NDXID
C local vars
      INTEGER NDX
      CHARACTER PRGSTR*(6)
      CHARACTER STR*(12)
      LOGICAL LOK
C data
      DATA PRGSTR /'RDDRST'/
      DATA STR    /'DISH,DISC'/
C begin
C get DISH, DISC
      CALL GIMME(IUNIT,PRGSTR,STR)
      CALL CHPREL(DISH,PRGSTR,STR)
      CALL CHPREL(DISC,PRGSTR,STR)
      CALL CHKMTY(PRGSTR,STR)
C
      IF (DISH .LE. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         PRINT 500,'DISH',DISH
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

      IF (DISC .LE. 0.0D0) THEN
         CALL FLAGLN(PRGSTR)
         PRINT 500,'DISC',DISC
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

      NDR = 0
C read until we get an END marker
      CALL GETLN(IUNIT,LOK)
      NDX = NDXID(FIOLIN)
C while ID == IERRID AND LOK do
 10   IF (NDX .EQ. IDUKN .AND. LOK) THEN
         NDR = NDR+1

         IF (NDR .GT. MAXNDR) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'too many distance restraints on file!'
            PRINT *,'maximum allowed: MAXNDR = ',MAXNDR
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

C read from the data line and check for correctness
         CALL CHPDSR(NRPT,PRGSTR,NDR,IDR1(NDR),JDR1(NDR),
     $        KDR1(NDR),LDR1(NDR),ICDR1(NDR))

         CALL CHPDSR(NRPT,PRGSTR,NDR,IDR2(NDR),JDR2(NDR),
     $        KDR2(NDR),LDR2(NDR),ICDR2(NDR))

C now R0 and W0
         CALL CHPREL(R0(NDR),PRGSTR,'R0')
         CALL CHPREL(W0(NDR),PRGSTR,'W0')
         CALL CHKMTY(PRGSTR,STR)

         CALL GETLN(IUNIT,LOK)
         NDX = NDXID(FIOLIN)
C end while
         GOTO 10
      ENDIF

      IF (.NOT. LOK) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,'read error while reading dist. rest. block!'
         PRINT *,'(probably END marker missing)'
         PRINT *
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

      IF (NDX .NE. IENDID) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,ENDEXP
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

 500  FORMAT(1X,'illegal value for ',A,' = ',F10.5)
C end rddrst
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR CHPDIH
C     SUBROUTINE CHPDIH(MIN,MAX,PRGSTR,STRLIM,NDXNUM,I,J,K,L,LOK)
C
C     A routine to get four integers specifying a dihedral angle.
C     LOK is returned with .TRUE. if the four integers I,J,K and L
C     were read successfully and MIN <= I,J,K,L, <= MAX .
COMMEND
      SUBROUTINE CHPDIH(MIN,MAX,PRGSTR,STRLIM,NDXNUM,I,J,K,L,LOK)
C args
      INTEGER MIN,MAX,NDXNUM,I,J,K,L
      LOGICAL LOK
      CHARACTER PRGSTR*(*)
      CHARACTER STRLIM*(*)
C local vars
      CHARACTER STR*(4)
      INTEGER NVAL(4)
      INTEGER NN
C data
      DATA STR /'IJKL'/
C begin
      LOK = .TRUE.
      DO 10 NN=1,4
         CALL CHPINT(NVAL(NN),PRGSTR,STR(NN:NN))
         IF (NVAL(NN) .LT. MIN .OR. NVAL(NN) .GT. MAX) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 500,STR(NN:NN),NDXNUM,NVAL(NN)
            PRINT 510,MIN,STRLIM,MAX
            LOK = .FALSE.
         ENDIF
 10   CONTINUE

      I = NVAL(1)
      J = NVAL(2)
      K = NVAL(3)
      L = NVAL(4)
 500  FORMAT (1X,'illegal value for ',A1,'(',I4,') =',I4)
 510  FORMAT (1X,'must be ',I3,'<= val <= ',A6,'=',I5)
C end chpdih
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR CHPDSR
C     SUBROUTINE CHPDSR(NRPT,PRGSTR,NDR,I,J,K,L,IC)
C
C     CHPDSR reads in and checks consistency of
C     I,J,K,L and IC which defines the atoms needed in L<DISRE>
C     to define pseudo/virtual atoms.
C
C     This is an auxiliary routine called by L<RDDRST> which reads
C     in a distance restraint specification block.
COMMEND
      SUBROUTINE CHPDSR(NRPT,PRGSTR,NDR,I,J,K,L,IC)
      INCLUDE 'disre.h'
      INCLUDE 'formats.h'
C args
      INTEGER I,J,K,L,IC,NRPT,NDR
      CHARACTER PRGSTR*(*)
C local vars
      LOGICAL LOK
      CHARACTER STR*(5)
C begin
      CALL CHPDIH(0,NRPT,PRGSTR,'NRPT',NDR,I,J,K,L,LOK)

      STR = 'HTYPE'
      CALL CHPINT(IC,PRGSTR,STR)
      IF (IC .LT. MINIHT .OR. IC .GT. MAXIHT) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMNII,STR,IC
         PRINT FMBTI,MINIHT,STR,MAXIHT
         LOK = .FALSE.
      ENDIF
      
      IF (.NOT. LOK) THEN
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

C I must always be a REAL*8 atom
      IF (I .LT. 1 .OR. I .GT. NRPT) THEN
         CALL FLAGLN(PRGSTR)
         PRINT 510,IC
         PRINT 520,I,J,K,L
         PRINT FMBTII,1,'I','NRPT',NRPT
         LOK = .FALSE.
      ENDIF

C now check for consistence
      IF (IC .EQ. IHTREL) THEN
C     real*8 HYDROGEN ATOM            (I)
         IF (J .NE. 0 .OR. K .NE. 0 .OR. L .NE. 0) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 510,IC
            PRINT 520,I,J,K,L
            PRINT *,'J,K and L must be 0 !'
            LOK = .FALSE.
         ENDIF

      ELSEIF (IC .EQ. IHLCH1) THEN
C     VIRTUAL H, ALIPHATIC CH1      (I,J,K,L)
         IF (J .EQ. 0 .OR. K .EQ. 0 .OR. L .EQ. 0) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 510,IC
            PRINT 520,I,J,K,L
            PRINT FMBTII,1,'I,J,K,L','NRPT',NRPT
            LOK = .FALSE.
         ENDIF
      ELSEIF (IC .EQ. IHRCH1 .OR. IC .EQ. IHLCH2 .OR. 
     $        IC .EQ. IVLCH2 .OR. IC .EQ. IH2CH3) THEN
C     VIRTUAL H, AROMATIC  CH1      (I,J,K)
C     PSEUDO  H, ALIPHATIC CH2      (I,J,K)
C     VIRTUAL H, ALIPHATIC CH2      (I,J,K)
C     PSEUDO  H, TWO CH3 GROUPS     (I,J,K) 
         IF (J .EQ. 0 .OR. K .EQ. 0 .OR. L .NE. 0) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 510,IC
            PRINT 520, I,J,K,L
            PRINT *,'J and K must be <> 0, L must be 0 !'
            LOK = .FALSE.
         ENDIF
      ELSEIF (IC .EQ. IH1CH3) THEN
C     PSEUDO  H, ONE CH3 GROUP      (I,J)
         IF (J .EQ. 0 .OR. K .NE. 0 .OR. L .NE. 0) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 510,IC
            PRINT 520, I,J,K,L
            PRINT *,'J must be <> 0, K and L must be 0 !'
            LOK = .FALSE.
         ENDIF
      ELSE
         PRINT *,PRGSTR,': unkown H type ',IC
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF
      
      IF (.NOT. LOK) THEN
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

 510  FORMAT (1X,'inconsistent I,J,K or L values for H type ',I2)
 520  FORMAT (1X,'I=',I4,/,1X,'J=',I4,/,1X,'K=',I4,/,1X,'L=',I4)
C end CHPDSR
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDDHRS
C     SUBROUTINE RDDHRS(IUNIT,NRPT,MAXNDL,NDLR,
C     $     IPLR,JPLR,KPLR,LPLR,CPLR,PDLR,ICPLR)
C
C     RDDHRS reads in a block of dihedral restraint specifications
C     from an open formatted file.
C     If the number of restraints in the block is larger than MAXNDL,
C     then an error is written to screen and the program is terminated.
C
C     NO CONVERSION OF THE VALUES READ IN IS PERFORMED IN THIS ROUTINE.
C
C     NDLR is returned with the number of distance restraints
C     read in.
COMMEND
      SUBROUTINE RDDHRS(IUNIT,NRPT,MAXNDL,NDLR,
     $     IPLR,JPLR,KPLR,LPLR,CPLR,PDLR,ICPLR)
C
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
C args
      INTEGER IUNIT,NRPT,MAXNDL,NDLR
      INTEGER IPLR(MAXNDL),JPLR(MAXNDL),KPLR(MAXNDL)
      INTEGER LPLR(MAXNDL),ICPLR(MAXNDL)
      real*8 CPLR(MAXNDL),PDLR(MAXNDL)
C external
      INTEGER NDXID
      EXTERNAL NDXID
C local vars
      INTEGER NDX
      LOGICAL LOK
      CHARACTER PRGSTR*(6)
      CHARACTER STR*(16)
C data
      DATA PRGSTR /'RDDHRS'/
      DATA STR    /'I,J,K,L, CP, PD'/
C begin
      NDLR = 0
C read until we get an END marker
      CALL GETLN(IUNIT,LOK)
      NDX = NDXID(FIOLIN)
C while ID = IDUKN AND LOK do
 10   IF (NDX .EQ. IDUKN .AND. LOK) THEN
         NDLR = NDLR+1

         IF (NDLR .GT. MAXNDL) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 500
            PRINT *,'maximum allowed: MAXNDL = ',MAXNDL
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

C--   read from the data line and check for correctness
         CALL CHPDIH(1,NRPT,PRGSTR,'NRPT',NDLR,
     $        IPLR(NDLR),JPLR(NDLR),KPLR(NDLR),LPLR(NDLR),LOK)

         CALL CHPREL(CPLR(NDLR),PRGSTR,STR)
         CALL CHPREL(PDLR(NDLR),PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)

C     NOTE: with normal dihedrals, the icp array points to a dihedral
C     type (defining the CP and PD values).
C     In dihedral restraints, however, we read in the CP and PD
C     values for each dihe. rest. _individually_, and so need to
C     set ICPLR below for the routine DIHANG to work.
C     end of NOTE

         ICPLR(NDLR) = NDLR

C--   end of reading one line

         CALL GETLN(IUNIT,LOK)
         NDX = NDXID(FIOLIN)
C end while
         GOTO 10
      ENDIF

      IF (.NOT. LOK) THEN
         CALL FLAGLN(PRGSTR)
         PRINT 510
         PRINT *,'(probably END marker missing)'
         PRINT *
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

      IF (NDX .NE. IENDID) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,ENDEXP
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF
 500  FORMAT (1X,'too many dihedral restraint specs on file!')
 510  FORMAT (1X,'read error while reading dihedral restraining ',
     $     'spec. block')
C end rddhrs
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDJVSP
C     SUBROUTINE RDJVSP(IUNIT,NRP,MAXNDL,NDLR,
C     $     IPLR,JPLR,KPLR,LPLR,
C     $     CPLR,PJR0,PSJR,A,B,C)
C
C     RDJVSP reads in a block of j-value restraint specifications
C     from an open formatted file.
C     If the number of restraints in the block
C     is larger than MAXNDL then an error is written to
C     screen and the program is terminated.
C
C     NDLR is returned with the number of distance restraints
C     read in.
COMMEND
      SUBROUTINE RDJVSP(IUNIT,NRP,MAXNDL,NDLR,
     $     IPLR,JPLR,KPLR,LPLR,
     $     CPLR,PJR0,PSJR,A,B,C)
C
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
C args
      INTEGER IUNIT,NRP,MAXNDL,NDLR
      INTEGER IPLR(MAXNDL),JPLR(MAXNDL),KPLR(MAXNDL)
      INTEGER LPLR(MAXNDL)
      real*8 CPLR(MAXNDL),PJR0(MAXNDL),PSJR(MAXNDL)
      real*8 A(MAXNDL),B(MAXNDL),C(MAXNDL)
C external
      INTEGER NDXID
      EXTERNAL NDXID
C local vars
      CHARACTER PRGSTR*(6)
      INTEGER NDX
      LOGICAL LOK
      DATA PRGSTR /'RDJVSP'/
C begin
      NDLR = 0
C read until we get an END marker
      CALL GETLN(IUNIT,LOK)
      NDX = NDXID(FIOLIN)
C while ID = IDUKN AND LOK do
 10   IF (NDX .EQ. IDUKN .AND. LOK) THEN
         NDLR = NDLR+1

         IF (NDLR .GT. MAXNDL) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 500
            PRINT *,'maximum allowed: MAXNDL = ',MAXNDL
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

C--   read from the data line and check for correctness
         CALL CHPDIH(1,NRP,PRGSTR,'NRP',NDLR,
     $        IPLR(NDLR),JPLR(NDLR),KPLR(NDLR),LPLR(NDLR),LOK)

         CALL CHPREL(CPLR(NDLR),PRGSTR,'CPLR')
         CALL CHPREL(PJR0(NDLR),PRGSTR,'PJR0')
         CALL CHPREL(PSJR(NDLR),PRGSTR,'PSJR')
         CALL CHPREL(A(NDLR),PRGSTR,'A')
         CALL CHPREL(B(NDLR),PRGSTR,'B')
         CALL CHPREL(C(NDLR),PRGSTR,'C')

C--   end of reading one line

         CALL GETLN(IUNIT,LOK)
         NDX = NDXID(FIOLIN)
C end while
         GOTO 10
      ENDIF

      IF (.NOT. LOK) THEN
         CALL FLAGLN(PRGSTR)
         PRINT 510
         PRINT *,'(probably END marker missing)'
         PRINT *
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

      IF (NDX .NE. IENDID) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,ENDEXP
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF
 500  FORMAT (1X,'too many j-value restraining specs on file!')
 510  FORMAT (1X,'read error while reading j-value restraining ',
     $     'spec. block')
C end rdjvsp
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDLESP
C     SUBROUTINE RDLESP(IUNIT,NRP,MAXDLE,NDLE
C     $     IPLR,JPLR,KPLR,LPLR)
C
C     RDLESP reads in a block of local elevation dihedral 
C     specifications from an open formatted file.
C     If the number of dihedrals  in the block
C     is larger than MAXDLE, an error is written to
C     screen and the program is terminated.
C
C     NDLE is returned with the number of LE dihedral specifications
C     read in.
COMMEND
      SUBROUTINE RDLESP(IUNIT,NRP,MAXDLE,NDLE,
     $     IPLE,JPLE,KPLE,LPLE)
C
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
C args
      INTEGER IUNIT,NRP,MAXDLE,NDLE
      INTEGER IPLE(MAXDLE),JPLE(MAXDLE),KPLE(MAXDLE)
      INTEGER LPLE(MAXDLE)
C external
      INTEGER NDXID
      EXTERNAL NDXID
C local vars
      CHARACTER PRGSTR*(6)
      INTEGER NDX
      LOGICAL LOK
      DATA PRGSTR /'RDLESP'/
C begin
      NDLE = 0
C read until we get an END marker
      CALL GETLN(IUNIT,LOK)
      NDX = NDXID(FIOLIN)
C while ID = IDUKN AND LOK do
 10   IF (NDX .EQ. IDUKN .AND. LOK) THEN
         NDLE = NDLE+1

         IF (NDLE .GT. MAXDLE) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 500
            PRINT *,'maximum allowed: MAXDLE = ',MAXDLE
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

C--   read from the data line and check for correctness
         CALL CHPDIH(1,NRP,PRGSTR,'NRP',NDLE,
     $        IPLE(NDLE),JPLE(NDLE),KPLE(NDLE),LPLE(NDLE),LOK)
         CALL CHKMTY(PRGSTR,'I, J, K, L ')

C--   end of reading one line

         CALL GETLN(IUNIT,LOK)
         NDX = NDXID(FIOLIN)
C end while
         GOTO 10
      ENDIF

      IF (.NOT. LOK) THEN
         CALL FLAGLN(PRGSTR)
         PRINT 510
         PRINT *,'(probably END marker missing)'
         PRINT *
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF

      IF (NDX .NE. IENDID) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,ENDEXP
         CALL STOPGM('BLOCKIO','FATAL ERROR')
      ENDIF
 500  FORMAT (1X,'too many local elevation dihedrals on file!')
 510  FORMAT (1X,'read error while reading local ',
     $     'elevation spec. block !')
C end rdlesp
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRTIME
C     SUBROUTINE WRTIME(IUNIT,LFORM,NSTEP,T)
C
C     Write NSTEP and T to a TIMESTEP block.
C     The block can be read using L<RDTIME>.
COMMEND
      SUBROUTINE WRTIME(IUNIT,LFORM,NSTEP,T)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NSTEP
      LOGICAL LFORM
      real*8 T
C local vars
      INTEGER IRES,NREC
      CHARACTER PRGSTR*(6)
      DATA PRGSTR /'WRTIME'/
C begin
      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(ITIMID))
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=500)NSTEP,T
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,'NSTEP,T',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
         NREC = 1
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAMEID(ITIMID), NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NSTEP,T
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'NSTEP,T',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
 500  FORMAT (I15,F15.9)
C end wrtime
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDTIME
C     SUBROUTINE RDTIME(IUNIT,LFORM,NREC,NSTEP,T)
C
C     Reads a timeblock written by L<WRTIME> from file.
COMMEND
      SUBROUTINE RDTIME(IUNIT,LFORM,NREC,NSTEP,T)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NREC,NSTEP
      LOGICAL LFORM
      real*8 T
C local vars
      INTEGER IRES
      CHARACTER PRGSTR*(6)
      DATA PRGSTR /'RDTIME'/
C begin
      IF (LFORM) THEN
         CALL GIMME(IUNIT,PRGSTR,'NSTEP,T')
         CALL CHPINT(NSTEP,PRGSTR,'NSTEP')
         CALL CHPREL(T,PRGSTR,'T')
C get the end marker
         CALL GTEND(IUNIT,PRGSTR)
      ELSE
         CALL CHNREC(PRGSTR,NREC,1)
C
         READ(UNIT = IUNIT, IOSTAT = IRES) NSTEP,T
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'NSTEP,T',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end rdtime
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRIIAV
C     SUBROUTINE WRIIAV(IUNIT,LFORM,NDR,RIIAVE)
C
C     Write a block of time averaged distance restraints
C     to an open file IUNIT. The block can be read using L<RDIIAV>.
COMMEND
      SUBROUTINE WRIIAV(IUNIT,LFORM,NDR,RIIAVE)
Cincludes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NDR
      LOGICAL LFORM
      real*8 RIIAVE(NDR)
C local vars
      INTEGER IRES,I,NREC
      CHARACTER PRGSTR*(6)
      DATA PRGSTR /'WRIIAV'/
C begin
      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IRIIAV))
C write RIIAVE
         DO 10 I=1,NDR
            WRITE(UNIT=IUNIT, IOSTAT=IRES, FMT=500)RIIAVE(I)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,'RIIAVE',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
            IF (MOD(I,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,I)
            ENDIF

 10      CONTINUE

C write end marker
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
         NREC = 1
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAMEID(IRIIAV), NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NDR,RIIAVE
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'RIIAVE',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
 500  FORMAT (E15.7)
C end wriiav
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDIIAV
C     SUBROUTINE RDIIAV(IUNIT,LFORM,NREC,NDR,RIIAVE)
C
C     reads a block of time averaged distance restraint
C     from an open file.
C     The block is written using L<WRIIAV>.
COMMEND
      SUBROUTINE RDIIAV(IUNIT,LFORM,NREC,NDR,RIIAVE)
Cincludes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NREC,NDR
      LOGICAL LFORM
      real*8 RIIAVE(NDR)
C local vars
      INTEGER IRES,I,NUMRED
      CHARACTER PRGSTR*(6)
      CHARACTER STR*(16)
      DATA PRGSTR /'RDIIAV'/
C begin
      IF (LFORM) THEN
C read RIIAVE
         STR = 'ENERGY VAL'
         DO I=1,NDR
            CALL GIMME(IUNIT,PRGSTR,STR)
            CALL CHPREL(RIIAVE(I),PRGSTR,STR)
            CALL CHKMTY(PRGSTR,STR)
         ENDDO

C get the end marker
         CALL GTEND(IUNIT,PRGSTR)
      ELSE
         CALL CHNREC(PRGSTR,NREC,1)
C
         READ(UNIT = IUNIT, IOSTAT = IRES) NUMRED,RIIAVE
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'RIIAVE',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         IF (NUMRED .NE. NDR) THEN
            PRINT 500,NDR,IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
 500  FORMAT (1X,'inconsistency in NDR (RIIAVE) on UNIT',I5)
C end rdiiav
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRJVAV
C     SUBROUTINE WRJVAV(IUNIT,LFORM,NJR,COSQAV,COSIAV)
C
C     Write a block of time averaged J-coupling constant
C     restraints to an open file IUNIT. The block can be read
C     using L<RDJVAV>.
COMMEND
      SUBROUTINE WRJVAV(IUNIT,LFORM,NJR,COSQAV,COSIAV)
Cincludes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NJR
      LOGICAL LFORM
      real*8 COSQAV(NJR),COSIAV(NJR)
C local vars
      INTEGER IRES,I,NREC
      CHARACTER PRGSTR*(6)
      CHARACTER VARSTR*(16)
      DATA PRGSTR /'WRJVAV'/
      DATA VARSTR /'COSQAV, COSIAV'/
C begin
      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(ICOSAV))
C write RJVAVE
         DO 10 I=1,NJR
            WRITE(UNIT=IUNIT, IOSTAT=IRES, FMT=500)COSQAV(I),COSIAV(I)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,VARSTR,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
            IF (MOD(I,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,I)
            ENDIF

 10      CONTINUE

C write end marker
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
         NREC = 1
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAMEID(ICOSAV), NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NJR,COSQAV,COSIAV
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,VARSTR,IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
 500  FORMAT (2E15.7)
C end wrjvav
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDJVAV
C     SUBROUTINE RDJVAV(IUNIT,LFORM,NREC,NJRT,COSQAV,COSIAV)
C
C     Reads a block of time averaged J-value restraints
C     from file. The block must contain NJRT values.
C     The block can be written using L<WRJVAV>.
COMMEND
      SUBROUTINE RDJVAV(IUNIT,LFORM,NREC,NJRT,COSQAV,COSIAV)
Cincludes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NREC,NJRT
      LOGICAL LFORM
      real*8 COSQAV(NJRT),COSIAV(NJRT)
C local vars
      INTEGER IRES,I,NUMRED
      CHARACTER PRGSTR*(6)
      CHARACTER STR*(16)
      DATA PRGSTR /'RDJVAV'/
C begin
      IF (LFORM) THEN
C read RJVAVE
         STR = 'ENERGY VAL'
         DO I=1,NJRT
            CALL GIMME(IUNIT,PRGSTR,STR)
            CALL CHPREL(COSQAV(I),PRGSTR,STR)
            CALL CHPREL(COSIAV(I),PRGSTR,STR)
            CALL CHKMTY(PRGSTR,STR)
         ENDDO

C get the end marker
         CALL GTEND(IUNIT,PRGSTR)
      ELSE
         CALL CHNREC(PRGSTR,NREC,1)
C
         READ(UNIT = IUNIT, IOSTAT = IRES) NUMRED,COSQAV,COSIAV
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'COSQAV,COSIAV',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         IF (NUMRED .NE. NJRT) THEN
            PRINT 500,IUNIT
            PRINT EXPFIL,NJRT
            PRINT BONFIL,NUMRED
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

      ENDIF
 500  FORMAT (1X,'inconsistency in NJR (COSQAV, COSIAV) on UNIT',I5)
C end rdjvav
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRLEME
C     SUBROUTINE WRLEME(IUNIT,LFORM,NDLE,NLECFG,NLEMEM,ILEMEM,NLEVST)
C
C     WRLEME writes the local elevation memory contained in ILEMEM
C     and NLEVST to file. No variables are modified.
C     The block can be read with L<RDLEME>.
C
C     NLEMEM  the number of configurations in the memory    
COMMEND
      SUBROUTINE WRLEME(IUNIT,LFORM,NDLE,NLECFG,NLEMEM,ILEMEM,NLEVST)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'coordsz.h'
C args
      INTEGER IUNIT
      INTEGER NDLE,NLECFG,NLEMEM
      INTEGER ILEMEM(MXLECF,MLECFG)
      INTEGER NLEVST(MLECFG)
      LOGICAL LFORM
C local params
      INTEGER NREC
      PARAMETER (NREC = 3)
C local vars
      INTEGER IRES,I,M
      CHARACTER PRGSTR*(6)
C data
      DATA PRGSTR /'WRLEME'/
C begin
      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(ILEMID))
         WRITE(IUNIT,IOSTAT=IRES,FMT=500)NDLE,NLECFG,NLEMEM
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'local elevation data',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         DO 10 I=1,NLEMEM
            WRITE(IUNIT,IOSTAT=IRES,FMT=500)
     $           NLEVST(I),(ILEMEM(M,I),M=1,NLECFG)
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,'local elevation data',IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
            IF (MOD(I,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,I)
            ENDIF            
 10      CONTINUE

C     write end marker
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
C binary write
         WRITE(UNIT=IUNIT, IOSTAT=IRES)NAMEID(ILEMID),NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT=IUNIT, IOSTAT=IRES)NDLE,NLECFG,NLEMEM
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'NDLE,NLECFG,NLEMEM',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT=IUNIT, IOSTAT=IRES)(NLEVST(I),I=1,NLEMEM)
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'NLEVST',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

         WRITE(UNIT=IUNIT, IOSTAT=IRES)
     $        ((ILEMEM(M,I),M=1,NLECFG),I=1,NLEMEM)
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'ILEMEM',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF

 500  FORMAT (8I10)
C end WRLEME
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDLEME
C      SUBROUTINE RDLEME(IUNIT,LFORM,NREC,
C     $     NDLE,NLECFG,NLEMEM,ILEMEM,NLEVST)
C
C     RDLEME reads in a local elevation memory block.
C     NDLE is checked with the value on file.
C     It is not modified.
C     NDLE must be in the range 1..MAXDLE. This is not
C     checked here, as it is checked on reading the
C     LE dihedral specifications.
C
C     NLECFG is checked with the value on file.
C     It is not modified.
C     NLECFG must be in the range 1..MXLECF
C     This is not checked here, as it is checked on
C     reading the LE dihedral specifications.
C
C     NLEMEM is read from file. This set to
C     the number of configurations read in.
C     NLEMEM will be in the range 1..MLECFG
COMMEND
      SUBROUTINE RDLEME(IUNIT,LFORM,NREC,
     $     NDLE,NLECFG,NLEMEM,ILEMEM,NLEVST)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'formats.h'
C args
      INTEGER IUNIT
      INTEGER NDLE,NLECFG,NLEMEM
      INTEGER ILEMEM(MXLECF,MLECFG)
      INTEGER NLEVST(MLECFG)
      LOGICAL LFORM
      INTEGER NREC
C local vars
      INTEGER IRES,I,M,ITMP1,ITMP2,NDX
      CHARACTER PRGSTR*(6)
      CHARACTER STR*(20)
C data
      DATA PRGSTR /'RDLEME'/
C begin
      IF (LFORM) THEN
         STR = 'NDLE,NLECFG,NLEMEM'
         CALL GIMME(IUNIT,PRGSTR,STR)
C get NDLE
         CALL CHPINT(ITMP1,PRGSTR,STR)
         IF (ITMP1 .NE. NDLE) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 500,'NDLE'

            PRINT EXPFIL,NDLE
            PRINT BONFIL,ITMP1
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

C get NLECFG
         CALL CHPINT(ITMP1,PRGSTR,STR)
         IF (ITMP1 .NE. NLECFG) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 500,'NLECFG'
            PRINT EXPFIL,NLECFG
            PRINT BONFIL,ITMP1
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

C get NLEMEM
         CALL CHPINT(NLEMEM,PRGSTR,STR)
         IF (NLEMEM .LT. 0 .OR. NLEMEM .GT. MLECFG) THEN
            CALL FLAGLN(PRGSTR)
            PRINT FMNII,'NLEMEM',NLEMEM
            PRINT FMBTII,0,'NLEMEM','MLECFG',MLECFG
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         CALL CHKMTY(PRGSTR,STR)

C now get the data
         I = 0
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C while loop
 10      IF (NDX .EQ. IDUKN .AND. I .LT. NLEMEM) THEN
            I = I+1
            CALL CHPINT(NLEVST(I),'NLEVST)I)',PRGSTR)
            DO 20 M=1,NLECFG
               CALL CHPINT(ILEMEM(M,I),'ILEMEM',PRGSTR)
 20         CONTINUE
            CALL CHKMTY(PRGSTR,'ILEMEM')
C
            CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
            GOTO 10
C end of while
         ENDIF

         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN(PRGSTR)
            PRINT *,'attempting to read NLEMEM = ',
     $           NLEMEM,'configurations'
            PRINT *,ENDEXP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ELSE
C binary read
         CALL CHNREC(PRGSTR,NREC,3)
C get NDLE,NLECFG,NLEMEM
         READ(UNIT=IUNIT, IOSTAT=IRES)ITMP1,ITMP2,NLEMEM
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'NDLE,NLECFG,NLEMEM',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF   
C check the values
         IF (ITMP1 .NE. NDLE) THEN
            PRINT 500,'NDLE'
            PRINT EXPFIL,NDLE
            PRINT BONFIL,ITMP1
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         IF (ITMP2 .NE. NLECFG) THEN
            PRINT 500,'NLECFG'
            PRINT EXPFIL,NLECFG
            PRINT BONFIL,ITMP2
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         IF (NLEMEM .LT. 0 .OR. NLEMEM .GT. MLECFG) THEN
            CALL FLAGLN(PRGSTR)
            PRINT FMNII,'NLEMEM',NLEMEM
            PRINT FMBTII,0,'NLEMEM','MLECFG',MLECFG
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

C get NLEVST
         READ(UNIT=IUNIT,IOSTAT=IRES)(NLEVST(I),I=1,NLEMEM)
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'NLEVST',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF

C get ILEMEM
         READ(UNIT=IUNIT, IOSTAT=IRES)
     $        ((ILEMEM(M,I),M=1,NLECFG),I=1,NLEMEM)
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'ILEMEM',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
 500  FORMAT (1X,'inconsistent input for ',A)
C end rdleme
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRP2M
C     SUBROUTINE WRP2M(IUNIT,LFORM,NAT,NPM,NSM,XX)
C
C     Writes a second moment block to file.
C     The block can be read using L<RDP2M>.
COMMEND
      SUBROUTINE WRP2M(IUNIT,LFORM,NAT,NPM,NSM,XX)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'coordsz.h'
C args
      INTEGER IUNIT,NAT,NPM,NSM
      LOGICAL LFORM
      real*8 XX(M2MAX,NAT)
C begin
      CALL DMPMXN(IUNIT,LFORM,IP2MID,'WRP2M',M2MAX,NAT,NPM,NSM,XX)
C end wrp2m
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDP2M
C     SUBROUTINE RDP2M(IUNIT,LFORM,NREC,NAT,XX)
C
C     Reads a second moment block from file.
C     The block can be written using L<WRP2M>.
COMMEND
      SUBROUTINE RDP2M(IUNIT,LFORM,NREC,NAT,XX)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'coordsz.h'
C args
      INTEGER IUNIT,NAT,NREC
      LOGICAL LFORM
      real*8 XX(M2MAX,NAT)
C begin
      CALL GETMXN(IUNIT,LFORM,NREC,'RDP2M',M2MAX,NAT,XX)
C end rdp2m
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRP3M
C     SUBROUTINE WRP3M(IUNIT,LFORM,NAT,NPM,NSM,XX)
C
C     Writes a second moment block to file.
C     The block can be read using L<RDP3M>.
COMMEND
      SUBROUTINE WRP3M(IUNIT,LFORM,NAT,NPM,NSM,XX)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'coordsz.h'
C args
      INTEGER IUNIT,NAT,NPM,NSM
      LOGICAL LFORM
      real*8 XX(M3MAX,NAT)
C begin
      CALL DMPMXN(IUNIT,LFORM,IP3MID,'WRP3M',M3MAX,NAT,NPM,NSM,XX)
C end wrp3m
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDP3M
C     SUBROUTINE RDP3M(IUNIT,LFORM,NREC,NAT,XX)
C
C     Reads a second moment block from file.
C     The block can be written using L<WRP3M>.
COMMEND
      SUBROUTINE RDP3M(IUNIT,LFORM,NREC,NAT,XX)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'coordsz.h'
C args
      INTEGER IUNIT,NAT,NREC
      LOGICAL LFORM
      real*8 XX(M3MAX,NAT)
C begin
      CALL GETMXN(IUNIT,LFORM,NREC,'RDP3M',M3MAX,NAT,XX)
C end rdp3m
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRP4M
C     SUBROUTINE WRP4M(IUNIT,LFORM,NAT,NPM,NSM,XX)
C
C     Writes a second moment block to file.
C     The block can be read using L<RDP4M>.
COMMEND
      SUBROUTINE WRP4M(IUNIT,LFORM,NAT,NPM,NSM,XX)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'coordsz.h'
C args
      INTEGER IUNIT,NAT,NPM,NSM
      LOGICAL LFORM
      real*8 XX(M3MAX,NAT)
C begin
      CALL DMPMXN(IUNIT,LFORM,IP4MID,'WRP4M',M4MAX,NAT,NPM,NSM,XX)
C end wrP4M
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDP4M
C     SUBROUTINE RDP4M(IUNIT,LFORM,NREC,NAT,XX)
C
C     Reads a second moment block from file.
C     The block can be written using L<WRP4M>.
COMMEND
      SUBROUTINE RDP4M(IUNIT,LFORM,NREC,NAT,XX)
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'coordsz.h'
C args
      INTEGER IUNIT,NAT,NREC
      LOGICAL LFORM
      real*8 XX(M3MAX,NAT)
C begin
      CALL GETMXN(IUNIT,LFORM,NREC,'RDP4M',M4MAX,NAT,XX)
C end rdp4m
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRP2MT
C     SUBROUTINE WRP2MT(IUNIT,LFORM,NAT,NPM,NSM,XTR)
C
C     Writes a trace of second moment block to file : XTR(1..NAT).
C     The block can be read using L<RDP2MT>.
COMMEND
      SUBROUTINE WRP2MT(IUNIT,LFORM,NAT,NPM,NSM,XTR)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NAT,NPM,NSM
      LOGICAL LFORM
      real*8 XTR(NAT)
C begin
      CALL DMPMXN(IUNIT,LFORM,IP2TID,'WRP2MT',1,NAT,NPM,NSM,XTR)
C end wrp2mt
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDP2MT
C     SUBROUTINE RDP2MT(IUNIT,LFORM,NREC,NAT,XTR)
C
C     Reads a trace of second moment block from file : XTR(1..NAT).
C     The block can be written using L<WRP2MT>.
COMMEND
      SUBROUTINE RDP2MT(IUNIT,LFORM,NREC,NAT,XTR)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NAT,NREC
      LOGICAL LFORM
      real*8 XTR(NAT)
C begin
      CALL GETMXN(IUNIT,LFORM,NREC,'RDP2MT',1,NAT,XTR)
C end rdp2mt
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRTCF
C     SUBROUTINE WRTCF(IUNIT,IMODE,LFORM,MAXNQV,NRQV,NQ,NTVC,Q,TC,FC,
C    $     NSPEC,WC,FAH,NDIS,DIS,QB,COMSTR)
C
C     Writes a either
C        IMODE=1, a time series of quantity block to file 
C                 TC(1..NRQV), Q(1..NRQV,NQ,1..NTVC), NQ=1,2
C        IMODE=2, a frequency distribution block to file
C                 QB(1..NDIS), DIS(1..NDIS)
C        IMODE=3, a time correlation function block to file
C                 TC(1..NRQV), FC(1..NRQV)
C        IMODE=4, a spectral density block to file
C                 WC(1..NSPEC), FAH(1..NSPEC)
C
C     The blocks (IMODE=1) can be read using L<RDTCF>.
COMMEND
      SUBROUTINE WRTCF(IUNIT,IMODE,LFORM,MAXNQV,NRQV,NQ,NTVC,Q,TC,FC,
     $     NSPEC,WC,FAH,NDIS,DIS,QB,COMSTR)
C includes
      INCLUDE 'dataid.h'
C local param statements
      INTEGER MYID
      INTEGER NREC
      PARAMETER (NREC = 3)
C args
      INTEGER IUNIT,IMODE,MAXNQV,NRQV,NQ,NTVC,NDIS,NSPEC
      LOGICAL LFORM
      real*8 Q(MAXNQV,2,3), TC(NRQV), FC(NRQV), WC(NSPEC), FAH(NSPEC),
     $     DIS(NDIS), QB(NDIS)
      CHARACTER*(*) COMSTR
C local vars
      INTEGER IRES,N,M
      CHARACTER PRGSTR*(5)
      CHARACTER*9 STR0
      CHARACTER*4 STR00
      CHARACTER*33 STR12
      CHARACTER*21 STR1
      CHARACTER*11 STR2

      CHARACTER*12 STR3
      CHARACTER*11 STR4
      CHARACTER*24 STR34

      CHARACTER*11 STR5
      CHARACTER*23 STR56

      CHARACTER*12 STR7
      CHARACTER*11 STR8
      CHARACTER*24 STR78
C data
      DATA PRGSTR /'WRTCF'/
      DATA STR0 /'NRQV,NTVC'/
      DATA STR78 /'WC(1..NRQV),FAH(1..NRQV)'/
      DATA STR8 /'WC(1..NRQV)'/
      DATA STR7 /'FAH(1..NRQV)'/
      DATA STR56 /'TC(1..NRQV),FC(1..NRQV)'/
      DATA STR5 /'FC(1..NRQV)'/
      DATA STR34 /'QB(1..NDIS),DIS(1..NDIS)'/
      DATA STR4 /'QB(1..NDIS)'/
      DATA STR3 /'DIS(1..NDIS)'/
      DATA STR2 /'TC(1..NRQV)'/
      DATA STR1 /'Q(1..NRQV,NQ,1..NTVC)'/
      DATA STR12 /'TC(1..NRQV),Q(1..NRQV,NQ,1..NTVC)'/
      DATA STR00 /'NDIS'/
C begin
      IF (IMODE.EQ.1) THEN
         MYID = IQTSID
      ELSEIF (IMODE.EQ.2 ) THEN
         MYID = IQDSID
      ELSEIF (IMODE.EQ.3 ) THEN
         MYID = IQTCID
      ELSEIF (IMODE.EQ.4 ) THEN
         MYID = IQSDID
      ELSE
         PRINT *,PRGSTR,' wrong IMODE = ',IMODE
      ENDIF
C
      IF (LFORM) THEN
C formatted write - write header
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(MYID))
         CALL WRCOMM(IUNIT,PRGSTR,COMSTR)
C
         IF (MYID.EQ.IQTSID) THEN
C write time series
            WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=500)
     $           (TC(N),(Q(N,NQ,M),M=1,NTVC),N=1,NRQV)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,STR12,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
         ELSEIF (MYID.EQ.IQDSID) THEN
C write frequency distribution
            WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=500)
     $           (QB(N),DIS(N),N=1,NDIS)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,STR34,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
         ELSEIF (MYID.EQ.IQTCID) THEN
C write time correlation function
            WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=500)
     $           (TC(N),FC(N),N=1,NRQV)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,STR56,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
         ELSEIF (MYID.EQ.IQSDID) THEN
C write spectral density
            WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=500)
     $           (WC(N),FAH(N),N=1,NSPEC)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,STR78,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
         ENDIF
C write end marker
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
C binary write - write header
         WRITE(UNIT=IUNIT, IOSTAT=IRES) NAMEID(MYID),NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,NAMEID(MYID),IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         IF (MYID.EQ.IQTSID) THEN
C write time series
            WRITE(UNIT=IUNIT, IOSTAT=IRES) NRQV,NTVC
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,STR0,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            WRITE(UNIT=IUNIT, IOSTAT=IRES)
     $           ((Q(N,NQ,M),N=1,NRQV),M=1,NTVC)
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,STR1,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            WRITE(UNIT=IUNIT, IOSTAT=IRES) (TC(N),N=1,NRQV)
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,STR2,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF          
C
         ELSEIF (MYID.EQ.IQDSID) THEN
C write frequency distribution
            WRITE(UNIT=IUNIT, IOSTAT=IRES) NDIS
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,STR00,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            WRITE(UNIT=IUNIT, IOSTAT=IRES) (DIS(N),N=1,NDIS)
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,STR3,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            WRITE(UNIT=IUNIT, IOSTAT=IRES) (QB(N),N=1,NDIS)
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,STR4,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
         ELSEIF (MYID.EQ.IQTCID) THEN
C write time correlation function
            WRITE(UNIT=IUNIT, IOSTAT=IRES) NRQV
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,STR0,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            WRITE(UNIT=IUNIT, IOSTAT=IRES) (FC(N),N=1,NRQV)
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,STR5,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            WRITE(UNIT=IUNIT, IOSTAT=IRES) (TC(N),N=1,NRQV)
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,STR2,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
         ELSEIF (MYID.EQ.IQSDID) THEN
C write spectral density
            WRITE(UNIT=IUNIT, IOSTAT=IRES) NRQV
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,STR0,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            WRITE(UNIT=IUNIT, IOSTAT=IRES) (FAH(N),N=1,NSPEC)
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,STR7,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            WRITE(UNIT=IUNIT, IOSTAT=IRES) (WC(N),N=1,NSPEC)
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,STR8,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
         ENDIF
      ENDIF
 500  FORMAT (4F20.9)
C end wrtcf
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDTCF
C     SUBROUTINE RDTCF(IUNIT,IMODE,LFORM,NREC,MAXNQV,
C    $     NRQV,NOFF1,NQ,NTVC,Q,TC)
C
C     Reads a time series of quantity block from file,
C     Q(1..NRQV,NQ,1..NTVC), NQ=1,2
C     The read values will be offset by NOFF1 in their first index
C     The block can be written using L<WRTCF>.
C
C     IMODE is currently unused - no program needs to read
C     the other blocks generated by WRTCF
C
COMMEND
      SUBROUTINE RDTCF(IUNIT,IMODE,LFORM,NREC,MAXNQV,
     $     NRQV,NOFF1,NQ,NTVC,Q,TC)
C includes
      INCLUDE 'dataid.h'
C params
C args
      INTEGER IUNIT,IMODE,NREC,MAXNQV,NRQV,NOFF1,NQ,NTVC
      LOGICAL LFORM
      real*8 Q(MAXNQV,2,3), TC(NRQV)
C local vars
      INTEGER NDX,IRES,MYID,N,M
      CHARACTER PRGSTR*(5)
      CHARACTER*9 STR0
      CHARACTER*33 STR12
      CHARACTER*21 STR1
      CHARACTER*11 STR2
C data
      DATA STR1 /'Q(1..NRQV,NQ,1..NTVC)'/
      DATA STR12 /'TC(1..NRQV),Q(1..NRQV,NQ,1..NTVC)'/
      DATA STR0 /'NRQV,NTVC'/
      DATA PRGSTR /'RDTCF'/
      DATA STR2 /'TC(1..NRQV)'/
C
C begin
C
      MYID = IQTSID
      IF (LFORM) THEN
C
C formatted read - read until we get an END or NRQ is exceeded
C
C QUANTIMESERIES
C
         CALL GIMME(IUNIT,PRGSTR,NAMEID(MYID))
         DO 10 N=NOFF1+1,NOFF1+NRQV
            CALL CHPREL(TC(N),PRGSTR,STR2)
            DO 20 M=1,NTVC
               CALL CHPREL(Q(N,NQ,M),PRGSTR,STR1)
 20         CONTINUE
 10      CONTINUE
         CALL GNDXID(IUNIT,PRGSTR,STR12,NDX)
         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN(PRGSTR)
            PRINT EXPFIL,NRQV
            PRINT *,ENDEXP
            PRINT *,'with NTVC= ',NTVC
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
      ELSE
C
C unformatted read - check NREC, read NRQV,NTVC and blocks
         CALL CHNREC(PRGSTR,NREC,3)
C
C QUANTIMESERIES
C
         READ(UNIT = IUNIT, IOSTAT = IRES) NRQV,NTVC
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,STR0,IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         READ(UNIT = IUNIT, IOSTAT = IRES) 
     $        ((Q(N,NQ,M),N=NOFF1+1,NOFF1+NRQV),M=1,NTVC)
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,STR1,IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         READ(UNIT = IUNIT, IOSTAT = IRES) 
     $        (TC(N),N=NOFF1+1,NOFF1+NRQV)
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,STR2,IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C end rdtcf
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRAVQ
C     SUBROUTINE WRAVQ(IUNIT,IMODE,LFORM,NRQ,AQ,AE,COMSTR)
C
C     Writes a either
C        IMODE=1, a quantity average block to file AQ(6,1,NRQ)
C        IMODE=2, a quantity energy average block to file AQ(6,2,NRQ)
C        IMODE=3, a quantity sum of energies block AE(6)
C     The block can be read using L<RDAVQ>.
COMMEND
      SUBROUTINE WRAVQ(IUNIT,IMODE,LFORM,NRQ,AQ,AE,COMSTR)
C includes
      INCLUDE 'dataid.h'
C local param statements
      INTEGER MYID
      INTEGER NREC
      PARAMETER (NREC = 2)
C args
      INTEGER IUNIT,IMODE,NRQ
      LOGICAL LFORM
      real*8 AQ(6,2,NRQ), AE(6)
      CHARACTER*(*) COMSTR
C local vars
      INTEGER MQ,K,IRES,N
      CHARACTER PRGSTR*(5)
      CHARACTER*3 STR0
      CHARACTER*11 STR1
      CHARACTER*5 STR2
C data
      DATA PRGSTR /'WRAVQ'/
      DATA STR0 /'NRQ'/
      DATA STR1 /'AQ(6,2,NRQ)'/
      DATA STR2 /'AE(6)'/
C begin
      IF (IMODE.EQ.1) THEN
         MYID = IQTAID
         MQ = 1
      ELSEIF (IMODE.EQ.2) THEN
         MYID = IQEAID
         MQ = 2
      ELSEIF (IMODE.EQ.3) THEN
         MYID = IQESID
         MQ = 0
      ELSE
         PRINT *,PRGSTR,' wrong IMODE = ',IMODE
      ENDIF
C
      IF (LFORM) THEN
C formatted write - write header
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(MYID))
         CALL WRCOMM(IUNIT,PRGSTR,COMSTR)
C write averages
         IF (MQ.NE.0) THEN
            WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=500)
     $           ((AQ(K,MQ,N),K=1,6),N=1,NRQ)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,STR1,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
         ELSE
            WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=500) (AE(K),K=1,6)
            IF (IRES .NE. 0) THEN
               PRINT FWFAIL,PRGSTR,STR2,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
         ENDIF
C write end marker
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
C binary write - write header
         WRITE(UNIT=IUNIT, IOSTAT=IRES) NAMEID(MYID),NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,NAMEID(MYID),IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         IF (MQ.NE.0) THEN
            WRITE(UNIT=IUNIT, IOSTAT=IRES) NRQ
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,STR0,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            WRITE(UNIT=IUNIT, IOSTAT=IRES)
     $           ((AQ(K,MQ,N),K=1,6),N=1,NRQ)
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,STR1,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
         ELSE
            WRITE(UNIT=IUNIT, IOSTAT=IRES) (AE(K),K=1,6)
            IF (IRES .NE. 0)THEN
               PRINT FWFAIL,PRGSTR,STR2,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
         ENDIF
      ENDIF
C increased precision, testing, Phil
 500  FORMAT (6E13.5)
C end wravq
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDAVQ
C     SUBROUTINE RDAVQ(IUNIT,IMODE,LFORM,NREC,NRQ,AQ,AE)
C
C     Reads either
C        IMODE=1, a quantity average block from file - AQ(6,1,NRQ)
C        IMODE=2, a quantity energy average block from file - AQ(6,2,NRQ)
C        IMODE=3, a quantity sum of energies block from file - AE(6)
C     The block can be written using L<WRAVQ>.
C
COMMEND
      SUBROUTINE RDAVQ(IUNIT,IMODE,LFORM,NREC,NRQ,AQ,AE)
C includes
      INCLUDE 'dataid.h'
C params
C args
      INTEGER IUNIT,IMODE,NREC,NRQ
      LOGICAL LFORM
      real*8 AQ,AE
      DIMENSION AQ(6,2,NRQ),AE(6)
C local vars
      INTEGER NDX,IRES,MYID,MQ,I,J
      CHARACTER PRGSTR*(5)
      CHARACTER*3 STR0
      CHARACTER*11 STR1
      CHARACTER*5 STR2
C data
      DATA PRGSTR /'RDAVQ'/
      DATA STR0 /'NRQ'/
      DATA STR1 /'AQ(6,2,NRQ)'/
      DATA STR2 /'AE(6)'/
C begin - select mode
C
      IF (IMODE.EQ.1) THEN
         MYID = IQTAID
         MQ = 1
      ELSEIF (IMODE.EQ.2) THEN
         MYID = IQEAID
         MQ = 2
      ELSEIF (IMODE.EQ.3) THEN
         MYID = IQESID
         MQ = 0
      ELSE
         PRINT *,PRGSTR,' wrong IMODE = ',IMODE
      ENDIF
C
      IF (LFORM) THEN
C
C formatted read - read until we get an END or NRQ is exceeded
C
         IF ( MQ.NE.0 ) THEN
C
C QUANTITYAVER, QUANENEAVER
C
            CALL GIMME(IUNIT,PRGSTR,NAMEID(MYID))
            DO I=1,NRQ
               CALL CHPREL(AQ(1,MQ,I),PRGSTR,STR1)
               CALL CHPREL(AQ(2,MQ,I),PRGSTR,STR1)
               CALL CHPREL(AQ(3,MQ,I),PRGSTR,STR1)
               CALL CHPREL(AQ(4,MQ,I),PRGSTR,STR1)
               CALL CHPREL(AQ(5,MQ,I),PRGSTR,STR1)
               CALL CHPREL(AQ(6,MQ,I),PRGSTR,STR1)
            ENDDO
            CALL GNDXID(IUNIT,PRGSTR,STR1,NDX)
            IF (NDX .NE. IENDID) THEN
               CALL FLAGLN(PRGSTR)
               PRINT EXPFIL,NRQ
               PRINT *,ENDEXP
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
C
         ELSE
C
C QUANSUMENEAVER
C
            CALL GIMME(IUNIT,PRGSTR,NAMEID(MYID))
            CALL CHPREL(AE(1),PRGSTR,STR2)
            CALL CHPREL(AE(2),PRGSTR,STR2)
            CALL CHPREL(AE(3),PRGSTR,STR2)
            CALL CHPREL(AE(4),PRGSTR,STR2)
            CALL CHPREL(AE(5),PRGSTR,STR2)
            CALL CHPREL(AE(6),PRGSTR,STR2)
            CALL GNDXID(IUNIT,PRGSTR,STR2,NDX)
            IF (NDX .NE. IENDID) THEN
               CALL FLAGLN(PRGSTR)
               PRINT EXPFIL,NRQ
               PRINT *,ENDEXP
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
         ENDIF
C
      ELSE
C
C unformatted read - check NREC, read NRQ and blocks
C
         CALL CHNREC(PRGSTR,NREC,2)
         IF ( MQ.NE.0 ) THEN
C QUANTITYAVER, QUANENEAVER
            READ(UNIT = IUNIT, IOSTAT = IRES) NRQ
            IF (IRES .NE. 0)THEN
               PRINT FRFAIL,PRGSTR,STR0,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
            READ(UNIT = IUNIT, IOSTAT = IRES)
     $           ((AQ(J,MQ,I),J=1,6),I=1,NRQ)
            IF (IRES .NE. 0)THEN
               PRINT FRFAIL,PRGSTR,STR1,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
         ELSE
C
C QUANSUMENEAVER
C
            READ(UNIT = IUNIT, IOSTAT = IRES) (AE(J),J=1,6)
            IF (IRES .NE. 0)THEN
               PRINT FRFAIL,PRGSTR,STR2,IUNIT
               CALL STOPGM('BLOCKIO','FATAL ERROR')
            ENDIF
         ENDIF
      ENDIF
C end rdavq
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRSVS
C     SUBROUTINE WRSVS(IUNIT,LFORM,NSVT,NMOL,NFRAME,NBIN,STATR,COMSTR)
C
C     Writes a solvent statistics average block to file.
C     The block can be read using L<RDSVS>.
COMMEND
      SUBROUTINE WRSVS(IUNIT,LFORM,NSVT,NMOL,NFRAME,NBIN,STATR,COMSTR)
C includes
      INCLUDE 'dataid.h'
C local param statements
      INTEGER MYID
      PARAMETER (MYID = ISVSID)
      INTEGER NREC
      PARAMETER (NREC = 2)
      INTEGER MAXSTA
      PARAMETER (MAXSTA=19)
C args
      INTEGER IUNIT,NSVT,NMOL,NFRAME,NBIN
      LOGICAL LFORM
      real*8 STATR
      DIMENSION STATR(MAXSTA,NBIN)
      CHARACTER*(*) COMSTR
C local vars
      INTEGER K,IRES,N
      CHARACTER PRGSTR*(5)
      DATA PRGSTR /'WRSVS'/
C begin
      IF (LFORM) THEN
C formatted write,
C write header
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(MYID))
         CALL WRCOMM(IUNIT,PRGSTR,COMSTR)
C write solvent molecule number in the list, molecule number in the 
C coordinates, number of frames, number of bins 
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=24)
     $        NSVT,NMOL,NFRAME,NBIN
C     write averages
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=25)
     $        ((STATR(K,N),K=1,MAXSTA),N=1,NBIN)
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,'STATR',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C write end marker
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
C binary write
C write header
         WRITE(UNIT=IUNIT, IOSTAT=IRES) NAMEID(MYID),NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C NSVT,NMOL,NFRAME,NBIN
         WRITE(UNIT=IUNIT, IOSTAT=IRES) NSVT,NMOL,NFRAME,NBIN
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'NSVT,NMOL,NFRAME,NBIN',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT=IUNIT, IOSTAT=IRES) STATR
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'STATR',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C     
 24   FORMAT(4I10)
 25   FORMAT(6E15.7)
C end wrsvs
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDSVS
C     SUBROUTINE RDSVS(IUNIT,LFORM,NREC,NSVT,NMOL,NFRAME,NBIN,STATR)
C
C     Reads a solvent statistics block from file.
C     The block can be written using L<WRSVS>.
COMMEND
      SUBROUTINE RDSVS(IUNIT,LFORM,NREC,NSVT,NMOL,NFRAME,NBIN,STATR)
C includes
      INCLUDE 'dataid.h'
C params
      INTEGER MAXSTA
      PARAMETER (MAXSTA=19)
C args
      INTEGER IUNIT,NREC,NSVT,NMOL,NFRAME,NBIN,IRES
      LOGICAL LFORM
      real*8 STATR
      DIMENSION STATR(MAXSTA,NBIN)
C local vars
      INTEGER NDX,NBINL,K,I
      CHARACTER*60 STR
      CHARACTER PRGSTR*(5)
      DATA PRGSTR /'RDSVS'/
C begin
      IF (LFORM) THEN
C formatted read: read until we get an END or
C NBIN is exceeded
         NBINL = 0
         STR = 'SOLVSTAT BLOCK'
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C read NSVT,NMOL,NFRAME,NBIN
         STR = 'NSVT,NMOL,NFRAME,NBIN'
         CALL CHPINT(NSVT,PRGSTR,STR)
         CALL CHPINT(NMOL,PRGSTR,STR)
         CALL CHPINT(NFRAME,PRGSTR,STR)
         CALL CHPINT(NBIN,PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)

C while no marker  AND still atoms to read
         STR = 'STATR'
         CALL GIMME(IUNIT,PRGSTR,STR)
         DO 10 I=1,NBIN
            DO 80 K=1,MAXSTA
               CALL CHPREL(STATR(K,I),PRGSTR,STR)
 80         CONTINUE
 10      CONTINUE
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)

         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 530,NBIN
            PRINT *,ENDEXP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ELSE
C unformatted reading
         CALL CHNREC(PRGSTR,NREC,2)
         READ(UNIT = IUNIT, IOSTAT = IRES) NSVT,NMOL,NFRAME,NBIN
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'NSVT,NMOL,NFRAME,NBIN',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         READ(UNIT = IUNIT, IOSTAT = IRES) STATR
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'data',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF

 530  FORMAT (1X,'attempting to read',I7,' solvent statistics bins')
C end rdsvs
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRDFS
C     SUBROUTINE WRDFS(IUNIT,LFORM,NSVT,NMOL,NFRAME,NCRTC,MAXCRT,XSVTD,
C           XSVTX,ROTR,TRANSR,XINITX,COMSTR)
C
C     Writes a diffusion statistics average block to file.
C     The block can be read using L<RDDFS>.
COMMEND
      SUBROUTINE WRDFS(IUNIT,LFORM,NSVT,NMOL,NFRAME,NCRTC,MAXCRT,XSVTD,
     $     XSVTX,ROTR,TRANSR,XINITX,COMSTR)
C includes
      INCLUDE 'dataid.h'
C local param statements
      INTEGER MYID
      PARAMETER (MYID = IDFSID)
      INTEGER NREC
      PARAMETER (NREC = 2)
      INTEGER MAXLEG
      PARAMETER (MAXLEG=3)
C args
      INTEGER IUNIT,NSVT,NMOL,NFRAME,NCRTC,MAXCRT
      LOGICAL LFORM
      real*8 XSVTD,XSVTX,ROTR,TRANSR,XINITX
      DIMENSION XSVTD(MAXCRT,3),XSVTX(MAXCRT,3),ROTR(MAXLEG,NCRTC),
     $     TRANSR(NCRTC),XINITX(3)
      CHARACTER*(*) COMSTR
C local vars
      INTEGER K,IRES,N
      CHARACTER*60 STR
      CHARACTER PRGSTR*(5)
C data
      DATA PRGSTR /'WRDFS'/
C begin
      IF (LFORM) THEN
C formatted write,
C write header
         STR = 'DIFFSTAT BLOCK'
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(MYID))
         CALL WRCOMM(IUNIT,PRGSTR,COMSTR)

C write solvent molecule number in the list, molecule number in the 
C coordinates, number of frames, number of frames to correlate
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=24)
     $        NSVT,NMOL,NFRAME,NCRTC
C write averages
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=25)
     $        ((XSVTD(K,N),K=1,NCRTC),N=1,3)
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,'XSVTD',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=25)
     $        ((XSVTX(K,N),K=1,NCRTC),N=1,3)
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,'XSVTX',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=25)
     $        ((ROTR(K,N),K=1,MAXLEG),N=1,NCRTC)
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,'ROTR',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=25)
     $        (TRANSR(K),K=1,NCRTC)
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,'TRANSR',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=25)
     $        (XINITX(K),K=1,3)
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,'XINITX',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C write end marker
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
C binary write
C write header
         WRITE(UNIT=IUNIT, IOSTAT=IRES) NAMEID(MYID),NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C NSVT,NMOL,NFRAME,NCRTC
         WRITE(UNIT=IUNIT, IOSTAT=IRES) NSVT,NMOL,NFRAME,NCRTC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'NSVT,NMOL,NFRAME,NCRTC',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT=IUNIT, IOSTAT=IRES) XSVTD,XSVTX,ROTR,
     $        TRANSR,XINITX
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'XSVTD, ETC',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C     
 24   FORMAT(4I10)
 25   FORMAT(6E15.7)
C end wrdfs
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDDFS
C     SUBROUTINE RDDFS(IUNIT,LFORM,NREC,NSVT,NMOL,NFRAME,MAXCRT,NCRTC,
C                      XSVTD,XSVTX,ROTR,TRANSR,XINITX)
C
C     Reads a diffusion statistics block from file.
C     The block can be written using L<WRDFS>.
COMMEND
      SUBROUTINE RDDFS(IUNIT,LFORM,NREC,NSVT,NMOL,NFRAME,MAXCRT,NCRTC,
     $     XSVTD,XSVTX,ROTR,TRANSR,XINITX)
C includes
      INCLUDE 'dataid.h'
C params
      INTEGER MAXLEG
      PARAMETER (MAXLEG=3)
C args
      INTEGER IUNIT,NREC,NSVT,NMOL,NFRAME,NCRTC,MAXCRT,IRES
      LOGICAL LFORM
      real*8 XSVTD,XSVTX,ROTR,TRANSR,XINITX
      DIMENSION XSVTD(MAXCRT,3),XSVTX(MAXCRT,3),ROTR(MAXLEG,NCRTC),
     $     TRANSR(NCRTC),XINITX(3)
C local vars
      INTEGER NDX,NCRTCL,K,I
      CHARACTER*60 STR
      CHARACTER PRGSTR*(5)
C data
      DATA PRGSTR /'RDDFS'/
C begin
      IF (LFORM) THEN
C formatted read: read until we get an END or
C NBIN is exceeded
         NCRTCL = 0
         STR = 'DIFFSTAT BLOCK'
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C read NSVT,NMOL,NFRAME,NBIN
         STR = 'NSVT,NMOL,NFRAME,NCRTC'
         CALL CHPINT(NSVT,PRGSTR,STR)
         CALL CHPINT(NMOL,PRGSTR,STR)
         CALL CHPINT(NFRAME,PRGSTR,STR)
         CALL CHPINT(NCRTC,PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)

C while no marker  AND still atoms to read
         STR = 'XSVTD'
         CALL GIMME(IUNIT,PRGSTR,STR)
         DO I=1,NCRTC
            CALL CHPREL(XSVTD(I,1),PRGSTR,STR)
         ENDDO
         DO I=1,NCRTC
            CALL CHPREL(XSVTD(I,2),PRGSTR,STR)
         ENDDO
         DO I=1,NCRTC
            CALL CHPREL(XSVTD(I,3),PRGSTR,STR)
         ENDDO

C
         STR = 'XSVTX'
         CALL GIMME(IUNIT,PRGSTR,STR)
         DO I=1,NCRTC
            CALL CHPREL(XSVTX(I,1),PRGSTR,STR)
         ENDDO
         DO I=1,NCRTC
            CALL CHPREL(XSVTX(I,2),PRGSTR,STR)
         ENDDO
         DO I=1,NCRTC
            CALL CHPREL(XSVTX(I,3),PRGSTR,STR)
         ENDDO
C
         STR = 'ROTR'
         CALL GIMME(IUNIT,PRGSTR,STR)
         DO I=1,NCRTC
            DO K=1,MAXLEG
               CALL CHPREL(ROTR(K,I),PRGSTR,STR)
            ENDDO
         ENDDO
C
         STR = 'TRANSR'
         CALL GIMME(IUNIT,PRGSTR,STR)
         DO 50 I=1,NCRTC
            CALL CHPREL(TRANSR(I),PRGSTR,STR)
 50      CONTINUE
C
         STR = 'XINITX'
         CALL GIMME(IUNIT,PRGSTR,STR)
         CALL CHPREL(XINITX(1),PRGSTR,STR)
         CALL CHPREL(XINITX(2),PRGSTR,STR)
         CALL CHPREL(XINITX(3),PRGSTR,STR)
C
         STR = 'DIFFSTAT BLOCK'
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)

         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 530,NCRTC
            PRINT *,ENDEXP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ELSE
C unformatted reading
         CALL CHNREC(PRGSTR,NREC,2)
         READ(UNIT = IUNIT, IOSTAT = IRES) NSVT,NMOL,NFRAME,NCRTC
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'NSVT,NMOL,NFRAME,NCRTC',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         READ(UNIT = IUNIT, IOSTAT = IRES) XSVTD,XSVTX,ROTR,
     $        TRANSR,XINITX
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'XSVTD, ETC.',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C
 530  FORMAT (1X,'attempting to read',I7,' solvent statistics bins')
C end rddfs
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRDMS
C     SUBROUTINE WRDMS(IUNIT,LFORM,NSVT,NMOL,NFRAME,NCRTC,MAXCRT,
C                      NBINC,MAXBNC,XMDIP,DIPC,COMSTR)
C
C     Writes a dipole moment statistics average block to file.
C     The block can be read using L<RDDMS>.
COMMEND
      SUBROUTINE WRDMS(IUNIT,LFORM,NSVT,NMOL,NFRAME,NCRTC,MAXCRT,
     $     NBINC,MAXBNC,XMDIP,DIPC,COMSTR)
C includes
      INCLUDE 'dataid.h'
C local param statements
      INTEGER MYID
      PARAMETER (MYID = IDMSID )
      INTEGER NREC
      PARAMETER (NREC = 2)
C args
      INTEGER IUNIT,NSVT,NMOL,NFRAME,NCRTC,MAXCRT,NBINC,MAXBNC
      LOGICAL LFORM
      real*8 XMDIP,DIPC
      DIMENSION XMDIP(MAXCRT,MAXBNC,3),DIPC(MAXCRT,NBINC)
      CHARACTER*(*) COMSTR
C local vars
      INTEGER K,IRES,N,I
      CHARACTER*60 STR
      CHARACTER PRGSTR*(5)
      DATA PRGSTR /'WRDMS'/
C begin
      IF (LFORM) THEN
C formatted write,
C write header
         STR = 'DIPMSTAT BLOCK'
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(MYID))
         CALL WRCOMM(IUNIT,PRGSTR,COMSTR)
C write solvent molecule number in the list, molecule number in the 
C coordinates, number of frames, number of frames to correlate
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=24)
     $        NSVT,NMOL,NFRAME,NCRTC,NBINC
C write averages
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=25)
     $        (((XMDIP(I,K,N),I=1,NCRTC),K=1,NBINC),N=1,3)
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,'XMDIP',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
         WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT=25)
     $        ((DIPC(I,K),I=1,NCRTC),K=1,NBINC)
         IF (IRES .NE. 0) THEN
            PRINT FWFAIL,PRGSTR,'DIPC',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C write end marker
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IENDID))
      ELSE
C binary write
C write header
         WRITE(UNIT=IUNIT, IOSTAT=IRES) NAMEID(MYID),NREC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'header',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C NSVT,NMOL,NFRAME,NCRTC,NBINC
         WRITE(UNIT=IUNIT, IOSTAT=IRES) NSVT,NMOL,NFRAME,NCRTC,NBINC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'NSVT,NMOL,NFRAME,NCRTC,NBINC',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         WRITE(UNIT=IUNIT, IOSTAT=IRES) XMDIP,DIPC
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'XMDIP,DIPC',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C     
 24   FORMAT(5I10)
 25   FORMAT(6E15.7)
C end wrdms
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RDDMS
C     SUBROUTINE RDDMS(IUNIT,LFORM,NREC,NSVT,NMOL,NFRAME,MAXCRT,NCRTC,
C                      NBINC,MAXBNC,XMDIP,DIPC)
C
C     Reads a dipole moment statistics block from file.
C     The block can be written using L<WRDMS>.
COMMEND
      SUBROUTINE RDDMS(IUNIT,LFORM,NREC,NSVT,NMOL,NFRAME,MAXCRT,NCRTC,
     $     NBINC,MAXBNC,XMDIP,DIPC)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NREC,NSVT,NMOL,NFRAME,NCRTC,MAXCRT,
     $     NBINC,MAXBNC,IRES
      LOGICAL LFORM
      real*8 XMDIP,DIPC
      DIMENSION XMDIP(MAXCRT,MAXBNC,3),DIPC(MAXCRT,NBINC)
C local vars
      INTEGER NDX,NCRTCL,K,I,J
      CHARACTER*60 STR
      CHARACTER PRGSTR*(5)
      DATA PRGSTR /'RDDFS'/
C begin
      IF (LFORM) THEN
C formatted read: read until we get an END or
C NBIN is exceeded
         NCRTCL = 0
         STR = 'DIPMSTAT BLOCK'
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
C read NSVT,NMOL,NFRAME,NCRTC,NBINC
         STR = 'NSVT,NMOL,NFRAME,NCRTC,NBINC'
         CALL CHPINT(NSVT,PRGSTR,STR)
         CALL CHPINT(NMOL,PRGSTR,STR)
         CALL CHPINT(NFRAME,PRGSTR,STR)
         CALL CHPINT(NCRTC,PRGSTR,STR)
         CALL CHPINT(NBINC,PRGSTR,STR)
         CALL CHKMTY(PRGSTR,STR)

C while no marker  AND still atoms to read
         STR = 'XMDIP'
         CALL GIMME(IUNIT,PRGSTR,STR)
         DO I=1,NBINC
            DO J=1,NCRTC
               CALL CHPREL(XMDIP(J,I,1),PRGSTR,STR)
            ENDDO
         ENDDO
         DO I=1,NBINC
            DO J=1,NCRTC
               CALL CHPREL(XMDIP(J,I,2),PRGSTR,STR)
            ENDDO
         ENDDO
         DO I=1,NBINC
            DO J=1,NCRTC
               CALL CHPREL(XMDIP(J,I,3),PRGSTR,STR)
            ENDDO
         ENDDO
C     
         STR = 'DIPC'
         CALL GIMME(IUNIT,PRGSTR,STR)
         DO J=1,NBINC
            DO I=1,NCRTC
               CALL CHPREL(DIPC(I,J),PRGSTR,STR)
            ENDDO
         ENDDO
C
         STR = 'DIPMSTAT BLOCK'
         CALL GNDXID(IUNIT,PRGSTR,STR,NDX)
         IF (NDX .NE. IENDID) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 530,NCRTC,NBINC
            PRINT *,ENDEXP
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ELSE
C unformatted reading
         CALL CHNREC(PRGSTR,NREC,2)
         READ(UNIT = IUNIT, IOSTAT = IRES) NSVT,NMOL,NFRAME,NCRTC,
     $        NBINC
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'NSVT,NMOL,NFRAME,NCRTC,NBINC',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
C
         READ(UNIT = IUNIT, IOSTAT = IRES) XMDIP,DIPC
         IF (IRES .NE. 0)THEN
            PRINT FRFAIL,PRGSTR,'XMDIP,DIPC',IUNIT
            CALL STOPGM('BLOCKIO','FATAL ERROR')
         ENDIF
      ENDIF
C
 530  FORMAT (1X,'attempting to read',I7,' correlation bins and ',
     $     I7,' radial statistics bins')
C end rddms
C     ==--------------------------------------------------------------==
      RETURN
      END

#include "inc_blockio1.h"
