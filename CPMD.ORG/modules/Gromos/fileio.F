Cmb - Revised on 16 Sept. 2005
C     simple file handling for GROMOS
C     Walter Scott ETH Zuerich
C
C------------------------------------------------------------
C     The module fileio.f contains some low level file handling
C     routines. All GROMOS files are opened using the
C     routines OPNFIL or DETOPN and closed using CLSFIL.
C     OPNFIL uses a routine OPNUNT. DETOPN uses a routine
C     DETUNT. These two last routines can be used as well
C     in exceptional circumstances, such as opening files
C     on non GROMOS units.
C     The use of DETOPN and OPNFIL is strongly urged, however.
C 
C     Routines for reading from formatted files with
C     comments (i.e. lines beginning with a # char are considered
C     to be comments) are the following:
C     RAWGET : the low level reading routine.
C              Only use when interested in whether a line is
C              a comment line or not.
C     GETLN  : returns the next line in the file which is not
C              a comment line
C     GIMME  : same as GETLN, but also calls FLAGLN if an
C              error occurs in the reading.
C
C     It is recommended to use GIMME whenever possible.
C
C     For printing out errors, it is nice to give the number
C     the error the line is on. The line number is initialized in
C     the routine INITLN, updated in the routine RAWGET, and can be
C     accessed by the routine LNNUM. However, indirect access
C     of the linenum through the routines FLAGLN or ERRLN
C     is strongly urged.
C     INITLN : resets the line number counter. Call this routine
C              after opening the file.
C     FLAGLN : writes out the subroutine name and line number
C              in the standard GROMOS format for errors.
C     ERRLN :  calls FLAGLN, writes a user specified error
C              description, and calls STOP.
C
C
C     Typical use of these routines would be:
C     (the following code is extracted from the routine RDTOPO)
C
C     open the file and initialize the line counter
C
C      CALL OPNFIL('TOPO','FORMATTED','OLD',IUNIT)
C      IF (IUNIT .LT. 0) THEN
C         PRINT *,'failed to open file!'
C         STOP
C      ENDIF
C
C     CALL INITLN(IUNIT)
C     PRGSTR is assigned the name of the subroutine doing the reading.
C     This is used for writing out errors and stays the same
C     throughout the subroutine
C     PRGSTR = 'RDTOPO'
C
C     In the following code, we have the following declaration:
C     INCLUDE 'fileio.h'
C 
C     Now start reading from the file....
C     We want to read an integer NRATT and a REAL*8 RBLA that
C     appear on a line together.
C     a) Assign the name of the variables we want to read
C     this is used for writing out errors.
C     STR = 'NRATT, RBLA'
C
C     b) Get the next data line from file.
C     the line is read into a character array FIOLIN
C     which is normally not accessed directly.
C     An error is written out if this read fails.
C     CALL GIMME(IUNIT,PRGSTR,STR)
C
C     c) We now have the data line, chop of the values
C     we want to read using CHPINT for integers and
C     CHPREL for reals. Both routines call ERRLN if
C     something goes wrong on reading.
C
C     CALL CHPINT(NRATT,PRGSTR,STR)
C     CALL CHPREL(RBLA,PRGSTR,STR)
C
C     d) Here we can check the values for correctness.
C     Use FLAGLN for writing the line number to screen.
C     E.g., for NRATT that could look like this:
C     IF (NRATT .LT. 0 .OR. NRATT .GT. MAXATT) THEN
C        CALL FLAGLN(PRGSTR)
C        PRINT *,'NRATT out of range!'
C        PRINT *,'read ',NRATT,'MAXATT=',MAXATT
C        STOP
C     ENDIF
C
C     repeat steps a)...d) for all other variable(s)
C
C     now close the file
C     CALL CLSFIL(IUNIT)
C
C------------------------------------------------------------


COMMSUBR  DETOPN
C     SUBROUTINE DETOPN(FSTR,IUNIT,LFORM)
C
C     DETOPN attempts to open an existing file and to
C     determine whether it is formatted or unformatted.
C     DETOPN first opens a file using L<OPNFIL>.
C     If OPNFIL returns a negative value in IUNIT
C     meaning that the open failed (e.g. no file there)
C     then DETOPN just returns that negative value 
C     in IUNIT.
C
C     If OPNFIL successfully opens a file, then
C     DETOPN attempts to determine its format.
C     If it fails to do this, it calls a STOP.
C
C     DETOPN returns IUNIT and LFORM on success.
C
C     Use this routine for READING configurations only.
C     For writing configurations, use OPNFIL specifying the
C     type of file (i.e. formatted, unformatted) required.
C
C     FSTR is a string determining what type of GROMOS
C        file to open -- such as 'XVI' for coordinate
C        files. See the routine L<OPNFIL> in fileio.f for
C        all possible values for FSTR.
C
C     IUNIT delivered with the unit opened.
C
C     LFORM delivered with 
C      .TRUE. if the file is formatted
C      .FALSE.if the file is unformatted ( binary)
COMMEND
      SUBROUTINE DETOPN(FSTR,IUNIT,LFORM)
C args
      CHARACTER FSTR*(*)
      INTEGER IUNIT
      LOGICAL LFORM
C includes
       INCLUDE 'units.h'
C local vars
      INTEGER NDX
      CHARACTER PRGSTR *(6)
      DATA PRGSTR /'DETOPN'/
C begin
C     determine the filetype
C     while not found and NDX < table size
      NDX = 1
 10   IF (NDX .LE. MAXNTS .AND. UNAME(NDX) .NE. FSTR) THEN
         NDX = NDX + 1
         GOTO 10
      ENDIF

      IF (NDX .GT. MAXNTS) THEN
         PRINT *,PRGSTR,': Error !'
         PRINT *,'Unknown FTNAME ''',FSTR,''''
         CALL STOPGM('DETOPN','FATAL ERROR')
      ELSE
         IUNIT = IUNUM(NDX)
      ENDIF

      CALL DETUNT(PRGSTR,IUNIT,LFORM)
C end detopn
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR DETUNT
C     SUBROUTINE DETUNT(PRGSTR,IUNIT,LFORM)
C
C     IN: PRGSTR
C     INOUT: IUNIT
C     OUT: LFORM
C
C     A low level routine to determine whether a unit IUNIT
C     is formatted or not. The file must start with a valid
C     GROMOS block header to be recognized.
C     This routine can be used when opening an existing file
C     for reading. 
C     PRGSTR: string of the calling routine. Used in printing
C     error messages.
C     IUNIT: the unit number to be checked.
C     On success, IUNIT is not changed.
C     On failure, IUNIT is returned with  a negative value
C
C     LFORM: .TRUE. if the file is formatted.
C     .FALSE. otherwise. 
C     The value of LFORM is undefined if IUNIT is negative.
C
C     Strategy used in this routine:
C     assume the file is formatted and attempt
C     to read the first block header.
C     If that fails, try reading the block header
C     unformatted.
C     If that fails, we give up -- probably the file
C     is unformatted but written by a machine with different
C     binary representation.
COMMEND
      SUBROUTINE DETUNT(PRGSTR,IUNIT,LFORM)
C args
      CHARACTER PRGSTR*(*)
      INTEGER IUNIT
      LOGICAL LFORM
C include
      INCLUDE 'dataid.h'
C local vars
      INTEGER IDBLK,NREC,IRES
C begin
C try formatted first
      LFORM = .TRUE.
      CALL OPNUNT(PRGSTR,'FORMATTED','OLD',IUNIT)
      IF (IUNIT .LT. 0) THEN
C the file doesnt exist: just return
C so we can print out something meaningful.
         RETURN 
      ENDIF
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .GT. IDERR) THEN
         REWIND(UNIT = IUNIT,IOSTAT = IRES)
         IF (IRES .NE. 0)THEN
            PRINT 500,IUNIT
            CALL STOPGM('FILEIO','FATAL ERROR')
         ENDIF
         CALL INITLN(IUNIT)
      ELSE
C the formatted reading failed:
C close and reopen unformatted
         CALL CLSFIL(IUNIT)
         LFORM = .FALSE.
         CALL OPNUNT(PRGSTR,'UNFORMATTED','OLD',IUNIT)
         CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
         IF (IDBLK .LE. IDERR) THEN
C something I can t read
            PRINT *,PRGSTR,': unrecognized file contents!'
            PRINT *,'UNIT = ',IUNIT
            PRINT *,
     $           'The file must start with a GROMOS Block',
     $           ' Type (E.g. TITLE) to be recognized'
            PRINT *
            IUNIT = -1
            RETURN
         ELSE
C this file is unformatted:
C rewind before returning
            REWIND(UNIT = IUNIT,IOSTAT = IRES)
            IF (IRES .NE. 0)THEN
               PRINT 500,IUNIT
               CALL STOPGM('FILEIO','FATAL ERROR')
            ENDIF
         ENDIF
      ENDIF
 500  FORMAT ('DETUNT: rewind failed on UNIT =',I5)
C end detunt
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR OPNUNT
C     SUBROUTINE OPNUNT(PRGSTR,FORMST, STATST, IUNIT)
C     IN: PRGSTR,FORMST, STATST
C     INOUT: IUNIT
C
C     The low level file opening routine.
C     It is attempted to open the unit IUNIT with the
C     format and status strings given. 
C     On a successful open, IUNIT is returned unchanged.
C     On a failed open, IUNIT is return with a negative value.
COMMEND
      SUBROUTINE OPNUNT(PRGSTR,FORMST, STATST, IUNIT)
C includes
       INCLUDE 'units.h'
C params
       CHARACTER PRGSTR*(*),FORMST*(*), STATST*(*)
       INTEGER IUNIT
C local vars
      INTEGER IRES
C      LOGICAL LEXIST
      CHARACTER FILEST*(12)
C begin

C     standard input requires special handling:
C     fortran unit 5 is always open and so we have to read from it,
C     even if there is a file fort.5 on disk.
C     if the fort.5 exists on disk, try to open it and read
C     standard input from there.
C     Otherwise, assume the user has piped something to the process.
      IF (IUNIT .EQ. ISTDIN) THEN
         RETURN
C        INQUIRE(FILE='fort.5',IOSTAT=IRES,EXIST=LEXIST)
C        IF (IRES .NE. 0) THEN
C            PRINT *,PRGSTR,': INQUIRE failed'
C            PRINT 500,IUNIT,IRES
C            CALL STOPGM('FILEIO','FATAL ERROR')
C         ENDIF
C         IF (.NOT. LEXIST) THEN
C user must have piped something and the file is already open
C            RETURN
C         ENDIF
      ENDIF

C standard F77 and F90:
C if STATUS = 'SCRATCH' there MAY NOT be a FILE=
C if STATUS = 'OLD', 'NEW', or 'REPLACE' there MUST be a FILE=
C this is usually not enforced (Sun, Microsoft, SGI), but is on
C Crays.
      IF (STATST .EQ. 'SCRATCH') THEN
         OPEN(IUNIT, IOSTAT=IRES, STATUS=STATST, FORM=FORMST)
      ELSE
C make my filename string
         IF (IUNIT .LT. 10) THEN
            WRITE(FILEST,IOSTAT=IRES,FMT='(''fort.'',I1)')IUNIT
         ELSE IF (IUNIT .LT. 100) THEN
            WRITE(FILEST,IOSTAT=IRES,FMT='(''fort.'',I2)')IUNIT
         ELSE IF (IUNIT .LT. 1000) THEN
            WRITE(FILEST,IOSTAT=IRES,FMT='(''fort.'',I3)')IUNIT
         ELSE IF (IUNIT .LT. 10000) THEN
            WRITE(FILEST,IOSTAT=IRES,FMT='(''fort.'',I4)')IUNIT
         ENDIF
         IF (IRES .NE. 0) THEN
            PRINT *,PRGSTR,': format error!'
            PRINT 500,IUNIT,IRES
            CALL STOPGM('FILEIO','FATAL ERROR')
         ENDIF
C try to open
         OPEN(IUNIT,IOSTAT=IRES, STATUS=STATST, FORM=FORMST,FILE=FILEST)
      ENDIF
      IF (IRES .NE. 0) THEN
         PRINT *,PRGSTR,': Error on call to OPEN !'
         PRINT *,'FORM: ',FORMST,' STATUS: ',STATST
         PRINT 500,IUNIT,IRES
C-----
C here try to diagnose what went wrong without using machine dependent
C value of IRES
         IF (STATST .EQ. 'NEW') THEN
C most probably the open failed because the file
C exists already
            PRINT *,'Does the file exist already ?'
         ELSEIF (STATST .EQ. 'OLD') THEN
C most probably the open failed because the file
C doesnt exist
            PRINT *,'Does the file exist at all ?'
         ENDIF
C add any other 'diagnostics' between  here....

C ... and here
         IUNIT = -1
      ENDIF

 500  FORMAT (1X,'OPEN (UNIT =',I3,') failed with IOSTAT = ',I6)
C end opnunt
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR OPNFIL
C     SUBROUTINE OPNFIL(FTNAME, FORMST, STATST, IUNIT)
C      IN : FTNAME, FORMST, STATST
C      OUT: IUNIT
C
C     OPNFIL opens a file whose TYPE is in FTNAME with the
C     attributes FORMST and STATST.
C
C     FORMST  determines whether the file should be
C             formatted or unformatted.
C             Valid vals for FORMST are: 'FORMATTED', 'UNFORMATTED'
C     STATST  determines the status of the file to be opened.
C             Valid vals for STATST are:
C             'NEW', 'OLD', 'SCRATCH' or 'UNKNOWN'.
C     IUNIT   delivered with the unit number opened if successful.
C             The unit is returned in IUNIT >0 0 on success.
C             On failure, the IUNIT is < 0 and
C             an error is printed.
C
C     example of use:
C     CALL OPNFIL('XVI','FORMATTED','OLD', MYFILE)
C        IF (MYFILE .LT. 0) THEN
C        PRINT *,'failed to open input coordinate file'
C        CALL STOPGM('FILEIO','FATAL ERROR')
C     ENDIF
C     here we can use MYFILE for IO
COMMEND
       SUBROUTINE OPNFIL(FTNAME, FORMST, STATST, IUNIT)
C includes
       INCLUDE 'units.h'
C params
       CHARACTER  FTNAME* (*),FORMST*(*), STATST*(*)
       INTEGER IUNIT
C local vars
      INTEGER NDX
      CHARACTER PRGSTR*(6)
C data
      DATA PRGSTR/'OPNFIL'/
C begin
C     determine the filetype
C     while not found and NDX < table size
      NDX = 1
 10   IF (NDX .LE. MAXNTS .AND. UNAME(NDX) .NE. FTNAME) THEN
         NDX = NDX + 1
         GOTO 10
      ENDIF

      IF (NDX .GT. MAXNTS) THEN
         PRINT *,PRGSTR,': Error !'
         PRINT *,'Unknown FTNAME ''',FTNAME,''''
         CALL STOPGM('FILEIO','FATAL ERROR')
      ELSE
         IUNIT = IUNUM(NDX)
      ENDIF

      CALL OPNUNT(PRGSTR,FORMST,STATST, IUNIT)
C end OPNFIL
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR CLSFIL
C     SUBROUTINE CLSFIL(IU)
C
C     CLSFIL closes a file previously opened with L<OPNFIL>.
COMMEND
      SUBROUTINE CLSFIL(IU)
      INTEGER IU
C begin
      CLOSE(IU)
C end clsfil
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR GOTOLN
C     SUBROUTINE GOTOLN(IUNIT,LNNUM)
C
C     A slow but portable way to set a formatted file
C     to a certain position.
C     The line in a file are numbered from 1 onwards.
C     Comment lines are also considered as lines.
C     The file is positioned such that a GETLN or GIMME
C     will read in the line number LNNUM.
COMMEND
      SUBROUTINE GOTOLN(IUNIT,LNNUM)
      INTEGER IUNIT,LNNUM
C local vars
      INTEGER IRES,I
C common
      INTEGER ILNUM,ISVUNT
      COMMON /SCOLIN/ILNUM,ISVUNT
C begin
      REWIND(IUNIT,IOSTAT=IRES)
      CALL INITLN(IUNIT)
      IF (IRES .NE. 0) THEN
         PRINT *,'GOTOLN: failed to rewind unit',IUNIT
         CALL STOPGM('FILEIO','FATAL ERROR')
      ENDIF
      DO I = 1,LNNUM
         CALL RAWGET(IUNIT,IRES)
         IF (IRES .NE. 0) THEN
            PRINT *,'GOTOLN: skipping failed for unit',IUNIT
            CALL STOPGM('FILEIO','FATAL ERROR')
         ENDIF
      ENDDO
C     end gotoln
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR INITLN
C     SUBROUTINE INITLN(IUNIT)
C
C     INITLN intializes the internal line counter for a
C     the unit IUNIT accessed by SUBR. L<LNNUM> used for
C     printing out error messages with the line of the error.
C
C     This routine should be called immediately after the file IUNIT
C     has been opened.
COMMEND
      SUBROUTINE INITLN(IUNIT)
      INTEGER IUNIT
C common
      INTEGER ILNUM,ISVUNT
      COMMON /SCOLIN/ILNUM,ISVUNT
C begin
      ILNUM = 0
      ISVUNT = IUNIT
C end
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR LNNUM
C     INTEGER FUNCTION LNNUM()
C
C     LNNUM returns the current value of the line number
C     this counter is initialized by the SUBR. L<INITLN>
C     and is incremented every time a line is read in by
C     SUBR. L<RAWGET>.
COMMEND

      INTEGER FUNCTION LNNUM()
C common
      INTEGER ILNUM,ISVUNT
      COMMON /SCOLIN/ILNUM,ISVUNT
C begin
      LNNUM = ILNUM
C end LNNUM
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR LNUNIT
C     INTEGER FUNCTION LNUNIT()
C
C     LNUNIT returns the value of the unit last passed
C     to INITLN
COMMEND
      INTEGER FUNCTION LNUNIT()
C common
      INTEGER ILNUM,ISVUNT
      COMMON /SCOLIN/ILNUM,ISVUNT
C begin
      LNUNIT = ISVUNT
C end LNUNIT
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RAWGET
C     SUBROUTINE RAWGET(IUNIT,IRES)
C
C     RAWGET is the lowest level reading routine
C     for GROMOS formatted files. It read a line of
C     90 chars regardless of whether the line contains
C     a comment (i.e. the line begins with a # char) or not,
C     and increments the line counter that is initialized by
C     SUBR. L<INITLN> and can be accessed by calling
C     SUBR. L<LNNUM>.
C
C     The line read is stored in FIOLIN in fileio.h
COMMEND 
      SUBROUTINE RAWGET(IUNIT,IRES)
C includes
      INCLUDE 'fileio.h'
C args
      INTEGER IUNIT,IRES
C common
      INTEGER ILNUM,ISVUNT
      COMMON /SCOLIN/ILNUM,ISVUNT
C begin
      READ(IUNIT, IOSTAT = IRES,FMT = '(A90)') FIOLIN
      ILNUM = ILNUM + 1
      NDXRI = 0
      NDXLE = 0
C end RAWGET
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR RESTLN
C     Reset the reading of a line, so that
C     we can read from the beginning of the line
C     again using the CHPXXX routines.
C     OMMEND
      SUBROUTINE RESTLN
      INCLUDE 'fileio.h'
C begin
      NDXRI = 0
      NDXLE = 0
C     ==--------------------------------------------------------------==
      RETURN
      END
      
COMMSUBR GETLN
C     SUBROUTINE GETLN(IUNIT,LOK)
C
C     GETLN returns the next line on the file
C     IUNIT which is not a comment, i.e. does not start
C     with a # character.
C
C     The line read in is delivered in FIOLIN which is in
C     a common block defined in 'fileio.h'
C     The logical variable LOK = .TRUE. iff the
C     reading was successful.
COMMEND
      SUBROUTINE GETLN(IUNIT,LOK)
C includes
      INCLUDE 'fileio.h'
C arguments
      INTEGER IUNIT
      LOGICAL LOK
C common
      INTEGER ILNUM,ISVUNT
      COMMON /SCOLIN/ILNUM,ISVUNT
C local vars
      LOGICAL LCOMM
      INTEGER IRES
      CHARACTER*(MAXDLN) BLANK
C data
      DATA BLANK / ' ' /
C begin
 10   CONTINUE
C     repeat until error or not blank line or not comment
C---------------
C     here we inline RAWGET instead of calling it for speed.
C     CALL RAWGET(IUNIT,IRES)
C---------------
      READ(IUNIT, IOSTAT = IRES,FMT = '(A90)') FIOLIN
      ILNUM = ILNUM + 1
      NDXRI = 0
      NDXLE = 0
C-----------------END OF RAWGET INLINE
         LCOMM = (FIOLIN(1:1) .EQ. '#') .OR. (FIOLIN .EQ. BLANK)
      IF ((IRES .EQ. 0) .AND. LCOMM) GOTO 10

      LOK = (IRES .EQ. 0)
C end getline
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR GIMME
C     SUBROUTINE GIMME(IUNIT,PROG,STR)
C
C     GIMME reads a line from the open file IUNIT
C     using L<GETLN>.
C     If an error occurs on reading, an error is written
C     out using L<FLAGLN>, the faulty line read in is
C     printed out, a line containing "Error reading STR"
C     is printed and a STOP is envoked.
C     STR my contain variable names that are expected by the
C     program.
COMMEND
      SUBROUTINE GIMME(IUNIT,PROG,STR)
C includes
      INCLUDE 'fileio.h'
C arguments
      INTEGER IUNIT
      CHARACTER *(*) PROG
      CHARACTER *(*) STR
C local vars
      LOGICAL LOK
C begin
      CALL GETLN(IUNIT,LOK)
      IF (.NOT. LOK) THEN
         CALL FLAGLN(PROG)
         PRINT '(A1,A75,A1)',''',FIOLIN,'''
         PRINT *,'Error reading ',STR
         CALL STOPGM('FILEIO','FATAL ERROR')
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR GTNREL
C      SUBROUTINE GTNREL(IUNIT,NWANT,X,NUMRED)
C
C     GTNREL attempts to read in NWANT reals from
C     the formatted file IUNIT into X (NWANT).
C     The actual number of reals read in is returned
C     in NUMRED. If the routine succeeded, NUMRED will be
C     delivered with NWANT. If not then NUMRED is returned with
C     a value .NE. NWANT.
C
C     Any number of comment lines AT THE BEGINNING of the
C     block of NWANT reals are skipped without an error, and
C     the read of reals is reattempted.
C
C     If the read fails, the file is positioned with BACKSPACE
C     to where it was before GNTREL was called.
C
C     NOTE: THIS ROUTINE ATTEMPTS TO UPDATE THE VARIABLE INDICATING
C     THE CURRENT LINE NUMBER (ILNUM). HOWEVER, AS THE REALS ARE READ
C     IN USING A LIST DIRECTED READ (FMT=*), MORE LINES CAN BE READ IN
C     THAN INDICATED. THIS HAPPENS, E.G. IN THE CASE OF
C     BLANK LINES. ALSO, WHEN A VARYING NUMBER OF REALS ARE ON ONE
C     LINE, ILNUM WILL BE WRONG.
C     ILNUM IS INCREMENTED, BY NWANT/3 ASSUMING THAT THERE ARE
C     THREE REALS ON ONE LINE.
C
COMMEND
      SUBROUTINE GTNREL(IUNIT,NWANT,X,NUMRED)
C includes
      INCLUDE 'fileio.h'
C args
      INTEGER IUNIT,NWANT,NUMRED
      real*8 X(NWANT)
C local vars
      INTEGER IRES
      LOGICAL LSTOP
Ccommon
      INTEGER ILNUM,ISVUNT
      COMMON /SCOLIN/ILNUM,ISVUNT
C begin
      NUMRED = -1
      LSTOP = .FALSE.
 10   IF (.NOT. LSTOP) THEN
         READ(UNIT=IUNIT,IOSTAT=IRES,FMT=*)X
         ILNUM = ILNUM + NWANT/3
         NDXRI = 0
         NDXLE = 0
         IF (IRES .NE. 0)THEN
C the read failed, try to do a backspace and read a string to see
C if its a comment line. If not the read fails.
            BACKSPACE(UNIT=IUNIT,IOSTAT=IRES)
            IF (IRES .NE. 0) THEN
               PRINT *,'GNTREL: BACKSPACE FAILED (1)!!!',IRES
               CALL STOPGM('FILEIO','FATAL ERROR')
            ENDIF
CC            PRINT *,'BACK 1'
            ILNUM = ILNUM - NWANT/3
C get a string. This is an inline version of RAWGET
            READ(IUNIT, IOSTAT = IRES,FMT = '(A90)') FIOLIN
            ILNUM = ILNUM + 1
C stop the loop if we have a read error or no comment line
            LSTOP = ((IRES .NE. 0) .OR. (FIOLIN(1:1) .NE. '#'))
            IF (LSTOP) THEN
               BACKSPACE(UNIT=IUNIT,IOSTAT=IRES)
               IF (IRES .NE. 0) THEN
                  PRINT *,'GNTREL: BACKSPACE FAILED (2)!!!',IRES
                  CALL STOPGM('FILEIO','FATAL ERROR')
               ENDIF
CC               PRINT *,'BACK 2',LSTOP
CC               PRINT *,FIOLIN
               ILNUM = ILNUM - 1
            ENDIF
CC            PRINT *,'COMM',FIOLIN
         ELSE
            NUMRED = NWANT
            LSTOP = .TRUE.
         ENDIF
         GOTO 10
      ENDIF
C end gtnrel
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR FLAGLN
C     SUBROUTINE FLAGLN(PROG)
C
C     FLAGLN writes out an error in a standard format.
C     All programs reading from formatted GROMOS files
C     should use this subroutine to flag errors.
COMMEND
      SUBROUTINE FLAGLN(PROG)
      CHARACTER *(*) PROG
C includes
      INCLUDE 'fileio.h'
C define external
      INTEGER LNNUM,LNUNIT
      EXTERNAL LNNUM,LNUNIT
C local vars
      INTEGER NLINE,NUNIT
C begin
      NLINE = LNNUM()
      NUNIT = LNUNIT()

      PRINT 50,PROG,NUNIT,NLINE

C     As an option, we can also write a line that emacs
C     can interprete for loading the file and moving
C     to the line with the error on it.
C
C      IF (NUNIT .LT. 10) THEN
C         PRINT 100,NUNIT,NLINE,PROG
C      ELSE
C         PRINT 101,NUNIT,NLINE,PROG
C      ENDIF

C now print what we have read on the input
      PRINT '(A1,A75,A1)','"',FIOLIN,'"'

 50   FORMAT (1X,A6,': Unit:',I4,' Line:',I10)
C     may not use a leading 1X here, or emacs wont parse
C     error list properly
C 100  FORMAT ('fort.',I1,':',I10,': Error in subroutine: ',A)
C 101  FORMAT ('fort.',I2,':',I10,': Error in subroutine: ',A)
C end flagln
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR ERRLN
C     SUBROUTINE ERRLN(PROG,STR)
C
C     ERRLN writes out an error using L<FLAGLN>,
C     then prints out STR which typically contains a
C     further explanation and then calls a STOP.
COMMEND
      SUBROUTINE ERRLN(PROG,STR)
C args
      CHARACTER *(*) PROG
      CHARACTER *(*) STR
C begin errln
      CALL FLAGLN(PROG)
      PRINT *,'Error formatting ',STR
      CALL STOPGM('FILEIO','FATAL ERROR')
C end errln
C     ==--------------------------------------------------------------==
      RETURN
      END

C return whether a line is empty or not
      LOGICAL FUNCTION LISMTY()
      INCLUDE 'fileio.h'
C local params
      INTEGER ITAB
      PARAMETER (ITAB = 9)
C local vars
      INTEGER LELOC
C begin
      LELOC = NDXRI+1
 10   IF (LELOC .LE. LEN(FIOLIN) ) THEN
         IF (FIOLIN(LELOC:LELOC) .EQ. ' '.OR.
     $        FIOLIN(LELOC:LELOC) .EQ. CHAR(ITAB)) THEN
            LELOC = LELOC + 1
            GOTO 10
         ENDIF
      ENDIF
      LISMTY =  (LELOC .GT. MAXDLN)
C     PRINT *,'LISMTY',NDXLE,LELOC,LISMTY
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR CHKMTY
C     SUBROUTINE CHKMTY(PRGSTR,STR)
C
C     If there are any non space characters left on the line
C     in the buffer, FLAGLN is called, a message is printed
C     and STOP is envoked.
C
C     This routine can be used to make sure all the
C     input on a line has been processed.
COMMEND
      SUBROUTINE CHKMTY(PRGSTR,VARSTR)
C args
      CHARACTER PRGSTR*(*)
      CHARACTER VARSTR*(*)
C externals
      LOGICAL LISMTY
      EXTERNAL LISMTY
C begin
      IF (.NOT. LISMTY()) THEN
         CALL FLAGLN(PRGSTR)
         PRINT *,'More input on line than expected!'
         PRINT *,VARSTR,'expected'
         CALL STOPGM('FILEIO','FATAL ERROR')
      ENDIF
C end chkmty
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR CHPSTR
C     SUBROUTINE CHPSTR(RETSTR,PRGSTR,STR)
C
C     CHPSTR chops a string from the current line and returns
C     it in RETSTR.
C     The string may not have any spaces in it.
C     Any preceding spaces are skipped.
C
C     This routine calls L<GETLN> if the last string read in is empty.
C     If GETLN returns LOK = .FALSE. indicating that a read error
C     (due to end of file) occurred, an error messages is printed and
C     STOP is envoked.
COMMEND
      SUBROUTINE CHPSTR(RETSTR,PRGSTR,STR)
      INCLUDE 'fileio.h'
C args
      CHARACTER RETSTR*(*)
      CHARACTER PRGSTR*(*)
      CHARACTER STR*(*)
C common
      INTEGER ILNUM,ISVUNT
      COMMON /SCOLIN/ILNUM,ISVUNT
C local params
C     the ASCII code for a TAB character
      INTEGER ITAB
      PARAMETER (ITAB = 9)
C externals
      LOGICAL LISMTY
      EXTERNAL LISMTY
C local vars
      LOGICAL LOK
      CHARACTER CHRTAB
C begin
      CHRTAB = CHAR(ITAB)
      NDXLE = NDXRI+1
      IF (LISMTY()) THEN
         CALL GETLN(ISVUNT,LOK)
         IF (.NOT. LOK) THEN
            PRINT *,PRGSTR,': end of file reached!'
            CALL STOPGM('FILEIO','FATAL ERROR')
         ENDIF
         NDXLE = NDXRI+1
      ENDIF

C advance left until no blanks or tabs
 10   IF (NDXLE .LE. MAXDLN .AND. 
     $     (FIOLIN(NDXLE:NDXLE) .EQ. ' ' .OR.
     $     FIOLIN(NDXLE:NDXLE) .EQ. CHRTAB)) THEN
C         PRINT *,NDXLE,ICHAR(FIOLIN(NDXLE:NDXLE))
         NDXLE = NDXLE + 1
         GOTO 10
      ENDIF
      IF (NDXLE .GT. MAXDLN) THEN
         CALL ERRLN(PRGSTR,STR)
      ENDIF

C advance right until blanks or tabs
      NDXRI = NDXLE
 20   IF (NDXRI .LE. MAXDLN .AND.
     $     (FIOLIN(NDXRI:NDXRI) .NE. ' ' .AND.
     $     FIOLIN(NDXRI:NDXRI) .NE. CHRTAB )) THEN
         NDXRI = NDXRI + 1
         GOTO 20
      ENDIF
      IF (NDXRI .GT. MAXDLN) THEN
         CALL ERRLN(PRGSTR,STR)
      ENDIF

      RETSTR = FIOLIN(NDXLE:NDXRI-1)
C      PRINT *,'CHPSTR',RETSTR
C end chpstr
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR SKPCHR
C     SUBROUTINE SKPCHR(ICHRS)
C
C     skip ICHRS characters from the line read in.
C     This routine can be used to discard the first
C     ICHRS characters on a line.
COMMEND
      SUBROUTINE SKPCHR(ICHRS)
      INCLUDE 'fileio.h'
C args
      INTEGER ICHRS
C begin
      IF (ICHRS .GE. MAXDLN) THEN
         PRINT *,'Error in SKPCHR!'
         PRINT *,'ICHRS too large!',ICHRS
         CALL STOPGM('FILEIO','FATAL ERROR')
      ENDIF
      NDXRI = ICHRS
C end skpchr
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR CHPINT
C     SUBROUTINE CHPINT(IRET,PRGSTR,STR)
C
C     CHPINT reads a integer value from the input stream.
C     The value read is returned in IRET.
C     If an error occurs, a message is printed to output
C     and STOP is envoked.
COMMEND
      SUBROUTINE CHPINT(IRET,PRGSTR,STR)
C includes
      INCLUDE 'fileio.h'
C args
      INTEGER IRET
      CHARACTER PRGSTR*(*)
      CHARACTER STR*(*)
C local vars
      INTEGER IRES,LENSTR
      CHARACTER FMTSTR*(8)
      CHARACTER STRVAL*(20)
C begin
      CALL CHPSTR(STRVAL,PRGSTR,STR)
      LENSTR = NDXRI-NDXLE

      IF (LENSTR .LE. FMIMAX) THEN
C try to read internally
         READ (STRVAL, IOSTAT = IRES, FMT = FMITAB(LENSTR))IRET
         IF (IRES .NE. 0) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 600,STR
            PRINT 610,STRVAL(1:LENSTR)
            PRINT 620,IRES,FMITAB(LENSTR)
            CALL STOPGM('FILEIO','FATAL ERROR')
         ENDIF
      ELSE
C create my format string
         WRITE(FMTSTR,IOSTAT = IRES,FMT = 500)LENSTR
         IF (IRES .NE. 0) THEN
            PRINT *,'CHPINT: ERROR IN FORMAT'
            CALL STOPGM('FILEIO','FATAL ERROR')
         ENDIF

C try to read internally
         READ (STRVAL, IOSTAT = IRES, FMT = FMTSTR)IRET
         IF (IRES .NE. 0) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 600,STR
            PRINT 610,STRVAL(1:LENSTR)
            PRINT 620,IRES,FMTSTR
            CALL STOPGM('FILEIO','FATAL ERROR')
         ENDIF
      ENDIF

 500  FORMAT ('(I',I3.3,')')
 600  FORMAT ('failed to read an INTEGER for ',A)
 610  FORMAT ('got ''',A,'''')
 620  FORMAT ('IRES =',I3,', FMTSTR = ',A)
C end CHPINT
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR CHPREL
C     SUBROUTINE CHPREL(RETVAL,PRGSTR,STR)
C
C     CHPREL reads a REAL*8 value from the input stream.
C     The value read is returned in RETVAL.
COMMEND

      SUBROUTINE CHPREL(RETVAL,PRGSTR,STR)
C includes
      INCLUDE 'fileio.h'
C local params
      INTEGER ITAB
      PARAMETER (ITAB = 9)
C args
      real*8 RETVAL
      CHARACTER PRGSTR*(*)
      CHARACTER STR*(*)
C local vars
      LOGICAL LDOT
      INTEGER I,IRES,LENSTR
      CHARACTER FMTSTR*(8)
      CHARACTER STRVAL*(40)
C begin
      CALL CHPSTR(STRVAL,PRGSTR,STR)
      LENSTR = NDXRI-NDXLE

C see if we have a dot ('.') in the string we have read in.
C if not, append one
      LDOT = .FALSE.
      I = 1
 10   IF (I .LE. LENSTR) THEN
         IF (.NOT. LDOT) THEN
            LDOT = (STRVAL(I:I) .EQ. '.')
            I = I + 1
            GOTO 10
         ENDIF
      ENDIF

      IF (.NOT. LDOT) THEN
         LENSTR = LENSTR + 1
         STRVAL(LENSTR:LENSTR) = '.'
      ENDIF

C     if we have a prestored format string, we use that;
C     otherwise we make our own format string
C     Note that list driven reading from internal files
C     is non standard.
      IF (LENSTR .LE. FMRMAX) THEN
C     try to read internally
         READ (STRVAL(1:LENSTR), IOSTAT=IRES, FMT=FMRTAB(LENSTR))RETVAL
         IF (IRES .NE. 0) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 600,STR
            PRINT 610,STRVAL(1:LENSTR)
            PRINT 620,IRES,FMRTAB(LENSTR)
            CALL STOPGM('FILEIO','FATAL ERROR')
         ENDIF
      ELSE
C create my format string
         WRITE(FMTSTR,IOSTAT = IRES,FMT = 500)LENSTR
         IF (IRES .NE. 0) THEN
            PRINT *,'CHPREL: ERROR IN FORMAT'
            CALL STOPGM('FILEIO','FATAL ERROR')
         ENDIF

C     try to read internally
         READ (STRVAL(1:LENSTR), IOSTAT=IRES, FMT=FMTSTR)RETVAL
         IF (IRES .NE. 0) THEN
            CALL FLAGLN(PRGSTR)
            PRINT 600,STR
            PRINT 610,STRVAL(1:LENSTR)
            PRINT 620,IRES,FMTSTR
            CALL STOPGM('FILEIO','FATAL ERROR')
         ENDIF
      ENDIF

 500  FORMAT ('(F',I3.3,'.1)')
 600  FORMAT ('failed to read a real*8 for ',A)
 610  FORMAT ('got ''',A,'''')
 620  FORMAT ('IRES =',I3,', FMTSTR = ',A)
C end CHPREL
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR CHPNRE
C     SUBROUTINE CHPNRE(NVAL,RETARR,PRGSTR,STR)
C
C     CHPREL reads NVAL REAL*8 values from the input stream.
C     The values read are returned in RETARR(NVAL).
COMMEND
      SUBROUTINE CHPNRE(NVAL,RETARR,PRGSTR,STR)
C includes
      INCLUDE 'fileio.h'
C local params
      INTEGER ITAB
      PARAMETER (ITAB = 9)
C args
      INTEGER NVAL
      real*8 RETARR(NVAL)
      CHARACTER PRGSTR*(*)
      CHARACTER STR*(*)
C local vars
      LOGICAL LDOT
      INTEGER IRES,LENSTR,I,II
      real*8 RETVAL
      CHARACTER FMTSTR*(8)
      CHARACTER STRVAL*(40)
C begin
C this could be optimized for speed.....
      DO 10 I=1,NVAL
         CALL CHPSTR(STRVAL,PRGSTR,STR)

         LENSTR = NDXRI-NDXLE
C see if we have a dot ('.') in the string we have read in.
C if not, append one
         LDOT = .FALSE.
         II = 1
 20      IF (II .LE. LENSTR) THEN
            IF (.NOT. LDOT) THEN
               LDOT = (STRVAL(II:II) .EQ. '.')
               II = II + 1
               GOTO 20
            ENDIF
         ENDIF
         
         IF (.NOT. LDOT) THEN
            LENSTR = LENSTR + 1
            STRVAL(LENSTR:LENSTR) = '.'
         ENDIF

         IF (LENSTR .LE. FMRMAX) THEN
C     try to read internally
            READ (STRVAL(1:LENSTR), IOSTAT=IRES,
     $           FMT=FMRTAB(LENSTR))RETVAL
            IF (IRES .NE. 0) THEN
               CALL FLAGLN(PRGSTR)
               PRINT 600,STR
               PRINT 610,STRVAL(1:LENSTR)
               PRINT 620,IRES,FMRTAB(LENSTR)
               CALL STOPGM('FILEIO','FATAL ERROR')
            ENDIF
         ELSE
C create my format string
            WRITE(FMTSTR,IOSTAT = IRES,FMT = 500)LENSTR
            IF (IRES .NE. 0) THEN
               PRINT *,'CHPREL: ERROR IN FORMAT'
               CALL STOPGM('FILEIO','FATAL ERROR')
            ENDIF
C     try to read internally
            READ (STRVAL(1:LENSTR), IOSTAT=IRES, FMT=FMTSTR)RETVAL
            IF (IRES .NE. 0) THEN
               CALL FLAGLN(PRGSTR)
               PRINT 600,STR
               PRINT 610,STRVAL(1:LENSTR)
               PRINT 620,IRES,FMTSTR
               CALL STOPGM('FILEIO','FATAL ERROR')
            ENDIF
         ENDIF

         RETARR(I) = RETVAL
 10   CONTINUE

 500  FORMAT ('(F',I3.3,'.1)')
 600  FORMAT ('failed to read a real*8 for ',A)
 610  FORMAT ('got ''',A,'''')
 620  FORMAT ('IRES =',I3,', FMTSTR = ',A)
C end CHPNRE
C     ==--------------------------------------------------------------==
      RETURN
      END

C---------------------------------
C     routines for formatted output of comments
C---------------------------------
COMMSUBR WRCOMM
C     SUBROUTINE WRCOMM(IUNIT,PRGSTR,STR)
C
C     write STR as a comment (i.e. with a preceding #)
C     to file IUNIT.
C     The open file must be formatted.
COMMEND
      SUBROUTINE WRCOMM(IUNIT,PRGSTR,STR)
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
      CHARACTER STR*(*)
C local vars
      INTEGER IRES,ILEN
      CHARACTER FMTSTR*(16)
C begin
      ILEN = LEN(STR)
C decrement until we find the first nonblank char
 10   IF (STR(ILEN:ILEN) .EQ. ' ' .AND. ILEN .GT. 1) THEN
         ILEN = ILEN -1
         GOTO 10
      ENDIF
C
      WRITE(FMTSTR,IOSTAT=IRES,FMT='(''(A1,1X,A'',I3.3,'')'')')ILEN
      IF (IRES .NE. 0) THEN
         PRINT 150,PRGSTR
         PRINT 160,STR
         CALL STOPGM('FILEIO','FATAL ERROR')
      ENDIF
      WRITE (IUNIT, IOSTAT=IRES, FMT = FMTSTR)'#',STR
      IF (IRES .NE. 0) THEN
         PRINT 250,PRGSTR
         PRINT *,FMTSTR
         CALL STOPGM('FILEIO','FATAL ERROR')
      ENDIF
 150  FORMAT (1X,'WRCOMM: format error in subroutine',A6)
 160  FORMAT (1X,'attempting to write:',A20)
 250  FORMAT (1X,'WRCOMM: write error in subroutine',A6)
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRSTR
C     SUBROUTINE WRSTR(IUNIT,PRGSTR,STR)
C
C     WRSTR writes STR to IUNIT.
C     It finds out the "true length" of the string by finding
C     the last non-blank character before writing out.
C     This prevents a lot of trailing blanks at the end of a line.
COMMEND
      SUBROUTINE WRSTR(IUNIT,PRGSTR,STR)
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
      CHARACTER STR*(*)
C local vars
      INTEGER IRES,ILEN
C      CHARACTER FMTSTR*(16)
C begin
      ILEN = LEN(STR)
C decrement until we find the first nonblack char
 10   IF (STR(ILEN:ILEN) .EQ. ' ' .AND. ILEN .GT. 1) THEN
         ILEN = ILEN -1
         GOTO 10
      ENDIF
C
      WRITE (UNIT=IUNIT, IOSTAT=IRES, FMT = 500)STR(1:ILEN)
      IF (IRES .NE. 0) THEN
         PRINT 150,PRGSTR
         PRINT 160,STR
         PRINT 170,IUNIT
         PRINT 180,IRES
         CALL STOPGM('FILEIO','FATAL ERROR')
      ENDIF
 150  FORMAT (1X,'WRSTR: write error in subroutine ',A6)
 160  FORMAT (1X,'attempting to write:',A)
 170  FORMAT (1X,'attempting to write to unit',I5)
 180  FORMAT (1X,'IOSTAT = ',I5)
 500  FORMAT (A)
C     ==--------------------------------------------------------------==
      RETURN
      END

COMMSUBR WRCINT
C     SUBROUTINE WRCINT(IUNIT,PRGSTR,K)
C
C     WRCINT writes an integer as a comment,i.e.
C     with a preceding # character.
COMMEND
      SUBROUTINE WRCINT(IUNIT,PRGSTR,K)
      INTEGER IUNIT
      CHARACTER PRGSTR*(*)
      INTEGER K
C local vars
      INTEGER IRES
C begin
      WRITE(UNIT=IUNIT,IOSTAT=IRES,FMT = '(A1,I10)')'#',K
      IF (IRES .NE. 0)THEN
         PRINT *,'WRCINT: WRITE ERROR IN',PRGSTR
         PRINT *,'ERROR WRITING TO UNIT ',IUNIT
         CALL STOPGM('FILEIO','FATAL ERROR')
      ENDIF
C end wrcint
C     ==--------------------------------------------------------------==
      RETURN
      END
