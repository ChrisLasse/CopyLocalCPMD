C     ==================================================================
Cmb - Revised on 16 Sept. 2005
COMMSUBR WRCFG
C      SUBROUTINE WRCFG(
C     $     IUNIT, LFORM,
C     $     NTIT,TITLE,
C     $     NAT,NPM,NSM,NDIM,NDR,NDJVT,
C     $     LWRX,X,
C     $     LWRV,V,
C     $     LWRB,BOX,
C     $     LWRXC,XC,
C     $     LWRSTC,SX,IG,
C     $     LWRPER,NDSZE,RLAM,DLSUM,RMU,DMSUM,
C     $     LWRDR,RIIAVE,
C     $     LWRJV,COSQAV,COSIAV,
C     $     LWRLE,NDLE,NLECFG,NLEMEM,ILEMEM,NLEVST)
C
C     WRCFG writes a configuration and other assorted data
C     needed for restarting an MD run to a file IUNIT
C     that must be open BEFORE calling this routine.
C
C     The file is NOT closed by this routine after writing.
C     This makes it easy for developers to add blocks to the end of
C     a startup file after having called WRCFG.
C     Likewise, their blocks can be read after calling L<RDCFG>.
C
C     WRCFG calls the following routines:
C     L<WRTIT>    : writes a title block
C     L<WRFMT>    : writes a format block
C     L<WRXPOS>   : writes a POSITION block (If LWRX = .TRUE.)
C     L<WRXVEL>   : writes a VELOCITY block (If LWRV = .TRUE.)
C     L<WRBOX>    : writes a BOX block (If LWRBOX = .TRUE.)
C     L<WRXC>     : writes a position re(con)straining reference
C                   position block (If LWRXC = .TRUE.)
C     L<WRIIAV>   : writes a distance retraining average block
C                   (If (LWRDR = .TRUE.)
C     L<WRSX>     : writes a stochastic integral block
C                   (If LWRSTC = .TRUE.)
C     L<WRPER>    : writes a perturbation block (If LWRPER = .TRUE.)
C     L<WRJVAV>   : writes a J-value average block (If LWRJV = .TRUE.)
C     L<WRLEME>   : writes a local elevation memory block
C                   (If LWRLE = .TRUE.)
COMMEND
C     ==--------------------------------------------------------------==
      SUBROUTINE WRCFG(
     $     IUNIT, LFORM,
     $     NTIT,TITLE,
     $     NAT,NPM,NSM,NDIM,NDRT,NDJVT,
     $     LWRX,X,
     $     LWRV,V,
     $     LWRB,BOX,
     $     LWRXC,XC,
     $     LWRSTC,SX,IG,
     $     LWRPER,NDSZE,RLAM,DLSUM,RMU,DMSUM,
     $     LWRDR,RIIAVE,
     $     LWRJV,COSQAV,COSIAV,
     $     LWRLE,NDLE,NLECFG,NLEMEM,ILEMEM,NLEVST)
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
C args
      INTEGER IUNIT,NAT,NPM,NSM,NDIM,NDRT,NDJVT,NTIT,NDSZE
      CHARACTER TITLE*(*)
      DIMENSION TITLE(NTIT)
      LOGICAL LFORM, LWRX,LWRV,LWRB,LWRXC,LWRDR
      LOGICAL LWRSTC,LWRPER,LWRJV,LWRLE
C
      real*8 X(NDIM*NAT),V(NDIM*NAT),BOX(NDIM),XC(NDIM*NAT)
      real*8 SX(NDIM*NAT)
      INTEGER IG
      real*8 RLAM,DLSUM(NDSZE)
      real*8 RMU, DMSUM(NDSZE)
      real*8 RIIAVE(NDRT)
      real*8 COSQAV(NDJVT),COSIAV(NDJVT)
      INTEGER NDLE,NLECFG,NLEMEM
      INTEGER ILEMEM(MXLECF,MLECFG)
      INTEGER NLEVST(MLECFG)
C local vars
      LOGICAL LWR4DX,LWR4DV
      PARAMETER (LWR4DX = .TRUE., LWR4DV = .TRUE.)
C begin
      CALL WRTIT(IUNIT,LFORM,NTIT,TITLE)

      CALL WRFMT(IUNIT,LFORM)

      IF (LWRX) THEN
         CALL WRXPOS(IUNIT,LFORM,LWR4DX,NAT,NPM,NSM,NDIM,X)
      ENDIF
C
      IF (LWRV) THEN
         CALL WRXVEL(IUNIT,LFORM,LWR4DV,NAT,NPM,NSM,NDIM,V)
      ENDIF
C
      IF (LWRB) THEN
         CALL WRBOX(IUNIT,LFORM,BOX)
      ENDIF
C
      IF (LWRXC) THEN
         CALL WRXC(IUNIT,LFORM,NAT,NPM,NSM,NDIM,XC)
      ENDIF

      IF (LWRDR) THEN
         CALL WRIIAV(IUNIT,LFORM,NDRT,RIIAVE)
      ENDIF

      IF (LWRSTC) THEN
         CALL WRSX(IUNIT,LFORM,LWR4DV,NAT,NPM,NSM,NDIM,SX,IG)
      ENDIF

      IF (LWRPER) THEN
         CALL WRPER(IUNIT,LFORM,NDSZE,RLAM,DLSUM,RMU,DMSUM)
      ENDIF

      IF (LWRJV) THEN
         CALL WRJVAV(IUNIT,LFORM,NDJVT,COSQAV,COSIAV)
      ENDIF

      IF (LWRLE) THEN
         CALL WRLEME(IUNIT,LFORM,NDLE,NLECFG,NLEMEM,ILEMEM,NLEVST)
      ENDIF
C end WRCFG
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
COMMSUBR RDCFG
C     SUBROUTINE RDCFG(IUNIT, LFORM,
C     $   MAXTIT,NTIT,TITLE,
C     $   LRD4DX,LRD4DV,
C     $   NAT,NDIM,NDRT,NDJVT,
C     $   LRDX,X,
C     $   LRDV,V,
C     $   LRDB,BOX,
C     $   LRDXC,XC,
C     $   LRDSTC,SX,IG,
C     $   LRDPER,NDSZE,RLAM,DLSUM,RMU,DMSUM,
C     $   LRDDIS,RIIAVE,
C     $   LRDJV,COSQAV,COSIAV,
C     $   LRDLE,NDLE,NLECFG,NLEMEM,ILEMEM,NLEVST)
C
C     RDCFG reads a configuration and other assorted
C     data needed to restart an MD run from a file IUNIT
C     that must be open BEFORE calling this routine.
C
C     The file is NOT closed by this routine after reading.
C     This makes it easy for developers to add blocks to the end of
C     a startup file and read them after calling RDCFG.
C     Likewise, their blocks can be written after calling L<WRCFG>.
C     
C     RDCFG calls the following routines:
C     L<RDBHDR>     : read a block header
C     L<RDTIT>      : read a title block
C     L<RDXPOS>     : read a position block (If LRDX = .TRUE.)
C     L<RDXRED>     : read a reduced position block (If LRDX = .TRUE.)
C     L<RDXVEL>     : read a velocity block (If LRDV = .TRUE.)
C     L<RDVRED>     : read a reduced velocity block (If LRDV = .TRUE.)
C     L<RDBOX>      : read the L<BOX> dimensions (If LRDB = .TRUE.)
C     L<RDXC>       : read position re(con)straints position
C                     block (If LRDXC = .TRUE.)
C     L<RDSX>       : read in stochastic integrals (If LRDSTC = .TRUE.)
C     L<RDPER>      : read in a perturbation block (If LRDPER = .TRUE.)
C     L<RDIIAV>     : read a distance restraint average block
C                     (If LRDDIS = .TRUE.)
C     L<RDJVAV>     : read a J-value average block (If LRDJV = .TRUE.)
C     L<RDLEME>     : read local memory block (If LRDLE = .TRUE.)
C
C     L<SKPBLK>     : skip over unwanted blocks on file
C     L<FLAGLN>     : write out line number in case of errors
COMMEND
C     ==--------------------------------------------------------------==
      SUBROUTINE RDCFG(
     $     IUNIT, LFORM,
     $     MAXTIT,NTIT,TITLE,
     $     LRD4DX,LRD4DV,
     $     NAT,NDIM,NDRT,NDJVT,
     $     LRDX,X,
     $     LRDV,V,
     $     LRDB,BOX,
     $     LRDXC,XC,
     $     LRDSTC,SX,IG,
     $     LRDPER,NDSZE,RLAM,DLSUM,RMU,DMSUM,
     $     LRDDIS,RIIAVE,
     $     LRDJV,COSQAV,COSIAV,
     $     LRDLE,NDLE,NLECFG,NLEMEM,ILEMEM,NLEVST)
C     ==--------------------------------------------------------------==
C include
      INCLUDE 'dataid.h'
      INCLUDE 'coordsz.h'
C args
      INTEGER IUNIT,NAT,NDIM,NDRT,IG,NDJVT,MAXTIT,NTIT,NDSZE
      CHARACTER TITLE*(*)
      DIMENSION TITLE(MAXTIT)
      LOGICAL LFORM, LRDX,LRDV,LRDB,LRDXC,LRDDIS,LRDSTC,LRDPER
      LOGICAL LRD4DX,LRD4DV,LRDJV,LRDLE
      real*8 X(NDIM*NAT),V(NDIM*NAT),BOX(NDIM),XC(NDIM*NAT)
      real*8 RIIAVE(NDRT)
      real*8 SX(NDIM*NAT)
      real*8 RLAM,DLSUM(NDSZE)
      real*8 RMU,DMSUM(NDSZE)
      real*8 COSQAV(NDJVT),COSIAV(NDJVT)

      INTEGER NDLE,NLECFG,NLEMEM
      INTEGER ILEMEM(MXLECF,MLECFG)
      INTEGER NLEVST(MLECFG)

C local vars
      LOGICAL LGOTIT,LGOTPO,LGOTVE,LGOTBO,LGOTXC,LGOTDR
      LOGICAL LGOTSX,LGOTPE,LGOTJV,LGOTLE
      LOGICAL LGTALL
      INTEGER IDBLK,NREC
      CHARACTER PRGSTR*(5)
C data
      DATA PRGSTR /'RDCFG'/
C begin
      LGOTIT = .FALSE.
      LGOTPO = .FALSE.
      LGOTVE = .FALSE.
      LGOTBO = .FALSE.
      LGOTXC = .FALSE.
      LGOTDR = .FALSE.
      LGOTSX = .FALSE.
      LGOTPE = .FALSE.
      LGOTJV = .FALSE.
      LGOTLE = .FALSE.

      LGTALL = .FALSE.


C first block must be a title
      write(6,*)IUNIT
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE. ITITID) THEN
         IF (LFORM) THEN
            CALL FLAGLN(PRGSTR)
         ENDIF
         PRINT FMEXP,PRGSTR,NAMEID(ITITID),IUNIT
         CALL STOPGM('GDCFG','FATAL ERROR')
      ENDIF
      CALL RDTIT(IUNIT,LFORM,NREC,MAXTIT,TITLE,NTIT)
      LGOTIT = .TRUE.


      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
C while not LGTALL and not an error reading
 10   IF (.NOT. LGTALL .AND. IDBLK .NE. IDERR) THEN
C pos
         IF (IDBLK .EQ. IPOSID) THEN
            IF (LRDX) THEN
               IF (LGOTPO) THEN
                  IF (LFORM) THEN
                     CALL FLAGLN(PRGSTR)
                  ENDIF
                  PRINT 550,PRGSTR,NAMEID(IDBLK)
                  CALL STOPGM('GDCFG','FATAL ERROR')
               ELSE
                  CALL RDXPOS(IUNIT,LFORM,LRD4DX,NREC,NAT,NDIM,X)
                  LGOTPO = .TRUE.
               ENDIF
            ELSE
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
               CALL SKPBLK(IUNIT,LFORM,NREC)
            ENDIF
C reduced position
         ELSEIF (IDBLK .EQ. IPORID) THEN
            IF (LRDX) THEN
               IF (LGOTPO) THEN
                  IF (LFORM) THEN
                     CALL FLAGLN(PRGSTR)
                  ENDIF
                  PRINT 550,PRGSTR,NAMEID(IDBLK)
                  CALL STOPGM('GDCFG','FATAL ERROR')
               ELSE
                  CALL RDXRED(IUNIT,LFORM,LRD4DX,NREC,NAT,NDIM,X)
                  LGOTPO = .TRUE.
               ENDIF
            ELSE
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
               CALL SKPBLK(IUNIT,LFORM,NREC)
            ENDIF
C vel
         ELSEIF (IDBLK .EQ. IVELID) THEN
            IF (LRDV) THEN
               IF (LGOTVE) THEN
                  IF (LFORM) THEN
                     CALL FLAGLN(PRGSTR)
                  ENDIF
                  PRINT 550,PRGSTR,NAMEID(IDBLK)
                  CALL STOPGM('GDCFG','FATAL ERROR')
               ELSE
                  CALL RDXVEL(IUNIT,LFORM,LRD4DV,NREC,NAT,NDIM,V)
                  LGOTVE = .TRUE.
               ENDIF
            ELSE
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
               CALL SKPBLK(IUNIT,LFORM,NREC)
            ENDIF
C velocity reduced
         ELSEIF (IDBLK .EQ. IVRDID) THEN
            IF (LRDV) THEN
               IF (LGOTVE) THEN
                  IF (LFORM) THEN
                     CALL FLAGLN(PRGSTR)
                  ENDIF
                  PRINT 550,PRGSTR,NAMEID(IDBLK)
                  CALL STOPGM('GDCFG','FATAL ERROR')
               ELSE
                  CALL RDVRED(IUNIT,LFORM,LRD4DV,NREC,NAT,NDIM,V)
                  LGOTVE = .TRUE.
               ENDIF
            ELSE
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
               CALL SKPBLK(IUNIT,LFORM,NREC)
            ENDIF
C box
         ELSEIF (IDBLK .EQ. IBOXID) THEN
            IF (LRDB) THEN
               IF (LGOTBO) THEN
                  IF (LFORM) THEN
                     CALL FLAGLN(PRGSTR)
                  ENDIF
                  PRINT 550,PRGSTR,NAMEID(IDBLK)
                  CALL STOPGM('GDCFG','FATAL ERROR')
               ELSE
                  CALL RDBOX(IUNIT,LFORM,NREC,BOX)
                  LGOTBO = .TRUE.
               ENDIF
            ELSE
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
               CALL SKPBLK(IUNIT,LFORM,NREC)
            ENDIF
C pos re(con)straints coordinates
         ELSEIF (IDBLK .EQ. IPRCID) THEN
            IF (LRDXC) THEN
               IF (LGOTXC) THEN
                  IF (LFORM) THEN
                     CALL FLAGLN(PRGSTR)
                  ENDIF
                  PRINT 550,PRGSTR,NAMEID(IDBLK)
                  CALL STOPGM('GDCFG','FATAL ERROR')
               ELSE
                  CALL RDXC(IUNIT,LFORM,NREC,NAT,NDIM,XC)
                  LGOTXC = .TRUE.
               ENDIF
            ELSE
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
               CALL SKPBLK(IUNIT,LFORM,NREC)
            ENDIF
C stochastic integrals
         ELSEIF (IDBLK .EQ. ISDIID) THEN
            IF (LRDSTC) THEN
               IF (LGOTSX) THEN
                  IF (LFORM) THEN
                     CALL FLAGLN(PRGSTR)
                  ENDIF
                  PRINT 550,PRGSTR,NAMEID(IDBLK)
                  CALL STOPGM('GDCFG','FATAL ERROR')
               ELSE
                  CALL RDSX(IUNIT,LFORM,.TRUE.,NREC,NAT,NDIM,SX,IG)
                  LGOTSX = .TRUE.
               ENDIF
            ELSE
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
               CALL SKPBLK(IUNIT,LFORM,NREC)
            ENDIF
C perturbation
         ELSEIF (IDBLK .EQ. IPRTID) THEN
            IF (LRDPER) THEN
               IF (LGOTPE) THEN
                  IF (LFORM) THEN
                     CALL FLAGLN(PRGSTR)
                  ENDIF
                  PRINT 550,PRGSTR,NAMEID(IDBLK)
                  CALL STOPGM('GDCFG','FATAL ERROR')
               ELSE
                  CALL RDPER(IUNIT,LFORM,NREC,
     $                 NDSZE,RLAM,DLSUM,RMU,DMSUM)
                  LGOTPE = .TRUE.
               ENDIF
            ELSE
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
               CALL SKPBLK(IUNIT,LFORM,NREC)
            ENDIF

C     NOTE: 4D coords will be read in by RDVEL, or RDPOS...
C     if the calling routine specified that with
C     LRD4DX,LRD4DV, = .TRUE. .
C     If LRD4DV is .FALSE., however, then we will encounter
C     these block types here. Just notify the user we have
C     seen them on file and not read them, then skip them...
C vel or pos in 4D
         ELSEIF (NDIM .GT. 3 .AND. 
     $           (((IDBLK .EQ. IP4DID .OR. IDBLK .EQ. IP4RID)
     $           .AND. LRD4DX) .OR.
     $           (IDBLK .EQ. IV4DID .OR. IDBLK .EQ. IV4RID)
     $           .AND. LRD4DV)) THEN
            IF (IDBLK .GE. MINIDT .AND. IDBLK .LE. MAXIDT) THEN
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
            ELSEIF (IDBLK .EQ. IDUKN) THEN
               PRINT 500,PRGSTR,IDLN,IUNIT
            ENDIF
            CALL SKPBLK(IUNIT,LFORM,NREC)

C distance restraint values
         ELSEIF (IDBLK .EQ. IRIIAV) THEN
            IF (LRDDIS) THEN
               IF (LGOTDR) THEN
                  IF (LFORM) THEN
                     CALL FLAGLN(PRGSTR)
                  ENDIF
                  PRINT 550,PRGSTR,NAMEID(IDBLK)
                  CALL STOPGM('GDCFG','FATAL ERROR')
               ELSE
                  CALL RDIIAV(IUNIT,LFORM,NREC,NDRT,RIIAVE)
                  LGOTDR = .TRUE.
               ENDIF
            ELSE
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
               CALL SKPBLK(IUNIT,LFORM,NREC)
            ENDIF

C j-val restraint values
         ELSEIF (IDBLK .EQ. ICOSAV) THEN
            IF (LRDJV) THEN
               IF (LGOTJV) THEN
                  IF (LFORM) THEN
                     CALL FLAGLN(PRGSTR)
                  ENDIF
                  PRINT 550,PRGSTR,NAMEID(IDBLK)
                  CALL STOPGM('GDCFG','FATAL ERROR')
               ELSE
                  CALL RDJVAV(IUNIT,LFORM,NREC,NDJVT,COSQAV,COSIAV)
                  LGOTJV = .TRUE.
               ENDIF
            ELSE
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
               CALL SKPBLK(IUNIT,LFORM,NREC)
            ENDIF
C local elevation memory
         ELSEIF (IDBLK .EQ. ILEMID) THEN
            IF (LRDLE) THEN
               IF (LGOTLE) THEN
                  IF (LFORM) THEN
                     CALL FLAGLN(PRGSTR)
                  ENDIF
                  PRINT 550,PRGSTR,NAMEID(IDBLK)
                  CALL STOPGM('GDCFG','FATAL ERROR')
               ELSE
                  CALL RDLEME(IUNIT,LFORM,NREC,
     $                 NDLE,NLECFG,NLEMEM,ILEMEM,NLEVST)
                  LGOTLE = .TRUE.
               ENDIF
            ELSE
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
               CALL SKPBLK(IUNIT,LFORM,NREC)
            ENDIF
         ELSE
C we have a block type we dont know about
            IF (IDBLK .GE. MINIDT .AND. IDBLK .LE. MAXIDT) THEN
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
            ELSEIF (IDBLK .EQ. IDUKN) THEN
               PRINT 500,PRGSTR,IDLN,IUNIT
            ENDIF
            CALL SKPBLK(IUNIT,LFORM,NREC)
         ENDIF

C check for termination conditions
         LGTALL = LGOTIT
         LGTALL = (LGTALL .AND. .NOT. (LRDX .AND. .NOT. LGOTPO))
         LGTALL = (LGTALL .AND. .NOT. (LRDV .AND. .NOT. LGOTVE))
         LGTALL = (LGTALL .AND. .NOT. (LRDB .AND. .NOT. LGOTBO))
         LGTALL = (LGTALL .AND. .NOT. (LRDXC .AND. .NOT. LGOTXC))
         LGTALL = (LGTALL .AND. .NOT. (LRDDIS .AND. .NOT. LGOTDR))
         LGTALL = (LGTALL .AND. .NOT. (LRDSTC .AND. .NOT. LGOTSX))
         LGTALL = (LGTALL .AND. .NOT. (LRDPER .AND. .NOT. LGOTPE))
         LGTALL = (LGTALL .AND. .NOT. (LRDJV .AND. .NOT. LGOTJV))
         LGTALL = (LGTALL .AND. .NOT. (LRDLE .AND. .NOT. LGOTLE))

         CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
C end while
         GOTO 10
      ENDIF

      IF (.NOT. LGTALL) THEN
         PRINT *,'RDCFG: error reading from UNIT ',IUNIT
         IF (LRDX .AND. .NOT. LGOTPO) THEN
            PRINT 560,'positions'
         ENDIF
         IF (LRDV .AND. .NOT. LGOTVE) THEN
            PRINT 560,'velocities'
         ENDIF
         IF (LRDB .AND. .NOT. LGOTBO) THEN
            PRINT 560,'box'
         ENDIF
         IF (LRDXC .AND. .NOT. LGOTXC) THEN
            PRINT 560,'position restraint positions'
         ENDIF
         IF (LRDDIS .AND. .NOT. LGOTDR) THEN
            PRINT 560,'time averaged distances'
         ENDIF
         IF (LRDLE .AND. .NOT. LGOTLE) THEN
            PRINT 560,'local elevation memory'
         ENDIF
         
         CALL STOPGM('GDCFG','FATAL ERROR')
      ENDIF

 500  FORMAT (1X,A,': BLOCK ',A18,'on UNIT ',I3,' SKIPPED')
 550  FORMAT (1X,A6,': ',A10,'block encountered twice!')
 560  FORMAT (1X,A,' required, but not on file!')
C end RDCFG
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
COMMSUBR RDXRST
C     SUBROUTINE RDXRST(NATTOT,NDIM,LRDCOR,LRDBFA)
C
C     RDXRST opens a file, reads and prints the title block,
C     then, if LRDCOR = .TRUE., reads a position re(con)straints
C     POSITION block into the array L<XC> defined in coordar.h
C     using subroutine L<RDXC>.
C
C     IF LRDBFA = .TRUE., then L<RDABF> is also called to read
C     in an isotropic atomic B factors block into the array
C     L<CXC> defined in coordar.h.
C
C     Note that in GROMOS, there are TWO block types
C     needed for position restraining.
C     a) a position restraints POSITION block
C        this defines the restraining reference positions
C        for EVERY atom in the system, i.e. this must contain a
C        complete set of coordinates.
C     b) a position restraints specification block
C        this defines which of the atoms is actually to
C        be re(con)strained.
C        This block contains only those atoms which are to be
C        re(con)strained.
C     An a) block can be in a separate file OR in the start
C     up configuration. This is because reference positions
C     can change e.g. when pressure coupling is on.
C     A b) block is always in a separate file for itself
C     and is therefore always read in by L<GTIRST>.
C
C     L<RDXRST> is only used if position re(con)straint positions
C     are on a separate file for themselves.
C     They can also be on a file with initial coordinates
C     and velocities, in which case they are read in by
C     subroutine L<RDCFG>.
C
C     NATTOT and NDIM are not modified.
COMMEND
C     ==--------------------------------------------------------------==
      SUBROUTINE RDXRST(NATTOT,NDIM,LRDCOR,LRDBFA)
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
C args
      INTEGER NATTOT,NDIM
      LOGICAL LRDCOR,LRDBFA
C local params
      INTEGER MAXTIT
      PARAMETER (MAXTIT = 10)
C local vars
      INTEGER IUNIT,IDBLK,NREC,NTIT
      LOGICAL LFORM
      LOGICAL LGCOR,LGBFA
      LOGICAL LOK
      CHARACTER PRGSTR*(6)
      CHARACTER TITLE*(MAXDLN)
      DIMENSION TITLE(MAXTIT)
C data
      DATA PRGSTR /'RDXRST'/
C begin
      IF (.NOT. LRDCOR .AND. .NOT. LRDBFA) THEN
         RETURN
      ENDIF

      CALL DETOPN('RSTCOR',IUNIT,LFORM)
      IF (IUNIT .LT. 0) THEN
         PRINT *,PRGSTR,': failed to open file for ',
     $      'reading position restraining coordinates',
     $      'and/or atomic B factors'
         CALL STOPGM('GDXRST','FATAL ERROR')
      ENDIF
      CALL INITLN(IUNIT)

C get title and print it
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE. ITITID ) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMEXP,PRGSTR,NAMEID(ITITID),IUNIT
         CALL STOPGM('GDXRST','FATAL ERROR')
      ENDIF
      CALL RDTIT(IUNIT,LFORM,NREC,MAXTIT,TITLE,NTIT)
      CALL PRTIT('POSITION RESTRAINING COORDINATES/BFA TITLE',
     $     NTIT,TITLE)

      LOK = .FALSE.
      LGCOR = .FALSE.
      LGBFA = .FALSE.

 10   IF (.NOT. LOK) THEN
         CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)

         IF (IDBLK .EQ. IPRCID) THEN
C     get coords into XC
            IF (LRDCOR) THEN
               IF (LGCOR) THEN
                  CALL FLAGLN(PRGSTR)
                  PRINT FMG2,PRGSTR,NAMEID(IDBLK),IUNIT
                  CALL STOPGM('GDXRST','FATAL ERROR')
               ELSE
                  CALL RDXC(IUNIT,LFORM,NREC,NATTOT,NDIM,XC)
                  LGCOR = .TRUE.
               ENDIF
            ELSE
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
               CALL SKPBLK(IUNIT,LFORM,NREC)
            ENDIF
         ELSEIF (IDBLK .EQ. IABFID) THEN
            IF (LRDBFA) THEN
               IF (LGBFA) THEN
                  CALL FLAGLN(PRGSTR)
                  PRINT FMG2,PRGSTR,NAMEID(IDBLK),IUNIT
                  CALL STOPGM('GDXRST','FATAL ERROR')
               ELSE
                  CALL RDABF(IUNIT,LFORM,NREC,NATTOT,CXC)
                  LGBFA = .TRUE.
               ENDIF
            ELSE
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
               CALL SKPBLK(IUNIT,LFORM,NREC)
            ENDIF
         ELSEIF (IDBLK .EQ. IDERR) THEN
            CALL FLAGLN(PRGSTR)
            IF (LRDBFA .AND. .NOT. LGBFA) THEN
               PRINT FMEXP,PRGSTR,NAMEID(IABFID),IUNIT
               PRINT 510,NAMEID(IABFID)
            ENDIF
            IF (LRDCOR .AND. .NOT. LGCOR) THEN
               PRINT FMEXP,PRGSTR,NAMEID(IPRCID),IUNIT
               PRINT 510,NAMEID(IPRCID)
            ENDIF
            CALL STOPGM('GDXRST','FATAL ERROR')
         ELSE
C some other block.. just skip
            IF (IDBLK .GE. MINIDT .AND. IDBLK .LE. MAXIDT) THEN
               PRINT 500,PRGSTR,NAMEID(IDBLK),IUNIT
            ELSEIF (IDBLK .EQ. IDUKN) THEN
               PRINT 500,PRGSTR,IDLN,IUNIT
            ENDIF
            CALL SKPBLK(IUNIT,LFORM,NREC)
         ENDIF

         LOK = (.NOT. LRDCOR .OR. (LRDCOR .AND. LGCOR))
         LOK = LOK .AND. (.NOT. LRDBFA .OR. (LRDBFA .AND. LGBFA))

         GOTO 10
C end of while loop
      ENDIF

      CALL CLSFIL(IUNIT)

 500  FORMAT (1X,A,': skipping block ',A,' on unit',I4)
 510  FORMAT (1X,'Looking for a ',A,'block, but none on file')
C end RDXRST
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
COMMSUBR GTIRST
C     SUBROUTINE GTIRST(NLIM,MAXXC,JRC,NRC)
C
C     GTIRST opens a file and reads and prints the title block.
C     It then reads in the indices of the atoms
C     to be re(con)strained into JRC using subroutine L<RDNDXC>.
C     The variables NRC and JRC are modified.
C     The variables NLIM and MAXXC are not modified.
C
C     The difference between reading CONstraining and REstraining
C     indices is in the allowed range of the indices:
C     for constraining:  1.. NRP    (the number of atoms in a solute molecule)
C     for restraining:   1.. NATTOT (the total number of atoms in the system)
C     The appropriate value should be passed in NLIM.
COMMEND
C     ==--------------------------------------------------------------==
      SUBROUTINE GTIRST(NLIM,MAXXC,JRC,NRC)
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
C args
      INTEGER NLIM,MAXXC,JRC(MAXXC),NRC
C local params
      INTEGER MAXTIT
      PARAMETER (MAXTIT = 10)
C local vars
      INTEGER IUNIT,IDBLK,NREC,NTIT
      LOGICAL LFORM
      CHARACTER PRGSTR*(6)
      CHARACTER TITLE*(MAXDLN)
      DIMENSION TITLE(MAXTIT)
C data
      DATA PRGSTR /'GTIRST'/
C begin
      CALL DETOPN('RSTNDX',IUNIT,LFORM)
      IF (IUNIT .LT. 0) THEN
         PRINT *,PRGSTR,': failed to open file for ',
     $        'reading atom sequence numbers for',
     $        ' position re(con)straining'
         CALL STOPGM('GTIRST','FATAL ERROR')
      ENDIF
      CALL INITLN(IUNIT)

C get title and print it
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE. ITITID ) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMEXP,PRGSTR,NAMEID(ITITID),IUNIT
         PRINT 520
         CALL STOPGM('GTIRST','FATAL ERROR')
      ENDIF
      CALL RDTIT(IUNIT,LFORM,NREC,MAXTIT,TITLE,NTIT)
      CALL PRTIT('POSITION RE(CON)STRAINING INDICES TITLE',
     $     NTIT,TITLE)

C get atom indices 
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE. IXPRID) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMEXP,PRGSTR,NAMEID(IXPRID),IUNIT
         PRINT 520
         CALL STOPGM('GTIRST','FATAL ERROR')
      ENDIF
      CALL RDNDXC(IUNIT,LFORM,NLIM,MAXXC,JRC,NRC)
      CALL CLSFIL(IUNIT)

      PRINT 510,NRC
      PRINT *

 510  FORMAT (1X,'TOTAL OF',I6,' POSITION RE(CON)STRAINED ATOMS')
 520  FORMAT (1X,'Failed to read atom sequence numbers for',
     $     ' position re(con)straining',/)
C end GTIRST
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
COMMSUBR GT4NDX
C     SUBROUTINE GT4NDX(NATTOT,C4D)
C
C     GT4NDX opens a file, reads and prints the title block
C     and then reads the information from file which
C     specifies the harmonic force constants in the
C     fourth dimension using subroutine L<RD4NDX>.
C     See variable L<C4D>.
C
C     NATTOT  specifies how many atoms are in the system as
C     a whole. NATTOT is not modified, but used to check
C     consistency of the file read in (there must be exactly
C     NATTOT values in the block read in).
C
C     C4D(NATTOT) is delivered with the numbers
C     read from file.
C
C     NOTE: NO SCALING OF THE C4D VALUES IS PERFORMED
C     IN THIS ROUTINE.
COMMEND
C     ==--------------------------------------------------------------==
      SUBROUTINE GT4NDX(NATTOT,C4D)
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
      INCLUDE 'coordsz.h'
C args
      INTEGER NATTOT
      real*8 C4D(NATTOT)
C local params
      INTEGER MAXTIT
      PARAMETER (MAXTIT = 10)
C local vars
      INTEGER IUNIT,IDBLK,NREC,NTIT
      LOGICAL LFORM
      CHARACTER PRGSTR*(6)
      CHARACTER TITLE*(MAXDLN)
      DIMENSION TITLE(MAXTIT)
C data
      DATA PRGSTR /'GT4NDX'/
C begin
      CALL DETOPN('4DNDX',IUNIT,LFORM)
      IF (IUNIT .LT. 0) THEN
         PRINT 500,PRGSTR
         PRINT *,'failed to open file'
         CALL STOPGM('GT4NDX','FATAL ERROR')
      ENDIF
      CALL INITLN(IUNIT)

C     get title and print it
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE. ITITID ) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMEXP,PRGSTR,NAMEID(ITITID),IUNIT
         CALL STOPGM('GT4NDX','FATAL ERROR')
      ENDIF
      CALL RDTIT(IUNIT,LFORM,NREC,MAXTIT,TITLE,NTIT)
      CALL PRTIT('4D ATOM SEQUENCE NUMBERS TITLE',NTIT,TITLE)

C     get atom indices
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE. IX4DID) THEN
         CALL FLAGLN(PRGSTR)
         PRINT 500,PRGSTR
         PRINT FMEXP,PRGSTR,NAMEID(IX4DID),IUNIT
         CALL STOPGM('GT4NDX','FATAL ERROR')
      ENDIF
      CALL RD4NDX(IUNIT,LFORM,NREC,NATTOT,C4D)
      CALL CLSFIL(IUNIT)

 500  FORMAT (1X,A6,': reading 4D atom sequence numbers')
 510  FORMAT (1X,A80)
C end gt4ndx
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
COMMSUBR GTDRST
C     SUBROUTINE GTDRST(NRPT)
C
C     Read in the distance restraint specification data.
C     NRPT : total number of solute atoms. Not modified.
C
C     GTDRST opens a file, reads and print the title block
C     and then reads in the distance restraint specification data
C     L<NDR>,
C     L<DISH>, L<DISC>,
C     L<IDR1>, L<JDR1, L<KDR1>, L<LDR1>, L<ICDR1>,
C     L<IDR2>, L<JDR2, L<KDR2>, L<LDR2>, L<ICDR2>,
C     L<R0>, L<W0>
C     using subroutine L<RDDRST>.
C     It then closes the file.
C
C     NRPT (= total number of solute atoms)
C     is not modified, but is used to check
C     L<IDR1>, L<JDR1, L<KDR1>, L<LDR1>, L<ICDR1>,
C     L<IDR2>, L<JDR2, L<KDR2>, L<LDR2> and L<ICDR2>
C     on reading.
COMMEND
C     ==--------------------------------------------------------------==
      SUBROUTINE GTDRST(NRPT)
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
C args
      INTEGER NRPT
C local params
      INTEGER MAXTIT
      PARAMETER (MAXTIT = 10)
C local vars
      INTEGER IUNIT,IDBLK,NREC,NTIT
      LOGICAL LFORM
      CHARACTER PRGSTR*(6)
      CHARACTER TITLE*(MAXDLN)
      DIMENSION TITLE(MAXTIT)
C data
      DATA PRGSTR /'GTDRST'/
C begin
      CALL DETOPN('RSTPAIR',IUNIT,LFORM)
      IF (IUNIT .LT. 0) THEN
         PRINT 500,PRGSTR
         PRINT *,'failed to open distance restraint file'
         CALL STOPGM('GTDRST','FATAL ERROR')
      ENDIF
      IF (.NOT. LFORM) THEN
         PRINT 500,PRGSTR
         PRINT *,'file must be formatted'
         CALL STOPGM('GTDRST','FATAL ERROR')
      ENDIF
      CALL INITLN(IUNIT)

C     get title and print it
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE. ITITID ) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMEXP,PRGSTR,NAMEID(ITITID),IUNIT
         CALL STOPGM('GTDRST','FATAL ERROR')
      ENDIF
      CALL RDTIT(IUNIT,LFORM,NREC,MAXTIT,TITLE,NTIT)
      CALL PRTIT('DISTANCE RESTRAINTS TITLE',NTIT,TITLE)

C     get disre atom indices
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE. IDRBLK) THEN
         CALL FLAGLN(PRGSTR)
         PRINT 500,PRGSTR
         PRINT FMEXP,PRGSTR,NAMEID(IDRBLK),IUNIT
         CALL STOPGM('GTDRST','FATAL ERROR')
      ENDIF

      CALL RDDRST(IUNIT,NRPT,MAXNDR,NDR,
     $     DISH,DISC,
     $     IDR1,JDR1,KDR1,LDR1,ICDR1,
     $     IDR2,JDR2,KDR2,LDR2,ICDR2,
     $     R0,W0)

      CALL CLSFIL(IUNIT)

      PRINT *,'READ ',NDR,' DISTANCE RESTRAINTS'
      PRINT *

 500  FORMAT (1X,A6,': reading distance restraint information')
C end gtdrst
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
COMMSUBR GTDHRS
C     SUBROUTINE GTDHRS(NRPT)
C
C     GTDHRS opens a file, reads and print the title block
C     and then reads in the dihedral restraint specification data
C     <NDLR>,
C     L<IPLR>, L<JPLR>, L<KPLR>, L<LPLR>,
C     L<CPLR>, L<PDLR> and L<ICPLRL> defined in coordar.h
C     using subroutine L<RDDHRS>.
C     It then closes the file.
C
C     NRPT (= total number of solute atoms)
C     is not modified, but is used to check
C     L<IPLR>, L<JPLR>, L<KPLR> and L<LPLR> on reading.
C
C     NO CONVERSION OF UNITS IS DONE IN THIS SUBROUTINE.
COMMEND
C     ==--------------------------------------------------------------==
      SUBROUTINE GTDHRS(NRPT)
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
C args
      INTEGER NRPT
C local params
      INTEGER MAXTIT
      PARAMETER (MAXTIT = 10)
C local vars
      INTEGER IUNIT,IDBLK,NREC,NTIT
      LOGICAL LFORM
      CHARACTER PRGSTR*(6)
      CHARACTER TITLE*(MAXDLN)
      DIMENSION TITLE(MAXTIT)
C data
      DATA PRGSTR /'GTDHRS'/
C begin
      CALL DETOPN('DIHESPEC',IUNIT,LFORM)
      IF (IUNIT .LT. 0) THEN
         PRINT 500,PRGSTR
         PRINT *,'failed to open dihedral restraint file'
         CALL STOPGM('GTDHRS','FATAL ERROR')
      ENDIF
      IF (.NOT. LFORM) THEN
         PRINT 500,PRGSTR
         PRINT *,'file must be formatted'
         CALL STOPGM('GTDHRS','FATAL ERROR')
      ENDIF
      CALL INITLN(IUNIT)

C     get title and print it
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE. ITITID ) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMEXP,PRGSTR,NAMEID(ITITID),IUNIT
         CALL STOPGM('GTDHRS','FATAL ERROR')
      ENDIF
      CALL RDTIT(IUNIT,LFORM,NREC,MAXTIT,TITLE,NTIT)
      CALL PRTIT('DIHEDRAL RESTRAINTS TITLE',NTIT,TITLE)

C     get dihedral restraint data
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE. IDHBLK) THEN
         CALL FLAGLN(PRGSTR)
         PRINT 500,PRGSTR
         PRINT FMEXP,PRGSTR,NAMEID(IDHBLK),IUNIT
         CALL STOPGM('GTDHRS','FATAL ERROR')
      ENDIF

      CALL RDDHRS(IUNIT,NRPT,MAXNDL,NDLR,
     $     IPLR,JPLR,KPLR,LPLR,CPLR,PDLR,ICPLR)

      CALL CLSFIL(IUNIT)

      PRINT *,'READ ',NDLR,' DIHEDRAL RESTRAINTS'
      PRINT *

 500  FORMAT (1X,A6,': reading dihedral restraint information')
C end GTDHRS
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
COMMSUBR GTLESP
C     SUBROUTINE GTLESP(NRP)
C
C     GTLESP opens a file, reads and print the title block on the file,
C     and then reads in the local elevation
C     dihedral specification data
C     L<NDLE>, L<IPLE>, L<JPLE>, L<KPLE>, L<LPLE>
C     which are defined in coordar.h using subroutine L<RDLESP>.
C     The file is then closed.
C
C     NRP is not modified, but used to check the contents of the
C     arrays L<IPLE>, L<JPLE>, L<KPLE>, L<LPLE> on reading.
COMMEND
C     ==--------------------------------------------------------------==
      SUBROUTINE GTLESP(NRP)
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
C args
      INTEGER NRP
C local params
      INTEGER MAXTIT
      PARAMETER (MAXTIT = 10)
C local vars
      INTEGER IUNIT,IDBLK,NREC,NTIT
      LOGICAL LFORM
      CHARACTER PRGSTR*(6)
      CHARACTER TITLE*(MAXDLN)
      DIMENSION TITLE(MAXTIT)
C data
      DATA PRGSTR /'GTLESP'/
C begin
      CALL DETOPN('LESPEC',IUNIT,LFORM)
      IF (IUNIT .LT. 0) THEN
         PRINT 500,PRGSTR
         PRINT *,'failed to open local elevation specification file'
         CALL STOPGM('GTLESP','FATAL ERROR')
      ENDIF
      IF (.NOT. LFORM) THEN
         PRINT 500,PRGSTR
         PRINT *,'file must be formatted'
         CALL STOPGM('GTLESP','FATAL ERROR')
      ENDIF
      CALL INITLN(IUNIT)

C     get title and print it
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE. ITITID ) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMEXP,PRGSTR,NAMEID(ITITID),IUNIT
         CALL STOPGM('GTLESP','FATAL ERROR')
      ENDIF
      CALL RDTIT(IUNIT,LFORM,NREC,MAXTIT,TITLE,NTIT)
      CALL PRTIT('LOCAL ELEVATION SPECIFICATIONS TITLE',
     $     NTIT,TITLE)

C     get local elevation data
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE.IXLEID ) THEN
         CALL FLAGLN(PRGSTR)
         PRINT 500,PRGSTR
         PRINT FMEXP,PRGSTR,NAMEID(IXLEID),IUNIT
         CALL STOPGM('GTLESP','FATAL ERROR')
      ENDIF

      CALL RDLESP(IUNIT,NRP,MAXDLE,NDLE,
     $     IPLE,JPLE,KPLE,LPLE)

      CALL CLSFIL(IUNIT)

C calculate here the number of INTS we need to store
C one configuration. See coordsz.h and coordar.h.
      NLECFG =  (NDLE + DIPINT -1)/DIPINT
C this error will only happen if the constants in
C coordsz.h have been misconfigured.
      IF (NLECFG .GT. MXLECF) THEN
         PRINT *,PRGSTR,': Too many LE dihedral angles specified!'
         PRINT *,'NLECFG =',NLECFG
         PRINT *,'MXLECF =',MXLECF
         CALL STOPGM('GTLESP','FATAL ERROR')
      ENDIF

      PRINT *,'READ ',NDLE,
     $     ' LOCAL ELEVATION DIHEDRAL SPECIFICATIONS'
      PRINT *

 500  FORMAT (1X,A6,': reading local elevation dihedral',
     $     ' spec. information')
C end GTLESP
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GTJVSP(NRP)
C     ==--------------------------------------------------------------==
C     SUBROUTINE GTJVSP(NRP)
C
C     GTJVSP opens a file, reads and prints a title block
C     and then reads in the j-value restraining
C     dihedral specification data into the arrays
C     L<NDJV>,
C     L<IPJV>, L<JPJV>, L<KPJV>, L<LPJV>,
C     L<CPJV>, L<PJR0>, L<PSJR>, L<AJV>, L<BJV>, L<CJV>
C     which are defined in coordar.h using subroutine L<RDJVSP>.
C     The file is then closed.
C     
C     NRP is not modified, but used to check the correctness
C     of the L<IPJV>, L<JPJV>, L<KPJV> and L<LPJV> arrays when
C     reading.
C
C     The L<PSJR> array is converted from degrees to radians.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
C args
      INTEGER NRP
C local params
      INTEGER MAXTIT
      PARAMETER (MAXTIT = 10)
      REAL*8 ONE
      PARAMETER (ONE = 1.0d0)
C local vars
      INTEGER IUNIT,IDBLK,NREC,NTIT,I
      LOGICAL LFORM
      REAL*8 CONV
      CHARACTER PRGSTR*(6)
      CHARACTER TITLE*(MAXDLN)
      DIMENSION TITLE(MAXTIT)
C data
      DATA PRGSTR /'GTJVSP'/
C begin
      CALL DETOPN('JVALSPEC',IUNIT,LFORM)
      IF (IUNIT .LT. 0) THEN
         PRINT 500,PRGSTR
         PRINT *,'failed to open j-val restraining specification file'
         CALL STOPGM('GTJVSP','FATAL ERROR')
      ENDIF
      IF (.NOT. LFORM) THEN
         PRINT 500,PRGSTR
         PRINT *,'file must be formatted'
         CALL STOPGM('GTJVSP','FATAL ERROR')
      ENDIF
      CALL INITLN(IUNIT)

C     get title and print it
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE. ITITID ) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMEXP,PRGSTR,NAMEID(ITITID),IUNIT
         CALL STOPGM('GTJVSP','FATAL ERROR')
      ENDIF
      CALL RDTIT(IUNIT,LFORM,NREC,MAXTIT,TITLE,NTIT)
      CALL PRTIT('J-VALUE RESTRAINING SPECIFICATIONS TITLE',
     $     NTIT,TITLE)

C     get local elevation data
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE.IXJVID) THEN
         CALL FLAGLN(PRGSTR)
         PRINT 500,PRGSTR
         PRINT FMEXP,PRGSTR,NAMEID(IXJVID),IUNIT
         CALL STOPGM('GTJVSP','FATAL ERROR')
      ENDIF

      CALL RDJVSP(IUNIT,NRP,MAXJ,NDJV,
     $     IPJV, JPJV, KPJV, LPJV,
     $     CPJV,PJR0,PSJR,AJV,BJV,CJV)

      CALL CLSFIL(IUNIT)

C now convert the PSJR array to radians
      CONV=45.0d0/DATAN(ONE)
!$OMP parallel do private(I)
      DO I=1,NDJV
         PSJR(I) = PSJR(I)/CONV
      ENDDO

      PRINT *,'READ ',NDJV,
     $     ' J-VALUE RESTRAINING SPECIFICATIONS'
      PRINT *

 500  FORMAT (1X,A6,': reading j-value restraining',
     $     ' spec. information')
C end GTJVSP
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE RDGAM(NATTOT,GAM)
C     ==--------------------------------------------------------------==
C     SUBROUTINE RDGAM(NATTOT,GAM)
C
C     Opens a file, reads and prints the title block, and then
C     reads in atomic friction coefficient values into
C     the array GAM using L<RDGARR>. The file is then closed.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'dataid.h'
      INCLUDE 'fileio.h'
      INCLUDE 'coordsz.h'
C args
      INTEGER NATTOT
      real*8 GAM(NATTOT)
C local params
      INTEGER MAXTIT
      PARAMETER (MAXTIT = 10)
C local vars
      INTEGER IUNIT,IDBLK,NREC,NTIT
      LOGICAL LFORM
      CHARACTER PRGSTR*(6)
      CHARACTER TITLE*(MAXDLN)
      DIMENSION TITLE(MAXTIT)
C data
      DATA PRGSTR / 'RDGAM'/
C begin
      CALL DETOPN('GAMMA',IUNIT,LFORM)
      IF (IUNIT .LT. 0) THEN
         PRINT 500,PRGSTR
         PRINT *,'failed to open GAMMA data file'
         CALL STOPGM('RDGAM','FATAL ERROR')
      ENDIF
      IF (.NOT. LFORM) THEN
         PRINT 500,PRGSTR
         PRINT *,'file must be formatted'
         CALL STOPGM('RDGAM','FATAL ERROR')
      ENDIF
      CALL INITLN(IUNIT)

C     get title and print it
      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE. ITITID ) THEN
         CALL FLAGLN(PRGSTR)
         PRINT FMEXP,PRGSTR,NAMEID(ITITID),IUNIT
         CALL STOPGM('RDGAM','FATAL ERROR')
      ENDIF
      CALL RDTIT(IUNIT,LFORM,NREC,MAXTIT,TITLE,NTIT)
      CALL PRTIT('ATOMIC FRICTION COEFFICIENTS TITLE',NTIT,TITLE)

      CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
      IF (IDBLK .NE. IXFRID) THEN
         CALL FLAGLN(PRGSTR)
         PRINT 500,PRGSTR
         PRINT FMEXP,PRGSTR,NAMEID(IXFRID),IUNIT
         CALL STOPGM('RDGAM','FATAL ERROR')
      ENDIF
C     now get the GAMMA array
      CALL RDGARR(IUNIT,LFORM,NREC,NATTOT,GAM)

      CALL CLSFIL(IUNIT)

 500  FORMAT (1X,A6,': reading atomic friction coefficients')
C end rdgam
C     ==--------------------------------------------------------------==
      RETURN
      END
