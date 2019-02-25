C hello emacs, this is -*- fortran -*-
      SUBROUTINE DMPRED2(IUNIT,LFORM,PRGSTR,IDBLK,NAT,NDIM,X,FRMT)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,IDBLK,NAT,NDIM
      LOGICAL LFORM
      CHARACTER PRGSTR*(*)
      real*8 X(NDIM*NAT)
C Phil
      CHARACTER*(*) FRMT
C '(3F15.9)' -> idem as DMPRED

C local vars
      INTEGER NREC,IRES,I,I3
C begin
      IF (LFORM) THEN
         CALL WRSTR(IUNIT,PRGSTR,NAMEID(IDBLK))
         I3 = 1
         DO 10 I=1,NAT
            WRITE (UNIT=IUNIT, IOSTAT=IRES,FMT=FRMT)
     $           X(I3),X(I3+1),X(I3+2)
            IF (IRES .NE. 0) THEN
               PRINT '(2A)',PRGSTR,' Error writing data!'
               STOP
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
            STOP
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES) NAT
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'NAT',IUNIT
            STOP
         ENDIF
C
         WRITE(UNIT = IUNIT, IOSTAT = IRES)
     $        ((X(NDIM*I+I3),I3=1,3),I=0,NAT-1)
         IF (IRES .NE. 0)THEN
            PRINT FWFAIL,PRGSTR,'coordinate data',IUNIT
            STOP
         ENDIF
      ENDIF
C     end DMPRED2
      END



COMMSUBR WRFRED
C     SUBROUTINE WRFRED(IUNIT,LFORM,LWR4,NAT,NDIM,X)
C
C     WRVRED writes forces to a file.
COMMEND

      SUBROUTINE WRFRED(IUNIT,LFORM,LWR4,NAT,NDIM,X)
C includes
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NAT,NDIM
      LOGICAL LFORM,LWR4
      real*8 X(NDIM*NAT)
C begin
      CALL DMPRED2(IUNIT,LFORM,'WRFRED',IFRDID,NAT,NDIM,X,'(3F20.11)')
      IF (LWR4 .AND. NDIM .GT. 3) THEN
         CALL DMP4RE(IUNIT,LFORM,'WRFRED',IF4RID,NAT,NDIM,X)
      ENDIF
C end WRFRED
      END

COMMSUBR RDVRED
C     SUBROUTINE RDFRED(IUNIT,LFORM,LDO4,NREC,NAT,NDIM,X)
C
C     from an open file.
COMMEND

      SUBROUTINE RDFRED(IUNIT,LFORM,LDO4,NREC,NAT,NDIM,X)
      INCLUDE 'dataid.h'
C args
      INTEGER IUNIT,NAT,NDIM,NREC
      LOGICAL LFORM,LDO4
      real*8 X(NDIM*NAT)
C local vars
      INTEGER IDBLK
      CHARACTER PRGSTR*(6)
      DATA PRGSTR /'RDFRED'/
C begin
      CALL GETRED(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
      IF (LDO4 .AND. NDIM .GT. 3) THEN
         CALL RDBHDR(IUNIT,LFORM,IDBLK,NREC)
         IF (IDBLK .EQ. IF4RID) THEN
            CALL GET4RE(IUNIT,LFORM,NREC,PRGSTR,NAT,NDIM,X)
         ELSE
            IF (LFORM) THEN
               CALL FLAGLN(PRGSTR)
            ENDIF
            PRINT *,NAMEID(IF4RID), 'expected!'
            STOP
         ENDIF
      ENDIF
C end RDFRED
      END

