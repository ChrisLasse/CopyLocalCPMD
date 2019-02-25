C     ==================================================================
      SUBROUTINE SCOCN2(LFIELD,MAX,STRA,STRB,STRRES)
C     ==--------------------------------------------------------------==
C string.f,v 1.6 1996/02/12 13:05:33 wscott Exp
Cmb - Revised on 15 Sept. 2005
C some string operations used in GROMOS
C
COMMSUBR SCOCN2
C     SUBROUTINE SCOCN2(LFIELD,MAX,STRA,STRB,STRRES)
C
C     SCOCN2 concatenates the the first LFIELD
C     characters of the strings STRA and STRB with
C     a hyphen ('-' character) between them.
C     The resulting string is written into STRRES.
C
C     LFIELD is the width of the field in num of chars.
C     MAX is the length of the strings A,B.
C     ==--------------------------------------------------------------==
C args
      INTEGER LFIELD,MAX
      CHARACTER STRA*(*)
      CHARACTER STRB*(*)
      CHARACTER STRRES*(*)

C local vars
      INTEGER LEA
C begin
      LEA = INDEX(STRA,' ')
C determine length of string A
      IF (LEA .EQ. 0) THEN
         LEA = MAX
      ELSE
         LEA = LEA - 1
      ENDIF
      STRRES = ' '
      IF (LEN(STRRES) .LT. 2*LFIELD+1) THEN
         PRINT *,'SCOCN2: string too short!'
         CALL STOPGM('SCOCN2','FATAL ERROR')
      ENDIF
C concat at the right pos
      STRRES(LFIELD-LEA+1:LFIELD) = STRA 
      STRRES(LFIELD+1:LFIELD+1) = '-'
      STRRES(LFIELD+2:2*LFIELD+1) = STRB
C end
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SCOCN3(LFIELD,MAX,STRA,STRB,STRC,STRRES)
C     ==--------------------------------------------------------------==
C LFIELD is the width of the field in num of chars 
C MAX is the length of the strings A,B,C
C args
      INTEGER LFIELD,MAX
      CHARACTER STRA*(*)
      CHARACTER STRB*(*)
      CHARACTER STRC*(*)
      CHARACTER STRRES*(*)

C local vars
      INTEGER MXNLE3
      INTEGER TOTLE3
      INTEGER LEA,LEB
      INTEGER LAF,LBF,LCF,LD1,LD2
C begin
C determine length of string A and B
      MXNLE3 = 3*MAX+2
      TOTLE3 = 3*LFIELD+2

      LEA = INDEX(STRA,' ')
      IF (LEA .EQ. 0) THEN
         LEA = MAX
      ELSE
         LEA = LEA - 1
      ENDIF

      LEB = INDEX(STRB,' ')
      IF (LEB .EQ. 0) THEN
         LEB = MAX
      ELSE
         LEB = LEB - 1
      ENDIF

      STRRES = ' '
      IF (LEN(STRRES) .LT. MXNLE3) THEN
         PRINT *,'SCOCN3: string too short!'
         CALL STOPGM('SCOCN3','FATAL ERROR')
      ENDIF

      LBF = (TOTLE3 - LEB)/2
      LD2 = LBF+LEB
      LCF  = LD2+1

      LD1 = LBF-1
      LAF  = LD1-LEA

      STRRES(LAF:LAF+LEA) = STRA
      STRRES(LD1:LD1)   = '-'
      STRRES(LBF:LBF+LEB) = STRB
      STRRES(LD2:LD2) = '-'
      STRRES(LCF:MXNLE3) = STRC
C end
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SCOCN4(LFIELD,MAX,STRA,STRB,STRC,STRD,STRRES)
C     ==--------------------------------------------------------------==
C args
      INTEGER LFIELD,MAX
      CHARACTER STRA*(*)
      CHARACTER STRB*(*)
      CHARACTER STRC*(*)
      CHARACTER STRD*(*)
      
      CHARACTER STRRES*(*)

C local vars
      INTEGER MXNLE4
      INTEGER TOTLE4
      INTEGER LEA,LEB,LEC
      INTEGER LAF,LBF,LCF,LDF
      INTEGER LD1,LD2,LD3
C begin
C determine length of string A and B
      TOTLE4 = 4*LFIELD+3
      MXNLE4 = 4*MAX+3
      IF (LEN(STRRES) .LT. MXNLE4) THEN
         PRINT *,'SCOCN4: string too short!'
         CALL STOPGM('SCOCN4','FATAL ERROR')
      ENDIF

      LEA = INDEX(STRA,' ')
      IF (LEA .EQ. 0) THEN
         LEA = MAX
      ELSE
         LEA = LEA - 1
      ENDIF

      LEB = INDEX(STRB,' ')
      IF (LEB .EQ. 0) THEN
         LEB = MAX
      ELSE
         LEB = LEB - 1
      ENDIF

      LEC = INDEX(STRC,' ')
      IF (LEC .EQ. 0) THEN
         LEC = MAX
      ELSE
         LEC = LEC - 1
      ENDIF

      STRRES = ' '

      LD2 = (2*LFIELD +2)
      LBF = LD2-LEB
      LD1 = LBF-1
      LAF = LD1-LEA

      LCF = LD2+1
      LD3 = LCF+LEC
      LDF = LD3+1

      STRRES(LAF:LAF+LEA) = STRA
      STRRES(LD1:LD1)   = '-'
      STRRES(LBF:LBF+LEB) = STRB
      STRRES(LD2:LD2) = '-'
      STRRES(LCF:LCF+LEC) = STRC
      STRRES(LD3:LD3) = '-'
      STRRES(LDF:MXNLE4) = STRD
C end scocn4
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SCOLFT(I,STRA)
C     ==--------------------------------------------------------------==
C format the string with left adjusted integer
C args
      INTEGER I
      CHARACTER STRA*(*)
C local params
      INTEGER MXFMT
      PARAMETER (MXFMT = 16)
C loca vars
      CHARACTER FSTR*(MXFMT)
      INTEGER IRES,IPOS,IBAR,ICHECK
C begin
      IF (I .GE. 0) THEN
         IPOS = 1
         ICHECK = I
      ELSE
         IPOS = 2
         ICHECK = -I
      ENDIF

      IBAR = 10
C while ICHECK >= IBAR
 10   IF (ICHECK .GE. IBAR) THEN
         IPOS = IPOS + 1
         IBAR = IBAR*10
C end while
         GOTO 10
      ENDIF
C
C AK: 2005/07/22.
C all Gromos files written through CPMD are not vital
C for the further execution of the program, so we just
C write a marker and continue. Gromos STOPs here.
      WRITE (FSTR,IOSTAT = IRES,FMT = '(''(I'',I1,'')'')')IPOS
      IF (IRES .NE. 0) THEN
         WRITE(6,*) ' SCOLFT| BLA ERROR!'
         WRITE(STRA,'(A)') '-#X#-'
         RETURN
      ENDIF
      WRITE (STRA,IOSTAT = IRES, FMT = FSTR)I
      IF (IRES .NE. 0)THEN
         WRITE(6,*) ' SCOLFT| FORMAT ERROR!'
         WRITE(STRA,'(A)') '-#%#-'
      ENDIF
C end SCOLFT
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SCOIN2(I,J,MAX,STRRES)
C     ==--------------------------------------------------------------==
C args
      INTEGER I,J,MAX
      CHARACTER STRRES*(*)
C params
      INTEGER MAXLOC
      PARAMETER (MAXLOC = 5)
C local vars
      CHARACTER STRA*(MAXLOC)
      CHARACTER STRB*(MAXLOC)
C begin
      CALL SCOLFT(I,STRA)
      CALL SCOLFT(J,STRB)
      CALL SCOCN2(MAX,MAXLOC,STRA,STRB,STRRES)
C end scoin3
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SCOIN3(I,J,K,STRRES)
C     ==--------------------------------------------------------------==
C args
      INTEGER I,J,K
      CHARACTER STRRES*(*)
C params
      INTEGER MAXLOC
      PARAMETER (MAXLOC = 5)
C local vars
      CHARACTER STRA*(MAXLOC)
      CHARACTER STRB*(MAXLOC)
      CHARACTER STRC*(MAXLOC)
C begin
      CALL SCOLFT(I,STRA)
      CALL SCOLFT(J,STRB)
      CALL SCOLFT(K,STRC)
      CALL SCOCN3(4,MAXLOC,STRA,STRB,STRC,STRRES)
C end scoin3
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SCOIN4(I,J,K,L,STRRES)
C     ==--------------------------------------------------------------==
C args
      INTEGER I,J,K,L
      CHARACTER STRRES*(*)
C local param
      INTEGER MAXLOC
      PARAMETER (MAXLOC = 5)
C local vars
      CHARACTER STRA*(MAXLOC)
      CHARACTER STRB*(MAXLOC)
      CHARACTER STRC*(MAXLOC)
      CHARACTER STRD*(MAXLOC)
C begin
      CALL SCOLFT(I,STRA)
      CALL SCOLFT(J,STRB)
      CALL SCOLFT(K,STRC)
      CALL SCOLFT(L,STRD)
      CALL SCOCN4(4,MAXLOC,STRA,STRB,STRC,STRD,STRRES)
C end scoin4
      END
C     ==--------------------------------------------------------------==
