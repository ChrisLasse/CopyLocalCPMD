C     ==================================================================
      SUBROUTINE WRTOPO(IUNIT)
C     ==--------------------------------------------------------------==
C     SUBROUTINE WRTOPO(IUNIT)
C
C     WRTOPO writes a topology to an open formatted file,
C     which it does not close afterwards.
C
C     Thus, typical usage would be:
C     CALL OPNFIL('TOPO','FORMATTED','NEW',IUNIT)
C     IF (IUNIT .LT.0) THEN
C        PRINT *,PRGSTR,': failed to open GROMOS95 topology file!'
C        STOP
C     ENDIF
C     CALL WRTOPO(IUNIT)
C     CALL CLSFIL(IUNIT)
C
C     The topology is assumed to be stored in the common blocks
C     defined in toposz.h> and topoar.h.
C     Some simple checking of consistency is made while
C     writing to file.If an error occurs, a message is printed
C     and STOP is called.
C
C     NOTE: this routine does NOT convert any of the variables
C     on writing, but writes the arrays as they are found.
C     Thus, for example, the calling routine must ensure that T0 is
C     in degrees, and not a COS of the angles as it is used in PROMD.
C     This can be accomplished by calling PWRTOP before writing
C     if necessary.
Cmb - Revised on 16 Sept. 2005
C     ==--------------------------------------------------------------==
C include the topology definition
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
      INCLUDE 'topblock.h'
C the maximum number of excluded atoms to write on a line after the atom
C information.
      INTEGER MAXEXC
      PARAMETER (MAXEXC = 4)
      INTEGER MAXEXL
      PARAMETER (MAXEXL = 12)
C args
      INTEGER IUNIT
C local vars
      REAL*8 PYE,RD2DEG,CGDIV,RDSGSQ
      INTEGER IRES,I,J,K,ILIM,ITMP,IATCG
      INTEGER KMIN,KMAX,KWRT
      LOGICAL LERR,LOK
      CHARACTER PRGSTR*(6)
      CHARACTER IACSTR*(40)
      CHARACTER TSTR*(76)
      DATA IACSTR /'integer (van der Waals) atom type code'/
      DATA PRGSTR /'WRTOPO'/
C begin
      PYE = 4.0D0*DATAN(1.d0)
      RD2DEG = 180.0D0/PYE
      RDSGSQ = RD2DEG*RD2DEG 
      LOK  = .TRUE.
C write a little header
      CALL WRCOMM(IUNIT,PRGSTR,'GROMOS TOPOLOGY FILE')
      CALL WRCOMM(IUNIT,PRGSTR,'WRTOPO version:')
      CALL WRCOMM(IUNIT,PRGSTR,
     .   'Id: wrtopo.f,v 1.19 1996/10/18 14:49:29 wscott Exp')
      CALL WRCOMM(IUNIT,PRGSTR,' ')
C title block
      CALL WRTIT(IUNIT,.TRUE.,NTPLNS,TOPTIT)

C units
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPUNT))
      CALL WRCOMM(IUNIT,PRGSTR,
     .'FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)')
      IF (FPEPSI .LE. 0.0D0)THEN
         PRINT *,PRGSTR,': invalid value for FPEPSI =',FPEPSI
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF
      WRITE (UNIT=IUNIT, IOSTAT=IRES, FMT=12) FPEPSI
      IF (IRES .NE. 0) THEN
         PRINT *,PRGSTR,': Error writing FPEPSI!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF
      CALL WRCOMM(IUNIT,PRGSTR,
     .     'HBAR: Planck''s constant HBAR = H/(2* PI)')
      WRITE (UNIT = IUNIT, IOSTAT = IRES, FMT = 12) HBAR
      IF (IRES .NE. 0) THEN
         PRINT *,PRGSTR,': Error writing HBAR!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C version
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPVER))
      WRITE (UNIT = IUNIT, IOSTAT = IRES, FMT = 13) TPVER
      IF (IRES .NE. 0) THEN
         PRINT *,PRGSTR,': Error writing TPVER!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C natname
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NATNAM))
C nratt
      CALL WRCOMM(IUNIT,PRGSTR,
     .     'NRATT: number of van der Waals atom types')
      WRITE (UNIT = IUNIT, IOSTAT = IRES, FMT = 14) NRATT
      IF (IRES .NE. 0) THEN
         PRINT *,PRGSTR,': Error writing NRATT!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF
      IF (NRATT .LT. 0 .OR. NRATT .GT. MAXATT) THEN
         PRINT *,PRGSTR,': NRATT out of range!'
         PRINT *,'WRITE ',NRATT,'MAXATT=',MAXATT
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C type
      CALL WRCOMM(IUNIT,PRGSTR,'FFTYPE: atom type names')
      DO 5 K=1,NRATT
         WRITE (IUNIT, IOSTAT = IRES, FMT = 15)FFTYPE(K)
         IF (IRES .NE. 0)THEN
            PRINT *,PRGSTR,': Error writing FFTYPE!'
            CALL STOPGM('WRTOPO','FATAL ERROR')
         ENDIF
C
         IF (MOD(K,10) .EQ. 0)THEN
            CALL WRCINT(IUNIT,PRGSTR,K)
         ENDIF
 5    CONTINUE
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C nraa2
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NRSNAM))
      IF (NRAA2 .LT. 0 .OR. NRAA2 .GT. MAXAA2) THEN
         PRINT *,PRGSTR,': NRAA2 out of range!'
         PRINT *,'WRITE ',NRAA2,'MAXAA2=',MAXAA2
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

      CALL WRCOMM(IUNIT,PRGSTR,
     .     'NRAA2: number of residues in a solute molecule')
      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NRAA2
      IF (IRES .NE. 0)THEN
         PRINT *,PRGSTR,': Error writing NRAA2!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C aanm
      CALL WRCOMM(IUNIT,PRGSTR,'AANM: residue names')
      DO 7 K=1,NRAA2
         WRITE (IUNIT, IOSTAT = IRES, FMT = 15) AANM(K)
         IF (IRES .NE. 0) THEN
            PRINT *,PRGSTR,': Error writing AANM!'
            CALL STOPGM('WRTOPO','FATAL ERROR')
         ENDIF
C
         IF (MOD(K,10) .EQ. 0)THEN
            CALL WRCINT(IUNIT,PRGSTR,K)
         ENDIF
 7    CONTINUE
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C solute atoms
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NSOLAT))
C nrp
      IF (NRP .LT. 0 .OR. NRP .GT. MAXNRP) THEN
         PRINT *,PRGSTR,': NRP out of range!'
         PRINT *,'WRITE ',NRP,'MAXNRP=',MAXNRP
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

      CALL WRCOMM(IUNIT,PRGSTR,'  NRP: number of solute atoms')
      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NRP
      IF (IRES .NE. 0)THEN
         PRINT *,PRGSTR,': Error writing NRP!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C winv,panm,mres,iac
      CALL WRCOMM(IUNIT,PRGSTR,' ATNM: atom number')
      CALL WRCOMM(IUNIT,PRGSTR,' MRES: residue number')
      CALL WRCOMM(IUNIT,PRGSTR,' PANM: atom name of solute atom')
      TSTR = '  IAC: '// IACSTR
      CALL WRCOMM(IUNIT,PRGSTR,TSTR)
      CALL WRCOMM(IUNIT,PRGSTR,' MASS: mass of solute atom')
      CALL WRCOMM(IUNIT,PRGSTR,'   CG: charge of solute atom')
      CALL WRCOMM(IUNIT,PRGSTR,'  CGC: charge group code (0 or 1)')
      CALL WRCOMM(IUNIT,PRGSTR,'  INE: number of excluded atoms')
      CALL WRCOMM(IUNIT,PRGSTR,'INE14: number of 1-4 interactions')

      WRITE (IUNIT, IOSTAT = IRES, FMT = '(2A5,A6,A4,2A11,2A4)')
     .   '#ATNM','MRES','PANM','IAC','MASS','CG','CGC','INE'

      WRITE (IUNIT, IOSTAT = IRES, FMT = '(A1,46X,A5)')'#','INE14'
C ----this here as a reminder
C EXCLUSIONS J OF ATOM I ARE POSITIONED AT    
C                     JSNE(KNE(I)+1),...,JSNE(KNE(I)+INE(I));
C                     IT IS ASSUMED THAT I.LE.J
C ----

      LERR = .FALSE.
      DO 9 K=1,NRP
C determine whether K terminates a cg
         ITMP = 1
C while
 111     IF (ITMP .LE. NCAG .AND. INC(ITMP) .LT. K) THEN
            ITMP = ITMP + 1
            GOTO 111
         ENDIF
C end while
         IF (INC(ITMP) .EQ. K) THEN
            IATCG = 1
         ELSE
            IATCG = 0
         ENDIF

C in order to avoid the line getting to long, we only write
C MAXEXC exclusions on this line...
C if there are more for this atom, i.e. INE > MAXEXC, then
C write them on the next line.

         KMIN = KNE(K)+1
         KWRT = INE(K)
         IF (KWRT .GT. MAXEXC) THEN
            KWRT = MAXEXC
         ENDIF
         KMAX = KNE(K)+KWRT

         CGDIV = DSQRT(FPEPSI)
         WRITE (IUNIT, IOSTAT = IRES,FMT = 17)
!    .        K,MRES(K),PANM(K),IAC(K),WMAS(K),CG(K),IATCG,
     .    K,MRES(K),PANM(K),IAC(K),WMAS(K),CG(K)/CGDIV,IATCG,  !QMMM
     .    INE(K),(JSNE(ITMP), ITMP = KMIN,KMAX)
         IF (IRES .NE. 0)THEN
            PRINT 165,PRGSTR
            CALL STOPGM('WRTOPO','FATAL ERROR')
         ENDIF
C now any additional exclusions if there are any
C while loop
 125     IF (KWRT .LT. INE(K)) THEN
            KWRT = MIN(INE(K),KWRT+MAXEXL)
            KMIN = KMAX+1
            KMAX = KNE(K)+KWRT
            WRITE (IUNIT,IOSTAT = IRES, FMT =14)
     .           (JSNE(ITMP), ITMP = KMIN,KMAX)
            IF (IRES .NE. 0)THEN
               PRINT 165,PRGSTR
               CALL STOPGM('WRTOPO','FATAL ERROR')
            ENDIF
C end of while loop
            GOTO 125
         ENDIF

C now 1-4 interactions
         KMIN = KNE14(K)+1
         KWRT = INE14(K)
         IF (KWRT .GT. MAXEXC) THEN
            KWRT = MAXEXC
         ENDIF
         KMAX = KNE14(K)+KWRT

         WRITE (IUNIT, IOSTAT = IRES,FMT = 18)
     .        INE14(K),(JSNE14(ITMP), ITMP = KMIN,KMAX)
         IF (IRES .NE. 0)THEN
            PRINT 165,PRGSTR
            CALL STOPGM('WRTOPO','FATAL ERROR')
         ENDIF

C now any additional 1-4s if there are any
C while loop
 127     IF (KWRT .LT. INE14(K)) THEN
            KWRT = MIN(INE14(K),KWRT+MAXEXL)
            KMIN = KMAX+1
            KMAX = KNE14(K)+KWRT
            WRITE (IUNIT,IOSTAT = IRES, FMT =14)
     .           (JSNE14(ITMP), ITMP = KMIN,KMAX)
            IF (IRES .NE. 0)THEN
               PRINT 165,PRGSTR
               CALL STOPGM('WRTOPO','FATAL ERROR')
            ENDIF
C end of while loop
            GOTO 127
         ENDIF
 9    CONTINUE

      IF (LERR) THEN
         LOK = .FALSE.
      ENDIF
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C nbty
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NBNDTY))
      IF (NBTY .LT. 0 .OR. NBTY .GT. MAXNBT) THEN
         PRINT *,PRGSTR,': NBTY out of range!'
         PRINT *,'WRITE ',NBTY,'MAXNBT=',MAXNBT
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' NBTY: number of covalent bond types')
      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NBTY
      IF (IRES .NE. 0)THEN
         PRINT *,PRGSTR,': Error writing NBTY!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C cb, b0
      CALL WRCOMM(IUNIT,PRGSTR,' CB: force constant')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' B0: bond length at minimum energy')
      WRITE (IUNIT, IOSTAT = IRES, FMT = '(A1,A14,A15)')
     .   '#','CB','B0'

      DO 150 K=1,NBTY
         WRITE (IUNIT, IOSTAT = IRES, FMT = 19)CB(K),B0(K)
         IF (IRES .NE. 0) THEN
            PRINT *,PRGSTR,': Error writing CB, B0!'
            CALL STOPGM('WRTOPO','FATAL ERROR')
         ENDIF
C
         IF (MOD(K,10) .EQ. 0)THEN
            CALL WRCINT(IUNIT,PRGSTR,K)
         ENDIF
 150  CONTINUE
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C nbonh
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NBNDHY))
      IF (NBONH .LT. 0 .OR. NBONH .GT. MAXBNH) THEN
         PRINT *,PRGSTR,': NBONH out of range!'
         PRINT *,'WRITE ',NBONH,'MAXBNH=',MAXBNH
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF
C
      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' NBONH: number of bonds involving H atoms in solute')
      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NBONH
      IF (IRES .NE. 0) THEN
         PRINT *,PRGSTR,': Error writing NBONH!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C ibh,jbh,icbh
      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' IBH, JBH: atom sequence numbers of atoms forming a bond')
      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' ICBH: bond type code')
      WRITE (IUNIT, IOSTAT = IRES, FMT = '(A1,A4,2A5)')
     .   '#','IBH','JBH','ICBH'

      CALL WBND(IUNIT,PRGSTR,NBONH,NRP,NBTY,IBH,JBH,ICBH,LOK)
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C nbon
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NBND))
      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' NBON: number of bonds NOT involving H atoms in solute')
      IF (NBON .LT. 0 .OR. NBON .GT. MAXBON)THEN
         PRINT *,PRGSTR,':NBON out of range!'
         PRINT *,'WRITE ',NBON,'MAXBON=',MAXBON
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF
      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NBON
      IF (IRES .NE. 0)THEN
         PRINT *,PRGSTR,': Error writing NBON'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF


C ib,jb,icb
      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' IB, JB: atom sequence numbers of atoms forming a bond')
      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' ICB: bond type code')
      WRITE (IUNIT, IOSTAT = IRES, FMT = '(A1,A4,2A5)')
     .   '#','IB','JB','ICB'

      CALL WBND(IUNIT,PRGSTR,NBON,NRP,NBTY,IB,JB,ICB,LOK)
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C bond angle types
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NBATY))
C ntty
      IF ( NTTY .LT. 0 .OR. NTTY .GT. MAXTTY) THEN
         PRINT *,PRGSTR,': NTTY out of range!'
         PRINT *,'WRITE ',NTTY,'MAXTTY=',MAXTTY
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' NTTY: number of bond angle types')
      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NTTY
      IF (IRES .NE. 0)THEN
         PRINT *,PRGSTR,': Error writing NTTY!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C ct,t0
      CALL WRCOMM(IUNIT,PRGSTR,' CT: force constant')
      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' T0: bond angle at minimum energy in degrees')
      WRITE (IUNIT, IOSTAT = IRES, FMT = '(A1,A14,A15)')
     .   '#','CT','T0'
      DO K=1,NTTY
!        WRITE (IUNIT, IOSTAT=IRES, FMT=19)CT(K),T0(K)
         WRITE (IUNIT, IOSTAT=IRES, FMT=19)CT(K),RD2DEG*DACOS(T0(K)) !QMMM
         IF (IRES .NE. 0) THEN
            PRINT *,PRGSTR,': Error writing CT,T0!'
            CALL STOPGM('WRTOPO','FATAL ERROR')
         ENDIF
C
         IF (MOD(K,10) .EQ. 0)THEN
            CALL WRCINT(IUNIT,PRGSTR,K)
         ENDIF
      ENDDO
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C ntheh
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NBAHY))
      IF (NTHEH .LT. 0 .OR. NTHEH .GT. MXQHEH) THEN
         PRINT *,PRGSTR,': NTHEH out of range!'
         PRINT *,'WRITE ',NTHEH,'MXQHEH=',MXQHEH
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' NTHEH: number of bond angles involving H atoms in solute')

      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NTHEH
      IF (IRES .NE. 0)THEN
         PRINT *,PRGSTR,': Error writing NTHEH'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C ith,jth,kth,icth
      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' ITH, JTH, KTH: atom sequence numbers ')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     '   of atoms forming a bond angle in solute')
      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' ICTH: bond angle type code')
      WRITE (IUNIT, IOSTAT = IRES, FMT = '(A1,A4,3A5)')
     .   '#','ITH','JTH','KTH','ICTH'
C
      CALL WRBAN(IUNIT,PRGSTR,NTHEH,NRP,NTTY,ITH,JTH,KTH,ICTH,LOK)
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C nthe
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NBAN))
      IF (NTHE .LT. 0 .OR. NTHE .GT. MAXTHE) THEN
         PRINT *,PRGSTR,': NTHE out of range!'
         PRINT *,'WRITE ',NTHE,'MAXTHE=',MAXTHE
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' NTHE: number of bond angles NOT ')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     '   involving H atoms in solute')
      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NTHE
      IF (IRES .NE. 0) THEN
         PRINT *,PRGSTR,': Error writing NTHE'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C it,jt,kt,ict
      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' IT, JT, KT: atom sequence numbers of atoms')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     '    forming a bond angle')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' ICT: bond angle type code')
      WRITE (IUNIT, IOSTAT = IRES, FMT = '(A1,A4,3A5)')
     .     '#','IT','JT','KT','ICT'
C
      CALL WRBAN(IUNIT,PRGSTR,NTHE,NRP,NTTY,IT,JT,KT,ICT,LOK)
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C nqty
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NIMDTY))
      IF (NQTY .LT. 0 .OR. NQTY .GT. MAXQTY) THEN
         PRINT *,PRGSTR,': NQTY out of range!'
         PRINT *,'WRITE ',NQTY,'MAXQTY=',MAXQTY
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' NQTY: number of improper dihedrals')
      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NQTY
      IF (IRES .NE. 0) THEN
         PRINT *,PRGSTR,': Error writing NQTY!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C cq,q0
      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' CQ: force constant of improper dihedral per degrees square')
      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' Q0: improper dihedral angle at minimum energy in degrees')
      WRITE (IUNIT, IOSTAT = IRES, FMT = '(A1,A14,A15)')
     .   '#','CQ','Q0'
C
      DO 52 K=1,NQTY
!        WRITE (IUNIT, IOSTAT=IRES, FMT=19)CQ(K),Q0(K)   
         WRITE (IUNIT, IOSTAT=IRES, FMT=19)CQ(K)/RDSGSQ,Q0(K)*RD2DEG !QMMM
         IF (IRES .NE. 0)THEN
            PRINT *,PRGSTR,': Error writing CQ, Q0!'
            CALL STOPGM('WRTOPO','FATAL ERROR')
         ENDIF
C
         IF (MOD(K,10) .EQ. 0)THEN
            CALL WRCINT(IUNIT,PRGSTR,K)
         ENDIF
 52   CONTINUE
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C nqhih
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NIMDH))
      IF (NQHIH .LT. 0 .OR. NQHIH .GT. MAXHIH) THEN
         PRINT *,PRGSTR,': NQHIH out of range!'
         PRINT *,'WRITE ',NQHIH,'MAXHIH=',MAXHIH
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF
      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' NQHIH: number of improper dihedrals')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     '   involving H atoms in the solute')

      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NQHIH
      IF (IRES .NE. 0)THEN
         PRINT *,PRGSTR,': Error writing NQHIH!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C iqh,jqh,kqh,lqh,icqh
      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' IQH,JQH,KQH,LQH: atom sequence numbers')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     '    of atoms forming an improper dihedral')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' ICQH: improper dihedral type code')
      WRITE (IUNIT, IOSTAT = IRES, FMT = '(A1,A4,4A5)')
     .     '#','IQH','JQH','KQH','LQH','ICQH'

      CALL WRQ(IUNIT,PRGSTR,NQHIH,NRP,NQTY,
     .     IQH,JQH,KQH,LQH,ICQH,LOK)
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))


C nqhi
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NIMD))
      IF (NQHI .LT. 0 .OR. NQHI .GT. MAXQHI) THEN
         PRINT *,PRGSTR,': NQHI out of range!'
         PRINT *,'WRITE ',NQHI,'MAXQHI=',MAXQHI
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF
      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' NQHI: number of improper dihedrals NOT')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     '   involving H atoms in solute')
      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NQHI
      IF (IRES .NE. 0)THEN
         PRINT *,PRGSTR,': Error writing NQHI!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C iq,jq,kq,lq,icq
      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' IQ,JQ,KQ,LQ: atom sequence numbers of atoms')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     '   forming an improper dihedral')
      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' ICQ: improper dihedral type code')
      WRITE (IUNIT, IOSTAT = IRES, FMT = '(A1,A4,4A5)')
     .   '#','IQ','JQ','KQ','LQ','ICQ'
C
      CALL WRQ(IUNIT,PRGSTR,NQHI,NRP,NQTY,
     .     IQ,JQ,KQ,LQ,ICQ,LOK)
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

Cnpty
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NDITY))
      IF (NPTY .LT. 0 .OR. NPTY .GT. MAXPTY) THEN
         PRINT *,PRGSTR,': NPTY out of range!'
         PRINT *,'WRITE ',NPTY,'MAXPTY=',MAXPTY
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF
      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' NPTY: number of dihedral types')
      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NPTY
      IF (IRES .NE. 0)THEN
         PRINT *,PRGSTR,': Error writing NPTY!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C cp,pd,np
      CALL WRCOMM(IUNIT,PRGSTR,' CP: force constant')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' PD: cosine of the phase shift')
      CALL WRCOMM(IUNIT,PRGSTR,' NP: multiplicity')
      WRITE (IUNIT, IOSTAT = IRES, FMT = '(A1,A9,A10,A5)')
     .     '#','CP','PD','NP'

      DO K=1,NPTY
         WRITE (IUNIT, IOSTAT = IRES, FMT = 20)
     .      CP(K),PD(K),NP(K)
         IF (IRES .NE. 0)THEN
            PRINT *,PRGSTR,': Error writing CP,PD,NP!'
            CALL STOPGM('WRTOPO','FATAL ERROR')
         ENDIF
         IF (MOD(K,10) .EQ. 0)THEN
            CALL WRCINT(IUNIT,PRGSTR,K)
         ENDIF
      ENDDO
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C nphih
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NDIHY))
      IF (NPHIH .LT. 0 .OR. NPHIH .GT. MXPHIH) THEN
         PRINT *,PRGSTR,': NPHIH out of range!'
         PRINT *,'WRITE ',NPHIH,'MXPHIH=',MXPHIH
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF
      CALL WRCOMM(IUNIT,PRGSTR,
     . ' NPHIH: number of dihedrals involving H atoms in solute')
      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NPHIH
      IF (IRES .NE. 0)THEN
         PRINT *,PRGSTR,': Error writing NPHIH!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C iph,jph,kph,lph,icph
      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' IPH, JPH, KPH, LPH: atom sequence numbers')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     '   of atoms forming a dihedral')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' ICPH: dihedral type code')
      WRITE (IUNIT, IOSTAT = IRES, FMT = '(A1,A4,4A5)')
     .     '#','IPH','JPH','KPH','LPH','ICPH'

      CALL WRQ(IUNIT,PRGSTR,NPHIH,NRP,NPTY,
     .     IPH,JPH,KPH,LPH,ICPH,LOK)
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C nphi
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NDIN))
      CALL WRCOMM(IUNIT,PRGSTR,
     . ' NPHI: number of dihedrals NOT involving H atoms in solute')
      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NPHI
      IF (IRES .NE. 0)THEN
         PRINT *,PRGSTR,': Error writing NPHI!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C ip,jp,kp,lp,icp
      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' IP, JP, KP, LP: atom sequence numbers')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     '    of atoms forming a dihedral')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' ICP: dihedral type code')
      WRITE (IUNIT, IOSTAT = IRES, FMT = '(A1,A4,4A5)')
     .     '#','IP','JP','KP','LP','ICP'

      CALL WRQ(IUNIT,PRGSTR,NPHI,NRP,NPTY,
     .     IP,JP,KP,LP,ICP,LOK)
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C nratt2
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NLJBLK))
      ILIM =  NRATT*(NRATT+1)/2
      IF (NRATT2 .NE. NRATT2) THEN
         PRINT *,PRGSTR,': inconsistent value for NRATT2'
         PRINT *,'when writing C12, C6, CS12, CS6 !'
         PRINT *,'WRITE',NRATT2,',expected ',ILIM
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF
      CALL WRCOMM(IUNIT,PRGSTR,
     .' NRATT2: number of LJ interaction types = NRATT*(NRATT+1)/2')
      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NRATT2
      IF (IRES .NE. 0)THEN
         PRINT *,PRGSTR,': Error writing NRATT2!'
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C c12,c6,cs12,cs6
      TSTR = ' IAC,JAC: '// IACSTR
      CALL WRCOMM(IUNIT,PRGSTR,TSTR)
      WRITE (TSTR,FMT=155)' C12: r**(-12)'
      CALL WRCOMM(IUNIT,PRGSTR,TSTR)
      WRITE (TSTR,FMT=155)'  C6: r**(-6)'
      CALL WRCOMM(IUNIT,PRGSTR,TSTR)
      WRITE (TSTR,FMT=156)'CS12: r**(-12)'
      CALL WRCOMM(IUNIT,PRGSTR,TSTR)
      WRITE (TSTR,FMT=156)' CS6: r**(-6)'
      CALL WRCOMM(IUNIT,PRGSTR,TSTR)

      WRITE(IUNIT,IOSTAT = IRES,FMT='(A1,A4,A5,4A15)')
     .     '#','IAC','JAC','C12','C6','CS12','CS6'

      K = 1
      DO I=1,NRATT
         DO J=1,I
C            WRITE (IUNIT, IOSTAT = IRES, FMT = '(I5,I5,4F10.5)')
            WRITE (IUNIT, IOSTAT = IRES, FMT = 16)
     .           J,I,C12(K),C6(K),CS12(K),CS6(K)
            IF (IRES .NE. 0)THEN
               PRINT *,PRGSTR,': Error writing C12,C6!'
               CALL STOPGM('WRTOPO','FATAL ERROR')
            ENDIF
            K = K+1
         ENDDO
         CALL WRCOMM(IUNIT,PRGSTR,' ')
      ENDDO
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C---------solvent part ---------------
C nram
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NSLVBK))
      IF (NRAM .LT. 0 .OR. NRAM .GT. MAXNRS) THEN
         PRINT *,PRGSTR,': NRAM out of range!'
         PRINT *,'WRITE ',NRAM,'MAXNRS=',MAXNRS
         CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF
      CALL WRCOMM(IUNIT,PRGSTR,
     .   ' NRAM: number of atoms per solvent molecule')
      WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NRAM
      IF (IRES .NE. 0)THEN
          PRINT *,PRGSTR,': Error writing NRAM!'
          CALL STOPGM('WRTOPO','FATAL ERROR')
      ENDIF

C winvs,anms,cgs,iacs
      CALL WRCOMM(IUNIT,PRGSTR,
     .     '    I: solvent atom sequence number')
      TSTR = ' IACS: '//IACSTR
      CALL WRCOMM(IUNIT,PRGSTR,TSTR)
      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' ANMS: atom name of solvent atom')
      CALL WRCOMM(IUNIT,PRGSTR,
     .     ' MASS: mass of solvent atom' )
      CALL WRCOMM(IUNIT,PRGSTR,
     .     '  CGS: charge of solvent atom')

      WRITE (IUNIT, IOSTAT = IRES, FMT = '(A1,A4,2A6,A9,A11)')
     .     '#','I','ANMS','IACS','MASS','CGS'

      CGDIV = DSQRT(FPEPSI)
      DO K=1,NRAM
         WRITE (IUNIT, IOSTAT = IRES,FMT = 21)
!    .        K,ANMS(K),IACS(K),WMASS(K),CGS(K)
     .        K,ANMS(K),IACS(K),WMASS(K),CGS(K)/CGDIV    !QMMM
         IF (IRES .NE. 0)THEN
            PRINT *,PRGSTR,': Error writing WINVS,ANMS,CGS,IACS!'
            CALL STOPGM('WRTOPO','FATAL ERROR')
         ENDIF
         IF (MOD(K,10) .EQ. 0) THEN
            CALL WRCINT(IUNIT,PRGSTR,K)
         ENDIF
      ENDDO
      CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))

C ncons
      IF (NCONS .EQ. 0) THEN
         CALL WRCOMM(IUNIT,PRGSTR,'no solvent constraints')
      ELSE
         CALL WRSTR(IUNIT,PRGSTR,TPNAME(NSLVCN))
         IF (NCONS .LT. 0 .OR. NCONS .GT. MXCONS)THEN
            PRINT *,PRGSTR,': NCONS out of range!'
            PRINT *,'WRITE ',NCONS,'MXCONS=',MXCONS
            CALL STOPGM('WRTOPO','FATAL ERROR')
         ENDIF
         CALL WRCOMM(IUNIT,PRGSTR,
     .        ' NCONS: number of constraints')
         WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NCONS
         IF (IRES .NE. 0)THEN
            PRINT *,PRGSTR,': Error writing NCONS!'
            CALL STOPGM('WRTOPO','FATAL ERROR')
         ENDIF

C     icons,jcons,cons
         CALL WRCOMM(IUNIT,PRGSTR,
     .     ' ICONS, JCONS: atom sequence numbers forming constraint')
         CALL WRCOMM(IUNIT,PRGSTR,
     .        '  CONS constraint length')
         WRITE (IUNIT, IOSTAT = IRES, FMT = '(A1,A5,A6,A15)')
     .        '#','ICONS','JCONS','CONS'

         DO K=1,NCONS
            WRITE (IUNIT, IOSTAT = IRES, FMT = 22)
!    .           ICONS(K),JCONS(K), CONS(K)
     .           ICONS(K),JCONS(K),DSQRT(CONS(K))        !QMMM
            IF (IRES .NE. 0)THEN
               PRINT *,PRGSTR,': Error writing ICONS,JCONS,CONS!'
               CALL STOPGM('WRTOPO','FATAL ERROR')
            ENDIF
            IF (MOD(K,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,K)
            ENDIF
         ENDDO
         CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))
      ENDIF

C-----path integral part
      IF (NPID .GT. 1) THEN
C     path integral block, values are checked in the invoking program
C     thus no checking here
         CALL WRSTR(IUNIT,PRGSTR,TPNAME(NPITBL))

         CALL WRCOMM(IUNIT,PRGSTR,
     .        ' NPID: number of discretizations per atom')
         WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NPID
         IF (IRES .NE. 0) THEN
            PRINT *,PRGSTR,': Error writing NPID!'
            CALL STOPGM('WRTOPO','FATAL ERROR')
         ENDIF

         CALL WRCOMM(IUNIT,PRGSTR,
     .        ' IPIC: bead codes (0 = classical atom)')
         DO K=1,NRP
            WRITE (IUNIT, IOSTAT = IRES, FMT = 14) IPIC(K)
            IF (IRES .NE. 0) THEN
               PRINT *,PRGSTR,': Error writing IPIC!'
               CALL STOPGM('WRTOPO','FATAL ERROR')
            ENDIF
            IF (MOD(K,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,K)
            ENDIF
         ENDDO
         CALL WRCOMM(IUNIT,PRGSTR,
     .      ' TPI: temperature as basis for the discretization')
         WRITE (IUNIT, IOSTAT = IRES, FMT = 19) TPI
         CALL WRCOMM(IUNIT,PRGSTR,
     .      ' NPIT: number of force constants between two pseudopart.')
         WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NPIT
         IF (IRES .NE. 0) THEN
            PRINT *,PRGSTR,': Error writing NPIT!'
            CALL STOPGM('WRTOPO','FATAL ERROR')
         ENDIF

         CALL WRCOMM(IUNIT,PRGSTR,
     .        ' CPI: force constants between two pseudoparticles')
         DO K=1,NPIT
            WRITE (IUNIT, IOSTAT = IRES, FMT = 19) CPI(K)
            IF (IRES .NE. 0) THEN
               PRINT *,PRGSTR,': Error writing CPI!'
               CALL STOPGM('WRTOPO','FATAL ERROR')
            ENDIF
            IF (MOD(K,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,K)
            ENDIF
         ENDDO
 
         CALL WRCOMM(IUNIT,PRGSTR,
     .    ' NPIB: number of harmonic springs between pseudoparticles')
         WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NPIB
         IF (IRES .NE. 0) THEN
            PRINT *,PRGSTR,': Error writing NPIB!'
            CALL STOPGM('WRTOPO','FATAL ERROR')
         ENDIF

         CALL WRCOMM(IUNIT,PRGSTR,
     .      ' IPIB, JPIB: atom sequence numbers forming a bond')
         CALL WRCOMM(IUNIT,PRGSTR,
     .      ' ICPIB: bond type code')
         CALL WRCOMM(IUNIT,PRGSTR,
     .      ' IPIB  JPIB  ICPIB')
         DO K=1,NPIB
            WRITE (IUNIT, IOSTAT = IRES, FMT = 14)
     .           IPIB(K), JPIB(K), ICPIB(K)
            IF (IRES .NE. 0) THEN
               PRINT *,PRGSTR,': Error writing IPIB,JPIB,ICPIB!'
               CALL STOPGM('WRTOPO','FATAL ERROR')
            ENDIF
            IF (MOD(K,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,K)
            ENDIF
         ENDDO

         CALL WRCOMM(IUNIT,PRGSTR,
     .      ' NPIA: number of atoms that are discretized')
         WRITE (IUNIT, IOSTAT = IRES, FMT = 14) NPIA
         IF (IRES .NE. 0) THEN
            PRINT *,PRGSTR,': Error writing NPIA!'
            CALL STOPGM('WRTOPO','FATAL ERROR')
         ENDIF
         CALL WRCOMM(IUNIT,PRGSTR,
     .      ' WMCL: classical mass of discretized atoms')
         DO K=1,NPIA
            WRITE (IUNIT, IOSTAT = IRES, FMT = 12)WMCL(K)
            IF (IRES .NE. 0) THEN
               PRINT *,PRGSTR,': Error writing WMCL!'
               CALL STOPGM('WRTOPO','FATAL ERROR')
            ENDIF
            IF (MOD(K,10) .EQ. 0) THEN
               CALL WRCINT(IUNIT,PRGSTR,K)
            ENDIF
         ENDDO

         CALL WRSTR(IUNIT,PRGSTR,TPNAME(NTPEND))
      ENDIF

      CALL WRCOMM(IUNIT,PRGSTR,'end of topology file')
      IF (LOK) THEN
         PRINT *,PRGSTR,': OK!'
      ELSE
         PRINT *
         PRINT *,PRGSTR,': ERRORS DETECTED IN TOPOLOGY WRITTEN'
         PRINT *
      ENDIF

C format statements (in order to save space needed for storing strings)
 155  FORMAT (A,' term in nonbonded interactions')
 156  FORMAT (A,' term in 1-4 nonbonded interactions')
 165  FORMAT (A,': Error writing solute atom data')

 12   FORMAT (E15.7)
 13   FORMAT (F3.1)
 14   FORMAT (12I7)
 15   FORMAT (A5)
 16   FORMAT (2I6,4E15.7)
 17   FORMAT (2I6,A6,I4,2F11.5,2I3,4I7)
 18   FORMAT (45X,I3,6I6)
 19   FORMAT (E15.7,F15.7)
 20   FORMAT (2F10.5,I6)
 21   FORMAT (I6,A6,I4,2F11.5)
 22   FORMAT (2I6,F15.7)
C end WRTOPO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE WBND(IUNIT,PRGSTR,NBONH,NRP,NBTY,IBH,JBH,ICBH,LOK)
C     ==--------------------------------------------------------------==
      INTEGER IUNIT
      LOGICAL LOK
      CHARACTER PRGSTR*(*)
C
      INTEGER NBONH,NRP,NBTY
      INTEGER IBH(NBONH),JBH(NBONH),ICBH(NBONH)
C local vars
      INTEGER K,I,IRES
      LOGICAL LERR
C begin
      DO K=1,NBONH
         WRITE (IUNIT, IOSTAT = IRES, FMT = '(3I7)')
     $      IBH(K),JBH(K),ICBH(K)
         IF (IRES .NE. 0 )THEN
            PRINT *,PRGSTR,': Error writing IBH,JBH,ICBH!'
            CALL STOPGM('WBND','FATAL ERROR')
         ENDIF
         IF (MOD(K,10) .EQ. 0)THEN
            CALL WRCINT(IUNIT,PRGSTR,K)
         ENDIF
      ENDDO

C check ibh,jbh,icbh
      LERR = .FALSE.
      DO 20 I=1,NBONH
         IF (IBH(I) .LT. 1 .OR. IBH(I) .GT. NRP) THEN
            LERR = .TRUE.
            PRINT *,PRGSTR,': illegal atom number in IBH!'
            PRINT *,'IBH(',I,') = ',IBH(I)
            PRINT *,'must be 0 < IBH(I) <= NRP = ',NRP
            PRINT *
         ENDIF

         IF (JBH(I) .LT. 1 .OR. JBH(I) .GT. NRP) THEN
            LERR = .TRUE.
            PRINT *,PRGSTR,': illegal atom sequence number in JBH!'
            PRINT *,'JBH(',I,') = ',JBH(I)
            PRINT *,'must be 0 < JBH(I) <= NRP = ',NRP
            PRINT *
         ENDIF

         IF (ICBH(I) .LT. 1 .OR. ICBH(I) .GT. NBTY) THEN
            LERR = .TRUE.
            PRINT *,PRGSTR,': illegal bond type in ICBH!'
            PRINT *,'ICBH(',I,') = ',ICBH(I)
            PRINT *,'must be 0 < ICBH(I) <= NBTY =',NBTY
            PRINT *
         ENDIF
 20   CONTINUE

      IF (LERR) THEN
         LOK = .FALSE.
      ENDIF
C end wbnd
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE WRBAN(IUNIT,PRGSTR,NTHEH,NRP,NTTY,
     $     ITH,JTH,KTH,ICTH,LOK)
C     ==--------------------------------------------------------------==
C args
      LOGICAL LOK
      INTEGER IUNIT,NTHEH,NRP,NTTY
      CHARACTER PRGSTR*(*)
      INTEGER ITH(NTHEH),JTH(NTHEH),KTH(NTHEH),ICTH(NTHEH)
C local vars
      LOGICAL LERR
      INTEGER K,IRES,I
C begin
      LERR = .FALSE.
      DO K=1,NTHEH
         WRITE (IUNIT, IOSTAT = IRES, FMT = '(4I7)')
     .      ITH(K),JTH(K),KTH(K),ICTH(K)
         IF (IRES .NE. 0)THEN
            PRINT *,PRGSTR,': Error writing ITH,JTH,KTH,ICTH'
            CALL STOPGM('WRBAN','FATAL ERROR')
         ENDIF
         IF (MOD(K,10) .EQ. 0)THEN
            CALL WRCINT(IUNIT,PRGSTR,K)
         ENDIF
      ENDDO

C check ith,jth,kth,icth
      DO 40 I = 1,NTHEH

         IF (ITH(I) .LT. 1 .OR. ITH(I) .GT. NRP) THEN
            LERR = .TRUE.
            PRINT *,PRGSTR,': illegal atom sequence number in ITH!'
            PRINT *,'ITH(',I,') = ',ITH(I)
            PRINT *,'must be 0 < ITH(I) <= NRP = ',NRP
            PRINT *
         ENDIF

         IF (JTH(I) .LT. 1 .OR. JTH(I) .GT. NRP) THEN
            LERR = .TRUE.
            PRINT *,PRGSTR,': illegal atom sequence number in JTH!'
            PRINT *,'JTH(',I,') = ',JTH(I)
            PRINT *,'must be 0 < JTH(I) <= NRP = ',NRP
            PRINT *
         ENDIF

         IF (KTH(I) .LT. 1 .OR. KTH(I) .GT. NRP) THEN
            LERR = .TRUE.
            PRINT *,PRGSTR,': illegal atom sequence number in KTH!'
            PRINT *,'KTH(',I,') = ',KTH(I)
            PRINT *,'must be 0 < KTH(I) <= NRP = ',NRP
            PRINT *
         ENDIF

         IF (ICTH(I) .LT. 1 .OR. ICTH(I) .GT. NTTY) THEN
            LERR = .TRUE.
            PRINT *,PRGSTR,': illegal bond type in ICTH!'
            PRINT *,'ICTH(',I,') = ',ICTH(I)
            PRINT *,'must be 0 < ICTH(I) <= NTTY =',NTTY
            PRINT *
         ENDIF
 40   CONTINUE

      IF (LERR) THEN
         LOK = .FALSE.
      ENDIF
C end wrban
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE  WRQ(IUNIT,PRGSTR,NQHIH,NRP,NQTY,
     $     IQH,JQH,KQH,LQH,ICQH,LOK)
C     ==--------------------------------------------------------------==
C args
      LOGICAL LOK
      INTEGER IUNIT,NQHIH,NRP,NQTY
      CHARACTER PRGSTR*(*)
      INTEGER IQH(NQHIH),JQH(NQHIH),KQH(NQHIH)
      INTEGER LQH(NQHIH),ICQH(NQHIH)
C local vars
      LOGICAL LERR
      INTEGER I,K,IRES
C begin
      LERR = .FALSE.

      DO K=1,NQHIH
         WRITE (IUNIT, IOSTAT = IRES, FMT = '(5I7)')
     .      IQH(K),JQH(K),KQH(K),LQH(K),ICQH(K)
         IF (IRES .NE. 0)THEN
            PRINT *,PRGSTR,': Error writing IQH,JQH,KQH,LQH,ICQH!'
            CALL STOPGM('WRQ','FATAL ERROR')
         ENDIF
         IF (MOD(K,10) .EQ. 0)THEN
            CALL WRCINT(IUNIT,PRGSTR,K)
         ENDIF
      ENDDO

C check iqh,jqh,kqh,lqh,icqh
      DO 60 I=1,NQHIH

         IF (IQH(I) .LT. 1 .OR. IQH(I) .GT. NRP) THEN
            LERR = .TRUE.
            PRINT *,PRGSTR,': illegal atom sequence number in IQH!'
            PRINT *,'IQH(',I,') = ',IQH(I)
            PRINT *,'must be 0 < IQH(I) <= NRP = ',NRP
            PRINT *
         ENDIF

         IF (JQH(I) .LT. 1 .OR. JQH(I) .GT. NRP) THEN
            LERR = .TRUE.
            PRINT *,PRGSTR,': illegal atom sequence number in JQH!'
            PRINT *,'JQH(',I,') = ',JQH(I)
            PRINT *,'must be 0 < JQH(I) <= NRP = ',NRP
            PRINT *
         ENDIF

         IF (KQH(I) .LT. 1 .OR. KQH(I) .GT. NRP) THEN
            LERR = .TRUE.
            PRINT *,PRGSTR,': illegal atom sequence number in KQH!'
            PRINT *,'KQH(',I,') = ',KQH(I)
            PRINT *,'must be 0 < KQH(I) <= NRP = ',NRP
            PRINT *
         ENDIF

         IF (LQH(I) .LT. 1 .OR. LQH(I) .GT. NRP) THEN
            LERR = .TRUE.
            PRINT *,PRGSTR,': illegal atom sequence number in LQH!'
            PRINT *,'LQH(',I,') = ',LQH(I)
            PRINT *,'must be 0 < LQH(I) <= NRP = ',NRP
            PRINT *
         ENDIF

         IF (ICQH(I) .LT. 1 .OR. ICQH(I) .GT. NQTY) THEN
            LERR = .TRUE.
            PRINT *,PRGSTR,': illegal bond type in ICQH!'
            PRINT *,'ICQH(',I,') = ',ICQH(I)
            PRINT *,'must be 0 < ICQH(I) <= NQTY =',NQTY
            PRINT *
         ENDIF
 60   CONTINUE

      IF (LERR) THEN
         LOK = .FALSE.
      ENDIF
C end WRQ
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PWRTOP
C     ==--------------------------------------------------------------==
C     Reverses the changes that RDTOPO does to the values
C     in the topology on reading.
C     Thus, calling RDTOPO, then PWRTOP and then WRTOPO
C     will produce the same topology file.
C     ==--------------------------------------------------------------==
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
C local vars
      INTEGER K
      REAL*8 RD2DEG,CGDIV,RDSGSQ
C begin
      RD2DEG = 45.0D0/DATAN(1.0D0)
      RDSGSQ = RD2DEG*RD2DEG

C bonds: no conversion needed

C bond angles:
C     T0 from COS(angle) --> angle in degrees
C     CT unchanged
!$OMP parallel do private(K)
      DO K=1,NTTY
         T0(K) = RD2DEG*DACOS(T0(K))
      ENDDO

C impropers:
C     Q0 from angles in rad --> angles in degrees
C     CQ from energy/ (rad^2) --> energy/(degrees^2)
!$OMP parallel do private(K)
      DO K=1,NQTY
         CQ(K) = CQ(K)/RDSGSQ
         Q0(K) = Q0(K)*RD2DEG
      ENDDO

C dihedrals: no conversion needed

C atom charges: CG from E/sqrt(FPEPSI) --> E
      CGDIV = DSQRT(FPEPSI)
!$OMP parallel do private(K)
      DO K=1,NRP
         CG(K)  = CG(K)/CGDIV
      ENDDO

!$OMP parallel do private(K)
      DO K=1,NRAM
         CGS(K)  = CGS(K)/CGDIV
      ENDDO

!$OMP parallel do private(K)
      DO K=1,NCONS
         CONS(K) = DSQRT(CONS(K))
      ENDDO

C end pwrtop
C     ==--------------------------------------------------------------==
      RETURN
      END
