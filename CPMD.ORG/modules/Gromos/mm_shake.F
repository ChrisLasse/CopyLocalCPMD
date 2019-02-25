C     ==================================================================
      SUBROUTINE mm_SHAKE(NATTOT,NR,NM,NCALCD,WINVSH,NSKIP,
     .     NCON,ICO,JCO,CONSTR,FCONSH,SHKTOL,XREF,XP,C4D,
     .     NCC,JRC,LEVERY,NITER,LSOLV)
C     ==--------------------------------------------------------------==
COMMSUBR mm_SHAKE W.F. VAN GUNSTEREN, GRONINGEN, AUG. 1986
C Id: shake.f,v 1.18 1995/12/21 12:27:01 wscott Exp
C     new version by Alan Mark and Walter Scott
Cmb - Revised on 16 Sept. 2005
C
C     SUBROUTINE mm_SHAKE(NATTOT,NR,NM,NCALCD,WINVSH,NSKIP,
C     $     NCON,ICO,JCO,CONSTR,FCONSH,
C     $     SHKTOL,XREF,XP,C4D,
C     $     NCC,JRC,LEVERY,NITER,LSOLV)
C
C     SHAKE will supply corrections to given coordinates XP, such
C     that for XP a list of constraints will be satisfied, each
C     geometrically within a specified tolerance SHKTOL.
C     The corrections are made along vectors derived from XREF which
C     is not modified.
C
C     When XP=XP(T+DT) results from a non constrained dynamics
C     step of length DT, starting at XREF=X(T), the corrections
C     are dynamically correct up to an order of DT**2.
C     Thus, the dynamic accuracy of SHAKE depends on both
C     SHKTOL and DT, whereas the geometrical accuracy only
C     depends on SHKTOL.
C
C     Description of shake: J. Comp. Phys. 23 (1977) 327.
C
C     Periodic boundary conditions can be taken into account
C     using L<BOX>,L<NTB> and L<BETA>.
C     More than one identical molecule can be treated.
C
C     If MAXIT or more steps appear necessary, or if XP
C     deviates from XREF (DT too large), a message is printed
C     and NITER is returned with zero.
C
C
C     NATTOT         total number of atoms in the X and XP arrays
C     NR             number of atoms per molecule
C     NM             number of molecules to shake.
C     WINVSH(NR)     inverse atomic masses
C     NSKIP          number of atoms in X and XP preceding the
C                    first atom of the first molecule
C     NCC            number of position constrained atoms
C     JRC(NCC)       index number of positionally constrained atoms.
C                    If NCC > 1 then NM must be 1.
C                    If atoms I or J are positionally constrained,
C                    the constraint IJ is not shaken.
C
C     NCON           number of elements in ICO,JCO and CONSTR
C     ICO, JCO(NCON) OF [1..NR]
C                    topology sequence numbers of atoms I and J
C                    forming constraints
C     CONSTR(NCON)   squares of the constraint distances
C     FCONSH(NCON)   sum of ACOR (see program code) for each constraint,
C                    used to determine the constraint force.
C                    This is used in free energy calculations.
C
C     SHKTOL         relative geometric tolerance for coord. resetting
C
C     XREF(NATTOT*NDIM) reference atom cartesian coordinates
C     XP(NATTOT*NDIM)atom cartesian coordinates that will be shaken
C
C     X and Y axes lie along box edges.
C     Corresponding atoms in X and XP must not be separated by more
C     than BOX/2 or BOX*SQRT(3)/4.
C
C     LEVERY        .TRUE.: quantities depending on L<BETA> are
C                           calculated I<every time> SHAKE is envoked.
C                   .FALSE.: quantities depending on L<BETA> are
C                           only calculated I<the first time>
C                           SHAKE is envoked.
C     NITER         delivered with the number of iterations performed.
C                   The value of zero is returned if an error occurs.
C
C     LSOLV         .TRUE.: solvent is being shaken.  Constraints do not 
C                           cross BOX.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'toposz.h'
      INCLUDE 'box.h'
C args
      INTEGER NATTOT,NR,NM,NSKIP,NCON,NITER,NCC,NCALCD
      INTEGER ICO(*),JCO(*),JRC(*)
      LOGICAL LEVERY,LSOLV
      REAL*8 WINVSH(*),CONSTR(*),FCONSH(*)
      REAL*8 XREF(*),XP(*),SHKTOL
      REAL*8 C4D(*)
C local params
      REAL*8 EPS
      PARAMETER (EPS = 1.0D-6)
      INTEGER MAXIT
      PARAMETER (MAXIT = 1000)
C commons

C----SHAKECOM -- be sure to modify this consistently within this file
      REAL*8 BOXOH,BOXOQ
      LOGICAL LMONO,LOCTO,LDOTRA
      COMMON /SHKCOM/BOXOH,BOXOQ,LMONO,LOCTO,LDOTRA
C---- 
C these are used for calculation the dimension of each constraint
      INTEGER NCDCO(MAXCON)
      LOGICAL LI4,LJ4
      INTEGER NCLCLO

C local vars
      LOGICAL LFIRST,LVAC
      SAVE LFIRST,LVAC

C these are 'I skip now' (skip a constraint in this shake iteration) and
C 'I skip next' ( skip a constraint in the next shake iteration)
      LOGICAL LSKNOW(MAXNRP),LSKNEX(MAXNRP)
      LOGICAL LFINI,LDOCOR
      INTEGER M,NMOL,NUMIT,NFIRST,NC,I,J,I3,J3,ICOFF,II
      REAL*8 BETAR,COSB,BOXHM
C
C     XPIJ: the distance between atoms I and J in the XP array
C     XREFIJ: the distance between atoms I and J in the XREF array
      REAL*8 XPIJ(MAXDIM),XREFIJ(MAXDIM)
      REAL*8 TOLER,DIFF,RPIJ2,BOXHM2,TOL2,RRPR,ACOR,XH

      CHARACTER PRGSTR*(5)
C     data
      DATA LFIRST/.TRUE./
      DATA PRGSTR /'SHAKE'/
C begin
      NITER = 0

      IF (LFIRST .OR. LEVERY) THEN
         LFIRST = .FALSE.
         
         LMONO = (NTB .GT. 0)
         LOCTO = (NTB. LT. 0)
         LVAC  = (NTB .EQ. 0)
         LDOTRA = .FALSE.

         IF (LMONO)THEN
C     we have a monoclinic box
            BETAR = ATAN(1.D0)*BETA/4.5D+1
            COSB = DCOS(BETAR)
            LDOTRA = (DABS(COSB).GE.1.D-4)
         ENDIF
      ENDIF

C Check if nothing to do!
      IF (NCON .EQ. 0)THEN
         NITER = 1
         RETURN
      ENDIF

      TOL2 = SHKTOL*2.0D0

C     TEST FOR THE MOST COMMON CASE 
C     ie: NCALCD is three we dont have a monoclinic box and
C     we are shaking solvent.
C     This loop could be inverted for vectorization.

      IF (LSOLV .AND. (.NOT. LDOTRA) .AND. NCALCD .EQ. NDRMAX) THEN
C     In this case skip
C     SHAKE FORCES
C     POSITION RESTRAINTS
C     BOX CHECKS (SOLVENT CANNOT CROSS BOX AND WE DO NOT
C     HAVE TO CALL TRACO)
C     4D CHECKS
C     INSTEAD OF SETTING NCLCLO TO THREE,WE JUST USE NDRMAX
C     AS THIS IS A PARAMETER, NOT A VARIABLE.
C     THIS SHOULD MAKE IT EASIER FOR COMPILERS TO UNROLL LOOPS

         DO 100 NMOL= 1,NM

C     the atom BEFORE the first atom of this molecule
            NFIRST = NSKIP + NR*(NMOL-1)
C     the coord offset

C     this (NDIM+1) less than the first coord in the array
            ICOFF = NDIM*(NFIRST-1)

            LFINI = .FALSE.
            NUMIT = 0
!$OMP parallel do private(M)
            DO M=1,NR
               LSKNEX(M) = .TRUE.
               LSKNOW(M) = .FALSE.
            ENDDO

C     while not finished (=convergence) do
 110        IF (.NOT. LFINI ) THEN

C     set LFINI to .TRUE.
C     if all constraints are skipped then LFINI will stay .TRUE.
               LFINI = .TRUE.
               DO 120 NC=1,NCON
                  I = ICO(NC)
                  J = JCO(NC)
C------begin of constraint block
                  IF (.NOT. (LSKNOW(I) .AND. LSKNOW(J))) THEN
                     TOLER = CONSTR(NC)
                     I3 = NDIM*I + ICOFF
                     J3 = NDIM*J + ICOFF
C     calc distance between I and J in XP array
                     RPIJ2 = 0.0D0
!$OMP parallel do private(M) reduction(+:RPIJ2)
                     DO M=1,NDRMAX
                        XPIJ(M) = XP(I3+M) - XP(J3+M)
                        RPIJ2 = RPIJ2 + XPIJ(M)*XPIJ(M)
                     ENDDO
                     DIFF = TOLER - RPIJ2
                     IF (DABS(DIFF) .GE. TOLER*TOL2) THEN
C     we have to shake
C     calc distance between I and J in XREF array and do the shake
                        RRPR = 0.0D0
!$OMP parallel do private(M) reduction(+:RRPR)
                        DO M=1,NDRMAX
                           XREFIJ(M) = XREF(I3+M) - XREF(J3+M)
                           RRPR = RRPR + XREFIJ(M)*XPIJ(M)
                        ENDDO

                        IF (RRPR .LT. TOLER*EPS) THEN
                           PRINT 500,
     .                          PRGSTR,RRPR,RPIJ2,DIFF,
     .                          NDIM,NDRMAX,NITER,NMOL,
     .                          NC,ICO(NC),JCO(NC),CONSTR(NC),
     .                          NSKIP,NFIRST,I,J,I3,J3
                        
                           PRINT 570
                           DO M=1,NDRMAX
                              PRINT 580,XP(I3+M),XP(J3+M),
     .                              XREF(I3+M),XREF(J3+M),
     .                              XREFIJ(M),XPIJ(M)
                           ENDDO
                           PRINT *
                        
                           NITER = 0
                           RETURN
                        ENDIF
                        ACOR=DIFF/(RRPR*(WINVSH(I)+WINVSH(J))*2.D0)
!$OMP parallel do private(M,XH)
                        DO M=1,NDRMAX
                           XH = XREFIJ(M)*ACOR
                           XP(I3+M) = XP(I3+M) + XH*WINVSH(I)
                           XP(J3+M) = XP(J3+M) - XH*WINVSH(J)
                        ENDDO

                        LSKNEX(I) = .FALSE.
                        LSKNEX(J) = .FALSE.
                        LFINI = .FALSE.
                     ENDIF
C
C     --- end of skip this constraint
                  ENDIF
C     ---- end of constraint block
C     
C     --- end of loop over NCON
 120           CONTINUE

               NUMIT = NUMIT + 1
               IF (NUMIT .GE. MAXIT) THEN
                  PRINT 550,PRGSTR,MAXIT,SHKTOL
                  NITER = 0
                  RETURN
               ENDIF
!$OMP parallel do private(M)
               DO M=1,NR
                  LSKNOW(M) = LSKNEX(M)
                  LSKNEX(M) = .TRUE.
               ENDDO

C     end of while
               GOTO 110
            ENDIF
      
            NITER = NITER + NUMIT

C     end of DO loop for each molecule
 100     CONTINUE

C     return the average NITER value
         NITER = NITER/NM

C   finished solvent return to main program
         RETURN
      ELSE
C -------------------------------------------------------------
C ----------else now do everything which is not solvent -------
C -------------------------------------------------------------

C check for position constraint consistency
         IF (NCC .NE. 0 .AND. NM .NE. 1) THEN
            PRINT *,PRGSTR,': NM must be 1 iff NCC > 0!'
            CALL STOPGM('mm_shake','FATAL ERROR')
         ENDIF

         IF (.NOT. LVAC) THEN
C     calc the minimum box dim
            BOXHM = BOXH(1)
            DO M=1,NDIM
               IF(BOXH(M) .LT. BOXHM) BOXHM = BOXH(M)
            ENDDO
            BOXHM2 = BOXHM*BOXHM
            IF (LOCTO) THEN
               BOXOH  = BOXH(1)
               BOXOQ  = BOX(1)*0.75D0
               BOXHM2 = BOXHM2*0.75D0
            ENDIF
         ENDIF

C Initialize the array to return shake forces.   This is done 
C outside the loop over molecules to save time.  Only one
C additional addition is then required.  Note it is not possible to
C keep track of individual molecules.  It does not effect the 
C calculation of free energy as the contributions from different 
C molecules are additive as implemented.
         CALL mm_AZZERO(FCONSH,NCON)

C if NCALCD is three, set all the ND* arrays to three
C here. If NCALCD = 4, then we have to calculate the ND arrays
C for each moledule separately

         IF (NCALCD .EQ. NDRMAX) THEN
!$OMP parallel do private(I)
            DO I=1,NCON
               NCDCO(I) = NDRMAX
            ENDDO
         ENDIF

C start loop over molecules 

         DO 300 NMOL= 1,NM

C the atom BEFORE the first atom of this molecule
            NFIRST = NSKIP + NR*(NMOL-1)
C the coord offset
C this (NDIM+1) less than the first coord in the array
            ICOFF = NDIM*(NFIRST-1)

C set the NCD arrays for this molecule if needed
            IF (NCALCD .NE. NDRMAX) THEN
               DO I=1,NCON
                  LI4 = (C4D(ICO(I)+NFIRST) .GE. 0.0D0)
                  LJ4 = (C4D(JCO(I)+NFIRST) .GE. 0.0D0)
                  IF (LI4 .AND. LJ4) THEN
                     NCDCO(I) = NDHMAX
                  ELSE
                     NCDCO(I) = NDRMAX
                  ENDIF
               ENDDO
            ENDIF

            LFINI = .FALSE.
            NUMIT = 0
!$OMP parallel do private(M) 
            DO M=1,NR
               LSKNEX(M) = .TRUE.
               LSKNOW(M) = .FALSE.
            ENDDO

C     while not finished (=convergence) do
 310        IF (.NOT. LFINI ) THEN

C     set LFINI to .TRUE.
C     if we skip all constraints in the following loops
C     then LFINI will stay .TRUE.
               LFINI = .TRUE.
            
C     for every iteration, set logical LSKNOW for each constrained
C     atom to .true.
!$OMP parallel do private(M) 
               DO M = 1, NCC
                  LSKNOW(JRC(M)) =  .TRUE.
               ENDDO
            
               DO NC=1,NCON
                  I = ICO(NC)
                  J = JCO(NC)
                  NCLCLO = NCDCO(NC)
C------begin of constraint block
                  IF (.NOT. (LSKNOW(I) .AND. LSKNOW(J))) THEN
                     TOLER = CONSTR(NC)
                     I3 = NDIM*I + ICOFF
                     J3 = NDIM*J + ICOFF
C     calc distance between I and J in XP array
                     RPIJ2 = 0.0D0
!$OMP parallel do private(M) reduction(+:RPIJ2)
                     DO M=1,NCLCLO
                        XPIJ(M) = XP(I3+M) - XP(J3+M)
                        RPIJ2   = RPIJ2 + XPIJ(M)*XPIJ(M)
                     ENDDO

                     LDOCOR = (.NOT. LVAC .AND. RPIJ2 .GE. BOXHM2)

C---  correct distance for PBC if necessary and recalc RPIJ2
                     IF (LDOCOR) THEN
                        CALL CORDST(LEVERY,NCLCLO,XPIJ)
                        RPIJ2 = 0.0D0
!$OMP parallel do private(M) reduction(+:RPIJ2)
                        DO M=1,NCLCLO
                           RPIJ2   = RPIJ2 + XPIJ(M)*XPIJ(M)
                        ENDDO
                     ENDIF

                     DIFF = TOLER - RPIJ2
                     IF (DABS(DIFF) .GE. TOLER*TOL2) THEN
C     we have to shake
C     calc distance between I and J in XREF array
!$OMP parallel do private(M)
                        DO M=1,NCLCLO
                           XREFIJ(M) = XREF(I3+M) - XREF(J3+M)
                        ENDDO
                        IF (LDOCOR) THEN
                           CALL CORDST(LEVERY,NCLCLO,XREFIJ)
                        ENDIF
C     --- do the shake
                        RRPR = 0.0D0
!$OMP parallel do private(M) reduction(+:RRPR)
                        DO M=1,NCLCLO
                           RRPR = RRPR + XREFIJ(M)*XPIJ(M)
                        ENDDO
                        IF (RRPR .LT. TOLER*EPS) THEN
                           PRINT 500,
     .                          PRGSTR,RRPR,RPIJ2,DIFF,
     .                          NDIM,NCLCLO,NITER,NMOL,
     .                          NC,ICO(NC),JCO(NC),CONSTR(NC),
     .                          NSKIP,NFIRST,I,J,I3,J3
                        
                           PRINT 570
                           DO M=1,NCLCLO
                              PRINT 580,XP(I3+M),XP(J3+M),XREF(I3+M),
     .                             XREF(J3+M),XREFIJ(M),XPIJ(M)
                           ENDDO
                           PRINT *
                        
                           NITER = 0
                           RETURN
                        ENDIF
                        ACOR = DIFF/(RRPR*(WINVSH(I)+WINVSH(J))*2.0D0)

C     Sum ACOR into FCONSH array for free energy calculations

                        FCONSH(NC) = FCONSH(NC) + ACOR
!$OMP parallel do private(M,XH)
                        DO M=1,NCLCLO
                           XH = XREFIJ(M)*ACOR
                           XP(I3+M) = XP(I3+M) + XH*WINVSH(I)
                           XP(J3+M) = XP(J3+M) - XH*WINVSH(J)
                        ENDDO

                        LSKNEX(I) = .FALSE.
                        LSKNEX(J) = .FALSE.
                        LFINI = .FALSE.
                     ENDIF
C     --- end of skip this constraint
                  ENDIF
C     ---- end of constraint block
C     
C     --- end of loop over NCON
               ENDDO


               NUMIT = NUMIT + 1
               IF (NUMIT .GE. MAXIT) THEN
                  PRINT 550,PRGSTR,MAXIT,SHKTOL
                  NITER = 0
                  RETURN
               ENDIF

!$OMP parallel do private(M) 
               DO M=1,NR
                  LSKNOW(M) = LSKNEX(M)
                  LSKNEX(M) = .TRUE.
               ENDDO
C     end of while

               GOTO 310
            ENDIF
            NITER = NITER + NUMIT
C     end of DO loop for each molecule
 300     CONTINUE

C     return the average NITER value
         NITER = NITER/NM
C     end non solvent shake
      ENDIF

C------format statements ----------
 500  FORMAT (1X,A6,
     . ': COORD RESET CANNOT BE ACCOMPLISHED, DEVIATION IS TOO LARGE',
     . 1X,'RRPR      =',F8.3,/
     . 1X,'RPIJ2     =',F8.3,/
     . 1X,'DIFF      =',F8.3,/
     . 1X,'NDIM      =',I5,/
     . 1X,'NCLCLO    =',I5,/
     . 1X,'NITER     =',I5,/
     . 1X,'NMOL      =',I5,/
     . 1X,'NC        =',I5,/
     . 1X,'IX(NC)    =',I5,/
     . 1X,'JX(NC)    =',I5,/
     . 1X,'CONSTR(NC)=',F8.3,/
     . 1X,'NSKIP     =',I5,/
     . 1X,'NFIRST    =',I5,/
     . 1X,'I         =',I5,/
     . 1X,'J         =',I5,/
     . 1X,'I3        =',I5,/
     . 1X,'J3        =',I5)

 550  FORMAT (1X,A6,
     .     ': COORDINATE RESETTING WAS NOT ACCOMPLISHED WITHIN'/
     .     1X,  'MAXITER = ',I8,' STEPS'/
     .     1X,  'SHKTOL  = ',F8.6/)

 570  FORMAT ('   XP(I3+M)','   XP(J3+M)',' XREF(I3+M)',
     $     ' XREF(J3+M)','  XREFIJ(M)','    XPIJ(M)')
 580  FORMAT (6F11.4)
C end shake
      END
C     ==================================================================
      SUBROUTINE CORDST(LEVERY,NCALCD,XPIJ)
C     ==--------------------------------------------------------------==
C Id: shake.f,v 1.18 1995/12/21 12:27:01 wscott Exp
C     subroutine CORDST corrects distances for periodic
C     boundary conditions.
C     It is an auxillary routine for L<SHAKE>.
C     ==--------------------------------------------------------------==
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'box.h'
C args
      LOGICAL LEVERY
      INTEGER NCALCD
      real*8 XPIJ(NDIM)

C----SHAKECOM -- be sure to modify this consistently within this file
      real*8 BOXOH,BOXOQ
      LOGICAL LMONO,LOCTO,LDOTRA
      COMMON /SHKCOM/BOXOH,BOXOQ,LMONO,LOCTO,LDOTRA
C---- 
C local vars
      INTEGER M
      real*8 ABSSUM
C begin
      IF (LMONO) THEN
         IF (LDOTRA) THEN
            CALL TRACO(1,0,XPIJ,BETA,1,LEVERY)
         ENDIF
C make sure XPIJ is shortest distance
!$OMP parallel do private(M)
         DO M=1,NCALCD
            IF (XPIJ(M) .GE. BOXH(M)) THEN
               XPIJ(M) = XPIJ(M) - BOX(M)
            ELSEIF (XPIJ(M) .LT. -BOXH(M)) THEN
               XPIJ(M) = XPIJ(M) + BOX(M)
            ENDIF
         ENDDO

         IF (LDOTRA) THEN
            CALL TRACO(1,0,XPIJ,BETA,-1,LEVERY)
         ENDIF
      ELSEIF (LOCTO ) THEN
C make sure XPIJ is shortest distance
!$OMP parallel do private(M)
         DO M=1,NCALCD
            IF (XPIJ(M) .GE. BOXH(M)) THEN
               XPIJ(M) = XPIJ(M) - BOX(M)
            ELSEIF (XPIJ(M) .LT. -BOXH(M)) THEN
               XPIJ(M) = XPIJ(M) + BOX(M)
            ENDIF
         ENDDO

         ABSSUM=DABS(XPIJ(1))+DABS(XPIJ(2))+DABS(XPIJ(3))

         IF (ABSSUM .GT. BOXOQ) THEN
!$OMP parallel do private(M)
            DO M = NDRMIN,NDRMAX
               XPIJ(M) = XPIJ(M) - SIGN(BOXOH,XPIJ(M))
            ENDDO
         ENDIF
      ENDIF
C end CORDST
C     ==--------------------------------------------------------------==
      RETURN
      END
