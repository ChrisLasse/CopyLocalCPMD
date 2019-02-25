C hello emacs, this is -*- fortran -*-
C      GROMOS UNIT DEFINITIONS
C

C     standard FORTRAN units
      INTEGER ISTDIN,ISTDOT
      PARAMETER (ISTDIN = 1005,ISTDOT = 6)

C-------define reserved unit numbers--------------
C     md run control file
      INTEGER IOMDCT
      PARAMETER (IOMDCT = ISTDIN)

C     final coords and velocities
      INTEGER IOXVE
      PARAMETER (IOXVE = 1011)

C     coordinates trajectory
      INTEGER IOTRJX
      PARAMETER (IOTRJX = 1012)

C     velocities trajectory
      INTEGER IOTRJV
      PARAMETER (IOTRJV = 1013)

C reserved for future use
      INTEGER IORSVD
      PARAMETER (IORSVD = 1014)

C     energies etc. trajectory
      INTEGER IOTRJE
      PARAMETER (IOTRJE = 1015)

C     free energies
      INTEGER IOTRJG
      PARAMETER (IOTRJG = 1016)



C input files start from unit 20 upwards
C     formatted GROMOS95 topology files
      INTEGER IOTOPO
      PARAMETER (IOTOPO = 1020)

C     initial coordinates and velocities
      INTEGER IOXVI
      PARAMETER (IOXVI = 1021)

C     reference coordinates for position re(con)straining
      INTEGER IOREST
      PARAMETER (IOREST = 1022)

C     sequence numbers of position re(con)strained atoms
      INTEGER IORSTA
      PARAMETER (IORSTA = 1023)

C     distance restrained atom pairs
      INTEGER IORSTP
      PARAMETER (IORSTP = 1024)

C     restrained dihedral specifications
      INTEGER IODIHE
      PARAMETER (IODIHE = 1025)

C     j-value restraining specifications
      INTEGER IOJVSP
      PARAMETER (IOJVSP = 1026)

C     local elevation dihedral specifications
      INTEGER IOLESP
      PARAMETER (IOLESP = 1027)

C     sequence numbers indicating atoms in 4D
      INTEGER IO4NDX
      PARAMETER (IO4NDX = 1028)

C     friction coefficients for SD
      INTEGER IOGAM
      PARAMETER (IOGAM = 1029)

C     data determining perturbation
      INTEGER IOPERT
      PARAMETER (IOPERT = 1030)

C     some conversion programs (GROMOS87 --> GROMOS95 file format)
C     read from unit 40
      INTEGER IOCNV
      PARAMETER (IOCNV = 1040)


C     numbers 50 to 59 are reserved for developers

      INTEGER IOD00,IOD01,IOD02,IOD03,IOD04
      INTEGER IOD05,IOD06,IOD07,IOD08,IOD09
      PARAMETER (IOD00 = 1050)
      PARAMETER (IOD01 = 1051)
      PARAMETER (IOD02 = 1052)
      PARAMETER (IOD03 = 1053)
      PARAMETER (IOD04 = 1054)
      PARAMETER (IOD05 = 1055)
      PARAMETER (IOD06 = 1056)
      PARAMETER (IOD07 = 1057)
      PARAMETER (IOD08 = 1058)
      PARAMETER (IOD09 = 1059)




C------consecutive numbers for table entries
C     The tables are used for mapping a string
C     to a unit number in SUBR. OPNFIL.
C     The table is scanned in a linear fashion which allows
C     the entries to be in any order.

      INTEGER JNMDCT
      INTEGER JNTRJX
      INTEGER JNTRJV
      INTEGER JNTRJE
      INTEGER JNTOPO
      INTEGER JNXVI
      INTEGER JNREST
      INTEGER JNRSTA
      INTEGER JNRSTP
      INTEGER JNPERT
      INTEGER JNDIHE
      INTEGER JNXVE
      INTEGER JND00,JND01,JND02,JND03,JND04
      INTEGER JND05,JND06,JND07,JND08,JND09
      INTEGER JN4COR,JN4NDX
      INTEGER JNCNV
      INTEGER JNGAM
      INTEGER JNJVSP
      INTEGER JNLESP
      INTEGER JNTRJG


      PARAMETER ( JNMDCT =  1001)
      PARAMETER ( JNTRJX =  1002)
      PARAMETER ( JNTRJV =  1003)
      PARAMETER ( JNTRJE =  1004)
      PARAMETER ( JNTOPO =  1005)
      PARAMETER ( JNXVI  =  1006)
      PARAMETER ( JNREST =  1007)
      PARAMETER ( JNRSTA =  1008)
      PARAMETER ( JNRSTP =  1009)
      PARAMETER ( JNPERT = 1010)
      PARAMETER ( JNDIHE = 1011)
      PARAMETER ( JNXVE  = 1012)
      PARAMETER ( JN4COR = 1013)
      PARAMETER ( JN4NDX = 1014)
      PARAMETER ( JNCNV  = 1015)
      PARAMETER ( JNGAM  = 1016)
      PARAMETER ( JNJVSP = 1017)
      PARAMETER ( JNLESP = 1018)
      PARAMETER ( JNTRJG = 1019)

C add other types values here.
C remember to change NLSRES accordingly!

      INTEGER NLSRES
      PARAMETER (NLSRES = JNTRJG)

      PARAMETER ( JND00  = NLSRES + 1)
      PARAMETER ( JND01  = NLSRES + 2)
      PARAMETER ( JND02  = NLSRES + 3)
      PARAMETER ( JND03  = NLSRES + 4)
      PARAMETER ( JND04  = NLSRES + 5)
      PARAMETER ( JND05  = NLSRES + 6)
      PARAMETER ( JND06  = NLSRES + 7)
      PARAMETER ( JND07  = NLSRES + 8)
      PARAMETER ( JND08  = NLSRES + 9)
      PARAMETER ( JND09  = NLSRES + 10)

C     the number of reserved units
      INTEGER MAXNTS
      PARAMETER (MAXNTS = JND09)


      
C     The arrays of these common blocks are initialized
C     in a blockdata statement in fileio.f

      INTEGER MAXFNM
      PARAMETER (MAXFNM = 8)

      CHARACTER*(MAXFNM) UNAME
      COMMON /UNCHAR/UNAME(MAXNTS)

      INTEGER IUNUM
      COMMON /UNINT/IUNUM(MAXNTS)
