Cmb Revised 21 August 2006
C GROMOS Molecular Dynamics, Stochastic Dynamics and Energy Minimisation
C main program

COMMPROG PROMD
C     PROGRAM PROMD performs an MD, SD or EM run for a number of
C     solute (protein) molecules plus solvent molecules in three or
C     four spatial dimensions.
C     The options possible, depending on the number of dimensions,
C     are found in the following table, where footnotes are marked with
C     (XX), e.g. (1) marks the first footnote.
COMMVERB
C-----------------------------------
C     Boundary Conditions    3D   4D
C
C     none (in vacuo run)   YES  YES
C     truncated octahedron  YES  YES
C     rectangular           YES  YES
C     monoclinic            YES   NO
C-----------------------------------
C     Coupling
C     Temperature coupling  YES  YES (1)
C     Virial calculation    YES   NO (2)
C     Pressure calculation  YES   NO (2)
C     Pressure coupling     YES   NO (2)
C-----------------------------------
C     Restraining
C
C     position restraining  YES  YES
C     position constraining YES  YES (3) 
C     distance restraining  YES  YES
C     dihedral restraining  YES  YES
C
C     Bond Constr. (SHAKE)  YES  YES (4)
C
C     The following subroutines can be used in 4 dimensions, however,
C     their effect is restricted to three dimensions. E.g. SHIAG
C     only puts the atoms into the box in 3 dimensions, the fourth dimension
C     is unchanged:
C     SHIAG, CENMAS
C     In addition, certain "very physical" forces or quantities are
C     not generalized to 4 dimensions:
C     Reaction field calulations
C     Friction calculations in L<FRIC> are only performed in 3D.
C
C     FOOTNOTES:
C     ( 1) Temperature coupling in four dimensions is to:
C          a: one, two or three baths in three dimensions
C          b: a separate bath for the whole (solute and solvent)
C          system in the fourth dimension
C     ( 2) The virial cannot be calculated in 4D, as a consequence
C          neither can the pressure be calculated, and so pressure
C          coupling is not possible in 4D.
C     ( 3) Although it is possible to use position constraining while
C          performing a simulation in 4D, a 4D atom may not be positionally
C          constrained.
C     ( 4) Performing a 4D simulation, if the user selects that bonds are
C          to be SHAKEN in 4D (using the NT4D array), then only those bonds
C          between atoms which are both in 4D are also shaken in 4D.
C          I.e. if one or both atoms of a bond is in 3D, then that bond
C          is SHAKEN in 3D.
COMMEND
C
C     A L<RLAM> dependent perturbation potential can be applied in
C     order to determine the free energy difference between a
C     state A and a state B of the system.
C     A L<MU> dependent perturbation can be applied in order to determine
C     the free energy difference between the system in three dimensions
C     and in four dimensions.
C
C     All common blocks are defined in include files which are organized
C     as follows:
COMMVERB
C     the topology                  toposz.h topoar.h
C     the coordinates (X,F, etc)    coordsz.h coordar.h
C     perturbation data             pertsz.h pertar.h
COMMEND
C     In all files, the PARAMETER statements defining the sizes of
C     arrays are in header files ending in *sz.h.
C     The actual arrays are defined in header files ending in *ar.h.
C     Thus, if one wishes to simulate a system which has more
C     atoms than defined in coordsz.h, then only that file
C     needs to be modified and GROMOS recompiled using the makefile
C     provided.
C
C     PROMD reads data from a number of fortran units (files)
C     depending on the values read in from unit 5 (standard input).
C
C     The following units are defined in the file units.h
COMMVERB
C--------------------------------------
C
C     Input units
C     -----------
C     MD run control file            IOMDCT (5) (standard input)
C
C     topology file                  IOTOPO (20)
C
C     data for a continuation run:
C     initial coordinates, velocities
C     initial position restraint positions
C     initial box size,
C     initial stochastic integral
C     values and random generator
C     seed. (for SD)
C     initial values for RLAM,DLSUM,
C     RMU,DMSUM
C     time-averaged distances for
C     restraining
C     time-averaged cos square and
C     cos for for J-value
C     restraining                    IOXVI (21)
C
C     friction coefficients for SD   IOGAM (29)
C
C     specification of dihedrals to
C     be used in J-value restraining IOJVSP (26)
C     
C     reference coordinates for
C     position re(con)straining      IOREST (22)
C
C     sequence numbers specifying
C     position re(con)strained atoms IORSTA (23)
C
C     distance restrained
C     atom pair specifications       IORSTP (24)
C
C     data determining a
C     perturbation                   IOPERT (30)
C
C     specification of restrained
C     dihedrals                      IODIHE (25)
C
C     local elevation dihedral
C     specifications                 IOLESP (27)
C
C     sequence numbers specifying    IO4NDX (28)
C     atoms in 4D                    
C--------------------------------------
C
C     Output units
C     ------------
C     promd output (printing)        ISTDOT (6 standard output)
C
C     coordinates trajectory         IOTRJX (12)
C
C     velocities trajectory          IOTRJV (13)
C     average restrained
C
C     energies trajectory            IOTRJE (15)
C
C     free energies trajectory       IOTRJG (16)
C
C     data allowing a subsequent continuation run
C     This data is typically read in on unit IOXVI (21)
C     by subsequent PROMD runs.
C     final coordinates, velocities,
C     reference positions, box size,
C     position restraint positions,
C     stochastic integral values,
C     time-averaged distances
C     time-averaged cos squared and
C     cos for j-value restraining
C                                    IOXVE (11)
COMMEND
C     L<PROMD> first prints out its RCS version number and a
C     welcome message, then calls L<RDTOPO> to read in the topology
C     file of the system to be simulated.
C
C     L<RDMD> is then called in order to read in the user provided
C     control data.
C     L<PRMD> is called in order to write the data
C     read by L<RDMD> to screen with comments 
C     about the type of simulation L<PROMD> is going to perform.
C
C     After checking that the system to be simulated will fit
C     into the arrays provided, and some checking of
C     data consistency between topology and PROMD control file
C     (L<CHKNSP> and L<CHKNRE>), additional files are can be read in
C     depending on the switches read in by L<RDMD>:
C
COMMLIST
CITM  L<GTDRST>
C        reads in distance restraint data
CITM L<GTCOOR>
C        reads in or generates initial coordinates L<X>, initial
C        velocities L<V>, the simulation box L<BOX>, position restraint
C        reference positions L<XC>, and L<RIIAVE> time averaged
C        distances for distance restraining.
C
C        Depending on the input switches defined, L<GTCOOR> calls
C        L<RDCFG> to read some of the quantities from file, and generates
C        others. Subroutine L<GAUSS> is called in order to produce
C        the initial velocities. If the simulation is performed in
C        4D, L<GT4NDX> is called in order to read in the L<C4D> array.
C        Subroutine L<RDXRST> can be called which reads in L<CXC>
C        (the atomic B factors) and/or the position restraint
C        reference positions L<XC>.
CITM L<GTDHRS>
C        reads in restrained dihedrals data
CITM L<RDPERT>
C        reads in data pertaining to a perturbation. L<PRPERT>
C        which prints the data to screen is also called in this case.
COMMEND
C
C     Subroutine L<SETCON> can be called in order to set up the
C     distance constraint data used by L<SHAKE> using
C     bond data in the topology.
C
C     If a minimisation is to be performed, L<RUNEM> is called
C     in other cases, i.e. for a molecular dynamics or stochastic
C     dynamics run, L<DOMDSD> is called which in turn calls L<RUNMD>.
C
C     The final configuration is written to file using L<WRCFG>.
COMMEND
COMMVERB
C     PROMD reads its control file from IOXVI (5) which consists of a
C     subset of the following blocks in alphabetic order.
C     On file they occur in logical order. The sequence number of the
C     logical order is given below in parentheses after the blockname
C     for each block.
CCC
C     BOUNDARY (7)
C     NTB, BOX(1), BOX(2), BOX(3), BETA, NRDBOX
C
C     NTB -2..2 define boundary conditions
C     <0 truncated octahedron periodic boundary conditions
C     >0 rectangular/monoclinic periodic boundary conditions
C     0 vacuum
C     abs(NTB) = 2 the virial (pressure) is calculated
C     BOX(1..3) > 0.0 periodic box dimensions
C     BETA > 0.0 angle of box (inclination between x- and z- axes)
C     in degrees
C     BETA must be 90.0 degrees if abs(NTB) = 2 or if NTB<0
C     NRDBOX 0,1 controls reading of BOX dimensions
C     0 use the box dimensions from the BOUNDARY block
C     1 read box dimensions from startup file (BOX) IOXVI (21)
C
CCC
C     CENTREOFMASS (11)
C     NDFMIN, NTCM, NSCM
C
C     NDFMIN >= 0 number of degrees of freedom to subtract from total
C     when calculating the system temperature
C     NTCM 0,1 controls initial centre of mass motion removal
C     0 no initial centre of mass motion removal
C     1 inital centre of mass motion is removed
C     NSCM >= 0 controls centre of mass motion removal during
C     simulation
C     0 no centre of mass motion removal
C     > 0 centre of mass motion removal every NSCM steps
C
CCC
C     DIHEREST (20)
C     NTDLR, CDLR
C
C     NTDLR 0,1,2 controls dihedral angle restraining
C     0 no dihedral restraining
C     1 restraining using force constant CDLR
C     2 restraining using force constant CDLR * CPLR
C     CDLR > 0.0 force constant in dihedral angle restraining function
C     specification of dihedrals to be restrained (DIHRESSPEC) is 
C     read from IODIHE (25)
CCC
C     DISTREST (19)
C     NTDR, CDIS, DR0, TAUDR, NRDDR
C
C     NTDR -2..2 controls distance restraining
C     0 no distance restraining
C     abs(NTDR)=1 restraining using force constant CDIS
C     abs(NTDR)=2 restraining using force constant CDIS * W0
C     NTDR<0 time-averaged restraining
C     NTDR>0 instantaneous restraining
C     CDIS >= 0.0 distance restraining force constant
C     DR0 >=0.0 distance offset in restraining function
C     TAUDR >= 0.0 coupling time for time averaging
C     NRDDR 0,1 controls reading of averages from startup file
C     0 do not read averages from IOXVI (21)
C     1 read averages from IOXVI (21)
C     specification of distance restraints (DISRESSPEC) is read
C     from IORSTP (24)
C
CCC
C     FORCE (15)
C     NTF array
C     NEGR NRE(1..NEGR)
C
C     NTF(1..10) 0,1 determines terms used in force calculation
C     NTF(I)=0 do not use term I
C     NTF(I)=1 use term I
C     NTF(1) bonds involving hydrogens
C     NTF(2) bonds not involving hydrogens
C     NTF(3) bond angles involving hydrogens
C     NTF(4) bond angles not involving hydrogens
C     NTF(5) improper dihedrals involving hydrogens
C     NTF(6) improper dihedrals not involving hydrogens
C     NTF(7) dihedrals involving hydrogens
C     NTF(8) dihedrals not involving hydrogens
C     NTF(9) nonbonded charge interactions
C     NTF(10) nonbonded interactions
C     NEGR >= 0 number of energy groups
C     NRE > 0 last atoms in each energy group
C     NRE values must be in ascending order
C
CCC
C     FOURDIM (23)
C     NT4DIM, CW4DA, TEMP4I, NDFMI4, NT4XI, NT4XO
C     NTT4, TEMP04, TAUT4
C     NTCW4D, CW4DB, TEMP0B
C     NTF(1..6)
C
C     NT4DIM 0,1,2 controls dimensionality of the simulation
C     0 perform the simulation in 3D
C     1 perform the simulation in 4D using CW4DA
C     2 perform the simulation in 4D using CW4DA * C4D
C     CW4DA >= 0.0 4thD restraining force constant
C     TEMP4I >= 0.0 initial 4thD temperature
C     NDFMI4 >= 0 degrees of freedom to remove from the 4th dimension
C     NT4XI 1..4 controls reading/generation of initial 4thD coordinates
C     (X4) and velocities (V4)
C     1 X4 set to zero, Maxwell distribution for V4
C     2 Boltzmann distribution for X4, Maxwell distribution for V4
C     3 read X4 from IOXVI (21), Maxwell distribution for V4
C     4 read X4 and V4 from IOXVI (21)
C     NT4XO 0,1 controls writing of 4thD data into trajectory files
C     0 write 3D coordinate and velocity trajectories
C     1 write 4D coordinate and velocity trajectories
C
C     NTT4 0,1 controls temperature coupling in the 4th dimension
C     TEMP04 >= 0.0 bath reference temperature
C     TAUT4 >= 0.0 coupling time
C     
C     NTCW4D 0,1 controls scaling of TEMP04 and CW4DA
C     0 constant reference temperature using TEMP04 and CW4DA
C     1 reference temperature is linearly changed from TEMP04 to TEMP0B and
C     force constant is changed from CW4DA to CW4DB over NSTLIM steps
C     CW4DB >= 0.0 final 4thD restraining force constant
C     TEMP0B >= 0.0 final bath reference temperature
C
C     NTF4(1..6) 0,1 controls dimensionality of force field terms
C     NTF4(I)=0 calculate I in 3D
C     NTF4(I)=1 calculate I in 4D
C     NTF4(1) bond interactions or SHAKE if specified
C     NTF4(2) bond angles
C     NTF4(3) dihedrals
C     NTF4(4) non-bonded interactions
C     NTF4(5) distance restraining
C     NTF4(6)=0 no 4thD harmonic restraining
C     NTF4(6)=1 4thD harmonic restraining
C     specification of 4D atoms (FOURDIMATOMSPEC) is read from IO4NDX (28)
C
CCC
C     J-VAL (21)
C     NTJR, NTJRH, CJR, TAUJR, NRDJR
C
C     NTJR -2..2 controls J-value restraining
C     0 no J-value restraining
C     abs(NTJR)=1 restraining using force constant 
C     CJR *(A cos^2 + B cos + C)
C     abs(NTJR)=2 restraining using force constant
C     CJR * CPJR(A cos^2 + B cos + C)
C     NTJR <0 use time-averaged distances in restraining
C     NTJR >0 use instantaneous distances in restraining
C     NTJRH 0,1 determines function used in restraining
C     0 use a half harmonic function
C     1 use a harmonic function
C     CJR >= 0.0 J-value restraining force constant
C     TAUJR >= 0.0 coupling time for time averaging
C     NRDJR 0,1 controls reading of averages from startup file IOXVI (21)
C     0 do not read averages from startup file IOXVI (21)
C     1 read averages from IOXVI (21)
C     specification of dihedrals to be restrained and A,B,C (JVALRESSPEC)
C     is read from IOJVSP (26)
C
CCC
C     LOCALELEVATION (22)
C     NTLE,CWLE,NRDLE
C
C     NTLE 0,1,2 controls local elevation search
C     0 no local elevation search
C     1 local elevation search using a Gaussian function
C     2 local elevation search using an inverse quadratic function
C     CWLE >= 0.0 local elevation force constant
C     NRDLE 0,1 controls reading of local elevation memory
C     0 do not read local elevation memory
C     0 read local elevation memory from startup file IOXVI (21)
C     specification of local elevation dihedrals (LEDIHSPEC) is read
C     from IOLESP (27)
C
CCC
C     LONGRANGE (17)
C     EPSRF, APPAK, RCRF
C
C     EPSRF =0.0,or >= 1.0 eps used in reaction-field calculation
C     APPAK >= 0.0 kappa used in reaction-field calculation
C     RCRF > 0.0 cut-off radius used in reaction-field calculation
C
CCC
C     MINIMISE (2)
C     NTEM, NCYC, DELE, DX0,DXM
C
C     NTEM 1,2 controls minimisation method
C     1 use steepest descent minimising
C     2 use conjugate gradient (CG) minimisation
C     NCYC > 0 resetting of search direction (CG)
C     DELE >= 0.0 energy threshold for stopping
C     DX0 >= 0.0 initial step size
C     DXM >= 0.0 maximum step size
C     the existence of this block on file means that PROMD will
C     perform an energy minimisation run
C
CCC   
C     PATHINT (25)
C     NTPI
C
C     NTPI 0,1 controls path integral calculation
C     0 no path integral calculation
C     1 perform path integral calculation
C     for path integral calculations, a path integral topology must
C     be specified
C
CCC   
C     PCOUPLE (10)
C     NTP, PRES0, COMP, TAUP
C
C     NTP 0,1,2 controls pressure coupling
C     0 no pressure coupling
C     1 isotropic pressure scaling
C     2 anisotropic pressure scaling (x-, y-, z- axes)
C     PRES0 > 0.0 reference pressure
C     COMP > 0.0 isothermal compressibility
C     TAUP > 0.0 coupling time
C     pressure coupling requires calculation of the virial 
C     (see BOUNDARY)
C
CCC
C     PERTURB (24)
C     NTG, NRDGL, RLAM, DLAMT, RMU, DMUT
C     ALPHLJ, ALPHC, NLAM, MMU
C     
C     NTG 0,1,2,3 controls free energy calculation
C     0 no free energy calculation is performed
C     1 calculate dVdLAMDA
C     2 calculate dVdMU
C     3 calculate both dVdLAMDA and dVdMU
C     NRDGL 0,1 controls reading of initial values
C     0 use values for LAMDA, MU given in IOMDCT (5)
C     1 read initial values from startup file IOXVI (21)
C     RLAM 0.0..1.0 (initial) value for LAMDA
C     DLAMT >= 0.0 LAMDA increase size per time unit
C     RMU 0.0..1.0 (initial) value for MU
C     DMUT >= 0.0 MU increase size per time unit
C     ALPHLJ >= 0.0 Lennard Jones soft core parameter
C     ALPHC >= 0.0 Coulomb soft core parameter
C     NLAM > 0 power dependence of LAMDA
C     MMU > 0 power dependence of MU
C     a perturbation topology file is read from IOPERT (30)
C
CCC
C     PLIST (16)
C     NTNB, NSNB, RCUTP, RCUTL
C
C     NTNB 0,1 controls pairlist construction before the first step
C     0 no pairlist is constructed before the first step
C     1 a pairlist is constructed before the first step
C     NSNB > 0 controls frequency (number of steps) a pairlist is
C     constructed
C     RCUTP > 0.0 cut-off used in pairlist construction
C     RCUTL > 0.0 cut-off used in long range interaction
C     
CCC
C     POSREST (18)
C     NTR, CHO, NRDRX
C
C     NTR 0,1,2,3 controls atom position re/constraining
C     0 no position re/constraining
C     1 restraining using force constant CHO
C     2 restraining using force constant CHO/atomic B factors
C     3 position constraining
C     position constraining is only possible with NPM = 1
C     CHO >= 0.0 position restraining force constant
C     NRDRX 0,1 controls reading of position re/constraint reference
C     positions
C     0 read reference coordinates and B factors from IOREST (22)
C     1 read reference coordinates from startup file IOXVI (21)
C     specification of the atoms to re/constrain (POSRESSPEC) is read
C     from IORSTA (23)
C     reference coordinates are in REFPOSITION blocks
C     B factors are in BFACTOR or BFACTORANISO blocks
C
CCC   
C     PRINT (12)
C     NTPR, NTPL, NTPP
C
C     NTPR >= 0 controls printing of energies to ISTDOT (6)
C     0 no printing out of energies
C     > 0 print out energies every NTPR steps
C     NTPL >= 0 controls centre of mass (c.o.m.) printing and
C     total energy least squares fitting
C     0 no printing out of c.o.m and total energy time line
C     > 0 print out c.o.m. and fit of total energy every NTPL steps
C     NTPP 0,1 controls dihedral angle transition monitoring
C     0 no dihedral angle transition monitoring
C     1 perform dihedral angle transition monitoring
C
CCC   
C     SHAKE (14)
C     NTC, TOL
C
C     NTC 1,2,3 controls use of SHAKE for the solute
C     1 no shake is performed (solute)
C     2 SHAKE the bonds involving hydrogens (solute)
C     3 SHAKE all bonds (solute)
C     TOL > 0.0 SHAKE tolerance (relative geometric precision), for
C     solute and solvent
C
CCC
C     START (5)
C     NTX, INIT, IG, TEMPI, HEAT, NTXO, BOLTZ
C
C     NTX 1..3 reading/generation of initial 3D coordinates (X),
C     3D velocities (V), and stochastic integrals (SX)
C     action                                  1    2    3  
C     read X from IOXVI (21)                  yes  yes  yes
C     set V to zero if TEMPI = 0.0            yes  no   no 
C     read V from IOXVI (21) if  TEMPI = 0.0  no   yes  yes
C     read SX from IOXVI (21)                 no   no   yes
C
C     INIT 1..4 startup configuration
C     action          1      2    3    4
C     shake X        yes     no   no   no
C     shake V        yes     yes  no   no
C     C.O.M removal  yes     yes  yes  no
C     if NTCM =1
C
C     IG any random number generator seed
C     TEMPI >= 0.0 initial temperature for velocities
C     HEAT any initial temperature scaling
C     > 0.0 scale all solute and solvent (3 and 4D) velocities by HEAT
C     < 0.0 scale all solute (3 and 4D) velocities by abs(HEAT)
C     NTXO 0,1 format of final configuration
C     0 binary
C     1 formatted
C     BOLTZ > 0.0 the Boltzmann constant
C     3D coordinates are in POSITION or POSITIONRED blocks
C     3D velocities are in VELOCITY or VELOCITYRED blocks
C     stochastic integrals are in STOCHINT blocks
C     the final configuration is written to IOXVE (11)
C
CCC
C     STEP (6)
C     NSTLIM, T, DT
C
C     NSTLIM > 0 number of steps to perform in the minimization or
C     simulation
C     T >= 0.0 time at beginning of simulation (MD and SD)
C     DT >= 0.0 time step in MD or SD simulations
C
CCC
C     STOCHASTIC (3)
C     NTFR, NSFR, NBREF, RCUTF, CFRIC, TEMPSD
C
C     NTFR 0..3 define atomic friction coefficients GAM
C     0 set GAM to 0.0
C     1 set GAM to CFRIC
C     2 set GAM to CFRIC*GAM0, where GAM0 is read from IOGAM (29)
C     3 calculate GAM using SUBR. FRIC
C
C     NSFR > 0 recalculate GAM every NSFR steps
C     NBREF > 0 number of neighbour atoms within
C     RCUTF distance for an atom to be considered buried
C     RCUTF >= 0.0 inter atom distance to consider when
C     calculating GAM
C     CFRIC >= 0.0 global weighting for GAM
C     TEMPSD >= 0.0 temperature of stochastic bath
C     NSFR, NBREF and RCUTF only have an effect if NTFR = 3
C     atomic friction coefficients are in FRICTIONSPEC blocks
C     the existence of this block on file means that PROMD will
C     perform an SD run
C
CCC
C     SUBMOLECULES (8)
C     NSPM, NSP(1.. NSPM)
C
C     NSPM > 0 number of (sub)molecules the solute consists of
C     NSP(I) > 0 atom sequence number in topology of last atom in
C     submolecule I
C
CCC
C     SYSTEM (4)
C     NPM, NSM
C
C     NPM >= 0 number of (identical) solute molecules
C     NSM >= 0 number of (identical) solvent molecules
C
CCC
C     TCOUPLE (9)
C     NTT(1)   TEMP0(1) TAUT(1)
C     NTT(2)   TEMP0(2) TAUT(2)
C     NTT(3)   TEMP0(3) TAUT(3)
C     
C     NTT(1) -3..3 controls temperature coupling of solute
C     internal and rotational degrees of freedom ( d.o.f.) in three dimensions
C     NTT(2) -3..3 controls temperature coupling of solute 
C     centre of mass translational d.o.f. in three dimensions
C     NTT(3) -3..3 controls temperature coupling of solvent d.o.f.
C     in three dimensions
C     0 no temperature coupling for set of d.o.f.
C     1 couple one set of d.o.f. to one bath
C     2,-2 couple two sets of d.o.f. to one bath (+ sign)
C     -3,3 couple three sets of d.o.f. to one bath (+ sign)
C     TEMP0 >= 0.0 bath reference temperature
C     TAUT >= 0.0 coupling time 
C
CCC
C     TITLE (1)
C     text
C
C     text consists of up to ten lines of arbitray text. A line
C     may not be longer than 80 characters in length
C
CCC
C     WRITE (13)
C     NTWX, NTWSE, NTWV, NTWE, NTWG, NTPW
C     
C     NTWX >= 0 controls writing of coordinate trajectory
C     0 no coordinate trajectory is written
C     >0 write solute and solvent coordinates every NTWX steps
C     <0 write solute coordinates every abs(NTWX) steps
C     NTWSE >= 0 selection criteria for trajectory writing
C     0 write normal coordinate trajectory
C     >0 write minimum energy trajectory
C
C     NTWV >= 0 controls writing of velocity trajectory
C     0 no velocity trajectory is written
C     >0 write solute and solvent velocities every NTWV steps
C     <0 write solute velocities every abs(NTWV) steps
C
C     NTWE >= 0 controls writing of energy trajectory
C     0 no energy trajectory is written
C     >0 write energy variables every NTWE steps
C
C     NTWG >= 0 controls writing of free energy trajectory
C     0 no free energy trajectory is written
C     >0 write free energy variables every NTWG steps
C     NTPW >= 0 controls format of trajectory files
C     0 write binary
C     1 write formatted
C     coordinates are written to IOTRJX (12) in POSRED blocks
C     velocities are written to IOTRJV (13) in VELRED blocks
C     energies are written to IOTRJE (15) or IOTRJX (12) in ENERGY blocks
C     volume and pressure quantities are written to IOTRJE (15) in
C     VOLUMEPRESSURE blocks
C     free energy quantities are written to IOTRJG (16) in FREEENERGY3D4
C     and FREEENERGYLAMBDA blocks
C
COMMEND
C     ==================================================================
      SUBROUTINE mm_setup(NRPT,NSOLV)
C     ==--------------------------------------------------------------==
C includes
C topology
      INCLUDE 'toposz.h'
      INCLUDE 'topoar.h'
C coords
      INCLUDE 'coordsz.h'
      INCLUDE 'coordar.h'
C force and energy arrays
      INCLUDE 'forcesz.h'
C the box data
      INCLUDE 'box.h'
C promd control file data
      INCLUDE 'md.h'
C unit mappings
      INCLUDE 'units.h'
c
      include '../system.h'
      INCLUDE '../parac.inc'
      include '../mm_input.inc'
      include '../fileopen.inc'

      character*120 fileop

C local parameters
C     standard output unit
      INTEGER ISTOUT
      PARAMETER (ISTOUT = 6)
      INTEGER MAXFLD
      PARAMETER (MAXFLD = 5)
      INTEGER MAXSTR
      PARAMETER (MAXSTR = 2*MAXFLD+1)
C local vars
      INTEGER IFNCFG,K,KFRST,KLST
      LOGICAL LFORMF,FERROR
      CHARACTER PRGSTR*(5)
      CHARACTER STR*(MAXSTR)
C
C     NATTOT total number of atoms
C      NSOLV number of solvent coordinates
C       NRPT number of solute atoms
C      NUM4DA number of atoms which are in 4D

      INTEGER NRPT,NSOLV,NATTOT,NUM4DA

C
C     DLSUM = cumulative value of the derivative of the perturbation
C             energy with respect to LAMBDA
C     DMSUM = cumulative value of the derivative of the perturbation
C             energy with respect to MU
      real*8 DLSUM(MXETBL),DMSUM(MXETBL)
      character(90)::my_char
C data statements
      DATA NUM4DA/0/
      DATA (DLSUM(K), K=1,MXETBL)/MXETBL* 0.0D0/
      DATA (DMSUM(K), K=1,MXETBL)/MXETBL* 0.0D0/

      DATA PRGSTR /'PROMD'/
#include "inc_promd4.h"
C begin PROMD
      call schifezza
      call schifezza1
      CALL FILEOPEN(ISTDIN,FILEGRIN,FO_OLD,FERROR)
      IF(FERROR) GOTO 114
      CALL FILEOPEN(IOTOPO,FILETOPO,FO_OLD,FERROR)
      IF(FERROR) GOTO 115
      CALL FILEOPEN(IOXVI,FILECONF,FO_OLD,FERROR)
      IF(FERROR) GOTO 116

#ifdef TIMING
      CALL INITIM()
      CALL PRNTIM()
      WRITE(6,*)'BEGIN PROMD'
#endif

      WRITE(6,*)'promd.f,v 1.55 1996/11/25 14:49:50 wscott Exp'
      WRITE(6,*)'Present revision: Tsukuba, 13 July 2005'
      WRITE(6,*)
      WRITE(6,*)'PROGRAM PROMD PERFORMS EITHER'
      WRITE(6,*)'  AN ENERGY MINIMISATION,'
      WRITE(6,*)'  A MOLECULAR DYNAMICS RUN,'
      WRITE(6,*)'  OR A STOCHASTIC DYNAMICS RUN'
      WRITE(6,*)'IN THREE OR FOUR DIMENSIONS FOR A NUMBER OF SOLUTE'
      WRITE(6,*)'(PROTEIN) MOLECULES PLUS SOLVENT MOLECULES, USING:'
      WRITE(6,*)
      WRITE(6,*)' A: SPECIFICATION DATA'
      WRITE(6,*)'-----------------------'
      WRITE(6,1000)
      WRITE(6,1100)
      WRITE(6,1200)
      WRITE(6,1300)
      WRITE(6,1400)
      WRITE(6,1500)
      WRITE(6,1600)
      WRITE(6,1700)
      WRITE(6,1800)
      WRITE(6,*)
      WRITE(6,*)' B: STARTUP DATA'
      WRITE(6,*)'-----------------'
      WRITE(6,1900)
      WRITE(6,2000)
      WRITE(6,2100)
      WRITE(6,*)

C*****READ MOLECULAR TOPOLOGY FILE, PRINT TITLE
      CALL RDTOPO
      WRITE(6,1010)
      CALL PRTIT('TOPOLOGY TITLE',NTPLNS,TOPTIT)
      CALL CHKCG(PRGSTR)

C*****READ DATA CHARACTERISING THE MD-RUN
      WRITE(6,1110)
      CALL RDMD
      CALL PRMD
      CALL GLBCHK(PRGSTR)

C check consistency of control file and topology file
C NRP = number of solute atoms per molecule
C NPM = number of solute (protein) molecules
C NRAM= number of solvent atoms per molecule
C NSM = number of solvent molecules
C NRPT = total number of solute atoms
 
      NRPT = NRP*NPM
      NSOLV= NRAM*NSM
      NATTOT= NRPT + NSOLV
cmb
      WRITE(6,'('' NRP    ='',I12,4x,''NPM ='',I6)')NRP,NPM
      WRITE(6,'('' NRPT   = NRP*NPM ='',I12)')NRPT
      WRITE(6,'('' NRAM   ='',I6,4x,''NSM ='',I12)')NRAM,NSM
      WRITE(6,'('' NSOLV  = NRAM*NSM ='',I12)')NSOLV
      WRITE(6,'('' NATTOT = NRPT + NSOLV ='',I12)')NATTOT
cmb
C check NRE(...) and NRP
      CALL CHKNRE(PRGSTR,NATTOT)

C check NSP
C     If the SUBMOLECULES input block was not on the input file,
C     (this is legal for energy minimising) NSPM will be zero.
C     However, if it is on file, it must be correct.. so check it.
      IF (NSPM .NE. 0 .OR. IDOPRO .NE. IDOEM) THEN
         CALL CHKNSP(PRGSTR,NDIM,NPM,NRP,NCAG,INC)
      ENDIF

C     check dihedral angle monitoring if applicable
      IF (NTPP .NE. NTPPNO .AND. MXDMON .LT. NPM*NPHI) THEN
         WRITE(6,*)PRGSTR,': too many dihedrals to monitor!'
         WRITE(6,550)NPM*NPHI
         WRITE(6,560)MXDMON
         CALL STOPGM('MM_SETUP','FATAL ERROR')
      ENDIF

CCCCC
C     Here follows the reading of any optional data files
CCCCC
      WRITE(6,1120)

C
C*****READ AND PRINT DATA DETERMINING THE PERTURBATION, WHEN REQUIRED
C do this BEFORE GTCOOR because we check 4D perturbation consistency
C there
      IF (NTG .NE. NTGOFF) THEN
         WRITE(6,1210)
!        CALL RDPERT
!        CALL PRPERT
C check consistency between topology and perturbation
!        CALL CHKPRT(PRGSTR)
      ENDIF


C*****READ DISTANCE RE(CON)STRAINT SPECIFICATION
      IF (NTR .NE. NTROFF) THEN
         WRITE(6,1310)
C check implementation restriction: constraining only possible
C with NPM = 1 (because of the modification of WMAS and WINV arrays)
         IF (NTR .EQ. NTRCON .AND. NPM .GT. 1) THEN
            WRITE(6,*)PRGSTR,': IMPLEMENTATION RESTRICTION!!'
            WRITE(6,*)'NPM must be 1 if position constraining',
     $           ' is specified'
            WRITE(6,*)'(This does not apply to position restraining)'
            CALL STOPGM('MM_SETUP','FATAL ERROR')
         ENDIF

C     now the indices of the position re(con)strained atoms
         IF (NTR .EQ. NTRCON) THEN
C     position constraining
            CALL GTIRST(NRP,MAXXC,JRC,NRCON)

C check that we dont perturb the mass of a positionally constrained atom
            IF (NTG .NE. NTGOFF) THEN
               CALL CHKPCN(PRGSTR,NRCON,JRC)
            ENDIF
         ELSE
C     position restraining
            CALL GTIRST(NATTOT,MAXXC,JRC,NRRST)
         ENDIF
      ENDIF

C*****READ DISTANCE RESTRAINT SPECIFICATION, PRINT TITLE, WHEN REQUIRED
C     do this BEFORE we call GTCOOR as GTCOOR also gets
C     distance restraint data for continuation run using
C     time averaged dist. rest.
      NDR = 0
      IF (NTDR .NE. NTDROF) THEN
         WRITE(6,1410)
         CALL GTDRST(NRPT)
         IF (ABS(NTDR) .EQ. NTDRCD) THEN
C     set W0 to CDIS
!$OMP parallel do private(K)
            DO K =1,NDR
               W0(K) = CDIS
            ENDDO
         ELSEIF (ABS(NTDR) .EQ. NTDRW0) THEN
C     scale W0 by CDIS
!$OMP parallel do private(K)
            DO K =1,NDR
               W0(K) = W0(K)*CDIS
            ENDDO
         ENDIF
      ENDIF

C now check we have enough space in RIIAVE for distance restraints
      IF (NPM*NDR .GT. MAXTDR) THEN
         WRITE(6,*)PRGSTR,'too many distance restraints specified!'
         WRITE(6,*)'Number of distance restraints per molecule:',NDR
         WRITE(6,*)'Number of (solute) molecules              :',NPM
         WRITE(6,*)'Total number of distance restraints       :',NPM*NDR
         WRITE(6,*)'Maximum allowed:                    MAXTDR=',MAXTDR
      ENDIF


C*****READ RESTRAINED DIHEDRALS, PRINT TITLE, WHEN REQUIRED
      IF (NTDLR .NE. NTDLRN) THEN
         CALL RDIHRS(PRGSTR,NRPT)
      ENDIF
      
C READ J-VALUE RESTRAINING SPECIFICATION
C     do this BEFORE we call GTCOOR as GTCOOR can also get
C     j-val data data for continuation run.
      NDJV = 0
      IF (NTJR .NE. NTJROF) THEN
         WRITE(6,1610)
         CALL GTJVSP(NRP)
      ENDIF


C READ LOCAL ELEVATION SPECIFICATION
C do this BEFORE GTCOOR because we could be getting LE data there
C for a continuation run
      IF (NTLE .NE. NTLEOF) THEN
         WRITE(6,1710)
         CALL GTLESP(NRP)
      ENDIF


C READ 4D ATOM INDICATOR LIST
      NUM4DA = 0
      IF (NT4DIM .NE. NT4OFF) THEN
         WRITE(6,1819)
         CALL GT4NDX(NATTOT,C4D)
C     count the number of atoms in 4D
         DO K=1,NATTOT
            IF (C4D(K) .GE. 0.0D0) THEN
               NUM4DA = NUM4DA + 1
            ENDIF
         ENDDO
         WRITE(6,580)NUM4DA
         WRITE(6,*)
      ENDIF

CCCCC
C     End of optional files
CCCCC

C*****READ COORDINATES AND VELOCITIES, PRINT TITLE,
C     ALSO READ IN IF NEEDED:
C     BOX DIMENSIONS
C     REFERENCE POSITIONS FOR RE(CON)STRAINING
C     TIME AVERAGED DIST RESTRAINING DATA
C     TIME AVERAGED J-VAL RESTRAINING DATA
C     LOCAL ELEVATION DATA
C     DLSUM,DMSUM

      WRITE(6,1130)
C check max number of atoms
      IF (NATTOT .GT. MAXNAT) THEN
         WRITE(6,*)PRGSTR,': too many atoms specified!'
         WRITE(6,500)NRPT
         WRITE(6,501)NSOLV
         WRITE(6,502)NATTOT
         WRITE(6,505)MAXNAT
         WRITE(6,*)
         CALL STOPGM('MM_SETUP','FATAL ERROR')
      ENDIF

C check coordinate array sizes
      IF (NDIM*NATTOT .GT. MAXXCO) THEN
         WRITE(6,*)PRGSTR,': too many atom coordinates specified!'
         WRITE(6,500)NRPT
         WRITE(6,501)NSOLV
         WRITE(6,502)NATTOT
         WRITE(6,506)NDIM
         WRITE(6,507)NDIM*NATTOT
         WRITE(6,510)MAXXCO
         CALL STOPGM('MM_SETUP','FATAL ERROR')
      ENDIF


C print out total number of atoms in the system
      WRITE(6,*)'SYSTEMBLOCK'
      WRITE(6,500)NRPT
      WRITE(6,501)NSOLV
      WRITE(6,502)NATTOT
      WRITE(6,520)NPM*NCAG
      WRITE(6,521)NSM
      WRITE(6,522)NSM+NPM*NCAG
      WRITE(6,*)'END SYSTEMBLOCK'
      WRITE(6,*)

#ifndef READ_COORD
      CALL GTCOOR(PRGSTR,NATTOT,MXEWRT,DLSUM,DMSUM)
#endif
#include "inc_promd3.h"

C*****BUILD THE CONSTRAINTS FROM BOND LENGTHS IF NEEDED,
      IF (NTC .NE. NTCNON .OR. 
     $     (INIT .EQ. INSHVX .OR. INIT .EQ. INSHV)) THEN
         CALL SETCON
      ENDIF

C*****SET ALL CHARGES TO 0.0 IF SO REQUIRED.
      IF (.NOT. LTF(NTFCG)) THEN
         CALL mm_AZZERO(CG,NRP)
         CALL mm_AZZERO(CGS,NRAM)
      ENDIF
      WRITE(6,1140)

C now print out the submolecule atoms
      IF (NSPM .GT. 1) THEN
         KFRST = 1
         DO K=1,NSPM
            KLST = NSP(K)
            CALL SCOIN2(KFRST,KLST,MAXFLD,STR)
            WRITE(6,570)K,STR
            KFRST = KLST + 1
         ENDDO
         WRITE(6,*)
      ENDIF

      LFORMF = (NTXO .EQ. NTXOFO)
C******WE HAVE SET UP ALL THE DATA NEEDED. NOW DECIDE WHAT TO DO
      IF (IDOPRO .EQ. IDOEM) THEN
C lets do EM
!        CALL RUNEM(NRPT,NATTOT,NUM4DA)
      ELSEIF (IDOPRO .EQ. IDOSD .OR. IDOPRO .EQ. IDOMD) THEN
C lets do MD or SD
         IF (NSTLIM .GT. 0) THEN
#ifdef READ_COORD
            CALL DOMDSD(PRGSTR,NRPT,NSOLV,NATTOT,NUM4DA,IFNCFG,
     .                  IURDCO,MXETBL,DLSUM,DMSUM)
#else
            CALL DOMDSD(PRGSTR,NRPT,NSOLV,NATTOT,NUM4DA,IFNCFG,
     .                  MXETBL,DLSUM,DMSUM)
#endif
         ENDIF
      ELSE
         WRITE(6,*)PRGSTR,': illegal value for IDOPRO =',IDOPRO
         CALL STOPGM('MM_SETUP','FATAL ERROR')
      ENDIF


      CALL FILECLOSE(ISTDIN)
      CALL FILECLOSE(IOTOPO)
      CALL FILECLOSE(IOXVI)
#include "inc_promd5.h"

 500  FORMAT (4X,'NUMBER OF SOLUTE ATOMS         :',I8)
 501  FORMAT (4X,'NUMBER OF SOLVENT ATOMS        :',I8)
 502  FORMAT (4X,'TOTAL NUMBER OF ATOMS          :',I8)
 505  FORMAT (4X,'MAXIMUM ALLOWED:         MAXNAT=',I8)
 506  FORMAT (4X,'DIMENSION SPECIFIED            :',I8)
 507  FORMAT (4X,'COORDINATE ARRAY SIZE NEEDED   :',I8)
 510  FORMAT (4X,'MAXIMUM ALLOWED:         MAXXCO=',I8)
C
 520  FORMAT (4X,'NUMBER OF SOLUTE CHARGE GROUPS :',I8)
 521  FORMAT (4X,'NUMBER OF SOLVENT CHARGE GROUPS:',I8)
 522  FORMAT (4X,'TOTAL NUMBER OF CHARGE GROUPS  :',I8)
C
 550  FORMAT (1X,'DIHEDRALS TO MONITOR           :',I8)
 560  FORMAT (1X,'MAXIMUM ALLOWED          MXDMON=',I8)

 570  FORMAT (1X,'SUBMOLECULE ',I3,': ATOMS ',A11)
 580  FORMAT (1X,'TOTAL OF',I7,' ATOMS IN 4D')

C titles as from 1000
 1010 FORMAT (//,1X,'1. M O L E C U L A R   T O P O L O G Y',/)
 1110 FORMAT (//,1X,'2. D A T A   C H A R A C T E R I Z I N G',
     $     '   T H E   RUN',/)
 1120 FORMAT (//,1X,'3. O P T I O N A L   D A T A   F I L E S',/)

 1130 FORMAT (//,1X,'4. A T O M I C   C O O R D I N A T E S   A N D',
     $     '   V E L O C I T I E S',/)

 1140 FORMAT (//,1X,'5. I N I T I A L I Z E   R U N',/)

 1000 FORMAT (1X,' 1. A MOLECULAR TOPOLOGY FILE')

 1100 FORMAT (1X,' 2. DATA CHARACTERIZING THE RUN')

 1200 FORMAT (1X,' 3. DATA DETERMINING THE FREE ENERGY PERTURBATION ',
     $     '(IF NTG>0)')

 1210 FORMAT (//,1X,'D A T A   D E T E R M I N I N G   T H E',
     $     '   F R E E',
     $     '   E N E R G Y   P E R T U R B A T I O N',/)

C all spec files
 1300 FORMAT (1X,' 4. RE(CON)STRAINED ATOM SEQUENCE NUMBERS (IF NTR>0)')
 1310 FORMAT (//,1X,'R E (C O N) S T R A I N E D   A T O M',
     $     '   S E Q U E N C E   N U M B E R S',/)

 1400 FORMAT (1X,' 5. DISTANCE RESTRAINT ATOM PAIRS (IF NTDR>0)')
 1410 FORMAT (//,1X,'D I S T A N C E   R E S T R A I N I N G',/)

 1500 FORMAT (1X,' 6. RESTRAINED DIHEDRALS (IF NTDLR>0)')

 1600 FORMAT (1X,' 7. J-VALUE RESTRAINING DIHEDRALS (IF NTJR >0)')
 1610 FORMAT (//,1X,'J - V A L U E   R E S T R A I N I N G',/)

 1700 FORMAT (1X,' 8. LOCAL ELEVATION DIHEDRALS (IF NTLE > 0)')
 1710 FORMAT (//,1X,'L O C A L   E L E V A T I O N',/)

 1800 FORMAT (1X,' 9. 4th DIMENSION ATOM INDICATORS (IF NT4DIM >0)')
 1819 FORMAT (//,1X,'4 t h   D I M E N S I O N   D A T A',/)

C all data files
 1900 FORMAT (1X,'10. INITIAL ATOMIC COORDINATES AND VELOCITIES')

 2000 FORMAT (1X,'11. ATOMIC FRICTION COEFFICIENTS (IF NTFR=2)')

 2100 FORMAT (1X,'12. REFERENCE ATOM COORDINATES FOR RE(CON)STRAINING',
     $     ' (IF NTR>0)')
      RETURN
114   CALL STOPGM('mm_setup','gromos input file not found')
115   CALL STOPGM('mm_setup','topology file not found')
116   CALL STOPGM('mm_setup','groomos coordinates file not found')
C end PROMD
      END
C     ==================================================================
      SUBROUTINE DOMDSD(PRGSTR,NRPT,NSOLV,NATTOT,NUM4DA,IFNCFG,
C     ==--------------------------------------------------------------==
C     DOMDSD opens the required trajectory files and performs
C     a molecular dynamics or stochastic runs by calling
C     L<RUNMD>.
C     It then closes the trajectory files.
COMMEND
#ifdef READ_COORD
     $     IURDCO,
#endif
     $     MXETBL,DLSUM,DMSUM)
C includes
      INCLUDE 'coordsz.h'
      INCLUDE 'md.h'
C args
      CHARACTER PRGSTR*(*)
      INTEGER NRPT,NSOLV,NATTOT,NUM4DA,IFNCFG,MXETBL
      real*8 DLSUM(MXETBL),DMSUM(MXETBL)
C local params
      INTEGER MAXSTR
      PARAMETER (MAXSTR = 12)
C local vars
      CHARACTER FRMSTR*(MAXSTR)
      CHARACTER STATST*(MAXSTR)
C LFORM is .TRUE. if trajectories are written formatted
      LOGICAL LFORM
      real*8 TIME
C units for writing trajectories
      INTEGER IUTRJX,IUTRJV,IUTRJE,IUTRJG
C data
C with STATST = UNKNOWN, any existing trajectory files
C will be overwritten. Setting this to 'NEW' will prevent this.
      DATA STATST /'UNKNOWN'/
      DATA IUTRJX /-1/
      DATA IUTRJV /-1/
      DATA IUTRJE /-1/
      DATA IUTRJG /-1/

#include "inc_promd1.h"

C begin
C open trajectory output files if necessary
      LFORM = (NTPW .EQ. NTPWFO)
      IF (LFORM) THEN
         FRMSTR = 'FORMATTED'
      ELSE
         FRMSTR = 'UNFORMATTED'
      ENDIF

      TIME = T

#include "inc_promd2.h"

      IF (IDOPRO .EQ. IDOSD) THEN
         WRITE(6,610)'STOCHASTIC'
      ELSE
         WRITE(6,610)'MOLECULAR'
      ENDIF

C*****PERFORM MD or SD-RUN
      IF (LTCM) THEN
         WRITE(6,*)
     $    'OVERALL TRANSLATIONAL AND ROTATIONAL MOTION IS STOPPED,'
         WRITE(6,*)
      ENDIF

CCC
C     uncomment this here for benchmarking
C     and then make the target promd.tim to link in the
C     timing routines
C     CALL INITIM
CCC

C FIXME: AK 2005/04/17 fix this preprocessor monstrosity
C                      and have separate functions instead.
      CALL mm_setup_dr(NRPT,NSOLV,NATTOT,NUM4DA,TIME,
     $     DLSUM,DMSUM,
     $     IUTRJX,IUTRJV,IUTRJE,IUTRJG,
#ifdef WRITE_FORCES
     $     IUTRJF,
#endif
#ifdef READ_COORD
     $     IURDCO,
#endif
     $     IFNCFG)

CCC
C uncomment this here for benchmarking
C     CALL PRNTIM
CCC

#ifdef TIMING
      CALL PRNTIM()
      WRITE(6,*)'END PROMD'
#endif

 610  FORMAT (1X,'PERFORMING A ',A,' DYNAMICS RUN',/)
C end domdsd
C     ==--------------------------------------------------------------==
      RETURN
      END
